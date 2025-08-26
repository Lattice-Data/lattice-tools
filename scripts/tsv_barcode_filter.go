package main

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"log"
	"math"
	"os"
	"strings"
	"sync"
)

// filtered barcode counter needs lock for worker writes, could use this for raw
// counts too, but lock not needed
type Counter struct {
	sync.RWMutex
	m map[string]int
}

// trying go idiomatic enum for no col magic number
type FragmentCols int

const (
	ColChrom FragmentCols = iota
	ColStart
	ColEnd
	ColBarcode
	ColreadSupport
)

type Stats struct {
	// matching key names in concatenator
	RawMatrix      string `json:"raw_matrix"`
	RawMin         int    `json:"raw_min"`
	RawMean        int    `json:"raw_mean"`
	RawRowCount    int    `json:"raw_row_num"`
	FilterMin      int    `json:"filt_min"`
	FilterMean     int    `json:"filt_mean"`
	FilterRowCount int    `json:"filt_row_num"`
	UniqueBarcodes int    `json:"unique_barcodes"`
}

func makeBarcodeMap() map[string]struct{} {
	barcodeSet := make(map[string]struct{})

	f, err := os.Open("barcodes.txt")
	if err != nil {
		log.Println("Error opening barcodes file:", err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		barcodeSet[line] = struct{}{}
	}
	if err := scanner.Err(); err != nil {
		log.Println("Error with scanner:", err)
	}
	return barcodeSet
}

func worker(label string, location string, barcodeSet *map[string]struct{}, counter *Counter, lines <-chan []string, results chan<- []string, wg *sync.WaitGroup) {
	defer wg.Done()
	var testBarcode string
	for line := range lines {
		if location == "prefix" {
			testBarcode = label + line[ColBarcode]
		} else {
			testBarcode = line[ColBarcode] + label
		}
		// map should pass by reference by default, trying out pointer dereference
		if _, exists := (*barcodeSet)[testBarcode]; exists {
			line[ColBarcode] = testBarcode
			// processedResult := fmt.Sprintf("Worker %d processed: %v", id, line)
			// fmt.Println(processedResult)
			counter.Lock()
			counter.m[testBarcode]++
			counter.Unlock()
			results <- line
		}
	}
}

func main() {
	args := os.Args
	fmt.Println(args)
	label := args[1]
	location := args[2]
	filePath := args[3]
	fileName := strings.Split(filePath, "/")[1]
	accession := strings.Split(fileName, "_")[0]

	// with gzip decompress on the fly, 2 filter workers seems about right
	// straight tsv file needs 3-4 workers for increased read speeds
	numWorkers := 2

	file, err := os.Open(filePath)
	if err != nil {
		log.Println("Error opening file:", err)
		return
	}
	defer file.Close()

	gzReader, err := gzip.NewReader(file)
	if err != nil {
		log.Println("Error with gzip reader: ", err)
		return
	}
	defer gzReader.Close()

	reader := csv.NewReader(gzReader)
	reader.Comma = '\t'
	reader.Comment = '#'

	outputFile := strings.Replace(filePath, "fragments.tsv.gz", "filtered_fragments.tsv", 1)
	output, err := os.Create(outputFile)
	if err != nil {
		log.Fatalf("failed to create output file: %v", err)
	}
	defer output.Close()

	writer := csv.NewWriter(output)
	writer.Comma = '\t'
	// flush here instead of in receiving channel
	defer writer.Flush()

	linesChan := make(chan []string, 1000)
	resultsChan := make(chan []string)
	var wg sync.WaitGroup
	var rawRowCount int
	filtBarcodeCounter := Counter{m: make(map[string]int)}
	barcodeSet := makeBarcodeMap()

	for i := 1; i <= numWorkers; i++ {
		wg.Add(1)
		go worker(label, location, &barcodeSet, &filtBarcodeCounter, linesChan, resultsChan, &wg)
	}

	// hashmap for proper raw stats
	rawBarcodeCounts := make(map[string]int)

	// read lines and send to workers
	go func() {
		for {
			record, err := reader.Read()
			if err != nil {
				if err.Error() == "EOF" {
					break
				}
				log.Println("Error reading record:", err)
				break
			}
			rawRowCount++
			barcode := record[ColBarcode]
			rawBarcodeCounts[barcode]++
			linesChan <- record
		}
		close(linesChan) // close channel when all lines read
	}()

	// collect the results
	go func() {
		for result := range resultsChan {
			if err := writer.Write(result); err != nil {
				log.Fatal("Line write error: ", err)
			}
		}
		if err := writer.Error(); err != nil {
			log.Fatal("Writer flush error: ", err)
		}
	}()

	wg.Wait()
	close(resultsChan)

	var filterRowCount int
	uniqueFiltBarcodes := len(filtBarcodeCounter.m)
	uniqueRawBarcodes := len(rawBarcodeCounts)

	filtMinVal := math.MaxInt
	for _, counts := range filtBarcodeCounter.m {
		filterRowCount += counts
		if counts < filtMinVal {
			filtMinVal = counts
		}
	}
	filtMean := filterRowCount / uniqueFiltBarcodes

	rawMinVal := math.MaxInt
	for _, counts := range rawBarcodeCounts {
		if counts < rawMinVal {
			rawMinVal = counts
		}
	}
	rawMean := rawRowCount / uniqueRawBarcodes

	stats := Stats{
		RawMatrix:      accession,
		RawMin:         rawMinVal,
		RawMean:        rawMean,
		RawRowCount:    rawRowCount,
		FilterMin:      filtMinVal,
		FilterMean:     filtMean,
		FilterRowCount: filterRowCount,
		UniqueBarcodes: uniqueFiltBarcodes,
	}
	jsonStats, err := json.Marshal(stats)
	if err != nil {
		log.Fatal("Error with json marshal", err)
	}

	// using this to parse results, current logging shows command line args
	// plus final results
	fmt.Println("||", string(jsonStats))

}

// linux build: GOOS=linux GOARCH=amd64 go build -o tsv_barcode_filter_linux_amd64 tsv_barcode_filter.go
