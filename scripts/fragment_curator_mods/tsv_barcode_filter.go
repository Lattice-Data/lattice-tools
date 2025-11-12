package main

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"runtime/pprof"
	"runtime/trace"
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

// [][]string size for results channel
const batchSize = 1_000_000

type Stats struct {
	// matching key names in fragment curator
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

// will keep worker function for now, not currently used
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

// using specific flags for normal inputs, parsing returns pointers
var label = flag.String("label", "", "label to prepend or append to barcode")
var location = flag.String("location", "", "prefix or suffix location or attaching label")
var filePath = flag.String("filepath", "", "path to raw fragment file")

// optional profiling flags
// if profiling in parallel, should modify prof and out files to contain above info to not overwrite
var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
var traceprofile = flag.String("traceprofile", "", "write trace execution to `file`")

func main() {
	flag.Parse()

	if *cpuprofile != "" {
		f, err := os.Create("./profiles/cpu.prof")
		if err != nil {
			log.Fatal("could not create CPU profile")
		}
		defer f.Close()
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}

	if *traceprofile != "" {
		t, err := os.Create("./profiles/trace.out")
		if err != nil {
			log.Fatal("could not create trace execution profile")
		}
		defer t.Close()
		if err := trace.Start(t); err != nil {
			log.Fatal("could not start trace: ", err)
		}
		defer trace.Stop()
	}

	fileName := strings.Split(*filePath, "/")[1]
	accession := strings.Split(fileName, "_")[0]

	file, err := os.Open(*filePath)
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

	outputFile := strings.Replace(*filePath, "fragments.tsv.gz", "filtered_fragments.tsv", 1)
	output, err := os.Create(outputFile)
	if err != nil {
		log.Fatalf("failed to create output file: %v", err)
	}
	defer output.Close()

	writer := csv.NewWriter(output)
	writer.Comma = '\t'
	// flush here instead of in receiving channel, could also defer flush in writer goroutine
	defer writer.Flush()

	resultsChan := make(chan [][]string, 100)
	var rawRowCount int
	filtBarcodeCounter := Counter{m: make(map[string]int)}
	barcodeSet := makeBarcodeMap()

	var wg sync.WaitGroup
	// one for reader, one for writer
	wg.Add(2)

	// hashmap for proper raw stats
	rawBarcodeCounts := make(map[string]int)

	// read lines, send in batch over channel
	go func() {
		defer wg.Done()
		var batch [][]string
		for {
			var testBarcode string
			record, err := reader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				log.Println("Error reading record:", err)
				break
			}

			if *location == "prefix" {
				testBarcode = *label + record[ColBarcode]
			} else {
				testBarcode = record[ColBarcode] + *label
			}

			if _, exists := barcodeSet[testBarcode]; exists {
				record[ColBarcode] = testBarcode
				filtBarcodeCounter.m[testBarcode]++
				batch = append(batch, record)

				if len(batch) == batchSize {
					resultsChan <- batch
					batch = nil
				}
			}
			rawRowCount++
			barcode := record[ColBarcode]
			rawBarcodeCounts[barcode]++
		}
		// if remaining batch < batchSize, need to send to channel
		if len(batch) > 0 {
			resultsChan <- batch
		}
		// go idiom for sender to close channel
		// channel will close after channel buffer is cleared
		close(resultsChan)
	}()

	// write lines to filtered file, wrtie batchSize
	go func() {
		defer wg.Done()
		for batch := range resultsChan {
			if err := writer.WriteAll(batch); err != nil {
				log.Fatal("write error: ", err)
			}
		}
		if err := writer.Error(); err != nil {
			log.Fatal("Writer flush error: ", err)
		}
	}()

	wg.Wait()

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

// macOS build: go build -o tsv_barcode_filter_macos_arm64 tsv_barcode_filter.go
// linux build: GOOS=linux GOARCH=amd64 go build -o tsv_barcode_filter_linux_amd64 tsv_barcode_filter.go
