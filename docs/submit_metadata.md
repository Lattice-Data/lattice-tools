Running submit_metadata.py to post new object and patch existing objects to the Lattice DB
---------------- 
To run this script, you should first set-up your environment according to the [lattice-tools instructions](../README.md)

This script takes in an Excel file with the metadatadata
This is a dryrun-default script, run with --update to enable patch and post

Defining object types
---------------- 
Name each sheet of the excel file the name of the object type you are using
with the format used on https://www.lattice-data.org/profiles/,
capitalization and underscores do not matter

The objects will be loaded in the order specified in encoded/src/loadxl.py ORDER

Use the --justtype argument to only submit one of the object types, even if your file contains more sheets
Ex: %(prog)s mydata.xsls --justtype Experiment

Use the --starttype argument to start at an object type and only submit the sequential objects

Spreadsheet formatting
---------------- 
The top row of each sheet should be the names of the fields specified in the schema.
The validity of these fields will be check during a 'dry-run' (i.e. if '--update' is not used)

If the first value in a row (column A) begins with '#', the entire row will be ignored
This is useful if you want rows with property descriptors (e.g. description, enum, linkTo)
These rows to ignore must be above metadata-to-submit rows or the Row number print-outs will be inaccurate

If the first cell (A1) is either empty or begins with 'schema_version=', the entire first column (A) will be ignored
This action happens after the above avoidance of rows starting with '#'

Any sheet with the name 'Cover Sheet' will be ignored, capitalization does not matter, but the space does

Objects with attachments (not tested)
---------------- 
To upload objects with attachments (Documents), have a column titled 'attachment'containing the
path and name of the file you wish to attach

Patch existing objects
---------------- 
If there is an idenitifying property (e.g. uuid, alias, accession, name, or term_id) in the sheet,
and an existing object is found using that identifier, it will ask if you want to PATCH that object
Use '--patchall' if you want to patch ALL objects in your document and ignore that message

Ontology terms
---------------- 
In most cases, a file does not need a separate OntologyTerm sheet. Instead, the field that linkTo:OntologyTerm
can be split into 2 fields: <property>.term_name and <property>.term_id. The submission script will pair these
together and submit the corresponding OntologyTerm, if it is not already in the DB. If it is in the DB, it will
confirm that term_name, etc. from the sheet match the DB metadata, and error if there is any mismatch.
Instead of a column for biosample_ontology with values that linkTo OntologyTerm (UBERON_0002113, CL_0000057),
you can submit 2 columns biosample_ontology.term_id (UBERON:0002113, CL:0000057) and biosample_ontology.term_name (kidney, fibroblast).
The above will not work for arrays of linkTo:OntologyTerm, like diseases. For that, a column labeled diseases is
expected and the values need to be comma-separated in the format to identify the OntologyTerm object to linkTo (MONDO_0005015, MONDO_0005565)

Embedded objects
---------------- 
Embedded objects are assumed to be of the format of dictionary objects or a list of dictionary objects

If you are submitting just one dictionary object...
Ex:

	 'plate_barcode_details': [
			{
				'barcode': 'ATGCCGCCG',
				'plate_location': 'A1'
			}
		]
Formatting in the document should be as follows for the above example:
| plate_barcode_details.barcode | plate_barcode_details.plate_location |
|:--|:--|
| ATGCCGCCG | A1 |

If you are submitting a list of multiple dictionary objects...
Ex:

	 'plate_barcode_details': [
			{
				'barcode': 'ATGCCGCCG',
				'plate_location': 'A1'
			},
			{
				'barcode': 'TGAAACGAC',
				'plate_location': 'A2'
			}
		]

An identifier (number or letter) should be appended to the property names w/ '-' in order to group them appropriately
| plate_barcode_details-1.barcode | plate_barcode_details-1.plate_location | plate_barcode_details-2.barcode | plate_barcode_details-2.plate_location |
|:--|:--|:--|:--|
| ATGCCGCCG| A1 | TGAAACGAC | A2 |

For more complex cases of embedded objects within embedded objects, the same logic applies but there will be two properties to append identifiers to
For example, these columns would be expected for a MatrixFile that has a layer with 2 filtering_cutoffs and another layer with 1 filter_cutoff

	layers-1.label
	layers-1.filtering_cutoffs-1.cutoff_value
	layers-1.filtering_cutoffs-1.cutoff_units
	layers-1.filtering_cutoffs-2.cutoff_value
	layers-1.filtering_cutoffs-2.cutoff_units
	layers-2.label
	layers-2.filtering_cutoffs-1.cutoff_value
	layers-2.filtering_cutoffs-1.cutoff_units
