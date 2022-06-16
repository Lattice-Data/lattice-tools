Running submit_metadata.py to post new object and patch existing objects to the Lattice DB
----------------
To run this script, you should first set-up your environment according to the [lattice-tools instructions](../README.md)

This script takes in a Google Sheet identifier with the metadata

This is a dryrun-default script, run with `--update` to enable patch and post

Defining object types
----------------
Name each tab of the Google sheet the name of the object type you are using with the format used in [loadxl.py ORDER]. This is also the format in the url when viewing schema on the portal (ex: https://www.lattice-data.org/profiles/human_postnatal_donor). Any tab that is not named after a valid schema type will be ignored.

The objects will be loaded in the order specified in [loadxl.py ORDER].

Use the `--justtype` argument to only submit one of the object types, even if your file contains more valid tabs. Example:
```
python submit_metadata.py -m local <sheet_id> --justtype raw_sequence_file
```

Use the `--starttype` argument to start at an object type and only submit the sequential objects.

Spreadsheet formatting
----------------
The top row of each sheet should be the names of the fields specified in the schema. The validity of these fields will be check during a "dry-run" (i.e. if `--update` is not used)

Any row, except for row 1, with the first value (column A) beginning with `#` and any column with the header value (row 1) beginning with `#` will be ignored. This is useful if you want rows with property descriptors (e.g. description, enum, linkTo) or columns with additional notes, annotations, etc.

Objects with attachments
----------------
To upload objects with attachments (**Documents**), have a column titled `attachment` containing the path and name of the file you wish to attach

Patch existing objects
----------------
To change or add metadata to existing objects, the sheet must include an identifying property (e.g. `uuid`, `alias`, `accession`, `name`, or `term_id`) and `--patchall` must be passed in order to idnetify the object to update.

Posting new objects and patching existing objects must be split into separate tasks. If `--patchall` is not used and existing objects are found to conflict based on identifying properties, an error will be reported.

Removing properties of existing objects
----------------
To remove values from existing objects, rather than updating them, use `--remove`. The first column in the sheet should be an identifying property. The other column(s) should be labeled with the property you wish to remove. It does not matter if values are filled in those columns, they will be ignored.

If a default value is set in the schema for a given property, removing that property will reset it to the default value.

Ontology terms
----------------
In most cases, a sheet does not need a separate ontology_term tab. Instead, the field in other objects that linkTo:OntologyTerm can be split into 2 fields: `<property>.term_name` and `<property>.term_id`. The submission script will pair these together and submit the corresponding OntologyTerm if it is not already in the DB. If it is in the DB, it will confirm that the `term_name` from the sheet match the DB metadata, and error if there is any mismatch. The submission is forgiving and will allow for `:` or `_` delimiter in the term_id.

For example, instead of a column for `biosample_ontology` with values that linkTo OntologyTerm (UBERON_0002113, CL_0000057), you can submit 2 columns `biosample_ontology.term_id` (UBERON:0002113, CL:0000057) and `biosample_ontology.term_name` (kidney, fibroblast).

The above will not work for arrays of linkTo:OntologyTerm, like `diseases`. For that, a column labeled `diseases` is expected and the values need to be comma-separated in the format to identify the OntologyTerm object to linkTo (MONDO_0005015, MONDO_0005565) and any terms not already represented in the DB will need to be defined on an ontology_term tab. The submission is forgiving and will allow for `:` or `_` delimiter in the term_id.

Embedded objects
----------------
Embedded objects are assumed to be of the format of dictionary objects or a list of dictionary objects

If you are submitting just one dictionary object...
```
'plate_barcode_details': [
	{
		'barcode': 'ATGCCGCCG',
		'plate_location': 'A1'
	}
]
```
Formatting in the document should be as follows for the above example:
| plate_barcode_details.barcode | plate_barcode_details.plate_location |
|:--|:--|
| ATGCCGCCG | A1 |

If you are submitting multiple dictionary objects...
```
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
```
An identifier (number or letter) should be appended to the property names with '-' in order to group them appropriately...
| plate_barcode_details-1.barcode | plate_barcode_details-1.plate_location | plate_barcode_details-2.barcode | plate_barcode_details-2.plate_location |
|:--|:--|:--|:--|
| ATGCCGCCG| A1 | TGAAACGAC | A2 |

For more complex cases of embedded objects within embedded objects, the same logic applies but there will be two properties to append identifiers to. For example, these columns would be expected for a MatrixFile that has a layer with 2 filtering_cutoffs and another layer with 1 filter_cutoff...
```
layers-1.label
layers-1.filtering_cutoffs-1.cutoff_value
layers-1.filtering_cutoffs-1.cutoff_units
layers-1.filtering_cutoffs-2.cutoff_value
layers-1.filtering_cutoffs-2.cutoff_units
layers-2.label
layers-2.filtering_cutoffs-1.cutoff_value
layers-2.filtering_cutoffs-1.cutoff_units
```

[loadxl.py ORDER]: https://github.com/Lattice-Data/encoded/blob/dev/src/encoded/loadxl.py#L15
