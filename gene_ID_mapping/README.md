# Ensembl gene ID mapping across versions

Symbols/names are often used to discuss genes in research, however, gene symbols cannot be relied upon to consistently refer to the same biological gene because the names are dynamic. Thus, Ensembl gene IDs in GENCODE reference annotations are used as stable identifiers. However, Ensembl gene IDs inadvertently change from GENCODE version-to-version. Through review of versiont-to-version changes and discussion with Ensembl, we have determined a mechanism to map such cases.

To identify mappings, each GENCODE annotation version from V19 until the one immediately preceding the target version was compared to the target version as described below.

A mapping was made from any 'old' gene ID to a 'target version' gene ID if the 'old' gene ID was not present in the target version and:
- the two genes had the same genomic coordinates **and/or**
- the two genes had the same transcript IDs associated with them

A mapping was then negated if either:
- the 'old' gene ID mapped to more than one different 'target version' gene ID amongst all annotation versions, **or**
- any GENCODE annotation version contained both the 'old' gene ID and the 'target version' gene ID

For each annotation V22 and after, the primary annotation was used. For V19-21, GENCODE had not yet released a primary annotation so the comprehensive annotation from ALL regions was used.

The coordinates of the V19 genes were converted from hg19 to GRCh38 coordinates via [liftOver](https://genome.ucsc.edu/FAQ/FAQdownloads.html#liftOver).

Target version 44 was used because it is the pinned annotation version for CELLxGENE Discover starting with [schema 5.0.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.0.0/schema.md#required-gene-annotations).

The output `gene_map_v44.json` was generated via running 
```
python gene_mapping.py -v 44
```

Future updates will include mouse annotation comparisons, but preliminary exploration indicated that there are significantly fewer cases that need to be mapped in this manner.
