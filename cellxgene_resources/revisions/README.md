# Revision specs

One JSON file per Collection revision, named after the Jira ticket. A curator writes it
after reading the paper; `curation_revision_perturbation.ipynb` reads it and does the
mechanical part. The code never invents a term or picks an obs column.

```json
{
  "collection_id": "<published Collection id>",
  "ticket": "CXG-###",
  "paper": "<doi, optional>",
  "note": "<what was done to which samples, optional>",
  "datasets": {
    "<dataset id>": {
      "source_column": "<author obs column that marks exposure>",
      "terms": {
        "<value in that column>": ["CHEBI:...", "EFO:...", "uniprot:...", "anti-uniprot:..."]
      }
    }
  }
}
```

* Datasets not listed are left as published.
* Observations whose source value is not listed get `"na"`.
* Several values may map to different term lists; a list with several terms is written
  sorted and `" || "`-delimited, as the schema requires.
* Terms are checked against schema 7.1.0 before any file is touched: allowed prefixes,
  ChEBI under chemical entity and outside the forbidden list, EFO limited to temperature,
  diet and its descendants, or fasting, UniProt accession format.
