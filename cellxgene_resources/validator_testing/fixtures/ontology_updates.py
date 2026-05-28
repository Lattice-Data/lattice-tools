from dataclasses import dataclass, fields
from typing import Any
from .valid_adatas import DELIMITER


@dataclass
class ColumnOntologies:
    """
    Class to hold valid ontologies per organism, check for inclusion in attirbute set
    """
    assay: set[str]
    cell_type: set[str]
    development_stage: set[str]
    disease: set[str]
    experimental_condition: set[str]
    organism: set[str]
    self_reported_ethnicity: set[str]
    sex: set[str]
    tissue: set[str]

    @property
    def full_ontology_column_mapping(self) -> dict[str, set[str]]:
        return {
            f"{attribute.name}_ontology_term_id": getattr(self, attribute.name)
            for attribute 
            in fields(self)
        }
    

# for multiple ontologies per species, use delimiter to list each ontology
ONTOLOGY_MAP = {
    "assay": "EFO",
    "cell_type": {
        "NCBITaxon:6239": "WBbt || CL",         # C. Elegans
        "NCBITaxon:7227": "FBbt || CL",         # Drosophila
        "NCBITaxon:7955": "ZFA || CL",          # Zebrafish
        "default": "CL",                        # For all other organisms, use CL
    },
    "development_stage": {
        "NCBITaxon:6239": "WBls",               # C. Elegans
        "NCBITaxon:7227": "FBdv",               # Drosophila
        "NCBITaxon:7955": "ZFS",                # Zebrafish
        "NCBITaxon:9606": "HsapDv",             # Human
        "NCBITaxon:10090": "MmusDv",            # Mouse
        "default": "MONDO",                     # For all other organisms, use MONDO
    },
    "experimental_condition": "EFO || CHEBI || uniprot || anti-uniprot",
    "disease": "MONDO || PATO",
    "organism": "NCBITaxon",
    "self_reported_ethnicity": {
        "NCBITaxon:9606": "AfPO || HANCESTRO", 
        "default": None
    },
    "sex": "PATO",
    "tissue": {
        "NCBITaxon:6239": "WBbt || UBERON || CL || CVCL",     # C. Elegans
        "NCBITaxon:7227": "FBbt || UBERON || CL || CVCL",     # Drosophila
        "NCBITaxon:7955": "ZFA || UBERON || CL || CVCL",      # Zebrafish
        "default": "UBERON || CL || CVCL",                    # For all other organisms, use UBERON
    },
}


def get_valid_ontology_mappings(organism: str, ontology_map: dict[str, Any] = ONTOLOGY_MAP) -> ColumnOntologies:
    """
    Take organism ontology term and return valid ontologies per field
    """
    result: dict[str, set[str]] = {}
    for ontology_column, column_mapping in ontology_map.items():
        match column_mapping:
            case str(column_mapping):
                result[ontology_column] = set(column_mapping.split(DELIMITER))
            case dict(column_mapping):
                if organism in column_mapping:
                    result[ontology_column] = set(column_mapping[organism].split(DELIMITER))
                else:
                    result[ontology_column] = set([column_mapping["default"]])

    return ColumnOntologies(**result)
