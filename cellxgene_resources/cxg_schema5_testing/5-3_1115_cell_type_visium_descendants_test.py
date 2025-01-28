"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell-curation/issues/1103
https://github.com/chanzuckerberg/single-cell-curation/pull/1115/
"""

import pytest
from fixtures.valid_adatas import validator_with_all_visiums


# ERROR: not added unless .validate_adata() called
CT_UNKNOWN_ERROR = "obs['cell_type_ontology_term_id'] must be 'unknown' when obs['assay_ontology_term_id'] is a descendant of 'EFO:0010961' (Visium Spatial Gene Expression) and uns['spatial']['is_single'] is True and in_tissue is 0."

@pytest.mark.parametrize(
    "assay_term,expected",
    (
        pytest.param("EFO:0010961", False, id="Visium Spatial Gene Expression"),
        pytest.param("EFO:0022858", True, id="Visium CytAssist Spatial Gene Expression V2"),
        pytest.param("EFO:0022860", True, id="Visium CytAssist Spatial Gene Expression, 11mm"),
        pytest.param("EFO:0022859", True, id="Visium CytAssist Spatial Gene Expression, 6.5mm"),
        pytest.param("EFO:0022857", True, id="Visium Spatial Gene Expression V1"),
        pytest.param("EFO:0030005", False, id="spatial transcriptomics by high-throughput sequencing"),
        pytest.param("EFO:0009922", False, id="10x 3' v3"),
    )
)
def test_visium_descendants(validator_with_all_visiums, assay_term, expected):
    validator = validator_with_all_visiums
    validator.adata.obs['assay_ontology_term_id'] = assay_term
    validator._is_visium_including_descendants() 
    assert validator.is_visium is expected


def test_visium_parent_with_error(validator_with_all_visiums):
    validator = validator_with_all_visiums
    validator.adata.obs['assay_ontology_term_id'] = "EFO:0010961"
    validator._is_visium_including_descendants() 
    assert validator.is_visium is False
    assert (
        "Invalid spatial assay. obs['assay_ontology_term_id'] must be a descendant of "
        "EFO:0010961 but NOT EFO:0010961 itself. "
    ) in validator.errors


parameters_visiums = (
    "assay_term",
    (
        pytest.param("EFO:0010961", id="Visium Spatial Gene Expression"),
        pytest.param("EFO:0022858", id="Visium CytAssist Spatial Gene Expression V2"),
        pytest.param("EFO:0022860", id="Visium CytAssist Spatial Gene Expression, 11mm"),
        pytest.param("EFO:0022859", id="Visium CytAssist Spatial Gene Expression, 6.5mm"),
        pytest.param("EFO:0022857", id="Visium Spatial Gene Expression V1"),
    )
)

@pytest.mark.parametrize(*parameters_visiums)
def test_in_tissue_zero_w_cell_type_partial(validator_with_all_visiums, assay_term):
    validator = validator_with_all_visiums
    random_cell_type = "CL:0001082"

    validator.adata.obs['assay_ontology_term_id'] = assay_term

    # set first index in obs to fail
    validator.adata.obs.iloc[0, validator.adata.obs.columns.get_loc("in_tissue")] = 0
    
    if random_cell_type not in validator.adata.obs["cell_type_ontology_term_id"].cat.categories:
        validator.adata.obs["cell_type_ontology_term_id"] = validator.adata.obs["cell_type_ontology_term_id"].cat.add_categories(random_cell_type)

    validator.adata.obs.iloc[0, validator.adata.obs.columns.get_loc("cell_type_ontology_term_id")] = random_cell_type

    # only partial test to check method works
    validator._validate_spatial_cell_type_ontology_term_id()

    assert CT_UNKNOWN_ERROR in validator.errors


@pytest.mark.parametrize(*parameters_visiums)
def test_in_tissue_zero_w_cell_type_full(validator_with_all_visiums, assay_term):
    my_keyword_dict ={
        "validator_with_all_visiums": "my file",
        "assay_term": "my file",
    }
    test_in_tissue_zero_w_cell_type_full(**my_keyword_dict)
    validator = validator_with_all_visiums
    random_cell_type = "CL:0001082"
    validator.adata.obs['assay_ontology_term_id'] = assay_term

    # set first index in obs to fail
    validator.adata.obs.iloc[0, validator.adata.obs.columns.get_loc("in_tissue")] = 0
    
    if random_cell_type not in validator.adata.obs["cell_type_ontology_term_id"].cat.categories:
        validator.adata.obs["cell_type_ontology_term_id"] = validator.adata.obs["cell_type_ontology_term_id"].cat.add_categories(random_cell_type)

    validator.adata.obs.iloc[0, validator.adata.obs.columns.get_loc("cell_type_ontology_term_id")] = random_cell_type

    validator.validate_adata()

    # full adata not yet valid due to beginning of 5.3.0 work
    assert not validator.is_valid

    # just check that specific error is present
    ERROR_WITH_PREFIX = "ERROR: " + CT_UNKNOWN_ERROR

    assert ERROR_WITH_PREFIX in validator.errors


# initial issue with categorical check, making sure fix can work with object, string, category, etc
@pytest.mark.parametrize(*parameters_visiums)
def test_assay_term_string(validator_with_all_visiums, assay_term):
    validator = validator_with_all_visiums
    validator.adata.obs["assay_ontology_term_id"] = assay_term
    assert validator.adata.obs["assay_ontology_term_id"].dtype == 'object'
    validator.validate_adata()
    assert not validator.is_valid
