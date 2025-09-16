import anndata as ad


# add more columns that neeed na values
NA_COLUMNS = [
    "development_stage_ontology_term_id",
    "sex_ontology_term_id",
]


def make_valid_cell_line_fixture(adata: ad.AnnData, na_columns: list[str]=NA_COLUMNS) -> ad.AnnData:
    """
    Helper function to create valid h5ad when tissue_type == "cell line"
    Schema 7.0 introduces several changes where these columns MUST BE na:
        - development_stage_ontology_term_id
        - sex_ontology_term_id
        - self_reported_ethnicity_ontology_term_id
        - donor_id
    Update this function as features land to single-cell-curation instead of 
    needing to update all the tests individually

    na_columns = list of obs columns to set to "na", default is NA_COLUMNS in this file

    :rtype AnnData with select NA_COLUMNS set to "na"
    """

    adata.obs["tissue_type"] = "cell line"
    adata.obs["tissue_type"] = adata.obs["tissue_type"].astype("category")
    # random Cellosuarus term to be vaild
    adata.obs['tissue_ontology_term_id'] = 'CVCL_2830'

    for column in na_columns:
        adata.obs[column] = "na"

    return adata
