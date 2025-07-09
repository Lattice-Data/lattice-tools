"""
QA testing for this issue:
https://github.com/chanzuckerberg/single-cell/issues/721
https://github.com/chanzuckerberg/single-cell-curation/pull/1095/

Incorporation of feedback on this PR:
https://github.com/chanzuckerberg/single-cell-curation/pull/1284

For fragment validation, need location of h5ad and fragment file, so need
to provide temp file location of files before calling process-fragment

tmpfile is magic item for pytest that will use tmp dir for test then breakdown/
delete after

current fixture setup will by default use all fixtures in h5ads list
can specify specific fixtures with this:
@pytest.mark.parametrize("yield_atac_h5ads", [{list specific h5ads}])
should apply this to other fixtures

"""

import anndata as ad
import gc
import pandas as pd
import numpy as np
import pytest
import pyarrow
from cellxgene_schema.atac_seq import (
    check_anndata_requires_fragment,
    process_fragment,
    human_chromosome_by_length,
    mouse_chromosome_by_length
)
from pandas._libs.parsers import STR_NA_VALUES
from pathlib import Path
from fixtures.create_fixtures import Organism
from fixtures.valid_adatas import (
    ATAC_H5ADS,
    FIXTURES_ROOT,
    read_h5ad,
    validator_with_adatas,
    yield_atac_fixture_data,
    yield_atac_h5ads,
    _to_anndata_file,
    to_temp_files
)


def test_mock(yield_atac_fixture_data, tmpdir):
    """
    get fixture
    modify
    call to_temp_files()
    run process_fragment
    assert
    """
    test_data = yield_atac_fixture_data

    temp_files = to_temp_files(test_data, tmpdir)
    results = process_fragment(**temp_files)

    assert results == []


@pytest.mark.parametrize("test_h5ads", ATAC_H5ADS)
class TestPairedRawCounts:
    def test_paired_no_raw(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["is_primary_data"] = True
        validator.adata.obs["assay_ontology_term_id"] = "EFO:0030059"
        validator.adata.obs["suspension_type"] = "nucleus"
        for col in ["assay_ontology_term_id", "suspension_type"]:
            validator.adata.obs[col] = validator.adata.obs[col].astype("category") 
        del validator.adata.raw
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            "ERROR: Raw data is missing: there is only a normalized matrix in X and no raw.X"
        ) in validator.errors

    def test_paired_raw_in_X(self, validator_with_adatas):
        validator = validator_with_adatas
        validator.adata.obs["is_primary_data"] = True
        validator.adata.obs["assay_ontology_term_id"] = "EFO:0030059"
        validator.adata.obs["suspension_type"] = "nucleus"
        for col in ["assay_ontology_term_id", "suspension_type"]:
            validator.adata.obs[col] = validator.adata.obs[col].astype("category") 
        validator.adata.X = validator.adata.raw.X
        del validator.adata.raw
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    @pytest.mark.parametrize("assay_term", ["EFO:0022045", "EFO:0030007", "EFO:0008925", "EFO:0008904"])
    def test_unpaired_no_raw(self, validator_with_adatas, assay_term):
        validator = validator_with_adatas
        validator.adata.obs["is_primary_data"] = True
        validator.adata.obs["assay_ontology_term_id"] = assay_term
        validator.adata.obs["suspension_type"] = "nucleus"
        for col in ["assay_ontology_term_id", "suspension_type"]:
            validator.adata.obs[col] = validator.adata.obs[col].astype("category") 
        del validator.adata.raw
        validator.validate_adata()
        assert validator.is_valid
        assert validator.errors == []

    @pytest.mark.parametrize("assay_term", ["EFO:0009922", "EFO:0009899", "EFO:0009294"])
    def test_non_atac_no_raw(self, validator_with_adatas, assay_term):
        validator = validator_with_adatas
        validator.adata.obs["is_primary_data"] = True
        validator.adata.obs["assay_ontology_term_id"] = assay_term
        for col in ["assay_ontology_term_id", "suspension_type"]:
            validator.adata.obs[col] = validator.adata.obs[col].astype("category") 
        del validator.adata.raw
        validator.validate_adata()
        assert not validator.is_valid
        assert (
            "ERROR: Raw data is missing: there is only a normalized matrix in X and no raw.X"
        ) in validator.errors
        

class TestIsPrimaryData:
    def test_all_false(self, yield_atac_fixture_data, tmpdir):
        test_data = yield_atac_fixture_data

        test_data.adata.obs["is_primary_data"] = False

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Anndata.obs.is_primary_data must all be True." in results

    def test_mix_true_false(self, yield_atac_fixture_data, tmpdir):
        test_data = yield_atac_fixture_data

        half_point = test_data.adata.shape[0] // 2
        half_index = test_data.adata.obs.iloc[half_point].name
        test_data.adata.obs.loc[half_index:, "is_primary_data"] = True
        test_data.adata.obs.loc[:half_index, "is_primary_data"] = False

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Anndata.obs.is_primary_data must all be True." in results


class TestAssayTerms:
    """
    through the cli, process_fragment() will first run check_anndata_requires_fragment() and return
    early if there are errors, but need to specifically call check_anndata_requires_fragment()
    for pytest

    output is kind of confusing:
        paired (10x multiome) will return False
        unpaired (other atac) will return True
        non-atac raises ValueErrors

    will try standard fixture to get h5ad and fragment file, only will use h5ad in tests
    """
    # only 10x multiome should return False as a paired assay that does not need fragment
    def test_assay_paired_is_false(self, yield_atac_fixture_data, tmpdir):
        test_data = yield_atac_fixture_data

        test_data.adata.obs["assay_ontology_term_id"] = "EFO:0030059"

        temp_adata = _to_anndata_file(test_data, tmpdir)
        results = check_anndata_requires_fragment(temp_adata)

        assert results is False

    # all of these are descendants of scATAC-seq, except for 10x multiome
    @pytest.mark.parametrize("assay_term", ["EFO:0022045", "EFO:0030007", "EFO:0008925", "EFO:0008904"])
    def test_assay_paired_is_true(self, yield_atac_fixture_data, tmpdir, assay_term):
        test_data = yield_atac_fixture_data

        test_data.adata.obs["assay_ontology_term_id"] = assay_term

        temp_adata = _to_anndata_file(test_data, tmpdir)
        results = check_anndata_requires_fragment(temp_adata)

        assert results is True

    @pytest.mark.parametrize("first_term", ["EFO:0022045", "EFO:0030007", "EFO:0008925", "EFO:0008904"])
    @pytest.mark.parametrize("second_term", ["EFO:0022045", "EFO:0030007", "EFO:0008925", "EFO:0008904"])
    def test_assay_mixed_paired_is_true(self, yield_atac_fixture_data, tmpdir, first_term, second_term):
        test_data = yield_atac_fixture_data
        
        half_point = test_data.adata.shape[0] // 2
        half_index = test_data.adata.obs.iloc[half_point].name

        del test_data.adata.obs["assay_ontology_term_id"]
        test_data.adata.obs.loc[half_index:, "assay_ontology_term_id"] = first_term
        test_data.adata.obs.loc[:half_index, "assay_ontology_term_id"] = second_term

        temp_adata = _to_anndata_file(test_data, tmpdir)
        results = check_anndata_requires_fragment(temp_adata)

        assert results is True

    @pytest.mark.parametrize("assay_term", ["EFO:0022045", "EFO:0030007", "EFO:0008925", "EFO:0008904"])
    def test_assay_paired_and_unpaired_error(self, yield_atac_fixture_data, tmpdir, assay_term):
        test_data = yield_atac_fixture_data

        half_point = test_data.adata.shape[0] // 2
        half_index = test_data.adata.obs.iloc[half_point].name

        del test_data.adata.obs["assay_ontology_term_id"]
        test_data.adata.obs.loc[half_index:, "assay_ontology_term_id"] = "EFO:0030059"
        test_data.adata.obs.loc[:half_index, "assay_ontology_term_id"] = assay_term

        temp_adata = _to_anndata_file(test_data, tmpdir)
        with pytest.raises(ValueError, match="Anndata.obs.assay_ontology_term_id has mixed paired and unpaired assay terms."):
            _ = check_anndata_requires_fragment(temp_adata)

    @pytest.mark.parametrize("first_term", ["EFO:0030059", "EFO:0022045", "EFO:0030007", "EFO:0008925", "EFO:0008904"])
    @pytest.mark.parametrize("second_term", ["EFO:0009922", "EFO:0008930", "EFO:0022858"])
    def test_assay_mixed_atac_non_atac_error(self, yield_atac_fixture_data, tmpdir, first_term, second_term):
        test_data = yield_atac_fixture_data
        
        half_point = test_data.adata.shape[0] // 2
        half_index = test_data.adata.obs.iloc[half_point].name

        del test_data.adata.obs["assay_ontology_term_id"]
        test_data.adata.obs.loc[half_index:, "assay_ontology_term_id"] = first_term
        test_data.adata.obs.loc[:half_index, "assay_ontology_term_id"] = second_term

        temp_adata = _to_anndata_file(test_data, tmpdir)
        with pytest.raises(ValueError, match="Anndata.obs.assay_ontology_term_id are not all descendants of EFO:0010891."):
            _ = check_anndata_requires_fragment(temp_adata)

    @pytest.mark.parametrize("assay_term", ["EFO:0009922", "EFO:0008930", "EFO:0022858"])
    def test_assay_not_atac_error(self, yield_atac_fixture_data, tmpdir, assay_term):
        test_data = yield_atac_fixture_data

        test_data.adata.obs["assay_ontology_term_id"] = assay_term

        temp_adata = _to_anndata_file(test_data, tmpdir)
        with pytest.raises(ValueError, match="Anndata.obs.assay_ontology_term_id are not all descendants of EFO:0010891."):
            _ = check_anndata_requires_fragment(temp_adata)


@pytest.mark.parametrize("organism_term", [organism.term_id for organism in Organism])
class TestOrganismTerms:
    def test_all_wrong_organism(self, yield_atac_fixture_data, tmpdir, organism_term):
        test_data = yield_atac_fixture_data

        test_data.adata.obs["organism_ontology_term_id"] = organism_term

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert (
            "Anndata.obs.organism_ontology_term_id must be one of "
            f"['NCBITaxon:9606', 'NCBITaxon:10090']. Got {organism_term}."
        ) in results


    def test_mixed_wrong_organism(self, yield_atac_fixture_data, tmpdir, organism_term):
        test_data = yield_atac_fixture_data

        original_org_term = test_data.adata.obs["organism_ontology_term_id"].unique()[0]
        half_point = test_data.adata.shape[0] // 2
        half_index = test_data.adata.obs.iloc[half_point].name

        test_data.adata.obs["organism_ontology_term_id"] = \
        test_data.adata.obs["organism_ontology_term_id"].cat.add_categories(organism_term)
        test_data.adata.obs.loc[half_index:, "organism_ontology_term_id"] = organism_term

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert (
            "Anndata.obs.organism_ontology_term_id must have exactly 1 unique value. "
            f"Found the following values:\n{original_org_term}\n\t{organism_term}"
        ) in results

    def test_mixed_human_mouse_organism(self, yield_atac_fixture_data, tmpdir, organism_term):
        # all methods in class expcet organism_term with class decorator, throw it away with _ assignment 
        _ = organism_term
        test_data = yield_atac_fixture_data

        organism_terms = ["NCBITaxon:9606", "NCBITaxon:10090"]

        for organism_term in organism_terms:
            if organism_term not in test_data.adata.obs["organism_ontology_term_id"].cat.categories:
                test_data.adata.obs["organism_ontology_term_id"] = \
                test_data.adata.obs["organism_ontology_term_id"].cat.add_categories(organism_term)

        half_point = test_data.adata.shape[0] // 2
        half_index = test_data.adata.obs.iloc[half_point].name
        test_data.adata.obs.loc[half_index:, "organism_ontology_term_id"] = organism_terms[0]
        test_data.adata.obs.loc[:half_index, "organism_ontology_term_id"] = organism_terms[1]

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert (
            "Anndata.obs.organism_ontology_term_id must have exactly 1 unique value. "
            f"Found the following values:\n{organism_terms[1]}\n\t{organism_terms[0]}"
        ) in results


class TestFragmentCol1Chr:
    @pytest.mark.parametrize("yield_atac_h5ads", ["valid_human.h5ad"])
    @pytest.mark.parametrize("other_chromosome", ["GL456210.1", "JH584295.1"])
    def test_mouse_chromosomes_in_human(self, yield_atac_fixture_data, tmpdir, other_chromosome):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[0, 0] = other_chromosome

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert (
            f"Chromosomes in the fragment do not match the organism(NCBITaxon:9606).\n{other_chromosome}"
        ) in results

    @pytest.mark.parametrize("yield_atac_h5ads", ["valid_mouse.h5ad"])
    @pytest.mark.parametrize("other_chromosome", ["KI270442.1", "GL000009.2"])
    def test_human_chromosomes_in_mouse(self, yield_atac_fixture_data, tmpdir, other_chromosome):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[0, 0] = other_chromosome

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert (
            f"Chromosomes in the fragment do not match the organism(NCBITaxon:10090).\n{other_chromosome}"
        ) in results


class TestFragmentCol2Start:
    # error message makes sense, will see if this is unclear in practice
    @pytest.mark.parametrize("yield_atac_h5ads", ["valid_human.h5ad"])
    @pytest.mark.parametrize(
        "chromosome_lengths", 
        [(chromosome, max) for chromosome, max in human_chromosome_by_length.items()]
    )
    @pytest.mark.parametrize("row_to_change", [0, -1])
    def test_over_max_chr_start_human(self, yield_atac_fixture_data, tmpdir, chromosome_lengths, row_to_change):
        test_data = yield_atac_fixture_data
        chromosome, length = chromosome_lengths

        test_data.fragment_df.iloc[row_to_change, 0] = chromosome
        test_data.fragment_df.iloc[row_to_change, 1] = length + 1
        test_data.fragment_df.iloc[row_to_change, 2] = length

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Stop coordinate must be greater than start coordinate." in results


    @pytest.mark.parametrize("yield_atac_h5ads", ["valid_human.h5ad"])
    @pytest.mark.parametrize(
        "chromosome_lengths", 
        [(chromosome, max) for chromosome, max in human_chromosome_by_length.items()]
    )
    @pytest.mark.parametrize("row_to_change", [0, -1])
    def test_over_max_chr_start_and_stop_human(self, yield_atac_fixture_data, tmpdir, chromosome_lengths, row_to_change):
        test_data = yield_atac_fixture_data
        chromosome, length = chromosome_lengths

        test_data.fragment_df.iloc[row_to_change, 0] = chromosome
        test_data.fragment_df.iloc[row_to_change, 1] = length + 1
        test_data.fragment_df.iloc[row_to_change, 2] = length + 1

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Stop coordinate must be greater than start coordinate." in results
        assert "Stop coordinate must be less than the chromosome length." in results

    @pytest.mark.parametrize("yield_atac_h5ads", ["valid_mouse.h5ad"])
    @pytest.mark.parametrize(
        "chromosome_lengths", 
        [(chromosome, max) for chromosome, max in mouse_chromosome_by_length.items()]
    )
    @pytest.mark.parametrize("row_to_change", [0, -1])
    def test_over_max_chr_start_mouse(self, yield_atac_fixture_data, tmpdir, chromosome_lengths, row_to_change):
        test_data = yield_atac_fixture_data
        chromosome, length = chromosome_lengths

        test_data.fragment_df.iloc[row_to_change, 0] = chromosome
        test_data.fragment_df.iloc[row_to_change, 1] = length + 1
        test_data.fragment_df.iloc[row_to_change, 2] = length

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Stop coordinate must be greater than start coordinate." in results


    @pytest.mark.parametrize("yield_atac_h5ads", ["valid_mouse.h5ad"])
    @pytest.mark.parametrize(
        "chromosome_lengths", 
        [(chromosome, max) for chromosome, max in mouse_chromosome_by_length.items()]
    )
    @pytest.mark.parametrize("row_to_change", [0, -1])
    def test_over_max_chr_start_and_stop_mouse(self, yield_atac_fixture_data, tmpdir, chromosome_lengths, row_to_change):
        test_data = yield_atac_fixture_data
        chromosome, length = chromosome_lengths

        test_data.fragment_df.iloc[row_to_change, 0] = chromosome
        test_data.fragment_df.iloc[row_to_change, 1] = length + 1
        test_data.fragment_df.iloc[row_to_change, 2] = length + 1

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Stop coordinate must be greater than start coordinate." in results
        assert "Stop coordinate must be less than the chromosome length." in results

    @pytest.mark.parametrize("yield_atac_h5ads", ["valid_human.h5ad"])
    @pytest.mark.parametrize(
        "chromosome", 
        [chromosome for chromosome in human_chromosome_by_length]
    )
    def test_start_under_one_human(self, yield_atac_fixture_data, tmpdir, chromosome):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[0, 0] = chromosome
        test_data.fragment_df.iloc[0, 1] = 0

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Start coordinate must be greater than 0." in results

    @pytest.mark.parametrize("yield_atac_h5ads", ["valid_mouse.h5ad"])
    @pytest.mark.parametrize(
        "chromosome", 
        [chromosome for chromosome in mouse_chromosome_by_length]
    )
    def test_start_under_one_mouse(self, yield_atac_fixture_data, tmpdir, chromosome):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[0, 0] = chromosome
        test_data.fragment_df.iloc[0, 1] = 0

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Start coordinate must be greater than 0." in results


class TestFragmentCol3Stop:
    def test_stop_equals_start(self, yield_atac_fixture_data, tmpdir):
        test_data = yield_atac_fixture_data

        start_value = test_data.fragment_df.iloc[0, 1]
        test_data.fragment_df.iloc[0, 2] = start_value

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Stop coordinate must be greater than start coordinate." in results

    def test_stop_less_than_start(self, yield_atac_fixture_data, tmpdir):
        test_data = yield_atac_fixture_data

        start_value = test_data.fragment_df.iloc[0, 1]
        test_data.fragment_df.iloc[0, 2] = start_value - 1

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Stop coordinate must be greater than start coordinate." in results

    @pytest.mark.parametrize("yield_atac_h5ads", ["valid_human.h5ad"])
    @pytest.mark.parametrize(
        "chromosome_lengths", 
        [(chromosome, max) for chromosome, max in human_chromosome_by_length.items()]
    )
    def test_stop_over_chr_max_human(self, yield_atac_fixture_data, tmpdir, chromosome_lengths):
        test_data = yield_atac_fixture_data
        chromosome, length = chromosome_lengths

        test_data.fragment_df.iloc[0, 0] = chromosome
        test_data.fragment_df.iloc[0, 2] = length + 1

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Stop coordinate must be less than the chromosome length." in results

    @pytest.mark.parametrize("yield_atac_h5ads", ["valid_mouse.h5ad"])
    @pytest.mark.parametrize(
        "chromosome_lengths", 
        [(chromosome, max) for chromosome, max in mouse_chromosome_by_length.items()]
    )
    def test_stop_over_chr_max_mouse(self, yield_atac_fixture_data, tmpdir, chromosome_lengths):
        test_data = yield_atac_fixture_data
        chromosome, length = chromosome_lengths

        test_data.fragment_df.iloc[0, 0] = chromosome
        test_data.fragment_df.iloc[0, 2] = length + 1

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Stop coordinate must be less than the chromosome length." in results

    def test_stop_equals_zero(self, yield_atac_fixture_data, tmpdir):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[0, 2] = 0

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Stop coordinate must be greater than start coordinate." in results


class TestFragmentCol4Barcodes:
    """
    Currently no check if a barcode is a null value; this will be valid
        validation looks for set equivalence of fragment barcodes and obs barcodes
        dask.count() ignores NA values, so not added to fragment barcode set

    Unclear if this will be an issue
    To curate, we read in the obs.index from a h5ad; AnnData does not allow null values
    as an obs index. The curation process might introduce a null value

    Might need to be explicit with read and write to tsv to keep strings literal 

    Some of the tests show that by default pandas will coerce certain stings to null;
    other similar strings will remain as their string literal
    """
    @pytest.mark.parametrize("row_to_add", [0, -1])
    def test_barcode_in_fragment_not_in_adata(self, yield_atac_fixture_data, tmpdir, row_to_add):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[row_to_add, 3] = "barcode not in obs"

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Barcodes don't match anndata.obs.index" in results

    def test_barcode_in_adata_not_in_fragment(self, yield_atac_fixture_data, tmpdir):
        test_data = yield_atac_fixture_data

        test_data.adata.obs.reset_index(inplace=True)
        test_data.adata.obs["index"][0] = "barcode not in fragment"
        test_data.adata.obs.set_index("index", inplace=True)

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Barcodes don't match anndata.obs.index" in results

    @pytest.mark.parametrize("row_to_add", [0, -1])
    def test_all_barcode_not_in_adata(self, yield_atac_fixture_data, tmpdir, row_to_add):
        test_data = yield_atac_fixture_data
        
        barcode_to_drop = test_data.fragment_df.iloc[row_to_add, 3]
        test_data.fragment_df = test_data.fragment_df[~test_data.fragment_df[3].isin([barcode_to_drop])]

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Barcodes don't match anndata.obs.index" in results

    def test_all_barcode_not_in_fragment(self, yield_atac_fixture_data, tmpdir):
        test_data = yield_atac_fixture_data
        
        num_cells = test_data.adata.shape[0]
        test_data.adata = test_data.adata[:num_cells - 1, :].copy()

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Barcodes don't match anndata.obs.index" in results

    @pytest.mark.parametrize("null_string", ["", "na", "null", "NA", "NaN", "Na", "none", "None", "pd.NA", "np.nan", "np.NaN"])
    def test_null_strings_in_obs_index(self, yield_atac_fixture_data, tmpdir, null_string):
        test_data = yield_atac_fixture_data
        
        test_data.adata.obs.reset_index(inplace=True)
        test_data.adata.obs["index"][0] = null_string
        test_data.adata.obs.set_index("index", inplace=True)

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Barcodes don't match anndata.obs.index" in results

    # null values in obs will remain string literals, so barcodes will not match
    @pytest.mark.parametrize("row_to_change", [0, -1])
    @pytest.mark.parametrize("null_string", ["", "na", "null", "NA", "NaN", "Na", "none", "None", "pd.NA", "np.nan", "np.NaN"])
    def test_null_strings_in_obs_and_fragment_barcode(self, yield_atac_fixture_data, tmpdir, null_string, row_to_change):
        test_data = yield_atac_fixture_data
        
        test_data.adata.obs.reset_index(inplace=True)
        test_data.adata.obs["index"][0] = null_string
        test_data.adata.obs.set_index("index", inplace=True)

        test_data.fragment_df.iloc[row_to_change, 3] = null_string

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Barcodes don't match anndata.obs.index" in results


class TestFragmentCol5ReadSupport:
    @pytest.mark.parametrize("values", [0, -1])
    @pytest.mark.parametrize("row_to_add", [0, -1])
    def test_less_than_one_in_read_support(self, yield_atac_fixture_data, tmpdir, values, row_to_add):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[row_to_add, 4] = values

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Read support must be greater than 0." in results


@pytest.mark.parametrize("row_to_change", [0, -1])
class TestFragmentColDtypes:
    """
    The error reporting for wrong dtypes in the fragment file is a little sloppy, but 
    probably managable to figure out what went wrong. 
    
    Overall, will need to see how this plays out with real datasets

    The validator will roughly catch exceptions in the error message but through the 
    cli, the traceback will also be dumped

    The chr column should be categorical; a pyarrow exception is raised during parquet conversion
    that will roughly point to 'chromosome' column parsing incorrectly

    The int columns raise a pandas exception that is better caught by the validator and lists
    more information about wrong dtype for a particular column

    The barcodes column will coerce any value to string or null:
        string literals will cause a validation error for mismatched barcodes
        strings coerced to null will be valid.
    """
    # col 0 and 3 need special tests due to dtype of category and str respectively
    @pytest.mark.parametrize("values", [2.3, True, "string", False])
    @pytest.mark.parametrize("column", [1, 2, 4])
    def test_different_dtypes_in_int_cols(self, yield_atac_fixture_data, tmpdir, values, column, row_to_change):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[row_to_change, column] = values

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)
        
        error_prefix = f"Error Parsing the fragment file. Check that columns match schema definition. Error:"
        assert any(s.startswith(error_prefix) for s in results)

    @pytest.mark.parametrize("values", [None, np.nan, pd.NA])
    @pytest.mark.parametrize("column", [1, 2, 4])
    def test_null_types_in_int_cols(self, yield_atac_fixture_data, tmpdir, values, column, row_to_change):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[row_to_change, column] = values

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        error_prefix = f"Error Parsing the fragment file. Check that columns match schema definition. Error:"
        assert any(s.startswith(error_prefix) for s in results)

    @pytest.mark.parametrize("values", ["", "na", "null", "NA", "NaN", ])
    @pytest.mark.parametrize("column", [1, 2])
    def test_null_strings_in_start_stop_cols(self, yield_atac_fixture_data, tmpdir, values, column, row_to_change):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[row_to_change, column] = values

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        # null strings produce variety of errors, harder to predict
        assert len(results) > 0

    @pytest.mark.parametrize("values", ["", "na", "null", "NA", "NaN", ])
    @pytest.mark.parametrize("column", [4])
    def test_null_strings_in_read_support_cols(self, yield_atac_fixture_data, tmpdir, values, column, row_to_change):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[row_to_change, column] = values

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        # null strings for read support now is apparently valid with pyarrow change, except for 'na'
        assert results == []

    @pytest.mark.parametrize("values", [2.3, True, "string", False, None, np.nan, pd.NA, "", "na", "null", "NA", "NaN"])
    @pytest.mark.parametrize("column", [0])
    def test_different_dtypes_in_chr_col(self, yield_atac_fixture_data, tmpdir, values, column, row_to_change):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[row_to_change, column] = values
        organism = test_data.adata.obs["organism_ontology_term_id"].unique()[0]

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        error_prefix = f"Chromosomes in the fragment do not match the organism({organism})"
        assert any(s.startswith(error_prefix) for s in results)

    # strings coerced to null will be valid
    # STR_NA_VALUES is an unordered set, xdist needs ordered interable to distribute tests to workers
    # sorted() returns sorted list to allow for parallel test running
    @pytest.mark.parametrize("null_string", sorted(STR_NA_VALUES))
    @pytest.mark.parametrize("column", [3])
    def test_null_strings_in_barcode_col_pass(self, yield_atac_fixture_data, tmpdir, null_string, column, row_to_change):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[row_to_change, column] = null_string

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Barcodes don't match anndata.obs.index" in results

    # strings that look like null but are kept literal will cause mismatched barcodes
    @pytest.mark.parametrize("null_string", ["Na", "na", "none", "pd.NA", "np.nan", "np.NaN"])
    @pytest.mark.parametrize("column", [3])
    def test_null_strings_in_barcode_col_fails(self, yield_atac_fixture_data, tmpdir, null_string, column, row_to_change):
        test_data = yield_atac_fixture_data

        test_data.fragment_df.iloc[row_to_change, column] = null_string

        temp_files = to_temp_files(test_data, tmpdir)
        results = process_fragment(**temp_files)

        assert "Barcodes don't match anndata.obs.index" in results
