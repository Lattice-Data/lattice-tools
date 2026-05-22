# ruff linter/formatter will remove unused imports, need comments below to keep them
# since pytest arugments/usage shadow imports
import pytest
import sys
from fixtures.valid_adatas import (
    ALL_H5ADS,                  # noqa: F401
    SPATIAL_H5ADS,              # noqa: F401
    NON_SPATIAL_H5ADS,          # noqa: F401
    MULTISPECIES_H5ADS,         # noqa: F401
    GUIDE_H5ADS,                # noqa: F401
    test_h5ads,                 # noqa: F401
    validator_with_adatas,      # noqa: F401
    yield_guide_validator,      # noqa: F401
    yield_guide_file_name,      # noqa: F401
    _to_anndata_file,           # noqa: F401
    to_temp_files               # noqa: F401
)


@pytest.mark.parametrize("yield_guide_file_name", GUIDE_H5ADS)
class TestGeneticPerturbationPasses:
    @pytest.fixture(autouse=True)
    def setup(self, yield_guide_validator):
        self.validator = yield_guide_validator


    def test_guide_h5ads_passes(self):
        assert "genetic_perturbations" in self.validator.adata.uns
        assert len(self.validator.adata.uns["genetic_perturbations"]) > 1
        assert "genetic_perturbation_id" in self.validator.adata.obs.columns
        assert "genetic_perturbation_strategy" in self.validator.adata.obs.columns
        assert "control" in self.validator.adata.obs["genetic_perturbation_strategy"].unique()

        print(self.validator.adata.obs[["genetic_perturbation_id", "genetic_perturbation_strategy"]].head(), file=sys.__stderr__)
        print(self.validator.adata.obs[["genetic_perturbation_id", "genetic_perturbation_strategy"]].tail(), file=sys.__stderr__)

        self.validator.validate_adata()
        assert self.validator.is_valid
        assert self.validator.errors == []
