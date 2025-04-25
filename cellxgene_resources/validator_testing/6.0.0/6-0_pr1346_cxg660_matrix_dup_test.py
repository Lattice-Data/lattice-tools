"""
QA testing for this issue: https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/gh/chanzuckerberg/single-cell-curation/1317
PR: https://github.com/chanzuckerberg/single-cell-curation/pull/1346
"""

import pytest
from fixtures.valid_adatas import (
    ALL_H5ADS,
    test_h5ads,
    validator_with_adatas
)





# should not pass:
# - duplicate a cell’s raw counts once (raw is in .X / raw is in raw.X)
# - duplicate many cell's raw counts (raw is in .X / raw is in raw.X)
# - duplicate a cell's raw counts many times (raw is in .X / raw is in raw.X)
#
#
#
#