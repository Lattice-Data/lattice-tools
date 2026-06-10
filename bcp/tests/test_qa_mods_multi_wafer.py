"""Unit tests for multi-wafer trimmer key helpers in qa_mods."""

from qa_mods import resolve_wafer_run_id, trimmer_failure_storage_key


class TestTrimmerFailureStorageKeyMod:
    def test_multiple_wafers_same_sublibrary_get_distinct_keys(self):
        base = "proj/order/LIB01/raw"
        keys = [
            f"{base}/441001-LIB01_GEX-Z5001-BC01_trimmer-failure_codes.csv",
            f"{base}/441002-LIB01_CRI-Z5002-BC02_trimmer-failure_codes.csv",
            f"{base}/441003-LIB01_GEX-Z5003-BC03_trimmer-failure_codes.csv",
        ]
        storage_keys = [trimmer_failure_storage_key(k)[0] for k in keys]
        assert storage_keys == ["441001", "441002", "441003"]
        assert len(set(storage_keys)) == 3

    def test_legacy_experiment_key_still_resolves_via_map(self):
        assert (
            resolve_wafer_run_id("order/LIB01", {"order/LIB01": "441001"}) == "441001"
        )
