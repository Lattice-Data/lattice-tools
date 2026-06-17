"""Tests for file_extract.retry."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

from botocore.exceptions import ClientError

from file_extract.retry import is_transient, retry_with_backoff


def _client_error(code: str, status: int | None = None) -> ClientError:
    response: dict = {"Error": {"Code": code, "Message": "boom"}}
    if status is not None:
        response["ResponseMetadata"] = {"HTTPStatusCode": status}
    return ClientError(response, "TestOp")


def test_is_transient_throttling() -> None:
    assert is_transient(_client_error("Throttling", 400))


def test_is_transient_503() -> None:
    assert is_transient(_client_error("ServiceUnavailable", 503))


def test_is_transient_no_such_key() -> None:
    assert not is_transient(_client_error("NoSuchKey", 404))


def test_is_transient_access_denied() -> None:
    assert not is_transient(_client_error("AccessDenied", 403))


def test_is_transient_runtime_error() -> None:
    assert not is_transient(RuntimeError("missing checksum"))


@patch("file_extract.retry.time.sleep")
def test_retry_succeeds_first_try(_sleep: MagicMock) -> None:
    result, err = retry_with_backoff(lambda x: x + 1, 1)
    assert result == 2
    assert err is None
    _sleep.assert_not_called()


@patch("file_extract.retry.time.sleep")
def test_retry_transient_then_success(_sleep: MagicMock) -> None:
    calls = {"n": 0}

    def flaky() -> str:
        calls["n"] += 1
        if calls["n"] < 2:
            raise _client_error("Throttling")
        return "ok"

    result, err = retry_with_backoff(flaky, retries=3)
    assert result == "ok"
    assert err is None
    assert _sleep.call_count == 1


@patch("file_extract.retry.time.sleep")
def test_retry_permanent_no_backoff(_sleep: MagicMock) -> None:
    def fail() -> None:
        raise _client_error("NoSuchKey", 404)

    result, err = retry_with_backoff(fail, retries=5)
    assert result is None
    assert "NoSuchKey" in (err or "")
    _sleep.assert_not_called()


@patch("file_extract.retry.time.sleep")
def test_retry_exhausts_attempts(_sleep: MagicMock) -> None:
    def always_throttle() -> None:
        raise _client_error("SlowDown")

    result, err = retry_with_backoff(always_throttle, retries=3)
    assert result is None
    assert err is not None
    assert _sleep.call_count == 2
