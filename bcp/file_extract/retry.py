from __future__ import annotations

import random
import time
from typing import Callable, TypeVar

from botocore.exceptions import (
    ClientError,
    ConnectTimeoutError,
    ConnectionClosedError,
    EndpointConnectionError,
    ReadTimeoutError,
)

from .constants import TRANSIENT_ERROR_CODES, TRANSIENT_HTTP_STATUSES

T = TypeVar("T")


def is_transient(exc: BaseException) -> bool:
    """True if the exception is worth retrying; False for deterministic failures."""
    if isinstance(exc, ClientError):
        err = exc.response.get("Error", {}) or {}
        code = str(err.get("Code", ""))
        status = (exc.response.get("ResponseMetadata", {}) or {}).get("HTTPStatusCode")
        if code in TRANSIENT_ERROR_CODES:
            return True
        if status in TRANSIENT_HTTP_STATUSES:
            return True
        return False
    if isinstance(
        exc,
        (
            EndpointConnectionError,
            ConnectTimeoutError,
            ReadTimeoutError,
            ConnectionClosedError,
        ),
    ):
        return True
    return False


def retry_with_backoff(
    fn: Callable[..., T],
    *args,
    retries: int = 5,
    base_delay: float = 1.0,
    **kwargs,
) -> tuple[T | None, str | None]:
    """Exponential-backoff retry for transient errors only.

    Returns (result, error_str). On success error_str is None.
    """
    last_err: str | None = None
    for attempt in range(retries):
        try:
            return fn(*args, **kwargs), None
        except Exception as e:
            last_err = str(e)
            if attempt < retries - 1 and is_transient(e):
                time.sleep(base_delay * (2**attempt) + random.uniform(0, 0.5))
            else:
                break
    return None, last_err
