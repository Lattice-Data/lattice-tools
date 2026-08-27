from functools import cache

from .constants import DEFAULT_MODE
from db2_flattener.gather.lattice import Connection


@cache
def get_connection(mode: str = DEFAULT_MODE) -> Connection:
    """One Connection per mode, built on first use"""
    return Connection(mode)


def configure_connection(mode: str = DEFAULT_MODE, **attributes) -> Connection:
    """Set attributes on connection if needed"""
    connection = get_connection(mode)
    for attribute, value in attributes.items():
        setattr(connection, attribute, value)
    return connection
