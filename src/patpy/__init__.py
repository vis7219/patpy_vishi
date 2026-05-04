from importlib.metadata import version

from . import datasets as dt, pl, pp, tl

__all__ = ["dt", "pl", "pp", "tl"]

__version__ = version("patpy")
