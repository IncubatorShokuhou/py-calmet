"""py-calmet: pure-NumPy CALMET-style diagnostic downscaling."""

__version__ = "0.1.0"

from .core.runner import run_calmet, CalmetResult
from .io.calmet_dat import read_calmet_dat, write_calmet_dat, write_calmet_netcdf, CalmetDataset

__all__ = [
    "__version__",
    "run_calmet",
    "CalmetResult",
    "read_calmet_dat",
    "write_calmet_dat",
    "write_calmet_netcdf",
    "CalmetDataset",
]
