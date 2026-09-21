from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
# Suite goldens are the real-WRF wrf_demo case only (synthetic domains removed).
GOLDENS = ROOT / "cases" / "wrf_demo" / "goldens"
INPUTS = GOLDENS / "inputs"


@pytest.fixture(scope="session")
def goldens_root() -> Path:
    return GOLDENS


@pytest.fixture(scope="session")
def inputs_dir() -> Path:
    return INPUTS
