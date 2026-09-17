from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
GOLDENS = ROOT / "cases" / "small_domain" / "goldens"
INPUTS = GOLDENS / "inputs"


@pytest.fixture(scope="session")
def goldens_root() -> Path:
    return GOLDENS


@pytest.fixture(scope="session")
def inputs_dir() -> Path:
    return INPUTS
