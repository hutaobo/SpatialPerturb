from pathlib import Path
import sys

import pytest


SRC_DIR = Path(__file__).resolve().parents[1] / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import spatialperturb as sp


@pytest.fixture()
def demo_adata():
    return sp.load_demo_dataset()


@pytest.fixture()
def demo_pair():
    return sp.load_demo_dataset(paired=True)


@pytest.fixture()
def demo_unannotated():
    return sp.load_demo_dataset(annotate=False)
