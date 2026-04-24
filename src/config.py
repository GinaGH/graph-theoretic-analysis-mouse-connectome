"""
Project-wide configuration for the mouse connectome analysis repository.
"""

from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
SRC_DIR = BASE_DIR / "src"
OUTPUT_DIR = BASE_DIR / "outputs"

NORMALIZED_CONNECTION_DENSITY_FILE = DATA_DIR / "normalized_connection_density.csv"

OUTPUT_DIR.mkdir(exist_ok=True)
