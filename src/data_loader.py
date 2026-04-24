"""
Data-loading utilities for Allen Institute mouse connectome analyses.

All scripts that use normalized_connection_density.csv should load it through
this module so the CSV parsing logic stays consistent across the repository.
See notes below concerning the ipsilateral and contralateral matrices. 
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from config import NORMALIZED_CONNECTION_DENSITY_FILE


@dataclass
class FullMatrixData:
    """
    Full normalized connection-density matrix used by the heatmap workflow.
    """

    technical_labels: list[str]
    column_labels: list[str]
    row_labels: list[str]
    matrix: pd.DataFrame
    matrix_np: np.ndarray


@dataclass
class SquareBlockData:
    """
    Square matrix block used by graph-theoretic analyses.
    """

    row_labels: list[str]
    column_labels: list[str]
    matrix: np.ndarray


def load_raw_normalized_connection_density(
    file_path: Path = NORMALIZED_CONNECTION_DENSITY_FILE,
    header=None,
) -> pd.DataFrame:
    """
    Load the raw normalized_connection_density.csv file.
    """
    print("Reading:", file_path)
    return pd.read_csv(file_path, header=header)


def load_full_normalized_matrix(
    file_path: Path = NORMALIZED_CONNECTION_DENSITY_FILE,
) -> FullMatrixData:
    """
    Load the full normalized connection-density matrix for heatmap generation.

    Expected CSV structure:
    - optional blank first row
    - one row of ipsi / contra labels
    - one row of biological target-region column labels
    - remaining rows are numeric matrix data
    - except:
    - first column contains source-region labels
    """
    df = pd.read_csv(file_path)
    df = df.rename(columns={df.columns[0]: "region"})

    start_row = 0
    if pd.isna(df.iloc[0, 0]):
        start_row = 1

    technical_labels = df.iloc[start_row, 1:].astype(str).tolist()
    column_labels = df.iloc[start_row + 1, 1:].astype(str).tolist()

    df_data = df.iloc[start_row + 2:].reset_index(drop=True)
    row_labels = df_data["region"].astype(str).tolist()

    matrix = (
        df_data.drop(columns=["region"])
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
    )

    return FullMatrixData(
        technical_labels=technical_labels,
        column_labels=column_labels,
        row_labels=row_labels,
        matrix=matrix,
        matrix_np=matrix.to_numpy(dtype=float),
    )


def split_ipsi_contra(matrix: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split a full matrix into ipsilateral and contralateral halves.
    The original data file (at the time of this posting)
    seems to have a 461x461 matrix for the ipsilateral side
    but missing a few columns for the contralateral side
    hence the hardcoded 461
    """
    n_cols = matrix.shape[1]
    half = 461  # if the data changes to two 461x461 matrices use: n_cols // 2

    ipsi = matrix.iloc[:, :half]
    contra = matrix.iloc[:, half:]

    return ipsi, contra


def load_ipsilateral_square_block(
    file_path: Path = NORMALIZED_CONNECTION_DENSITY_FILE,
    expected_shape: tuple[int, int] = (461, 461),
) -> SquareBlockData:
    """
    Load the strict ipsilateral square block used by degree analysis.

    Expected raw CSV layout:
    - row 0: ipsi/contra labels 
    - row 1: biological target-region column labels
    - row 2 onward: numeric data rows
    - column 0: row/source-region labels
    - columns 1:462: ipsilateral square block
    """
    df = load_raw_normalized_connection_density(file_path=file_path, header=None)

    column_labels = df.iloc[1, 1:462].astype(str).tolist()
    row_labels = df.iloc[2:, 0].astype(str).tolist()
    raw_matrix = df.iloc[2:, 1:462]

    numeric_matrix = raw_matrix.apply(pd.to_numeric, errors="coerce")

    bad_mask = numeric_matrix.isna() & raw_matrix.notna()
    if bad_mask.any().any():
        bad_positions = np.argwhere(bad_mask.to_numpy())

        examples = []
        for r, c in bad_positions[:10]:
            examples.append(
                f"(csv_row={r + 2}, csv_col={c + 1}, value={raw_matrix.iat[r, c]!r})"
            )

        raise ValueError(
            "Non-numeric data found in ipsilateral matrix. Examples: "
            + ", ".join(examples)
        )

    matrix = numeric_matrix.to_numpy(dtype=float)

    if matrix.shape != expected_shape:
        raise ValueError(f"Expected matrix shape {expected_shape}, got {matrix.shape}")

    if len(row_labels) != expected_shape[0]:
        raise ValueError(
            f"Expected {expected_shape[0]} row labels, got {len(row_labels)}"
        )

    if len(column_labels) != expected_shape[1]:
        raise ValueError(
            f"Expected {expected_shape[1]} column labels, got {len(column_labels)}"
        )

    return SquareBlockData(
        row_labels=row_labels,
        column_labels=column_labels,
        matrix=matrix,
    )
