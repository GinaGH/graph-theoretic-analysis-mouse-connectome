"""
Ipsilateral degree analysis for the Allen Institute mouse connectome data.

This script:
- loads the strict 461x461 ipsilateral block
- thresholds the weighted matrix into an adjacency matrix
- computes in-degree and out-degree
- saves a CSV and several visualizations
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


try:
    from .config import OUTPUT_DIR
    from .data_loader import load_ipsilateral_square_block
except ImportError:
    from config import OUTPUT_DIR
    from data_loader import load_ipsilateral_square_block


THRESHOLD = 0.0001  # threshold to remove noise from the data


def compute_degrees(
    matrix: np.ndarray,
    threshold: float = THRESHOLD,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert a weighted matrix to a binary adjacency matrix based on the selected threshold and compute degrees.

    A_ij is interpreted as an edge i -> j.
    """
    adjacency = (matrix > threshold).astype(int)

    out_degree = adjacency.sum(axis=1)
    in_degree = adjacency.sum(axis=0)

    return out_degree, in_degree


def save_degree_csv(
    labels: list[str],
    out_degree: np.ndarray,
    in_degree: np.ndarray,
    filename: str = "ipsi_degree.csv",
) -> Path:
    """
    Save degree results to outputs/ipsi_degree.csv.
    """
    degree_df = pd.DataFrame(
        {
            "region": labels,
            "out_degree": out_degree,
            "in_degree": in_degree,
        }
    )

    csv_path = OUTPUT_DIR / filename
    degree_df.to_csv(csv_path, index=False)

    return csv_path


def save_histogram(
    values: np.ndarray,
    title: str,
    xlabel: str,
    filename: str,
) -> None:
    """
    Save a histogram plot.
    """
    plt.figure(figsize=(8, 5))
    plt.hist(values, bins=30)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / filename, dpi=220, bbox_inches="tight")
    plt.show()


def save_top_bar(
    labels: list[str],
    values: np.ndarray,
    title: str,
    filename: str,
    top_n: int = 20,
) -> None:
    """
    Save a bar chart of the top regions by degree.
    """
    order = np.argsort(values)[::-1][:top_n]
    top_labels = [labels[i] for i in order]
    top_values = values[order]

    plt.figure(figsize=(12, 6))
    plt.bar(range(len(top_values)), top_values)
    plt.xticks(range(len(top_values)), top_labels, rotation=90, fontsize=8)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / filename, dpi=220, bbox_inches="tight")
    plt.show()


def save_in_vs_out_scatter(
    out_degree: np.ndarray,
    in_degree: np.ndarray,
) -> None:
    """
    Save scatter plot comparing out-degree and in-degree.
    """
    plt.figure(figsize=(6, 6))
    plt.scatter(out_degree, in_degree)
    plt.xlabel("Out-degree")
    plt.ylabel("In-degree")
    plt.title("Out vs In Degree (Ipsilateral)")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "in_vs_out_scatter.png", dpi=220, bbox_inches="tight")
    plt.show()


def main() -> None:
    """
    Run the full ipsilateral degree-analysis workflow.
    """
    OUTPUT_DIR.mkdir(exist_ok=True)

    data = load_ipsilateral_square_block()

    row_labels = data.row_labels
    column_labels = data.column_labels
    matrix = data.matrix

    if set(row_labels) != set(column_labels):
        print("Mismatched Labels! Row and column label sets differ!")

    print("Matrix shape:", matrix.shape)

    out_degree, in_degree = compute_degrees(matrix, threshold=THRESHOLD)

    csv_path = save_degree_csv(row_labels, out_degree, in_degree)
    print("Saved:", csv_path)

    save_histogram(
        out_degree,
        "Ipsilateral Out-Degree Distribution",
        "Out-degree",
        "ipsi_out_degree_hist.png",
    )

    save_histogram(
        in_degree,
        "Ipsilateral In-Degree Distribution",
        "In-degree",
        "ipsi_in_degree_hist.png",
    )

    save_top_bar(
        row_labels,
        out_degree,
        "Top 20 Regions by Out-Degree",
        "top_out_degree.png",
    )

    save_top_bar(
        row_labels,
        in_degree,
        "Top 20 Regions by In-Degree",
        "top_in_degree.png",
    )

    save_in_vs_out_scatter(out_degree, in_degree)

    print("Saved plots to:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
