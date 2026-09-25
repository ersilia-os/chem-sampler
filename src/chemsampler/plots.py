"""Visualization utilities for optimization results."""

import pandas as pd


def plot_optimization_trajectory(
    summary_csv: str, annotator_id: str | None = None, output_path: str | None = None
):
    """
    Plot annotator scores over rounds from a summary.csv file.

    Parameters
    ----------
    summary_csv : str
        Path to summary.csv from a hill_climb run.
    annotator_id : str, optional
        Which annotator score to plot. If None, plots the primary active annotator
        (the one from the final row, or the first if all None).
    output_path : str, optional
        Save plot to this file (e.g., "trajectory.png"). If None, displays the plot.
    """
    import matplotlib.pyplot as plt

    summary = pd.read_csv(summary_csv)

    if summary.empty:
        print("Summary is empty, nothing to plot.")
        return

    if annotator_id is None:
        active_ids = summary["active_annotator_id"].dropna().unique()
        if len(active_ids) > 0:
            annotator_id = active_ids[-1]
        else:
            print("No active_annotator_id found in summary.")
            return

    if annotator_id not in summary.columns:
        print(
            f"Annotator '{annotator_id}' not found in summary columns: {summary.columns.tolist()}"
        )
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    rounds = summary["round"]
    scores = summary[annotator_id]

    ax.plot(rounds, scores, marker="o", linewidth=2, markersize=6, label=annotator_id)
    ax.set_xlabel("Round", fontsize=12)
    ax.set_ylabel(f"{annotator_id} Score", fontsize=12)
    ax.set_title(f"Optimization Trajectory: {annotator_id}", fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"Plot saved to {output_path}")
    else:
        plt.show()
