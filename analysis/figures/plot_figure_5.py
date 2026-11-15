"""Draw Figure 5 from the frozen paired Network results."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

from plot_common import renderer_args

COHORT_COLORS = {"GSE128537": "#286b8b", "GSE179328": "#c88332"}


def main():
    args = renderer_args(5, __doc__)
    data = pd.read_csv(args.inputs / "NETWORK_PAIRED_CONTRASTS.csv")
    # The archived Network panels inherited these text sizes from the three-tier plot.
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "axes.titlesize": 11, "axes.labelsize": 9.5,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
    })
    endpoints = [
        ("primary_singleton", "A  Primary singleton endpoint"),
        ("secondary_allowed_set", "B  Separate allowed-set endpoint"),
    ]
    entries = [(cohort, metric) for cohort in COHORT_COLORS for metric in ("Top1", "Top3")]
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4), layout="constrained")
    for axis, (endpoint, title) in zip(axes, endpoints, strict=True):
        for i, (cohort, metric) in enumerate(entries):
            row = data[
                data.cohort.eq(cohort) & data.endpoint.eq(endpoint)
                & data.route.eq("floor_clip") & data.metric.eq(metric)
            ].iloc[0]
            value, low, high = 100*row.paired_difference, 100*row.CI95_low, 100*row.CI95_high
            axis.errorbar(value, 3-i, xerr=[[value-low], [high-value]], fmt="o", capsize=4,
                          color=COHORT_COLORS[cohort])
            axis.text(value, 3-i+0.18, f"{value:+.2f} [{low:+.2f}, {high:+.2f}]", ha="center", fontsize=9)
        axis.axvline(0, color="#666", ls="--", lw=1)
        axis.set_ylim(-0.5, 3.7)
        axis.set_yticks([3, 2, 1, 0], [f"{cohort} · {metric}" for cohort, metric in entries])
        axis.set_title(title, loc="left", fontweight="bold")
        axis.set_xlabel("BrainTrace minus standard scaling (percentage points)")
        axis.grid(axis="x", alpha=0.15)
    axes[0].set_xlim(-7, 11)
    axes[1].set_xlim(-30, 18)
    fig.suptitle(
        "External Network comparison: paired donor differences and descriptive 95% intervals",
        fontsize=12, fontweight="bold",
    )
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
