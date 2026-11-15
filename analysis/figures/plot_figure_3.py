"""Draw Figure 3 from frozen GTEx donor contrasts and regional Top1 recall."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_common import renderer_args

METHODS = ("BrainTrace_RANK", "Pearson_centroid", "RF500")
LABELS = {"BrainTrace_RANK": "BrainTrace RANK", "Pearson_centroid": "Pearson centroid", "RF500": "RF500"}
COLORS = {"BrainTrace_RANK": "#28648a", "Pearson_centroid": "#d28a32", "RF500": "#7f5a9f"}


def recall_heatmap(axis, data):
    regions = data.region.drop_duplicates().tolist()
    matrix = np.array([
        [float(data[(data.region == region) & (data.method == method)].Top1.iloc[0]) * 100
         for method in METHODS]
        for region in regions
    ])
    image = axis.imshow(matrix, cmap="Blues", vmin=0, vmax=100, aspect="auto")
    axis.set_xticks(range(3), [LABELS[method] for method in METHODS], rotation=18, ha="right")
    axis.set_yticks(range(len(regions)), [
        region.replace("Brain - ", "").replace(" (basal ganglia)", "").replace(
            "Frontal Cortex (BA9)", "Frontal cortex (BA9)"
        ) for region in regions
    ])
    for i, j in np.ndindex(matrix.shape):
        value = matrix[i, j]
        axis.text(j, i, f"{value:.1f}", ha="center", va="center",
                  color="white" if value >= 60 else "#142534")
    axis.set_title("B  Regional Top1 recall (%)", loc="left", fontweight="bold", pad=14)
    return image


def main():
    args = renderer_args(3, __doc__)
    contrasts = pd.read_csv(args.inputs / "GTEX_PAIRED_CONTRASTS.csv")
    regions = pd.read_csv(args.inputs / "GTEX_REGIONAL_RECALL.csv")
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })

    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.4), gridspec_kw={"width_ratios": [1.12, 1]}, layout="constrained")
    compare = axes[0]
    order = [(method, metric) for metric in ("Top1", "Top3") for method in METHODS[1:]]
    for y, (method, metric) in enumerate(order):
        row = contrasts[(contrasts.comparator == method) & (contrasts.metric == metric)].iloc[0]
        value = row.paired_difference * 100
        low = row.CI95_low * 100
        high = row.CI95_high * 100
        compare.errorbar(value, 3-y, xerr=[[value-low], [high-value]], fmt="o",
                         color=COLORS[method], capsize=4, markersize=7)
        compare.text(value, 3-y+0.18, f"{value:+.2f} [{low:+.2f}, {high:+.2f}]", ha="center", fontsize=9)
    compare.axvline(0, color="#757575", lw=1, ls="--")
    compare.set_yticks([3, 2, 1, 0], [f"{LABELS[method]} · {metric}" for method, metric in order])
    compare.set_xlim(-32, 11)
    compare.set_ylim(-0.5, 3.65)
    compare.set_xlabel("Comparator − BrainTrace (percentage points)")
    compare.set_title("A  GTEx primary: paired donor contrasts", loc="left", fontweight="bold", pad=14)
    compare.grid(axis="x", alpha=0.15)

    image = recall_heatmap(axes[1], regions)
    fig.colorbar(image, ax=axes[1], fraction=0.04, pad=0.04, label="Top1 (%)")
    fig.suptitle("GTEx: 961 samples from 345 donors; fixed RANK adaptation", fontsize=13, fontweight="bold")
    fig.savefig(args.output, dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
