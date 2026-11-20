"""Draw Figure 4 from saved Top1 intersections for GTEx HIP and HPA frontal profiles."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_common import renderer_args


def main():
    args = renderer_args(4, __doc__)
    data = pd.read_csv(args.inputs / "HUMAN_TOP1_INTERSECTIONS.csv")
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42,
        "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
    })
    categories = [
        ("both_hit", "Both hit", "#547f98"),
        ("BrainTrace_only_hit", "BrainTrace only", "#7bb5b0"),
        ("comparator_only_hit", "Comparator only", "#ddac58"),
        ("both_miss", "Both miss", "#b86e70"),
    ]
    panels = [
        ("GTEx", "GTEx hippocampus: 255 samples / donors"),
        ("HPA", "HPA frontal subregions: 55 samples"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.9), sharey=True)
    for axis, (dataset, title) in zip(axes, panels, strict=True):
        rows = data[data.dataset == dataset].set_index("comparator").loc[["Pearson_centroid", "RF500"]]
        left = np.zeros(2)
        for field, label, color in categories:
            counts = rows[field].to_numpy()
            percentages = 100 * counts / rows.samples.to_numpy()
            axis.barh([0, 1], percentages, left=left, color=color, label=label, height=0.55)
            for i, (fraction, count) in enumerate(zip(percentages, counts, strict=True)):
                if fraction >= 3:
                    axis.text(
                        left[i]+fraction/2, i, str(count), ha="center", va="center", fontsize=10,
                        color="white" if field in ("both_hit", "both_miss") else "#17212b",
                    )
            left += percentages
        axis.set_yticks([0, 1], ["Pearson", "RF500"])
        axis.set_xlim(0, 100)
        axis.set_xlabel("Share of matched samples (%)")
        axis.set_title(title, fontsize=11, pad=14, loc="left")
        axis.invert_yaxis()
    axes[1].tick_params(axis="y", labelleft=True)
    fig.suptitle("Shared and method-specific Top1 errors", x=0.04, ha="left", fontweight="bold", fontsize=14)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.04))
    fig.text(
        0.04, 0.015, "Numbers inside bars are sample counts. HPA has no donor links; no donor-level inference is made.",
        fontsize=9, color="#424c59",
    )
    fig.subplots_adjust(top=0.78, bottom=0.28, wspace=0.16, left=0.10, right=0.98)
    fig.savefig(args.output, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
