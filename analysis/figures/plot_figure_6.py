"""Draw Figure 6 from the frozen pooled Network recall table."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_common import renderer_args

COHORTS = ("GSE128537", "GSE179328")
ROUTES = ("baseline", "floor_clip", "floor_only", "clip_only")
DISPLAY = ["Standard\nscaling", "BrainTrace\nfloor+clip", "Floor only", "Clip only"]
SHORT_LABELS = {
    "Cingulate gyrus": "Cingulate",
    "Frontal (agranular frontal motor areas)": "Motor",
    "Hippocampal formation": "Hippocampal",
    "Lateral Prefrontal Cortex": "Lateral PFC",
    "Occipital/Temporal": "Occipital/temporal",
    "Operculum/Insula": "Operculum/insula",
    "Orbitomedial Prefrontal Cortex (OMPFC)": "Orbitomedial PFC",
    "Parietal, and Parieto-occipital region": "Parietal",
    "Subcortical": "Subcortical",
    "Temporal": "Temporal",
}


def main():
    args = renderer_args(6, __doc__)
    data = pd.read_csv(args.inputs / "NETWORK_CLASS_RECALL.csv")
    # Match the text sizes inherited by the archived Network panels.
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "axes.titlesize": 11, "axes.labelsize": 9.5,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
    })
    classes = sorted(data.truth.unique())
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 7.7), layout="constrained")
    cmap = plt.get_cmap("Blues").copy()
    cmap.set_bad("#e3e3e3")
    for axis, cohort in zip(axes, COHORTS, strict=True):
        table = data[data.cohort.eq(cohort)].pivot(index="truth", columns="route", values="pooled_Top1")
        matrix = table.reindex(index=classes, columns=list(ROUTES)).to_numpy() * 100
        image = axis.imshow(np.ma.masked_invalid(matrix), vmin=0, vmax=100, cmap=cmap, aspect="auto")
        axis.set_yticks(range(10), [SHORT_LABELS[label] for label in classes])
        axis.set_xticks(range(4), DISPLAY, rotation=0, ha="center")
        for i, j in np.ndindex(matrix.shape):
            value = matrix[i, j]
            axis.text(
                j, i, "NA" if np.isnan(value) else f"{value:.1f}", ha="center", va="center",
                fontsize=9, color="white" if np.isfinite(value) and value >= 60 else "#20303b",
            )
        axis.set_title(cohort, loc="left", fontweight="bold")
    fig.colorbar(image, ax=axes, fraction=0.025, pad=0.025, label="Pooled within-class Top1 recall (%)")
    fig.suptitle("Primary Network recall: BrainTrace and preprocessing comparators", fontsize=12, fontweight="bold")
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
