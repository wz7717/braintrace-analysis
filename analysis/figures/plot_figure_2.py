"""Draw Figure 2 from the frozen RAW-fixed confusion counts and recall intervals."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from plot_common import read_rows, renderer_args


def main():
    args = renderer_args(2, __doc__)
    metrics = read_rows(args.inputs / "RAW_CLASS_RECALL.csv")
    counts = read_rows(args.inputs / "RAW_CONFUSION_COUNTS.csv")
    classes = [row["class_id"] for row in metrics]
    labels = [row["label"] for row in metrics]
    predictions = classes + ["FAIL"]
    matrix = np.zeros((len(classes), len(predictions)), dtype=int)
    for row in counts:
        matrix[classes.index(row["truth"]), predictions.index(row["prediction"])] = int(row["n"])

    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, (confusion, recall) = plt.subplots(1, 2, figsize=(16.5, 6.2), gridspec_kw={"width_ratios": [1.15, 1]})
    confusion.imshow(matrix, cmap="Blues", vmin=0, vmax=128, aspect="auto")
    for truth, prediction in np.ndindex(matrix.shape):
        count = matrix[truth, prediction]
        confusion.text(
            prediction, truth, str(count), ha="center", va="center", fontsize=9.4,
            color="white" if count > 75 else "#243746",
        )
    confusion.set_xticks(range(10), labels + ["FAIL"], rotation=55, ha="right", fontsize=9)
    confusion.set_yticks(range(9), labels, fontsize=10)
    confusion.set_xlabel("Predicted region")
    confusion.set_ylabel("True region")
    confusion.set_title("A  Human RAW out-of-fold predictions", loc="left", fontsize=13, weight="bold", pad=13)

    for i, row in enumerate(metrics):
        value = float(row["donor_balanced_recall"]) * 100
        low = float(row["lower_95"]) * 100
        high = float(row["upper_95"]) * 100
        recall.plot([low, high], [i, i], color="#25788a", lw=2)
        recall.plot(value, i, "o", color="#25788a", ms=6)
        label = f'{value:.1f}%   {row["truth_n"]}/{row["truth_donors_n"]}'
        recall.text(101, i, label, va="center", fontsize=9.5)
    recall.set_yticks(range(9), labels, fontsize=10)
    recall.set_ylim(8.6, -0.6)
    recall.set_xlim(70, 125)
    recall.set_xticks([70, 80, 90, 100])
    recall.grid(axis="x", alpha=0.2)
    recall.set_axisbelow(True)
    recall.set_xlabel("Donor-balanced regional recall (%)")
    recall.set_title("B  Regional recall and 95% donor intervals", loc="left", fontsize=13, weight="bold", pad=13)
    recall.text(101, -0.42, "Recall   n/donors", fontsize=9.2)
    recall.spines["bottom"].set_bounds(70, 100)
    for side in ("top", "right"):
        recall.spines[side].set_visible(False)
    fig.subplots_adjust(left=0.115, right=0.98, bottom=0.25, top=0.88, wspace=0.5)
    fig.savefig(args.output, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
