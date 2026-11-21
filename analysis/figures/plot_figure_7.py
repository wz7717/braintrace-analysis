"""Draw Figure 7 from saved three-tier results and their original evaluation sets."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

from plot_common import read_rows, renderer_args
from plot_raster_label import save_labeled_figure


def one(rows, **keys):
    row, = [row for row in rows if all(row[key] == value for key, value in keys.items())]
    return row


def main():
    args = renderer_args(7, __doc__)
    source = read_rows(args.inputs / "SOURCE_THREE_TIER.csv")
    external = read_rows(args.inputs / "EXTERNAL_THREE_TIER.csv")
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9.5, 'axes.titlesize': 11, 'axes.labelsize': 9.5, 'xtick.labelsize': 9, 'ytick.labelsize': 9, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
    blue = '#226E92'
    gray = '#777777'
    fig = plt.figure(figsize=(10.2, 5.15), layout='constrained')
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.3], width_ratios=[1, 1.05])
    ax = fig.add_subplot(gs[0, :])
    ax.set_axis_off()
    ax.set_xlim(-0.006, 1.01)
    ax.set_ylim(0, 1)
    ax.text(0, 1.0, 'A  BrainTrace macaque workflow', ha='left', va='top', weight='bold', fontsize=11)

    def box(x, y, w, h, txt, fill):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.006,rounding_size=0.008', linewidth=0.8, edgecolor='#6A7780', facecolor=fill))
        ax.text(x + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=9.3)

    def arrow(a, b):
        ax.annotate('', xy=b, xytext=a, arrowprops={'arrowstyle': '->', 'lw': 1.1, 'color': '#444444'})
    input_label = "Align counts to source genes"
    box(0.005, 0.28, 0.2, 0.42, f'Individual profile\nFull-library ln(1 + CPM)\n{input_label}', '#F0F2F4')
    box(0.25, 0.5, 0.27, 0.34, 'Network\n5,000 genes; floor + clip\nFrozen logistic classifier', '#E8F1F6')
    box(0.565, 0.5, 0.17, 0.34, 'Top3 Network\ncandidate regions', '#E8F1F6')
    box(0.785, 0.56, 0.21, 0.28, 'Group\nPearson over 200 genes\nMaximum member score', '#EAF3EF')
    box(0.785, 0.12, 0.21, 0.28, 'Exact\n0.25 z(r50) + 0.75 z(r100)\nFine-region ranking', '#EAF3EF')
    arrow((0.205, 0.59), (0.25, 0.67))
    arrow((0.52, 0.67), (0.565, 0.67))
    arrow((0.735, 0.67), (0.785, 0.7))
    ax.plot([0.755, 0.755], [0.26, 0.67], c='#444444', lw=1.1)
    arrow((0.755, 0.26), (0.785, 0.26))
    ax.plot([0.105, 0.105, 0.715], [0.28, 0.18, 0.18], c='#444444', lw=1.0)
    arrow((0.715, 0.18), (0.715, 0.48))
    ax.text(0.4, 0.23, 'Raw lower-tier inputs and source references', ha='center', fontsize=9.1)
    ax.text(0.4, 0.08, 'No Network score multiplier; Group does not gate Exact', ha='center', fontsize=9.1)
    ax = fig.add_subplot(gs[1, 0])
    tiers = ['Network', 'Group', 'Exact']
    y = np.arange(3)
    for arm, color, offset in [('original_three_tier', gray, -0.12), ('floor_clip_hybrid', blue, 0.12)]:
        for i, tier in enumerate(tiers):
            r = one(source, arm=arm, tier=tier)
            vals = []
            for metric, marker in [('Top1', 'o'), ('Top3', 's')]:
                v = float(r[metric]) * 100
                vals.append(v)
                lo = float(r[metric + '_CI95_low']) * 100
                hi = float(r[metric + '_CI95_high']) * 100
                ax.errorbar(v, i + offset, xerr=[[v - lo], [hi - v]], fmt=marker, color=color, markersize=4.5, capsize=2, elinewidth=0.8)
            ax.plot(vals, [i + offset] * 2, color=color, alpha=0.25, lw=1)
    ax.set_yticks(y, tiers)
    ax.set_ylim(2.55, -0.65)
    ax.set_xlim(0, 103)
    ax.set_xlabel('Donor-averaged accuracy (%)')
    ax.grid(axis='x', alpha=0.16)
    ax.set_title('B  Source holdouts\n819 specimens / 9 donors / 110 fine labels', loc='left')
    ax.legend(handles=[Line2D([], [], color=gray, label='VSD-Pearson'), Line2D([], [], color=blue, label='BrainTrace'), Line2D([], [], color='black', marker='o', ls='', label='Top1'), Line2D([], [], color='black', marker='s', ls='', label='Top3')], loc='upper left', bbox_to_anchor=(0, -0.27), ncol=4, frameon=False, fontsize=8.5, handlelength=1.3, columnspacing=0.9)
    ax = fig.add_subplot(gs[1, 1])
    labels = []
    for j, cohort in enumerate(['GSE128537', 'GSE179328']):
        for k, tier in enumerate(tiers):
            r = one(external, cohort=cohort, tier=tier, evaluation_set='primary_singleton' if tier == 'Network' else 'primary_unique')
            pos = j * 3 + k
            labels.append(f"{(cohort if k == 0 else '')}\n{tier}: {r['samples']} / {r['donors']}" if k == 0 else f"{tier}: {r['samples']} / {r['donors']}")
            for metric, marker in [('Top1', 'o'), ('Top3', 's')]:
                v = float(r[metric]) * 100
                lo = float(r[metric + '_CI95_low']) * 100
                hi = float(r[metric + '_CI95_high']) * 100
                ax.errorbar(v, pos, xerr=[[v - lo], [hi - v]], fmt=marker, color=blue, capsize=2, markersize=4.5, elinewidth=0.8)
    ax.set_yticks(range(6), labels, fontsize=8.5)
    ax.set_ylim(5.6, -0.75)
    ax.set_xlim(0, 103)
    ax.axhline(2.5, color='#CCCCCC', lw=0.7)
    ax.grid(axis='x', alpha=0.16)
    ax.set_xlabel('Donor-averaged accuracy (%)')
    ax.set_title('C  BrainTrace external results\nDistinct Network and fine evaluation sets', loc='left')
    ax.text(0, -0.3, 'Labels show specimens / donors; fine sets contain 5 labels each.', transform=ax.transAxes, fontsize=8.4)
    save_labeled_figure(
        fig, args.output, input_label, font_name="DejaVuSans.ttf", font_size=39,
        bounds=(221, 294, 798, 336), background="#f0f2f4", position=(510, 296),
        anchor="mt", dpi=300,
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
