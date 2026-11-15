"""Draw Figure 1 as a vector schematic and rasterize the manuscript PNG."""

import argparse
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


def draw_workflow_pdf(path, regular_font, bold_font):
    pdfmetrics.registerFont(TTFont('Arial', str(regular_font)))
    pdfmetrics.registerFont(TTFont('ArialBold', str(bold_font)))
    width, height = (720, 405)
    c = canvas.Canvas(str(path), pagesize=(width, height), pageCompression=1)
    c.setTitle('BrainTrace reference construction and single-sample inference')
    c.setAuthor('')
    navy = '#14253f'
    blue = '#0752a3'
    purple = '#692396'

    def rect(x, y, w, h, stroke, fill, r=4):
        c.setStrokeColor(HexColor(stroke))
        c.setFillColor(HexColor(fill))
        c.setLineWidth(0.6)
        c.roundRect(x, height - y - h, w, h, r, fill=1, stroke=1)

    def text(x, y, s, size=10, bold=False, color=navy, anchor='left'):
        c.setFont('ArialBold' if bold else 'Arial', size)
        c.setFillColor(HexColor(color))
        getattr(c, {'left': 'drawString', 'center': 'drawCentredString', 'right': 'drawRightString'}[anchor])(x, height - y, s)

    def block(x, y, w, h, lines, color, fill, size=10, bold=False):
        rect(x, y, w, h, color, fill)
        gap = size * 1.22
        start = y + h / 2 - (len(lines) - 1) * gap / 2 + size * 0.34
        for i, line in enumerate(lines):
            text(x + w / 2, start + i * gap, line, size, bold, color, anchor='center')

    def arrow(x1, y1, x2, y2, color):
        c.setStrokeColor(HexColor(color))
        c.setFillColor(HexColor(color))
        c.setLineWidth(0.8)
        c.line(x1, height - y1, x2, height - y2)
        p = c.beginPath()
        if x1 == x2:
            p.moveTo(x2, height - y2)
            p.lineTo(x2 - 2.8, height - y2 + 5)
            p.lineTo(x2 + 2.8, height - y2 + 5)
        else:
            p.moveTo(x2, height - y2)
            p.lineTo(x2 - 5, height - y2 - 2.8)
            p.lineTo(x2 - 5, height - y2 + 2.8)
        p.close()
        c.drawPath(p, fill=1, stroke=0)
    text(12, 21, 'BrainTrace reference construction and single-sample inference', 16, True)
    for y, label, col, fill in [(34, 'Human', blue, '#edf5fc'), (210, 'Macaque', purple, '#f5eff9')]:
        rect(10, y, 700, 166, col, '#ffffff', 5)
        rect(10, y, 64, 166, col, fill, 5)
        text(42, y + 88, label, 12, True, col, anchor='center')
        text(84, y + 16, 'Reference construction', 11, True, col)
        text(223, y + 16, 'Within training splits', 10, color=navy)
        xs = [84, 240, 405, 575]
        ws = [136, 145, 150, 125]
        ty = y + 26
        if label == 'Human':
            trains = [['UKBEC', 'microarrays'], ['Reference preprocessing', 'RAW / RANK'], ['200 genes + scaling', 'Logistic regression'], ['Frozen human', 'reference']]
            queries = [['One compatible', 'CEL'], ['Apply reference', 'preprocessing'], ['Apply panel, scales', 'and classifier'], ['9 regional', 'candidates']]
        else:
            trains = [['Cynomolgus', 'counts'], ['Full-library', 'ln(1 + CPM)'], ['5,000 genes + scaling', 'Logistic regression'], ['Frozen macaque', 'reference']]
            queries = [['One compatible', 'count profile'], ['Normalize and', 'align genes'], ['Apply floor+clip', 'and classifier'], ['10 anatomical', 'candidates']]
        for i, (x, w, lines) in enumerate(zip(xs, ws, trains)):
            block(x, ty, w, 37, lines, col, fill, bold=i == 3)
            if i < 3:
                arrow(x + w + 2, ty + 18.5, xs[i + 1] - 4, ty + 18.5, col)
        if label == 'Macaque':
            text(480, ty + 47, 'floor+clip Network classifier', 10, color=col, anchor='center')
        fy = y + 79
        c.setStrokeColor(HexColor(col))
        c.setLineWidth(0.8)
        c.line(637.5, height - (ty + 37), 637.5, height - fy)
        c.line(637.5, height - fy, 480, height - fy)
        arrow(480, fy, 480, y + 100, col)
        if label == 'Human':
            c.line(480, height - fy, 312.5, height - fy)
            arrow(312.5, fy, 312.5, y + 100, col)
        rect(489, fy - 3, 105, 13, '#ffffff', '#ffffff', 0)
        text(541.5, fy + 7, 'Fixed parameters', 10, True, col, anchor='center')
        text(84, y + 91, 'Single-sample inference', 11, True, col)
        for i, (x, w, lines) in enumerate(zip(xs, ws, queries)):
            block(x, y + 103, w, 36, lines, col, fill, bold=i == 3)
            if i < 3:
                arrow(x + w + 2, y + 121, xs[i + 1] - 4, y + 121, col)
        note = 'GTEx / HPA: explicit RANK adapters' if label == 'Human' else 'Rhesus: orthology mapping and missing-feature adaptation'
        text(88, y + 156, note, 10, color=col, bold=label == 'Human')
    text(width / 2, 395, 'Separate species models     |     Top1: leading candidate     |     Top3: three candidates', 10, anchor='center')
    c.showPage()
    c.save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--font-regular", type=Path, required=True, help="Arial regular font on this host")
    parser.add_argument("--font-bold", type=Path, required=True, help="Arial bold font on this host")
    parser.add_argument("--pdftoppm", default="pdftoppm", help="Poppler pdftoppm executable")
    parser.add_argument("--output", type=Path, required=True, help="Destination PNG path")
    args = parser.parse_args()
    if args.output.suffix.lower() != ".png":
        parser.error("--output must end in .png")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="braintrace_figure_1_") as workspace:
        pdf = Path(workspace) / "workflow.pdf"
        draw_workflow_pdf(pdf, args.font_regular, args.font_bold)
        subprocess.run(
            [args.pdftoppm, "-singlefile", "-r", "600", "-png", str(pdf), str(args.output.with_suffix(""))],
            check=True,
        )


if __name__ == "__main__":
    main()
