"""Raster typography for the two labels set separately in the v38 figures."""

from io import BytesIO
from pathlib import Path

import matplotlib
from PIL import Image, ImageDraw, ImageFont, PngImagePlugin


def save_labeled_figure(
    figure, output, text, *, font_name, font_size, bounds, background, position,
    anchor, **savefig_options,
):
    """Save vector labels natively and set the PNG label at its print position."""
    if output.suffix.lower() != ".png":
        figure.savefig(output, **savefig_options)
        return

    stream = BytesIO()
    figure.savefig(stream, format="png", **savefig_options)
    stream.seek(0)
    font_path = Path(matplotlib.get_data_path()) / "fonts" / "ttf" / font_name
    font = ImageFont.truetype(str(font_path), font_size)
    with Image.open(stream) as raster:
        pen = ImageDraw.Draw(raster)
        pen.rectangle(bounds, fill=background)
        pen.text(position, text, font=font, anchor=anchor, fill="black")
        metadata = PngImagePlugin.PngInfo()
        for key, value in raster.info.items():
            if isinstance(value, str):
                metadata.add_text(key, value)
        raster.save(output, format="PNG", dpi=raster.info["dpi"], pnginfo=metadata)
