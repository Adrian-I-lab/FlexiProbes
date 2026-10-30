#!/usr/bin/env python3
"""Generate the flexiprobes app icon in every format the build needs.

Outputs (under packaging/assets/ unless --out is given):
  • icon.png          1024px master (also copied to the GUI package for the
                      runtime window icon)
  • icon.ico          Windows (PNG-embedded multi-size ICO)
  • icon.icns         macOS (built via iconutil when available)

The artwork is drawn with Qt (no Pillow needed). Re-run after changing the
design; commit the regenerated assets.

Usage:
    python packaging/make_icons.py            # needs PySide6 installed
"""
from __future__ import annotations

import os
import struct
import subprocess
import sys
from pathlib import Path

# Headless rendering — no display required.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QBuffer, QByteArray, QIODevice, QPointF, Qt  # noqa: E402
from PySide6.QtGui import (  # noqa: E402
    QBrush, QColor, QImage, QPainter, QPainterPath, QPen,
)
from PySide6.QtWidgets import QApplication  # noqa: E402

MASTER = 1024
ICO_SIZES = (16, 32, 48, 64, 128, 256)
ICNS_SIZES = (16, 32, 64, 128, 256, 512, 1024)


# Logo geometry on a 100-unit canvas, y pointing down. Two bent arms are the
# halves of a split probe; the gap between them is the ligation site.
NAVY, WHITE, CYAN = "#0B3C5D", "#FFFFFF", "#00AEEF"
ARM_WIDTH = 12
LEFT_ARM = ((8.4, 37.68), (24.0, 53.28), (38.8, 53.28))
RIGHT_ARM = ((61.2, 53.28), (76.0, 53.28), (91.6, 37.68))


def render_master(size: int = MASTER) -> QImage:
    """Draw the flexiprobes logo: a white and a cyan probe half on a navy rounded square."""
    img = QImage(size, size, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(Qt.GlobalColor.transparent)
    p = QPainter(img)
    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    unit = size / 100.0

    radius = size * 0.22
    bg = QPainterPath()
    bg.addRoundedRect(0, 0, size, size, radius, radius)
    p.fillPath(bg, QBrush(QColor(NAVY)))

    for points, colour in ((LEFT_ARM, WHITE), (RIGHT_ARM, CYAN)):
        pen = QPen(QColor(colour))
        pen.setWidthF(ARM_WIDTH * unit)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        p.setPen(pen)
        path = QPainterPath()
        path.moveTo(QPointF(points[0][0] * unit, points[0][1] * unit))
        for x, y in points[1:]:
            path.lineTo(QPointF(x * unit, y * unit))
        p.drawPath(path)

    p.end()
    return img


def _png_bytes(img: QImage, size: int) -> bytes:
    scaled = img.scaled(
        size, size,
        Qt.AspectRatioMode.KeepAspectRatio,
        Qt.TransformationMode.SmoothTransformation,
    )
    ba = QByteArray()
    buf = QBuffer(ba)
    buf.open(QIODevice.OpenModeFlag.WriteOnly)
    scaled.save(buf, "PNG")
    buf.close()
    return bytes(ba)


def write_ico(master: QImage, path: Path) -> None:
    """Assemble a PNG-embedded .ico (Vista+ format) by hand — no Pillow."""
    entries, blobs = [], []
    offset = 6 + 16 * len(ICO_SIZES)
    for size in ICO_SIZES:
        png = _png_bytes(master, size)
        w = h = 0 if size >= 256 else size  # 0 means 256 in ICO
        entries.append(struct.pack("<BBBBHHII", w, h, 0, 0, 1, 32, len(png), offset))
        blobs.append(png)
        offset += len(png)
    with open(path, "wb") as f:
        f.write(struct.pack("<HHH", 0, 1, len(ICO_SIZES)))  # ICONDIR
        for e in entries:
            f.write(e)
        for b in blobs:
            f.write(b)


def write_icns(master: QImage, path: Path) -> bool:
    """Build a .icns via macOS iconutil. Returns False if unavailable."""
    import shutil
    iconutil = shutil.which("iconutil")
    if not iconutil:
        return False
    iconset = path.with_suffix(".iconset")
    iconset.mkdir(parents=True, exist_ok=True)
    name = {16: "16x16", 32: "32x32", 128: "128x128", 256: "256x256", 512: "512x512"}
    for size in ICNS_SIZES:
        png = _png_bytes(master, size)
        if size in name:
            (iconset / f"icon_{name[size]}.png").write_bytes(png)
        if size // 2 in name:                      # @2x retina variant
            (iconset / f"icon_{name[size // 2]}@2x.png").write_bytes(png)
    subprocess.run([iconutil, "-c", "icns", str(iconset), "-o", str(path)],
                   check=True)
    for child in iconset.iterdir():
        child.unlink()
    iconset.rmdir()
    return True


def main() -> None:
    here = Path(__file__).resolve().parent
    out = here / "assets"
    out.mkdir(parents=True, exist_ok=True)

    # A QApplication must exist before any QImage/QPainter work; we only need
    # the side effect of constructing it, not a reference.
    QApplication.instance() or QApplication(sys.argv)
    master = render_master()

    png_path = out / "icon.png"
    master.save(str(png_path), "PNG")
    print(f"OK  {png_path}")

    ico_path = out / "icon.ico"
    write_ico(master, ico_path)
    print(f"OK  {ico_path}")

    icns_path = out / "icon.icns"
    if write_icns(master, icns_path):
        print(f"OK  {icns_path}")
    else:
        print("SKIP icon.icns (iconutil not found — build on macOS to get it)")

    # Runtime window icon: ship the PNG inside the GUI package.
    gui_assets = here.parent / "src" / "flexprobe_gui" / "assets"
    gui_assets.mkdir(parents=True, exist_ok=True)
    master.save(str(gui_assets / "icon.png"), "PNG")
    print(f"OK  {gui_assets / 'icon.png'}")


if __name__ == "__main__":
    main()
