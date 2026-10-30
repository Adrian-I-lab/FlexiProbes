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

import math
import os
import struct
import subprocess
import sys
from pathlib import Path

# Headless rendering — no display required.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QBuffer, QByteArray, QIODevice, QPointF, Qt  # noqa: E402
from PySide6.QtGui import (  # noqa: E402
    QBrush, QColor, QImage, QLinearGradient, QPainter, QPainterPath, QPen,
)
from PySide6.QtWidgets import QApplication  # noqa: E402

MASTER = 1024
ICO_SIZES = (16, 32, 48, 64, 128, 256)
ICNS_SIZES = (16, 32, 64, 128, 256, 512, 1024)


def render_master(size: int = MASTER) -> QImage:
    """Draw a DNA double-helix on a gradient rounded square."""
    img = QImage(size, size, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(Qt.GlobalColor.transparent)
    p = QPainter(img)
    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)

    # Rounded-square background with a teal→indigo gradient.
    grad = QLinearGradient(0, 0, size, size)
    grad.setColorAt(0.0, QColor("#14b8a6"))
    grad.setColorAt(1.0, QColor("#4f46e5"))
    radius = size * 0.22
    bg = QPainterPath()
    bg.addRoundedRect(0, 0, size, size, radius, radius)
    p.fillPath(bg, QBrush(grad))

    # Double helix.
    cx = size / 2.0
    amp = size * 0.17
    top, bottom = size * 0.2, size * 0.8
    span = bottom - top
    turns = 2.5
    n = 240

    def strand(phase: float) -> QPainterPath:
        path = QPainterPath()
        for i in range(n + 1):
            frac = i / n
            y = top + frac * span
            ang = frac * turns * 2 * math.pi + phase
            x = cx + amp * math.sin(ang)
            if i == 0:
                path.moveTo(QPointF(x, y))
            else:
                path.lineTo(QPointF(x, y))
        return path

    white = QColor(255, 255, 255)
    rung_pen = QPen(QColor(255, 255, 255, 150))
    rung_pen.setWidthF(size * 0.022)
    rung_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    p.setPen(rung_pen)
    rungs = 11
    for j in range(1, rungs):
        frac = j / rungs
        y = top + frac * span
        ang = frac * turns * 2 * math.pi
        xa = cx + amp * math.sin(ang)
        xb = cx - amp * math.sin(ang)
        p.drawLine(QPointF(xa, y), QPointF(xb, y))

    strand_pen = QPen(white)
    strand_pen.setWidthF(size * 0.055)
    strand_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    strand_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    p.setPen(strand_pen)
    p.drawPath(strand(0.0))
    p.drawPath(strand(math.pi))

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
