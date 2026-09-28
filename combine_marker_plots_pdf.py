"""
Combines every PNG in full_recording_marker_plots\ into a single
multi-page PDF (one plot per page, landscape, to match the wide aspect
ratio of the plots), so you can scroll through all subjects in one file
instead of opening 50 separate images.
"""
import glob
import os

from fpdf import FPDF
from PIL import Image

from config import OUT_DIR
PLOTS_DIR = os.path.join(OUT_DIR, 'full_recording_marker_plots')
OUT_PATH = os.path.join(PLOTS_DIR, 'all_subjects_markers.pdf')

pngs = sorted(glob.glob(os.path.join(PLOTS_DIR, '*.png')))
print(f'Found {len(pngs)} plots to combine.')

pdf = FPDF(orientation='L', format='A4')  # landscape - matches the wide plot shape
pdf.set_auto_page_break(auto=False)

page_w, page_h = 297, 210  # A4 landscape, mm
margin = 8

for png_path in pngs:
    with Image.open(png_path) as img:
        img_w, img_h = img.size
    aspect = img_h / img_w

    avail_w = page_w - 2 * margin
    avail_h = page_h - 2 * margin
    draw_w = avail_w
    draw_h = draw_w * aspect
    if draw_h > avail_h:
        draw_h = avail_h
        draw_w = draw_h / aspect

    x = (page_w - draw_w) / 2
    y = (page_h - draw_h) / 2

    pdf.add_page()
    pdf.image(png_path, x=x, y=y, w=draw_w, h=draw_h)

pdf.output(OUT_PATH)
print(f'Saved combined PDF: {OUT_PATH}')

# remove the individual PNGs now that they're combined - keep only the PDF
removed = 0
for png_path in pngs:
    os.remove(png_path)
    removed += 1
print(f'Removed {removed} individual PNG files - only the combined PDF remains.')
