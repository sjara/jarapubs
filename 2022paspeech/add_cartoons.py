"""
Add the cartoons (in SVG) to the SVG with the plots.

Run as:
python add_cartoons.py FIG_NUMBER
 or
python add_cartoons.py FIG_NUMBER pdf
(where FIG_NUMBER is the figure number, e.g.: 2)
"""

import sys
import os
import svgutils.compose as svgcomp
from jaratoolbox import settings
import studyparams

if len(sys.argv) < 2:
    print('Usages:\n python add_cartoons.py FIG_NUMBER \n python add_cartoons.py FIG_NUMBER pdf')
    sys.exit(1)
if (len(sys.argv) == 3) and (sys.argv[2]=='pdf'):
    make_pdf = True
else:
    make_pdf = False

plots_path = settings.TEMP_OUTPUT_PATH
cartoons_path = settings.TEMP_OUTPUT_PATH
output_path = settings.TEMP_OUTPUT_PATH

fig = int(sys.argv[1])
if fig==1:
    plots_file = 'plots_behavior.svg'
    cartoons_file = ['cartoon_behavior.svg']
    cartoons_scale = [2.5]
    cartoons_pos = [(53, 0)]
elif fig==2:
    plots_file = 'plots_speech_responsiveness.svg'
    cartoons_file = ['headfixed_mouse_ephys.svg',
                     'histology.svg']
    cartoons_scale = [1.2, 0.45]
    cartoons_pos = [(8, 10), (120, 5)]
elif fig==100:
    plots_file = ''
    cartoons_file = []
    cartoons_scale = []
    cartoons_pos = []

ref_svg = svgcomp.SVG(os.path.join(plots_path, plots_file))
cartoons_svg = []
for indf, onefile in enumerate(cartoons_file):
    fullfile = os.path.join(cartoons_path, onefile)
    oneSVG = svgcomp.SVG(fullfile).scale(cartoons_scale[indf]).move(*cartoons_pos[indf])
    cartoons_svg.append(oneSVG)
output_svg = os.path.join(output_path, plots_file.replace('plots', 'figure'))

# For some reason svgutils doesn't read the correct sizes
UNITS_FACTOR = 0.8  # From svgutils px to pt

width = UNITS_FACTOR * ref_svg.width
height = UNITS_FACTOR * ref_svg.height
svgcomp.Figure(
    f'{width}pt',
    f'{height}pt',
    ref_svg,
    *cartoons_svg,
).save(output_svg)

print('Saved merged SVG to:', output_svg)

# -- Convert SVG to PDF --
if make_pdf:
    output_pdf = output_svg.replace('svg', 'pdf')
    os.system(f'inkscape --export-type="pdf" --export-filename={output_pdf} {output_svg}')
    print('Saved PDF to:', output_pdf)


