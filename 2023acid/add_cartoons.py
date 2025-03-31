"""
Add the cartoons (in SVG) to the SVG with the plots.

Run as:
python add_cartoons.py FIG_NUMBER
 or
python add_cartoons.py FIG_NUMBER pdf
(where FIG_NUMBER is the figure number, e.g.: 2)

This script requires installing 'svgutils' 
See ~/src/svg_utils for a version without requirements to install with:
pip3 install -e ./

NOTES:
- It looks like the cartoon file cannot be in units of 'inches'. Gives error "KeyError: 'in'
- The scaling is messed up.
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
cartoon_path = '/home/sjara/papers/2023acid/figures/' #settings.TEMP_OUTPUT_PATH
output_path = settings.TEMP_OUTPUT_PATH

fig = int(sys.argv[1])
if fig==1:
    plots_file = 'plots_overall_firing.svg'
    cartoon_files = ['cartoon_injection.svg']
    cartoon_scale = [2.8]
    cartoon_pos = [(16, 28)]
elif fig==2:
    plots_file = 'plots_freqtuning.svg'
    cartoon_files = ['cartoon_tuningcurve.svg']
    cartoon_scale = [2.8]
    cartoon_pos = [(658, 10)]
elif fig==3:
    plots_file = 'plots_oddball_enhancement.svg'
    cartoon_files = [ 'cartoon_oddball_chord.svg', 'cartoon_oddball_FM.svg']
    cartoon_scale = [3.2, 3.2]
    cartoon_pos = [(16, 24), (16, 218)]
elif fig==100:
    plots_file = ''
    cartoon_files = []
    cartoon_scale = []
    cartoon_pos = []

ref_svg = svgcomp.SVG(os.path.join(plots_path, plots_file))
cartoon_svg = []
for indf, onefile in enumerate(cartoon_files):
    fullfile = os.path.join(cartoon_path, onefile)
    oneSVG = svgcomp.SVG(fullfile).scale(cartoon_scale[indf]).move(*cartoon_pos[indf])
    cartoon_svg.append(oneSVG)
output_svg = os.path.join(output_path, plots_file.replace('plots', 'figure'))

# For some reason svgutils doesn't read the correct sizes
UNITS_FACTOR = 0.8  # From svgutils px to pt

width = UNITS_FACTOR * ref_svg.width
height = UNITS_FACTOR * ref_svg.height
svgcomp.Figure(
    f'{width}pt',
    f'{height}pt',
    ref_svg,
    *cartoon_svg,
).save(output_svg)

print('Saved merged SVG to:', output_svg)

# -- Convert SVG to PDF --
if make_pdf:
    output_pdf = output_svg.replace('svg', 'pdf')
    os.system(f'inkscape --export-type="pdf" --export-filename={output_pdf} {output_svg}')
    print('Saved PDF to:', output_pdf)


