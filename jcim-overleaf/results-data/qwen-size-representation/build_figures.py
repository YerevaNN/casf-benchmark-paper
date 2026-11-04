"""Reproduce the two approved Qwen comparisons from archived measurements."""
from pathlib import Path
import shutil
import subprocess
import sys

here = Path(__file__).resolve().parent
subprocess.run([sys.executable, str(here / 'plot_main_two_row.py')], check=True)
subprocess.run([sys.executable, str(here / 'plot_two_row.py'), '--simplified'], check=True)
artwork = here.parents[1] / 'figures'
for source, destination in [
    ('08-main-size-and-representation', 'figure-qwen-size-representation'),
    ('07-representation-and-size-simplified', 'figure-qwen-molecule-effects'),
]:
    for extension in ['pdf', 'svg', 'png']:
        shutil.copy2(here / f'{source}.{extension}', artwork / f'{destination}.{extension}')
