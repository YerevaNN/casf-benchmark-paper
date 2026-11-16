# Main Figure 1 and Supporting Figure S8

Render both figures from the current archived measurements:

```bash
python jcim-overleaf/results-data/diversity-combined/render.py
python jcim-overleaf/results-data/diversity-combined/verify.py
```

Main Figure 1 combines all ten PLINDER methods' cluster counts and MAT-R
with the eleven CASF core/ref rows (including ChEMBL3D-PB) at the ChEMBL-count
target. Panel axes differ; CASF retains its matched stored-ensemble guides.
Supporting Figure S8 shows the seven threshold curves previously in the
PLINDER figure: RDKit raw and six learned methods. No scores, method selection,
cutoffs, budgets, or aggregations changed. Qwen uses FSQ v2 step41054 and
PLINDER torsion uses the current independent raw/minimized runs.

`plinder-points.csv` and `casf-points.csv` record plotted coordinates.
`provenance.json` records source, renderer and PDF/PNG/SVG hashes. The shared
appendix includes S8 after S7, preserving existing supporting figure numbers.
