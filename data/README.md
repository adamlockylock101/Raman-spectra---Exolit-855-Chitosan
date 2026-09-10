Real spectra go here as `.txt` files: two columns, wavenumber then intensity.

`python3 analyse_raman.py` with no arguments analyses every `.txt` in this
directory. The loader auto-detects the delimiter (tab, comma, semicolon or
whitespace), skips comment and instrument-preamble lines, ignores columns beyond
the first two, and handles descending wavenumber order.

Filenames containing `treated` or `coated` are auto-labelled as the treated
sample; those containing `control`, `ctrl` or `reference` as the control. That
labelling is what produces the treated-vs-control comparison table. Override it
with `--label`, given once per file in argument order:

```bash
python3 analyse_raman.py a.txt b.txt \
    --label "Ch/Ex" \
    --label "Control"
```

This directory is empty in the repository. Runnable example spectra — with
known ground truth — live in `examples/synthetic/`, and the README's
thirty-second demo uses those.
