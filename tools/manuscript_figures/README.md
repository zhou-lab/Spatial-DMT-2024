# Data-driven manuscript figures

`build_data_figures.py` generates quantitative manuscript Figures 2-5 from
existing revised pipeline outputs and audited published/comparison tables.
Conceptual Figure 1 remains in `results/manuscript_figures/generated_schematics`.

Run with the project environment:

```bash
micromamba run -n qc_env python tools/manuscript_figures/build_data_figures.py \
  --repo-root /home/suc3/Lab/Spatial-DMT-2024 \
  --output-dir /home/suc3/Lab/Spatial-DMT-2024/results/manuscript_figures/data_driven
```

Each figure is exported as PDF, SVG, PNG, and LZW-compressed TIFF. Tidy metric
and provenance TSV files are written beside the figures. Missing values remain
missing and technical/control runs are excluded.

Figure 2 currently summarizes real diagnostic classifications across completed
runs. It does not claim controlled classification accuracy because a labeled
synthetic or manually adjudicated validation set has not yet been generated.

Figure 4 uses the audited E6.5 same-input comparison available now. A matched
per-pixel concordance panel is intentionally omitted until equivalent
per-pixel outputs are available from both processing implementations.
