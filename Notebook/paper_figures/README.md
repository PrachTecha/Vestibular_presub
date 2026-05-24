# Vestibular scRNA-seq — Code for Figures 1 & 2

Reproducible notebooks for the figure panels listed below.

> *Deep scRNAseq reveals novel heterogeneity among adult VGNs and VHCs.*

## Figure 1 — VGN + VHC heterogeneity

| Figure | Description | Notebook |
| --- | --- | --- |
| Fig. 1d | Feature plots of VHC type and striolar/extrastriolar marker genes | [`Fig1d.ipynb`](Fig1d.ipynb) |
| Fig. 1h | Matrix map (DEG + SCENIC GRN) across VHC and VGN cell types | [`Fig1h.ipynb`](Fig1h.ipynb) |
| Suppl. Fig. 1e | Self-marker dot plots of VHC and VGN subtypes | [`SupplFig1e.ipynb`](SupplFig1e.ipynb) |

> Suppl. Fig. 1e is a **self-marker / staircase dotplot** where each cluster (e.g. `HC_Foxp2`) is paired with the gene it is named after (Foxp2). `Slc17a7` is added on the VGN side as the shared Vglut1 marker. Cluster labels are renamed in Illustrator: `HC_*` → `VHC-I.*` / `VHC-II.*`, `neuron_*` → `VGN-*`.

## Figure 2 — VHC functional and spatial diversity, vulnerability index

| Figure | Description | Notebook |
| --- | --- | --- |
| Fig. 2f | Dot plot of mech/voltage-gated ion channels across VHC subtypes | [`Fig2f.ipynb`](Fig2f.ipynb) |
| Fig. 2g | Metabolic profiles of VHC subtypes (scCellFie radial plot) | [`Fig2g.ipynb`](Fig2g.ipynb) |
| Fig. 2h | Schematic of vulnerability index workflow (Illustrator only — no code) | — |
| Fig. 2i–j | Vulnerability index validation in bulk + scRNA-seq IHC/OHC | [`Fig2hij.ipynb`](Fig2hij.ipynb) |
| Fig. 2k–l | Vulnerability index across VHC subtypes | [`Fig2kl.ipynb`](Fig2kl.ipynb) |

## Figure 4 — VGN hierarchical molecular profile

| Figure | Description | Notebook |
| --- | --- | --- |
| Fig. 4a | Metacluster-marker staircase dot plot (4 genes × 13 VGN subtypes) | [`Fig4a.ipynb`](Fig4a.ipynb) |
| Fig. 4c | Biological processes per metacluster (balanced consensus NMF + GO:BP) | [`Fig4c.ipynb`](Fig4c.ipynb) |

> Fig 4a uses the original **4-MC** marker panel `[Lypd1, Sall3, Vcan, Lmo3]`. Fig 4c uses the consolidated **3-MC** scheme (Sall3+ folded into the other three) — this is intentional, see Fig4a.ipynb header for the rationale.

## Environment

Tested with Python 3.10. Key packages:

- `scanpy`, `anndata`, `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`
- `loompy`, `PyComplexHeatmap` — Fig 1h
- `openpyxl` — reads `.xlsx` inputs
- `sccellfie` — Fig 2g

Tested versions in [`requirements.txt`](requirements.txt).

## File layout

```
github_upload/
├── README.md               # this file
├── requirements.txt
├── .gitignore
├── _gene_sets.py           # shared metabolic / proteostasis panels
├── _dotplot.py             # shared marker_dotmap_simple helper
├── Fig1d.ipynb
├── Fig1h.ipynb
├── SupplFig1e.ipynb
├── Fig2f.ipynb
├── Fig2g.ipynb
├── Fig2hij.ipynb
├── Fig2kl.ipynb
├── Fig4a.ipynb
└── Fig4c.ipynb
```

Each notebook writes PDFs and CSVs to a local `figures/` subdirectory (created on first run).
