"""Shared dot-plot helper used by Fig 2f and Suppl Fig 1e.

`marker_dotmap_simple` reproduces the dot-plot style used throughout the
paper: dot size = fraction of cells with expression > `expression_cutoff`,
dot colour = mean expression z-scored per gene (`standard_scale='var'`).
Optional hierarchical clustering of rows (genes) and/or columns (clusters).

Colour palette for the cluster colour bar at the bottom is read from
`adata.uns['<groupby>_colors']` (scanpy convention).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy import sparse
from scipy.cluster.hierarchy import dendrogram, linkage


def marker_dotmap_simple(
    adata,
    genes,
    groupby='cell_type',
    color_map='Blues',
    width_per_column=0.5,
    height_per_row=0.3,
    use_raw=False,
    standard_scale='var',
    expression_cutoff=0.0,
    fontsize=12,
    spines=False,
    show_rownames=True,
    show_colnames=True,
    rows_cluster=False,
    cols_cluster=False,
    cluster_method='average',
    cluster_metric='euclidean',
    show_dendrogram=True,
    dendro_width=0.5,
    dendro_height=0.5,
    min_frac=None,
):
    cell_groups = adata.obs[groupby].cat.categories
    matrix = adata.raw.X if (use_raw and adata.raw is not None) else adata.X

    means, fracs = [], []
    for gene in genes:
        gene_means, gene_fracs = [], []
        gene_idx = adata.var_names.get_loc(gene)
        for group in cell_groups:
            mask = adata.obs[groupby].values == group
            gene_data = matrix[mask][:, gene_idx]
            if sparse.issparse(gene_data):
                gene_data = gene_data.toarray().flatten()
            else:
                gene_data = np.asarray(gene_data).flatten()
            gene_means.append(np.mean(gene_data))
            gene_fracs.append(np.mean(gene_data > expression_cutoff))
        means.append(gene_means)
        fracs.append(gene_fracs)
    means, fracs = np.array(means), np.array(fracs)

    if min_frac is not None:
        keep = fracs.max(axis=1) >= min_frac
        if not keep.any():
            raise ValueError(f'No genes passed min_frac = {min_frac}.')
        means, fracs = means[keep, :], fracs[keep, :]
        genes = [g for g, k in zip(genes, keep) if k]

    if standard_scale == 'var':
        means = (means - means.mean(axis=1, keepdims=True)) / means.std(axis=1, keepdims=True)
    elif standard_scale == 'group':
        means = (means - means.mean(axis=0, keepdims=True)) / means.std(axis=0, keepdims=True)
    means_clean = np.nan_to_num(means, nan=0.0)

    row_linkage = col_linkage = None
    if rows_cluster:
        row_linkage = linkage(means_clean, method=cluster_method, metric=cluster_metric)
        order = dendrogram(row_linkage, no_plot=True)['leaves']
        means, fracs, means_clean = means[order, :], fracs[order, :], means_clean[order, :]
        genes = [genes[i] for i in order]
    if cols_cluster:
        col_linkage = linkage(means_clean.T, method=cluster_method, metric=cluster_metric)
        order = dendrogram(col_linkage, no_plot=True)['leaves']
        means, fracs = means[:, order], fracs[:, order]
        cell_groups = [cell_groups[i] for i in order]

    main_w, main_h = len(cell_groups) * width_per_column, len(genes) * height_per_row
    col_dendro_h = dendro_height if (cols_cluster and show_dendrogram) else 0
    row_dendro_w = dendro_width if (rows_cluster and show_dendrogram) else 0
    fig_w, fig_h = row_dendro_w + main_w + 3, col_dendro_h + 0.5 + main_h + 2
    fig = plt.figure(figsize=(fig_w, fig_h))

    main_left = 0.10 + row_dendro_w / fig_w
    main_bottom = 0.15 / fig_h
    main_w_norm = main_w / fig_w * 0.6
    main_h_norm = main_h / fig_h * 0.7
    ax_main = fig.add_axes([main_left, main_bottom, main_w_norm, main_h_norm])

    if rows_cluster and show_dendrogram:
        ax_rd = fig.add_axes([
            main_left - row_dendro_w / fig_w - 0.01, main_bottom,
            row_dendro_w / fig_w, main_h_norm,
        ])
        dendrogram(row_linkage, orientation='left', ax=ax_rd,
                   color_threshold=0, above_threshold_color='black')
        ax_rd.set_xticks([]); ax_rd.set_yticks([])
        for s in ax_rd.spines.values():
            s.set_visible(False)
        ax_rd.set_ylim(ax_main.get_ylim())

    if cols_cluster and show_dendrogram:
        ax_cd = fig.add_axes([
            main_left, main_bottom + main_h_norm + 0.01,
            main_w_norm, col_dendro_h / fig_h,
        ])
        dendrogram(col_linkage, orientation='top', ax=ax_cd,
                   color_threshold=0, above_threshold_color='black')
        ax_cd.set_xticks([]); ax_cd.set_yticks([])
        for s in ax_cd.spines.values():
            s.set_visible(False)
        xlim = ax_cd.get_xlim()
        scale = (xlim[1] - xlim[0]) / len(cell_groups)
        offset = 5 - scale / 2
        ax_cd.set_xlim(offset - 0.5 * scale, offset + (len(cell_groups) - 0.5) * scale)

    ax_main.grid(True, linestyle='-', alpha=0.2, color='gray')
    for i in range(len(genes)):
        for j in range(len(cell_groups)):
            size = fracs[i, j] * 200
            if size > 0:
                ax_main.scatter(j, i, s=size, c=means[i, j], cmap=color_map,
                                vmin=0, vmax=1, edgecolor='gray', linewidth=0.3)

    if show_colnames:
        ax_main.set_xticks(range(len(cell_groups)))
        ax_main.set_xticklabels(cell_groups, rotation=90, ha='right', fontsize=fontsize - 2)
    if show_rownames:
        ax_main.yaxis.tick_right()
        ax_main.yaxis.set_label_position('right')
        ax_main.set_yticks(range(len(genes)))
        ax_main.set_yticklabels(genes, fontsize=fontsize - 2)
    ax_main.tick_params(axis='both', length=0)
    if not spines:
        for s in ax_main.spines.values():
            s.set_visible(False)

    color_key = f'{groupby}_colors'
    if color_key in adata.uns:
        color_dict = dict(zip(adata.obs[groupby].cat.categories, adata.uns[color_key]))
        cell_colors = [color_dict.get(g, '#808080') for g in cell_groups]
    else:
        cell_colors = plt.cm.tab20(np.linspace(0, 1, len(cell_groups)))
    color_bar_y = -1.0 if not (cols_cluster and show_dendrogram) else -0.8
    for j, color in enumerate(cell_colors):
        ax_main.add_patch(Rectangle((j - 0.5, color_bar_y), 1, 0.4,
                                    facecolor=color, edgecolor='none', clip_on=False))
    ax_main.set_xlim(-0.5, len(cell_groups) - 0.5)
    y_upper = -1.2 if not (cols_cluster and show_dendrogram) else -0.8
    ax_main.set_ylim(len(genes) - 0.5, y_upper)
    return fig
