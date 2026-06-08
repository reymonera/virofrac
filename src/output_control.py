import src.frac as frac
from src.utils import GlobalTimer
from itertools import combinations
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
from tqdm import tqdm
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
from concurrent.futures import ThreadPoolExecutor, as_completed
from ete3 import Tree
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import pcoa

# This function manages the partial OTU tables that will
# be used as an input when compairing each pair.
def get_partial_otu_tables(otu_df):
    
    otu_header = otu_df.columns[1:].tolist()
    partial_otu_tables = []

    for couple in combinations(otu_header, 2):

        i, j = couple
        partial_df = otu_df[[otu_df.columns[0], i, j]]
        partial_otu_tables.append(partial_df)

    return partial_otu_tables

# This function outputs the lines that are part of the
# partial OTU tables.
def get_input_otu_lines(partial_otu_tables):
    otu_table_lines_list = []

    for partial_table in partial_otu_tables:
        lines = []
        lines.append('\t'.join(partial_table.columns))
        data_array = partial_table.values.astype(str)
        lines.extend(['\t'.join(row) for row in data_array])

        otu_table_lines_list.append(lines)

    return otu_table_lines_list

# This function is managing the comparisons necessary
# for the matrix output. For this, the function first
# calculates the size of the matrix and then builds
# the matrix based on the selected distance.
def get_frac_matrix_output(otu_df, tree, distance_type):
    sample_names = otu_df.columns[1:].tolist()
    n = len(sample_names)
    total_pairs = n * (n - 1) // 2
    
    matrix = np.zeros((n, n))
    
    for (i, j) in tqdm(combinations(range(n), 2), 
                       total=total_pairs,
                       desc=f"[{GlobalTimer.elapsed():7.2f}s] Brewing UniFrac",
                       unit="pairs"):
        
        col_i = sample_names[i]
        col_j = sample_names[j]
        sub_df = otu_df[[otu_df.columns[0], col_i, col_j]]
        
        lines = ['\t'.join(sub_df.columns)]
        data_array = sub_df.values.astype(str)
        lines.extend(['\t'.join(row) for row in data_array])
        
        if distance_type == 'normalized weighted unifrac':
            distance = frac.get_normalized_weighted_unifrac(tree, lines)
        elif distance_type == 'unnormalized weighted unifrac':
            distance = frac.get_unnormalized_weighted_unifrac(tree, lines)
        else:
            distance = frac.get_unweighted_unifrac(tree, lines)
        
        matrix[i, j] = distance
        matrix[j, i] = distance
    
    return pd.DataFrame(matrix, index=sample_names, columns=sample_names)

# This function is managing the comparisons necessary
# for the matrix output. For this, the function first
# calculates the size of the matrix and then builds
# the matrix based on the selected distance.
def get_net_frac_matrix_output(network, distance_type, count_table):
    count_table = count_table.set_index(count_table.columns[0])
    sample_names = count_table.columns.tolist()
    n = len(sample_names)
    total_pairs = n * (n - 1) // 2
    matrix = np.zeros((n, n))

    # Precompute once
    if distance_type == 'spectral clustering':
        edge_data = frac.put_edge_community_data(network)

    for (i, j) in tqdm(combinations(range(n), 2),
                        total=total_pairs,
                        desc=f"[{GlobalTimer.elapsed():7.2f}s] Brewing UniFrac",
                        unit="pairs"):
        sample_i = sample_names[i]
        sample_j = sample_names[j]

        if distance_type == 'normalized weighted unifrac':
            distance = frac.get_weighted_netunifrac(network, sample_i, sample_j)
        elif distance_type == 'spectral clustering':
            distance = frac.get_based_spectral_clustering(edge_data, sample_i, sample_j)
        else:
            distance = frac.get_unweighted_netunifrac(network, sample_i, sample_j)

        matrix[i, j] = distance
        matrix[j, i] = distance

    return pd.DataFrame(matrix, index=sample_names, columns=sample_names)

# This function saves the matrix in a tab delimited file
# that should be available after the run of this pipeline.
def get_dataframe_from_matrix(matrix, otu_df):
    sample_names = otu_df.columns[1:].tolist()
    matrix_as_df = pd.DataFrame(matrix, index=sample_names, columns=sample_names)

    GlobalTimer.log("Saving matrix as a dataframe...")

    matrix_as_df.to_csv('matrix_as_dataframe.tsv', sep='\t')

    return matrix_as_df

# This function controls the heatmap gradient for the
# final plot. It can accept 2 hexacodes as input for
# a custom color scale. Palettes are also available.
def get_color_gradient_for_heatmap(color_gradient):
    if not color_gradient:
        GlobalTimer.log("Using default colormap gradient: 'coolwarm'")
        return 'coolwarm'
    
    if ',' in color_gradient:
        colors = [c.strip() for c in color_gradient.split(',')]
        return LinearSegmentedColormap.from_list('custom', colors)
    
    if isinstance(color_gradient, str):
        if color_gradient.startswith('#'):
            GlobalTimer.log("Color pair not specified, using white in the custom gradiet...")
            return LinearSegmentedColormap.from_list('custom', ['#FFFFFF', color_gradient])
        else:
            GlobalTimer.log("Using predefined colormap")
            return color_gradient

# ESTA FUNCIÓN ESTÁ DE RESIDUO PERO PUEDE QUE FUNCIONE
def get_plot_network_output(network):
    """Plot the network and save to file."""
    import networkx as nx
    
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(network, k=0.5, iterations=50, seed=42)
    edges = network.edges(data=True)
    weights = [d.get('weight', 1.0) for _, _, d in edges]
    
    nx.draw_networkx_nodes(network, pos, node_size=20, node_color='steelblue', alpha=0.8)
    nx.draw_networkx_edges(network, pos, width=0.5, alpha=0.5, 
                           edge_color=weights, edge_cmap=plt.cm.coolwarm)
    
    plt.title(f"ANI Network ({network.number_of_nodes()} nodes, {network.number_of_edges()} edges)")
    plt.axis('off')
    plt.tight_layout()
    plt.savefig("virofrac_network_plot.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Network plot saved to: virofrac_network_plot.png")


# --- PERMUTATION TESTS BASED ON CLUSTERING-SIGNIFICANCE-TEST GITHUB REPO ---
# Script by Yasas1994 at https://github.com/Yasas1994/Clustering-significance-test

# This function builds a newick tree from node and labels.
# It will return a formatted version of this.
def build_newick(node, labels):
    if node.is_leaf():
        return labels[node.id]
    left = build_newick(node.left, labels)
    right = build_newick(node.right, labels)
    dist = node.dist

    return f"({left}:{dist},{right}:{dist})"

# This function converts the matrix to a proper
# ete3 tree.
def scipy_linkage_to_ete3(linkage_matrix, labels):
    root, _ = to_tree(linkage_matrix, rd=True)
    newick = build_newick(root, labels) + ";"
    tree = Tree(newick, format=1)

    return tree

# SIMULATE() FUNCTION FROM CLUSTERING-SIGNIFICANCE-TEST GITHUB REPO
def simulate(annot, column):
    annot = annot.copy()
    annot['tmp'] = annot[column].sample(frac=1).values
    tmp2 = annot.groupby('cluster')['tmp'].value_counts().reset_index()
    tmp2.columns = ['cluster', 'tmp', 'counts']
    sum2 = tmp2.groupby('cluster')['counts'].sum()
    max2 = tmp2.groupby('cluster')['counts'].max()

    return sum(max2) / sum(sum2)

# RUN_PERMUTATION_TEST() FUNCTION BASED ON CLUSTERING-SIGNIFICANCE-TEST GITHUB REPO
def run_permutation_test(tree, metadata, column, replicates=1000, min_leaves=3, workers=10):
    clust_n = 0
    clust_leaf = []
    for branch in tree.traverse("postorder"):
        if not branch.is_leaf():
            tmp_ = []
            for leaf in branch.get_leaves():
                tmp_.append([clust_n, leaf.name])
            if len(tmp_) > min_leaves:
                clust_leaf.extend(tmp_)
                clust_n += 1

    if not clust_leaf:
        GlobalTimer.log(f"ERROR: No clusters found for permutation test on {column}.")
        return 1.0

    clust_leaf = pd.DataFrame(clust_leaf, columns=['cluster', 'leaf_lab'])
    annot = pd.merge(
        clust_leaf,
        metadata[[column]].reset_index(),
        left_on='leaf_lab',
        right_on=metadata.index.name or metadata.reset_index().columns[0]
    )

    # Observed global purity
    tmp = annot.groupby('cluster')[column].value_counts().reset_index()
    tmp.columns = ['cluster', column, 'counts']
    sum_cluster = tmp.groupby('cluster')['counts'].sum()
    max_per_cluster = tmp.groupby('cluster')['counts'].max()
    observed_purity = sum(max_per_cluster) / sum(sum_cluster)

    simulated = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(simulate, annot, column) for _ in range(replicates)]
        for future in as_completed(futures):
            simulated.append(future.result())

    p_value = sum(s >= observed_purity for s in simulated) / replicates

    return p_value


# Adds global p-value under color strips
# one for each legend column.
def add_significance_labels_to_heatmap(g, row_linkage, labels, metadata, legend_columns, replicates=1000):
    tree = scipy_linkage_to_ete3(row_linkage, labels)

    pos = g.ax_col_colors.get_position()
    bar_height = (pos.y1 - pos.y0) / len(legend_columns)

    for col_idx, legend_col in enumerate(legend_columns):
        if legend_col not in metadata.columns:
            continue

        GlobalTimer.log(f"Running permutation test for {legend_col}...")
        p_value = run_permutation_test(tree, metadata, legend_col, replicates=replicates)
        GlobalTimer.log(f"Permutation test for {legend_col}: p={p_value} (n={replicates} replicates)")
        star = '*' if p_value < 0.05 else ''
        label = f'p={p_value:.3f}{star}'

        x_pos = pos.x1 + 0.01
        y_pos = pos.y1 - (col_idx + 0.5) * bar_height

        g.figure.text(
            x_pos, y_pos,
            label,
            fontsize=8,
            va='center',
            ha='left',
            rotation=0
        )

# --- END OF PERMUTATION FUNCTIONS ---

# This function controls the heatmap output. There are two 
# options: When there is no metadata input and when there 
# is metadata input. A legend and a color column needs to 
# be specified for a correct labelling.
def get_heatmap_output(matrix, otu_df, distance_type, color_gradient, metadata_file=None, legend_columns=None, color_columns=None):
    gradient = get_color_gradient_for_heatmap(color_gradient)
    matrix_df = get_dataframe_from_matrix(matrix, otu_df)
    dist_condensed = squareform(matrix)
    
    row_linkage = linkage(dist_condensed, method='average')
    col_linkage = row_linkage

    if not metadata_file or not legend_columns:
        g = sns.clustermap(
            matrix,
            row_linkage=row_linkage,
            col_linkage=col_linkage,
            cmap=gradient,
            vmin=0,
            vmax=1,
            annot=False,
            fmt='.3f',
            figsize=(15, 13),
            cbar_kws={
                'label': 'UniFrac Distance',
                'orientation': 'horizontal'
            },
            linewidths=0.5,
            linecolor='white',
            xticklabels=matrix_df.columns,
            yticklabels=matrix_df.index,
            dendrogram_ratio=(0.08, 0.08)
        )

        g.ax_heatmap.set_xlabel('Samples', fontsize=12)
        g.ax_heatmap.set_ylabel('Samples', fontsize=12)
        g.figure.suptitle(
            f'{(distance_type).title()} Distance Matrix with Hierarchical Clustering',
            fontsize=14,
            x=0.45,
            y=0.96
        )
        
        plt.savefig('virofrac_heatmap.png', dpi=300, bbox_inches='tight')
        GlobalTimer.log("✓ Heatmap ready! Showing the final plot...")
        plt.show()
        GlobalTimer.log("✓ Heatmap ready! Showing the final plot...")
        return g
    
    with open(metadata_file, 'r') as f:
        first_line = f.readline()
        if '\t' in first_line:
            sep = '\t'
        elif ',' in first_line:
            sep = ','
        else:
            sep = None
    
    metadata = pd.read_csv(metadata_file, sep=sep) if sep else pd.read_csv(metadata_file, sep=None, engine='python')

    sample_col = 'sample' if 'sample' in metadata.columns.str.lower() else metadata.columns[0]
    metadata = metadata.set_index(sample_col)
    metadata = metadata.loc[matrix_df.index]

    colors_combined = pd.DataFrame(index=matrix_df.index)
    legend_handles = []

    for idx, legend_col in enumerate(legend_columns):
        if legend_col not in metadata.columns:
            GlobalTimer.log("⚠ Legend column not found, skipping")
            continue
    
        values = metadata[legend_col]
        unique_values = sorted(values.dropna().unique())

        if color_columns and idx < len(color_columns):
            color_col = color_columns[idx]
            
            if color_col in metadata.columns:
                color_map = {}
                for val in unique_values:
                    mask = metadata[legend_col] == val
                    colors_for_val = metadata.loc[mask, color_col].dropna().unique()
                    if len(colors_for_val) > 0:
                        color_map[val] = colors_for_val[0]
                    else:
                        fallback_idx = len(color_map) % 10
                        color_map[val] = sns.color_palette('tab10')[fallback_idx]
            else:
                GlobalTimer.log("⚠ Color column not found, using automatic colors")
                palette = sns.color_palette('tab10', n_colors=len(unique_values))
                color_map = dict(zip(unique_values, palette))
        else:
            palette = sns.color_palette('tab10', n_colors=len(unique_values))
            color_map = dict(zip(unique_values, palette))
    
        values_filled = values.fillna('__NA__')
        color_map['__NA__'] = '#FFFFFF'
        row_colors = values_filled.map(color_map)
        colors_combined[legend_col] = row_colors

        if idx < len(legend_columns) - 1:
            colors_combined[f'_space_{idx}'] = '#FFFFFF'
        
        legend_items = [Patch(facecolor=color_map[val], label=str(val)) 
                       for val in unique_values]
        legend_handles.append((legend_col, legend_items))

    if len(colors_combined.columns) == 0:
        GlobalTimer.log("⚠ No valid columns, generating simple heatmap")
        return get_heatmap_output(matrix, otu_df, distance_type, None, None, None)

    g = sns.clustermap(
        matrix,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        row_colors=colors_combined,
        col_colors=colors_combined,
        cmap=gradient,
        vmin=0,
        vmax=1,
        annot=False,
        figsize=(15, 13),
        cbar_kws={'label': 'UniFrac Distance'},
        linewidths=0,
        xticklabels=matrix_df.columns,
        yticklabels=matrix_df.index,
        dendrogram_ratio=(0.08, 0.08),
        cbar_pos=(0.02, 0.92, 0.05, 0.15)
    )
    
    g.ax_row_colors.set_xticklabels([])
    g.ax_col_colors.set_yticklabels([])

    for spine in g.ax_row_colors.spines.values():
        spine.set_visible(False)

    for spine in g.ax_col_colors.spines.values():
        spine.set_visible(False)

    g.ax_row_colors.tick_params(left=False, right=False, top=False, bottom=False)
    g.ax_col_colors.tick_params(left=False, right=False, top=False, bottom=False)
    g.figure.subplots_adjust(right=0.8)

    legend_ax = g.figure.add_axes([0.85, 0.3, 0.15, 0.4])
    legend_ax.axis('off')

    y_start = 0.7
    legend_height = 0.25

    for idx, (col_name, handles) in enumerate(legend_handles):
        y_pos = y_start - (idx * legend_height)
        legend_ax = g.figure.add_axes([0.85, y_pos, 0.12, legend_height])
        legend_ax.axis('off')
        
        legend_ax.legend(
            handles=handles,
            title=col_name,
            frameon=True,
            fontsize=9,
            loc='upper left',
            title_fontsize=10
        )

    g.ax_heatmap.set_xlabel('Samples', fontsize=12)
    g.ax_heatmap.set_ylabel('Samples', fontsize=12)
    
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, ha='right', fontsize=6)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=6)
    
    g.figure.subplots_adjust(top=0.93) 
    g.figure.suptitle(
        f'{(distance_type).title()} Distance Matrix with Hierarchical Clustering',
        fontsize=14,
        x=0.45,
        y=0.96
    )

    # Llamando a add_significance_labels_to_heatmap para agregar
    # p-values debajo de cada barra de colores
    labels = list(matrix_df.index)
    add_significance_labels_to_heatmap(g, row_linkage, labels, metadata, legend_columns)

    plt.savefig('virofrac_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig('virofrac_heatmap.svg', format='svg', bbox_inches='tight')

    # Necesito poner esto como opcón.
    #GlobalTimer.log("✓ Heatmap ready! Showing the final plot...")
    #plt.show()
    GlobalTimer.log("✓ Heatmap saved in virofrac_heatmap.png")
    
    return g

# These functions control the PCoA output. There are two functions,
# one specified for performing the permutation tests and the other
# one that actually does the plot. The PCoA can be implemented
# alongside an EnvFit test, in which vectors are drawn on the
# PCoA plot.

def get_envfit(vals, x, y, n_perm=999, seed=42):
    r_x = np.corrcoef(vals, x)[0, 1]
    r_y = np.corrcoef(vals, y)[0, 1]
    r2_obs = r_x**2 + r_y**2
    rng = np.random.default_rng(seed)
    r2_perm = np.array([
        np.corrcoef(v := rng.permutation(vals), x)[0, 1]**2 + np.corrcoef(v, y)[0, 1]**2
        for _ in range(n_perm) 
    ])
    pval = (np.sum(r2_perm >= r2_obs) + 1) / (n_perm + 1)

    return r_x, r_y, r2_obs, pval

def get_pcoa_output(matrix, otu_df, distance_type, metadata_file=None, legend_columns=None, color_columns=None, env_columns=None):
    env_r2_min=0.20 # Minimum R2 for vectors
    env_pval_max=0.05 # Minimum p-value for vectors
    n_perm=999 # Might take time in big datasets
    
    if not metadata_file or not legend_columns:
        GlobalTimer.log("⚠ WARNING: PCoA requires metadata and at least one legend column. Skipping.")
        return None
    
    # Loading matrix...
    matrix_df = get_dataframe_from_matrix(matrix, otu_df)
    sample_names = list(matrix_df.index)

    with open(metadata_file, 'r') as f:
        first_line = f.readline()
    sep = '\t' if '\t' in first_line else (',' if ',' in first_line else None)
    metadata = (
        pd.read_csv(metadata_file, sep=sep)
        if sep
        else pd.read_csv(metadata_file, sep=None, engine='python')
    )
    sample_col = (
        metadata.columns[metadata.columns.str.lower() == 'sample'][0]
        if any(metadata.columns.str.lower() == 'sample')
        else metadata.columns[0]
    )
    metadata = metadata.set_index(sample_col).loc[sample_names]

    # PCoA
    dm = DistanceMatrix(matrix, ids=sample_names)
    pcoa_results = pcoa(dm)
    pcoa_coords = pcoa_results.samples
    explained   = pcoa_results.proportion_explained
    pc1_var = round(explained.iloc[0] * 100, 1)
    pc2_var = round(explained.iloc[1] * 100, 1)

    # EnvFit
    envfit_vectors = {}

    if env_columns:
        missing = [c for c in env_columns if c not in metadata.columns]
        if missing:
            GlobalTimer.log(f"⚠ env_columns not found in metadata and will be skipped: {missing}")

        for col in env_columns:
            if col not in metadata.columns:
                continue
            vals_series = metadata[col].dropna()
            common = vals_series.index.intersection(pcoa_coords.index)
            if len(common) < 3:
                GlobalTimer.log(f"⚠ '{col}' has fewer than 3 common samples. Skipping.")
                continue

            vals = vals_series.loc[common].values
            x    = pcoa_coords.loc[common, 'PC1'].values
            y    = pcoa_coords.loc[common, 'PC2'].values

            r_x, r_y, r2, pval = get_envfit(vals, x, y, n_perm=n_perm)

            if r2 < env_r2_min:
                GlobalTimer.log(f"  envfit '{col}': r²={r2:.3f} < {env_r2_min} → skipped (weak)")
                continue
            if pval > env_pval_max:
                GlobalTimer.log(f"  envfit '{col}': p={pval:.3f} > {env_pval_max} → skipped (not significant)")
                continue

            envfit_vectors[col] = (r_x, r_y, r2, pval)
            GlobalTimer.log(f"  envfit '{col}': r²={r2:.3f}, p={pval:.3f} ✓")

    # Colors
    color_maps = {}
    for idx, legend_col in enumerate(legend_columns):
        if legend_col not in metadata.columns:
            GlobalTimer.log(f"⚠ Legend column '{legend_col}' not found, skipping.")
            continue

        values       = metadata[legend_col]
        unique_values = sorted(values.dropna().unique())

        if color_columns and idx < len(color_columns):
            color_col = color_columns[idx]
            if color_col in metadata.columns:
                color_map = {}
                for val in unique_values:
                    mask = metadata[legend_col] == val
                    cols_for_val = metadata.loc[mask, color_col].dropna().unique()
                    color_map[val] = cols_for_val[0] if len(cols_for_val) > 0 else sns.color_palette('tab10')[len(color_map) % 10]
            else:
                GlobalTimer.log(f"⚠ Color column '{color_col}' not found, using automatic colors.")
                color_map = dict(zip(unique_values, sns.color_palette('tab10', n_colors=len(unique_values))))
        else:
            color_map = dict(zip(unique_values, sns.color_palette('tab10', n_colors=len(unique_values))))

        color_maps[legend_col] = color_map

    valid_legend_cols = list(color_maps.keys())
    if not valid_legend_cols:
        GlobalTimer.log("⚠ No valid legend columns found. Skipping PCoA.")
        return None

    # Plot
    n_plots = len(valid_legend_cols)
    fig, axes = plt.subplots(1, n_plots, figsize=(7 * n_plots, 6), squeeze=False)

    # Vector being projected at 0.75%
    coord_range = max(
        pcoa_coords['PC1'].abs().max(),
        pcoa_coords['PC2'].abs().max()
    )
    arrow_scale = coord_range * 0.75

    for col_idx, legend_col in enumerate(valid_legend_cols):
        ax        = axes[0][col_idx]
        color_map = color_maps[legend_col]
        values    = metadata[legend_col]

        # Scatter per group
        for val in sorted(color_map.keys()):
            mask = values == val
            samples_in_group = values[mask].index.intersection(pcoa_coords.index)
            if samples_in_group.empty:
                continue
            ax.scatter(
                pcoa_coords.loc[samples_in_group, 'PC1'],
                pcoa_coords.loc[samples_in_group, 'PC2'],
                c=[color_map[val]], label=str(val),
                s=70, edgecolors='white', linewidths=0.5, zorder=3,
            )

        # Muestras sin valor → gris
        na_samples = values[values.isna()].index.intersection(pcoa_coords.index)
        if not na_samples.empty:
            ax.scatter(
                pcoa_coords.loc[na_samples, 'PC1'],
                pcoa_coords.loc[na_samples, 'PC2'],
                c='#CCCCCC', label='NA',
                s=70, edgecolors='white', linewidths=0.5, zorder=3,
            )

        # Envfit Arrows
        for env_col, (r_x, r_y, r2, pval) in envfit_vectors.items():
            length = np.sqrt(r2)
            dx = (r_x / length) * arrow_scale * length
            dy = (r_y / length) * arrow_scale * length

            ax.annotate(
                '', xy=(dx, dy), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='#333333', lw=1.5, mutation_scale=12),
                zorder=5,
            )
            label = f"{env_col}\nr²={r2:.2f}, p={pval:.3f}"
            ax.text(
                dx * 1.12, dy * 1.12, label,
                fontsize=7.5, ha='center', va='center',
                color='#333333', fontweight='bold', zorder=6,
            )

        ax.axhline(0, color='grey', linewidth=0.5, linestyle='--', zorder=1)
        ax.axvline(0, color='grey', linewidth=0.5, linestyle='--', zorder=1)
        ax.set_xlabel(f'PC1 ({pc1_var}%)', fontsize=11)
        ax.set_ylabel(f'PC2 ({pc2_var}%)', fontsize=11)
        ax.set_title(
            f'{distance_type.title()} PCoA\nColored by: {legend_col}',
            fontsize=12, pad=10,
        )
        ax.legend(
            title=legend_col, bbox_to_anchor=(1.02, 1), loc='upper left',
            fontsize=8, title_fontsize=9, frameon=True,
        )
        ax.set_aspect('equal', adjustable='datalim')

    fig.suptitle(f'{distance_type.title()} Distance Matrix — PCoA', fontsize=14, y=1.02)
    fig.tight_layout()
    plt.savefig('virofrac_pcoa.png', dpi=300, bbox_inches='tight')
    #plt.savefig('virofrac_pcoa.svg', format='svg', bbox_inches='tight')
    GlobalTimer.log("✓ PCoA saved → virofrac_pcoa.png / virofrac_pcoa.svg")

    return fig