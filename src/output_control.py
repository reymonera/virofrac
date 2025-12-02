import src.frac as frac
from itertools import combinations
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from tqdm import tqdm
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
from src.utils import GlobalTimer

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
def get_frac_matrix_output(otu_df, tree, unifrac_type):
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
        
        if unifrac_type == 'normalized weighted unifrac':
            distance = frac.getNormalizedWeightedUniFrac(tree, lines)
        elif unifrac_type == 'unnormalized weighted unifrac':
            distance = frac.getUnnormalizedWeightedUniFrac(tree, lines)
        else:
            distance = frac.getUnweightedUnifrac(tree, lines)
        
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
        #print("Using default colormap gradient: 'coolwarm'")
        GlobalTimer.log("Using default colormap gradient: 'coolwarm'")
        return 'coolwarm'
    
    if ',' in color_gradient:
        colors = [c.strip() for c in color_gradient.split(',')]
        return LinearSegmentedColormap.from_list('custom', colors)
    
    if isinstance(color_gradient, str):
        if color_gradient.startswith('#'):
            #print('Color pair not specified, using white in the custom gradiet...')
            GlobalTimer.log("Color pair not specified, using white in the custom gradiet...")
            return LinearSegmentedColormap.from_list('custom', ['#FFFFFF', color_gradient])
        else:
            #print(f"Using predefined colormap: '{color_gradient}'")
            GlobalTimer.log("Using predefined colormap")
            return color_gradient
    

# This function controls the heatmap output. There are two 
# options: When there is no metadata input and when there 
# is metadata input. A legend and a color column needs to 
# be specified for a correct labelling.
def get_heatmap_output(matrix, otu_df, unifrac_type, color_gradient, metadata_file=None, legend_columns=None, color_columns=None):
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
            f'{(unifrac_type).title()} Distance Matrix with Hierarchical Clustering',
            fontsize=14,
            x=0.45,
            y=0.96
        )
        
        plt.savefig('virofrac_heatmap.png', dpi=300, bbox_inches='tight')
        GlobalTimer.log("✅ Heatmap ready! Showing the final plot...")
        #print("✅ Heatmap ready! Showing the final plot...")
        plt.show()
        GlobalTimer.log("✅ Heatmap ready! Showing the final plot...")
        #print(f"✅ Heatmap saved in virofrac_heatmap.png")
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
            GlobalTimer.log("⚠️ Legend column not found, skipping")
            #print("⚠️ Legend column not found, skipping")
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
                GlobalTimer.log("⚠️ Color column not found, using automatic colors")
                #print("⚠️ Color column not found, using automatic colors")
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
        GlobalTimer.log("⚠️ No valid columns, generating simple heatmap")
        #print("⚠️ No valid columns, generating simple heatmap")
        return get_heatmap_output(matrix, otu_df, unifrac_type, None, None, None)

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
        f'{(unifrac_type).title()} Distance Matrix with Hierarchical Clustering',
        fontsize=14,
        x=0.45,
        y=0.96
    )
    
    # Output files being handled here.
    plt.savefig('virofrac_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig('virofrac_heatmap.svg', format='svg', bbox_inches='tight')

    GlobalTimer.log("✅ Heatmap ready! Showing the final plot...")
    #print("✅ Heatmap ready! Showing the final plot...")
    plt.show()
    GlobalTimer.log("✅ Heatmap saved in virofrac_heatmap.png")
    #print("✅ Heatmap saved in virofrac_heatmap.png")
    
    return g


