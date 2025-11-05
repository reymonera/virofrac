import frac
from itertools import combinations
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

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

# def get_input_otu_lines(partial_otu_tables):
#     """Convierte DataFrames a líneas usando NumPy (muy rápido)"""
#     otu_table_lines_list = []
    
#     print(f"Convirtiendo {len(partial_otu_tables)} tablas parciales...")
#     for i, partial_table in enumerate(partial_otu_tables):
#         if i % 2000 == 0:
#             print(f"  Tabla {i}/{len(partial_otu_tables)}...")
        
#         lines = []
        
#         # Header
#         lines.append('\t'.join(partial_table.columns))
        
#         # Datos con NumPy (mucho más rápido que iterrows)
#         data_array = partial_table.values.astype(str)
#         lines.extend(['\t'.join(row) for row in data_array])
        
#         otu_table_lines_list.append(lines)
    
#     return otu_table_lines_list

# This function is managing the comparisons necessary
# for the matrix output. For this, the function first
# calculates the size of the matrix and then builds
# the matrix based on the selected distance.
def get_frac_matrix_output(otu_df, tree, unifrac_type):
    sample_names = otu_df.columns[1:].tolist()
    n = len(sample_names)

    matrix = np.zeros((n, n))

    for (i,j) in combinations(range(n), 2):
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

# def get_frac_matrix_output(otu_df, tree, unifrac_type):
#     from itertools import combinations
    
#     print(f"Calculando matriz de distancias...")
#     sample_names = otu_df.columns[1:].tolist()
#     n = len(sample_names)
    
#     print(f"Total de pares a procesar: {n*(n-1)//2}")
    
#     matriz = np.zeros((n, n))
    
#     # Procesar cada par directamente sin almacenar todas las tablas
#     for idx, (i, j) in enumerate(combinations(range(n), 2)):
#         if idx % 500 == 0:
#             print(f"  Par {idx}/{n*(n-1)//2}...")
        
#         # Crear sub-tabla solo para este par
#         col_i = sample_names[i]
#         col_j = sample_names[j]
#         sub_df = otu_df[[otu_df.columns[0], col_i, col_j]]
        
#         # Convertir a líneas directamente
#         lines = ['\t'.join(sub_df.columns)]
#         data_array = sub_df.values.astype(str)
#         lines.extend(['\t'.join(row) for row in data_array])
        
#         # Calcular distancia
#         if unifrac_type == 'normalized weighted unifrac':
#             distance = frac.getNormalizedWeightedUniFrac(tree, lines)
#         elif unifrac_type == 'unnormalized weighted unifrac':
#             distance = frac.getUnnormalizedWeightedUniFrac(tree, lines)
#         else:  # unweighted
#             distance = frac.getUnweightedUnifrac(tree, lines)
        
#         matriz[i, j] = distance
#         matriz[j, i] = distance
    
#     return pd.DataFrame(matriz, index=sample_names, columns=sample_names)

def get_dataframe_from_matrix(matrix, otu_df):
    sample_names = otu_df.columns[1:].tolist()
    matrix_as_df = pd.DataFrame(matrix, index=sample_names, columns=sample_names)

    matrix_as_df.to_csv('matrix_as_dataframe.tsv', sep='\t')

    return matrix_as_df

# AGREGAR STRIP DE COLOR STRIP DE COLOR
def get_heatmap_output(matrix, otu_df, unifrac_type):
    matrix_df = get_dataframe_from_matrix(matrix, otu_df)
    dist_condensed = squareform(matrix)
    
    row_linkage = linkage(dist_condensed, method='average')
    col_linkage = row_linkage 
    
    g = sns.clustermap(
        matrix,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cmap='coolwarm',
        vmin=0,
        vmax=1,
        annot=False,
        fmt='.3f',
        figsize=(10, 8),
        cbar_kws={
            'label': 'UniFrac Distance',
            'orientation': 'horizontal',
            },
        linewidths=0.5,
        linecolor='white',
        xticklabels=matrix_df.columns,
        yticklabels=matrix_df.index 
    )
    
    g.ax_heatmap.set_xlabel('Samples', fontsize=12)
    g.ax_heatmap.set_ylabel('Samples', fontsize=12)
    plt.suptitle(f'{(unifrac_type).title()} Distance Matrix with Hierarchical Clustering',
             fontsize=14, y=0.95)
    
    plt.tight_layout()
    plt.savefig('virofrac_heatmap.png', dpi=300, bbox_inches='tight')
    plt.show()

    return g