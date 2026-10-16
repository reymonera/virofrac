from ete3 import Tree
#import networkx as nx
import numpy as np
from decimal import Decimal, getcontext

# Gives a precise number of decimals for the final
# result. 20 digits of precision are determined.
getcontext().prec = 20

####################################################
#           PHYLOGENETIC - BASED                   #
####################################################

# This function processes the lines from the OTU
# table that is included in the input files. Mainly
# it works by determining community assignment per
# taxa. It ignores the header.
def get_processed_otu_table(lines):
    taxon_communities = {}
    header = None
    
    for i, line in enumerate(lines):
        columns = line.split('\t')
        
        if i == 0:
            header = columns
            continue
        
        taxon_name = columns[0]
        count_community_1 = float(columns[1])
        count_community_2 = float(columns[2])
        
        if count_community_1 > 0 and count_community_2 > 0:
            taxon_communities[taxon_name] = (header[1], header[2])
        elif count_community_1 > 0 and count_community_2 == 0:
            taxon_communities[taxon_name] = (header[1],)
        elif count_community_2 > 0:
            taxon_communities[taxon_name] = (header[2],)
    
    return taxon_communities

# This function puts a new feature on the tree. In this
# case, it refers to the community feature. The output
# is a tree with leaves and the featured community
# they belong to.

# def put_community_feature_on_tree(tree, lines):
#     taxon_communities = get_processed_otu_table(lines)
    
#     for leaf in tree.iter_leaves():
#         if leaf.name in taxon_communities:
#             leaf.community = taxon_communities[leaf.name]
    
#     return tree

def put_community_feature_on_tree(tree, lines):
    taxon_communities = get_processed_otu_table(lines)
    for leaf in tree.iter_leaves():
        leaf.community = taxon_communities.get(leaf.name, ())
    return tree

# This function gets all the branches that are unique.
# It generates a Python list that will be filled with
# the lengths of said unique branches. This function
# is not considering the root of the tree.
def get_unique_branches(tree):
    unique_branches = []
    community_list = []

    # for node in tree.traverse():
    #     if node.is_root(): 
    #         continue

    for node in tree.traverse():
        if node.is_root(): 
            continue
        descendants = node.get_leaves()
        for leaf in descendants:
            community_list.append(leaf.community)
        if len(set(community_list)) == 1 and len(community_list[0]) == 1:
            unique_branches.append(node)
        community_list = []
    return(unique_branches)

# getSumUniqueBranchLengths just sums all the 
# unique lengths that were in the previous
# created list.
def get_sum_unique_branch_lengths(tree):
    unique_branches = get_unique_branches(tree)
    unique_lengths = []
    for node in unique_branches:
        unique_lengths.append(node.dist)
    return(sum(unique_lengths))

# getSumAllBranchLengths just sums all the lengths.
# This function does not consider the root.
def get_sum_all_branch_lengths(tree):
    all_lengths = []

    for node in tree.traverse():
        if node.is_root():
            continue
        all_lengths.append(node.dist)

    return(sum(all_lengths))

# getCommunityDictionaries is a function that will
# create dictionaries for the communities that are
# being compared with each other.
def get_community_dictionaries(lines):
    first_community = {}
    second_community = {}
    community_dictionaries = []
    for i, line in enumerate(lines):

        if i == 0:
            continue
        
        columns = line.split('\t')
        
        first_community[columns[0]] = float(columns[1])
        second_community[columns[0]] = float(columns[2])
    
    community_dictionaries.append(first_community)
    community_dictionaries.append(second_community)

    return community_dictionaries

# getCommunityProportionsDictionaries is a function
# that is assigning relative proportions per taxa 
# in each community. 
def get_community_proportions_dictionaries(lines):
    communities_proportions = []
    community_dictionaries = get_community_dictionaries(lines)

    first_community = community_dictionaries[0]
    second_community = community_dictionaries[1]
    
    all_first_community = sum(first_community.values())
    all_second_community = sum(second_community.values())
    
    first_community.update((x, y/all_first_community) for x, y in first_community.items())
    second_community.update((x, y/all_second_community) for x, y in second_community.items())
    
    communities_proportions.append(first_community)
    communities_proportions.append(second_community)
    
    return communities_proportions

# This function outputs the Unweighted Unifrac value
# between 2 communities.
def get_unweighted_unifrac(community_tree, lines):
    tree = put_community_feature_on_tree(community_tree, lines)

    sum_unique_branches = Decimal(get_sum_unique_branch_lengths(tree))
    sum_all_branches = Decimal(get_sum_all_branch_lengths(tree))

    return sum_unique_branches/sum_all_branches

# This function outputs the Normalized Weighted Unifrac 
# value between 2 communities.
def get_unnormalized_weighted_unifrac(community_tree, lines):
    tree = put_community_feature_on_tree(community_tree, lines)

    first_community = get_community_proportions_dictionaries(lines)[0]
    second_community = get_community_proportions_dictionaries(lines)[1]

    sum_weighted_branches = []
    
    for node in tree.traverse():
        if node.is_root():
            continue
            
        descendants = node.get_leaves()
        sum_first_proportions = []
        sum_second_proportions = []
        
        for leaf in descendants:
            sum_first_proportions.append(first_community[leaf.name])
            sum_second_proportions.append(second_community[leaf.name])
        
        p_i = sum(sum_first_proportions)
        q_i = sum(sum_second_proportions)
        
        sum_weighted_branches.append(node.dist * abs(p_i - q_i))

    return(Decimal(sum(sum_weighted_branches)))

# This function outputs the Normalized Weighted Unifrac 
# value between 2 communities.
def get_normalized_weighted_unifrac(community_tree, lines):
    tree = put_community_feature_on_tree(community_tree, lines)

    first_community = get_community_proportions_dictionaries(lines)[0]
    second_community = get_community_proportions_dictionaries(lines)[1]

    sum_weighted_branches = []
    sum_all_weighted_branches = []
    
    for node in tree.traverse():
        if node.is_root():
            continue
            
        descendants = node.get_leaves()
        sum_first_proportions = []
        sum_second_proportions = []
        
        for leaf in descendants:
            sum_first_proportions.append(first_community[leaf.name])
            sum_second_proportions.append(second_community[leaf.name])
        
        p_i = sum(sum_first_proportions)
        q_i = sum(sum_second_proportions)
        
        sum_weighted_branches.append(node.dist * abs(p_i - q_i))
        sum_all_weighted_branches.append(node.dist * (p_i + q_i))
    
    return(Decimal(sum(sum_weighted_branches)/sum(sum_all_weighted_branches)))

####################################################
#              NETWORK - BASED
####################################################

# Extracts the edges of the A∪B subnetwork for one pair: both nodes with abundance > 0 in A or B
# (the same filter as before). Only looks at edges leaving A∪B nodes, never the whole network.
# Returns, per relevant edge, its weight and the membership/abundance of both nodes in A and B.
def get_pair_edges(network, column_A, column_B):
    abundance_A = network.abundance[:, column_A].toarray().ravel()
    abundance_B = network.abundance[:, column_B].toarray().ravel()
    in_A = abundance_A > 0
    in_B = abundance_B > 0
    in_union = in_A | in_B

    candidate_edges = network.edges_by_source[np.flatnonzero(in_union)].indices
    relevant_edges = candidate_edges[in_union[network.edge_target[candidate_edges]]]

    node1 = network.edge_source[relevant_edges]
    node2 = network.edge_target[relevant_edges]

    return {
        'weight': network.edge_weight[relevant_edges],
        'node1': node1,
        'node2': node2,
        'n1_in_A': in_A[node1],
        'n1_in_B': in_B[node1],
        'n2_in_A': in_A[node2],
        'n2_in_B': in_B[node2],
        'node1_A': abundance_A[node1],
        'node1_B': abundance_B[node1],
        'node2_A': abundance_A[node2],
        'node2_B': abundance_B[node2],
        'abundance_A': abundance_A,
        'abundance_B': abundance_B,
    }

# def get_distances_from_edges(edge_list):
#     distance_list = [1 - x for x in edge_list]
#     return distance_list

# def get_similarities_from_edges(edge_list):
#     similarity_list = [x for x in edge_list]
#     return similarity_list

# def get_all_weights_from_edges(network, community_A, community_B):            
#     all_weights = []
#     for node1, node2, weight in network.edges(data='weight'):
#         n1_in_A = community_A in network.nodes[node1]['communities']
#         n1_in_B = community_B in network.nodes[node1]['communities']
#         n2_in_A = community_A in network.nodes[node2]['communities']
#         n2_in_B = community_B in network.nodes[node2]['communities']
        
#         # Edge is relevant if BOTH nodes belong to at least one of the communities
#         if (n1_in_A or n1_in_B) and (n2_in_A or n2_in_B):
#             all_weights.append(weight)
    
#     return all_weights

# # Get volume for any imput community, which is always the first one
# def get_volume_of_first_sample(network, community_A, community_B):            
#     edge_weights_of_sample_A = []
#     for node1, node2, weight in network.edges(data='weight'):
#         n1_in_A = community_A in network.nodes[node1]['communities']
#         n1_in_B = community_B in network.nodes[node1]['communities']
#         n2_in_A = community_A in network.nodes[node2]['communities']
#         n2_in_B = community_B in network.nodes[node2]['communities']
        
#         # Edge must be relevant (both nodes belong to at least one community)
#         if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
#             continue
        
#         # Volume of A = edges where at least one node belongs to A
#         if n1_in_A or n2_in_A:
#             edge_weights_of_sample_A.append(weight)
    
#     return sum(edge_weights_of_sample_A)

# def get_monochromatic_edges(network, community_A, community_B):
#     monochromatic_weights = []
#     for node1, node2, weight in network.edges(data='weight'):
#         n1_in_A = community_A in network.nodes[node1]['communities']
#         n1_in_B = community_B in network.nodes[node1]['communities']
#         n2_in_A = community_A in network.nodes[node2]['communities']
#         n2_in_B = community_B in network.nodes[node2]['communities']
        
#         # First: Edge must be relevant (same filter as all_edges)
#         if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
#             continue
        
#         # Second: Check if monochromatic
#         # Both nodes in A only (neither in B)
#         both_A_only = (n1_in_A and n2_in_A and not n1_in_B and not n2_in_B)
#         # Both nodes in B only (neither in A)
#         both_B_only = (n1_in_B and n2_in_B and not n1_in_A and not n2_in_A)
        
#         #if (n1_in_A and n2_in_A) or (n1_in_B and n2_in_B):
#         if both_A_only or both_B_only:
#             monochromatic_weights.append(weight)
    
#     return monochromatic_weights

# def get_bichromatic_edges(network, community_A, community_B):
#     bichromatic_weights = []
#     for node1, node2, weight in network.edges(data='weight'):
#         n1_in_A = community_A in network.nodes[node1]['communities']
#         n1_in_B = community_B in network.nodes[node1]['communities']
#         n2_in_A = community_A in network.nodes[node2]['communities']
#         n2_in_B = community_B in network.nodes[node2]['communities']
        
#         # First: Edge must be relevant (same filter as all_edges)
#         if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
#             continue

#         both_in_both = (n1_in_A and n1_in_B and n2_in_A and n2_in_B)
#         original_bichromatic = (n1_in_A != n2_in_A or n1_in_B != n2_in_B)
        
#         if original_bichromatic or both_in_both:
#             bichromatic_weights.append(weight)
    
#     return bichromatic_weights

# def get_total_for_first_community(network, community_A, community_B):
#     node_abundances = {}
    
#     for node1, node2 in network.edges():
#         n1_in_A = community_A in network.nodes[node1]['communities']
#         n1_in_B = community_B in network.nodes[node1]['communities']
#         n2_in_A = community_A in network.nodes[node2]['communities']
#         n2_in_B = community_B in network.nodes[node2]['communities']
        
#         if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
#             continue
        
#         node_abundances[node1] = network.nodes[node1]['communities'].get(community_A, 0)
#         node_abundances[node2] = network.nodes[node2]['communities'].get(community_A, 0)
    
#     return sum(node_abundances.values())

# def get_unweighted_netunifrac(network, community_A, community_B):
#     monochromatic_edges = get_monochromatic_edges(network, community_A, community_B)
#     all_edges = get_all_weights_from_edges(network, community_A, community_B)

#     monochromatic_distances = get_distances_from_edges(monochromatic_edges)
#     all_distances = get_distances_from_edges(all_edges)

#     sum_monochromatic_distances = sum(monochromatic_distances)
#     sum_all_distances = sum(all_distances)

#     if sum_all_distances == 0:
#         # No shared structure = maximum distance (completely different communities)
#         return 1.0

#     # NetUniFrac
#     unweighted_netunifrac = sum_monochromatic_distances/sum_all_distances

#     return unweighted_netunifrac

# def get_weighted_netunifrac(network, community_A, community_B):
#     list_num = []
#     list_dem = []

#     total_community_A = get_total_for_first_community(network, community_A, community_B)
#     total_community_B = get_total_for_first_community(network, community_B, community_A)

#     for node1, node2, weight in network.edges(data='weight'):
#         n1_in_A = community_A in network.nodes[node1]['communities']
#         n1_in_B = community_B in network.nodes[node1]['communities']
#         n2_in_A = community_A in network.nodes[node2]['communities']
#         n2_in_B = community_B in network.nodes[node2]['communities']

#         node1_A = network.nodes[node1]['communities'].get(community_A, 0)
#         node1_B = network.nodes[node1]['communities'].get(community_B, 0)
#         node2_A = network.nodes[node2]['communities'].get(community_A, 0)
#         node2_B = network.nodes[node2]['communities'].get(community_B, 0)

#         if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
#             continue

#     # Go through all distances
#     # Multiply them by the proportions of node wieghts
#         if total_community_A == 0 or total_community_B == 0:
#             continue
#         else:
#             # The (1 - weight) here is doing the distance transformation, instead of using the function assigned
#             # for this transformation as in the get_weighted_netunifrac() function
#             num = (1 - weight)*abs((node1_A + node2_A)/total_community_A - (node1_B + node2_B)/total_community_B)
#             dem = (1 - weight)*abs((node1_A + node2_A)/total_community_A + (node1_B + node2_B)/total_community_B)
#         # Save this on a list
#             list_num.append(num)
#             list_dem.append(dem)
#     # Return the division
#     # Handle empty lists
#     if sum(list_dem) == 0:
#         weighted_netunifrac = 1
#     else:
#         weighted_netunifrac = sum(list_num) / sum(list_dem)
    
#     return weighted_netunifrac

# Weighted but with monochromatic edges
#def get_weighted_netunifrac(network, community_A, community_B):
    # list_num = []
    # list_dem = []
    
    # total_community_A = get_total_for_first_community(network, community_A, community_B)
    # total_community_B = get_total_for_first_community(network, community_B, community_A)
    
    # if total_community_A == 0 or total_community_B == 0:
    #     return 1.0
    
    # for node1, node2, weight in network.edges(data='weight'):
    #     n1_in_A = community_A in network.nodes[node1]['communities']
    #     n1_in_B = community_B in network.nodes[node1]['communities']
    #     n2_in_A = community_A in network.nodes[node2]['communities']
    #     n2_in_B = community_B in network.nodes[node2]['communities']
        
    #     if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
    #         continue
        
    #     node1_A = network.nodes[node1]['communities'].get(community_A, 0)
    #     node1_B = network.nodes[node1]['communities'].get(community_B, 0)
    #     node2_A = network.nodes[node2]['communities'].get(community_A, 0)
    #     node2_B = network.nodes[node2]['communities'].get(community_B, 0)
        
    #     prop_diff = abs((node1_A + node2_A)/total_community_A - 
    #                     (node1_B + node2_B)/total_community_B)
    #     prop_sum  = abs((node1_A + node2_A)/total_community_A + 
    #                     (node1_B + node2_B)/total_community_B)
        
    #     # Monochromatic check inline
    #     both_A_only = (n1_in_A and n2_in_A and not n1_in_B and not n2_in_B)
    #     both_B_only = (n1_in_B and n2_in_B and not n1_in_A and not n2_in_A)
    #     is_mono     = both_A_only or both_B_only
        
    #     if is_mono:
    #         list_num.append((1 - weight) * prop_diff)
        
    #     list_dem.append((1 - weight) * prop_sum)
    
    # if sum(list_dem) == 0:
    #     return 1.0
    
    # return sum(list_num) / sum(list_dem)
    
# def put_edge_community_data(network):
#     edge_data = []
#     for node1, node2, weight in network.edges(data='weight'):
#         communities_1 = network.nodes[node1].get('communities', {})
#         communities_2 = network.nodes[node2].get('communities', {})
#         edge_data.append((communities_1, communities_2, weight))

#     return edge_data

# This function performs the Spectre distance. It
# requires similarity instead of distance as an input. If the denominator
# is equal to 0, then it will automatically return 1.
# def get_based_spectral_clustering(edge_data, community_A, community_B):
#     sum_bichromatic = 0.0
#     volume_A = 0.0
#     volume_B = 0.0

#     for communities_1, communities_2, weight in edge_data:
#         n1_in_A = community_A in communities_1
#         n1_in_B = community_B in communities_1
#         n2_in_A = community_A in communities_2
#         n2_in_B = community_B in communities_2

#         # Edge must be relevant
#         if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
#             continue

#         # Volume A: at least one node in A
#         if n1_in_A or n2_in_A:
#             volume_A += weight

#         # Volume B: at least one node in B
#         if n1_in_B or n2_in_B:
#             volume_B += weight

#         # Bichromatic: nodes differ in membership, or both in both
#         both_in_both = (n1_in_A and n1_in_B and n2_in_A and n2_in_B)
#         original_bichromatic = (n1_in_A != n2_in_A or n1_in_B != n2_in_B)
#         if original_bichromatic or both_in_both:
#             sum_bichromatic += weight

#     if volume_A == 0 or volume_B == 0:
#         return 1.0

#     return 1 - ((sum_bichromatic / volume_A + sum_bichromatic / volume_B) * 0.5)

# def get_based_fuzzy_spectral_clustering(edge_data, community_A, community_B):
    # total_A = 0.0
    # total_B = 0.0
    # for communities_1, communities_2, weight in edge_data:
    #     for communities in [communities_1, communities_2]:
    #         if community_A in communities:
    #             total_A += communities[community_A]
    #         if community_B in communities:
    #             total_B += communities[community_B]

    # if total_A == 0 or total_B == 0:
    #     return 1.0

    # sum_bichromatic = 0.0
    # volume_A = 0.0
    # volume_B = 0.0

    # for communities_1, communities_2, weight in edge_data:
    #     n1_in_A = community_A in communities_1
    #     n1_in_B = community_B in communities_1
    #     n2_in_A = community_A in communities_2
    #     n2_in_B = community_B in communities_2

    #     if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
    #         continue

    #     node1_A = communities_1.get(community_A, 0)
    #     node1_B = communities_1.get(community_B, 0)
    #     node2_A = communities_2.get(community_A, 0)
    #     node2_B = communities_2.get(community_B, 0)

    #     delta1 = node1_A - node1_B
    #     delta2 = node2_A - node2_B
    #     factor = weight * ((1 - delta1 * delta2) / 2)

    #     # Volume A: weighted by max proportion of A across both nodes
    #     if n1_in_A or n2_in_A:
    #         sum1 = node1_A + node1_B
    #         sum2 = node2_A + node2_B
    #         prop1_A = node1_A / sum1 if sum1 > 0 else 0.0
    #         prop2_A = node2_A / sum2 if sum2 > 0 else 0.0
    #         volume_A += weight * max(prop1_A, prop2_A)

    #     # Volume B: weighted by max proportion of B across both nodes
    #     if n1_in_B or n2_in_B:
    #         sum1 = node1_A + node1_B
    #         sum2 = node2_A + node2_B
    #         prop1_B = node1_B / sum1 if sum1 > 0 else 0.0
    #         prop2_B = node2_B / sum2 if sum2 > 0 else 0.0
    #         volume_B += weight * max(prop1_B, prop2_B)

    #     both_in_both   = (n1_in_A and n1_in_B and n2_in_A and n2_in_B)
    #     is_bichromatic = (n1_in_A != n2_in_A or n1_in_B != n2_in_B) or both_in_both

    #     if is_bichromatic:
    #         sum_bichromatic += factor

    # if volume_A == 0 or volume_B == 0:
    #     return 1.0
    
    # fs_distance = 1 - 0.5 * (sum_bichromatic / volume_A + sum_bichromatic / volume_B)
    # return 1.0 if fs_distance < 0 else fs_distance
    # #return max(0.0, fs_distance)

    # #return 1 - 0.5 * (sum_bichromatic / volume_A + sum_bichromatic / volume_B)

# Monochromatic edges: both nodes only in A, or both only in B.
def get_monochromatic_mask(pair_edges):
    n1_in_A, n1_in_B = pair_edges['n1_in_A'], pair_edges['n1_in_B']
    n2_in_A, n2_in_B = pair_edges['n2_in_A'], pair_edges['n2_in_B']

    both_A_only = n1_in_A & n2_in_A & ~n1_in_B & ~n2_in_B
    both_B_only = n1_in_B & n2_in_B & ~n1_in_A & ~n2_in_A

    return both_A_only | both_B_only

# Bichromatic edges: the nodes differ in membership, or both are in both samples.
def get_bichromatic_mask(pair_edges):
    n1_in_A, n1_in_B = pair_edges['n1_in_A'], pair_edges['n1_in_B']
    n2_in_A, n2_in_B = pair_edges['n2_in_A'], pair_edges['n2_in_B']

    both_in_both = n1_in_A & n1_in_B & n2_in_A & n2_in_B
    original_bichromatic = (n1_in_A != n2_in_A) | (n1_in_B != n2_in_B)

    return original_bichromatic | both_in_both

# Unweighted NetUniFrac: sum of monochromatic distances / sum of all distances (distance = 1 - weight).
# If the denominator is 0, returns 1.0.
def get_unweighted_netunifrac(pair_edges):
    distances = 1 - pair_edges['weight']
    sum_all_distances = distances.sum()

    if sum_all_distances == 0:
        return 1.0

    monochromatic = get_monochromatic_mask(pair_edges)

    return float(monochromatic.sum() / len(distances)) #this is the count version, idk what to do with this.
    #return float(distances[monochromatic].sum() / sum_all_distances)

# Weighted NetUniFrac. Totals = abundance over the nodes of the A∪B subnetwork, each node once.
# Each edge contributes (1 - weight) * |pA - pB| to the numerator and (1 - weight) * |pA + pB| to the denominator.
def get_weighted_netunifrac(pair_edges):
    subnetwork_nodes = np.unique(np.concatenate([pair_edges['node1'], pair_edges['node2']]))
    total_community_A = pair_edges['abundance_A'][subnetwork_nodes].sum()
    total_community_B = pair_edges['abundance_B'][subnetwork_nodes].sum()

    if total_community_A == 0 or total_community_B == 0:
        return 1

    distances = 1 - pair_edges['weight']
    proportion_A = (pair_edges['node1_A'] + pair_edges['node2_A']) / total_community_A
    proportion_B = (pair_edges['node1_B'] + pair_edges['node2_B']) / total_community_B

    sum_num = (distances * np.abs(proportion_A - proportion_B)).sum()
    sum_dem = (distances * np.abs(proportion_A + proportion_B)).sum()

    if sum_dem == 0:
        return 1

    return float(sum_num / sum_dem)

# Spectre distance (uses similarity). Volumes = weights of edges touching A (or B);
# bichromatic = weights of bichromatic edges. If either volume is 0, returns 1.0.
def get_based_spectral_clustering(pair_edges):
    weight = pair_edges['weight']

    volume_A = weight[pair_edges['n1_in_A'] | pair_edges['n2_in_A']].sum()
    volume_B = weight[pair_edges['n1_in_B'] | pair_edges['n2_in_B']].sum()

    if volume_A == 0 or volume_B == 0:
        return 1.0

    sum_bichromatic = weight[get_bichromatic_mask(pair_edges)].sum()

    return float(1 - ((sum_bichromatic / volume_A + sum_bichromatic / volume_B) * 0.5))

# Fuzzy spectre: like spectre, but each edge is weighted by the fuzzy factor (delta with raw counts)
# and the volumes by each node's maximum proportion of A (or B). If the result is < 0, returns 1.0.
# The old pre-check over the whole network is dropped: it implied volume = 0, which already returns 1.0.
def get_based_fuzzy_spectral_clustering(pair_edges):
    weight = pair_edges['weight']
    node1_A, node1_B = pair_edges['node1_A'], pair_edges['node1_B']
    node2_A, node2_B = pair_edges['node2_A'], pair_edges['node2_B']

    delta1 = node1_A - node1_B
    delta2 = node2_A - node2_B
    factor = weight * ((1 - delta1 * delta2) / 2)

    sum1 = node1_A + node1_B
    sum2 = node2_A + node2_B
    prop1_A = np.divide(node1_A, sum1, out=np.zeros_like(sum1), where=sum1 > 0)
    prop2_A = np.divide(node2_A, sum2, out=np.zeros_like(sum2), where=sum2 > 0)
    prop1_B = np.divide(node1_B, sum1, out=np.zeros_like(sum1), where=sum1 > 0)
    prop2_B = np.divide(node2_B, sum2, out=np.zeros_like(sum2), where=sum2 > 0)

    touches_A = pair_edges['n1_in_A'] | pair_edges['n2_in_A']
    touches_B = pair_edges['n1_in_B'] | pair_edges['n2_in_B']
    volume_A = (weight * np.maximum(prop1_A, prop2_A))[touches_A].sum()
    volume_B = (weight * np.maximum(prop1_B, prop2_B))[touches_B].sum()

    if volume_A == 0 or volume_B == 0:
        return 1.0

    sum_bichromatic = factor[get_bichromatic_mask(pair_edges)].sum()
    fs_distance = 1 - 0.5 * (sum_bichromatic / volume_A + sum_bichromatic / volume_B)

    return 1.0 if fs_distance < 0 else float(fs_distance)