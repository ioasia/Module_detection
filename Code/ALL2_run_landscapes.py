from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import pickle as pkl
import matplotlib.colors as mcolors
import statistics
import random
from tqdm import tqdm
import time

df = pd.read_csv('~/GitHub/Module_detection/Data/ALL_cell_line_merged_igraph_predictions_sig.txt', sep='\t')    
edges = df[['Protein1', 'Protein2']]# probably not needed, but nice to have dataframe with only edges
proteomics_quant = pd.read_csv('~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Ioannis Siavelis files - Protein_landscape/ALL_cell_lines/ALL_cell_line_proteomics_median_centered.txt',
                               delimiter = '\t')

proteins = proteomics_quant.index.tolist()
samples = list(proteomics_quant.columns)

# Create an undirected graph
G_all = nx.from_pandas_edgelist(edges, source='Protein1', target='Protein2')
G = G_all.subgraph(proteins).copy()
# Function to find neighbors with lower function values
def sign(x):
    return (x > 0) - (x < 0)
#  option=""): #potential option not to flow to differently signed nodes
def lower_neighbors(node, G, func):
    return [nbr for nbr in G.neighbors(node) if func[nbr] < func[node]]
def lower_neighbors_same_sign(node, G, func):
    return [nbr for nbr in G.neighbors(node) if (func[nbr] < func[node]) & ((sign(func[node]) - sign(func[nbr])) == 0)]
#             func[nbr] < func[node] & sign(func[node] == sign(func[nbr])]
largest_cc_nodes = max(nx.connected_components(G), key=len)
largest_cc = G.subgraph(largest_cc_nodes).copy()

#################### run basin finder for all samples #####################

# sample = samples[17]
all_landscapes = {}
all_blocked = {}
all_flows = {}
for sample in tqdm(samples):
    print(sample) # choose a sample
    proteomics_quant_sample = proteomics_quant[sample].dropna()

    
    func = -proteomics_quant_sample #option to take negative
    func = func.to_dict()
    
    # find non-degenerate critical nodes and attraction basins
    
    # Step 1. Get protein names from your sample dataset
    proteins_in_sample = set(proteomics_quant_sample.index)
    
    # Step 2. Filter graph to keep only those proteins
    subG = G.subgraph(proteins_in_sample).copy()
    
    # Step 3. Find the largest connected component
    largest_cc_nodes = max(nx.connected_components(subG), key=len)
    
    # Step 4. Create the graph of the largest component
    largest_subG = subG.subgraph(largest_cc_nodes).copy()


    G_curr = largest_subG.copy()
    crit = {}
    blocked = {}
    k = 0
    flow = {}
    while G_curr.number_of_nodes() > 0:
        crit_dict = {}
        flow_dict = {}
        boundary = []
        blocked[k] = []
        for prot in sorted(func, key=func.get):
    #         print(prot)
            if prot in G_curr.nodes:
    #             print(prot)
                if (len(lower_neighbors(prot, G_curr, func)) == 0):# option to restrict only to driver genes:  & (prot in driver_genes)
                    crit_dict[prot] = [prot]
                elif len([key for key, value_list in crit_dict.items() if set(lower_neighbors_same_sign(prot, G_curr, func)) & set(value_list)]) == 0:
                    blocked[k].append(prot)
                elif len([key for key, value_list in crit_dict.items() if set(lower_neighbors_same_sign(prot, G_curr, func)) & set(value_list)]) == 1:
                    crit_dict[[key for key, value_list in crit_dict.items() if set(lower_neighbors_same_sign(prot, G_curr, func)) & set(value_list)][0]].append(prot)
                elif len([key for key, value_list in crit_dict.items() if set(lower_neighbors_same_sign(prot, G_curr, func)) & set(value_list)]) > 1:
                    boundary.append(prot)   
                    flow_dict[prot] = [key for key, value_list in crit_dict.items() if set(lower_neighbors_same_sign(prot, G_curr, func)) & set(value_list)]
        crit[k] = crit_dict     
        flow[k] = flow_dict
        k += 1
        G_curr = G_curr.subgraph(boundary)
#         G_curr.remove_nodes_from([elem for lst in crit_dict.values() for elem in lst])
        print('number of nodes: ',len(G_curr.nodes()))
    all_landscapes[sample] = crit
    all_blocked[sample] = blocked
    all_flows[sample] = flow
# Save
code_dir = Path(__file__).resolve().parent          # .../Module_detection/Code
repo_dir = code_dir.parent                          # .../Module_detection

outpath = repo_dir / "Data" / "ALL2_sign_change_condition_negated_basins.pkl"

with open(outpath, 'wb') as f:
    pkl.dump(all_landscapes, f)

outpath = repo_dir / "Data" / "ALL2_sign_change_condition_negated_blocked_prots.pkl"

with open(outpath, 'wb') as g: 
    pkl.dump(all_blocked, g)

outpath = repo_dir / "Data" / "ALL2_sign_change_condition_negated_flows.pkl"

with open(outpath, 'wb') as g: 
    pkl.dump(all_flows, g)


