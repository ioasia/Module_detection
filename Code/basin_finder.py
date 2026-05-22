"""
Basin of attraction finder for protein abundance landscapes on a network.

For each sample, protein abundances are used as a scalar function on the
network graph. Proteins "flow" toward local minima; the set of proteins
that reach the same minimum forms one module (basin of attraction).

Run twice — once with and once without --no-negate — to capture modules
defined by both high and low abundance.

Usage:
    python Code/basin_finder.py \\
        --network=Data/my_dataset_igraph_predictions_sig.txt \\
        --proteomics=Data/my_proteomics.txt \\
        --output-dir=Data \\
        --prefix=my_dataset

    python Code/basin_finder.py ... --no-negate   # positive-abundance basins

Required arguments:
    --network       TSV with Protein1 / Protein2 columns (network edge list)
    --proteomics    TSV, proteins x samples, first column = protein names

Optional arguments:
    --output-dir    directory for output pickle files (default: Data)
    --prefix        prefix for output filenames (default: dataset)
    --no-negate     use raw abundances instead of negated (default: negate)

Outputs (in output-dir):
    {prefix}_basins_neg.pkl   basins from negated abundances  (default run)
    {prefix}_basins.pkl       basins from raw abundances      (--no-negate run)
"""

from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
import pickle as pkl
import argparse
from tqdm import tqdm

# ---- argument parsing -------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--network",     required=True,
                    help="TSV edge list with Protein1 / Protein2 columns")
parser.add_argument("--proteomics",  required=True,
                    help="TSV abundance matrix (proteins x samples, first col = protein names)")
parser.add_argument("--output-dir",  default="Data",
                    help="Directory for output pickle files (default: Data)")
parser.add_argument("--prefix",      default="dataset",
                    help="Prefix for output filenames (default: dataset)")
parser.add_argument("--no-negate",   action="store_true",
                    help="Use raw abundances (default: negate so local minima = high abundance)")
args = parser.parse_args()

negate = not args.no_negate

# ---- load data --------------------------------------------------------------
df = pd.read_csv(args.network, sep='\t')
edges = df[['Protein1', 'Protein2']]

proteomics_quant = pd.read_csv(args.proteomics, sep='\t', index_col=0)
proteins = proteomics_quant.index.tolist()
samples  = list(proteomics_quant.columns)

# ---- build graph filtered to quantified proteins ----------------------------
G_all = nx.from_pandas_edgelist(edges, source='Protein1', target='Protein2')
G     = G_all.subgraph(proteins).copy()

# ---- basin finder helpers ---------------------------------------------------
def sign(x):
    return (x > 0) - (x < 0)

def lower_neighbors(node, G, func):
    return [nbr for nbr in G.neighbors(node) if func[nbr] < func[node]]

def lower_neighbors_same_sign(node, G, func):
    return [nbr for nbr in G.neighbors(node)
            if (func[nbr] < func[node]) and (sign(func[node]) == sign(func[nbr]))]

def basin_finder(G, func):
    """
    Iteratively assigns each node to the basin of its local minimum.

    Returns
    -------
    crit    : dict  iteration -> {critical_node: [basin members]}
    blocked : dict  iteration -> [nodes that could not be assigned]
    """
    G_curr  = G.copy()
    crit    = {}
    blocked = {}
    k = 0
    while G_curr.number_of_nodes() > 0:
        crit_dict = {'boundary': []}
        blocked[k] = []
        for prot in sorted(func, key=func.get):
            if prot not in G_curr.nodes:
                continue
            matching = [key for key, val_list in crit_dict.items()
                        if set(lower_neighbors_same_sign(prot, G_curr, func)) & set(val_list)]
            if len(lower_neighbors(prot, G_curr, func)) == 0:
                crit_dict[prot] = [prot]
            elif len(matching) == 0:
                blocked[k].append(prot)
            elif len(matching) == 1:
                crit_dict[matching[0]].append(prot)
            else:
                crit_dict['boundary'].append(prot)
        crit[k] = crit_dict
        k += 1
        G_curr = G_curr.subgraph(crit_dict['boundary'])
        print('number of nodes: ', len(G_curr.nodes()))
    return crit, blocked

# ---- run basin finder for all samples ---------------------------------------
all_landscapes = {}
all_blocked    = {}

for sample in tqdm(samples):
    print(sample)
    proteomics_sample = proteomics_quant[sample].dropna()

    func = -proteomics_sample if negate else proteomics_sample
    func = func.to_dict()

    proteins_in_sample = set(proteomics_sample.index)
    subG = G.subgraph(proteins_in_sample).copy()
    largest_cc_nodes = max(nx.connected_components(subG), key=len)
    largest_subG     = subG.subgraph(largest_cc_nodes).copy()

    crit, blocked = basin_finder(largest_subG, func)
    all_landscapes[sample] = crit
    all_blocked[sample]    = blocked

# ---- save -------------------------------------------------------------------
output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

suffix   = "_neg" if negate else ""
outpath  = output_dir / f"{args.prefix}_basins{suffix}.pkl"
outpath2 = output_dir / f"{args.prefix}_blocked{suffix}.pkl"

with open(outpath, 'wb') as f:
    pkl.dump(all_landscapes, f)

with open(outpath2, 'wb') as f:
    pkl.dump(all_blocked, f)

print(f"Basins saved to {outpath}")
