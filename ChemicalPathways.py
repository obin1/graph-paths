import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from ChemicalCase import ChemicalCase


# Load in data -- note pandas is cleaner than CSV reader in this case
# If not using this repository, make sure you have the edge list 
# as downloaded from https://github.com/obin1/graph-paths in the same folder
df = pd.read_csv("gckpp_EdgeList.csv")

# Create DiGraph
B = nx.DiGraph()

# Add edges
B.add_edges_from([(df["from"][i],df["to"][i]) for i in range(0,len(df["to"]))])

# load in a data example: assumed local, though a full path can be specified
filename = "Amazon_L1_20180101_2100.txt"
amazon_surface = ChemicalCase(filename)

# Loop over all reaction nodes
for i in range(1,914):
    # For all edges coming in and out of a reaction node, set to half the rxn timescale 
    for u,v in B.in_edges("R"+str(i)):
        B[u][v]['timescale'] = 0.5 / (amazon_surface.reaction_rates[i-1] + 1e-20) # 1e-20 to avoid divide by 0
    for u,v in B.out_edges("R"+str(i)):
        B[u][v]['timescale'] = 0.5 / (amazon_surface.reaction_rates[i-1] + 1e-20)

# get a list of species nodes from B.nodes
# species nodes are all nodes that don't have the form "R" and then an integer
rxn_node_list = ['R' + str(i) for i in range(1,914)]
spc_node_list = [node for node in list(B.nodes) if node not in rxn_node_list]


# Then use Lucas's nested for loop approach, but this time with the weighted call
# And just find shortest paths from species to species
shortest_paths = []   # Loop through all nodes in the graph as target nodes
shortest_timescales = []
count = 0
for source_node in spc_node_list: # Ensure the source node is not the same as the target node
    for target_node in spc_node_list: # Check if there exists a path from the source node to the target node
        if source_node != target_node and nx.has_path(B, source_node, target_node):
            count += 1
            # delete last count in the same print line
            print('\rCalculating shortest path #', str(count), end='')
            path = nx.shortest_path(B, source=source_node, target=target_node,weight='timescale')
            timescale = nx.path_weight(B, path,weight='timescale')
            shortest_paths.append(path)
            shortest_timescales.append(timescale)

sns.kdeplot(shortest_timescales, log_scale=True)
plt.xlabel("Chemical Pathway Timescale [seconds/molec/cc]")
plt.savefig("AmazonShortestPaths.pdf")

            