import random
import networkx as nx
import csv
import os

FILENAME1 = "graph_query.csv"
FILENAME2 = "graph_target.csv"
LABELS = 10

# Function to save graph to CSV
def save_graph_to_csv(graph, filename):
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Number of vertices " + str(graph.number_of_nodes())])
        writer.writerow(['u', 'v', 'l1', 'l2'])  # Header
        for edge in graph.edges():
            u, v = edge
            l1 = graph.nodes[u]['label']
            l2 = graph.nodes[v]['label']
            writer.writerow([u, v, l1, l2])

# Graph initialization
g1 = nx.gnp_random_graph(550, 0.55, 42)
g2 = nx.gnp_random_graph(550, 0.55, 42)

for node in g1.nodes():
    label = random.randint(0, LABELS - 1)
    g1.nodes[node]["label"] = label
    g2.nodes[node]["label"] = label

# Save graphs to CSV files
dirname = os.path.dirname(__file__)

save_graph_to_csv(g1, os.path.join(dirname, FILENAME1))
save_graph_to_csv(g2, os.path.join(dirname, FILENAME2)) 

