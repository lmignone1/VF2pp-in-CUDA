import networkx as nx
import csv, os

FILENAME1 = "graph_query_example.csv"
FILENAME2 = "graph_target_example.csv"

def save_graph_to_csv(graph : nx, filename):
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
g1 = nx.Graph()
g2 = nx.Graph()

# Aggiunta dei nodi e degli archi al primo grafo
g1.add_edges_from([(0, 1), (0, 2), (1, 2), (1, 3)])     # ogni id del nodo è stato decrementato di 1

# Aggiunta dei nodi e degli archi al secondo grafo
g2.add_edges_from([(0,1), (0,2), (1,2), (2,3)])
# Verifica dell'isomorfismo


# VF2++ initialization
g1.nodes[0]["label"] = 0
g1.nodes[1]["label"] = 2
g1.nodes[2]["label"] = 1
g1.nodes[3]["label"] = 1

# G2.nodes[0]["label"] = 2
# G2.nodes[1]["label"] = 0
g2.nodes[0]["label"] = 1
g2.nodes[1]["label"] = 0
g2.nodes[2]["label"] = 2
g2.nodes[3]["label"] = 1


G1_labels = nx.get_node_attributes(g1, "label")
G2_labels = nx.get_node_attributes(g2, "label")

# print(G1_labels)
# print(G2_labels)

dirname = os.path.dirname(__file__)

save_graph_to_csv(g1, os.path.join(dirname, FILENAME1))
save_graph_to_csv(g2, os.path.join(dirname, FILENAME2)) 


# VF2++
# t0 = time.time()
# m = vf2pp_isomorphism(g1, g2, node_label="label")
# print('mapping')
# print(m)

# print(f"VF2++ elapsed time: {time.time() - t0}")

