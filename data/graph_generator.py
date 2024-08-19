# MIGNONE LORENZO 0622701866 L.MIGNONE@STUDENTI.UNISA.IT
# Course: High Performance Computing 2022/2023
# Lecturer: Francesco Moscato	fmoscato@unisa.it

# Copyright (C) 2024 - All Rights Reserved

# This file is part of VF2pp-in-CUDA.

# VF2pp-in-CUDA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# VF2pp-in-CUDA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with VF2pp-in-CUDA.  If not, see <http://www.gnu.org/licenses/>.

import random
import networkx as nx
import csv
import os

class GraphGenerator:
    '''
    This class is used to generate random graphs with different densities and sizes. The graphs are saved in CSV format 
    where each row represents an edge between two vertices. The first column is the source vertex, the second column 
    is the target vertex, the third column is the label of the source vertex and the fourth column is the label of the target vertex.
    ''' 
    def __init__(self, nodes, labels):
        self._nodes = nodes
        self._labels = labels
        self._path = os.path.dirname(__file__)
    
    # Function to save graph to CSV
    def _save_graph_to_csv(self, graph : nx, filename):
        filename = os.path.join(self._path, filename)

        with open(filename, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Number of vertices " + str(graph.number_of_nodes())])
            writer.writerow(['u', 'v', 'l1', 'l2'])  # Header
            for edge in graph.edges():
                u, v = edge
                l1 = graph.nodes[u]['label']
                l2 = graph.nodes[v]['label']
                writer.writerow([u, v, l1, l2])

    # Generate dense graph if p_dense is near 1 and sparse graph if p_dense is near 0. Density is near p_dense
    def generate_sparse_dense_graph(self, p_dense):
        kind_graph = 0 if p_dense < 0.5 else 1

        for v in self._nodes:
            g1 = nx.gnp_random_graph(v, p_dense, 42)
            g2 = nx.gnp_random_graph(v, p_dense, 42)

            for node in g1.nodes():
                label = random.randint(0, self._labels - 1)
                g1.nodes[node]["label"] = label
                g2.nodes[node]["label"] = label

            self._save_graph_to_csv(g1, f"graph_query_{v}_{kind_graph}.csv")
            self._save_graph_to_csv(g2, f"graph_target_{v}_{kind_graph}.csv")
    
    # Generate a complete graph
    def generate_complete_graph(self):
        for v in self._nodes:
            g1 = nx.complete_graph(v)
            g2 = nx.complete_graph(v)

            for node in g1.nodes():
                label = random.randint(0, self._labels - 1)
                g1.nodes[node]["label"] = label
                g2.nodes[node]["label"] = label

            self._save_graph_to_csv(g1, f"graph_query_{v}_2.csv")
            self._save_graph_to_csv(g2, f"graph_target_{v}_2.csv")
    
    # Generate a random graph
    def generate_random_graph(self):
        for v in self._nodes:
            g1 = nx.gnp_random_graph(v, 0.55, 42)
            g2 = nx.gnp_random_graph(v, 0.55, 42)

            for node in g1.nodes():
                label = random.randint(0, self._labels - 1)
                g1.nodes[node]["label"] = label
                g2.nodes[node]["label"] = label

            self._save_graph_to_csv(g1, f"graph_query_{v}_3.csv")
            self._save_graph_to_csv(g2, f"graph_target_{v}_3.csv")
    
if __name__ == "__main__":
    generator = GraphGenerator([3000, 5000, 8000, 10000], 10)
    generator.generate_random_graph()
    generator.generate_sparse_dense_graph(0.2)
    generator.generate_sparse_dense_graph(0.7)
    generator.generate_complete_graph()


        

        

    