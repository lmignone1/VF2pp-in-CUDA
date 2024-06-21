#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int* matrix;
    int numVertices;
    int* nodesToLabel; 
    int** labelToNodes;
    int* labelsCardinalities;
    int* degrees;
    int maxDegree;
} Graph;

typedef struct {
    int *mapping1;  
    int *mapping2;  
    int *T1;        
    int *T2;        
    int* T1_out;     
    int* T2_out;
} State;

bool set_intersection_count(int *set1, int size1, int *set2, int size2) {
    int count = 0;
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            if (set1[i] == set2[j]) {
                count++;
                break;
            }
        }
    }
    return count;
}

bool contains_same_labels(int *neighbors1, int *labels1, int size1, int *neighbors2, int *labels2, int size2) {
    int *unique_labels1 = malloc(size1 * sizeof(int));
    int *unique_labels2 = malloc(size2 * sizeof(int));
    int count1 = 0, count2 = 0;

    for (int i = 0; i < size1; i++) {
        bool found = false;
        for (int j = 0; j < count1; j++) {
            if (unique_labels1[j] == labels1[neighbors1[i]]) {
                found = true;
                break;
            }
        }
        if (!found) {
            unique_labels1[count1++] = labels1[neighbors1[i]];
        }
    }

    for (int i = 0; i < size2; i++) {
        bool found = false;
        for (int j = 0; j < count2; j++) {
            if (unique_labels2[j] == labels2[neighbors2[i]]) {
                found = true;
                break;
            }
        }
        if (!found) {
            unique_labels2[count2++] = labels2[neighbors2[i]];
        }
    }

    bool result = (count1 == count2);
    for (int i = 0; i < count1 && result; i++) {
        bool found = false;
        for (int j = 0; j < count2; j++) {
            if (unique_labels1[i] == unique_labels2[j]) {
                found = true;
                break;
            }
        }
        if (!found) {
            result = false;
        }
    }

    free(unique_labels1);
    free(unique_labels2);
    return result;
}

bool _cut_PT(int u, int v, Graph *G1, Graph *G2, int *G1_labels, int *G2_labels, State *state) {
    int *neighbors_u = &G1->matrix[u * G1->numVertices];
    int *neighbors_v = &G2->matrix[v * G2->numVertices];

    // Check if neighbors of u and v have the same labels
    if (!contains_same_labels(neighbors_u, G1_labels, G1->numVertices, neighbors_v, G2_labels, G2->numVertices)) {
        return true;
    }

    // Iterate through unique labels and perform the intersection checks
    for (int i = 0; i < G1->numVertices; i++) {
        if (neighbors_u[i]) {
            int label = G1_labels[i];
            int *G1_nbh = &G1->labelToNodes[label][0];
            int *G2_nbh = &G2->labelToNodes[label][0];
            int G1_nbh_size = G1->labelsCardinalities[label];
            int G2_nbh_size = G2->labelsCardinalities[label];

            if (set_intersection_count(state->T1, G1_nbh_size, G1_nbh, G1_nbh_size) != set_intersection_count(state->T2, G2_nbh_size, G2_nbh, G2_nbh_size)) {
                return true;
            }
            if (set_intersection_count(state->T1_out, G1_nbh_size, G1_nbh, G1_nbh_size) != set_intersection_count(state->T2_out, G2_nbh_size, G2_nbh, G2_nbh_size)) {
                return true;
            }
        }
    }

    return false;
}


// --------------------------------------

#include <stdio.h>
#include <stdlib.h>

// Assume these structures are defined as provided in the question
typedef struct {
    int* matrix;
    int numVertices;
    int* nodesToLabel; 
    int** labelToNodes; // alternativa: funzione
    int* labelsCardinalities;
    int* degrees;
    int maxDegree;
} Graph;

typedef struct {
    int *mapping1;  // mapping from query to target
    int *mapping2;  // mapping from target to query
    int *T1;        // Ti contains uncovered neighbors of covered nodes from Gi, i.e. nodes that are not in the mapping, but are neighbors of nodes that are.
    int *T2;        
    int* T1_out;     //Ti_out contains all the nodes from Gi, that are neither in the mapping nor in Ti
    int* T2_out;
} State;

// Function to update Ti and Ti_out
void update_Tinout(int new_node1, int new_node2, Graph* graph1, Graph* graph2, State* state) {
    int numVertices1 = graph1->numVertices;
    int numVertices2 = graph2->numVertices;

    // Update T1 and T1_out for graph1
    for (int i = 0; i < numVertices1; i++) {
        if (graph1->matrix[new_node1 * numVertices1 + i] && state->mapping1[i] == -1) {
            state->T1[i] = 1;
            state->T1_out[i] = 0;
        }
    }
    state->T1[new_node1] = 0;
    state->T1_out[new_node1] = 0;

    // Update T2 and T2_out for graph2
    for (int i = 0; i < numVertices2; i++) {
        if (graph2->matrix[new_node2 * numVertices2 + i] && state->mapping2[i] == -1) {
            state->T2[i] = 1;
            state->T2_out[i] = 0;
        }
    }
    state->T2[new_node2] = 0;
    state->T2_out[new_node2] = 0;
}

int main() {
    // Example usage
    Graph graph1, graph2;
    State state;
    
    // Initialize graphs and state...
    
    int new_node1 = 0; // Example node from graph1
    int new_node2 = 1; // Example node from graph2
    
    // Call the function to update T1, T2, T1_out, T2_out
    update_Tinout(new_node1, new_node2, &graph1, &graph2, &state);
    
    return 0;
}







#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int* matrix;
    int numVertices;
    int* nodesToLabel; 
    int** labelToNodes;
    int* labelsCardinalities;
    int* degrees;
    int maxDegree;
} Graph;

typedef struct {
    int *mapping1;
    int *mapping2;
    int *T1;
    int *T2;
    int* T1_out;
    int* T2_out;
} State;

State* createState(Graph* g1, Graph* g2) {
    State* s = (State*)malloc(sizeof(State));
    
    s->mapping1 = (int*)malloc(g1->numVertices * sizeof(int));
    s->mapping2 = (int*)malloc(g2->numVertices * sizeof(int));
    s->T1 = (int*)malloc(g1->numVertices * sizeof(int));
    s->T2 = (int*)malloc(g2->numVertices * sizeof(int));
    s->T1_out = (int*)malloc(g1->numVertices * sizeof(int));
    s->T2_out = (int*)malloc(g2->numVertices * sizeof(int));

    if (s == NULL || s->mapping1 == NULL || s->mapping2 == NULL || s->T1 == NULL || s->T2 == NULL || s->T1_out == NULL || s->T2_out == NULL) {
        printf("Error allocating memory in createState\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < g1->numVertices; i++) {
        s->mapping1[i] = -1;
        s->T1[i] = -1;
        s->T1_out[i] = 1;

        s->mapping2[i] = -1;
        s->T2[i] = -1;
        s->T2_out[i] = 1;
    }

    return s;
}

// -------------------------------------------

void restore_Tinout(int popped_node1, int popped_node2, Graph* g1, Graph* g2, State* state) {
    int* mapping = state->mapping1;
    int* reverse_mapping = state->mapping2;
    int* T1 = state->T1;
    int* T1_out = state->T1_out;
    int* T2 = state->T2;
    int* T2_out = state->T2_out;

    int is_added = 0;

    for (int i = 0; i < g1->numVertices; i++) {
        if (g1->matrix[popped_node1 * g1->numVertices + i] == 1) {
            int nbr = i;
            if (mapping[nbr] != -1) {
                T1[popped_node1] = 1; 
                is_added = 1;
            } else {
                int has_covered_neighbor = 0;
                for (int j = 0; j < g1->numVertices; j++) {
                    if (g1->matrix[nbr * g1->numVertices + j] == 1 && mapping[j] != -1) {
                        has_covered_neighbor = 1;
                        break;
                    }
                }
                if (!has_covered_neighbor) {
                    T1[nbr] = -1;
                    T1_out[nbr] = 1;
                }
            }
        }
    }

    if (!is_added) {
        T1_out[popped_node1] = 1;
    }

    is_added = 0;

    for (int i = 0; i < g2->numVertices; i++) {
        if (g2->matrix[popped_node2 * g2->numVertices + i] == 1) {
            int nbr = i;
            if (reverse_mapping[nbr] != -1) {
                T2[popped_node2] = 1;
                is_added = 1;
            } else {
                int has_covered_neighbor = 0;
                for (int j = 0; j < g2->numVertices; j++) {
                    if (g2->matrix[nbr * g2->numVertices + j] == 1 && reverse_mapping[j] != -1) {
                        has_covered_neighbor = 1;
                        break;
                    }
                }
                if (!has_covered_neighbor) {
                    T2[nbr] = -1;
                    T2_out[nbr] = 1;
                }
            }
        }
    }

    if (!is_added) {
        T2_out[popped_node2] = 1;
    }
}

int main() {
    // Example usage
    Graph g1, g2;
    // Initialize g1 and g2 properly here
    State* state = createState(&g1, &g2);

    int popped_node1 = 0; // Example node
    int popped_node2 = 1; // Example node

    restore_Tinout(popped_node1, popped_node2, &g1, &g2, state);

    // Free allocated memory
    free(state->mapping1);
    free(state->mapping2);
    free(state->T1);
    free(state->T2);
    free(state->T1_out);
    free(state->T2_out);
    free(state);

    return 0;
}







