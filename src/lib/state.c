#include "state.h"

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

void freeState(State* s) {
    free(s->mapping1);
    free(s->mapping2);
    free(s->T1);
    free(s->T2);
    free(s->T1_out);
    free(s->T2_out);
    free(s);
    s = NULL;
}

void printState(State* s, int numVertices) {
    printf("Mapping 1\n");
    for (int i = 0; i < numVertices; i++) {
        printf("%d ", s->mapping1[i]);
    }

    printf("\nMapping 2\n");
    for (int i = 0; i < numVertices; i++) {
        printf("%d ", s->mapping2[i]);
    }

    printf("\nT1\n");
    for (int i = 0; i < numVertices; i++) {
        if(s->T1[i] != -1) {
            printf("%d ", i);
        }
        // printf("%d ", s->T1[i]);
    }

    printf("\nT2\n");
    for (int i = 0; i < numVertices; i++) {
        if(s->T2[i] != -1) {
            printf("%d ", i);
        }
        // printf("%d ", s->T2[i]);
    }

    printf("\nT1_out\n");
    for (int i = 0; i < numVertices; i++) {
        if(s->T1_out[i] != -1) {
            printf("%d ", i);
        }
        // printf("%d ", s->T1_out[i]);
    }

    printf("\nT2_out\n");
    for (int i = 0; i < numVertices; i++) {
        if(s->T2_out[i] != -1) {
            printf("%d ", i);
        }
        // printf("%d ", s->T2_out[i]);
    }
}

bool isMappingFull(Graph* g, State* state) {
    for(int i = 0; i < g->numVertices; i++) {
        if(state->mapping1[i] == -1) {
            return false;
        }
    }
    return true;
}

void updateState(Graph* g1, Graph* g2, State* state, int node, int candidate) {
    for(int adjVertex = 0; adjVertex < g1->numVertices; adjVertex++) {
        
        if(g1->matrix[node * g1->numVertices + adjVertex] == 1 && state->mapping1[adjVertex] == -1) {
            state->T1[adjVertex] = 1;
            state->T1_out[adjVertex] = -1;
        }

        if(g2->matrix[candidate * g2->numVertices + adjVertex] == 1 && state->mapping2[adjVertex] == -1) {
            state->T2[adjVertex] = 1;
            state->T2_out[adjVertex] = -1;
        }
    }

    state->T1[node] = -1;
    state->T1_out[node] = -1;
    state->T2[candidate] = -1;
    state->T2_out[candidate] = -1;

    // printf("\n");
    // printf("Update state\n");
    // printState(state, g1->numVertices);
    // printf("\n");
}

void restoreState(Graph* g1, Graph* g2, State* state, int node, int candidate) {
    bool isAdded = false;
    for(int adjVertex = 0; adjVertex < g1->numVertices; adjVertex++) {
        if(g1->matrix[node * g1->numVertices + adjVertex] == 1) {
            
            if(state->mapping1[adjVertex] != -1) {
                state->T1[node] = 1;
                isAdded = true;     
            }
            else {
                bool hasCoveredNeighbor = false;
                for(int adjVertex2 = 0; adjVertex2 < g1->numVertices; adjVertex2++) {
                    if(g1->matrix[adjVertex * g1->numVertices + adjVertex2] == 1 && state->mapping1[adjVertex2] != -1) {
                        hasCoveredNeighbor = true;
                        break;
                    }
                }

                if(!hasCoveredNeighbor) {
                    state->T1[adjVertex] = -1;
                    state->T1_out[adjVertex] = 1;
                }
            }
        }
    }

    if(!isAdded) {
        state->T1_out[node] = 1;
    }

    isAdded = false;    
    for(int adjVertex = 0; adjVertex < g2->numVertices; adjVertex++) {
        if(g2->matrix[candidate * g2->numVertices + adjVertex] == 1) {
            
            if(state->mapping2[adjVertex] != -1) {
                state->T2[candidate] = 1;
                isAdded = true;
            }
            else {
                bool hasCoveredNeighbor = false;
                for(int adjVertex2 = 0; adjVertex2 < g2->numVertices; adjVertex2++) {
                    if(g2->matrix[adjVertex * g2->numVertices + adjVertex2] == 1 && state->mapping2[adjVertex2] != -1) {
                        hasCoveredNeighbor = true;
                        break;
                    }
                }

                if(!hasCoveredNeighbor) {
                    state->T2[adjVertex] = -1;
                    state->T2_out[adjVertex] = 1;
                }
            }
        }
    }
    
    if(!isAdded) {
        state->T2_out[candidate] = 1;
    }

    // printf("\n");
    // printf("Restore state\n");
    // printState(state, g1->numVertices);
    // printf("\n");

}