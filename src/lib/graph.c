/*
* MIGNONE LORENZO 0622701866 L.MIGNONE@STUDENTI.UNISA.IT
* Course: High Performance Computing 2022/2023
* Lecturer: Francesco Moscato	fmoscato@unisa.it
* 
* Copyright (C) 2024 - All Rights Reserved
* 
* This file is part of VF2pp-in-CUDA.
* 
* VF2pp-in-CUDA is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* VF2pp-in-CUDA is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with VF2pp-in-CUDA.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "graph.h"

void initGraph(Graph* g) {
    g->matrix = (int*)malloc(g->numVertices * g->numVertices * sizeof(int));
    g->nodesToLabel = (int*)malloc(g->numVertices * sizeof(int));
    g->labelsCardinalities = (int*)malloc(LABELS * sizeof(int));
    g->labelToNodes = (int**)malloc(LABELS * sizeof(int*));
    g->degrees = (int*)malloc(g->numVertices * sizeof(int));

    if (g->nodesToLabel == NULL || g->labelsCardinalities == NULL || g->labelToNodes == NULL || g->degrees == NULL || g->matrix == NULL) {
        printf("Error allocating memory in initGraph\n");
        exit(EXIT_FAILURE);
    }
    
    for (int vertex = 0; vertex < g->numVertices; vertex++) {
        g->nodesToLabel[vertex] = -1;
        g->degrees[vertex] = 0;

        for (int adjVertex = 0; adjVertex < g->numVertices; adjVertex++) {
            g->matrix[vertex * g->numVertices + adjVertex] = 0;
        }
    }

    for (int label = 0; label < LABELS; label++) {
        g->labelsCardinalities[label] = 0;
        g->labelToNodes[label] = (int*)malloc(g->numVertices * sizeof(int));
    }
}

void setLabel(Graph* g, int node, int label) {
    if (g->nodesToLabel[node] == -1) {
        g->nodesToLabel[node] = label;
        g->labelsCardinalities[label]++;
        g->labelToNodes[label][g->labelsCardinalities[label] - 1] = node;
    }
}

void addEdge(Graph* g, int src, int target) {
    g->matrix[src * g->numVertices + target] = 1;
    g->matrix[target * g->numVertices + src] = 1;
    g->degrees[src]++;
    g->degrees[target]++;
}

Graph* createGraph() {
    Graph* g = (Graph*)malloc(sizeof(Graph));

    if (g == NULL) {
        printf("Error allocating memory in createGraph\n");
        exit(EXIT_FAILURE);
    }

    g->matrix = NULL;
    g->numVertices = 0;
    g->nodesToLabel = NULL;
    g->labelsCardinalities = NULL;
    g->degrees = NULL;
    g->labelToNodes = NULL;
    return g;
}

Graph* readGraph(char* path) {
    int src, target, srcLabel, targetLabel;
    Graph* g = createGraph();
   
    FILE* f = fopen(path, "r");
    if (f == NULL) {
        printf("Error opening file\n");
        exit(EXIT_FAILURE);
    }

    char line[128];
    fgets(line, sizeof(line), f);
    sscanf(line, "%*s%*s%*s%d", &g->numVertices);
    fgets(line, sizeof(line), f); // skip the header

    initGraph(g);

    while (fgets(line, sizeof(line), f)) {
        sscanf(line, "%d,%d,%d,%d", &src, &target, &srcLabel, &targetLabel);
        addEdge(g, src, target);
        setLabel(g, src, srcLabel);
        setLabel(g, target, targetLabel);
    }
    
    fclose(f);

    for(int label = 0; label < LABELS; label++) {
        g->labelToNodes[label] = (int*)realloc(g->labelToNodes[label], g->labelsCardinalities[label] * sizeof(int));
    }
    
    return g;
}

void printGraph(Graph* g) {
    for (int i = 0; i < g->numVertices; i++) {
        for (int j = 0; j < g->numVertices; j++) {
            printf("%d ", g->matrix[i * g->numVertices + j]);
        }
        printf("\tVertex %d, label %d, degree %d\n", i, g->nodesToLabel[i], g->degrees[i]);
    }

    printf("\nCardinalities\n");
    for (int i = 0; i < LABELS; i++) {
        printf("Label %d: %d\n", i, g->labelsCardinalities[i]);
    }

    for(int i = 0; i < LABELS; i++) {
       printf("\nLabel %d\n", i);
       for(int j = 0; j < g->labelsCardinalities[i]; j++) {
           printf("%d ", g->labelToNodes[i][j]);
       }
    }
}

void freeGraph(Graph* g) {
    for(int i = 0; i < LABELS; i++) {
        free(g->labelToNodes[i]);
    }
    free(g->labelToNodes);
    free(g->matrix);
    free(g->nodesToLabel);
    free(g->labelsCardinalities);
    free(g->degrees);
    free(g);
    g = NULL;
}