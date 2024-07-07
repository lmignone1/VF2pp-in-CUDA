#ifndef GRAPH_H
#define GRAPH_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#define LABELS 10

typedef struct {
    int* matrix;
    int numVertices;
    int* nodesToLabel; 
    int** labelToNodes; 
    int* labelsCardinalities;
    int* degrees;
} Graph;

void initGraph(Graph*);
Graph* createGraph();
void addEdge(Graph*, int, int);
Graph* readGraph(char*);
void printGraph(Graph*);
void freeGraph(Graph*);
void setLabel(Graph*, int, int);

#ifdef __cplusplus
}
#endif

#endif
