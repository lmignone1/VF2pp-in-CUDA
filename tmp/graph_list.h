#ifndef GRAPH_LIST_H
#define GRAPH_LIST_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct Node {
    int vertex;
    struct Node* next;
} Node;

typedef struct AdjList {
    Node* head;
    int label;
} AdjList;

typedef struct Graph {
    int numVertices;
    AdjList* nodesList;
} Graph;

Node* createNode(int);
Graph* createGraph(int);
void addEdge(Graph*, int, int, int, int);
void printGraph(Graph*);
void freeGraph(Graph*);
Graph* readGraph(char*);

#ifdef __cplusplus
}
#endif

#endif
