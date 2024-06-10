#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct Node {
    int vertex;
    int label;  // forse non serve
    struct Node* next;
} Node;

typedef struct AdjList {
    Node* head;
} AdjList;

typedef struct Graph {
    int numVertices;
    AdjList* nodesList;
} Graph;

Node* createNode(int, int);
Graph* createGraph(int);
void addEdge(Graph*, int, int, int, int);
void printGraph(Graph*);
void freeGraph(Graph*);
Graph* readGraph(char*);

#endif 
