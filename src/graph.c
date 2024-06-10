#include "graph.h"

Node* createNode(int vertex, int label) {
    Node* node = (Node*)malloc(sizeof(Node));
    node->vertex = vertex;
    node->label = label;
    node->next = NULL;
    return node;
}

Graph* createGraph(int numVertices) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->numVertices = numVertices;
    graph->nodesList = (AdjList*)malloc(numVertices * sizeof(AdjList));
    
    for (int i = 0; i < numVertices; ++i) {
        graph->nodesList[i].head = NULL;
    }
    
    return graph;
}

void addEdge(Graph* graph, int src, int target, int srcLabel, int targetLabel) {
    Node* node = createNode(target, targetLabel);
    node->next = graph->nodesList[src].head;
    graph->nodesList[src].head = node;
    
    node = createNode(src, srcLabel);
    node->next = graph->nodesList[target].head;
    graph->nodesList[target].head = node;
}

void printGraph(Graph* graph) {
    for (int v = 0; v < graph->numVertices; ++v) {
        Node* currentNode = graph->nodesList[v].head;
        printf("\n // Adjency list of vertex %d\n head", v);
        while (currentNode) {
            printf(" -> %d", currentNode->vertex);
            currentNode = currentNode->next;
        }
        printf("\n");
    }
}

void freeGraph(Graph* graph) {
    for (int i = 0; i < graph->numVertices; ++i) {
        Node* currentNode = graph->nodesList[i].head;
        while (currentNode) {
            Node* tmp = currentNode;
            currentNode = currentNode->next;
            free(tmp);
        }
    }
    free(graph->nodesList);
    free(graph);
}

Graph* readGraph(char* path) {
    int src, target, srcLabel, targetLabel, numVertices;

    FILE* f = fopen(path, "r");
    if (f == NULL) {
        printf("Error opening file\n");
        exit(EXIT_FAILURE);
    }

    char line[128];
    fgets(line, sizeof(line), f);
    sscanf(line, "%*s%*s%*s%d", &numVertices);
    fgets(line, sizeof(line), f); // skip the header

    Graph* graph = createGraph(numVertices);

    while (fgets(line, sizeof(line), f)) {
        sscanf(line, "%d,%d,%d,%d", &src, &target, &srcLabel, &targetLabel);
        addEdge(graph, src, target, srcLabel, targetLabel);
    }

    fclose(f);
    return graph;
}



