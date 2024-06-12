#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define FILENAME_QUERY "../data/graph_query_example.csv"
#define FILENAME_TARGET "../data/graph_target_example.csv"
#define LABELS 10

typedef struct {
    int* matrix;
    int numVertices;
    int* nodesToLabel;  //manca labelToNodes che potrebbe essere implementato come una funzione a cui passo la label e amen (oppure linked list)
    int* labelsCardinalities;
    int* degrees;
} Graph;

typedef struct {
    int *mapping1;  // mapping from query to target
    int *mapping2;  // mapping from target to query
    int *T1;        // Ti contains uncovered neighbors of covered nodes from Gi, i.e. nodes that are not in the mapping, but are neighbors of nodes that are.
    int *T2;        
    int* T1_out;     //Ti_out contains all the nodes from Gi, that are neither in the mapping nor in Ti
    int* T2_out;
} State;

/***** GRAPH PROTOTYPES *****/
int* createAdjMatrix(int);
void initLabels(Graph*);
Graph* createGraph();
void addEdge(Graph*, int, int);
Graph* readGraph(char*);
void printGraph(Graph*);
void freeGraph(Graph*);
void setLabel(Graph*, int, int);
int* labelsCardinalities(int);
int* degrees(Graph*);

/***** STATE PROTOTYPES *****/
State* createState(Graph*, Graph*);
void freeState(State*);
void printState(State*, int);

/***** VF2++ PROTOTYPES *****/
void vf2pp(Graph*, Graph*, State*);

int main() {

    Graph* g1 = readGraph(FILENAME_QUERY);
    Graph* g2 = readGraph(FILENAME_TARGET);
    State* s = createState(g1, g2);
    
    printState(s, g1->numVertices);
    
    printf("Graph query:\n");
    printGraph(g1);
    printf("\nGraph target:\n");
    printGraph(g2);

    freeGraph(g1);
    freeGraph(g2);
    freeState(s);

    printf("fine");
    return 0;
}

/***** GRAPH FUNCTIONS *****/

int* createAdjMatrix(int numVertices) {
    int* matrix = (int*)malloc(numVertices * numVertices * sizeof(int));
    
    if (matrix == NULL) {
        printf("Error allocating memory in createAdjMatrix\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < numVertices; i++) {
        for (int j = 0; j < numVertices; j++) {
            matrix[i * numVertices + j] = 0;
        }
    }
    
    return matrix;
}

void initLabels(Graph* g) {
    g->nodesToLabel = (int*)malloc(g->numVertices * sizeof(int));
    g->labelsCardinalities = (int*)malloc(LABELS * sizeof(int));

    if (g->nodesToLabel == NULL || g->labelsCardinalities == NULL) {
        printf("Error allocating memory in initLabels\n");
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < g->numVertices; i++) {
        g->nodesToLabel[i] = -1;
    }

    for (int i = 0; i < LABELS; i++) {
        g->labelsCardinalities[i] = 0;
    }
}

void setLabel(Graph* g, int node, int label) {
    if (g->nodesToLabel[node] == -1) {
        g->nodesToLabel[node] = label;
        g->labelsCardinalities[label]++;
    }
}

void addEdge(Graph* g, int src, int target) {
    g->matrix[src * g->numVertices + target] = 1;
    g->matrix[target * g->numVertices + src] = 1;
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

    g->matrix = createAdjMatrix(g->numVertices);
    initLabels(g);


    while (fgets(line, sizeof(line), f)) {
        printf("%s", line);
        sscanf(line, "%d,%d,%d,%d", &src, &target, &srcLabel, &targetLabel);
        addEdge(g, src, target);
        setLabel(g, src, srcLabel);
        setLabel(g, target, targetLabel);
    }
    
    fclose(f);
    
    g->degrees = degrees(g);

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
}

void freeGraph(Graph* g) {
    free(g->matrix);
    free(g->nodesToLabel);
    free(g->labelsCardinalities);
    free(g->degrees);
    free(g);
    g = NULL;
}

int* degrees(Graph* g) {
    int* degrees = (int*)malloc(g->numVertices * sizeof(int));
    
    if (degrees == NULL) {
        printf("Error allocating memory in degrees\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < g->numVertices; i++) {
        for (int j = 0; j < g->numVertices; j++) {
            degrees[i] += g->matrix[i * g->numVertices + j];
        }
    }
    
    return degrees;
}

/***** STATE FUNCTIONS *****/
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
        s->T1[i] = 1;
        s->T1_out[i] = -1;

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
        printf("%d ", s->T1[i]);
    }

    printf("\nT2\n");
    for (int i = 0; i < numVertices; i++) {
        printf("%d ", s->T2[i]);
    }

    printf("\nT1_out\n");
    for (int i = 0; i < numVertices; i++) {
        printf("%d ", s->T1_out[i]);
    }

    printf("\nT2_out\n");
    for (int i = 0; i < numVertices; i++) {
        printf("%d ", s->T2_out[i]);
    }
}

/***** VF2++ FUNCTIONS *****/