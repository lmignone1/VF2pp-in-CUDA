#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define FILENAME_QUERY "../data/graph_query.csv"
#define FILENAME_TARGET "../data/graph_target.csv"
#define LABELS 10
#define INF 99999

/***** STRUCTS *****/
typedef struct {
    int* matrix;
    int numVertices;
    int* nodesToLabel; 
    int** labelToNodes; // alternativa: funzione
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

typedef struct {
    int* data;
    int head;
    int tail;
    int capacity;
} Queue;

/***** GRAPH PROTOTYPES *****/
int* createAdjMatrix(int);
void initGraph(Graph*);
Graph* createGraph();
void addEdge(Graph*, int, int);
Graph* readGraph(char*);
void printGraph(Graph*);
void freeGraph(Graph*);
void setLabel(Graph*, int, int);
int* labelsCardinalities(int);

/***** STATE PROTOTYPES *****/
State* createState(Graph*, Graph*);
void freeState(State*);
void printState(State*, int);

/***** VF2++ PROTOTYPES *****/
void vf2pp(Graph*, Graph*, State*);
bool checkGraphProperties(Graph*, Graph*);
bool checkSequenceDegree(int*, int*, int);
int compare(const void*, const void*);
int* ordering(Graph*, Graph*);
int* copyArray(int*, int);
int* bfs(Graph*, int, int*);
int findLevelNodes(Graph*, int*, int*, int);
void removeElementArray(int*, int*, int);
void processDepth(int*, int*, Graph*, int*, int*, int*, int*, int*);

/***** Queue FUNCTIONS *****/
Queue* createQueue(int);
bool isQueueEmpty(Queue*);
void resizeQueue(Queue*);
void enqueue(Queue*, int);
int dequeue(Queue*);
void freeQueue(Queue*);

int main() {
    printf("%s\n", FILENAME_QUERY);
    Graph* g1 = readGraph(FILENAME_QUERY);
    printf("%s\n", FILENAME_TARGET);
    Graph* g2 = readGraph(FILENAME_TARGET);
    State* s = createState(g1, g2);
    vf2pp(g1, g2, s);
    freeGraph(g1);
    freeGraph(g2);
    freeState(s);

    return EXIT_SUCCESS;
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

void initGraph(Graph* g) {
    g->nodesToLabel = (int*)malloc(g->numVertices * sizeof(int));
    g->labelsCardinalities = (int*)malloc(LABELS * sizeof(int));
    g->labelToNodes = (int**)malloc(LABELS * sizeof(int*));
    g->degrees = (int*)malloc(g->numVertices * sizeof(int));

    if (g->nodesToLabel == NULL || g->labelsCardinalities == NULL || g->labelToNodes == NULL || g->degrees == NULL) {
        printf("Error allocating memory in initGraph\n");
        exit(EXIT_FAILURE);
    }
    
    for (int vertex = 0; vertex < g->numVertices; vertex++) {
        g->nodesToLabel[vertex] = -1;
        g->degrees[vertex] = 0;
    }

    for (int label = 0; label < LABELS; label++) {
        g->labelsCardinalities[label] = 0;
        g->labelToNodes[label] = malloc(g->numVertices * sizeof(int));
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

    g->matrix = createAdjMatrix(g->numVertices);
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
    for(int i = 0; i < g->numVertices; i++) {
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

// int* degrees(Graph* g) {
//     int* degrees = (int*)malloc(g->numVertices * sizeof(int));
    
//     if (degrees == NULL) {
//         printf("Error allocating memory in degrees\n");
//         exit(EXIT_FAILURE);
//     }

//     for (int i = 0; i < g->numVertices; i++) {
//         for (int j = 0; j < g->numVertices; j++) {
//             degrees[i] += g->matrix[i * g->numVertices + j];
//         }
//     }
    
//     return degrees;
// }

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
bool checkGraphProperties(Graph* g1, Graph* g2) {
    if (g1->numVertices != g2->numVertices || g1->numVertices == 0 || g2->numVertices == 0) {
        return false;
    }

    if (!checkSequenceDegree(g1->degrees, g2->degrees, g1->numVertices)) {
        return false;
    }

    for(int i = 0; i < LABELS; i++) {
        if (g1->labelsCardinalities[i] != g2->labelsCardinalities[i]) {
            return false;
        }
    }

    return true;
}

int* copyArray(int* arr, int size) {
    int* copy = (int*)malloc(size * sizeof(int));
    
    if (copy == NULL) {
        printf("Error allocating memory in copyArray\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < size; i++) {
        copy[i] = arr[i];
    }

    return copy;
}

bool checkSequenceDegree(int* degree1, int* degree2, int size) {
    int* tmp1 = copyArray(degree1, size);
    int* tmp2 = copyArray(degree2, size);

    qsort(tmp1, size, sizeof(int), compare);
    qsort(tmp2, size, sizeof(int), compare);
   
    bool ret = true;
    for (int i = 0; i < size; i++) {
        if (tmp1[i] != tmp2[i]) {
            ret = false;
            break;
        }
    }

    free(tmp1);
    free(tmp2);
    return ret;
}

int compare(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

int* ordering(Graph* g1, Graph* g2) {
    int* order = (int*)malloc(g1->numVertices * sizeof(int));   // order of the nodes of g1
    int* labelRarity = copyArray(g1->labelsCardinalities, LABELS);  // FORSE POSSIAMO EVITARE LA COPIA E USARE DIRETTAMENTE g1->labelsCardinalities
    int* connectivityG1 = (int*)malloc(g1->numVertices * sizeof(int));   // number of neighbors already ordered
    int* V1Unordered = (int*)malloc(g1->numVertices * sizeof(int));    // V1Unordered[i] = -1 if node has not been ordered yet else 1 

    if (order == NULL || labelRarity == NULL || connectivityG1 == NULL || V1Unordered == NULL) {
        printf("Error allocating memory in ordering\n");
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < g1->numVertices; i++) {
        connectivityG1[i] = 0;
        V1Unordered[i] = -1;
    }

    int order_index = 0;
    while (order_index < g1->numVertices) {
        int maxRarity = INF;
        int maxNode = -1;
        for (int vertex = 0; vertex < g1->numVertices; vertex++) {
            if (V1Unordered[vertex] == -1) {
                int rarity = labelRarity[g1->nodesToLabel[vertex]];
                if (rarity < maxRarity) {
                    maxRarity = rarity;
                    maxNode = vertex;
                } else if (rarity == maxRarity) {
                    if (g1->degrees[vertex] > g1->degrees[maxNode]) {
                        maxNode = vertex;
                    }
                }
            }
        }
        
        int maxDepth = 0;
        int* levels = bfs(g1, maxNode, &maxDepth);
        int* levelNodes = (int*)malloc((g1->numVertices) * sizeof(int));

        for (int depth = 0; depth <= maxDepth; depth++) {
            int levelSize = findLevelNodes(g1, levels, levelNodes, depth);
            processDepth(order, &order_index, g1, connectivityG1, labelRarity, V1Unordered, levelNodes, &levelSize);
        }
        free(levels);
        free(levelNodes);
    }           
    free(labelRarity);
    free(connectivityG1);
    free(V1Unordered);

    for(int i = 0; i < g1->numVertices; i++) {
        printf("%d ", order[i]);
    }
    return order;
}

int findLevelNodes(Graph* g, int* levels, int* levelNodes, int depth) {
    int i = 0;
    for(int vertex = 0; vertex < g->numVertices; vertex++) {
        if (levels[vertex] == depth) {
            levelNodes[i] = vertex;
            i++;
        }
    }
    return i;
}

int* bfs(Graph* g, int root, int* maxDepth) {
    int* levels = (int*)malloc(g->numVertices * sizeof(int));

    if (levels == NULL) {
        printf("Error allocating memory in bfs\n");
        exit(EXIT_FAILURE);
    }

    for(int vertex = 0; vertex < g->numVertices; vertex++) {
        levels[vertex] = -1;
    }

    Queue* queue = createQueue(g->numVertices/2);
    enqueue(queue, root);
    levels[root] = 0;

    *maxDepth = 0;
    while (!isQueueEmpty(queue)) {
        int node = dequeue(queue);
        for (int adjVertex = 0; adjVertex < g->numVertices; adjVertex++) {
            if (g->matrix[node * g->numVertices + adjVertex] == 1 && levels[adjVertex] == -1) {
                enqueue(queue, adjVertex);
                int depth = levels[node] + 1;
                levels[adjVertex] = depth;
                if (depth > *maxDepth) {
                    *maxDepth = depth;
                }
            }
        }
    }
    freeQueue(queue);
    return levels;
}

void removeElementArray(int* arr, int* size, int element) {
    for (int i = 0; i < *(size); i++) {
        if (arr[i] == element) {
            for (int j = i; j < *(size) - 1; j++) {
                arr[j] = arr[j + 1];
            }
            break;
        }
    }
    *size = *size - 1;
}

void processDepth(int* order, int* order_index, Graph* g, int* connectivityG1, int* labelRarity, int* V1Unordered, int* levelNodes, int* levelSize) {
    while (*levelSize > 0) {
        int maxConnectivity = -INF;  
        int maxDegree = -INF;
        int maxRarity = INF;
        int nextNode = -1;

        for(int i = 0; i < *levelSize; i++) {
            int vertex = levelNodes[i];
            int conn = connectivityG1[vertex];
            if (conn > maxConnectivity) {
                maxConnectivity = conn;
                maxDegree = -INF; 
            }
            if (conn == maxConnectivity) {
                int degree = g->degrees[vertex];
                if (degree > maxDegree) {
                    maxDegree = degree;
                    maxRarity = INF;
                }
                if (degree == maxDegree) {
                    int rarity = labelRarity[g->nodesToLabel[vertex]];
                    if(rarity < maxRarity) {
                        maxRarity = rarity;
                        nextNode = vertex;
                    }
                }
            }
        }
        
        order[*order_index] = nextNode;
        *order_index = *order_index + 1;    
        
        for (int adjVertex = 0; adjVertex < g->numVertices; adjVertex++) {
            if (g->matrix[nextNode * g->numVertices + adjVertex] == 1) {
                connectivityG1[adjVertex]++;
            }
        }
        
        removeElementArray(levelNodes, levelSize, nextNode);
        labelRarity[g->nodesToLabel[nextNode]]--;
        V1Unordered[nextNode] = 1;
    }
}

int* match(int node, Graph* g2, State* state) {
    
}

void vf2pp(Graph* g1, Graph* g2, State* state) {
    if (!checkGraphProperties(g1, g2)) {
        printf("Graphs are not isomorphic\n");
        return;
    }

    int* order = ordering(g1, g2);
    
    free(order);

    printf("Graphs are isomorphic\n");
}

/***** QUEUE FUNCTIONS *****/
Queue* createQueue(int capacity) {
    Queue* q = (Queue*)malloc(sizeof(Queue));
    q->data = (int*)malloc(capacity * sizeof(int));
    q->head = -1;
    q->tail = -1;
    q->capacity = capacity;

    if (q == NULL || q->data == NULL) {
        printf("Error allocating memory in createQueue\n");
        exit(EXIT_FAILURE);
    }

    return q;
}

bool isQueueEmpty(Queue* q) {
    return q->head == -1;
}

void resizeQueue(Queue* q) {
    q->capacity *= 2;
    q->data = (int*)realloc(q->data, q->capacity * sizeof(int));

    if (q->data == NULL) {
        printf("Error reallocating memory in resizeQueue\n");
        exit(EXIT_FAILURE);
    }
}

void enqueue(Queue* q, int item) {
    if (q->tail == q->capacity - 1) {
        resizeQueue(q);
    }
    if (isQueueEmpty(q)) {
        q->head = 0;
    }
    q->tail++;
    q->data[q->tail] = item;
}

int dequeue(Queue* q) {
    if (isQueueEmpty(q)) {
        return -1;
    } else {
        int item = q->data[q->head];
        q->head++;
        if (q->head > q->tail) {
            q->head = -1;
            q->tail = -1;
        }
        return item;
    }
}

void freeQueue(Queue* q) {
    free(q->data);
    free(q);
    q = NULL;
}


