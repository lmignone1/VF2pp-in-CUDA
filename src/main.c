#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define FILENAME_QUERY "../data/graph_query_example.csv"
#define FILENAME_TARGET "../data/graph_target_example.csv"
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
    int maxDegree;
} Graph;

typedef struct {
    int *mapping1;  // mapping from query to target
    int *mapping2;  // mapping from target to query
    int *T1;        // Ti contains uncovered neighbors of covered nodes from Gi, i.e. nodes that are not in the mapping, but are neighbors of nodes that are.
    int *T2;        
    int* T1_out;     //Ti_out contains all the nodes from Gi, that are neither in the mapping nor in Ti. Cioe nodi che non sono in mapping e non sono vicini di nodi coperti
    int* T2_out;
} State;

typedef struct {
    int* data;
    int head;
    int tail;
    int capacity;
} Queue;

typedef struct {
    int vertex;
    int* candidates;
    int sizeCandidates;
    int candidateIndex;
} Info;

typedef struct StackNode {
    Info* info;
    struct StackNode* next;
} StackNode;

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
bool isMappingFull(Graph*, State*);
void updateState(Graph*, Graph*, State*, int, int);
void restoreState(Graph*, Graph*, State*, int, int);

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
int* findCoveredNeighbors(Graph*, State*, int, int*);
int* findCandidates(Graph*, Graph*, State*, int, int*);
int intersectionCount(int*, int, int*);
int* findNodesOfLabel(int*, Graph*, int, int, int*);
bool cutISO(Graph*, Graph*, State*, int, int);
int* findNeighbors(Graph*, int, int*);
bool checkLabels(Graph*, Graph*, int*, int*, int, int, int*);

/***** QUEUE PROTOTYPES *****/
Queue* createQueue(int);
bool isQueueEmpty(Queue*);
void resizeQueue(Queue*);
void enqueue(Queue*, int);
int dequeue(Queue*);
void freeQueue(Queue*);

/***** STACK PROTOTYPES *****/
Info* createInfo(int* candidates, int sizeCandidates, int vertex);
StackNode* createStackNode(Info*);
void push(StackNode**, int*, int, int);
Info* pop(StackNode**);
bool isStackEmpty(StackNode*);
void freeStack(StackNode*);
void printStack(StackNode*);
void printInfo(Info*);
void freeInfo(Info*);
StackNode* createStack();
Info* peek(StackNode**);

int main() {
    printf("%s\n", FILENAME_QUERY);
    Graph* g1 = readGraph(FILENAME_QUERY);
    printf("%s\n", FILENAME_TARGET);
    Graph* g2 = readGraph(FILENAME_TARGET);
    State* s = createState(g1, g2);
    
    vf2pp(g1, g2, s);
    
    printf("\nMapping\n");
    for(int i = 0; i < g1->numVertices; i++) {
        printf("%d -> %d\n", i, s->mapping1[i]);
    }

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
    int degree = g->degrees[src] > g->degrees[target] ? g->degrees[src] : g->degrees[target];   // si possono rimuovere forse
    g->maxDegree = degree > g->maxDegree ? degree : g->maxDegree;   // si possono rimuovere forse
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
    g->maxDegree = 0;
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
}

void restoreState(Graph* g1, Graph* g2, State* state, int node, int candidate) {
    bool isAdded = false;
    for(int adjVertex = 0; adjVertex < g1->numVertices; adjVertex++) {
        if(g1->matrix[node * g1->numVertices + adjVertex] == 1) {
            
            if(state->mapping1[adjVertex] == 1) {
                state->T1[node] = 1;
                isAdded = true;     // verificare se l operazione legata a questo falg possa essere spostata qui
            }
            else {
                bool hasCoveredNeighbor = false;
                for(int adjVertex2 = 0; adjVertex2 < g1->numVertices; adjVertex2++) {
                    if(g1->matrix[adjVertex * g1->numVertices + adjVertex2] == 1 && state->mapping1[adjVertex2] == 1) {
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
            
            if(state->mapping2[adjVertex] == 1) {
                state->T2[candidate] = 1;
                isAdded = true;
            }
            else {
                bool hasCoveredNeighbor = false;
                for(int adjVertex2 = 0; adjVertex2 < g2->numVertices; adjVertex2++) {
                    if(g2->matrix[adjVertex * g2->numVertices + adjVertex2] == 1 && state->mapping2[adjVertex2] == 1) {
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

    // for(int i = 0; i < g1->numVertices; i++) {
    //     printf("%d ", order[i]);
    // }
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

int* findCandidates(Graph* g1, Graph* g2, State* state, int node, int* sizeCandidates) {  
    // void findCandidates(int u, GraphParams *graph_params, StateParams *state_params, int *candidates, int *candidate_count) {
    // Graph *G1 = graph_params->G1;
    // Graph *G2 = graph_params->G2;
    // int *G1_labels = graph_params->G1_labels;
    // int **nodes_of_G2Labels = graph_params->nodes_of_G2Labels;
    // int **G2_nodes_of_degree = graph_params->G2_nodes_of_degree;

    // int *mapping = state_params->mapping;
    // int *reverseMapping = state_params->reverseMapping;
    // bool *T2_out = state_params->T2_out;
    // ------------------------
    int coveredNeighborsSize = 0;
    int *coveredNeighbors = findCoveredNeighbors(g1, state, node, &coveredNeighborsSize);

    printf("CoveredNeighbors: ");
    for(int i = 0; i < coveredNeighborsSize; i++) {
        printf("%d ", coveredNeighbors[i]);
    }

    if(coveredNeighborsSize == 0) {
        // A
        int label = g1->nodesToLabel[node];
        int maxSizeCandidates = g2->labelsCardinalities[label];
        // int* candidates = copyArray(g2->labelToNodes[label], sizeCandidates);
        int* candidates = malloc(maxSizeCandidates * sizeof(int));

        *sizeCandidates = 0;
        for(int i = 0; i < maxSizeCandidates; i++) {
            int vertex = g2->labelToNodes[label][i];  // g2->labelToNodes[label] is list of candidates
            if(g2->degrees[vertex] == g1->degrees[node] && state->T2_out[vertex] == 1 && state->mapping2[vertex] != 1) {
                candidates[*sizeCandidates] = vertex;
                *sizeCandidates = *sizeCandidates + 1;
            }
        }
        free(coveredNeighbors);
    
        printf("\nCandidates:");
        for(int i = 0; i < *sizeCandidates; i++) {
            printf("%d ", candidates[i]);
        }
        return candidates;
    }
    
    // B
    printf("entro in B");
    int* commonNodes = (int*)malloc(g2->numVertices * sizeof(int));
    for(int i = 0; i < g2->numVertices; i++) {
        commonNodes[i] = 1;  
    }

    int count = 0;
    for(int i = 0; i < coveredNeighborsSize; i++) {
        int nbrG1 = coveredNeighbors[i];
        int mappedG2 = state->mapping1[nbrG1];

        for(int adjVertex = 0; adjVertex < g2->numVertices; adjVertex++) {
            if(g2->matrix[mappedG2 * g2->numVertices + adjVertex] == 0) {
                commonNodes[adjVertex] = 0;
                count++;
            }
        }
    }

    int commonNodesSize = g2->numVertices - count;
    int* candidates = (int*)malloc(commonNodesSize * sizeof(int));

    int label = g1->nodesToLabel[node];
    int maxSizeCandidates = g2->labelsCardinalities[label];
    *sizeCandidates = 0;
    for(int i = 0; i < maxSizeCandidates; i++) {
        int vertex = g2->labelToNodes[label][i];
        if(commonNodes[vertex] == 1 && g2->degrees[vertex] == g1->degrees[node] && state->mapping2[vertex] != -1) {
            candidates[*sizeCandidates] = vertex;
            *sizeCandidates = *sizeCandidates + 1;
        }
    }
    
    // int label = g1->nodesToLabel[node];
    // for(int i = 0; i < commonNodesSize; i++) {
    //     int* verticies = g2->labelToNodes[label];



    //     if(commonNodes[i] && g2->nodesToLabel[i] == label && g2->degrees[i] == g1->degrees[node] && state->T2_out[i] == 1 && state->mapping2[i] == -1) {
    //         candidates[i] = i;
    //     }
    // }


    
    // int* candidates = (int*)malloc(maxSizeCandidates * sizeof(int));


    free(coveredNeighbors);
    free(commonNodes);
    printf("CommonNodes: ");
    for(int i = 0; i < *sizeCandidates; i++) {
        printf("%d ", candidates[i]);
    }
    return candidates;



    // if (covered_neighbors[0] == -1) {
    //     // A
    //     for (int i = 0; nodes_of_G2Labels[G1_labels[u]][i] != -1; i++) {
    //         int node = nodes_of_G2Labels[G1_labels[u]][i];
    //         if (G2_nodes_of_degree[G1->degrees[u]][node] &&
    //             T2_out[node] && reverseMapping[node] == -1) {
    //             if (!isMultiGraph(G1) || G1->matrix[u][u] == G2->matrix[node][node]) {
    //                 candidates[*candidate_count] = node;
    //                 (*candidate_count)++;
    //             }
    //         }
    //     }
    // } else {
    //     //B
    //     int common_nodes[G2->numVertices];
    //     int common_count = 0;
    //     for (int i = 0; i < G2->numVertices; i++) {
    //         common_nodes[i] = 1;
    //     }
    //     for (int i = 0; covered_neighbors[i] != -1; i++) {
    //         int nbr1 = covered_neighbors[i];
    //         for (int j = 0; j < G2->numVertices; j++) {
    //             if (!G2->matrix[mapping[nbr1]][j]) {
    //                 common_nodes[j] = 0;
    //             }
    //         }
    //     }

    //     for (int i = 0; i < G2->numVertices; i++) {
    //         if (common_nodes[i] &&
    //             G2_nodes_of_degree[G1->degrees[u]][i] &&
    //             nodes_of_G2Labels[G1_labels[u]][i] &&
    //             reverseMapping[i] == -1) {
    //             if (!isMultiGraph(G1) || G1->matrix[u][u] == G2->matrix[i][i]) {
    //                 candidates[*candidate_count] = i;
    //                 (*candidate_count)++;
    //             }
    //         }
    //     }
    // }
    // free(covered_neighbors);
}

int* findCoveredNeighbors(Graph* g, State* state, int node, int* size) {
    int* coveredNeighbors = (int*)malloc(g->degrees[node] * sizeof(int));   // possono essere al massimo maxDegreen (non vero)
    *size = 0;

    if(coveredNeighbors == NULL) {
        printf("Error allocating memory in findCoveredNeighbors\n");
        exit(EXIT_FAILURE);
    }

        for (int adjVertex = 0; adjVertex < g->numVertices; adjVertex++) {
        if (g->matrix[node * g->numVertices + adjVertex] == 1 && state->mapping1[adjVertex] != -1) {
            coveredNeighbors[*size] = adjVertex;
                *size = *size + 1;
            }
        }
    
    return coveredNeighbors;
}

void vf2pp(Graph* g1, Graph* g2, State* state) {
    if (!checkGraphProperties(g1, g2)) {
        printf("Graphs are not isomorphic\n");
        return;
    }

    int* order = ordering(g1, g2);
    
    int sizeCandidates = 0;
    int* candidates = findCandidates(g1, g2, state, order[0], &sizeCandidates);
    
    StackNode* stack = createStack(); 
    push(&stack, candidates, sizeCandidates, order[0]);
    int matchingNode = 1;
    while (!isStackEmpty(stack)) {
        Info* info = peek(&stack);
        bool isMatch = false;

        for(int i = info->candidateIndex; i < info->sizeCandidates; i++) {
            int candidate = info->candidates[i];
            info->candidateIndex = i + 1;

            if(!cutISO(g1, g2, state, info->vertex, candidate)) {
                
                if(isMappingFull(g1, state)) {
                    freeStack(stack);
                    free(order);
                    printf("\nGraphs are isomorphic\n");
                    return;
                }

                state->mapping1[info->vertex] = candidate;
                state->mapping2[candidate] = info->vertex;
                updateState(g1, g2, state, info->vertex, candidate);
                candidates = findCandidates(g1, g2, state, order[matchingNode], &sizeCandidates);
                push(&stack, candidates, sizeCandidates, order[matchingNode]);
                matchingNode++;
                isMatch = true;
                break;
            }
        }

        // no more candidates
        if(!isMatch) {
            Info* tmp = pop(&stack);
            freeInfo(tmp);
            matchingNode--;

            // backtracking
            if(!isStackEmpty(stack)) {
                Info* prevInfo = peek(&stack);
                int candidate = state->mapping1[prevInfo->vertex];
                state->mapping1[prevInfo->vertex] = -1;
                state->mapping2[candidate] = -1;
                restoreState(g1, g2, state, prevInfo->vertex, candidate);
            }
        }
    }
    free(order);
    freeStack(stack);
}

bool cutISO(Graph* g1, Graph* g2, State* state, int node1, int node2){
    int nbrSize1 = 0, nbrSize2 = 0;
    int* neighbors1 = findNeighbors(g1, node1, &nbrSize1);
    int* neighbors2 = findNeighbors(g2, node2, &nbrSize2);

    int* labelsNbr = (int*)malloc(LABELS * sizeof(int));

    if(labelsNbr == NULL) {
        printf("Error allocating memory in cutISO\n");
        exit(EXIT_FAILURE);
    }

    if(!checkLabels(g1, g2, neighbors1, neighbors2, nbrSize1, nbrSize2, labelsNbr)) {   // verifica almeno una corrispondenza G1 a G2
        free(neighbors1);
        free(neighbors2);
        free(labelsNbr);
        return true;
    }

    // se sono qui vuol dire che le label sono uguali (hanno tutti una corrispondenza)
    for(int label = 0; label < LABELS; label++) {
        if (labelsNbr[label] == 1) {    // se è effettivamente una label del vicinato
            int G1NodesOfLabelSize = 0;
            int* G1NodesOfLabel = findNodesOfLabel(neighbors1, g1, label, nbrSize1, &G1NodesOfLabelSize);
            int G2NodesOfLabelSize = 0;
            int* G2NodesOfLabel = findNodesOfLabel(neighbors2, g2, label, nbrSize2, &G2NodesOfLabelSize);

            int count1 = intersectionCount(G1NodesOfLabel, G1NodesOfLabelSize, state->T1);
            int count2 = intersectionCount(G2NodesOfLabel, G2NodesOfLabelSize, state->T2);
            int count3 = intersectionCount(G1NodesOfLabel, G1NodesOfLabelSize, state->T1_out);
            int count4 = intersectionCount(G2NodesOfLabel, G2NodesOfLabelSize, state->T2_out);

            free(G1NodesOfLabel);
            free(G2NodesOfLabel);

            if(count1 != count2 || count3 != count4) {
                free(neighbors1);
                free(neighbors2);
                free(labelsNbr);
                return true;
            }
        }
    }
    free(neighbors1);
    free(neighbors2);
    free(labelsNbr);
    return false;
}

int intersectionCount(int *arr1, int size1, int *arr2) {
    int count = 0;
    for(int i = 0; i < size1; i++) {
        int vertex = arr1[i];
        if(arr2[vertex] == 1) {
            count++;
        }
    }
    return count;
}

int* findNodesOfLabel(int* neighbors, Graph* g, int label, int maxSize, int* size) {
    int* nodes = (int*)malloc(maxSize * sizeof(int));
    
    if (nodes == NULL) {
        printf("Error allocating memory in findNodeOfLabel\n");
        exit(EXIT_FAILURE);
    }

    *size = 0;
    for(int i = 0; i < maxSize; i++) {
        int v = neighbors[i];
        if (g->nodesToLabel[v] == label) {
            nodes[*size] = v;
            *size = *size + 1;
        }
    }

    return nodes;
}

bool checkLabels(Graph*g1, Graph*g2, int* neighbors1, int* neighbors2, int nbr1Size, int nbr2Size, int* labelsNbr) {
    for(int i = 0; i < LABELS; i++) {
        labelsNbr[i] = 0;
    }

    for(int i = 0; i < nbr1Size; i++) {
        int nbr1 = neighbors1[i];
        int labelNbr1 = g1->nodesToLabel[nbr1];
        bool found = false;

        for(int j = 0; j < nbr2Size; j++) {
            int nbr2 = neighbors2[j];
            int labelNbr2 = g2->nodesToLabel[nbr2];

            if(labelNbr1 == labelNbr2)
                found = true;
                labelsNbr[labelNbr1] = 1;
                break;
        }

        if (!found) {
            return false;
        }
    }
    // *labelSize = 0;
    // for(int label = 0; label < LABELS; label++) {
    //     if(labelsNbr[label] == 1) {
    //         *labelSize = *labelSize + 1;
    //     }
    // }
    return true;
}

int* findNeighbors(Graph* g, int node, int* size) {
    int* neighbors = (int*)malloc(g->degrees[node] * sizeof(int));   // prima era maxDegree (in realta no). Se corretto puoi togleire il calcolo del maxDegree
    
    if(neighbors == NULL) {
        printf("Error allocating memory in findNeighbors\n");
        exit(EXIT_FAILURE);
    }

    for (int adjVertex = 0; adjVertex < g->numVertices; adjVertex++) {
        if (g->matrix[node * g->numVertices + adjVertex] == 1) {
            neighbors[*size] = adjVertex;
            *size = *size + 1;
        }
    }
    
    return neighbors;
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

/***** STACK FUNCTIONS *****/
Info* createInfo(int* candidates, int sizeCandidates, int vertex) {
    Info* info = (Info*)malloc(sizeof(Info));
    if (info == NULL) {
        printf("Error allocating memory in createInfo\n");
        exit(EXIT_FAILURE);
    }
    info->vertex = vertex;
    info->candidates = candidates;
    info->sizeCandidates = sizeCandidates;
    info->candidateIndex = 0;
    return info;
}

StackNode* createStackNode(Info* info) {
    StackNode* node = (StackNode*)malloc(sizeof(StackNode));
    if (node == NULL) {
        printf("Error allocating memory in createStackNode\n");
        exit(EXIT_FAILURE);
    }
    node->info = info;
    node->next = NULL;
    return node;
}

void push(StackNode** top, int* candidates, int sizeCandidates, int vertex) {
    Info* info = createInfo(candidates, sizeCandidates, vertex);
    StackNode* node = createStackNode(info);
    node->next = *top;
    *top = node;
}

Info* pop(StackNode** top) {
    if (isStackEmpty(*top)) {
        printf("Stack is empty, cannot pop\n");
        return NULL;
    }
    StackNode* node = *top;
    Info* info = node->info;
    *top = node->next;
    free(node);
    return info;
}

bool isStackEmpty(StackNode* top) {
    return top == NULL;
}

void freeStack(StackNode* top) {
    while (!isStackEmpty(top)) {
        StackNode* node = top;
        top = top->next;
        freeInfo(node->info);
        free(node);
    }
}

void printStack(StackNode* top) {
    StackNode* current = top;
    while (current != NULL) {
        printInfo(current->info);
        current = current->next;
    }
}

void printInfo(Info* info) {
    printf("Vertex: %d\n", info->vertex);
    printf("Index seen: %d\n", info->candidateIndex);
    printf("Candidates: ");
    for (int i = 0; i < info->sizeCandidates; i++) {
        printf("%d ", info->candidates[i]);
    }
    printf("\n");
}

void freeInfo(Info* info) {
    free(info->candidates);
    free(info);
}

StackNode* createStack() {
    return NULL;
}

Info* peek(StackNode** top) {
    return (*top)->info;
}