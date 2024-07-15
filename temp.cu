#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <cuda_runtime.h>
#include <limits.h>
#include <string.h>

#define FILENAME_QUERY "data/graph_query_10000.csv"
#define FILENAME_TARGET "data/graph_target_10000.csv"
#define STREAMS 6
#define LABELS 10
#define INF 99999
#define MAXCONSTMEM 16000

int blockSize;
__constant__ int constMem[MAXCONSTMEM];
__constant__ int constMemVar;

#define CUDA_CHECK_ERROR(err)           \
    if (err != cudaSuccess) {            \
        printf("CUDA error: %s\n", cudaGetErrorString(err)); \
        printf("Error in file: %s, line: %i\n", __FILE__, __LINE__); \
        exit(EXIT_FAILURE);              \
    }

/***** STRUCTS *****/
typedef struct {
    int* matrix;
    int numVertices;
    int* nodesToLabel;
    int** labelToNodes;
    int* labelsCardinalities;
    int* degrees;
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
void initGraphGPU(Graph*);
Graph* createGraph();
void addEdge(Graph*, int, int);
Graph* readGraph(char*);
void printGraph(Graph*);
void freeGraph(Graph*);
void setLabel(Graph*, int, int);

Graph* graphGPU(Graph*);
void freeGraphGPU(Graph*);

/***** STATE PROTOTYPES *****/
State* createStateGPU(Graph*, Graph*, cudaStream_t*);
void freeStateGPU(State*);
void printState(State*, int);
void updateStateGPU(Graph*, Graph*, State*, int, int, cudaStream_t*);
void restoreStateGPU(Graph*, Graph*, State*, int, int, cudaStream_t*);

/***** VF2++ PROTOTYPES *****/
void vf2ppGPU(Graph*, Graph*, State*, Graph*, Graph*, cudaStream_t*);
bool checkGraphPropertiesGPU(Graph*, Graph*, Graph*, Graph*, cudaStream_t*);
int compare(const void*, const void*);
int* orderingGPU(Graph*, Graph*, cudaStream_t*, Graph*);
void findRootGPU(Graph*, int*, int*, int*, int*, cudaStream_t*);
void processDepthGPU(Graph*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, cudaStream_t*, Graph*);
int* findCandidatesGPU(Graph*, Graph*, State*, int, int*, Graph*, Graph*, cudaStream_t*);
bool cutISOGPU(Graph*, Graph*, State*, int, int, Graph*, Graph*, cudaStream_t*);

/***** STACK PROTOTYPES *****/
Info* createInfo(int* candidates, int sizeCandidates, int vertex);
StackNode* createStackNode(Info*);
void push(StackNode**, Info*);
Info* pop(StackNode**);
bool isStackEmpty(StackNode*);
void freeStack(StackNode*);
void printStack(StackNode*);
void printInfo(Info*);
void freeInfo(Info*);
StackNode* createStack();
Info* peek(StackNode**);

int main() {
    blockSize = 256;

    Graph* h_g1 = readGraph(FILENAME_QUERY);
    Graph* h_g2 = readGraph(FILENAME_TARGET);
    Graph* d_g1 = graphGPU(h_g1);
    Graph* d_g2 = graphGPU(h_g2);

    cudaStream_t streams[STREAMS];

    for(int i = 0; i < STREAMS; i++) {
        cudaStreamCreate(&streams[i]);
    }

    State* d_state = createStateGPU(d_g1, d_g2, streams);

    vf2ppGPU(d_g1, d_g2, d_state, h_g1, h_g2, streams);

    int* mapping1 = (int*)malloc(d_g1->numVertices * sizeof(int));

    if(mapping1 == NULL) {
        printf("Error allocating memory in main\n");
        exit(EXIT_FAILURE);
    }

    CUDA_CHECK_ERROR(cudaMemcpy(mapping1, d_state->mapping1, d_g1->numVertices * sizeof(int), cudaMemcpyDeviceToHost));

    printf("Mapping\n");
    for(int i = 0; i < d_g1->numVertices; i++) {
        printf("%d -> %d\n", i, mapping1[i]);
    }

    for(int i = 0; i < STREAMS; i++) {
        cudaStreamDestroy(streams[i]);
    }

    freeGraph(h_g1);
    freeGraph(h_g2);
    freeGraphGPU(d_g1);
    freeGraphGPU(d_g2);

    freeStateGPU(d_state);

    return EXIT_SUCCESS;
}

__global__ void initArrayKernel(int* d_array, int size, int value) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < size) {
        d_array[idx] = value;
    }
}

__global__ void initMatrixKernel(int* d_matrix, int V, int value) {
    int col = threadIdx.y + blockIdx.y * blockDim.y;
    int row = threadIdx.x + blockIdx.x * blockDim.x;

    if (row < V && col < V) {
        d_matrix[row * V + col] = value;
    }
}

/***** GRAPH FUNCTIONS *****/
void initGraphGPU(Graph* g) {
    g->matrix = (int*)malloc(g->numVertices * g->numVertices * sizeof(int));
    g->nodesToLabel = (int*)malloc(g->numVertices * sizeof(int));
    g->labelsCardinalities = (int*)malloc(LABELS * sizeof(int));
    g->labelToNodes = (int**)malloc(LABELS * sizeof(int*));
    g->degrees = (int*)malloc(g->numVertices * sizeof(int));

    if (g->nodesToLabel == NULL || g->labelsCardinalities == NULL || g->labelToNodes == NULL || g->degrees == NULL || g->matrix == NULL) {
        printf("Error allocating memory in initGraph\n");
        exit(EXIT_FAILURE);
    }

    int *d_nodesToLabel, *d_degrees, *d_matrix;

    cudaStream_t stream1, stream2, stream3;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);
    cudaStreamCreate(&stream3);

    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_nodesToLabel, g->numVertices * sizeof(int)));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_degrees, g->numVertices * sizeof(int)));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_matrix, g->numVertices * g->numVertices * sizeof(int)));

    int gridSize = (g->numVertices + blockSize - 1) / blockSize;

    initArrayKernel<<<gridSize, blockSize, 0, stream1>>>(d_nodesToLabel, g->numVertices, -1);
    initArrayKernel<<<gridSize, blockSize, 0, stream2>>>(d_degrees, g->numVertices, 0);

    int gridSizeX = (g->numVertices + blockSize - 1) / blockSize;
    int gridSizeY = (g->numVertices + blockSize - 1) / blockSize;
    dim3 gridSizeM(gridSizeX, gridSizeY);

    dim3 blockSizeM(blockSize, blockSize);
    initMatrixKernel<<<gridSizeM, blockSizeM, 0, stream3>>>(d_matrix, g->numVertices, 0);

    for (int label = 0; label < LABELS; label++) {
        g->labelsCardinalities[label] = 0;
        g->labelToNodes[label] = (int*)malloc(g->numVertices * sizeof(int));
    }

    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);
    cudaStreamSynchronize(stream3);

    CUDA_CHECK_ERROR(cudaMemcpy(g->nodesToLabel, d_nodesToLabel, g->numVertices * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK_ERROR(cudaMemcpy(g->degrees, d_degrees, g->numVertices * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK_ERROR(cudaMemcpy(g->matrix, d_matrix, g->numVertices * g->numVertices * sizeof(int), cudaMemcpyDeviceToHost));

    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);

    cudaFree(d_nodesToLabel);
    cudaFree(d_degrees);
    cudaFree(d_matrix);
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

    initGraphGPU(g);

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

Graph* graphGPU(Graph* h_g) {
    Graph* d_g = createGraph();

    size_t sizeMatrix = h_g->numVertices * h_g->numVertices * sizeof(int);
    size_t size1 = h_g->numVertices * sizeof(int);
    size_t size2 = LABELS * sizeof(int);

    d_g->numVertices = h_g->numVertices;
    d_g->labelToNodes = (int**)malloc(LABELS * sizeof(int*));   // it is a vector on host wich contains pointers to vectors on device

    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_g->matrix, sizeMatrix));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_g->nodesToLabel, size1));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_g->labelsCardinalities, size2));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_g->degrees, size1));

    for(int label = 0; label < LABELS; label++) {
        CUDA_CHECK_ERROR(cudaMalloc((void**)&d_g->labelToNodes[label], h_g->labelsCardinalities[label] * sizeof(int))); // each vector on device has a size equal to the cardinality of the label
    }

    CUDA_CHECK_ERROR(cudaMemcpy(d_g->matrix, h_g->matrix, sizeMatrix, cudaMemcpyHostToDevice));
    CUDA_CHECK_ERROR(cudaMemcpy(d_g->nodesToLabel, h_g->nodesToLabel, size1, cudaMemcpyHostToDevice));
    CUDA_CHECK_ERROR(cudaMemcpy(d_g->labelsCardinalities, h_g->labelsCardinalities, size2, cudaMemcpyHostToDevice));
    CUDA_CHECK_ERROR(cudaMemcpy(d_g->degrees, h_g->degrees, size1, cudaMemcpyHostToDevice));

    for(int label = 0; label < LABELS; label++) {
        CUDA_CHECK_ERROR(cudaMemcpy(d_g->labelToNodes[label], h_g->labelToNodes[label], h_g->labelsCardinalities[label] * sizeof(int), cudaMemcpyHostToDevice));
    }

    return d_g;
}

void freeGraphGPU(Graph* g) {
    for(int label = 0; label < LABELS; label++) {
        cudaFree(g->labelToNodes[label]);
    }
    free(g->labelToNodes);
    cudaFree(g->matrix);
    cudaFree(g->nodesToLabel);
    cudaFree(g->labelsCardinalities);
    cudaFree(g->degrees);
    free(g);
    g = NULL;
}

/***** STATE FUNCTIONS *****/
State* createStateGPU(Graph* g1, Graph* g2, cudaStream_t* streams) {
    State* s = (State*)malloc(sizeof(State));

    if (s == NULL) {
        printf("Error allocating memory in createStateGPU\n");
        exit(EXIT_FAILURE);
    }

    size_t size1 = g1->numVertices * sizeof(int);
    size_t size2 = g2->numVertices * sizeof(int);

    CUDA_CHECK_ERROR(cudaMalloc((void**)&s->mapping1, size1));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&s->mapping2, size2));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&s->T1, size1));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&s->T2, size2));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&s->T1_out, size1));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&s->T2_out, size2));

    int gridSize1 = (g1->numVertices + blockSize - 1) / blockSize;
    int gridSize2 = (g2->numVertices + blockSize - 1) / blockSize;

    cudaStream_t stream1 = streams[0];
    cudaStream_t stream2 = streams[1];
    cudaStream_t stream3 = streams[2];

    initArrayKernel<<<gridSize1, blockSize, 0, stream1>>>(s->mapping1, g1->numVertices, -1);
    initArrayKernel<<<gridSize1, blockSize, 0, stream2>>>(s->T1, g1->numVertices, -1);
    initArrayKernel<<<gridSize1, blockSize, 0, stream3>>>(s->T1_out, g1->numVertices, 1);

    initArrayKernel<<<gridSize2, blockSize, 0, stream1>>>(s->mapping2, g2->numVertices, -1);
    initArrayKernel<<<gridSize2, blockSize, 0, stream2>>>(s->T2, g2->numVertices, -1);
    initArrayKernel<<<gridSize2, blockSize, 0, stream3>>>(s->T2_out, g2->numVertices, 1);

    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);
    cudaStreamSynchronize(stream3);

    return s;
}

void freeStateGPU(State* s) {
    cudaFree(s->mapping1);
    cudaFree(s->mapping2);
    cudaFree(s->T1);
    cudaFree(s->T2);
    cudaFree(s->T1_out);
    cudaFree(s->T2_out);
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

__global__ void updateStateKernel(int* matrix, int V, int* mapping, int* T, int* T_out, int node1, int node2) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if( idx < V) {
        if(idx == node1) {
            mapping[idx] = node2;
        }

        __syncthreads();

        if(matrix[node1 * V + idx] == 1 && mapping[idx] == -1) {
            T[idx] = 1;
            T_out[idx] = -1;
        }

        __syncthreads();

        if(idx == node1) {
            T[idx] = -1;
            T_out[idx] = -1;
        }
    }
}

void updateStateGPU(Graph* d_g1, Graph* d_g2, State* d_state, int node, int candidate, cudaStream_t* streams) {
    cudaStream_t stream1 = streams[0];
    cudaStream_t stream2 = streams[1];

    int gridSize = (d_g1->numVertices + blockSize - 1) / blockSize;
    updateStateKernel<<<gridSize, blockSize, 0, stream1>>>(d_g1->matrix, d_g1->numVertices, d_state->mapping1, d_state->T1, d_state->T1_out, node, candidate);

    gridSize = (d_g2->numVertices + blockSize - 1) / blockSize;
    updateStateKernel<<<gridSize, blockSize, 0, stream2>>>(d_g2->matrix, d_g2->numVertices, d_state->mapping2, d_state->T2, d_state->T2_out, candidate, node);
}

__global__ void restoreStateKernel(int* matrix, int V, int node, int* T, int* T_out, int offset, int usage, int* mapping) {  // offset introduced because of the two streams
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= V) {
        return;
    }

    int isAdded = 0;

    if(matrix[node * V + idx] == 1) {

        if(usage == 1) {

            if(constMem[offset + idx] != -1) {  // constMem contains mapping1 if offset is 0, mapping2 if offset is V
                atomicExch(&T[node], 1);
                isAdded = 1;
            }
            else {
                int hasCoveredNeighbor = 0;
                for(int adjVertex2 = 0; adjVertex2 < V; adjVertex2++) {
                    if(matrix[idx * V + adjVertex2] == 1 && constMem[offset + adjVertex2] != -1) {   // constMem contains mapping1 if offset is 0, mapping2 if offset is V
                        hasCoveredNeighbor = 1;
                        break;
                    }
                }

                if(hasCoveredNeighbor == 0) {
                    T[idx] = -1;
                    T_out[idx] = 1;
                }
            }

        } else {

            if(mapping[idx] != -1) { 
                atomicExch(&T[node], 1);
                isAdded = 1;
            }
            else {
                int hasCoveredNeighbor = 0;
                for(int adjVertex2 = 0; adjVertex2 < V; adjVertex2++) {
                    if(matrix[idx * V + adjVertex2] == 1 && mapping[adjVertex2] != -1) {   
                        hasCoveredNeighbor = 1;
                        break;
                    }
                }

                if(hasCoveredNeighbor == 0) {
                    T[idx] = -1;
                    T_out[idx] = 1;
                }
            }

        }
    }

    if(isAdded == 0) {
        atomicExch(&T_out[node], 1);
    }
}

void restoreStateGPU(Graph* d_g1, Graph* d_g2, State* d_state, int node, int candidate, cudaStream_t* streams) {
    cudaStream_t stream1 = streams[0];
    cudaStream_t stream2 = streams[1];

    size_t size1 = d_g1->numVertices * sizeof(int);
    size_t size2 = d_g2->numVertices * sizeof(int);
    int value = -1;

    int useConstMem = d_g1->numVertices + d_g2->numVertices < MAXCONSTMEM;

    CUDA_CHECK_ERROR(cudaMemcpyAsync(d_state->mapping1 + node, &value, sizeof(int), cudaMemcpyHostToDevice, stream1));
    CUDA_CHECK_ERROR(cudaMemcpyAsync(d_state->mapping2 + candidate, &value, sizeof(int), cudaMemcpyHostToDevice, stream2));

    if(useConstMem) {
        CUDA_CHECK_ERROR(cudaMemcpyToSymbolAsync(constMem, d_state->mapping1, size1, 0, cudaMemcpyDeviceToDevice, stream1));
        CUDA_CHECK_ERROR(cudaMemcpyToSymbolAsync(constMem, d_state->mapping2, size2, size1, cudaMemcpyDeviceToDevice, stream2));
    }
    
    int gridSize = (d_g1->numVertices + blockSize - 1) / blockSize;
    restoreStateKernel<<<gridSize, blockSize, 0, stream1>>>(d_g1->matrix, d_g1->numVertices, node, d_state->T1, d_state->T1_out, 0, useConstMem, d_state->mapping1);
    
    gridSize = (d_g2->numVertices + blockSize - 1) / blockSize;
    restoreStateKernel<<<gridSize, blockSize, 0, stream2>>>(d_g2->matrix, d_g2->numVertices, candidate, d_state->T2, d_state->T2_out, d_g1->numVertices, useConstMem, d_state->mapping2);

    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);
}

/***** VF2++ FUNCTIONS *****/
__global__ void equalKernel(int* arr1, int* arr2, int size, int* ret) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < size) {
        if (arr1[idx] != arr2[idx]) {
            *ret = 0;
        }
    }
}

int compare(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

bool checkGraphPropertiesGPU(Graph* h_g1, Graph* h_g2, Graph* d_g1, Graph* d_g2, cudaStream_t* streams) {
    if (h_g1->numVertices != h_g2->numVertices || h_g1->numVertices == 0 || h_g2->numVertices == 0) {
        return false;
    }

    int h_ret1 = 1, h_ret2 = 1;
    int* d_ret1, *d_ret2;

    cudaStream_t stream1 = streams[0];
    cudaStream_t stream2 = streams[1];

    // first check: the cardinalities of the labels must be the same
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_ret1, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMemcpyAsync(d_ret1, &h_ret1, sizeof(int), cudaMemcpyHostToDevice, stream1));

    int gridSize1 = (LABELS + blockSize - 1) / blockSize;
    equalKernel<<<gridSize1, blockSize, 0, stream1>>>(d_g1->labelsCardinalities, d_g2->labelsCardinalities, LABELS, d_ret1);

    // second check: the sequence of the degrees must be the same
    int size = h_g1->numVertices;
    int *d_tmp1, *d_tmp2;

    int* tmp1 = (int*)malloc(size * sizeof(int));
    int* tmp2 = (int*)malloc(size * sizeof(int));

    memcpy(tmp1, h_g1->degrees, size * sizeof(int));
    memcpy(tmp2, h_g2->degrees, size * sizeof(int));
    
    qsort(tmp1, size, sizeof(int), compare);
    qsort(tmp2, size, sizeof(int), compare);

    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_tmp1, size * sizeof(int)));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_tmp2, size * sizeof(int)));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_ret2, sizeof(int)));

    CUDA_CHECK_ERROR(cudaMemcpyAsync(d_tmp1, tmp1, size * sizeof(int), cudaMemcpyHostToDevice, stream2));
    CUDA_CHECK_ERROR(cudaMemcpyAsync(d_tmp2, tmp2, size * sizeof(int), cudaMemcpyHostToDevice, stream2));
    CUDA_CHECK_ERROR(cudaMemcpyAsync(d_ret2, &h_ret2, sizeof(int), cudaMemcpyHostToDevice, stream2));
   
    int gridSize2 = (size + blockSize - 1) / blockSize;
    equalKernel<<<gridSize2, blockSize, 0, stream2>>>(d_tmp1, d_tmp2, size, d_ret2);

    CUDA_CHECK_ERROR(cudaMemcpyAsync(&h_ret1, d_ret1, sizeof(int), cudaMemcpyDeviceToHost, stream1));
    CUDA_CHECK_ERROR(cudaMemcpyAsync(&h_ret2, d_ret2, sizeof(int), cudaMemcpyDeviceToHost, stream2));

    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);

    cudaFree(d_ret1);
    cudaFree(d_ret2);
    cudaFree(d_tmp2);
    cudaFree(d_tmp1);
    
    free(tmp1);
    free(tmp2);

    return h_ret1 && h_ret2;
}

__global__ void bfsKernel(int* matrix, int V, int* levels, int* d_done, int depth) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    __shared__ int s_done;
    extern __shared__ int s_levels[];

    if(threadIdx.x == 0) {
        s_done = 0;
    }

    // it copies the whole levels array in shared memory
    // example: thread 0 loads levels[0] and levels[4] if blockSize is 4
    for(int i = threadIdx.x; i < V; i += blockDim.x) {
      s_levels[i] = levels[i];
    }

    __syncthreads();

    // levels is used as visited too
    if(idx < V && s_levels[idx] == depth) {    // it blocks all thread with size greater than V and all threads not at the current depth

        for(int adjVertex = 0; adjVertex < V; adjVertex++) {
            if(matrix[idx * V + adjVertex] == 1 && s_levels[adjVertex] == -1) {
                atomicExch(&levels[adjVertex], depth + 1);
                s_done = 1;
            }
        }
    }

    __syncthreads();

    if(threadIdx.x == 0) {
        atomicExch(d_done, s_done);
    }
}

__global__ void maxRarityConstMemKernel(int V, int* d_nodesToLabel, int* d_maxRarity, int* d_is_good) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ int s_data[];

    if(idx < V && d_is_good[idx]) {
        s_data[threadIdx.x] = constMem[d_nodesToLabel[idx]];    // constMem contains labelRarity
    } else {
        s_data[threadIdx.x] = INF;
    }

    __syncthreads();

    // parallel reduction: at each step, the block is halved and each thread computes the min of its value with the value of the other one at distance s in
    // the same block
    for(unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if(threadIdx.x < s) {
            s_data[threadIdx.x] = min(s_data[threadIdx.x], s_data[threadIdx.x + s]);   // At the end of the loop, the min value is in sdata[0]
        }
        __syncthreads();
    }

    // each thread 0 of each block computes the global min
    if(threadIdx.x == 0) {
        atomicMin(d_maxRarity, s_data[0]);
    }
}

__global__ void maxRarityConstMemFilterKernel(int V, int* d_nodesToLabel, int* d_maxRarity, int* d_is_good) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx < V) {
        if(constMem[d_nodesToLabel[idx]] != *d_maxRarity) {     // constMem contains labelRarity
            d_is_good[idx] = 0;
        }
    }
}

__global__ void maxDegreeKernel(int V, int* d_degrees, int* d_maxDegree, int* d_is_good) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ int s_data[];

    if(idx < V && d_is_good[idx]) {
        s_data[threadIdx.x] = d_degrees[idx];
    }
    else {
        s_data[threadIdx.x] = -INF;
    }

    __syncthreads();

    for(unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if(threadIdx.x < s) {
            s_data[threadIdx.x] = max(s_data[threadIdx.x], s_data[threadIdx.x + s]);
        }
        __syncthreads();
    }

    if(threadIdx.x == 0) {
        atomicMax(d_maxDegree, s_data[0]);
    }
}

__global__ void maxDegreeFilterKernel(int V, int* d_degrees, int* d_maxDegree, int* d_is_good) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx < V) {
        if(d_degrees[idx] != *d_maxDegree) {
            d_is_good[idx] = 0;
        }
    }
}

__global__ void findNodeKernel(int V, int* is_good, int* d_node) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < V && is_good[idx]) {
        atomicExch(d_node, idx);
    }
}

void findRootGPU(Graph* d_g, int* d_root, int* d_is_good, int* d_maxRarity, int* d_maxDegree, cudaStream_t* streams) {
    cudaStream_t stream1 = streams[0];
    cudaStream_t stream2 = streams[1];

    int h_maxRarity = INF, h_maxDegree = -INF;
    int gridSize = (d_g->numVertices + blockSize - 1) / blockSize;
    size_t size2 = LABELS * sizeof(int);

    initArrayKernel<<<gridSize, blockSize, 0, stream1>>>(d_is_good, d_g->numVertices, 1);
    CUDA_CHECK_ERROR(cudaMemcpyAsync(d_maxRarity, &h_maxRarity, sizeof(int), cudaMemcpyHostToDevice, stream2));
    CUDA_CHECK_ERROR(cudaMemcpyAsync(d_maxDegree, &h_maxDegree, sizeof(int), cudaMemcpyHostToDevice, stream1));
    CUDA_CHECK_ERROR(cudaMemcpyToSymbolAsync(constMem, d_g->labelsCardinalities, size2, 0, cudaMemcpyDeviceToDevice, stream2)); // constMem contains labelRarity (d_g1->labelsCardinalities)

    size_t sharedMemSize = blockSize * sizeof(int); // each block has a shared memory of size blockSize

    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);

    maxRarityConstMemKernel<<<gridSize, blockSize, sharedMemSize, stream1>>>(d_g->numVertices, d_g->nodesToLabel, d_maxRarity, d_is_good);
    maxRarityConstMemFilterKernel<<<gridSize, blockSize, 0, stream1>>>(d_g->numVertices, d_g->nodesToLabel, d_maxRarity, d_is_good);
    maxDegreeKernel<<<gridSize, blockSize, sharedMemSize, stream1>>>(d_g->numVertices, d_g->degrees, d_maxDegree, d_is_good);
    maxDegreeFilterKernel<<<gridSize, blockSize, 0, stream1>>>(d_g->numVertices, d_g->degrees, d_maxDegree, d_is_good);
    findNodeKernel<<<gridSize, blockSize, 0, stream1>>>(d_g->numVertices, d_is_good, d_root);
}

__global__ void findLevelNodesKernel(int* levels, int depth, int* levelNodes, int V, int* levelSize) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    __shared__ int s_levelSize;

    if(threadIdx.x == 0) {
        s_levelSize = 0;
    }

    __syncthreads();

    if(idx < V && levels[idx] == depth) {
        atomicAdd(&s_levelSize, 1);
        levelNodes[idx] = 1;
    }

    __syncthreads();

    if(threadIdx.x == 0) {
        atomicAdd(levelSize, s_levelSize);
    }
}

__global__ void initBfsKernel(int* levels, int V, int* root, int depth) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(idx < V) {
        if(idx == *root) {
            levels[idx] = depth;
        }
        else {
            levels[idx] = -1;
        }
    }
}

// __global__ void printArrayKernel(int* arr, int V) {
//     int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     if( idx == 0)
//       for(int i = 0; i < V; i++)
//         printf("%d ", arr[i]);
// }

// __global__ void printVarKernel(int* d) {
//     printf(" .%d. ", *d);
// }

int* orderingGPU(Graph* d_g1, Graph* d_g2, cudaStream_t* streams, Graph* h_g1) {
    // findRootGPU
    int *d_root, *d_is_good, *d_maxRarity, *d_maxDegree;

    // BFS
    int* d_levels, *d_done;
    int* h_done;
    int depth = 0;

    // findLevelNodesKernel
    int *d_levelNodes, *d_levelSize;

    // processDepth
    int *d_V1Unordered, *d_labelRarity, *d_connectivityG1, *d_maxConnectivity;  // d_root is reused as d_nextNode
    int* order = (int*)malloc(d_g1->numVertices * sizeof(int));   // order of the nodes of g1
    int order_index = 0;

    size_t size1 = d_g1->numVertices * sizeof(int);
    size_t size2 = LABELS * sizeof(int);
    int gridSize = (d_g1->numVertices + blockSize - 1) / blockSize;

    cudaStream_t stream1 = streams[0];
    cudaStream_t stream2 = streams[1];
    cudaStream_t stream3 = streams[2];
    cudaStream_t stream4 = streams[3];
    
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_root, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_is_good, size1));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_maxRarity, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_maxDegree, sizeof(int)));
    
    findRootGPU(d_g1, d_root, d_is_good, d_maxRarity, d_maxDegree, streams);

    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_levels, size1));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_done, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_done, sizeof(int))); // pinned memory for faster host-device communication

    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_levelNodes, size1));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_levelSize, sizeof(int)));
    
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_V1Unordered, size1));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_labelRarity, size2));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_connectivityG1, size1));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_maxConnectivity, sizeof(int)));

    initArrayKernel<<<gridSize, blockSize, 0, stream2>>>(d_V1Unordered, d_g1->numVertices, -1); // stream2 empty because already sync in findRootGPU
    CUDA_CHECK_ERROR(cudaMemcpyAsync(d_labelRarity, d_g1->labelsCardinalities, size2, cudaMemcpyDeviceToDevice, stream3));
    CUDA_CHECK_ERROR(cudaMemsetAsync(d_connectivityG1, 0, size1, stream4));

    initBfsKernel<<<gridSize, blockSize, 0, stream1>>>(d_levels, d_g1->numVertices, d_root, depth); // executed after findRootGPU because of same stream
    
    size_t sharedMemSize = size1;

    do {
        *h_done = 0;
        CUDA_CHECK_ERROR(cudaMemsetAsync(d_done, 0, sizeof(int), stream1));

        bfsKernel<<<gridSize, blockSize, sharedMemSize, stream1>>>(d_g1->matrix, d_g1->numVertices, d_levels, d_done, depth);

        CUDA_CHECK_ERROR(cudaMemcpyAsync(h_done, d_done, sizeof(int), cudaMemcpyDeviceToHost, stream1));
        depth++;
       
        cudaStreamSynchronize(stream1);
    } while(*h_done);
  
    cudaStreamSynchronize(stream2); // bfs at stream1 is implicitly synchronized 
    cudaStreamSynchronize(stream3);
    cudaStreamSynchronize(stream4);
    
    for (int d = 0; d < depth; d++) {
        CUDA_CHECK_ERROR(cudaMemsetAsync(d_levelNodes, 0, size1, stream2));
        CUDA_CHECK_ERROR(cudaMemsetAsync(d_levelSize, 0, sizeof(int), stream1));

        cudaStreamSynchronize(stream2);
        
        findLevelNodesKernel<<<gridSize, blockSize, 0, stream1>>>(d_levels, d, d_levelNodes, d_g1->numVertices, d_levelSize);
        
        processDepthGPU(d_g1, order, &order_index, d_connectivityG1, d_labelRarity, d_V1Unordered, d_levelNodes, d_levelSize, 
                        d_is_good, d_maxRarity, d_maxDegree, d_maxConnectivity, d_root, streams, h_g1);
    }

    // findRootGPU
    cudaFree(d_root);
    cudaFree(d_is_good);
    cudaFree(d_maxRarity);
    cudaFree(d_maxDegree);

    //BFS
    cudaFreeHost(h_done);
    cudaFree(d_done);
    cudaFree(d_levels);
    
    // findLevelNodesKernel
    cudaFree(d_levelNodes);
    cudaFree(d_levelSize);

    // processDepth
    cudaFree(d_V1Unordered);
    cudaFree(d_labelRarity);
    cudaFree(d_connectivityG1);
    cudaFree(d_maxConnectivity);

    return order;
}

__global__ void maxConnectivityKernel(int* d_connectivityG1, int* d_maxConnectivity, int* d_is_good, int V) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ int s_data[];

    if(idx < V && d_is_good[idx]) {         // d_is_good already contains the information about the nodes of the current level
        int conn = d_connectivityG1[idx];
        s_data[threadIdx.x] = conn;
    } else {
        s_data[threadIdx.x] = -INF;
    }

    __syncthreads();

    for(unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if(threadIdx.x < s) {
            s_data[threadIdx.x] = max(s_data[threadIdx.x], s_data[threadIdx.x + s]);
        }
        __syncthreads();
    }

    if(threadIdx.x == 0) {
        atomicMax(d_maxConnectivity, s_data[0]);
    }
}

__global__ void maxConnectivityFilterKernel(int* d_connectivityG1, int* d_maxConnectivity, int* d_is_good, int V) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx < V) {
        if(d_connectivityG1[idx] != *d_maxConnectivity) {
            d_is_good[idx] = 0;
        }
    }
}

__global__ void unorderedFilterKernel(int* d_is_good, int* d_V1Unordered, int* d_levelNodes, int V) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx >= V) {
        return;
    }

    if(d_levelNodes[idx]) {
        if(d_V1Unordered[idx] == 1) {
            d_is_good[idx] = 0;
        }
    } else {
        d_is_good[idx] = 0;
    }
}

__global__ void updateConnKernel(int* matrix, int V, int* connectivity, int node) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx < V) {
        if(matrix[node * V + idx] == 1) {
            atomicAdd(&connectivity[idx], 1);
        }
    }
}

// we need to update labelRarity so we cannot use constMem in order to avoid overhead
__global__ void maxRarityKernel(int V, int* d_labelRarity, int* d_nodesToLabel, int* d_maxRarity, int* d_is_good) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ int s_data[];

    if(idx < V && d_is_good[idx]) {
        s_data[threadIdx.x] = d_labelRarity[d_nodesToLabel[idx]];
    } else {
        s_data[threadIdx.x] = INF;
    }

    __syncthreads();

    // parallel reduction: at each step, the block is halved and each thread computes the min of its value with the value of the other one at distance s in
    // the same block
    for(unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if(threadIdx.x < s) {
            s_data[threadIdx.x] = min(s_data[threadIdx.x], s_data[threadIdx.x + s]);   // At the end of the loop, the min value is in sdata[0]
        }
        __syncthreads();
    }

    // each thread 0 of each block computes the global min
    if(threadIdx.x == 0) {
        atomicMin(d_maxRarity, s_data[0]);
    }
}

__global__ void maxRarityFilterKernel(int V, int* d_labelRarity, int* d_nodesToLabel, int* d_maxRarity, int* d_is_good) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx < V) {
        if(d_labelRarity[d_nodesToLabel[idx]] != *d_maxRarity) {
            d_is_good[idx] = 0;
        }
    }
}

void processDepthGPU(Graph* d_g, int* order, int* order_index, int* d_connectivityG1, int* d_labelRarity, int* d_V1Unordered, int* d_levelNodes,
                    int* d_levelSize, int* d_is_good, int* d_maxRarity, int* d_maxDegree, int* d_maxConnectivity, int* d_nextNode, cudaStream_t* streams,
                    Graph* h_g) {
    
    size_t size1 = d_g->numVertices * sizeof(int);
    size_t size2 = LABELS * sizeof(int);

    int h_levelSize;
    int* h_nodesToLabel = h_g->nodesToLabel;

    cudaStream_t stream1 = streams[0];
    cudaStream_t stream2 = streams[1];
    cudaStream_t stream3 = streams[2];
    cudaStream_t stream4 = streams[3];

    CUDA_CHECK_ERROR(cudaMemcpyAsync(&h_levelSize, d_levelSize, sizeof(int), cudaMemcpyDeviceToHost, stream1));

    // pinned memory for faster host-device communication
    int *h_maxRarity, *h_maxDegree, *h_maxConnectivity, *h_nextNode; 
    int* h_labelRarity, *h_V1Unordered;

    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_maxRarity, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_maxDegree, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_maxConnectivity, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_nextNode, sizeof(int)));

    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_labelRarity, size2));
    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_V1Unordered, size1));       

    *h_maxRarity = INF; *h_maxDegree = -INF; *h_maxConnectivity = -INF;

    int gridSize = (d_g->numVertices + blockSize - 1) / blockSize;
    size_t sharedMemSize = blockSize * sizeof(int);

    cudaStreamSynchronize(stream1);

    while(h_levelSize > 0) {
        initArrayKernel<<<gridSize, blockSize, 0, stream1>>>(d_is_good, d_g->numVertices, 1);
        CUDA_CHECK_ERROR(cudaMemcpyAsync(d_maxRarity, h_maxRarity, sizeof(int), cudaMemcpyHostToDevice, stream2));
        CUDA_CHECK_ERROR(cudaMemcpyAsync(d_maxDegree, h_maxDegree, sizeof(int), cudaMemcpyHostToDevice, stream3));
        CUDA_CHECK_ERROR(cudaMemcpyAsync(d_maxConnectivity, h_maxConnectivity, sizeof(int), cudaMemcpyHostToDevice, stream4));

        cudaStreamSynchronize(stream2);
        cudaStreamSynchronize(stream3);
        cudaStreamSynchronize(stream4);

        unorderedFilterKernel<<<gridSize, blockSize, 0, stream1>>>(d_is_good, d_V1Unordered, d_levelNodes, d_g->numVertices);
        maxConnectivityKernel<<<gridSize, blockSize, sharedMemSize, stream1>>>(d_connectivityG1, d_maxConnectivity, d_is_good, d_g->numVertices);
        maxConnectivityFilterKernel<<<gridSize, blockSize, 0, stream1>>>(d_connectivityG1, d_maxConnectivity, d_is_good, d_g->numVertices);
        maxDegreeKernel<<<gridSize, blockSize, sharedMemSize, stream1>>>(d_g->numVertices, d_g->degrees, d_maxDegree, d_is_good);   
        maxDegreeFilterKernel<<<gridSize, blockSize, 0, stream1>>>(d_g->numVertices, d_g->degrees, d_maxDegree, d_is_good);
        maxRarityKernel<<<gridSize, blockSize, sharedMemSize, stream1>>>(d_g->numVertices, d_labelRarity, d_g->nodesToLabel, d_maxRarity, d_is_good);
        maxRarityFilterKernel<<<gridSize, blockSize, 0, stream1>>>(d_g->numVertices, d_labelRarity, d_g->nodesToLabel, d_maxRarity, d_is_good);
        findNodeKernel<<<gridSize, blockSize, 0, stream1>>>(d_g->numVertices, d_is_good, d_nextNode);
        
        CUDA_CHECK_ERROR(cudaMemcpyAsync(h_nextNode, d_nextNode, sizeof(int), cudaMemcpyDeviceToHost, stream1));

        CUDA_CHECK_ERROR(cudaMemcpyAsync(h_labelRarity, d_labelRarity, size2, cudaMemcpyDeviceToHost, stream2));        
        CUDA_CHECK_ERROR(cudaMemcpyAsync(h_V1Unordered, d_V1Unordered, size1, cudaMemcpyDeviceToHost, stream3)); 
        
        cudaStreamSynchronize(stream1);
        
        updateConnKernel<<<gridSize, blockSize, 0, stream1>>>(d_g->matrix, d_g->numVertices, d_connectivityG1, *h_nextNode);
        order[(*order_index)++] = *h_nextNode;
        h_levelSize--;

        cudaStreamSynchronize(stream2);
        cudaStreamSynchronize(stream3);

        h_labelRarity[h_nodesToLabel[*h_nextNode]]--;
        h_V1Unordered[*h_nextNode] = 1;

        CUDA_CHECK_ERROR(cudaMemcpyAsync(d_labelRarity, h_labelRarity, size2, cudaMemcpyHostToDevice, stream2));
        CUDA_CHECK_ERROR(cudaMemcpyAsync(d_V1Unordered, h_V1Unordered, size1, cudaMemcpyHostToDevice, stream3));

        cudaStreamSynchronize(stream1);
        cudaStreamSynchronize(stream2);
        cudaStreamSynchronize(stream3);
    }

    cudaFreeHost(h_maxRarity);
    cudaFreeHost(h_maxDegree);
    cudaFreeHost(h_maxConnectivity);
    cudaFreeHost(h_nextNode);
    cudaFreeHost(h_labelRarity);
    cudaFreeHost(h_V1Unordered);
}

__global__ void findCoveredNeighborsKernel(int* matrix1, int* mapping1, int node, int* coveredNeighbors, int* coveredNeighborsSize,
                                            int numVertices) {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx >= numVertices)
        return;

    if(matrix1[node * numVertices + idx] == 1 && mapping1[idx] != -1) {
        int index = atomicAdd(coveredNeighborsSize, 1);
        coveredNeighbors[index] = idx;
    }
}

__global__ void findCandidatesKernel(int g1_label, int maxSizeCandidates, int* g2_vertexList, int g1_degree, int* g2_degrees, int* T2_out, 
                                    int* mapping2, int* candidates, int* candidateSize, int g2_numVertices, int* commonNodes, int* g2_matrix, 
                                    int* g2_nodesToLabel, int offset, int useConstMem, int* coveredNeighbors, int* mapping1) {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(constMemVar == 0) {  // constMemVar contains coveredNeighborsSize

        if(idx < maxSizeCandidates) {
            int vertex = g2_vertexList[idx];  // g2_labelToNodes[label]

            if(g2_degrees[vertex] == g1_degree && T2_out[vertex] == 1 && mapping2[vertex] == -1) {
                int index = atomicAdd(candidateSize, 1);
                candidates[index] = vertex;
            }
        }
    }
    else {

        if(idx < g2_numVertices) {
           commonNodes[idx] = 1;

            if(useConstMem) {
                
                for (int i = 0; i < constMemVar; i++) {
                    int nbrG1 = constMem[i];        // constMem contains coveredNeighbors with offset 0
                    int mappedG2 = constMem[offset + nbrG1]; // constMem contains mapping1 with offset degrees[node]
                    if (g2_matrix[mappedG2 * g2_numVertices + idx] == 0) {
                        commonNodes[idx] = 0;
                    }
                }

            } else {

                for (int i = 0; i < constMemVar; i++) {
                    int nbrG1 = coveredNeighbors[i];
                    int mappedG2 = mapping1[nbrG1];
                    if (g2_matrix[mappedG2 * g2_numVertices + idx] == 0) {
                        commonNodes[idx] = 0;
                    }
                }

            }

            if (commonNodes[idx] && mapping2[idx] != -1) {
                commonNodes[idx] = 0;
            }

            if (commonNodes[idx] && g2_degrees[idx] != g1_degree) {
                commonNodes[idx] = 0;
            }

            if (commonNodes[idx] && g2_nodesToLabel[idx] != g1_label) {
                commonNodes[idx] = 0;
            }

            if (commonNodes[idx] == 1) {
                int index = atomicAdd(candidateSize, 1);
                candidates[index] = idx;
            }
        }
    }
}

int* findCandidatesGPU(Graph* d_g1, Graph* d_g2, State* d_state, int node, int* sizeCandidates, Graph* h_g1, Graph* h_g2, cudaStream_t* streams) {
    cudaStream_t stream1 = streams[0];
    cudaStream_t stream2 = streams[1];

    size_t neighSize = h_g1->degrees[node] * sizeof(int);
    size_t size2 = d_g2->numVertices * sizeof(int);
    size_t size1 = d_g1->numVertices * sizeof(int);

    int useConstMem = d_g1->numVertices + h_g1->degrees[node] < MAXCONSTMEM;

    // writes on constMem must be sequential
    if(useConstMem) {
        CUDA_CHECK_ERROR(cudaMemcpyToSymbolAsync(constMem, d_state->mapping1, size1, neighSize, cudaMemcpyDeviceToDevice, stream2));    // neighSize is the offset
    }
    
    int* d_coveredNeighbors, *d_coveredNeighborsSize;
    
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_coveredNeighbors, neighSize));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_coveredNeighborsSize, sizeof(int)));

    CUDA_CHECK_ERROR(cudaMemsetAsync(d_coveredNeighborsSize, 0, sizeof(int), stream1));
    
    int gridSize = (d_g1->numVertices + blockSize - 1) / blockSize;

    findCoveredNeighborsKernel<<<gridSize, blockSize, 0, stream1>>>(d_g1->matrix, d_state->mapping1, node, d_coveredNeighbors, 
                                d_coveredNeighborsSize, d_g1->numVertices);

    int g1_label = h_g1->nodesToLabel[node];
    int g1_degree = h_g1->degrees[node];
    int maxSizeCandidates = h_g2->labelsCardinalities[g1_label];
    
    int* d_candidates, *d_candidateSize, *d_commonNodes;

    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_candidates, size2));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_candidateSize, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_commonNodes, size2));

    CUDA_CHECK_ERROR(cudaMemsetAsync(d_candidateSize, 0, sizeof(int), stream2));
    
    gridSize = (d_g2->numVertices + blockSize - 1) / blockSize;
    
    cudaStreamSynchronize(stream2); // transfer of constMem in the stream2 must be finished before the second write on constMem

    if(useConstMem) {
        CUDA_CHECK_ERROR(cudaMemcpyToSymbolAsync(constMem, d_coveredNeighbors, neighSize, 0, cudaMemcpyDeviceToDevice, stream1));   // it is queued before the kernel call on stream1 so implicit synchronization
    }
    
    CUDA_CHECK_ERROR(cudaMemcpyToSymbolAsync(constMemVar, d_coveredNeighborsSize, sizeof(int), 0, cudaMemcpyDeviceToDevice, stream1)); // constMemVar contains coveredNeighborsSize

    findCandidatesKernel<<<gridSize, blockSize, 0, stream1>>>(g1_label, maxSizeCandidates, d_g2->labelToNodes[g1_label], g1_degree, 
                            d_g2->degrees, d_state->T2_out, d_state->mapping2, d_candidates, d_candidateSize, d_g2->numVertices, 
                            d_commonNodes, d_g2->matrix, d_g2->nodesToLabel, h_g1->degrees[node], useConstMem, );

    int* candidates = (int*)malloc(maxSizeCandidates * sizeof(int));
    
    cudaStreamSynchronize(stream1);

    CUDA_CHECK_ERROR(cudaMemcpyAsync(candidates, d_candidates, maxSizeCandidates * sizeof(int), cudaMemcpyDeviceToHost, stream2));
    CUDA_CHECK_ERROR(cudaMemcpyAsync(sizeCandidates, d_candidateSize, sizeof(int), cudaMemcpyDeviceToHost, stream1));

    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);

    cudaFree(d_coveredNeighbors);
    cudaFree(d_coveredNeighborsSize);

    cudaFree(d_candidates);
    cudaFree(d_candidateSize);
    cudaFree(d_commonNodes);

    return candidates;
}


void vf2ppGPU(Graph* d_g1, Graph* d_g2, State* d_state, Graph* h_g1, Graph* h_g2, cudaStream_t* streams) {
    if (!checkGraphPropertiesGPU(h_g1, h_g2, d_g1, d_g2, streams)) {
        return;
    }

    int* order = orderingGPU(d_g1, d_g2, streams, h_g1);

     //printf("Order:\t");
     //for(int i = 0; i < h_g1->numVertices; i++) {
     //   printf("%d ", order[i]);
     //}
     //printf("\n");

    int sizeCandidates = 0;
    int* candidates = findCandidatesGPU(d_g1, d_g2, d_state, order[0], &sizeCandidates, h_g1, h_g2, streams);

    StackNode* stack = createStack();
    Info* info = createInfo(candidates, sizeCandidates, order[0]);
    push(&stack, info);

    int matchingNode = 1;
    while (!isStackEmpty(stack)) {
        Info* info = peek(&stack);
        bool isMatch = false;

        // printInfo(info);

        for(int i = info->candidateIndex; i < info->sizeCandidates; i++) {
            int candidate = info->candidates[i];
            info->candidateIndex = i + 1;

            int ret = cutISOGPU(d_g1, d_g2, d_state, info->vertex, candidate, h_g1, h_g2, streams);

            // printf("CutISO: %d\n", ret);
            // printf("\n");

            if(!ret) {
                // printf("\nMatch %d -> %d\n", info->vertex, candidate);

                updateStateGPU(d_g1, d_g2, d_state, info->vertex, candidate, streams);

                if(matchingNode >= d_g1->numVertices) {
                    freeStack(stack);
                    free(order);
                    printf("Graphs are isomorphic\n");
                    cudaStreamSynchronize(streams[0]);  // wait for updates on the state
                    cudaStreamSynchronize(streams[1]);
                    return;
                }

                cudaStreamSynchronize(streams[0]);
                cudaStreamSynchronize(streams[1]);

                candidates = findCandidatesGPU(d_g1, d_g2, d_state, order[matchingNode], &sizeCandidates, h_g1, h_g2, streams);
                Info* info = createInfo(candidates, sizeCandidates, order[matchingNode]);
                push(&stack, info);
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
                int candidate = prevInfo->candidates[prevInfo->candidateIndex - 1];
                restoreStateGPU(d_g1, d_g2, d_state, prevInfo->vertex, candidate, streams);
            }
        }
    }
    free(order);
    freeStack(stack);   
}

__global__ void findNeighborsKernel(int* matrix, int node, int* neighbors, int* size, int numVertices) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx < numVertices) {
        if(matrix[node * numVertices + idx] == 1) {
            int index = atomicAdd(size, 1);
            neighbors[index] = idx;
        }
    }
}

__global__ void checkLabelsKernel(int* neighbors1, int nbrSize1, int nbrSize2, int* labelsNbr, int numVertices,
    int* g1_nodesToLabel, int* d_result, int offset, int usage, int* neighbors2, int* g2_nodesToLabel) {   // idx is the id of the thread in the grid

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx < nbrSize1) {
        int nbr1 = neighbors1[idx];
        int labelNbr1 = g1_nodesToLabel[nbr1];
        bool found = false;

        if(usage) {

            for(int i = 0; i < nbrSize2; i++) {
                int nbr2 = constMem[i];     // constMem contains the neighbors2 with offset 0
                if(labelNbr1 == constMem[offset + nbr2]) {  // constMem contains the g2_nodesToLabel with offset degrees[node]
                    found = true;
                    labelsNbr[labelNbr1] = 1;
                    break;
                }
            }

        } else {

            for(int i = 0; i < *nbrSize2; i++) {
                int nbr2 = neighbors2[i];
                if(labelNbr1 == g2_nodesToLabel[nbr2]) {
                    found = true;
                    labelsNbr[labelNbr1] = 1;
                    break;
                }
            }
        }

        if(!found) {
            atomicExch(d_result, 0);   // d_result is initialized to 1 by default
        }
    }
}

__global__ void findNodesOfLabelKernel(int* neighbors, int* g_nodesToLabel, int label, int maxSize, int* size, int* nodes) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx < maxSize) {
        int vertex = neighbors[idx];

        if(g_nodesToLabel[vertex] == label) {
            int index = atomicAdd(size, 1);   // atomicAdd returns the index respect to global memory (so can't be used shared memory)
            nodes[index] = vertex;
        }
    }
}

__global__ void intersectionCountKernel(int* nodes, int* size, int* stateSet, int* count) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;    // idx is the id of the thread in the grid among all threads
    __shared__ int localCount;  // each block has own localCount variable (shared memory)

    if(threadIdx.x == 0) {   // only one thread in the block initializes the localCount
        localCount = 0;
    }

    __syncthreads();

    if(idx < *size) {
        int vertex = nodes[idx];

        if(stateSet[vertex] == 1) {
            atomicAdd(&localCount, 1);
        }
    }

    __syncthreads();

    if(threadIdx.x == 0) {  // only the thread with id 0 of each block updates the global size
        atomicAdd(count, localCount);
    }
}

bool cutISOGPU(Graph* d_g1, Graph* d_g2, State* d_state, int node1, int node2, Graph* h_g1, Graph* h_g2, cudaStream_t* streams) {
    cudaStream_t stream1 = streams[0];
    cudaStream_t stream2 = streams[1];
    cudaStream_t stream3 = streams[2];
    cudaStream_t stream4 = streams[3];
    cudaStream_t stream5 = streams[4];
    cudaStream_t stream6 = streams[5];

    int nbrSize1 = h_g1->degrees[node1];
    int nbrSize2 = h_g2->degrees[node2];

    size_t nbrSize1_bytes = nbrSize1 * sizeof(int);
    size_t nbrSize2_bytes = nbrSize2 * sizeof(int);
    size_t size1 = d_g1->numVertices * sizeof(int);
    size_t size2 = d_g2->numVertices * sizeof(int);
    size_t labelSize = LABELS * sizeof(int);

    int useConstMem = d_g2->numVertices + nbrSize2 < MAXCONSTMEM;

    if(useConstMem) {
        CUDA_CHECK_ERROR(cudaMemcpyToSymbolAsync(constMem, d_g2->nodesToLabel, size2, nbrSize2_bytes, cudaMemcpyDeviceToDevice, stream3)); // constMem contains g2_nodesToLabel
    }

    int* d_neighbors1, *d_neighbors2;
    int* d_nbrSize1, *d_nbrSize2;

    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_neighbors2, nbrSize2_bytes));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_nbrSize2, sizeof(int)));

    CUDA_CHECK_ERROR(cudaMemsetAsync(d_nbrSize2, 0, sizeof(int), stream2));
   
    int gridSize = (d_g2->numVertices + blockSize - 1) / blockSize;

    findNeighborsKernel<<<gridSize, blockSize, 0, stream2>>>(d_g2->matrix, node2, d_neighbors2, d_nbrSize2, d_g2->numVertices);
    
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_neighbors1, nbrSize1_bytes));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_nbrSize1, sizeof(int)));
    
    CUDA_CHECK_ERROR(cudaMemsetAsync(d_nbrSize1, 0, sizeof(int), stream1));
    
    gridSize = (d_g1->numVertices + blockSize - 1) / blockSize;

    findNeighborsKernel<<<gridSize, blockSize, 0, stream1>>>(d_g1->matrix, node1, d_neighbors1, d_nbrSize1, d_g1->numVertices);

    int* d_labelsNbr;
    int* d_result;
    int result = 1;                               

    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_labelsNbr, labelSize));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_result, sizeof(int)));
    
    CUDA_CHECK_ERROR(cudaMemsetAsync(d_labelsNbr, 0, labelSize, stream4));
    CUDA_CHECK_ERROR(cudaMemcpyAsync(d_result, &result, sizeof(int), cudaMemcpyHostToDevice, stream4));

    gridSize = (nbrSize1 + blockSize - 1) / blockSize;

    cudaStreamSynchronize(stream3); // the first write on constMem must be finished before the second write on constMem
    
    if(useConstMem) {
        CUDA_CHECK_ERROR(cudaMemcpyToSymbolAsync(constMem, d_neighbors2, nbrSize2_bytes, 0, cudaMemcpyDeviceToDevice, stream2)); // constMem contains neighbors2. it is queued after findNeighborsKernel on stream2 
    }

    cudaStreamSynchronize(stream4);
    cudaStreamSynchronize(stream2);

    checkLabelsKernel<<<gridSize, blockSize, 0, stream1>>>(d_neighbors1, nbrSize1, nbrSize2, d_labelsNbr, d_g1->numVertices, 
                        d_g1->nodesToLabel, d_result, nbrSize2, useConstMem, d_neighbors2, d_g2->nodesToLabel);

    int* labelsNbr = (int*)malloc(labelSize);
    CUDA_CHECK_ERROR(cudaMemcpyAsync(&result, d_result, sizeof(int), cudaMemcpyDeviceToHost, stream1));
    CUDA_CHECK_ERROR(cudaMemcpyAsync(labelsNbr, d_labelsNbr, labelSize, cudaMemcpyDeviceToHost, stream1));

    cudaStreamSynchronize(stream1);

    if(result == 0) {
        cudaFree(d_neighbors1);
        cudaFree(d_neighbors2);
        cudaFree(d_nbrSize1);
        cudaFree(d_nbrSize2);
        cudaFree(d_labelsNbr);
        cudaFree(d_result);
        free(labelsNbr);
        return true;
    }

    // pinned memory for faster host-device communication
    int *h_size1, *h_size2;
    int *h_count1, *h_count2, *h_count3, *h_count4;
    int *d_count1, *d_count2, *d_count3, *d_count4;

    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_size1, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_size2, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_count1, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_count2, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_count3, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMallocHost((void**)&h_count4, sizeof(int)));

    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_count1, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_count2, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_count3, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_count4, sizeof(int)));

    int *d_nodes_g1, *d_size_g1;
    int *d_nodes_g2, *d_size_g2;

    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_nodes_g1, nbrSize1_bytes));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_size_g1, sizeof(int)));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_nodes_g2, nbrSize2_bytes));
    CUDA_CHECK_ERROR(cudaMalloc((void**)&d_size_g2, sizeof(int)));

    int gridSize1 = (nbrSize1 + blockSize - 1) / blockSize;
    int gridSize2 = (nbrSize2 + blockSize - 1) / blockSize;
   
    bool ret = false;
    for(int label = 0; label < LABELS; label++) {
        if(labelsNbr[label] == 1) {

            CUDA_CHECK_ERROR(cudaMemsetAsync(d_size_g1, 0, sizeof(int), stream1));
            CUDA_CHECK_ERROR(cudaMemsetAsync(d_size_g2, 0, sizeof(int), stream2));

            findNodesOfLabelKernel<<<gridSize1, blockSize, 0, stream1>>>(d_neighbors1, d_g1->nodesToLabel, label, nbrSize1, d_size_g1, d_nodes_g1);
            findNodesOfLabelKernel<<<gridSize2, blockSize, 0, stream2>>>(d_neighbors2, d_g2->nodesToLabel, label, nbrSize2, d_size_g2, d_nodes_g2);

            CUDA_CHECK_ERROR(cudaMemcpyAsync(h_size1, d_size_g1, sizeof(int), cudaMemcpyDeviceToHost, stream1));
            CUDA_CHECK_ERROR(cudaMemcpyAsync(h_size2, d_size_g2, sizeof(int), cudaMemcpyDeviceToHost, stream2));

            CUDA_CHECK_ERROR(cudaMemsetAsync(d_count1, 0, sizeof(int), stream3));
            CUDA_CHECK_ERROR(cudaMemsetAsync(d_count2, 0, sizeof(int), stream4));
            CUDA_CHECK_ERROR(cudaMemsetAsync(d_count3, 0, sizeof(int), stream5));
            CUDA_CHECK_ERROR(cudaMemsetAsync(d_count4, 0, sizeof(int), stream6));

            cudaStreamSynchronize(stream1);
            cudaStreamSynchronize(stream2);

            gridSize = (*h_size1 + blockSize - 1) / blockSize;
            intersectionCountKernel<<<gridSize, blockSize, 0, stream3>>>(d_nodes_g1, d_size_g1, d_state->T1, d_count1);
            intersectionCountKernel<<<gridSize, blockSize, 0, stream5>>>(d_nodes_g1, d_size_g1, d_state->T1_out, d_count3);

            CUDA_CHECK_ERROR(cudaMemcpyAsync(h_count1, d_count1, sizeof(int), cudaMemcpyDeviceToHost, stream3));
            CUDA_CHECK_ERROR(cudaMemcpyAsync(h_count3, d_count3, sizeof(int), cudaMemcpyDeviceToHost, stream5));
            
            gridSize = (*h_size2 + blockSize - 1) / blockSize;
            intersectionCountKernel<<<gridSize, blockSize, 0, stream4>>>(d_nodes_g2, d_size_g2, d_state->T2, d_count2);
            intersectionCountKernel<<<gridSize, blockSize, 0, stream6>>>(d_nodes_g2, d_size_g2, d_state->T2_out, d_count4);

            CUDA_CHECK_ERROR(cudaMemcpyAsync(h_count2, d_count2, sizeof(int), cudaMemcpyDeviceToHost, stream4));
            CUDA_CHECK_ERROR(cudaMemcpyAsync(h_count4, d_count4, sizeof(int), cudaMemcpyDeviceToHost, stream6));
           
            cudaStreamSynchronize(stream3);
            cudaStreamSynchronize(stream4);
            cudaStreamSynchronize(stream5);
            cudaStreamSynchronize(stream6);

            if(*h_count1 != *h_count2 || *h_count3 != *h_count4) {
                ret = true;
                break;
            }
        }
    }

    cudaFree(d_neighbors1);
    cudaFree(d_neighbors2);
    cudaFree(d_nbrSize1);
    cudaFree(d_nbrSize2);
    cudaFree(d_labelsNbr);
    cudaFree(d_result);
    free(labelsNbr);

    cudaFreeHost(h_size1);
    cudaFreeHost(h_size2);
    cudaFreeHost(h_count1);
    cudaFreeHost(h_count2);
    cudaFreeHost(h_count3);
    cudaFreeHost(h_count4);
    
    cudaFree(d_count1);
    cudaFree(d_count2);
    cudaFree(d_count3);
    cudaFree(d_count4);

    cudaFree(d_nodes_g1);
    cudaFree(d_size_g1);
    cudaFree(d_nodes_g2);
    cudaFree(d_size_g2);

    return ret;
}

/***** STACK FUNCTIONS *****/
Info* createInfo(int* candidates, int sizeCandidates, int vertex) {
    Info* info = (Info*)malloc(sizeof(Info));
    if (info == NULL) {
        printf("Error allocating memory in createInfo\n");
        exit(EXIT_FAILURE);
    }
    info->vertex = vertex;
    info->candidates = (int*)realloc(candidates, sizeCandidates * sizeof(int));
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

void push(StackNode** top, Info* info) {
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
    printf("\nVertex: %d\n", info->vertex);
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