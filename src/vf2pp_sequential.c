#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "lib/stack.h"
#include "lib/queue.h"
#include "lib/graph.h"
#include "lib/state.h"

#define PATH_QUERY "data/graph_query_%d.csv"
#define PATH_TARGET "data/graph_target_%d.csv"
#define INF 99999

void vf2pp(Graph*, Graph*, State*);
bool checkGraphProperties(Graph*, Graph*);
bool checkSequenceDegree(int*, int*, int);
int compare(const void*, const void*);
int* ordering(Graph*, Graph*);
int* bfs(Graph*, int, int*);
int findLevelNodes(Graph*, int*, int*, int);
void processDepth(int*, int*, Graph*, int*, int*, int*, int*, int);
int* findCoveredNeighbors(Graph*, State*, int, int*);
int* findCandidates(Graph*, Graph*, State*, int, int*);
int intersectionCount(int*, int, int*);
int* findNodesOfLabel(int*, Graph*, int, int, int*);
bool cutISO(Graph*, Graph*, State*, int, int);
int* findNeighbors(Graph*, int, int*);
bool checkLabels(Graph*, Graph*, int*, int*, int, int, int*);

int main(int argc, char* argv[]) {
    int V = atoi(argv[1]); // number of vertices
    int D = atoi(argv[2]);  // degree
    
    char path1[256];
    char path2[256];
    sprintf(path1, PATH_QUERY, V);
    sprintf(path2, PATH_TARGET, V);

    Graph* g1 = readGraph(path1);
    Graph* g2 = readGraph(path2);
    State* s = createState(g1, g2);

    // printf("\n");
    // printState(s, g1->numVertices);
    // printf("\n");
    
    vf2pp(g1, g2, s);
    
    printf("Mapping\n");
    for(int i = 0; i < g1->numVertices; i++) {
        printf("%d -> %d\n", i, s->mapping1[i]);
    }

    freeGraph(g1);
    freeGraph(g2);
    freeState(s);

    return EXIT_SUCCESS;
}

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

bool checkSequenceDegree(int* degree1, int* degree2, int size) {
    int* tmp1 = (int*)malloc(size * sizeof(int));
    int* tmp2 = (int*)malloc(size * sizeof(int));

    memcpy(tmp1, degree1, size * sizeof(int));
    memcpy(tmp2, degree2, size * sizeof(int));

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
    int* labelRarity = (int*)malloc(LABELS * sizeof(int)); 
    int* connectivityG1 = (int*)malloc(g1->numVertices * sizeof(int));   // number of neighbors already ordered
    int* V1Unordered = (int*)malloc(g1->numVertices * sizeof(int));    // V1Unordered[i] = -1 if node has not been ordered yet else 1 

    if (order == NULL || labelRarity == NULL || connectivityG1 == NULL || V1Unordered == NULL) {
        printf("Error allocating memory in ordering\n");
        exit(EXIT_FAILURE);
    }

    memcpy(labelRarity, g1->labelsCardinalities, LABELS * sizeof(int)); 

    for(int i = 0; i < g1->numVertices; i++) {
        connectivityG1[i] = 0;
        V1Unordered[i] = -1;
    }

    int maxRarity = INF;
    int maxNode = -1;
    for (int vertex = 0; vertex < g1->numVertices; vertex++) {
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
    
    int maxDepth = 0;
    int* levels = bfs(g1, maxNode, &maxDepth);
    int* levelNodes = (int*)malloc((g1->numVertices) * sizeof(int));
    int order_index = 0;

    for (int depth = 0; depth <= maxDepth; depth++) {
        int levelSize = findLevelNodes(g1, levels, levelNodes, depth);
        processDepth(order, &order_index, g1, connectivityG1, labelRarity, V1Unordered, levelNodes, levelSize);
    }
    free(levels);
    free(levelNodes);
            
    free(labelRarity);
    free(connectivityG1);
    free(V1Unordered);

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

void processDepth(int* order, int* order_index, Graph* g, int* connectivityG1, int* labelRarity, int* V1Unordered, int* levelNodes, int levelSize) {
    int l = levelSize;
    
    while (levelSize > 0) {
        int maxConnectivity = -INF;  
        int maxDegree = -INF;
        int maxRarity = INF;
        int nextNode = -1;
    
        for(int i = 0; i < l; i++) {
            int vertex = levelNodes[i];

            if (V1Unordered[vertex] == 1) {
                continue;
            }

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
        
        levelSize--;
        labelRarity[g->nodesToLabel[nextNode]]--;
        V1Unordered[nextNode] = 1;
    }
}

int* findCandidates(Graph* g1, Graph* g2, State* state, int node, int* sizeCandidates) {  
    int coveredNeighborsSize = 0;
    int *coveredNeighbors = findCoveredNeighbors(g1, state, node, &coveredNeighborsSize);

    // printf("CoveredNeighbors with size %d of node %d: ", coveredNeighborsSize, node);
    // for(int i = 0; i < coveredNeighborsSize; i++) {
    //     printf("%d ", coveredNeighbors[i]);
    // }
    // printf("\n");

    if(coveredNeighborsSize == 0) {
        int label = g1->nodesToLabel[node];
        int maxSizeCandidates = g2->labelsCardinalities[label];
        int* candidates = malloc(maxSizeCandidates * sizeof(int));

        *sizeCandidates = 0;
        for(int i = 0; i < maxSizeCandidates; i++) {
            int vertex = g2->labelToNodes[label][i];  // g2->labelToNodes[label] is a list of candidates
            if(g2->degrees[vertex] == g1->degrees[node] && state->T2_out[vertex] == 1 && state->mapping2[vertex] == -1) {
                candidates[*sizeCandidates] = vertex;
                *sizeCandidates = *sizeCandidates + 1;
            }
        }
        free(coveredNeighbors);
    
        // printf("\nCandidates for node %d: ", node);
        // for(int i = 0; i < *sizeCandidates; i++) {
        //     printf("%d ", candidates[i]);
        // }
        // printf("\n");

        return candidates;
    }
    
    int* commonNodes = (int*)malloc(g2->numVertices * sizeof(int));
    for(int i = 0; i < g2->numVertices; i++) {
        commonNodes[i] = 1;  
    }
    
    int count = 0;
    for(int i = 0; i < coveredNeighborsSize; i++) {
        int nbrG1 = coveredNeighbors[i];
        int mappedG2 = state->mapping1[nbrG1];

        for(int adjVertex = 0; adjVertex < g2->numVertices; adjVertex++) {
            if(g2->matrix[mappedG2 * g2->numVertices + adjVertex] == 0 && commonNodes[adjVertex] != 0) {    // avoids duplicate counts
                commonNodes[adjVertex] = 0;
                count++;
            }
        }
    }
    
    int commonNodesSize = g2->numVertices - count;

    // printf("first common nodes size is %d and are\n", commonNodesSize);
    // for(int i = 0; i < g2->numVertices; i++) {
    //     printf("%d ", commonNodes[i]);
    // }
    // printf("\n");

    count = 0;
    for(int vertex = 0; vertex < g2->numVertices; vertex++) {
       
        if(commonNodes[vertex] == 1 && state->mapping2[vertex] != -1) {
            commonNodes[vertex] = 0;
            count++;
        }

    }

    commonNodesSize = commonNodesSize - count;

    // printf("second common nodes size is %d and are\n", commonNodesSize);
    // for(int i = 0; i < g2->numVertices; i++) {
    //     printf("%d ", commonNodes[i]);
    // }
    // printf("\n");

    int degree = g1->degrees[node];
    count = 0;
    for(int vertex = 0; vertex < g2->numVertices; vertex++) {
        if(commonNodes[vertex] == 1 && g2->degrees[vertex] != degree) {
            commonNodes[vertex] = 0;
            count++;
        }
    }

    commonNodesSize = commonNodesSize - count;

    // printf("third common nodes size is %d and are\n", commonNodesSize);
    // for(int i = 0; i < g2->numVertices; i++) {
    //     printf("%d ", commonNodes[i]);
    // }
    // printf("\n");

    int label = g1->nodesToLabel[node];
    count = 0;
    for(int vertex = 0; vertex < g2->numVertices; vertex++) {
        if(commonNodes[vertex] == 1 && g2->nodesToLabel[vertex] != label) {
            commonNodes[vertex] = 0;
            count++;
        }
    }

    commonNodesSize = commonNodesSize - count;
    
    // printf("fourth common nodes size is %d and are\n", commonNodesSize);
    // for(int i = 0; i < g2->numVertices; i++) {
    //     printf("%d ", commonNodes[i]);
    // }
    // printf("\n");
    
    int* candidates = (int*)malloc(commonNodesSize * sizeof(int));
    *sizeCandidates = 0;
    for(int i = 0; i < g2->numVertices; i++) {
        if(commonNodes[i] == 1) {
            candidates[*sizeCandidates] = i;
            *sizeCandidates = *sizeCandidates + 1;
        }
    }

    // printf("CommonNodes for node %d: ", node);
    // for(int i = 0; i < *sizeCandidates; i++) {
    //     printf("%d ", candidates[i]);
    // }
    // printf("\n");

    free(coveredNeighbors);
    free(commonNodes);
    return candidates;
}

int* findCoveredNeighbors(Graph* g, State* state, int node, int* size) {
    int* coveredNeighbors = (int*)malloc(g->degrees[node] * sizeof(int));   
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
        return;
    }

    int* order = ordering(g1, g2);
    
    // printf("Order:\t");
    // for(int i = 0; i < g1->numVertices; i++) {
    //     printf("%d ", order[i]);
    // }
    // printf("\n");

    int sizeCandidates = 0;
    int* candidates = findCandidates(g1, g2, state, order[0], &sizeCandidates);
    
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

            int ret = cutISO(g1, g2, state, info->vertex, candidate);

            // printf("CutISO: %d\n", ret);
            // printf("\n");

            if(!ret) {
                
                // printf("\nMatch %d -> %d\n", info->vertex, candidate);
                
                state->mapping1[info->vertex] = candidate;
                state->mapping2[candidate] = info->vertex;

                if(isMappingFull(g1, state)) {
                    freeStack(stack);
                    free(order);
                    printf("Graphs are isomorphic\n");
                    return;
                }

                updateState(g1, g2, state, info->vertex, candidate);
                candidates = findCandidates(g1, g2, state, order[matchingNode], &sizeCandidates);
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
    // printf("\nCutISO node1 %d e node2 %d\n", node1, node2);

    int nbrSize1 = 0, nbrSize2 = 0;
    int* neighbors1 = findNeighbors(g1, node1, &nbrSize1);

    // printf("Neighbors1 of node1 %d: ", node1);
    // for(int i = 0; i < nbrSize1; i++) {
    //     printf("%d ", neighbors1[i]);
    // }
    // printf("\n");

    int* neighbors2 = findNeighbors(g2, node2, &nbrSize2);

    // printf("Neighbors2 of node2 %d: ", node2);
    // for(int i = 0; i < nbrSize2; i++) {
    //     printf("%d ", neighbors2[i]);
    // }
    // printf("\n");

    int* labelsNbr = (int*)malloc(LABELS * sizeof(int));

    if(labelsNbr == NULL) {
        printf("Error allocating memory in cutISO\n");
        exit(EXIT_FAILURE);
    }

    bool result = checkLabels(g1, g2, neighbors1, neighbors2, nbrSize1, nbrSize2, labelsNbr); // check if labels are the same
    
    // printf("CheckLabels: %d\n", result);
    // printf("LabelsNbr: ");
    // for(int i = 0; i < LABELS; i++) {
    //     if(labelsNbr[i] == 1) {
    //         printf("%d ", i);
    //     }
    //     // printf("%d ", labelsNbr[i]);
    // }
    // printf("\n");

    if(!result) {   // if labels are different, cut
        free(neighbors1);
        free(neighbors2);
        free(labelsNbr);
        return true;
    }

    // in this case labels are the same
    // check if the number of neighbors with the same label in G1 and G2 are the same
    for(int label = 0; label < LABELS; label++) {
        if (labelsNbr[label] == 1) {    // if the label effectively the same
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

            if(labelNbr1 == labelNbr2) {
                found = true;
                labelsNbr[labelNbr1] = 1;
                break;
            }
        }

        if (!found) {
            return false;
        }
    }

    return true;
}

int* findNeighbors(Graph* g, int node, int* size) {
    int* neighbors = (int*)malloc(g->degrees[node] * sizeof(int));  
    
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