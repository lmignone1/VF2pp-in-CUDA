#ifndef STATE_H
#define STATE_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "graph.h"

typedef struct {
    int *mapping1;  // mapping from query to target
    int *mapping2;  // mapping from target to query
    int *T1;        // Ti contains uncovered neighbors of covered nodes from Gi, i.e. nodes that are not in the mapping, but are neighbors of nodes that are.
    int *T2;        
    int* T1_out;     //Ti_out contains all the nodes from Gi, that are neither in the mapping nor in Ti. Cioe nodi che non sono in mapping e sono vicini di nodi non coperti
    int* T2_out;
} State;

State* createState(Graph*, Graph*);
void freeState(State*);
void printState(State*, int);
bool isMappingFull(Graph*, State*);
void updateState(Graph*, Graph*, State*, int, int);
void restoreState(Graph*, Graph*, State*, int, int);

#ifdef __cplusplus
}
#endif

#endif
