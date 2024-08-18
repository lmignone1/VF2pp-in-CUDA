/*
* MIGNONE LORENZO 0622701866 L.MIGNONE@STUDENTI.UNISA.IT
* Course: High Performance Computing 2022/2023
* Lecturer: Francesco Moscato	fmoscato@unisa.it
* 
* Copyright (C) 2024 - All Rights Reserved
* 
* This file is part of VF2pp-in-CUDA.
* 
* VF2pp-in-CUDA is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* VF2pp-in-CUDA is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with VF2pp-in-CUDA.  If not, see <http://www.gnu.org/licenses/>.
*/

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
