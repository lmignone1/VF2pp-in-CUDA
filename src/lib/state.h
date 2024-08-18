
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
    int* T1_out;     //Ti_out contains all the nodes from Gi, that are neither in the mapping nor in Ti, i.e. nodes that are not in the mapping and are neighbors of uncovered nodes
    int* T2_out;
} State;

State* createState(Graph*, Graph*);
void freeState(State*);
void printState(State*, int);
void updateState(Graph*, Graph*, State*, int, int);
void restoreState(Graph*, Graph*, State*, int, int);

#ifdef __cplusplus
}
#endif

#endif
