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

#ifndef STACK_H
#define STACK_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

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

#ifdef __cplusplus
}
#endif

#endif
