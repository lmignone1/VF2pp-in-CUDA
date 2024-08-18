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

#include "stack.h"

Info* createInfo(int* candidates, int sizeCandidates, int vertex) {
    Info* info = (Info*)malloc(sizeof(Info));
    // if (info == NULL) {
    //     printf("Error allocating memory in createInfo\n");
    //     exit(EXIT_FAILURE);
    // }
    info->vertex = vertex;
    info->candidates = (int*)realloc(candidates, sizeCandidates * sizeof(int));
    info->sizeCandidates = sizeCandidates;
    info->candidateIndex = 0;
    return info;
}

StackNode* createStackNode(Info* info) {
    StackNode* node = (StackNode*)malloc(sizeof(StackNode));
    // if (node == NULL) {
    //     printf("Error allocating memory in createStackNode\n");
    //     exit(EXIT_FAILURE);
    // }
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
    // if (isStackEmpty(*top)) {
    //     printf("Stack is empty, cannot pop\n");
    //     return NULL;
    // }
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