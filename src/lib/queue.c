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

#include "queue.h"

Queue* createQueue(int capacity) {
    Queue* q = (Queue*)malloc(sizeof(Queue));
    q->data = (int*)malloc(capacity * sizeof(int));
    q->head = -1;
    q->tail = -1;
    q->capacity = capacity;

    // if (q == NULL || q->data == NULL) {
    //     printf("Error allocating memory in createQueue\n");
    //     exit(EXIT_FAILURE);
    // }

    return q;
}

bool isQueueEmpty(Queue* q) {
    return q->head == -1;
}

void resizeQueue(Queue* q) {
    q->capacity *= 2;
    q->data = (int*)realloc(q->data, q->capacity * sizeof(int));

    // if (q->data == NULL) {
    //     printf("Error reallocating memory in resizeQueue\n");
    //     exit(EXIT_FAILURE);
    // }
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