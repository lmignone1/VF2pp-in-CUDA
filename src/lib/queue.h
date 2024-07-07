#ifndef QUEUE_H
#define QUEUE_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct {
    int* data;
    int head;
    int tail;
    int capacity;
} Queue;

Queue* createQueue(int);
bool isQueueEmpty(Queue*);
void resizeQueue(Queue*);
void enqueue(Queue*, int);
int dequeue(Queue*);
void freeQueue(Queue*);

#ifdef __cplusplus
}
#endif

#endif
