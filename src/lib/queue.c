#include "queue.h"

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