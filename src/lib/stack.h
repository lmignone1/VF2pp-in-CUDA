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
