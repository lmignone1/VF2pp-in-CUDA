#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "graph_list.h"

#define FILENAME_QUERY "../data/graph_query_example.csv"
#define FILENAME_TARGET "../data/graph_target_example.csv"

int main() {
    Node* node = createNode(1);
    printf("Node: %d\n", node->vertex);
    return EXIT_SUCCESS;
}