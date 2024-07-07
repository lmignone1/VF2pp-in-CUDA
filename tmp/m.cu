#include <stdio.h>
#include <cuda_runtime.h>
#include "graph_list.h"

// Dichiarazione della funzione esterna
extern "C" Node *createNode(int vertex);
// Kernel CUDA semplice
__global__ void add(int *a, int *b, int *c) {
    int index = threadIdx.x;
    c[index] = a[index] + b[index];
}

int main(void) {
    int a_host[5] = {1, 2, 3, 4, 5};
    int b_host[5] = {10, 20, 30, 40, 50};
    int c_host[5];

    int *a_device, *b_device, *c_device;
    int size = 5 * sizeof(int);

    // Allocazione della memoria sul dispositivo
    cudaMalloc((void **)&a_device, size);
    cudaMalloc((void **)&b_device, size);
    cudaMalloc((void **)&c_device, size);

    // Copia dei dati dalla memoria host alla memoria device
    cudaMemcpy(a_device, a_host, size, cudaMemcpyHostToDevice);
    cudaMemcpy(b_device, b_host, size, cudaMemcpyHostToDevice);

    // Lancio del kernel CUDA
    add<<<1, 5>>>(a_device, b_device, c_device);

    // Copia dei risultati dalla memoria device alla memoria host
    cudaMemcpy(c_host, c_device, size, cudaMemcpyDeviceToHost);

    // Stampa dei risultati
    printf("Risultati dell'addizione nel kernel CUDA:\n");
    for (int i = 0; i < 5; i++) {
        printf("c_host[%d] = %d\n", i, c_host[i]);
    }

    // Utilizzo della funzione esterna
    Node *node = createNode(10);
    printf("Node: %d\n", node->vertex);

    // Libera la memoria allocata sul dispositivo
    cudaFree(a_device);
    cudaFree(b_device);
    cudaFree(c_device);

    return 0;
}
