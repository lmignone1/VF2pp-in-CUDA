#include <stdio.h>
#include <cuda_runtime.h>

// Definizione della struttura dati state
typedef struct {
    float *vector1;
    float *vector2;
    float *vector3;
    int** vector4;
} state;

// Funzione per controllare eventuali errori CUDA
void checkCudaError(cudaError_t err, const char *msg) {
    if (err != cudaSuccess) {
        fprintf(stderr, "Errore CUDA: %s: %s\n", msg, cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

// Kernel CUDA per riempire i vettori
__global__ void fillVectors(float *vector1, float *vector2, float *vector3, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        vector1[idx] = idx * 1.0f;
        vector2[idx] = idx * 2.0f;
        vector3[idx] = idx * 3.0f;
    }
}

__global__ void fillVector4(int* v, int N, int value) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        v[idx] = value;
    }
}

int main() {
    // Dimensione dei vettori
    int N = 100;
    int size = N * sizeof(float);

    // Allocazione dinamica della struttura state
    state *s;
    s = (state *)malloc(sizeof(state));
    if (s == NULL) {
        fprintf(stderr, "Errore di allocazione della memoria per la struttura state\n");
        return EXIT_FAILURE;
    }

    // Allocazione dei vettori su device (GPU)
    checkCudaError(cudaMalloc((void**)&(s->vector1), size), "Allocazione vector1");
    checkCudaError(cudaMalloc((void**)&(s->vector2), size), "Allocazione vector2");
    checkCudaError(cudaMalloc((void**)&(s->vector3), size), "Allocazione vector3");

    s->vector4 = (int**)malloc(3 * sizeof(int*));
    for (int i = 0; i < 3; i++) {
        checkCudaError(cudaMalloc((void**)&(s->vector4[i]), 10*sizeof(int)), "Allocazione vector4");
    }

    // Assicurarsi che l'allocazione sia avvenuta correttamente
    if (s->vector1 == NULL || s->vector2 == NULL || s->vector3 == NULL) {
        fprintf(stderr, "Errore di allocazione della memoria per i vettori\n");
        return EXIT_FAILURE;
    }

    for(int i = 0; i < 3; i++) {
        int* v = s->vector4[i];
        cudaMemset(v, 0, 10*sizeof(int));
    }

    // Lanciare il kernel per riempire i vettori
    int threadsPerBlock = 256;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    fillVectors<<<blocksPerGrid, threadsPerBlock>>>(s->vector1, s->vector2, s->vector3, N);
    checkCudaError(cudaGetLastError(), "Lancio del kernel fillVectors");

    for(int i = 0; i < 3; i++) {
        int* v = s->vector4[i];
        fillVector4<<<(10 + threadsPerBlock - 1) / threadsPerBlock, threadsPerBlock>>>(v, 10, i+1);
    }
   
    checkCudaError(cudaDeviceSynchronize(), "Sincronizzazione del dispositivo");

    // Copia i dati dai vettori di dispositivo alla memoria host per la verifica
    float *h_vector1 = (float *)malloc(size);
    float *h_vector2 = (float *)malloc(size);
    float *h_vector3 = (float *)malloc(size);
    int* h_vector4 = (int*)malloc(10*sizeof(int));

    if (h_vector1 == NULL || h_vector2 == NULL || h_vector3 == NULL) {
        fprintf(stderr, "Errore di allocazione della memoria per i vettori host\n");
        return EXIT_FAILURE;
    }

    checkCudaError(cudaMemcpy(h_vector1, s->vector1, size, cudaMemcpyDeviceToHost), "Copia vector1 da device a host");
    checkCudaError(cudaMemcpy(h_vector2, s->vector2, size, cudaMemcpyDeviceToHost), "Copia vector2 da device a host");
    checkCudaError(cudaMemcpy(h_vector3, s->vector3, size, cudaMemcpyDeviceToHost), "Copia vector3 da device a host");

    for(int i = 0; i < 3; i++) {
        int* v = s->vector4[i];
        checkCudaError(cudaMemcpy(h_vector4, v, 10*sizeof(int), cudaMemcpyDeviceToHost), "Copia vector4 da device a host");
        printf("Primi 10 elementi di vector4:\n");
        for (int i = 0; i < 10; i++) {
            printf("%d ", h_vector4[i]);
        }
        printf("\n");
    }

    // Stampa i primi 10 elementi di ogni vettore per la verifica
    printf("Primi 10 elementi di vector1:\n");
    for (int i = 0; i < 10; i++) {
        printf("%f ", h_vector1[i]);
    }
    printf("\n");

    printf("Primi 10 elementi di vector2:\n");
    for (int i = 0; i < 10; i++) {
        printf("%f ", h_vector2[i]);
    }
    printf("\n");

    printf("Primi 10 elementi di vector3:\n");
    for (int i = 0; i < 10; i++) {
        printf("%f ", h_vector3[i]);
    }
    printf("\n");

    // Deallocazione della memoria su device
    checkCudaError(cudaFree(s->vector1), "Deallocazione vector1");
    checkCudaError(cudaFree(s->vector2), "Deallocazione vector2");
    checkCudaError(cudaFree(s->vector3), "Deallocazione vector3");

    // Deallocazione della memoria della struttura state
    free(s);
    // Deallocazione della memoria host
    free(h_vector1);
    free(h_vector2);
    free(h_vector3);

    return 0;
}
