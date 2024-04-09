#include "helpers.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
    int rank, size, N = 200;
    double *local_x, *local_y, *x, *y, *A, **local_A;
    MPI_Comm commNumerator, commDenominator;

    // Inițializarea MPI
    initializeMPI(&argc, &argv, &rank, &size);

    // Alocarea memoriei pentru vectorii și matricea completă în procesul 0
    if (rank == 0) {
        x = (double *)malloc(N * sizeof(double));
        y = (double *)malloc(N * sizeof(double));
        A = (double *)malloc(N * N * sizeof(double));
    }

    // Alocarea memoriei pentru segmentele locale ale vectorilor și matricei
    local_x = (double *)malloc(N / size * sizeof(double));
    local_y = (double *)malloc(N / size * sizeof(double));
    local_A = (double **)malloc(N / size * sizeof(double*));
    for (int i = 0; i < N / size; i++) {
        local_A[i] = (double *)malloc(N * sizeof(double));
    }

    // Citirea și distribuția datelor
    if (rank == 0) {
        // Citirea datelor complete
        readVectors("x.dat", "y.dat", x, y, N, rank, size);
        readMatrix("mat.dat", &A, N, rank, size);
    }

    // Distribuția datelor către toate procesele
    distributeData(local_x, local_y, local_A, x, y, A, N, N/size, rank, size);

    // Crearea a două comunicatoare pentru numărător și numitor
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &commNumerator);
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &commDenominator);

    // Calculul parțial al numărătorului și numitorului
    double partialNumerator = calculateNumerator(local_x, (const double **)local_A, local_y, N / size, N);
    double partialDenominator = calculateDenominator(local_x, local_y, N / size);

    // Agregarea rezultatelor parțiale
    double totalNumerator = aggregateNumerator(partialNumerator, rank, commNumerator);
    double totalDenominator = aggregateDenominator(partialDenominator, rank, commDenominator);

    if (rank == 0) {
        double weightedAverage = totalNumerator / totalDenominator;
        printf("Media ponderata calculata: %lf\n", weightedAverage);
    }

    free(local_x);
    free(local_y);
    for (int i = 0; i < N / size; i++) {
        free(local_A[i]);
    }
    free(local_A);
    if (rank == 0) {
        free(x);
        free(y);
        free(A);
    }
    finalizeMPI();

    return 0;
}
