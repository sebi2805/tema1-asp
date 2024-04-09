#include "helpers.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void initializeMPI(int *argc, char ***argv, int *rank, int *size) {
    MPI_Init(argc, argv);  
    MPI_Comm_rank(MPI_COMM_WORLD, rank); 
    MPI_Comm_size(MPI_COMM_WORLD, size); 
}

void finalizeMPI() {
    MPI_Finalize();  
}

void readVectors(const char *xFile, const char *yFile, double *local_x, double *local_y, int N, int rank, int size) {
    int local_N = N / size; // Asumăm că N este divizibil exact cu numărul de procese pentru simplitate
    MPI_File fh_x, fh_y;
    MPI_Status status;

    MPI_Offset offset_x = rank * local_N * sizeof(double);
    MPI_Offset offset_y = rank * local_N * sizeof(double);

    MPI_File_open(MPI_COMM_SELF, xFile, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_x);
    MPI_File_open(MPI_COMM_SELF, yFile, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_y);

    MPI_File_read_at(fh_x, offset_x, local_x, local_N, MPI_DOUBLE, &status);
    MPI_File_read_at(fh_y, offset_y, local_y, local_N, MPI_DOUBLE, &status);

    MPI_File_close(&fh_x);
    MPI_File_close(&fh_y);
}

void readMatrix(const char *matrixFile, double **local_A, int N, int rank, int size) {
    int local_N = N / size; // Numărul de rânduri din matrice alocat fiecărui proces
    MPI_File fh;
    MPI_Status status;

    MPI_Offset offset = rank * local_N * N * sizeof(double); // fiecare rând are N elemente

    *local_A = (double *)malloc(local_N * N * sizeof(double));
    if (*local_A == NULL) {
        fprintf(stderr, "Eroare la alocarea memoriei pentru blocul local al matricei\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    MPI_File_open(MPI_COMM_SELF, matrixFile, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_at(fh, offset, *local_A, local_N * N, MPI_DOUBLE, &status);
    MPI_File_close(&fh);
}

void distributeData(double *local_x, double *local_y, double **local_A, const double *x, const double *y, const double *A, int N, int local_N, int rank, int size) {
    // int local_N = N / MPI_Comm_size(group_comm, &size); // Presupunem divizibilitate exactă pentru simplitate
    
    if (rank == 0) {
        double *x = (double *)malloc(N * sizeof(double));
        double *y = (double *)malloc(N * sizeof(double));
        double *A = (double *)malloc(N * N * sizeof(double));
        readVectors("x.dat", "y.dat", x, y, N, rank, size);
        readMatrix("mat.dat", A, N, rank, size);

        // MPI_Scatter(x, local_N, MPI_DOUBLE, local_x, local_N, MPI_DOUBLE, 0, group_comm);
        MPI_Scatter(y, local_N, MPI_DOUBLE, local_y, local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(y, local_N, MPI_DOUBLE, local_y, local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (int i = 0; i < size; i++) {
            if (i == 0) {
                memcpy(*local_A, A, local_N * N * sizeof(double));
            } else {
                MPI_Send(&A[i * local_N * N], local_N * N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }
        }

        free(x);
        free(y);
        free(A);
    } else {
        MPI_Scatter(NULL, local_N, MPI_DOUBLE, local_x, local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(NULL, local_N, MPI_DOUBLE, local_y, local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Recv(*local_A, local_N * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

double calculateNumerator(const double *local_x, const double **local_A, const double *local_y, int local_N, int N) {
    double sum = 0.0;
    for (int i = 0; i < local_N; i++) {
        for (int j = 0; j < N; j++) {
            sum += local_x[i] * local_A[i][j] * local_y[j];
        }
    }
    return sum;
}


double calculateDenominator(const double *local_x, const double *local_y, int local_N) {
    double sum = 0.0;
    for (int i = 0; i < local_N; i++) {
        sum += local_x[i] * local_y[i];
    }
    return sum;
}

double aggregateNumerator(double partialNumerator, int rank, MPI_Comm commNumerator) {
    double totalNumerator = 0.0;
    MPI_Reduce(&partialNumerator, &totalNumerator, 1, MPI_DOUBLE, MPI_SUM, 0, commNumerator);
    return (rank == 0) ? totalNumerator : 0;
}

double aggregateDenominator(double partialDenominator, int rank, MPI_Comm commDenominator) {
    double totalDenominator = 0.0;
    MPI_Reduce(&partialDenominator, &totalDenominator, 1, MPI_DOUBLE, MPI_SUM, 0, commDenominator);
    return (rank == 0) ? totalDenominator : 0;
}

