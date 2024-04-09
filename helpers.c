#include "helpers.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void initializeMPI(int *argc, char ***argv, int *rank, int *size) {
    MPI_Init(argc, argv); // Inițializează mediul MPI
    MPI_Comm_rank(MPI_COMM_WORLD, rank); // Obține rangul procesului curent
    MPI_Comm_size(MPI_COMM_WORLD, size); // Obține numărul total de procese
}

void finalizeMPI() {
    MPI_Finalize();  
}

void readVectors(const char *xFile, const char *yFile, double *local_x, double *local_y, int N, int rank, int size) {
    int local_N = N / size; // Asumăm că N este divizibil exact cu numărul de procese pentru simplitate
    MPI_File fh_x, fh_y;
    MPI_Status status;

    // Calculul offset-ului pentru fiecare proces
    MPI_Offset offset_x = rank * local_N * sizeof(double);
    MPI_Offset offset_y = rank * local_N * sizeof(double);

    // Deschiderea fișierelor pentru citire
    MPI_File_open(MPI_COMM_SELF, xFile, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_x);
    MPI_File_open(MPI_COMM_SELF, yFile, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_y);

    // Citirea porțiunilor de vector
    MPI_File_read_at(fh_x, offset_x, local_x, local_N, MPI_DOUBLE, &status);
    MPI_File_read_at(fh_y, offset_y, local_y, local_N, MPI_DOUBLE, &status);

    // Închiderea fișierelor
    MPI_File_close(&fh_x);
    MPI_File_close(&fh_y);
}

void readMatrix(const char *matrixFile, double **local_A, int N, int rank, int size) {
    int local_N = N / size; // Numărul de rânduri din matrice alocat fiecărui proces
    MPI_File fh;
    MPI_Status status;

    // Calculăm dimensiunea fiecărui bloc local și offset-ul de la începutul fișierului
    MPI_Offset offset = rank * local_N * N * sizeof(double); // fiecare rând are N elemente

    // Alocăm spațiu pentru blocul local al matricei
    *local_A = (double *)malloc(local_N * N * sizeof(double));
    if (*local_A == NULL) {
        fprintf(stderr, "Eroare la alocarea memoriei pentru blocul local al matricei\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // Deschidem fișierul pentru citire
    MPI_File_open(MPI_COMM_SELF, matrixFile, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    // Citim blocul local al matricei
    MPI_File_read_at(fh, offset, *local_A, local_N * N, MPI_DOUBLE, &status);

    // Închidem fișierul
    MPI_File_close(&fh);
}

void distributeData(double *local_x, double *local_y, double **local_A, int N, int rank, MPI_Comm group_comm) {
    int local_N = N / MPI_Comm_size(group_comm, &size); // Presupunem divizibilitate exactă pentru simplitate
    
    if (rank == 0) {
        // Procesul de rang 0 din grup citește întreaga dată
        double *x = (double *)malloc(N * sizeof(double));
        double *y = (double *)malloc(N * sizeof(double));
        double *A = (double *)malloc(N * N * sizeof(double));
        // Presupunem că funcțiile readFullVectors și readFullMatrix sunt implementate pentru a citi datele complete
        readFullVectors("x.dat", "y.dat", x, y, N);
        readFullMatrix("mat.dat", A, N);

        // Distribuie vectorii x și y către toate procesele din grup
        MPI_Scatter(x, local_N, MPI_DOUBLE, local_x, local_N, MPI_DOUBLE, 0, group_comm);
        MPI_Scatter(y, local_N, MPI_DOUBLE, local_y, local_N, MPI_DOUBLE, 0, group_comm);

        // Pentru matrice, poate fi nevoie de o abordare custom, deoarece MPI_Scatter nu suportă direct matrici 2D
        // Aici, trimitem fiecare segment al matricei către procesul corespunzător
        for (int i = 0; i < size; i++) {
            if (i == 0) {
                // Procesul 0 deja are datele sale, așa că doar copiem direct
                memcpy(*local_A, A, local_N * N * sizeof(double));
            } else {
                // Trimite segmentul corespunzător matricei procesului i
                MPI_Send(&A[i * local_N * N], local_N * N, MPI_DOUBLE, i, 0, group_comm);
            }
        }

        // Eliberează memoria alocată pentru datele complete, dacă este necesar
        free(x);
        free(y);
        free(A);
    } else {
        // Celelalte procese primesc segmentele lor de date
        MPI_Scatter(NULL, local_N, MPI_DOUBLE, local_x, local_N, MPI_DOUBLE, 0, group_comm);
        MPI_Scatter(NULL, local_N, MPI_DOUBLE, local_y, local_N, MPI_DOUBLE, 0, group_comm);

        // Primeste segmentul matricei de la procesul 0
        MPI_Recv(*local_A, local_N * N, MPI_DOUBLE, 0, 0, group_comm, MPI_STATUS_IGNORE);
    }
}

double calculateNumerator(const double *local_x, const double *local_A, const double *local_y, int local_N, int N) {
    double sum = 0.0;
    for (int i = 0; i < local_N; i++) {
        for (int j = 0; j < N; j++) {
            sum += local_x[i] * local_A[i * N + j] * local_y[j];
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

