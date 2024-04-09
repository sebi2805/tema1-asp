#ifndef HELPERS_H
#define HELPERS_H

#include <mpi.h>

typedef struct
{
    int N;  
} DataDims;

void initializeMPI(int *argc, char ***argv, int *rank, int *size);
void readVectors(const char *xFile, const char *yFile, double *local_x, double *local_y, int N, int rank, int size);
void readMatrix(const char *matrixFile, double **local_A, int N, int rank, int size);
void distributeData(double *local_x, double *local_y, double **local_A, const double *x, const double *y, const double *A, int N, int local_N, int rank, int size);
double calculateNumerator(const double *local_x, const double **local_A, const double *local_y, int local_N, int N);
double calculateDenominator(const double *local_x, const double *local_y, int local_N);
double aggregateDenominator(double partialDenominator, int rank, MPI_Comm commDenominator);
double aggregateNumerator(double partialNumerator, int rank, MPI_Comm commNumerator);
void finalizeMPI();

#endif // HELPERS_H
