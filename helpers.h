#ifndef HELPERS_H
#define HELPERS_H

#include <mpi.h>

// Structura pentru a stoca dimensiunile vectorilor și matricii
typedef struct
{
    int N; // Dimensiunea N pentru vectori și matrice
} DataDims;

// Funcție pentru inițializarea mediului MPI și determinarea rangului și dimensiunii comunicatorului
void initializeMPI(int *argc, char ***argv, int *rank, int *size);

// Funcție pentru citirea vectorilor x și y din fișiere
void readVectors(const char *xFile, const char *yFile, double *local_x, double *local_y, int N, int rank, int size) {

// Funcție pentru citirea matricei A din fișier
void readMatrix(const char *matrixFile, double **A, int N);

// Funcție pentru distribuirea segmentelor de vectori și matrice către toate procesele
void distributeData(double *local_x, double *local_y, double **local_A, const double *x, const double *y, const double *A, int N, int rank, int size);

// Funcție pentru calculul parțial al numărătorului AVG
double calculateNumerator(const double *local_x, const double **local_A, const double *local_y, int local_N);

// Funcție pentru calculul parțial al numitorului AVG
double calculateDenominator(const double *local_x, const double *local_y, int local_N);

double aggregateDenominator(double partialDenominator, int rank, MPI_Comm commDenominator) {
double aggregateNumerator(double partialNumerator, int rank, MPI_Comm commNumerator) {

// Funcție pentru finalizarea MPI
void finalizeMPI();

#endif // HELPERS_H
