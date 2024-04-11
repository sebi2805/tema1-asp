#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double calculNumarator(int N, double *x, double *y, double **A);
double calculNumitor(int N, double *x, double *y);

int main() {
  struct timespec start, end;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);
    int N = 200;
    // Alocare dinamică pentru doi vectori de tip double
    double *vectorX = (double*)calloc(N, sizeof(double));
    double *vectorY = (double*)calloc(N, sizeof(double));

    // Alocare dinamică pentru o matrice NxN de tip double
    double **matrice = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        matrice[i] = (double*)calloc(N, sizeof(double));
    }

    FILE *fx = fopen("x.dat", "r");
    FILE *fy = fopen("y.dat", "r");
    FILE *fmat = fopen("mat.dat", "r");

    // Verificăm dacă fișierele s-au deschis corect
    if (fx == NULL || fy == NULL || fmat == NULL) {
        printf("Eroare la deschiderea fișierelor.\n");
        // Eliberăm memoria și închidem fișierele deschise
        free(vectorX);
        free(vectorY);
        for (int i = 0; i < N; i++) {
            free(matrice[i]);
        }
        free(matrice);
        if (fx != NULL) fclose(fx);
        if (fy != NULL) fclose(fy);
        if (fmat != NULL) fclose(fmat);
        return 1;
    }

    // Citirea vectorului X
    for (int i = 0; i < N; i++) {
        fscanf(fx, "%lf", &vectorX[i]);
    }

    // Citirea vectorului Y
    for (int i = 0; i < N; i++) {
        fscanf(fy, "%lf", &vectorY[i]);
    }

    // Citirea matricei
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fscanf(fmat, "%lf", &matrice[i][j]);
        }
    }

     // Calculăm numărătorul și numitorul
    double numarator = calculNumarator(N, vectorX, vectorY, matrice);
    double numitor = calculNumitor(N, vectorX, vectorY);

    // Verificăm dacă numitorul este diferit de zero înainte de a calcula AVG
    if (numitor != 0) {
        printf("Numaratorul este: %lf\n", numarator);
        printf("Numitorul este: %lf\n", numitor);
        double AVG = numarator / numitor;
        printf("Valoarea lui AVG este: %lf\n", AVG);
    } else {
        printf("Eroare: Numitorul este zero, nu se poate calcula AVG.\n");
    }

    // Închiderea fișierelor
    fclose(fx);
    fclose(fy);
    fclose(fmat);

    // Eliberarea memoriei alocate
    free(vectorX);
    free(vectorY);
    for (int i = 0; i < N; i++) {
        free(matrice[i]);
    }
    free(matrice);
    clock_gettime(CLOCK_MONOTONIC, &end);

    elapsed = end.tv_sec - start.tv_sec;
    elapsed += (end.tv_nsec - start.tv_nsec) / 1000000000.0; // Conversie nanosecunde în secunde

    printf("Timpul de execuție: %.9f secunde\n", elapsed);

    return 0;
}

double calculNumarator(int N, double *x, double *y, double **A) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            sum += x[i] * A[i][j] * y[j];
        }
    }
    return sum;
}

double calculNumitor(int N, double *x, double *y) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += x[i] * y[i];
    }
    return sum;
}