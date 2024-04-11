#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 200

int main(int argc, char *argv[]) {
    int rank, size;
    double **mat = NULL;
    double *data = NULL;
    double *recvbuf = NULL;
    double local_sum = 0.0, total_sum = 0.0;
    int rows_per_proc, i, j;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    rows_per_proc = N / size;

    // Alocă memorie pentru toată matricea în procesul 0
    if (rank == 0) {
        data = (double *)malloc(N * N * sizeof(double));
        mat = (double **)malloc(N * sizeof(double *));
        for (i = 0; i < N; i++) {
            mat[i] = &data[i * N];
        }

        // Procesul 0 citește matricea din fișier
        FILE *file = fopen("mat.dat", "r");
        if (file == NULL) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                fscanf(file, "%lf", &mat[i][j]);
            }
        }
        fclose(file);
    }

    // Alocă memorie pentru partea de matrice primită de fiecare proces
    recvbuf = (double *)malloc(rows_per_proc * N * sizeof(double));

    // Distribuie partea de matrice către toate procesele
    MPI_Scatter(data, rows_per_proc * N, MPI_DOUBLE, recvbuf, rows_per_proc * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calculează suma pentru partea primită
    for (i = 0; i < rows_per_proc * N; i++) {
        local_sum += recvbuf[i];
    }

    // Adună sumele locale în suma totală în procesul cu rank-ul 0
    MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Procesul cu rank-ul 0 afișează suma totală
    if (rank == 0) {
        printf("Suma totală a elementelor matricei este: %f\n", total_sum);
    }

    // Eliberează memoria și finalizează MPI
    if (rank == 0) {
        free(data);
        free(mat);
    }
    free(recvbuf);

    MPI_Finalize();

    return 0;
}
