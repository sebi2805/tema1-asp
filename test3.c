#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 200

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double *mat_data = NULL, *vec_data = NULL;
    double *sub_mat = NULL, *sub_vec = NULL;
    int rows_per_proc = N / size;
    int elements_per_proc = N / size; // Pentru vector

    if (rank == 0) {
        // Alocare și citire pentru matrice
        mat_data = (double *)malloc(N * N * sizeof(double));
        FILE *file_mat = fopen("mat.dat", "r");
        if (!file_mat) {
            fprintf(stderr, "Eroare la deschiderea fișierului pentru matrice.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (int i = 0; i < N * N; i++) {
            fscanf(file_mat, "%lf", &mat_data[i]);
        }
        fclose(file_mat);

        // Alocare și citire pentru vector
        vec_data = (double *)malloc(N * sizeof(double));
        FILE *file_vec = fopen("x.dat", "r");
        if (!file_vec) {
            fprintf(stderr, "Eroare la deschiderea fișierului pentru vector.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (int i = 0; i < N; i++) {
            fscanf(file_vec, "%lf", &vec_data[i]);
        }
        fclose(file_vec);
    }

    // Alocare pentru submatrice și subvector pentru fiecare proces
    sub_mat = (double *)malloc(rows_per_proc * N * sizeof(double));
    sub_vec = (double *)malloc(elements_per_proc * sizeof(double));

    // Distribuirea matricei
    MPI_Scatter(mat_data, rows_per_proc * N, MPI_DOUBLE, sub_mat, rows_per_proc * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // Distribuirea vectorului
    MPI_Scatter(vec_data, elements_per_proc, MPI_DOUBLE, sub_vec, elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // De aici, fiecare proces poate opera pe sub_mat și sub_vec
    double local_sum_product = 0.0, total_sum_product = 0.0;
    // Calculează suma pentru submatrice
    for (int i = 0; i < rows_per_proc; i++) {
        for (int i = 0; i < elements_per_proc; i++) {
            local_sum_product += sub_vec[i] * sub_mat[i * elements_per_proc + i];
    }
    }

    // Reducerea sumelor vectorilor la procesul cu rank-ul 0
    MPI_Reduce(&local_sum_product, &total_sum_product, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Reducerea sumelor matricelor la procesul cu rank-ul 0
    MPI_Reduce(&local_sum_product, &total_sum_product, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Procesul cu rank-ul 0 afișează sumele totale
    if (rank == 0) {
        printf("Suma totală a elementelor vectorului este: %f\n", total_sum_product);
    }
    // Eliberare memorie și finalizare
    if (rank == 0) {
        free(mat_data);
        free(vec_data);
    }
    free(sub_mat);
    free(sub_vec);

    MPI_Finalize();
    return 0;
}
