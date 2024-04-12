#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 200

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double *mat_data = NULL, *vec_x_data = NULL, *vec_y_data = NULL;
    double *sub_mat = NULL, *sub_x_vec = NULL, *sub_y_vec = NULL;
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
        vec_x_data = (double *)malloc(N * sizeof(double));
        vec_y_data = (double *)malloc(N * sizeof(double));
        FILE *file_x_vec = fopen("x.dat", "r");
        FILE *file_y_vec = fopen("y.dat", "r");
        if (!file_x_vec||!file_y_vec) {
            fprintf(stderr, "Eroare la deschiderea fișierului pentru vector.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (int i = 0; i < N; i++) {
            fscanf(file_x_vec, "%lf", &vec_x_data[i]);
            fscanf(file_y_vec, "%lf", &vec_y_data[i]);
        }
        fclose(file_x_vec);
        fclose(file_y_vec);
    }

    // Alocare pentru submatrice și subvector pentru fiecare proces
    sub_mat = (double *)malloc(rows_per_proc * N * sizeof(double));
    sub_x_vec = (double *)malloc(elements_per_proc * sizeof(double));
    sub_y_vec = (double *)malloc(elements_per_proc * sizeof(double));

    // Distribuirea matricei
    MPI_Scatter(mat_data, rows_per_proc * N, MPI_DOUBLE, sub_mat, rows_per_proc * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // Distribuirea vectorului
    MPI_Scatter(vec_x_data, elements_per_proc, MPI_DOUBLE, sub_x_vec, elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(vec_y_data, N, MPI_DOUBLE, sub_y_vec, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    // De aici, fiecare proces poate opera pe sub_mat și sub_vec
    double local_sum_product = 0.0, total_sum_product = 0.0;
    double local_sum_product_2 = 0.0, total_sum_product_2 = 0.0;
    // Calculează suma pentru submatrice
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < elements_per_proc; i++) {   
            local_sum_product += sub_x_vec[i] * sub_mat[i * N + j] * sub_y_vec[j];
        }
    } 

    for (int i = 0; i < rows_per_proc; i++) {
        local_sum_product_2 += sub_x_vec[i] * sub_y_vec[i];
    }

    // Reducerea sumelor vectorilor la procesul cu rank-ul 0
    MPI_Reduce(&local_sum_product_2, &total_sum_product_2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Reducerea sumelor matricelor la procesul cu rank-ul 0
    MPI_Reduce(&local_sum_product, &total_sum_product, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Procesul cu rank-ul 0 afișează sumele totale
    if (rank == 0) {
        printf("Suma totală a elementelor vectorului este: %f\n", total_sum_product);
        printf("Suma totală a elementelor vectorului este: %f\n", total_sum_product_2);
        printf("Suma totală a elementelor vectorului este: %f\n", total_sum_product/total_sum_product_2);

    }
    // Eliberare memorie și finalizare
    if (rank == 0) {
        free(mat_data);
        free(vec_x_data);
        free(vec_y_data);
    }
    free(sub_mat);
    free(sub_x_vec);
    free(sub_y_vec);

    MPI_Finalize();
    return 0;
}
