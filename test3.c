#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 200

void compute_nominator(double *sub_mat, double *sub_x_vec, double *sub_y_vec, int elements_per_proc, double *local_sum_product) {
    *local_sum_product = 0.0;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < elements_per_proc; i++) {
            *local_sum_product += sub_x_vec[i] * sub_mat[i * N + j] * sub_y_vec[j];
        }
    }
}

void compute_denominator(double *sub_x_vec, double *sub_y_vec, int rows_per_proc, double *local_sum_product_2) {
    *local_sum_product_2 = 0.0;
    for (int i = 0; i < rows_per_proc; i++) {
        *local_sum_product_2 += sub_x_vec[i] * sub_y_vec[i];
    }
}

void compute_all(double *local_nominator, double *local_denominator, double *global_nominator, double *global_denominator, MPI_Comm comm) {
    MPI_Reduce(local_nominator, global_nominator, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(local_denominator, global_denominator, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
        printf("Total sum of matrix-vector products: %f\n", *global_nominator);
        printf("Total sum of vector elements: %f\n", *global_denominator);
        if (*global_denominator != 0) { // Check for divide-by-zero
            printf("Fraction result: %f\n", *global_nominator / *global_denominator);
        } else {
            printf("Cannot compute fraction due to division by zero.\n");
        }
    }
}

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
    MPI_Scatter(vec_x_data, elements_per_proc, MPI_DOUBLE, sub_x_vec, elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(vec_y_data, N, MPI_DOUBLE, sub_y_vec, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);



    // Calculează suma pentru submatrice
    double local_nominator = 0.0, global_nominator = 0.0;
    double local_denominator = 0.0, global_denominator = 0.0;

    compute_nominator(sub_mat, sub_x_vec, sub_y_vec, elements_per_proc, &local_nominator);
    compute_denominator(sub_x_vec, sub_y_vec, rows_per_proc, &local_denominator);
    compute_all(&local_nominator, &local_denominator, &global_nominator, &global_denominator, MPI_COMM_WORLD);
   
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
