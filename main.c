#include <mpi.h>

#include <stdio.h>

#include <stdlib.h>

#define N 200

void compute_nominator(double * sub_mat, double * sub_x_vec, double * sub_y_vec, int elements_per_proc, double * local_nominator) {
    * local_nominator = 0.0;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < elements_per_proc; i++) {
            * local_nominator += sub_x_vec[i] * sub_mat[i * N + j] * sub_y_vec[j];
        }
    }
}

void compute_denominator(double * sub_x_vec, double * sub_y_vec, int rows_per_proc, double * local_denominator) {
    * local_denominator = 0.0;
    for (int i = 0; i < rows_per_proc; i++) {
        * local_denominator += sub_x_vec[i] * sub_y_vec[i];
    }
}


void finalize_and_cleanup(int group, int group_rank, double *mat_data, double *vec_x_data, double *vec_y_data, double *sub_mat, double *sub_x_vec, double *sub_y_vec, MPI_Comm *group_comm) {
    // Clean up dynamically allocated memory
    if (group_rank == 0) {
        if (group == 0) {
            free(mat_data); // Free matrix data if it exists
        }
        free(vec_x_data); // Free x vector data
        free(vec_y_data); // Free y vector data
    }
    if (group == 0) {
        free(sub_mat); // Free sub matrix data if in group 0
    }
    free(sub_x_vec); // Free sub x vector data
    free(sub_y_vec); // Free sub y vector data

    // Free MPI communicator and finalize MPI
    MPI_Comm_free(group_comm);
    MPI_Finalize();
}

void send_results_to_master(int group_rank, int group, double global_nominator, double global_denominator, MPI_Comm comm) {
    if (group_rank == 0) {
        if (group == 0) {
            MPI_Send(&global_nominator, 1, MPI_DOUBLE, 0, 28, comm);
        } else if (group == 1) {
            MPI_Send(&global_denominator, 1, MPI_DOUBLE, 0, 5, comm);
        }
    }
}

void receive_and_compute(int rank, MPI_Comm comm) {
    if (rank == 0) {
        MPI_Status status;
        double nominator, denominator;

        MPI_Recv(&nominator, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 28, comm, &status);
        MPI_Recv(&denominator, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 5, comm, &status);

        if (denominator != 0) {
            printf("AVG: %f\n", nominator / denominator);
        } else {
            printf("Cannot compute average due to division by zero.\n");
        }
    }
}

void initialize_mpi_and_data(int argc, char **argv, int *rank, int *size, int *group, MPI_Comm *group_comm,
                             int *group_rank, int *group_size, double **mat_data, double **vec_x_data, double **vec_y_data,
                             double **sub_mat, double **sub_x_vec, double **sub_y_vec, int *rows_per_proc, int *elements_per_proc) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
    MPI_Comm_size(MPI_COMM_WORLD, size);

    int split_point = *size / 2;
    *group = *rank < split_point ? 0 : 1;

    MPI_Comm_split(MPI_COMM_WORLD, *group, *rank, group_comm);

    MPI_Comm_rank(*group_comm, group_rank);
    MPI_Comm_size(*group_comm, group_size);

    *mat_data = NULL;
    *vec_x_data = NULL;
    *vec_y_data = NULL;
    *sub_mat = NULL;
    *sub_x_vec = NULL;
    *sub_y_vec = NULL;

    *rows_per_proc = N / *group_size;
    *elements_per_proc = N / *group_size;
}

void initialize_matrix(int group, int group_rank, MPI_Comm group_comm, double **mat_data, double **sub_mat, int rows_per_proc) {
    if (group == 0) {
        // Allocate memory for the entire matrix only in group 0
        *mat_data = (double *) malloc(N * N * sizeof(double));
        if (*mat_data == NULL) {
            fprintf(stderr, "Failed to allocate memory for matrix data.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        if (group_rank == 0) {
            // Load matrix data from file only at the leader of group 0
            FILE *file_mat = fopen("mat.dat", "r");
            if (!file_mat) {
                fprintf(stderr, "Error opening file for matrix.\n");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            for (int i = 0; i < N * N; i++) {
                if (fscanf(file_mat, "%lf", &(*mat_data)[i]) != 1) {
                    fprintf(stderr, "Error reading matrix data from file.\n");
                    fclose(file_mat);
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                }
            }
            fclose(file_mat);
        }

        // Allocate memory for the sub-matrix that each process will handle
        *sub_mat = (double *) malloc(rows_per_proc * N * sizeof(double));
        if (*sub_mat == NULL) {
            fprintf(stderr, "Failed to allocate memory for sub-matrix data.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // Scatter the matrix data across processes in group 0
        MPI_Scatter(*mat_data, rows_per_proc * N, MPI_DOUBLE, *sub_mat, rows_per_proc * N, MPI_DOUBLE, 0, group_comm);
    }
}

void initialize_vectors(int group, int group_rank, MPI_Comm group_comm, int elements_per_proc,
                        double **vec_x_data, double **vec_y_data, double **sub_x_vec, double **sub_y_vec) {
    // Allocate memory for vectors
    *vec_x_data = (double *) malloc(N * sizeof(double));
    *vec_y_data = (double *) malloc(N * sizeof(double));
    if (*vec_x_data == NULL || *vec_y_data == NULL) {
        fprintf(stderr, "Failed to allocate memory for vectors.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (group_rank == 0) {
        // Load vector data from files only at the leader of group 0
        FILE *file_x_vec = fopen("x.dat", "r");
        FILE *file_y_vec = fopen("y.dat", "r");
        if (!file_x_vec || !file_y_vec) {
            fprintf(stderr, "Error opening files for vectors.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        for (int i = 0; i < N; i++) {
            if (fscanf(file_x_vec, "%lf", &(*vec_x_data)[i]) != 1 || fscanf(file_y_vec, "%lf", &(*vec_y_data)[i]) != 1) {
                fprintf(stderr, "Error reading vector data from files.\n");
                fclose(file_x_vec);
                fclose(file_y_vec);
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }
        fclose(file_x_vec);
        fclose(file_y_vec);
    }

    // Allocate memory for sub-vectors that each process will handle
    *sub_x_vec = (double *) malloc(elements_per_proc * sizeof(double));
    if (*sub_x_vec == NULL) {
        fprintf(stderr, "Failed to allocate memory for sub-vector X.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Scatter the x vector data across all processes
    MPI_Scatter(*vec_x_data, elements_per_proc, MPI_DOUBLE, *sub_x_vec, elements_per_proc, MPI_DOUBLE, 0, group_comm);

    // Handle y vector based on the group
    if (group == 0) {
        // Broadcast y vector data within group 0
        MPI_Bcast(*vec_y_data, N, MPI_DOUBLE, 0, group_comm);
    } else if (group == 1) {
        // Allocate and scatter y vector data within group 1
        *sub_y_vec = (double *) malloc(elements_per_proc * sizeof(double));
        if (*sub_y_vec == NULL) {
            fprintf(stderr, "Failed to allocate memory for sub-vector Y.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        MPI_Scatter(*vec_y_data, elements_per_proc, MPI_DOUBLE, *sub_y_vec, elements_per_proc, MPI_DOUBLE, 0, group_comm);
    }
}

int main(int argc, char * argv[]) {
    int rank, size, group;
    MPI_Comm group_comm;
    int group_rank, group_size;
    double *mat_data, *vec_x_data, *vec_y_data;
    double *sub_mat, *sub_x_vec, *sub_y_vec;
    int rows_per_proc, elements_per_proc;

    // Initialize MPI and data
    initialize_mpi_and_data(argc, argv, &rank, &size, &group, &group_comm, &group_rank, &group_size,
                            &mat_data, &vec_x_data, &vec_y_data, &sub_mat, &sub_x_vec, &sub_y_vec,
                            &rows_per_proc, &elements_per_proc);
    double start_time = MPI_Wtime();

   initialize_matrix(group, group_rank, group_comm, &mat_data, &sub_mat, rows_per_proc);

    initialize_vectors(group, group_rank, group_comm, elements_per_proc, &vec_x_data, &vec_y_data, &sub_x_vec, &sub_y_vec);

    double global_nominator = 0.0, global_denominator = 0.0;
    if (group == 0) {
        double local_nominator = 0.0;
        compute_nominator(sub_mat, sub_x_vec, vec_y_data, elements_per_proc, & local_nominator);
        MPI_Reduce( & local_nominator, & global_nominator, 1, MPI_DOUBLE, MPI_SUM, 0, group_comm);
    } else if (group == 1) {
        double local_denominator = 0.0;
        compute_denominator(sub_x_vec, sub_y_vec, rows_per_proc, & local_denominator);
        MPI_Reduce( & local_denominator, & global_denominator, 1, MPI_DOUBLE, MPI_SUM, 0, group_comm);
    }

    // Sending part
send_results_to_master(group_rank, group, global_nominator, global_denominator, MPI_COMM_WORLD);

// Receiving and computing part
receive_and_compute(rank, MPI_COMM_WORLD);

 double end_time = MPI_Wtime();

    // Calculul și afișarea duratei
    double elapsed_time = end_time - start_time;
    printf("Elapsed time: %f seconds\n", elapsed_time);
  finalize_and_cleanup(group, group_rank, mat_data, vec_x_data, vec_y_data, sub_mat, sub_x_vec, sub_y_vec, &group_comm);
    return 0;
}