#include "helpers.h"

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

void finalize_and_cleanup(int group, int group_rank, double * mat_data, double * vec_x_data, double * vec_y_data, double * sub_mat, double * sub_x_vec, double * sub_y_vec, MPI_Comm * low_comm, MPI_Comm * high_comm) {
    if (group_rank == 0) {
        if (group == 0) {
            free(mat_data); 
        }
        free(vec_x_data); 
        free(vec_y_data);  
    }
    if (group == 0) {
        free(sub_mat);  
    }
    free(sub_x_vec); 
    free(sub_y_vec);  
    if(group==0){
        MPI_Comm_free(high_comm);
    }
    else if(group==1){
        MPI_Comm_free(low_comm);
    }
    MPI_Finalize();
}

void send_results_to_master(int group_rank, int group, double global_nominator, double global_denominator, MPI_Comm comm) {
    if (group_rank == 0) {
        if (group == 0) {
            MPI_Send( & global_nominator, 1, MPI_DOUBLE, 0, 28, comm);
        } else if (group == 1) {
            MPI_Send( & global_denominator, 1, MPI_DOUBLE, 0, 5, comm);
        }
    }
}

void receive_and_compute(int rank, MPI_Comm comm) {
    if (rank == 0) {
        MPI_Status status;
        double nominator, denominator;

        MPI_Recv( & nominator, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 28, comm, & status);
        MPI_Recv( & denominator, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 5, comm, & status);

        if (denominator != 0) {
            FILE * file = fopen("avg.dat", "w");
            if (file != NULL) {
                fprintf(file, "AVG: %lf\n", nominator / denominator);
                fclose(file);
            } else {
                printf("Eroare: Nu am putut deschide fisierul.\n");
            }
        }
    }
}

void initialize_mpi_and_data(int argc, char **argv, int *rank, int *size,
                             MPI_Comm *low_comm, MPI_Comm *high_comm,
                             int *group_rank, int *group_size, int * group,
                             double **mat_data, double **vec_x_data, double **vec_y_data,
                             double **sub_mat, double **sub_x_vec, double **sub_y_vec,
                             int *rows_per_proc, int *elements_per_proc) {
    // Initializarea MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
    MPI_Comm_size(MPI_COMM_WORLD, size);

    // Calcularea split point-ului
    int split_point = *size / 2;
    * group = * rank < split_point ? 0 : 1;
    // Determinarea grupului si a rank-ului in grup
    if (*rank < split_point) {
        // High group
        MPI_Comm_split(MPI_COMM_WORLD, 1, *rank, high_comm);
        MPI_Comm_rank(*high_comm, group_rank);
        MPI_Comm_size(*high_comm, group_size);
    } else {
          // Low group
        MPI_Comm_split(MPI_COMM_WORLD, 0, *rank, low_comm);
        MPI_Comm_rank(*low_comm, group_rank);
        MPI_Comm_size(*low_comm, group_size);
    }

    *mat_data = (double *)malloc((*rows_per_proc) * N * sizeof(double));
    *vec_x_data = (double *)malloc(N * sizeof(double));
    *vec_y_data = (double *)malloc((*rows_per_proc) * sizeof(double));
    *sub_mat = (double *)malloc((*rows_per_proc) * N * sizeof(double));
    *sub_x_vec = (double *)malloc(N * sizeof(double));
    *sub_y_vec = (double *)malloc((*rows_per_proc) * sizeof(double));

    // if (!*mat_data || !*vec_x_data || !*vec_y_data || !*sub_mat || !*sub_x_vec || !*sub_y_vec) {
    //     MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    // }

    *rows_per_proc = N / *group_size;
    *elements_per_proc = N / *group_size;
}


void initialize_matrix(int group, int group_rank, MPI_Comm group_comm, double ** mat_data, double ** sub_mat, int rows_per_proc) {
    if (group == 0) {
        // Doar in grupul 0 se citeste matricea pentru doar el are nevoie de ea 
        * mat_data = (double * ) malloc(N * N * sizeof(double));
        if ( * mat_data == NULL) {
            fprintf(stderr, "Nu am putut aloca memorie matriceo.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        if (group_rank == 0) {
            FILE * file_mat = fopen("mat.dat", "r");
            if (!file_mat) {
                fprintf(stderr, "Eroare: Nu am putut deschide fisierul pentru matrice.\n");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            for (int i = 0; i < N * N; i++) {
                if (fscanf(file_mat, "%lf", & ( * mat_data)[i]) != 1) {
                    fprintf(stderr, "Eroare: Nu am citi fisierul.\n");
                    fclose(file_mat);
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                }
            }
            fclose(file_mat);
        }

        // Aloc matrice pentru submatricele fiecarui proces
        * sub_mat = (double * ) malloc(rows_per_proc * N * sizeof(double));
        if ( * sub_mat == NULL) {
            fprintf(stderr, "Eroarea nu am putut aloca memorie pentru submatrice.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // Impartirea matricei pe procese
        MPI_Scatter( * mat_data, rows_per_proc * N, MPI_DOUBLE, * sub_mat, rows_per_proc * N, MPI_DOUBLE, 0, group_comm);
    }
}

void initialize_vectors(int group, int group_rank, MPI_Comm low_comm, MPI_Comm high_comm, int elements_per_proc,
    double ** vec_x_data, double ** vec_y_data, double ** sub_x_vec, double ** sub_y_vec) {

    * vec_x_data = (double * ) malloc(N * sizeof(double));
    * vec_y_data = (double * ) malloc(N * sizeof(double));
    if ( * vec_x_data == NULL || * vec_y_data == NULL) {
        fprintf(stderr, "Failed to allocate memory for vectors.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (group_rank == 0) {
        // Citesc datele vectorilor doar la liderul grupului respectiv
        FILE * file_x_vec = fopen("x.dat", "r");
        FILE * file_y_vec = fopen("y.dat", "r");
        if (!file_x_vec || !file_y_vec) {
            fprintf(stderr, "Error opening files for vectors.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        for (int i = 0; i < N; i++) {
            if (fscanf(file_x_vec, "%lf", & ( * vec_x_data)[i]) != 1 || fscanf(file_y_vec, "%lf", & ( * vec_y_data)[i]) != 1) {
                fprintf(stderr, "Error reading vector data from files.\n");
                fclose(file_x_vec);
                fclose(file_y_vec);
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }
        fclose(file_x_vec);
        fclose(file_y_vec);
    }

    * sub_x_vec = (double * ) malloc(elements_per_proc * sizeof(double));
    if ( * sub_x_vec == NULL) {
        fprintf(stderr, "Failed to allocate memory for sub-vector X.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if(group==0){
        MPI_Scatter( * vec_x_data, elements_per_proc, MPI_DOUBLE, * sub_x_vec, elements_per_proc, MPI_DOUBLE, 0,  high_comm);
    } else if(group==1)
        {MPI_Scatter( * vec_x_data, elements_per_proc, MPI_DOUBLE, * sub_x_vec, elements_per_proc, MPI_DOUBLE, 0, low_comm);
    }

    if (group == 0) {
        // Dau broadcast la vectorul y doar in grupul 0 pentru ca am nevoie de toate coloanele matricei 
        MPI_Bcast( * vec_y_data, N, MPI_DOUBLE, 0, high_comm);
    } else if (group == 1) {
        
        // In grupul 1, fiecare proces primeste un sub-vector y pentru ca pot sa fac inmultirea cu vectorul x
        * sub_y_vec = (double * ) malloc(elements_per_proc * sizeof(double));
        if ( * sub_y_vec == NULL) {
            fprintf(stderr, "Failed to allocate memory for sub-vector Y.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        MPI_Scatter( * vec_y_data, elements_per_proc, MPI_DOUBLE, * sub_y_vec, elements_per_proc, MPI_DOUBLE, 0, low_comm);
    }
}
