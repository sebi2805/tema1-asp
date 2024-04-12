#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "helpers.h"


int main(int argc, char * argv[]) {
    int rank, size, group;
    MPI_Comm group_comm;
    int group_rank, group_size;
    double * mat_data, * vec_x_data, * vec_y_data;
    double * sub_mat, * sub_x_vec, * sub_y_vec;
    int rows_per_proc, elements_per_proc;

    // Initialize MPI and data
    initialize_mpi_and_data(argc, argv, & rank, & size, & group, & group_comm, & group_rank, & group_size, &
        mat_data, & vec_x_data, & vec_y_data, & sub_mat, & sub_x_vec, & sub_y_vec, &
        rows_per_proc, & elements_per_proc);
    double start_time = MPI_Wtime();

    initialize_matrix(group, group_rank, group_comm, & mat_data, & sub_mat, rows_per_proc);

    initialize_vectors(group, group_rank, group_comm, elements_per_proc, & vec_x_data, & vec_y_data, & sub_x_vec, & sub_y_vec);

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
    finalize_and_cleanup(group, group_rank, mat_data, vec_x_data, vec_y_data, sub_mat, sub_x_vec, sub_y_vec, & group_comm);
    return 0;
}