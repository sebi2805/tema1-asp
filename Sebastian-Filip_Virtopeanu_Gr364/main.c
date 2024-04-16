#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "helpers.h"

// helpers.h contine prototipurile functiilor folosite in main.c
// helpers.c contine implementarea functiilor folosite in main.c

// abordarea de parelelizare a codului este bazata pe impartirea datelor si a sarcinilor intre doua grupuri de procese
// astfel pentru numarator impart matricea in N/group_size linii si vectorul x in N/group_size elemente
// am nevoie de intreg y pentru ca fac suma pe toate coloanele matricei, as putea sa fac o preprocesare a matricei inmultite
// cu vectorul y, tot paralel, dar asta ar insemna sa fac mai multe operatii de comunicare si ar creste complexitatea programului

// group_rank = 0 citestea matricea si o imprastie in N/group_size linii, pastrand si el o submatrice
// fiecare proces din grup face inmultirea x[i] * A[i][j] * y[j] si trimite suma la liderul grupului

// pentru numitor impart vectorul x in N/group_size elemente si vectorul y in N/group_size elemente
// am access la a face inmultirea celor 2 paralel si sa fac suma pe toate elementele
// fiecare proces din grup face inmultirea x[i] * y[j] si trimite suma la liderul grupului


// procesul de rank 0 din comunicatorul global asteapta raspunsurile celor 2 grupuri si scrie in fisierul avg.dat
// la final in terminal sunt afisati timpii de executie
int main(int argc, char * argv[]) {
    int rank, size, group;
    MPI_Comm low_comm, high_comm;
    int group_rank, group_size;
    double *mat_data, *vec_x_data, *vec_y_data;
    double *sub_mat, *sub_x_vec, *sub_y_vec;
    int rows_per_proc, elements_per_proc;

    initialize_mpi_and_data(argc, argv, &rank, &size,
                            &low_comm, &high_comm, &group_rank, &group_size, &group,
                            &mat_data, &vec_x_data, &vec_y_data, 
                            &sub_mat, &sub_x_vec, &sub_y_vec,
                            &rows_per_proc, &elements_per_proc);
    double start_time = MPI_Wtime();

    initialize_matrix(group, group_rank, high_comm, & mat_data, & sub_mat, rows_per_proc);

    initialize_vectors(group, group_rank, low_comm, high_comm, elements_per_proc, & vec_x_data, & vec_y_data, & sub_x_vec, & sub_y_vec);

    double global_nominator = 0.0, global_denominator = 0.0;
    if (group == 0) {
        double local_nominator = 0.0;
        compute_nominator(sub_mat, sub_x_vec, vec_y_data, elements_per_proc, & local_nominator);
        MPI_Reduce( & local_nominator, & global_nominator, 1, MPI_DOUBLE, MPI_SUM, 0, high_comm);
    } else if (group == 1) {
        double local_denominator = 0.0;
        compute_denominator(sub_x_vec, sub_y_vec, rows_per_proc, & local_denominator);
        MPI_Reduce( & local_denominator, & global_denominator, 1, MPI_DOUBLE, MPI_SUM, 0, low_comm);
    }

    send_results_to_master(group_rank, group, global_nominator, global_denominator, MPI_COMM_WORLD);

    receive_and_compute(rank, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();

    // Calculul și afișarea duratei
    double elapsed_time = end_time - start_time;
    printf("Elapsed time: %f seconds\n", elapsed_time);
    finalize_and_cleanup(group, group_rank, mat_data, vec_x_data, vec_y_data, sub_mat, sub_x_vec, sub_y_vec, & low_comm, & high_comm);
    return 0;
}