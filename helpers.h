#ifndef MPI_HELPERS_H
#define MPI_HELPERS_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 200

// Calculeaza numaratorul pentru o operatie specifica folosind o submatrice si doua subvectori. Aceasta functie este responsabila pentru multiplicarea elementelor matricei cu elementele corespunzatoare din vectorul x si y, sumand rezultatele pentru a obtine numaratorul local. Elementele procesate sunt determinate de numarul de elemente pe proces.
void compute_nominator(double *sub_mat, double *sub_x_vec, double *sub_y_vec, int elements_per_proc, double *local_nominator);

// Calculeaza numitorul folosind subvectorii asignati fiecarui proces. Functia sumeaza produsul elementelor corespunzatoare din subvectorii x si y. Numarul de randuri per proces determina cate elemente sunt procesate de aceasta functie.
void compute_denominator(double *sub_x_vec, double *sub_y_vec, int rows_per_proc, double *local_denominator);

// Functie de finalizare si curatare care elibereaza memoria alocata si inchide comunicatorii MPI deschise. Aceasta functie asigura ca toate resursele sunt corect eliberate inainte de terminarea programului, prevenind scurgerile de memorie.
void finalize_and_cleanup(int group, int group_rank, double *mat_data, double *vec_x_data, double *vec_y_data, double *sub_mat, double *sub_x_vec, double *sub_y_vec, MPI_Comm *group_comm);

// Trimite rezultatele calculului (numaratorul si numitorul global) catre procesul master (rank 0) din comunicatorul global. Functia este folosita in contextul unde fiecare grup de procese calculeaza o parte a datelor si apoi trimite rezultatele catre un proces central pentru agregare sau prelucrare ulterioara.
void send_results_to_master(int group_rank, int group, double global_nominator, double global_denominator, MPI_Comm comm);

// Receives the nominator and denominator from other processes and computes the final result. If the denominator is not zero, the function writes the computed average to a file. This function is typically executed by the master process to aggregate and finalize the computation.
void receive_and_compute(int rank, MPI_Comm comm);

// Initializeaza mediul MPI si datele necesare. Aceasta functie seteaza rangul si dimensiunea comunicatorului, grupurile, si distribuie datele matricei si vectorilor catre procesele corespunzatoare.
void initialize_mpi_and_data(int argc, char **argv, int *rank, int *size, int *group, MPI_Comm *group_comm, int *group_rank, int *group_size, double **mat_data, double **vec_x_data, double **vec_y_data, double **sub_mat, double **sub_x_vec, double **sub_y_vec, int *rows_per_proc, int *elements_per_proc);

// Initializeaza matricea si distribuie submatricile catre procese in functie de grupul si rangul acestora in grup. Aceasta functie este responsabila pentru incarcarea datelor matricei dintr-un fisier sau o sursa externa si pentru impartirea acestora intre procese.
void initialize_matrix(int group, int group_rank, MPI_Comm group_comm, double **mat_data, double **sub_mat, int rows_per_proc);

// Initializeaza vectorii x si y si distribuie subvectorii catre procese. Similar cu initializarea matricei, aceasta functie incarca datele din vectori si le distribuie proceselor pe baza rangului si grupului din MPI.
void initialize_vectors(int group, int group_rank, MPI_Comm group_comm, int elements_per_proc, double **vec_x_data, double **vec_y_data, double **sub_x_vec, double **sub_y_vec);

#endif // MPI_HELPERS_H
