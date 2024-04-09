#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

double calculate_partial_numerator(int N, double *x, double *y, double **A) {
    double partial_sum = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            partial_sum += x[i] * A[i][j] * y[j];
        }
    }
    return partial_sum;
}

// Funcția pentru calculul sumei parțiale pentru numitor
double calculate_partial_denominator(int N, double *x, double *y) {
    double partial_sum = 0.0;
    for (int i = 0; i < N; i++) {
        partial_sum += x[i] * y[i];
    }
    return partial_sum;
}

void read_data_from_file(const char *file_name, double *data, int size) {
    FILE *file = fopen(file_name, "r");
    if (!file) {
        fprintf(stderr, "Eroare la deschiderea fisierului %s\n", file_name);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < size; i++) {
        if (fscanf(file, "%lf\n", &data[i]) != 1) {
            fprintf(stderr, "Eroare la citirea datelor din fisierul %s\n", file_name);
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
}

// Functia pentru a citi matricea de date dintr-un fisier
void read_matrix_from_file(const char *file_name, double **matrix, int size) {
    FILE *file = fopen(file_name, "r");
    if (!file) {
        fprintf(stderr, "Eroare la deschiderea fisierului %s\n", file_name);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (fscanf(file, "%lf", &matrix[i][j]) != 1) {
                fprintf(stderr, "Eroare la citirea datelor din fisierul %s\n", file_name);
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }
    }

    fclose(file);
}



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    double start = MPI_Wtime();
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_size < 2) {
        fprintf(stderr, "Sunt necesare cel puțin 2 procese.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Comm group_comm;
    int group = world_rank % 2;

    MPI_Comm_split(MPI_COMM_WORLD, group, world_rank, &group_comm);

    int group_rank, group_size;
    MPI_Comm_rank(group_comm, &group_rank);
    MPI_Comm_size(group_comm, &group_size);

    double *dataX = NULL, *dataY = NULL, **dataA = NULL;
    double my_partial_sum = 0.0, total_sum = 0.0;
    int N = 200; // Presupunem că dimensiunea este 200

    if (group_rank == 0) {
        dataX = (double*)malloc(N * sizeof(double));
        dataY = (double*)malloc(N * sizeof(double));
        if (group == 0) {
            dataA = (double**)malloc(N * sizeof(double*));
            for (int i = 0; i < N; i++) {
                dataA[i] = (double*)malloc(N * sizeof(double));
            }
        }

        // Presupunem că aceste funcții sunt implementate pentru a citi datele
        read_data_from_file("x.dat", dataX, N);
        read_data_from_file("y.dat", dataY, N);
        if (group == 0) {
            read_matrix_from_file("mat.dat", dataA, N);
        }

        for (int i = 1; i < group_size; i++) {
            MPI_Send(dataX, N, MPI_DOUBLE, i, 0, group_comm);
            MPI_Send(dataY, N, MPI_DOUBLE, i, 0, group_comm);
            if (group == 0) {
                for (int j = 0; j < N; j++) {
                    MPI_Send(dataA[j], N, MPI_DOUBLE, i, 0, group_comm);
                }
            }
        }
    } else {
        dataX = (double*)malloc(N * sizeof(double));
        dataY = (double*)malloc(N * sizeof(double));
        if (group == 0) {
            dataA = (double**)malloc(N * sizeof(double*));
            for (int i = 0; i < N; i++) {
                dataA[i] = (double*)malloc(N * sizeof(double));
            }
        }

        MPI_Recv(dataX, N, MPI_DOUBLE, 0, 0, group_comm, MPI_STATUS_IGNORE);
        MPI_Recv(dataY, N, MPI_DOUBLE, 0, 0, group_comm, MPI_STATUS_IGNORE);
        if (group == 0) {
            for (int i = 0; i < N; i++) {
                MPI_Recv(dataA[i], N, MPI_DOUBLE, 0, 0, group_comm, MPI_STATUS_IGNORE);
            }
        }
    }

    my_partial_sum = (group == 0) ? calculate_partial_numerator(N, dataX, dataY, dataA) :
                                    calculate_partial_denominator(N, dataX, dataY);

    MPI_Reduce(&my_partial_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, group_comm);

  double numerator_total = 0.0, denominator_total = 0.0;

if (group_rank == 0) {
    if (world_rank == 0) {
        // Procesul cu rank 0 global primește sumele de la ambele grupuri
        if (group == 0) {
            numerator_total = total_sum; // Aceasta este suma totală a numărătorului
            MPI_Recv(&denominator_total, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            denominator_total = total_sum; // Aceasta este suma totală a numitorului
            MPI_Recv(&numerator_total, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Calculăm AVG
        if (denominator_total != 0) {
            double AVG = numerator_total / denominator_total;
            printf("Valoarea AVG este: %f\n", AVG);
        } else {
            fprintf(stderr, "Eroare: Numitorul este zero.\n");
        }
    } else {
        // Procesele cu rank 0 din fiecare grup trimit suma totală la procesul cu rank 0 global
        MPI_Send(&total_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}

    // Eliberarea memoriei
    free(dataX);
    free(dataY);
    if (group == 0) {
        for (int i = 0; i < N; i++) {
            free(dataA[i]);
        }
        free(dataA);
    }
double end = MPI_Wtime();
printf("Durata: %f secunde\n", end - start);
    MPI_Comm_free(&group_comm);
    MPI_Finalize();
    return 0;
}

