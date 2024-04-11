#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int N = 200;
    int elements_per_proc = N / size;

    // Alocăm memorie pentru ambele vectori compleți și subvectorii corespunzători
    double* x = (double*)malloc(N * sizeof(double));
    double* y = (double*)malloc(N * sizeof(double));
    double* sub_x = (double*)malloc(elements_per_proc * sizeof(double));
    double* sub_y = (double*)malloc(elements_per_proc * sizeof(double));

    if (rank == 0) {
        // Procesul root citește datele din ambele fișiere
        FILE* file_x = fopen("x.dat", "r");
        FILE* file_y = fopen("y.dat", "r");
        if (!file_x || !file_y) {
            fprintf(stderr, "Eroare la deschiderea fișierelor\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (int i = 0; i < N; i++) {
            fscanf(file_x, "%lf\n", &x[i]);
            fscanf(file_y, "%lf\n", &y[i]);
        }

        fclose(file_x);
        fclose(file_y);
    }

    // Distribuie părți egale din vectorii x și y către toate procesele
    MPI_Scatter(x, elements_per_proc, MPI_DOUBLE, sub_x, elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(y, elements_per_proc, MPI_DOUBLE, sub_y, elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calculează suma elementelor din subvectorii x și y
    double local_sum_x = 0, local_sum_y = 0;
    for (int i = 0; i < elements_per_proc; i++) {
        local_sum_x += sub_x[i];
        local_sum_y += sub_y[i];
    }

     // Adună toate sumele locale pentru a obține suma totală pentru ambele vectori
    double total_sum_x = 0, total_sum_y = 0;
    MPI_Reduce(&local_sum_x, &total_sum_x, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sum_y, &total_sum_y, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Afișează sumele totale în procesul root
    if (rank == 0) {
        printf("Suma totala a elementelor vectorului x este: %lf\n", total_sum_x);
        printf("Suma totala a elementelor vectorului y este: %lf\n", total_sum_y);
    }

    // Eliberăm memoria alocată
    free(x);
    free(y);
    free(sub_x);
    free(sub_y);

    MPI_Finalize();
    return 0;
}
