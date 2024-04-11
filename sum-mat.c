#include <stdio.h>
#include <stdlib.h>

#define N 200

int main() {
    FILE *file;
    double matrix[N][N];
    int i, j;
    double sum = 0;

    // Deschide fișierul
    file = fopen("mat.dat", "r");
    if (file == NULL) {
        perror("Eroare la deschiderea fișierului");
        return 1;
    }

    // Citește valorile din fișier și le stochează în matrice
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            // Presupunem că valorile sunt separate prin spațiu sau sfârșit de linie
            if (fscanf(file, "%lf", &matrix[i][j]) != 1) {
                perror("Eroare la citirea din fișier");
                fclose(file);
                return 1;
            }
            sum += matrix[i][j];
        }
    }

    // Închide fișierul
    fclose(file);

    // Afișează suma
    printf("Suma elementelor matricei este: %lf\n", sum);

    return 0;
}
