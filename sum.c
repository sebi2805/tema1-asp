#include <stdio.h>
#include <stdlib.h>

int main() {
    FILE *file;
    float value, sum = 0.0;

    // Încercăm să deschidem fișierul
    file = fopen("x.dat", "r");
    if (file == NULL) {
        perror("Eroare la deschiderea fișierului");
        return EXIT_FAILURE;
    }

    // Citim valorile float și calculăm suma
    while (fscanf(file, "%f", &value) == 1) {
        sum += value;
    }

    // Închidem fișierul
    fclose(file);

    // Afișăm suma
    printf("Suma valorilor este: %f\n", sum);

    return 0;
}
