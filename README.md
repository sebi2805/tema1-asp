# Proiect MPI - Calculul Paralel al Mediei Ponderate

## Descriere

Proiectul implementează un algoritm de calcul paralel pentru determinarea mediei ponderate a produselor unor vectori și o matrice, folosind MPI. Problema se axează pe calculul următoarei cantități:

\[ AVG = \frac{\sum*{i,j} x_i A*{ij} y*j}{\sum*{i} x_i y_i} \]

unde `x` și `y` sunt vectori cu `N = 200` elemente, iar `A` este o matrice `N x N`. Datele sunt stocate în fișierele `x.dat`, `y.dat`, și `mat.dat` în format tabelar.

## Procesul de Calcul Paralel

### Inițializarea MPI și Crearea Grupurilor

- Inițializarea mediului MPI.
- Determinarea numărului de procese și a rangului fiecărui proces.
- Împărțirea proceselor în două grupuri pentru calculul numărătorului și numitorului, respectiv.

### Citirea și Distribuirea Datelor

- Procesele de rang 0 din fiecare grup citesc datele din fișierele de intrare.
- Distribuția datelor către celelalte procese din grupuri se face folosind operații de comunicare colective.

### Calculul Paralel

- Fiecare grup calculează în paralel partea sa din expresie, utilizând operații colective pentru sumarea parțială.

### Comunicarea și Agregarea Rezultatelor

- Rezultatele parțiale sunt combinate de procesul de rang 0 folosind `MPI_Reduce()`.
- Rezultatul final este trimis către procesul principal și stocat într-un fișier de ieșire.

### Finalizarea Programului MPI

- După finalizarea calculului și stocarea rezultatului, mediul MPI este închis.

## Sarcini și Structura de Fișiere

- `main.c`: Inițializarea MPI, împărțirea în grupuri și coordonarea proceselor.
- `makefile`: Pentru compilarea separată a modulelor și a întregului program.
- `mat.dat` : Fisier sursa care contine matricea
- `x.dat` : Fisier sursa care contine vectorul X
- `y.dat` : Fisier sursa care contine vectorul Y
- `helpers.h` : fisierul header care contine diferite functii ajutatoarea
- `helpers.c` : implementarea functiilor ajutatoarea
