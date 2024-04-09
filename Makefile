# Numele programului executabil
TARGET = mpi_program

# Compilatorul MPI
MPICC = mpicc

# Opțiuni de compilare, de exemplu -Wall pentru toate avertismentele
CFLAGS = -Wall

# Fișierele sursă
SOURCES = main.c helpers.c

# Fișierele obiect generate din fișierele sursă
OBJECTS = $(SOURCES:.c=.o)

# Regula implicită folosită atunci când rulați doar `make`
all: $(TARGET)

# Cum să construiască ținta finală
$(TARGET): $(OBJECTS)
	$(MPICC) $(CFLAGS) -o $@ $^

# Regulă generică pentru fișierele obiect
%.o: %.c
	$(MPICC) $(CFLAGS) -c $<

# Curăță proiectul
clean:
	rm -f $(TARGET) $(OBJECTS)

# Reguli suplimentare pentru a preveni conflicte cu fișierele existente
.PHONY: all clean
 