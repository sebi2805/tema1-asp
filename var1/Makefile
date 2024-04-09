# Numele executabilului
TARGET = prog

# Compilatorul: mpicc pentru MPI
CC = mpicc

# Opțiuni de compilare, de exemplu: -Wall pentru toate warning-urile
CFLAGS = -Wall

# .c files
SRC = main.c helpers.c

# .h files
HDR = helpers.h

# Obiect files
OBJ = $(SRC:.c=.o)

# Regula implicită
all: $(TARGET)

# Linking the program
$(TARGET): $(OBJ)
	$(CC) -o $@ $^

# Compiling every .c to .o
%.o: %.c $(HDR)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean-up
clean:
	rm -f $(OBJ) $(TARGET)
