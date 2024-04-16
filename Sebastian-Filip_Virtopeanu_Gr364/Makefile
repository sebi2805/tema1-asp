CC=mpicc

CFLAGS=-Wall -O2

INCLUDES=

LFLAGS=

LIBS=

SRCS=main.c helpers.c

OBJS=$(SRCS:.c=.o)

MAIN=program

.PHONY: depend clean

all:    $(MAIN)
	@echo  "Program has been compiled."

$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE  
