# Makefile for elsrtm using MPI

CC = mpicc      #Open MPI Wrapper compiler

LFLAGS = -lm  #math library flag

SRC =    SYTOEP.c                \
         SPOL.c                  \
         DPSS.c                  \
         Util.c                  \
         main.c

elsrtm:   $(SRC)
	$(CC) $(SRC) -o dpss $(LFLAGS)

