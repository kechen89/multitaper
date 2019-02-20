# Makefile for elsrtm using MPI

CC = mpicc      #Open MPI Wrapper compiler

LFLAGS = -lm  #math library flag

SRC =    DPSS.c                  \
         Util.c                  \
         main.c

elsrtm:   $(SRC)
	$(CC) $(SRC) -o dpss $(LFLAGS)

