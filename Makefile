CC=gcc
FLAGS = -g -Wall
TARGET = code
CFLAGS = -fopenmp 
LDFLAGS = -lm -O4
SPARSE = csparse

all: code.o csparse.o blas.h
	$(CC) $(CFLAGS) $(FLAGS) $(TARGET).o $(SPARSE).o -o $(TARGET) $(LDFLAGS)
csparse.o: csparse.c csparse.h
	$(CC) $(CFLAGS) $(FLAGS) -c $(SPARSE).c -o $(SPARSE).o $(LDFLAGS)
code.o: code.c 
	$(CC) $(CFLAGS) $(FLAGS) -c $(TARGET).c -o $(TARGET).o $(LDFLAGS)
clean:
	rm $(TARGET) $(SPARSE).o $(TARGET).o -f
