CC=gcc
MCC=mpicc
CFLAGS=-Wall -g -std=c99
SRC=./src
INCLUDE=./include
BIN=./bin

#all: 
#	$(CC) $(CFLAGS) -I$(INCLUDE) $(SRC)/single.c -o $(BIN)/single
all: 
	$(MCC) $(CFLAGS) -I$(INCLUDE) $(SRC)/single_MPI.c -o $(BIN)/single_MPI

clean:
	rm  $(BIN)/*
