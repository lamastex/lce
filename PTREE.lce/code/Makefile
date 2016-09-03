# Raaz made this makefile from Bob's compile.csh for ease
CC=cc
#DEBUG=-DNDEBUG
OPT=#-O2
CFLAGS=-Wall -W -ansi -pedantic  $(OPT) $(DEBUG)

TARGET=ptree#purely Bob's original PTREE
SRC=allocate.c comb.c debug.c file.c main.c memory.c path.c pattern.c prob.c sample.c setup.c store.c utility.c

all: $(SRC)
	cc $(CFLAGS) -o $(TARGET) $(SRC) -lm
clean:
	rm -f $(TARGET) 
