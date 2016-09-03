# Raaz made this makefile from Bob's compile.csh for ease
CC=cc
#DEBUG=-DNDEBUG
OPT=#-O2
CFLAGS=-Wall -W -ansi -pedantic  $(OPT) $(DEBUG)

TARGET=simulate#purely Bob's original PTREE
SRC=SIMULATE.C

all: $(SRC)
	cc $(CFLAGS) -o $(TARGET) $(SRC) 
clean:
	rm -f $(TARGET) 
