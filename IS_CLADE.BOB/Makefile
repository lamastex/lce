# Raaz made this makefile from Bob's compile.csh for ease
# This is purely Bob's original code
CC=gcc
#DEBUG=-DNDEBUG
OPT=#-O
CFLAGS=-Wall -W -ansi -pedantic  $(OPT) $(DEBUG)
TARGETS=constant exponential exp_flat
SRC=is.c time.c uniform.c

all: $(SRC)
	$(CC) -O -o constant -DCONSTANT is.c uniform.c -lm
	$(CC) -O -o exponential -DEXPONENTIAL is.c uniform.c -lm
	$(CC) -O -o exp_flat -DEXP_FLAT is.c uniform.c -lm
clean:
	rm -f $(TARGETS) 
