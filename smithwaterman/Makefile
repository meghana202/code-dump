# Makefile for genomikon demo programs

LIB = -L../genomikon -lgenomikon -lm
INC = -I../genomikon
CC = gcc
CFLAGS = -O2 -Wall -Werror 

sw: sw.o align.o dpmatrix.o 
	$(CC) $(CFLAGS) $(INC) -o sw sw.o align.o dpmatrix.o  $(LIB)

sw.o: sw.c align.c dpmatrix.c
	$(CC) $(CFLAGS) $(INC) -c sw.c
	
align.o: align.c align.h dpmatrix.h 
	$(CC) $(CFLAGS) $(INC) -c align.c
	
dpmatrix.o: dpmatrix.c dpmatrix.h
	$(CC) $(CFLAGS) $(INC) -c dpmatrix.c
	
clean:
	rm -f *.o sw

test: $(APP1)
	./sw a.fa b.fa
