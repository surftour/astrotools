#
#  Makefile for code to analyze output of P-Gadget-2 simulations 
#

OBJS = 	readsnap.o allocate.o printhead.o system.o setup_units.o \
		localvalues.o cooling.o tree.o overhead.o 
CC = /usr/bin/gcc
#CC = /usr/bin/cc
F77 = f77 
INCL = overhead.h proto.h allvars.h nrutil.h 
LIB = -lm
#CFLAGS  =  -shared -m32
#CFLAGS =  -xarch=v9
#INC = 

sim :  ${OBJS} ${INC} main.c
	${CC} -c ${CFLAGS} main.c 
	${F77} -o $@ main.o ${OBJS} ${LIB} 

dustmap :  ${OBJS} ${INC} makedustmap.c 
	${CC} ${CFLAGS} makedustmap.c -o $@ ${OBJS} ${LIB}

readsnap.o : readsnap.c
	${CC} -c ${CFLAGS} readsnap.c 

allocate.o : allocate.c
	${CC} -c ${CFLAGS} allocate.c 

printhead.o : printhead.c
	${CC} -c ${CFLAGS} printhead.c 

system.o : system.c
	${CC} -c ${CFLAGS} system.c 

setup_units.o : setup_units.c
	${CC} -c ${CFLAGS} setup_units.c 

cooling.o : cooling.c
	${CC} -c ${CFLAGS} cooling.c 

tree.o : tree.c
	${CC} -c ${CFLAGS} tree.c 

localvalues.o : localvalues.c
	${CC} -c ${CFLAGS} localvalues.c 

clean: 
	/bin/rm -f core *.o sim

