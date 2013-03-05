CC       =  g++            # sets the C-compiler

OPTIMIZE =   -shared -m32 


OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2)

LIB = stellar_spectrum

OBJS   = main.o stellar_spectrum.o

INCL   = stellar_spectrum.h 

CFLAGS =   $(OPTIONS) -I/usr/include

LIBS   =  -lm -L/usr/local/gsl32/lib -lgsl

$(LIB): $(OBJS) Makefile
	$(CC) $(CFLAGS) -o  $(LIB)  $(OBJS) $(LIBS)


$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

