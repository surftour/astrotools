CC       =  cc            # sets the C-compiler

OPTIMIZE =  -shared -m32  #-O2 #-Wall 


OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2)

EXEC   = colors


OBJS   = main.o \
	 pl.o  interp.o

INCL   = pl.h

CFLAGS =   $(OPTIONS)

LIBS   =  -lm

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)   -o  $(EXEC)  


$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

