CC       =  cc            # sets the C-compiler

OPTIMIZE =  #-O2 #-Wall 


OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2)

EXEC   = colors


OBJS   = main.o \
	 luminosity.o

INCL   = luminosity.h

CFLAGS =  -xarch=v9  $(OPTIONS)

LIBS   =  -lm

$(EXEC): $(OBJS) 
	ld  -G $(OBJS) $(LIBS)  -o $(EXEC)
#$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  


$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

