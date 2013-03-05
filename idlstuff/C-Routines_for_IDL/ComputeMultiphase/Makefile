CC       =  gcc            # sets the C-compiler

OPTIMIZE =  -shared -m32  #-O2 #-Wall 
#OPT1     =  -xarch=v9
OPT2     =  -DCOOLING

OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2)

EXEC   = mphase


OBJS   = main.o \
	 cooling.o 

INCL   = allvars.h cooling.h proto.h prototypes.h

CFLAGS =   $(OPTIONS)

LIBS   =   -lm 

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $(EXEC)



$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

