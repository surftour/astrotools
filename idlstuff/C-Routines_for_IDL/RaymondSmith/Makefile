CC       =  gcc            # sets the C-compiler

OPTIMIZE =   -shared -m32 #-xarch=v9 #-O2 #-Wall 


OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2)

EXEC   = raymond_smith 

OBJS   = main.o emissivity_callable.o

INCL   = 

CFLAGS =   $(OPTIONS) -I/usr/lib/include

LIBS   =  -lm -L/usr/lib

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o  $(EXEC)


$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

