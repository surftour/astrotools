EXEC   = Test

OPTIMIZE =   


OBJS   = main.o udmget.o xspexrav_callable.o

CC     = gcc -I/sw/lib/include

INCL   = cfortran.h 

LIBS   = -lm -lg2c -lstdc++ -L/sw/lib


$(EXEC): $(OBJS) 
	 $(CC) $(OBJS) $(LIBS) -o $(EXEC)  

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

