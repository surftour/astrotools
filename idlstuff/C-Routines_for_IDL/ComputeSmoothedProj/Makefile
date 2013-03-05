CC       =  cc            # sets the C-compiler

OPTIMIZE =   -shared -m32 #-xarch=v9 #-O2 #-Wall 


OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2)

LIB = smooth.so 


OBJS   = main.o \
         begrun.o \
	 endrun.o global.o \
	 b_init.o \
	 restart.o  io.o \
	 allvars.o \
       	 cooling.o  ewald.o \
	 system.o allocate.o \
	 b_read_ic.o \
         forcetree.o bdensity.o


INCL   = allvars.h proto.h forcetree.h cooling.h

CFLAGS =   $(OPTIONS)

LIBS   =  -lm

$(LIB): $(OBJS) Makefile
	$(CC) $(CFLAGS) -o  $(LIB)  $(OBJS)


$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

