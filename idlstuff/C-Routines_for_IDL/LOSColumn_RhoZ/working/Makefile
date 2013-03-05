
EXEC   = getnh.so

SRCS   = idlmakedustmap.c allocate.c printhead.c system.c setup_units.c \
	 localvalues.c cooling.c tree.c readsnap.c overhead.c


OBJS   = $(SRCS:.c=.o)
INCL   = overhead.h   \
	 nrutil.h


CFLAGS =  -shared -m32

LIBS   =   

CC     =  cc 

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(SRCS) 



clean:
	rm $(OBJS) $(EXEC)










