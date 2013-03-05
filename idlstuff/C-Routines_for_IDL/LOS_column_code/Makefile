#OPT    += -DDIRECT_RAY_INTEGRATION

EXEC   = getnh.so

SRCS   = idlmakedustmap.c allocate.c printhead.c system.c setup_units.c \
	 localvalues.c cooling.c tree.c readsnap.c overhead.c \
	 ngbtree3d.c nrutil.c selectb.c no_tree.c


OBJS   = $(SRCS:.c=.o)
INCL   = overhead.h   \
	 nrutil.h


CFLAGS =  -shared -m32 $(OPT)
#CFLAGS =  -bundle $(OPT)

LIBS   =   

CC     =  gcc 

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(SRCS) 



clean:
	rm $(OBJS) $(EXEC)










