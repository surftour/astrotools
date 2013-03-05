#OPT    += -DDIRECT_RAY_INTEGRATION

EXEC   = getnh.so

SRCS   = idlmakedustmap.c allocate.c printhead.c system.c setup_units.c \
	 localvalues.c cooling.c tree.c readsnap.c overhead.c \
	 ngbtree3d.c nrutil.c selectb.c no_tree.c


OBJS   = $(SRCS:.c=.o)
INCL   = overhead.h   \
	 nrutil.h




CFLAGS =  -O3  -fPIC  # -m32 $(OPT)

LIBS   =   

CC     =  cc 

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) -shared $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(SRCS) 



clean:
	rm $(OBJS) $(EXEC)










