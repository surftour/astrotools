OPT    += -DUSE_FULL_NEIGHBOR_CALC  # use the full neighbor calculation rather than tree
				     # this is more reliable but slower

EXEC   = getnh.so

SRCS   = idlmakedustmap.c allocate.c printhead.c system.c setup_units.c \
	 localvalues.c cooling.c tree.c readsnap.c overhead.c no_tree.c \
	 ngbtree3d.c nrutil.c selectb.c


OBJS   = $(SRCS:.c=.o)
INCL   = overhead.h   \
	 nrutil.h


CFLAGS =  -shared -m32 $(OPT)

LIBS   =   

CC     =  cc 

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(SRCS) 



clean:
	rm $(OBJS) $(EXEC)










