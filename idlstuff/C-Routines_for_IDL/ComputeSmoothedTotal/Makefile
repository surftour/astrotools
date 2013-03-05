
EXEC   = smtotal.so

SRCS   = main.c ngbtree2d.c  \
	 nrutil.c selectb.c sort_flt.c sort2_flt.c 


OBJS   = $(SRCS:.c=.o)
INCL   = allvars.h   \
	 nrutil.h


CFLAGS =  -shared -m32

LIBS   =   

CC     =  cc 

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(INCL) 



clean:
	rm $(OBJS) $(EXEC)










