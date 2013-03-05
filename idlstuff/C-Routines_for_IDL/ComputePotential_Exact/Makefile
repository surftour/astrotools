EXEC   = CalcPotential

SRCS   = main.c 


OBJS   = $(SRCS:.c=.o)
INCL   =

.KEEP_STATE:

CFLAGS =  -g   -fPIC

LIBS   =   

CC     =  cc


$(EXEC): $(OBJS) 
	ld -shared $(OBJS) $(LIBS)  -o potential.so


$(OBJS): $(INCL) 



clean:
	rm -rf $(OBJS)   potential.so










