#OPT    += -DDIRECT_RAY_INTEGRATION

EXEC   = attenuate.so

SRCS   = attenuate.c


OBJS   = $(SRCS:.c=.o)
INCL   = 


CFLAGS =  -shared -m32 $(OPT)
#CFLAGS =  -bundle $(OPT)

LIBS   =   

CC     =  gcc 

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(SRCS) 



clean:
	rm $(OBJS) $(EXEC)










