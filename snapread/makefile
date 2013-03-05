EXEC   = snapread 

SRCS   = read.c loadsnap.c globvars.c

OBJS   = $(SRCS:.c=.o)
INCL   = globvars.h 

.KEEP_STATE:

OPT    =  -DHDF5  -DH5_USE_16_API   # -DTJ_VERSION

CFLAGS = $(OPT) -g  -O3  -I/usr/local/hdf5/include

LIBS   =  -lm  -lhdf5  -L/usr/local/hdf5/lib  -L/usr/local/lib

CC     =  cc      # gcc , or cc on NERSC


$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) -o $(EPATH)$(EXEC)  $(LIBS)

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)

