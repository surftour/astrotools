EXEC   = MakePics

SRCS   = main.c ngbtree2d.c  \
	 nrsrc/nrutil.c nrsrc/selectb.c nrsrc/sort_flt.c nrsrc/sort2_flt.c 


OBJS   = $(SRCS:.c=.o)
INCL   = allvars.h   \
	 nrsrc/nrutil.h

.KEEP_STATE:

CFLAGS =  -O3 -fPIC

LIBS   =   

CC     =  cc


$(EXEC): $(OBJS) 
	ld -shared $(OBJS) $(LIBS)  -o nearestval.so


$(OBJS): $(INCL) 




clearn:
	rm -rf $(OBJS) nearestval.so









