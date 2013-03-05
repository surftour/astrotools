EXEC   = isophotfit.so

OBJS   = main.o  nrsrc/nrutil.o  nrsrc/zriddr.o


INCL   =  proto.h nrsrc/nrutil.h

CFLAGS =  -m32  -shared


LIBS   =  -lm

CC     =  cc


$(EXEC): $(OBJS) 
	 $(CC)  $(CFLAGS) $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC)


 










