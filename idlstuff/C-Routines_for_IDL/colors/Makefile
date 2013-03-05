
EXEC   = colors

OBJS   = main.o colors.o

INCL   = colors.h 


CFLAGS =  -O3 -g  -fPIC  #  -Wall   -I/usr/include

LIBS   =  -lm

CC     =   c++   # g++

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -shared $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

