EXEC   = slicer.so

OBJS   = main.o 


INCL   = proto.h

CFLAGS =  -m32 -O3 -g -fPIC -Wall 
LNKCMD =  ld -m elf_i386  -L/usr/lib -L/usr/local/lib  -shared  




#SYSTYPE="Regatta"


ifeq ($(SYSTYPE),"Regatta")
CFLAGS =  -O2   $(OPTIONS)
LNKCMD =   cc   -bM:SRE -bnoentry -bE:export.txt
endif




LIBS   =   -lm

CC     =  cc

$(EXEC): $(OBJS)
	 $(LNKCMD)  $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(INCL)

clean:
	 rm -f $(OBJS) $(EXEC)




