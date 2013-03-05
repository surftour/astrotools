EXEC   = HsmlAndProject.so

OBJS   = main.o tree.o peano.o allvars.o

#OPTIONS += -DPERIODIC


INCL   = proto.h allvars.h

CFLAGS =  -m32 -O3 $(OPTIONS) -g -fPIC -Wall
LNKCMD =  ld -shared -m elf_i386  -L/usr/lib -L/usr/local/lib


#SYSTYPE="Regatta"


ifeq ($(SYSTYPE),"Regatta")
CFLAGS =  -O3   $(OPTIONS)
LNKCMD =   cc   -bM:SRE -bnoentry -bE:export.txt
endif




LIBS   =  -lm

CC     =  cc

$(EXEC): $(OBJS)
	 $(LNKCMD)  $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(INCL)

clean:
	 rm -f $(OBJS) $(EXEC)




