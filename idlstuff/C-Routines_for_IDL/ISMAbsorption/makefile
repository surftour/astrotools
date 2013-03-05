EXEC   = ism_absorption

OPTIMIZE =  -O3  -fPIC # -m32


OBJS   = main.o interstellar_absorption.o \
	 manual_gsl_accel.o manual_gsl_interp.o \
	 manual_gsl_linear.o manual_gsl_bsearch.o

CC     = cc

INCL   = interstellar_absorption.h  manual_gsl_interp.h

CFLAGS = $(OPTIMIZE)

LIBS   = -lm 


$(EXEC): $(OBJS) 
	 $(CC) $(CFLAGS) -shared $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

