      INTEGER FUNCTION LENACT(CBUF)
      CHARACTER CBUF*(*)
C---
C Function to return the active length of a character string, not
C counting any trailing blanks.  N.B. an all blank string will
C return ZERO as the length.
C---
C CBUF    I    String whose length is to be measured.
C---
C 1988-Jun-13 - Standard Fortran version [AFT]
C---
      INTEGER   I
C---
      DO 190 I=LEN(CBUF),1,-1
         IF(CBUF(I:I).NE.' ') THEN
            LENACT=I
            RETURN
         END IF
  190 CONTINUE
      LENACT=0
      RETURN
      END
