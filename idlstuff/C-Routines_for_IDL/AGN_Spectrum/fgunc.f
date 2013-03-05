c Output routines for the FUNCTION object.

c Data types
c  CHARACTER*4        solar            The Solar abundance table in use.
c  CHARACTER*4        xsect            The photoelectric cross-sections in use.
c  CHARACTER*128      datdir           Directory for data files used.
c  REAL               solfil           Solar abundance table read in from a file
c  CHARACTER*128      strpar(2,MXSTPR) An array for model string parameters
c  INTEGER            nstrpr           Number of model string parameters

c Output routines
c       fgsolr          function returns the Solar abundance table in use.
c       fgxsct          function returns cross-sections in use.
c       fgabnd          function returns the Solar abundance for the input elt.
c       fgdatd          function returns the data file directory.
c       fgmstr          function returns the value for the model string parameter
c       fgnmst          subroutine returns name and value for the ith model 
c                       string parameter.

c *****************************************************************************
      FUNCTION fgsolr()

      INCLUDE 'func.inc'

      CHARACTER fgsolr*4

c Returns the Solar abundance table in use
c Arguments :
c      fgsolr   C*4    r: The Solar abundance table

      fgsolr = solar

      RETURN
      END

c *****************************************************************************
      FUNCTION fgxsct()

      INCLUDE 'func.inc'

      CHARACTER fgxsct*4

c Returns the photoelectric cross-sections in use.
c Arguments :
c      fgxsct   C*4    r: The cross-sections

      fgxsct = xsect

      RETURN
      END

c *****************************************************************************
      FUNCTION fgabnd(element)

      INCLUDE 'func.inc'

      REAL      fgabnd
      CHARACTER element*2

c Returns the Solar abundance for the input element
c Arguments :
c      fgabnd   R    r: The Solar abundance for the input element

      REAL        feld(NELTS), angr(NELTS), aneb(NELTS), grsa(NELTS)
      REAL        wilm(NELTS), lodd(NELTS)
      CHARACTER*2 elts(NELTS)

      INTEGER i

      CHARACTER fgsolr*4
      EXTERNAL  fgsolr

c Solar abundance tables. 'feld' is from Feldman, U., 1992. Physica Scripta,
c 46, 202 (abundances not included in this paper are set to grsa). 
c 'angr' is from Anders, E. & Grevesse, N., 1989. Geochimica and
c Cosmochimica Acta 53, 197. 'aneb' is from Anders, E. & Ebihara, 1982. 
c Geochimica and Cosmochimica Acta 46, 2363. 'grsa' is from Grevesse and
c Sauval, 1998, Space Science Reviews, 85, 161. 'wilm' is from Wilms, Allen,
c and McCray, 2000, ApJ 542, 914 (abundances are set to zero for those 
c elements not included in the paper). 'lodd' are the solar photospheric
c abundances from Lodders, K, 2003 ApJ 591, 1220, 'file' is from a file 
c read in by the user.

      DATA elts/'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 
     &          'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 
     &          'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 
     &          'Ni', 'Cu', 'Zn'/
      DATA feld/1.00e+00, 9.77e-02, 1.26e-11, 2.51e-11, 3.55e-10, 
     &          3.98e-04, 1.00e-04, 8.51e-04, 3.63e-08, 1.29e-04, 
     &          2.14e-06, 3.80e-05, 2.95e-06, 3.55e-05, 2.82e-07, 
     &          1.62e-05, 1.88e-07, 4.47e-06, 1.32e-07, 2.29e-06, 
     &          1.48e-09, 1.05e-07, 1.00e-08, 4.84e-07, 4.84e-07, 
     &          3.24e-05, 8.60e-08, 1.78e-06, 1.62e-08, 3.98e-08/
      DATA angr/1.00e+00, 9.77e-02, 1.45e-11, 1.41e-11, 3.98e-10, 
     &          3.63e-04, 1.12e-04, 8.51e-04, 3.63e-08, 1.23e-04, 
     &          2.14e-06, 3.80e-05, 2.95e-06, 3.55e-05, 2.82e-07, 
     &          1.62e-05, 3.16e-07, 3.63e-06, 1.32e-07, 2.29e-06, 
     &          1.26e-09, 9.77e-08, 1.00e-08, 4.68e-07, 2.45e-07, 
     &          4.68e-05, 8.32e-08, 1.78e-06, 1.62e-08, 3.98e-08/
      DATA aneb/1.00e+00, 8.01e-02, 2.19e-09, 2.87e-11, 8.82e-10, 
     &          4.45e-04, 9.12e-05, 7.39e-04, 3.10e-08, 1.38e-04, 
     &          2.10e-06, 3.95e-05, 3.12e-06, 3.68e-05, 3.82e-07, 
     &          1.89e-05, 1.93e-07, 3.82e-06, 1.39e-07, 2.25e-06, 
     &          1.24e-09, 8.82e-08, 1.08e-08, 4.93e-07, 3.50e-07, 
     &          3.31e-05, 8.27e-08, 1.81e-06, 1.89e-08, 4.63e-08/
      DATA grsa/1.00e+00, 8.51e-02, 1.26e-11, 2.51e-11, 3.55e-10, 
     &          3.31e-04, 8.32e-05, 6.76e-04, 3.63e-08, 1.20e-04, 
     &          2.14e-06, 3.80e-05, 2.95e-06, 3.55e-05, 2.82e-07, 
     &          2.14e-05, 3.16e-07, 2.51e-06, 1.32e-07, 2.29e-06, 
     &          1.48e-09, 1.05e-07, 1.00e-08, 4.68e-07, 2.45e-07, 
     &          3.16e-05, 8.32e-08, 1.78e-06, 1.62e-08, 3.98e-08/
      DATA wilm/1.00e+00, 9.77e-02, 0.00e+00, 0.00e+00, 0.00e+00,
     &          2.40e-04, 7.59e-05, 4.90e-04, 0.00e+00, 8.71e-05, 
     &          1.45e-06, 2.51e-05, 2.14e-06, 1.86e-05, 2.63e-07,
     &          1.23e-05, 1.32e-07, 2.57e-06, 0.00e+00, 1.58e-06, 
     &          0.00e+00, 6.46e-08, 0.00e+00, 3.24e-07, 2.19e-07, 
     &          2.69e-05, 8.32e-08, 1.12e-06, 0.00e+00, 0.00e+00/
      DATA lodd/1.00e+00, 7.92e-02, 1.90e-09, 2.57e-11, 6.03e-10, 
     &          2.45e-04, 6.76e-05, 4.90e-04, 2.88e-08, 7.41e-05, 
     &          1.99e-06, 3.55e-05, 2.88e-06, 3.47e-05, 2.88e-07, 
     &          1.55e-05, 1.82e-07, 3.55e-06, 1.29e-07, 2.19e-06, 
     &          1.17e-09, 8.32e-08, 1.00e-08, 4.47e-07, 3.16e-07, 
     &          2.95e-05, 8.13e-08, 1.66e-06, 1.82e-08, 4.27e-08/


      fgabnd = 0.

      DO i = 1, NELTS
         IF ( elts(i) .EQ. element ) THEN
            IF ( fgsolr() .EQ. 'feld' ) THEN
               fgabnd = feld(i)
            ELSEIF ( fgsolr() .EQ. 'angr' ) THEN
               fgabnd = angr(i)
            ELSEIF ( fgsolr() .EQ. 'aneb' ) THEN
               fgabnd = aneb(i) 
            ELSEIF ( fgsolr() .EQ. 'grsa' ) THEN
               fgabnd = grsa(i)
            ELSEIF ( fgsolr() .EQ. 'wilm' ) THEN
               fgabnd = wilm(i)
            ELSEIF ( fgsolr() .EQ. 'lodd' ) THEN
               fgabnd = lodd(i)
            ELSEIF ( fgsolr() .EQ. 'file' ) THEN
               fgabnd = solfil(i)
            ENDIF
            RETURN
         ENDIF
      ENDDO

      RETURN
      END

c *****************************************************************************
      FUNCTION fgdatd()

      INCLUDE 'func.inc'

      CHARACTER fgdatd*(*)

c Returns the data file directory in use
c Arguments :
c      fgdatd   C*(*)    r: The data file directory

      fgdatd = datdir

      RETURN
      END

c *****************************************************************************
      FUNCTION fgmstr(cvalue)

      INCLUDE 'func.inc'

      CHARACTER fgmstr*(*), cvalue*(*)

c Returns the value corresponding to the input model string parameter name.
c Arguments :
c      cvalue   C*(*)    r: The name
c      fgmstr   C*(*)    r: The value

      INTEGER i

      fgmstr = ' '
      CALL upc(cvalue)
      DO i = 1, nstrpr
         IF ( cvalue .EQ. strpar(1,i) ) THEN
            fgmstr = strpar(2,i)
            RETURN
         ENDIF
      ENDDO

      RETURN
      END

c *****************************************************************************
      SUBROUTINE fgnmst(i, cvalue1, cvalue2)

      INCLUDE 'func.inc'

      INTEGER i
      CHARACTER cvalue1*(*), cvalue2*(*)

c Returns the name and value for the ith string parameter.
c Arguments :
c      i         I        i: The parameter to return
c      cvalue1   C*(*)    r: The name
c      cvalue2   C*(*)    r: The value

      cvalue1 = ' '
      cvalue2 = ' '
      IF ( i .LT. 1 .OR. i .GT. nstrpr ) RETURN

      cvalue1 = strpar(1,i)
      cvalue2 = strpar(2,i)

      RETURN
      END
