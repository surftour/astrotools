
c Include file for FUNCTION object

c  CHARACTER*4       solar            The Solar abundance table in use.
c  CHARACTER*4       xsect            The photoelectric cross-sections in use.
c  CHARACTER*128     datdir           Directory for data files used.
c  REAL              solfil(NELTS)    Solar abundance table read in from a file
c  CHARACTER*128     strpar(2,MXSTPR) An array for model string parameters
c  INTEGER           nstrpr           Number of model string parameters

      INTEGER MXSTPR
      PARAMETER (MXSTPR=50)
      INTEGER NELTS
      PARAMETER (NELTS=30)

      INTEGER   nstrpr
      REAL      solfil(NELTS)
      CHARACTER datdir*128, solar*4, xsect*4
      CHARACTER strpar(2,MXSTPR)*128

      COMMON /funcco1/ solfil, nstrpr
      COMMON /funcco2/ solar, xsect, datdir, strpar

