PRO readmatrix,a,name,columns=Columns,lines=Lines,omitLines=OmitLines,ReadLines=ReadLines,SILENT=silent, $
              maxcolumns=maxcolumns,quiet=quiet,lastLines=lastLines
IF NOT(IS_DEF(name)) THEN BEGIN
   print,"SYNTAX: readmatrix,mymatrix,filemane,columns=Columns,lines=Lines,"
   print,"                   omitLines=OmitLines,ReadLines=ReadLines,SILENT=silent,"
   print,"                   maxcolumns=maxcolumns,quiet=quiet"
   return
END
tmp_preparefile,name
IF N_ELEMENTS(maxcolumns) EQ 0 THEN BEGIN
   maxcolumns=1000
ENDIF
IF N_ELEMENTS(name) EQ 0 THEN BEGIN
  name=PICKFILE(/NOCONFIRM,/READ,/MUST_EXIST,FILTER='hrd.plot')
  IF name EQ '' THEN RETURN
ENDIF

IF N_ELEMENTS(silent) EQ 0 THEN silent=0
IF IS_DEF(quiet) then silent=1

IF N_ELEMENTS(lastLines) EQ 0 THEN lastLines=0

IF silent EQ 0 THEN BEGIN
    PRINT,'Bestimme Datenformat...'
ENDIF
IF N_ELEMENTS(omitLines) EQ 0 THEN omitLines=0
x=''
OPENR,unit,name,/GET_LUN
FOR i=1,omitLines+1 DO READF,unit,x
x=STRCOMPRESS(STRTRIM(x,2))
columns=0L
i=0L
WHILE (i NE -1) DO BEGIN
  i=strpos(x,' ',i)
  IF (i ne -1) THEN i=i+1
  columns=columns+1
  IF(columns GE maxcolumns) THEN i=-1
ENDWHILE
IF N_ELEMENTS(ReadLines) EQ 0 THEN BEGIN
  lines=1L
  WHILE (NOT EOF(unit)) DO BEGIN
    READF,unit,x
    lines=lines+1L
  ENDWHILE
ENDIF ELSE BEGIN
  lines = ReadLines
ENDELSE 
lines = lines + lastLines
CLOSE,unit
IF silent EQ 0 THEN BEGIN
    PRINT,'column:',columns,' lines:',lines
ENDIF
OPENR,unit,name
a=DBLARR(columns,lines)
FOR i=1,omitLines DO READF,unit,x
IF maxcolumns EQ 1000 THEN BEGIN
    READF,unit,a
ENDIF ELSE BEGIN
   bbb=fltarr(columns)
   FOR i=1L,lines DO BEGIN
      READF,unit,bbb
      a(*,i-1)=bbb
   ENDFOR
ENDELSE
CLOSE,unit
FREE_LUN,unit
a=TRANSPOSE(a)
IF silent EQ 0 THEN BEGIN
    PRINT,'Daten engelesen.'
ENDIF
RETURN
END

