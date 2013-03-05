;+
; NAME:
;	VAR_TYPE
;
; PURPOSE:
;	This function returns the IDL code of the variable type.
;
; CATEGORY:
;	Miscellaneous
;
; CALLING SEQUENCE:
;	Result = VAR_TYPE( invar )
;
; INPUTS:
;	Invar:  The variable to have its type returned.
;
; KEYWORD PARAMETERS:
;	HELP:  If set the function prints the name of the variable
;	       type to screen.  Default is no printing.
;	TEXT:  If set the function returns a text string instead of a number.
;
; OUTPUTS:
;	Result:  The IDL code of the variable type.  See the HELP
;	         option section of the function for interpretation
;	         of the code, or use the HELP keyword.
;
; PROCEDURE:
;	This function reads the variable type index from the SIZE
;	function.
;
; EXAMPLE:
;	Define a floating point number.
;	  x = 1.2
;	Find out its variable type.
;	  result = var_type( x, /help )
;
; MODIFICATION HISTORY:
; 	Written by:	Edward C. Wiebe, 2000-01-21.
;	Modified:	Daithi A. Stone, 2000-06-29 (changed behaviour
;			of HELP keyword).
;       Modified:       Edward C. Wiebe, 2001-05-08 (added text keyword)
;-

;***********************************************************************

FUNCTION VAR_TYPE, Invar             $
                 , HELP=helpopt      $
                 , TEXT=text
 
;***********************************************************************
; Determine Variable Type

  siz = size(invar)
  type = siz[ siz[0]+1 ]

  names = ['Undefined','Byte','Integer','Longword integer' $
          ,'Floating point','Double-precision floating'    $
          ,'Complex floating','String','Structure'         $
          ,'Double-precision complex floating'             $
          ,'Pointer','Object reference']

; If HELP is set
  if (Keyword_Set(helpopt)) then begin
    Print, names[type]
  endif
  
  if (Keyword_Set(text)) then type=names[type]

;***********************************************************************
;The End

  return, type
END
