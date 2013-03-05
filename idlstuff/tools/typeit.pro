;+
; NAME:
; TYPEIT
;
; PURPOSE:
; Finds the type class of a variable.
;
; CATEGORY:
; Programming.
;
; CALLING SEQUENCE:
; Result = TYPEIT(X)
;
; INPUTS:
;    X
; Arbitrary, doesn't even need to be defined.
;
; OUTPUTS:
; Returns the type of X as a long integer, in the (0,11) range.
;  0 Undefined
;  1 Byte
;  2 Integer
;  3 Long integer
;  4 Float
;  5 Double precision
;  6 Complex number
;  7 String
;  8 Structure
;  9   Double complex
;  10 Pointer
;  11 Object reference
;
; PROCEDURE:
; Extracts information from the SIZE function.
;
; EXAMPLE:
; To find the type class of a variable:
; IDL> print,TYPEIT(7)
;   2
; IDL> print,TYPEIT(7D)
;   5
; IDL> print,TYPEIT('7')
;   7
;
; MODIFICATION HISTORY:
; Created 15-JUL-1991 by Mati Meron, University of Chicago.
; Modified 7 November 1997 by Paul Krummel,
; CSIRO Division of Atmospheric Research.
; Added in help message and expanded header
; info (Pointers and Object references).
;
;-
Function Typeit, x, help=help
;
on_error,2
if keyword_set(help) then begin
   doc_library,'TYPEIT'
endif
;
; ++++
    dum = size(x)
    return, dum(dum(0) + 1)
;
; ++++
;
end

