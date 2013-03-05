;+
; NAME:
; FILE_SIZE
;
; PURPOSE:
; This function finds the number of bytes in an ASCII data file.
; It should be platform independent (well Windows and UNIX at least!).
;
; CATEGORY:
; Read/Write OR Input/Output.
;
; CALLING SEQUENCE:
;
; Result = FILE_SIZE(File_name)
;
; INPUTS:
; File_name: The name of the file to find the number of bytes in.
;
; OUTPUTS:
; This function returns the number of bytes in a file.
;
; PROCEDURE:
; Uses fstat to find information about the opened unit number.
;
; EXAMPLE:
; To find the size in bytes of the file test.dat enter:
; IDL> out=FILE_SIZE('test.dat')
;
; MODIFICATION HISTORY:
; Copyright    R.Bauer 2. Jan. 1996
; The idea to use fstat instead of spawn ls -l was given by Phil Williams.
;
; Modified by Paul Krummel, 12 February 1997, CSIRO Division of Atmospheric
; Research. Changed error messages to english and modified them. Added complete
; header information and usage information (help keyword).
; Added some more comments and a check to see if filename is a string.
;-

FUNCTION FILE_SIZE, filename, help=help
;
; =====>> HELP
;
on_error,2
if (N_PARAMS(0) lt 1) or keyword_set(help) then begin
   doc_library,'FILE_SIZE'
   if N_PARAMS(0) ne 1 and not keyword_set(help) then $
               message,'Incorrect number of parameters, see above for usage.'
   return,-1
;
   ioerr:
   message,'Error reading file, '+filename+', does not exist' ,/inform
   return,-1
   filerr:
   message,'Filename must be of type string' ,/inform
   return,-1
endif
;
; ++++
; Check the filename to see if it is a string, if not display error message.
if typeit(filename) ne 7L then goto, filerr
;
; ++++
; Open the file name, if it does not exist, display error message.
openr, lun, filename, /get_lun, error=err
if err ne 0 then goto, ioerr
;
; ++++
; Use the fstat function to find information about the opened unit.
stats = fstat(lun)
free_lun, lun
;
; ++++
; Return the file size!
return, stats.size
;
; ++++
end

