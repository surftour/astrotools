;========================================
;
;   Read the O'Sullivan data
;
;========================================
; This is the O'Sullivan, Forbes & Ponman
; 2001, MNRAS, 328, 461 data:
; A catalogue and analysis of L_X of early-type galaxies 
;
pro read_data_osullivan_01, loglb, loglx, ttype


; this first file has spaces in the name which 
; screws things up, use  the second one.
;osullivanfile= '/home/tcox/osullivan_01_all.txt'
osullivanfile= '/n/home/tcox/Documents/Data_E_xray/osullivan_01_all2.txt'


;
;  Format of this file is:
;#
;#
;#
;#  Name		D	Log LB	 	Log LX	Source	T
;# 		(Mpc)	(LB)	 	(erg s1)	 	 
;ESO10114	30.12	9.93*	<	41.02	B	3.0
;ESO1074	38.89	10.22	<	40.94	B	4.0
; etc...
;

spawn, "wc "+osullivanfile,result
lines=long(result)
datalines=lines(0)-5
loglb= fltarr(datalines)
loglx= fltarr(datalines)
ttype= fltarr(datalines)

openr, 1, osullivanfile
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
	readf, 1, junk
	;print, junk
	tempjunk= strsplit(junk,/extract,count=count)
	name= tempjunk(0)
	;distance= float(tempjunk(1))
	loglb(i)= float(tempjunk(2))   ; doesn't need a trap for asterisk
	if count eq 7 then begin
		idontknow= tempjunk(3)
		loglx(i)= float(tempjunk(4))
		source= tempjunk(5)
		ttype(i)= float(tempjunk(6))
	endif else begin
		loglx(i)= float(tempjunk(3))
		source= tempjunk(4)
		ttype(i)= float(tempjunk(5))
	endelse

endfor

close, 1



end




