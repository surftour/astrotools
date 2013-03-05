;========================================
;
;   Read the O'Sullivan data
;
;========================================
; This is the O'Sullivan, Forbes & Ponman
; 2001, MNRAS, 328, 461 data:
; A catalogue and analysis of L_X of early-type galaxies 
;
pro read_osullivan_01, loglb, loglx, ttype


; this first file has spaces in the name which 
; screws things up, use  the second one.
;osullivanfile= '/home/tcox/osullivan_01_all.txt'

osullivanfile= '/home/tcox/OSullivan/osullivan_01_all2.txt'


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




;========================================
;========================================
; This is the O'Sullivan, Ponman & Collins
; 2003, MNRAS, 340, 1375 data:
; X-ray scaling properties of early-type galaxies
;
pro read_osullivan_03, loglx, loglb, tempx


; first we'll read the '01 data, to load the L_B
; values
osullivanfile01= '/home/tcox/OSullivan/osullivan_01_all2.txt'

spawn, "wc "+osullivanfile01,result
lines=long(result)
datalines01=lines(0)-5
name01= strarr(datalines01)
loglb01= fltarr(datalines01)
loglx01= fltarr(datalines01)
ttype01= fltarr(datalines01)

openr, 1, osullivanfile01
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
	tempjunk= strsplit(junk,/extract,count=count)
	name01(i)= tempjunk(0)
	loglb01(i)= float(tempjunk(2))   ; doesn't need a trap for asterisk
	if count eq 7 then begin
		loglx01(i)= float(tempjunk(4))
		ttype01(i)= float(tempjunk(6))
	endif else begin
		loglx01(i)= float(tempjunk(3))
		ttype01(i)= float(tempjunk(5))
	endelse

endfor

close, 1





; next we'll actually read the 03 data, and fix
; the values to send back
osullivanfile03= '/home/tcox/OSullivan/osullivan_03_lxfixed.txt'

spawn, "wc "+osullivanfile03,result
lines=long(result)
datalines03=lines(0)-5
loglb= fltarr(datalines03)
loglx= fltarr(datalines03)
tempx= fltarr(datalines03)

openr, 1, osullivanfile03
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
        tempjunk= strsplit(junk,/extract,count=count)
        name= tempjunk(0)
	idx= where(name eq name01)
	if idx(0) ne -1 then begin
		loglb(i)= loglb01(idx)
		;nh(i)= float(tempjunk(1))                 ; 10^21 cm^2
		tempx(i)= float(tempjunk(2))               ; keV
		;metallicity(i)= float(tempjunk(3))        ; in solar
		loglx(i)= float(tempjunk(4))
	endif else begin
		print, "PROBLEM: can't find ",name
	endelse

endfor

close, 1




end




;========================================
;========================================
; This is the O'Sullivan, Ponman & Collins
; 2003, MNRAS, 340, 1375 data (informational):
; X-ray scaling properties of early-type galaxies
;
pro read_osullivan_03b, sigma


; we read the 03 data, note that this is mainly
; informational data, so sigma, R_e, etc.  L_x, T_x,
; are all enclosed above.
osullivanfile03= '/home/tcox/OSullivan/osullivan_03_infofixed.txt'

spawn, "wc "+osullivanfile03,result
lines=long(result)
datalines03=lines(0)-5
sigma= fltarr(datalines03)

openr, 1, osullivanfile03
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
        tempjunk= strsplit(junk,/extract,count=count)
        name= tempjunk(0)
	sigma(i)= float(tempjunk(1))                 ; km sec^-1
	;distance(i)= float(tempjunk(3))             ; Mpc
	;re(i)= float(tempjunk(4))                   ; arcmin
	;ttype(i)= float(tempjunk(5))
	;environment(i)= float(tempjunk(6))

endfor

close, 1




end







;========================================
;========================================
; Mulchaey & Jeltema (2010) data
;
;
;
pro read_mulchaey_10, lk, lx, upperlimit

datafile= '/home/tcox/Documents/Data_E_xray/Mulchaey_2010_data_v2.txt'

spawn, "wc "+datafile,result
lines=long(result)
;datalines=lines(0)-5
datalines=lines(0)-2
lk= fltarr(datalines)
lx= fltarr(datalines)
upperlimit= fltarr(datalines)

openr, 1, datafile
junk=''

; read the header
readf, 1, junk
readf, 1, junk
;readf, 1, junk
;readf, 1, junk
;readf, 1, junk

; read the data
for i=0,datalines-1 do begin
        readf, 1, junk
        tempjunk= strsplit(junk,/extract,count=count)
        name= tempjunk(0)
        lk(i)= float(tempjunk(7))
        lx(i)= float(tempjunk(8))     ; x 10^39 erg/s
        if tempjunk(3) eq "UL" then upperlimit(i)= 1 else upperlimit= 0

endfor

close, 1


end









;========================================
;========================================
; Jeltema et al. (2008) data
;
; their own data and that of Sun et al.
;
;
pro read_jeltema_08, lk, lx

datafile= '/home/tcox/Documents/Data_E_xray/JeltemaSun_ClusterGroup_LxLkdata.txt'

spawn, "wc "+datafile,result
lines=long(result)
;datalines=lines(0)-5
datalines=lines(0)-1
lk= fltarr(datalines)
lx= fltarr(datalines)

openr, 1, datafile
junk=''

; read the header
readf, 1, junk
;readf, 1, junk
;readf, 1, junk
;readf, 1, junk
;readf, 1, junk

; read the data
for i=0,datalines-1 do begin
        readf, 1, junk
        tempjunk= strsplit(junk,/extract,count=count)
        name= tempjunk(0)
        lk(i)= float(tempjunk(1))
        lx(i)= float(tempjunk(2))     ; x 10^39 erg/s

endfor

close, 1


end







