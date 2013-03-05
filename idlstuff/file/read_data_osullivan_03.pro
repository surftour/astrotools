;========================================
;========================================
; This is the O'Sullivan, Ponman & Collins
; 2003, MNRAS, 340, 1375 data:
; X-ray scaling properties of early-type galaxies
;
pro read_data_osullivan_03, loglx, loglb, tempx


; first we'll read the '01 data, to load the L_B
; values
osullivanfile01= '/n/home/tcox/Documents/Data_E_xray/osullivan_01_all2.txt'

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
osullivanfile03= '/n/home/tcox/Documents/Data_E_xray/osullivan_03_lxfixed.txt'

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

