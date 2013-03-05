;========================================
;========================================
; Mulchaey & Jeltema (2010) data
;
;
;
pro read_data_mulchaey_10, lk, lx, upperlimit

datafile= '/n/home/tcox/Documents/Data_E_xray/Mulchaey_2010_data_v2.txt'

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
        if strmid(tempjunk(8),0,1) eq "<" then lx(i)= float(strmid(tempjunk(8),1)) else lx(i)= float(tempjunk(8))     ; x 10^39 erg/s
        if tempjunk(3) eq "UL" then upperlimit(i)= 1 else upperlimit(i)= 0

	lx(i)= 39.0 + alog10(lx(i))

endfor

close, 1


end




