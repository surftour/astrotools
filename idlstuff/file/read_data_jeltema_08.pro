;========================================
;========================================
; Jeltema et al. (2008) data
;
; their own data and that of Sun et al.
;
;
pro read_data_jeltema_08, lk, lx

datafile= '/n/home/tcox/Documents/Data_E_xray/JeltemaSun_ClusterGroup_LxLkdata.txt'

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







