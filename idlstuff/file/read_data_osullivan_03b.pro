;========================================
;========================================
; This is the O'Sullivan, Ponman & Collins
; 2003, MNRAS, 340, 1375 data (informational):
; X-ray scaling properties of early-type galaxies
;
pro read_data_osullivan_03b, sigma


; we read the 03 data, note that this is mainly
; informational data, so sigma, R_e, etc.  L_x, T_x,
; are all enclosed above.
osullivanfile03= '/n/home/tcox/Documents/Data_E_xray/osullivan_03_infofixed.txt'

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








