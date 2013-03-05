;========================================
;
;   Read the Bender data
;
;========================================




pro read_bender_data, filename, ellip, mb, vmax, sigma, $
				no_neg_trap=no_neg_trap


;
;  Format of this file is:
; -------------------------
;#
;#
;#
;# 		Ellip	M_b		V_max	Sigma
;#					(km/s)	(km/s)
;NGC205		0.40	-15.7		20	50
;M32		0.20	-15.5		15	60
;Fornax		0.33	-12.6		1.7	6
; etc...
;

spawn, "wc "+filename,result
lines=long(result)
datalines=lines(0)-5
ellip= fltarr(datalines)
mb= fltarr(datalines)
vmax= fltarr(datalines)
sigma= fltarr(datalines)

openr, 1, filename
print, "opening: ", filename
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
        ellip(i)= float(tempjunk(1))
        mb(i)= float(tempjunk(2))
        vmax(i)= float(tempjunk(3))
        sigma(i)= float(tempjunk(4))
endfor

close, 1


; trap for -1 values
if not keyword_set(no_neg_trap) then begin
   idx= where((vmax lt 0.0) or (sigma lt 0.0))
   if idx(0) ne -1 then begin
	vmax(idx)= -1.0
	sigma(idx)= 1.0
	print,"Number of data points: ", n_elements(ellip)
	print,"Fixing ", n_elements(idx), " for negative values."
   endif
endif


end







;================================================================================






; This is the Bender (1988) data:
; 
; 
;
pro oplot_bender_vsig_data, junk, addkey=addkey



; --------

benderfile= '/n/home/tcox/Documents/Data_RotvAnisSupport/bender88_ellipticals.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=2, symsize=1.0, color= 0
; --------


benderfile= '/n/home/tcox/Documents/Data_RotvAnisSupport/bender90_dwarfe_ellipticals.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=1, symsize=1.0, color= 0
; --------


benderfile= '/n/home/tcox/Documents/Data_RotvAnisSupport/bender90_lowlum_ellipticals.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=1, symsize=1.0, color= 0
; --------


;---------------------------------------------------------
; normal coordinates
x0=0.25
y0=0.93

; data coordinates
x0d= 0.05
y0d= 1.58

if keyword_set(addkey) then begin
	xyouts, x0+.01, y0-0.07, 'Bender (1988)', color= 0, size=0.8, charthick=2.0, /normal
	oplot, [x0d], [y0d-0.06], psym=2, color= 0

	xyouts, x0+.01, y0-0.11, 'Bender & Nieto (1990)', color= 0, size=0.8, charthick=2.0, /normal
	oplot, [x0d], [y0d-0.14], psym=1, color= 0
endif


;---------------------------------------------------------





end



;================================================================================



