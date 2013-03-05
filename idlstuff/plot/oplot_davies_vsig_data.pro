;========================================
;
;   Read the Davies data
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

; This is the Davies et al. (1983) data:
; 
; 
;
pro oplot_davies_vsig_data, junk, addkey=addkey

; this has the exact same file format
; as the bender data, so we'll use the 
; identical procedure
;
;


benderfile= '/n/home/tcox/Documents/Data_RotvAnisSupport/davies83_ellipticals.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma

; low luminosity ellipticals
symsize= 1.0
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
idx= where(mb gt -20.5)
oplot, ellip(idx), vmax(idx)/sigma(idx), psym=8, color= 0
; --------

; high luminosity ellipticals
symsize= 0.8
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
idx= where(mb le -20.5)
oplot, ellip(idx), vmax(idx)/sigma(idx), psym=8, color= 0
; --------


; bulges
benderfile= '/n/home/tcox/Documents/Data_RotvAnisSupport/davies83_bulges.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=7, symsize=0.8, color= 0, thick=2.0
; --------


; ----------------------------------------------------------
; normal coordinates
x0=0.25
y0=0.93

; data coordinates
x0d= 0.05
y0d= 1.61

if keyword_set(addkey) then begin
	xyouts, x0+.01, y0-0.03, 'Davies et al. (1983)', color= 0, size=0.8, charthick=2.0, /normal
	symsize=1.0
	usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	oplot, [x0d-0.025], [y0d], psym=8, color= 0
	xyouts, x0-.0615, y0-0.03, ',', color= 0, size=0.8, charthick=2.0, /normal
	symsize= 0.8
	usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
	oplot, [x0d], [y0d], psym=8, color= 0
	xyouts, x0-.0325, y0-0.03, ',', color= 0, size=0.8, charthick=2.0, /normal
	oplot, [x0d+0.025], [y0d], psym=7, color= 0,thick=2.0
endif
; ----------------------------------------------------------



end







;================================================================================


