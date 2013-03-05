;========================================
;
;   Read the deZeeuw data
;
;========================================



;
; ----------------------------------------------------------
pro read_dezeeuw_data, filename, ellip, mb, vmax, sigma


;
;  Format of this file is:
; -------------------------
;#  "Field"
;#  Elliptical, Lenticular, and Spiral data from de Zeeuw et al. (2002) - SAURON II
;#
;#Galaxy Type            T-type  Vsys    d_m     M_b     (B-V)   Mg2     R_e     e       mu_e    V_max   Sigma   Gamma   M_BH
;#                               (km/s)  (mag)   (mag)   (mag)   mag     (arcsec)        (mag)   (km/s)  (km/s)          a-b(10^c) Msol
;Ellipticals
;NGC821  E6?             4.2     1742    31.86   20.44   1.020   0.316   50      0.32    22.02   91      208     -1      3.0-7.0(7)
;NGC2699 E:              5.0     1825    31.83   18.85   0.980   0.282   -1      0.06    -1      -1      -1      -1      -1 
;NGC2768 E6:             3.1     1324    31.66   21.15   0.960   0.276   64      0.42    21.94   148     188     -1      -1 
;NGC2974 E4              3.6     1983    31.93   20.32   1.005   0.305   24      0.39    20.74   207     229     -1      -1 
; etc...
;#
;#
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
        ellip(i)= float(tempjunk(9))
        mb(i)= float(tempjunk(5))
        vmax(i)= float(tempjunk(11))
        sigma(i)= float(tempjunk(12))
endfor

close, 1


end










; This is the de Zeeuw et al. (2002) SAURON data:
; 
; 
;
pro oplot_dezeeuw_vsig_data, junk, addkey=addkey

; this first file has spaces in the name which 
; screws things up, use  the second one.
;osullivanfile= '/n/home/tcox/osullivan_01_all.txt'

benderfile= '/n/home/tcox/Documents/Data_RotvAnisSupport/deZeeuw02_cluster_es.txt'
read_dezeeuw_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

;benderfile= '/n/home/tcox/Documents/Data_RotvAnisSupport/deZeeuw02_cluster_len.txt'
;read_dezeeuw_data, benderfile, ellip, mb, vmax, sigma
;oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

;benderfile= '/n/home/tcox/Documents/Data_RotvAnisSupport/deZeeuw02_cluster_sp.txt'
;read_dezeeuw_data, benderfile, ellip, mb, vmax, sigma
;oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0


benderfile= '/n/home/tcox/Documents/Data_RotvAnisSupport/deZeeuw02_field_es.txt'
read_dezeeuw_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

;benderfile= '/n/home/tcox/Documents/Data_RotvAnisSupport/deZeeuw02_field_len.txt'
;read_dezeeuw_data, benderfile, ellip, mb, vmax, sigma
;oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

;benderfile= '/n/home/tcox/Documents/Data_RotvAnisSupport/deZeeuw02_field_sp.txt'
;read_dezeeuw_data, benderfile, ellip, mb, vmax, sigma
;oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

; normal coordinates
x0=0.25
y0=0.93

; data coordinates
x0d= 0.05
y0d= 1.58

if keyword_set(addkey) then begin
	xyouts, x0+.01, y0-0.15, 'de Zeeuw (2002)', color= 0, size=0.8, charthick=2.0, /normal
	oplot, [x0d], [y0d-0.23], psym=5, color= 0
endif


;---------------------------------------------------------

end



