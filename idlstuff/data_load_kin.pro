;================================================================================
;
;
;     V/sigma  Data
;  -------------------
;
;
;
;
;
;
;================================================================================





;========================================
;
;   Read the Bender data
;
;========================================

; This is the Bender (1988) data:
; 
; 
;
pro readandplot_bender_data, junk, addkey=addkey

; this first file has spaces in the name which 
; screws things up, use  the second one.
;osullivanfile= '/home/tcox/osullivan_01_all.txt'


; --------

benderfile= '/home/tcox/RotvAnisSupport/bender88_ellipticals.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=2, symsize=1.0, color= 0
; --------


benderfile= '/home/tcox/RotvAnisSupport/bender90_dwarfe_ellipticals.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=1, symsize=1.0, color= 0
; --------


benderfile= '/home/tcox/RotvAnisSupport/bender90_lowlum_ellipticals.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=1, symsize=1.0, color= 0
; --------


;---------------------------------------------------------
; normal coordinates
x0=0.30
y0=0.957

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






;========================================
;
;   Read the Davies data
;
;========================================

; This is the Davies et al. (1983) data:
; 
; 
;
pro readandplot_davies_data, junk, addkey=addkey

; this has the exact same file format
; as the bender data, so we'll use the 
; identical procedure
;
;


benderfile= '/home/tcox/RotvAnisSupport/davies83_ellipticals.txt'
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
benderfile= '/home/tcox/RotvAnisSupport/davies83_bulges.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=7, symsize=0.8, color= 0, thick=2.0
; --------


; ----------------------------------------------------------
; normal coordinates
x0=0.30
y0=0.957

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






pro read_bender_5, filename, ellip, mb, vmax_min, vmax_maj, sigma, $
				no_neg_trap=no_neg_trap


;
;  Format of this file is:
; -------------------------
;#   5x of these
;NGC205		0.40	-15.7		20	50    0.0
; etc...
;

spawn, "wc "+filename,result
lines=long(result)
datalines=lines(0)-5
ellip= fltarr(datalines)
mb= fltarr(datalines)
vmax_maj= fltarr(datalines)
vmax_min= fltarr(datalines)
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
        vmax_min(i)= float(tempjunk(3))
        vmax_maj(i)= float(tempjunk(4))
        sigma(i)= float(tempjunk(5))
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








;========================================
;
;   Read the deZeeuw data
;
;========================================

; This is the de Zeeuw et al. (2002) SAURON data:
; 
; 
;
pro readandplot_deZeeuw_data, junk, addkey=addkey

; this first file has spaces in the name which 
; screws things up, use  the second one.
;osullivanfile= '/home/tcox/osullivan_01_all.txt'

benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_cluster_es.txt'
read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

;benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_cluster_len.txt'
;read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
;oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

;benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_cluster_sp.txt'
;read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
;oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0


benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_field_es.txt'
read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

;benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_field_len.txt'
;read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
;oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

;benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_field_sp.txt'
;read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
;oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

; normal coordinates
x0=0.30
y0=0.957

; data coordinates
x0d= 0.05
y0d= 1.58

if keyword_set(addkey) then begin
	xyouts, x0+.01, y0-0.15, 'de Zeeuw (2002)', color= 0, size=0.8, charthick=2.0, /normal
	oplot, [x0d], [y0d-0.23], psym=5, color= 0
endif


;---------------------------------------------------------

end





;
; ----------------------------------------------------------
pro read_deZeeuw_data, filename, ellip, mb, vmax, sigma


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








;
; ----------------------------------------------------------
pro read_bbf_data, filename, mb, vsigst, a4diva


;
; # E's, cE's, dE's and bulges with known internal kinematics and/or Mg_2
; # ----------------------------------------------------------------------
; #
; # obj       l     b    Type     D   S  lgsig  S  lg re  SB_e     M_T   S  A_B   lg vds*   Mg2  S  B-Vo  S   a4   lg rc  R.C.
; # (1)      (2)   (3)    (4)    (5) (6)  (7)  (8)  (9)   (10)    (11) (12) (13)   (14)    (15)(16) (17)(18) (19)   (20)  (21)
; N0315   124.6 -32.5   LA    107.2 1  2.546  1  1.486  22.36  -23.61  1  0.26  -1.046   0.283 1  1.00  1 -0.31   9.999  3 - giant
; N0584   149.8 -67.6   E4     39.6 1  2.337  1  0.724  20.44  -21.72  1  0.13   0.190   0.283 1  0.93  1  1.50  -0.609  3   E's
;
;

spawn, "wc "+filename,result
lines=long(result)
;datalines=lines(0)-5
datalines=114               ; set this manually
ellip= fltarr(datalines)
mb= fltarr(datalines)
vmax= fltarr(datalines)
sigma= fltarr(datalines)
vsigst= fltarr(datalines)
a4diva= fltarr(datalines)

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

	if name eq 'stop' then break

        ;ellip(i)= float(tempjunk(x))
        mb(i)= float(tempjunk(10))
        ;vmax(i)= float(tempjunk(x))
        ;sigma(i)= float(tempjunk(x))
	vsigst(i)= 10^(float(tempjunk(13)))
	a4diva(i)= float(tempjunk(18))
endfor

close, 1


end







;
; ----------------------------------------------------------
pro read_bender_current, name, ellip, a4a, vmaj, vmin, sig, pK, dPA, M_t, SB_e, L_R, L_X



filename='/home/tcox/RotvAnisSupport/bender_current_data.txt'
headerlen= 6
;#  TABLE : e_kinxr.dat
;#  upper limits for velocities are given by negative values,
;#  no data available or a4 not classifiable: 999; if dPA = 999: galaxy
;#  re-classified as `face-on' SB0, if a4>4 then a4:=4;  pK = pec.core kinem., strong mi.ax.rot (no SB0)
;#   NGC 1-b/a   a4/a  Vmaj Vmin <sig> pK dPA   M_t   SB_e   L_R    L_X
;#                     km/s km/s km/s     deg   mag  mag/""  W/Hz   erg/s
;  N067A 0.53   3.70  999  999   999  0   1  -99.99 -99.99  99.99  99.99
;  N227  0.36   1.30  999  999   999  0   8  -21.90  21.33  99.99  99.99
;  N315  0.25  -0.31   30  -30   290  0   2  -23.61  22.36  24.31  41.93


spawn, "wc "+filename,result
lines=long(result)
datalines=lines(0)-headerlen
;datalines=114               ; set this manually

name= strarr(datalines)

ellip= fltarr(datalines)     ; 1
a4a= fltarr(datalines)       ; 2
vmaj= fltarr(datalines)      ; 3
vmin= fltarr(datalines)      ; 4
sig= fltarr(datalines)       ; 5
pK= fltarr(datalines)        ; 6
dPA= fltarr(datalines)       ; 7
M_t= fltarr(datalines)       ; 8
SB_e= fltarr(datalines)      ; 9
L_R= fltarr(datalines)       ; 10
L_X= fltarr(datalines)       ; 11

openr, 1, filename
print, "opening: ", filename
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,(lines(0)-headerlen-1) do begin
        readf, 1, junk
        ;print, junk
        tempjunk= strsplit(junk,/extract,count=count)

        name(i)= tempjunk(0)
	if name(i) eq 'stop' then break

        ellip(i)= float(tempjunk(1))
        a4a(i)= float(tempjunk(2))
        vmaj(i)= float(tempjunk(3))
        vmin(i)= float(tempjunk(4))
        sig(i)= float(tempjunk(5))
	pK(i)= float(tempjunk(6))
	dPA(i)= float(tempjunk(7))
	M_t(i)= float(tempjunk(8))
	SB_e(i)= float(tempjunk(9))
	L_R(i)= float(tempjunk(10))
	L_X(i)= float(tempjunk(11))
endfor

close, 1


end







