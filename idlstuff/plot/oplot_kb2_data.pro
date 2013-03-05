


;
; ----------------------------------------------------------
pro read_bender_current, name, ellip, a4a, vmaj, vmin, sig, pK, dPA, M_t, SB_e, L_R, L_X



filename='/n/home/tcox/Documents/Data_RotvAnisSupport/bender_current_data.txt'
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







pro oplot_kb2_data, junk

        ; --------------
        ; read b,s & g '94
        ;bbffile= '/home/tcox/RotvAnisSupport/bender94_minoraxis.txt'
        ;read_bender_5, bbffile, ellip, mb, vmin, vmaj, a4diva, /no_neg_trap

        ;mu= vmin / sqrt(vmin*vmin + vmaj*vmaj)

        ;idx= where(a4diva lt 3.5)
        ;if idx(0) ne -1 then begin
        ;       a4diva= a4diva(idx)
        ;       mu= mu(idx)
        ;endif

        ;f = (2.0*!pi/16.0)*findgen(17)
        ;usersym,0.7*cos(f),0.7*sin(f),/fill
        ;oplot, a4diva, mu, psym=8, color=0, thick=3.0



        ; -------------------------
        ; read bender current data
        read_bender_current, name, ellip, a4a, vmaj, vmin, sig, pK, dPA, M_t, SB_e, L_R, L_X

        ; make sure there is a major axis velocity
        idx= where(vmaj lt 800.0)
        if idx(0) ne -1 then begin
                vmaj= vmaj(idx)
                vmin= vmin(idx)
                a4a= a4a(idx)
        endif

        ; make sure there is a minor axis velocity
        idx= where(vmin lt 800.0)
        if idx(0) ne -1 then begin
                vmaj= vmaj(idx)
                vmin= vmin(idx)
                a4a= a4a(idx)
        endif

        ; ditch if both velocities are lower limits
        ;idx= where((vmin gt 0.0) and (vmaj gt 0.0))
        ;if idx(0) ne -1 then begin
        ;       vmaj= vmaj(idx)
        ;       vmin= vmin(idx)
        ;       a4a= a4a(idx)
        ;endif

        ; make lower limits for vmaj
        idx= where(vmaj lt 0.0)
        if idx(0) ne -1 then begin
           vmaj(idx)= -1.0 * vmaj(idx)
           ;for i=0,n_elements(idx)-1 do begin
        ;       vmaj(idx(i))= randomu(seed)*10.0
           ;endfor
        endif

        ; make lower limits for vmin
        idx= where(vmin lt 0.0)
        if idx(0) ne -1 then begin
            vmin(idx)= -1.0 * vmin(idx)
           ;for i=0,n_elements(idx)-1 do begin
        ;       vmin(idx(i))= -1.0*randomu(seed)*vmin(idx(i))
           ;endfor
        endif

        mu= vmin / sqrt(vmin*vmin + vmaj*vmaj)

        idx= where(a4a lt 50.0)
        if idx(0) ne -1 then begin
                a4a= a4a(idx)
                mu= mu(idx)
        endif

        f = (2.0*!pi/16.0)*findgen(17)
        usersym,0.7*cos(f),0.7*sin(f),/fill
        oplot, a4a, mu, psym=8, color=0, thick=3.0


end


