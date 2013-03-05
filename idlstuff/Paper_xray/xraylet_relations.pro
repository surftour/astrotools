


; ----------------------------
;  Read hotgas.txt file
; ----------------------------
pro read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, $
				gas_tot, gas_hot, gas_cold, gas_sf

hgasfile= frun+'/hotgas.txt'

spawn, "wc "+hgasfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then hgas_data= fltarr(9,lines)

openr, 1, hgasfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, hgas_data
close, 1


time= hgas_data[0,*]
temp_keV_X= hgas_data[1,*]
temp_keV= hgas_data[2,*]
temp_K= hgas_data[3,*]
entropy= hgas_data[4,*]
gas_tot= hgas_data[5,*]
gas_hot= hgas_data[6,*]
gas_cold= hgas_data[7,*]
gas_sf= hgas_data[8,*]


end




; ----------------------------
;  Read xrays.txt file
; ----------------------------
pro read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, $
				xray, xray_sf, $
				xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

hgasfile= frun+'/xrays.txt'

spawn, "wc "+hgasfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then hgas_data= fltarr(11,lines)

openr, 1, hgasfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, hgas_data
close, 1


time= hgas_data[0,*]
temp_keV_X= hgas_data[1,*]

; with new xrays.txt file (i.e. it has X-ray emission weighted Z)
z_X= hgas_data[2,*]
mass_xraygas= hgas_data[3,*]
entropy= hgas_data[4,*]
xray= hgas_data[5,*]
xray_sf= hgas_data[6,*]
xray_rs_s= hgas_data[7,*]
xray_rs_h= hgas_data[8,*]
xray_rs0_s= hgas_data[9,*]
xray_rs0_h= hgas_data[10,*]

end





; -------------------------------------------------------------------------------------------






;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------
;
;       L_x   -   T_x  
;
;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------

pro plot_Lx_Tx, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_Lx_Tx, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='LxTx.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
yaxistitle = "Log L!DX!N (ergs s!E-1!N)"
xaxistitle = "Temperature (keV)"
;xmax = 4.0
xmax = 3.0
;xmin = 0.10
xmin = 0.08
;xmin = 0.04
ymax = 43.5
;ymin = 38.6
ymin = 38.5


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
	/xlog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata





; vc3
;----------------------------
;readandplot_lx_and_else, "/raid4/tcox/vc3bvc3b", 3, /tx, yevolfac=2.0/3.0
;readandplot_lx_and_else, "/raid4/tcox/vc3bvc3b", 1, /tx, yevolfac=-1.0

;readandplot_lx_and_else, "/raid4/tcox/vc3bvc3b_no", 4, /tx, yevolfac=7.0/6.0
;readandplot_lx_and_else, "/raid4/tcox/vc3bvc3b_no", 3, /tx, yevolfac=-1.0

readandplot_lx_and_else, "/raid4/tcox/ds/d0e", 1, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d1e", 1, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d2e", 1, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d3e7", 1, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d4e", 1, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d5e", 1, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d6e", 1, /tx, yevolfac=-1

readandplot_lx_and_else, "/raid4/tcox/ds/d0h", 3, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d1h", 3, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d2h", 3, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d3h7", 3, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d4h", 3, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d5h", 3, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d6h", 3, /tx, yevolfac=-1

readandplot_lx_and_else, "/raid4/tcox/ds/d0k", 6, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d1k", 6, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d2k", 6, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d3k7", 6, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d4k", 6, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d5k", 6, /tx, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d6k", 6, /tx, yevolfac=-1






; extras
; ---------
;xyouts, 0.2, 0.92, 'bulge', /normal, charthick=3.0, size=1.33, color= 0


; Put actual data on there
; -------------------------
read_osullivan_03, loglx, loglb, tempx
;symsize= 0.2
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;oplot, tempx, loglx, psym=8, color= 0
oplot, tempx, loglx, psym=7, color= 0   ;, symsize=0.5


; a little key
;x0=0.06  &  x1= 0.28      ; first two for xmin/xmax= 0.04/4.0
;y0=42.6  &  y1= 43.0
x0=0.115  &  x1= 0.35       ;  for xmin/xmax = 0.08/3.0
y0=42.6  &  y1= 43.0
;xyouts, 0.7, 0.28, 'OFP (2001)', color= 0, charthick=2.0, /normal
;oplot, [x0+0.13], [y0+0.75], psym=8, color= 0
xyouts, 0.31, 0.83, 'OPC (2003)', color= 0, charthick=2.0, /normal
oplot, [x0+0.017], [y0+0.20], psym=7, color= 0
oplot, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], psym=-3, color=0



; Lx_T fit
; ----------
x= [0.01,100.0]
y= 42.45 + 4.8*(alog10(x))

oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0

; draw arrow for merger time
;timemerge=1.05/0.7
;arrow, timemerge, 36.8, timemerge, 37.5, COLOR=0, THICK=3.0, hthick=3.0, /data


; helpful info
; --------------

; bh seed mass
;xyouts, 0.13, 42.6, 'BH seed mass', color= 0, charthick=1.2, /data
;arrow, 0.13, 42.5, 0.13, 42.0, COLOR=0, THICK=3.0, hthick=3.0, /data


;xyouts, 0.22, 0.22, 'Without discrete source contribution', /normal, color= 0, size=1.4


device, /close




end






;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------
;
;       L_x   -   L_b  
;
;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------

pro plot_Lx_Lb, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_Lx_Lb, junk"
	print, " "
	print, " "
	return
endif

filename='LxLb.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
yaxistitle = "!6Log L!DX!N (ergs s!E-1!N)"
xaxistitle = "!6Log L!DB!N (L!D!9n!6!N)"
xmax = 12.5
xmin = 8.5
ymax = 43.5
ymin = 37.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
	;/xlog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;readandplot_lx_and_else, "/raid4/tcox/vc3bvc3b", 3, /lb, yevolfac=0.66
;readandplot_lx_and_else, "/raid4/tcox/vc3bvc3b", 1, /lb, yevolfac=-1

;readandplot_lx_and_else, "/raid4/tcox/vc3bvc3b_no", 4, /lb, yevolfac=7.0/6.0
;readandplot_lx_and_else, "/raid4/tcox/vc3bvc3b_no", 3, /lb, yevolfac=-1

readandplot_lx_and_else, "/raid4/tcox/ds/d0e", 1, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d1e", 1, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d2e", 1, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d3e7", 1, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d4e", 1, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d5e", 1, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d6e", 1, /lb, yevolfac=-1

readandplot_lx_and_else, "/raid4/tcox/ds/d0h", 3, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d1h", 3, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d2h", 3, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d3h7", 3, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d4h", 3, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d5h", 3, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d6h", 3, /lb, yevolfac=-1

readandplot_lx_and_else, "/raid4/tcox/ds/d0k", 6, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d1k", 6, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d2k", 6, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d3k7", 6, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d4k", 6, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d5k", 6, /lb, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d6k", 6, /lb, yevolfac=-1




; extras
; ---------
;xyouts, 0.2, 0.92, 'bulge', /normal, charthick=3.0, size=1.33, color= 0


; Put actual data on there
; -------------------------
read_osullivan_01, loglb, loglx, ttype
symsize= 0.2
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
oplot, loglb, loglx, psym=8, color= 0
;oplot, loglb, loglx, psym=7, color= 0, symsize= 0.5

read_osullivan_03, loglx, loglb, tempx
oplot, loglb, loglx, psym=7, color= 0


; a little key
x0=10.9  &  x1= 12.15
y0=37.4  &  y1= 38.4
xyouts, 0.7, 0.28, 'OFP (2001)', color= 0, charthick=2.0, /normal
oplot, [x0+0.13], [y0+0.75], psym=8, color= 0
xyouts, 0.7, 0.22, 'OPC (2003)', color= 0, charthick=2.0, /normal
oplot, [x0+0.13], [y0+0.25], psym=7, color= 0
oplot, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], psym=-3, color=0





; Lx_Lb fit
; ----------
; fit for all E's (OFP '01)
;x= [1.0,10.0,100.0]
;y= 17.98 + 2.17*x
;oplot, x, y, psym=-3, linestyle= 3, thick=2.0, color= 0

; bright E galaxy slope (OPC '03)
x= [1.0,10.0,100.0]
;y= 11.9 + 2.7*x     ; this is the OPC fit, but B lums are for h=0.5
y= 12.5 + 2.7*x     ; corrected for this plot because we use B lums from OFP, h=0.75
oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0



; discrete source Lx  (Ciotti et al. 1991)
x= [1.0,10.0,100.0]
y= 29.45 + x
oplot, x, y, psym=-3, linestyle= 1, thick=2.0, color= 0


; helpful info
oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0



; discrete source Lx  (Ciotti et al. 1991)
x= [1.0,10.0,100.0]
y= 29.45 + x
oplot, x, y, psym=-3, linestyle= 1, thick=2.0, color= 0


; helpful info
; --------------

; bh seed mass
;xyouts, 8.8, 42.6, 'BH seed mass', color= 0, charthick=1.2, /data
;arrow, 8.8, 42.5, 8.8, 42.0, COLOR=0, THICK=3.0, hthick=3.0, /data


;xyouts, 0.22, 0.22, 'Without discrete source contribution', /normal, color= 0, size=1.4

device, /close


end




;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------
;
;       Sigma   -   T_x  
;
;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------

pro plot_sigma_Tx, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_sigma_Tx, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='SigmaTx.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
yaxistitle = "!7r!3 (km s!E-1!N)"
xaxistitle = "Temperature (keV)"
xmax = 3.0
xmin = 0.08
ymax = 450.0
;ymin = 38.6
ymin = 90.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
	/xlog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;readandplot_lx_and_else, "/raid4/tcox/vc3bvc3b", 1, /sig, yevolfac=-1
;readandplot_lx_and_else, "/raid4/tcox/vc3bvc3b_no", 3, /sig, yevolfac=-1

readandplot_lx_and_else, "/raid4/tcox/ds/d0e", 1, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d1e", 1, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d2e", 1, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d3e7", 1, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d4e", 1, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d5e", 1, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d6e", 1, /sig, yevolfac=-1

readandplot_lx_and_else, "/raid4/tcox/ds/d0h", 3, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d1h", 3, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d2h", 3, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d3h7", 3, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d4h", 3, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d5h", 3, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d6h", 3, /sig, yevolfac=-1

readandplot_lx_and_else, "/raid4/tcox/ds/d0k", 6, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d1k", 6, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d2k", 6, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d3k7", 6, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d4k", 6, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d5k", 6, /sig, yevolfac=-1
readandplot_lx_and_else, "/raid4/tcox/ds/d6k", 6, /sig, yevolfac=-1





; extras
; ---------
;xyouts, 0.2, 0.92, 'bulge', /normal, charthick=3.0, size=1.33, color= 0


; Put actual data on there
; -------------------------
read_osullivan_03, loglx, loglb, tempx
read_osullivan_03b, sigma
;symsize= 0.2
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;oplot, tempx, loglx, psym=8, color= 0
oplot, tempx, sigma, psym=7, color= 0   ;, symsize=0.5


; a little key
x0=0.11  &  x1= 0.33
y0=335.0  &  y1= 390.0
;xyouts, 0.7, 0.28, 'OFP (2001)', color= 0, charthick=2.0, /normal
;oplot, [x0+0.13], [y0+0.75], psym=8, color= 0
xyouts, 0.31, 0.83, 'OPC (2003)', color= 0, charthick=2.0, /normal
oplot, [x0+0.020], [y0+25.0], psym=7, color= 0
oplot, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], psym=-3, color=0



; Sigma_T Beta_spec=1 (I think?)
; -------------------
x= [0.001,10.0]
y= 10^(0.480426*alog10(x) + alog10(300.0))

oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0



; helpful info
; --------------

; bh seed mass
;xyouts, 0.13, 42.6, 'BH seed mass', color= 0, charthick=1.2, /data
;arrow, 0.13, 42.5, 0.13, 42.0, COLOR=0, THICK=3.0, hthick=3.0, /data


;xyouts, 0.22, 0.22, 'Without discrete source contribution', /normal, color= 0, size=1.4


device, /close




end










; genereric plotting program
; ----------------------------
pro readandplot_lx_and_else, frun, pointselection, msg, y0, $
				tx=tx, $
				lb=lb, $
				sig=sig, $
				yevolfac=yevolfac

    thispsym= -3
    thiscolor= 150   ;red 

    ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
    read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

    ; get last index
    lstidx= n_elements(time)-1  ; take last

    ; default y-axis
    ;yval= xray[lstidx]
    yval= xray_rs_s[lstidx]

    ; T_x
    ; ----
    if keyword_set(tx) or keyword_set(sig) then begin
	xval= temp_keV_X[lstidx]
    endif

    ; L_b
    ; -----
    read_colors_file, frun, time,  mags
    lstidx= n_elements(time)-1
    Mb= mags[2,*]
    Lumb= -0.4*(Mb-5.51)

    if keyword_set(lb) then begin
	xval= Lumb[lstidx]
    endif

    ; Sigma
    ; ------
    if keyword_set(sig) then begin
	read_sigma_file, frun, time, Asigxy, Asigxz, Asigyz, Asigavg, Asigerr, $
                                Bsigxy, Bsigxz, Bsigyz, Bsigavg, Bsigerr, $
                                gsigavg, gsigerr
	slstidx= n_elements(time)-1
	yval= Asigavg[slstidx]
    endif

    ; add point source component to L_x
    ;L_x_discrete = 29.5 + Lumb[lstidx]
    ;yval= yval - 38.0
    ;L_x_discrete= L_x_discrete - 38.0
    ;yval= (10^(yval)) + (10^(L_x_discrete))
    ;yval= alog10(yval) + 38.0


;  open (or filled) square
; --------------------------
if pointselection eq 1 then begin
        symsize= 1.5
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=6.0
        symsel= 8
        symcolor= 50
endif

;  x's
; -----
if pointselection eq 2 then begin
        symsel= 7
        symcolor= 100
endif

;  open circle
; --------------
if pointselection eq 3 then begin
        symsize= 1.5
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        ;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
        symsel= 8
        symcolor= 150
endif

;  open triangle
; ----------------
if pointselection eq 4 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
        usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
        symsel= 8
        symcolor= 0
endif

;  asterisk
; ----------
if pointselection eq 5 then begin
        symsel= 2
        symcolor= 200
endif

;  filled square
; --------------------------
if pointselection eq 6 then begin
        symsize= 1.6
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 8
        symcolor= 100
endif

;  diamond
; ----------
if pointselection eq 7 then begin
        symsel= 4 
        symcolor= 120
endif

;  triangle (open)
; -----------------
if pointselection eq 8 then begin
        symsel= 5 
        symcolor= 180
endif

;  plus sign
; ----------
if pointselection eq 9 then begin
        symsel= 1 
        symcolor= 170
endif


print, xval, yval

    oplot, [xval], [yval], thick=3.0, psym=symsel, color= symcolor

    ;oplot, [2.2, 2.4], [41.48,41.48], thick=3.0, psym=-2, color= 50
    if not keyword_set(msg) then msg= ' '
    if not keyword_set(y0) then y0= 0.0
    xyouts, 0.73, y0, msg, /normal, charthick=3.0, size=1.33, color=thiscolor

if not keyword_set(sig) and yevolfac gt 0 then begin
   ; draw arrow for pt. source contribution
   L_x_discrete = 29.5 + Lumb[lstidx]
   yval_wpts= yval - 38.0
   L_x_discrete= L_x_discrete - 38.0
   yval_wpts= (10^(yval_wpts)) + (10^(L_x_discrete))
   yval_wpts= alog10(yval_wpts) + 38.0
   arrow, [xval], [yval], [xval], [yval_wpts], COLOR=100, THICK=3.0, hthick=3.0, /data

   yval=yval_wpts

   ; draw arrow for 5 Gyr evolution
   xevolfac=0.0
   if not keyword_set(yevolfac) then yevolfac= 0.0
   yevolfac=alog10(1+yevolfac)
   if keyword_set(lb) then xevolfac=-1.0*alog10(2.0)
   arrow, [xval], [yval], [xval+xevolfac], [yval+yevolfac], COLOR=symcolor, THICK=3.0, hthick=3.0, /data

endif


end




; genereric plotting program
;  - only this one does time
;    evolution
; ----------------------------
pro readandplot_lx_and_else_time, frun, msg, y0, $
				tx=tx, $
				lb=lb

    thispsym= -3
    thiscolor= 150   ;red 

    ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
    read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h


    ;yval= xray
    yval= xray_rs_s

    ; T_x
    ; ----
    if keyword_set(tx) then begin
	xval= temp_keV_X
    endif

    ; L_b
    ; -----
    read_colors_file, frun, time,  mags
    lstidx= n_elements(time)-1
    Mb= mags[2,*]
    Lumb= -0.4*(Mb-5.51)

    if keyword_set(lb) then begin
	xval= Lumb
    endif

    ; add point source component to L_x
    ; ------------------------------------
    ;L_x_discrete = 29.5 + Lumb
    ;yval= yval - 38.0
    ;L_x_discrete= L_x_discrete - 38.0
    ;yval= (10^(yval)) + (10^(L_x_discrete))
    ;yval= alog10(yval) + 38.0


    oplot, xval, yval, thick=3.0, psym=thispsym, color= thiscolor

    ;oplot, [2.2, 2.4], [41.48,41.48], thick=3.0, psym=-2, color= 50
    if not keyword_set(msg) then msg= ' '
    if not keyword_set(y0) then y0= 0.0
    xyouts, 0.73, y0, msg, /normal, charthick=3.0, size=1.33, color=thiscolor

end



; -----------------------------------------------------------------------------------




;=========================================================
;
;  Use Raymond & Smith ('77, I think) to determine
;  the X-ray luminosity.
;
;
;=========================================================

pro load_raymondsmith_xraylum, hard_xray_lum, soft_xray_lum, $
					zero_metallicity=zero_metallicity


ngas= long(fload_npart(0))
;Ttime= float(fload_time(1))
Redshift= double(0.0)

GasMetallicity = (1/0.02)*fload_gas_metallicity(1)
GasNe = fload_gas_nume(1)
GasMass = 1.0e+10*fload_gas_mass(1)  ; in solar masses
GasTemp = float(fload_gas_temperature(1))   ; in K (for some reason this comes back in double)
GasHsml = fload_gas_hsml(1)


rho= fload_gas_rho(1)

; units
;hub= 0.7
hub= 1.0
ProtonMass=         1.6726e-24
UnitLength_in_cm=   3.085678d+21
UnitMass_in_g=      1.989d+43
UnitDensity_in_cgs= UnitMass_in_g/(UnitLength_in_cm^3)
print, UnitDensity_in_cgs

; hydrogen number density (nH cm-3)
GasHIRho= rho*0.76*(hub*hub)*UnitDensity_in_cgs/ProtonMass
GasHIRho= float(GasHIRho)

; electron number density (ne cm-3)
GasNeRho= float(GasNe*GasHIRho)

idx=where((rho lt 0.000854924) and (GasTemp ge 1.0e+5))
if idx(0) ne -1 then begin
    print, n_elements(idx), " out of ",ngas," particles will have diffuse X-ray emission."
    ngas= n_elements(idx)
    GasMetallicity= GasMetallicity(idx)
    GasNeRho= GasNeRho(idx)
    GasHIRho= GasHIRho(idx)
    GasMass= GasMass(idx)
    GasTemp= GasTemp(idx)
    GasHsml= GasHsml(idx)
endif else begin
    ngas= 0
endelse


; for testing purposes
;    ngas= 10L
;    GasMetallicity= GasMetallicity(5600:5609)
;    GasNeRho= GasNeRho(5600:5609)
;    GasHIRho= GasHIRho(5600:5609)
;    GasMass= GasMass(5600:5609)
;    GasTemp= GasTemp(5600:5609)
;    GasHsml= GasHsml(5600:5609)
;

if keyword_set(zero_metallicity) then begin
	GasMetallicity(*)= 0.0
endif


;
; new incarnation of code doesn't
; like the zero metallicities
;
idx=where(GasMetallicity le 0.0)
if idx(0) ne -1 then begin
    GasMetallicity(idx)= 1.0e-5
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;       Find luminosities
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if(ngas gt 0) then begin

        soft_xray_lum = fltarr(ngas)
        hard_xray_lum = fltarr(ngas)

        ;S = CALL_EXTERNAL('/home/brant/code/idl/raymond_smith/raymond_smith', $
        S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/RaymondSmith/raymond_smith', $
                'raymond_smith', $
                Ngas, $
                Redshift, $
                GasMetallicity, $
                GasNeRho, $
		GasHIRho, $
                GasMass, $
                GasTemp, $
                GasHsml, $
                soft_xray_lum, $
                hard_xray_lum)

        ;LUMINOSITIES are in h^-1 units

endif else begin
        print,'No gas, no x-ray luminosity.'
	hard_xray_lum= [0]
	soft_xray_lum= [0]
endelse



; brant returns this in solar luminosities
; multiply by 3.989e33 to get ergs/sec


hard_xray_lum = hard_xray_lum*3.989d33
soft_xray_lum = soft_xray_lum*3.989d33
print, "Total hard= ",total(hard_xray_lum)," erg sec-1"
print, "Total soft= ",total(soft_xray_lum)," erg sec-1"


end












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




