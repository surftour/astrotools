; =================================================
;
;  SFR Properties
;
; ==================================================
;
pro get_sfrps_fruns, fruns, subdir=subdir, sfrps

nfruns= n_elements(fruns)

sfrps= fltarr(nfruns)

for i= 0, nfruns-1 do begin

        frun= '/raid4/tcox/'+subdir+'/'+fruns[i]

	;  get sfr properties
	; ---------------------
	fload_sfr_data, frun, spSFRmax, spSFRavg, aSFRmax, aSFRavg

        ;print, fruns[i], aSFRmax, spSFRavg

        sfrps[i]= aSFRmax/spSFRavg
endfor


end




; -------------------------------------------------------------------------



function get_sfrps, fruns, subdir=subdir

nfruns= n_elements(fruns)

sfrps= fltarr(nfruns)
    
for i= 0, nfruns-1 do begin

        frun= '/raid4/tcox/'+subdir+'/'+fruns[i]

        ;  get sfr properties
        ; ---------------------
        fload_sfr_data, frun, spSFRmax, spSFRavg, aSFRmax, aSFRavg

        print, fruns[i], aSFRmax, spSFRavg
    
        sfrps[i]= aSFRmax/spSFRavg
endfor
    
return, sfrps


end




; -------------------------------------------------------------------------
;
; spSFRavg= spiral phase average SFR
; aSFRmax= active phase maximum SFR
;
;
pro fload_sfr_data, frun, spSFRmax, spSFRavg, aSFRmax, aSFRavg, $
			isolated=isolated

    ;   get SFR rate from txt
    ;---------------------------
    open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass

    ;---------------------------
    if keyword_set(isolated) then begin
	atmin= 0.1
	atmax= max(sfrtime)
    endif else begin

	t_merg= fload_bh_mergertime(frun)
	if t_merg le 0.0 then t_merg= max(sfrtime)*0.5

	; define active phase
	atmin= t_merg-0.3
	atmax= atmin+1.0
    endelse

    idx= where(sfrtime lt atmin)
    spsfr= sfrsfr(idx)

    idx= where((sfrtime ge atmin) and (sfrtime lt atmax))
    asfr= sfrsfr(idx)

    ;---------------------------
    spSFRmax = max(spsfr)
    spSFRavg= mean(spsfr)

    aSFRmax= max(asfr)
    aSFRavg= mean(asfr)

end



; --------------------------------------------------






;========================================================================================
;
;
;
;          --------------
;          |            |
; SFR      |            |
; property |            |
;          |            |
;          |            |
;          --------------
;              Mass
;
;
pro sfrps1, junk


if not keyword_set(junk) then begin
	print, " "
	print, " sfrps1, junk"
	print, " "
	print, "    .run sfr_multi "
	print, " "
	print, " "
	return
endif

filename='sfrps.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



xaxistitle = "!6Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10

yaxistitle = "!6SFR!Dstarburst!N / SFR!Dprior!N "
ymax = 250.0
;ymax = 25.0
;ymax = 15.0
;ymax = 10.0
;ymin = 7.0
;ymin = 6.0
ymin = 2.5e-2


; ---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]

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
        ytickformat='exp_label', $
        xtickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; ---------------------------


; 5% gas, e orbit
;--------------------
fruns= ['z0e', 'z1e', 'z2e', 'z3e', 'z4e', 'z5e', 'z6e']
stmass_5e= [8.1668, 23.81, 70.727, 190.48, 529.709, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
sfrps_5e= get_sfrps(fruns, subdir='zs')
pointt= 7
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_5e, sfrps_5e, psym=-thispsym, color=thiscolor, thick=3.0


; 20% gas, e orbit
;--------------------
fruns= ['e0e', 'e1e', 'e2e', 'e3e', 'e4e', 'e5e', 'e6e']
stmass_20e= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
sfrps_20e= get_sfrps(fruns, subdir='es')
pointt= 6
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_20e, sfrps_20e, psym=-thispsym, color=thiscolor, thick=3.0


; 40% gas, e orbit
;--------------------
fruns=["d0e2_q", "d1e2_q", "d2e2_q", "d3e7", "d4e2_q", "d5e2_q", "d6e2_q"]
stmass_40e= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
sfrps_40e= get_sfrps(fruns, subdir='ds')
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_40e, sfrps_40e, psym=-thispsym, color=thiscolor, thick=3.0
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl, '40e_q', size=1.2, color=thiscolor, /data, charthick=4.0


; 40% gas, h orbit
;--------------------
fruns= ['d0h2_q', 'd1h2_q', 'd2h2_q', 'd3h7', 'd4h2_q', 'd5h2_q', 'd6h2_q']
stmass_40h= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
sfrps_40h= get_sfrps(fruns, subdir='ds')
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_40h, sfrps_40h, psym=-thispsym, color=thiscolor, thick=3.0


; 80% gas, e orbit
;--------------------
fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
stmass_80e= [8.1668, 23.81, 70.727, 190.48, 529.71, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
sfrps_80e= get_sfrps(fruns, subdir='bs')
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_80e, sfrps_80e, psym=-thispsym, color=thiscolor, thick=3.0


; 80% gas, h orbit
;--------------------
fruns=['b0h', 'b1h', 'b2h', 'b3h', 'b4h', 'b5h', 'b6h']
stmass_80h= [8.1668, 23.81, 70.727, 190.48, 529.71, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
sfrps_80h= get_sfrps(fruns, subdir='bs')
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_80h, sfrps_80h, psym=-thispsym, color=thiscolor, thick=3.0






;  extras
; ---------
x=[xmin,xmax]
y=[1.0,1.0]
oplot, x, y, psym=-3, linestyle=2, color=0


;  done
; ---------
device, /close




end







;==========================================================================
;
;
;     1 line with error bars
;
;
;          --------------
;          |            |
; SFR      |            |
; property |            |
;          |            |
;          |            |
;          --------------
;              Mass
;
;
pro sfrps2, junk


if not keyword_set(junk) then begin
	print, " "
	print, " sfrps2, junk"
	print, " "
	print, "    .run sfr_multi "
	print, " "
	print, " "
	return
endif

filename='sfrps.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



xaxistitle = "!6Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10

;yaxistitle = "!6SFR!Dstarburst!N / SFR!Dprior!N "
yaxistitle = "Log (!6SFR!Dstarburst!N / SFR!Dprior!N)"
ymax = alog10(200.0)
;ymax = 15.0
;ymax = 10.0
;ymin = 7.0
;ymin = 6.0
ymin = alog10(2.5e-2)


; ---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]

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
        xtickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; ---------------------------------------------
totmass= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass

f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f)

;fruns= ['z0e', 'z1e', 'z2e', 'z3e', 'z4e', 'z5e', 'z6e']
;get_sfrps_fruns, fruns, subdir='zs', sfrps_5e
;;oplot, totmass, sfrps_5e, psym=8, color= 0, thick=1.0, symsize= 1.0
;oplot, totmass, alog10(sfrps_5e), psym=8, color= 0, thick=1.0, symsize= 1.0

;fruns= ['e0e', 'e1e', 'e2e', 'e3e', 'e4e', 'e5e', 'e6e']
;get_sfrps_fruns, fruns, subdir='es', sfrps_20e
;;oplot, totmass, sfrps_20e, psym=8, color= 0, thick=1.0, symsize= 1.0
;oplot, totmass, alog10(sfrps_20e), psym=8, color= 0, thick=1.0, symsize= 1.0

fruns=["d0e2_q", "d1e2_q", "d2e2_q", "d3e7", "d4e2_q", "d5e2_q", "d6e2_q"]
get_sfrps_fruns, fruns, subdir='ds', sfrps_40e
;oplot, totmass, sfrps_40e, psym=8, color= 0, thick=1.0, symsize= 1.0
oplot, totmass, alog10(sfrps_40e), psym=8, color= 0, thick=1.0, symsize= 1.0

;fruns= ['d0h2_q', 'd1h2_q', 'd2h2_q', 'd3h7', 'd4h2_q', 'd5h2_q', 'd6h2_q']
;get_sfrps_fruns, fruns, subdir='ds', sfrps_40h
;;oplot, totmass, sfrps_40h, psym=8, color= 0, thick=1.0, symsize= 1.0
;oplot, totmass, alog10(sfrps_40h), psym=8, color= 0, thick=1.0, symsize= 1.0

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
get_sfrps_fruns, fruns, subdir='bs', sfrps_80e
;oplot, totmass, sfrps_80e, psym=8, color= 0, thick=1.0, symsize= 1.0
oplot, totmass, alog10(sfrps_80e), psym=8, color= 0, thick=1.0, symsize= 1.0

;fruns=['b0h', 'b1h', 'b2h', 'b3h', 'b4h', 'b5h', 'b6h']
;get_sfrps_fruns, fruns, subdir='bs', sfrps_80h
;;oplot, totmass, sfrps_80h, psym=8, color= 0, thick=1.0, symsize= 1.0
;oplot, totmass, alog10(sfrps_80h), psym=8, color= 0, thick=1.0, symsize= 1.0


ni= n_elements(totmass)
sfrps_avg= fltarr(ni)
sfrps_err= fltarr(ni)
for i=0, ni-1 do begin
   ;x= [sfrps_40e[i], sfrps_40h[i], sfrps_80e[i], sfrps_80h[i], sfrps_20e[i], sfrps_5e[i]]
   x= [sfrps_40e[i], sfrps_80e[i]]
   x= alog10(x)
   sfrps_avg[i]= mean(x)
   sfrps_err[i]= sqrt(variance(x))
endfor

f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill

oplot, totmass, sfrps_avg, psym=-8, color= 0, thick=3.0, symsize= 1.5
oploterror, totmass, sfrps_avg, sfrps_err, errcolor=0, color=0, thick= 3.0, errthick= 3.0

ylbl= -0.7
xpt1= 1e+11
xpt2= 2e+11
xpt3= 3e+11
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-8, color=0, thick=3.0, symsize= 1.1
xyouts, xpt3, ylbl-0.01, 'no SB winds (wBH)', size=1.2, color=thiscolor, /data, charthick=4.0






; ---------------------------------------------
winddir= 'sb10_mass'
thiscolor= 100
thispsym= 2

fruns= ['d0e_no', 'd1e_no', 'd2e_no', 'd3e_no', 'd4e_no', 'd5e_no', 'd6e_no']
get_sfrps_fruns, fruns, subdir=winddir, sfrps_d
;oplot, totmass, sfrps_d, psym=thispsym, color= thiscolor, thick=1.0, symsize= 1.0
oplot, totmass, alog10(sfrps_d), psym=thispsym, color= thiscolor, thick=1.0, symsize= 1.0

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
get_sfrps_fruns, fruns, subdir=winddir, sfrps_b
;oplot, totmass, sfrps_b, psym=thispsym, color= thiscolor, thick=1.0, symsize= 1.0
oplot, totmass, alog10(sfrps_b), psym=thispsym, color= thiscolor, thick=1.0, symsize= 1.0


ni= n_elements(totmass)
sfrps_avg= fltarr(ni)
sfrps_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [sfrps_d[i], sfrps_b[i]]
   x= alog10(x)
   sfrps_avg[i]= mean(x)
   sfrps_err[i]= sqrt(variance(x))
endfor

oplot, totmass, sfrps_avg, psym=-thispsym, color= thiscolor, thick=3.0, symsize= 1.5
oploterror, totmass, sfrps_avg, sfrps_err, errcolor=thiscolor, color=thiscolor, thick= 3.0, errthick= 3.0

ylbl= -0.9
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0, symsize= 1.1
xyouts, xpt3, ylbl-0.01, 'sb10', size=1.2, color=thiscolor, /data, charthick=4.0




; ---------------------------------------------
winddir= 'sb8_mass'
thiscolor= 50
thispsym= 4

fruns= ['d0e_no', 'd1e_no', 'd2e_no', 'd3e_no', 'd4e_no', 'd5e_no', 'd6e_no']
get_sfrps_fruns, fruns, subdir=winddir, sfrps_d
;oplot, totmass, sfrps_d, psym=thispsym, color= thiscolor, thick=1.0, symsize= 1.0
oplot, totmass, alog10(sfrps_d), psym=thispsym, color= thiscolor, thick=1.0, symsize= 1.0

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
get_sfrps_fruns, fruns, subdir=winddir, sfrps_b
;oplot, totmass, sfrps_b, psym=thispsym, color= thiscolor, thick=1.0, symsize= 1.0
oplot, totmass, alog10(sfrps_b), psym=thispsym, color= thiscolor, thick=1.0, symsize= 1.0


ni= n_elements(totmass)
sfrps_avg= fltarr(ni)
sfrps_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [sfrps_d[i], sfrps_b[i]]
   x= alog10(x)
   sfrps_avg[i]= mean(x)
   sfrps_err[i]= sqrt(variance(x))
endfor

oplot, totmass, sfrps_avg, psym=-thispsym, color= thiscolor, thick=3.0, symsize= 1.5
oploterror, totmass, sfrps_avg, sfrps_err, errcolor=thiscolor, color=thiscolor, thick= 3.0, errthick= 3.0

ylbl= -1.1
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0, symsize= 1.1
xyouts, xpt3, ylbl-0.01, 'sb8', size=1.2, color=thiscolor, /data, charthick=4.0





; ---------------------------------------------
winddir= 'sb13_mass'
thiscolor= 150
thispsym= 6

fruns= ['d0e_no', 'd1e_no', 'd2e_no', 'd3e_no', 'd4e_no', 'd5e_no', 'd6e_no']
get_sfrps_fruns, fruns, subdir=winddir, sfrps_d
;oplot, totmass, sfrps_d, psym=thispsym, color= thiscolor, thick=1.0, symsize= 1.0
oplot, totmass, alog10(sfrps_d), psym=thispsym, color= thiscolor, thick=1.0, symsize= 1.0

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
get_sfrps_fruns, fruns, subdir=winddir, sfrps_b
;oplot, totmass, sfrps_b, psym=thispsym, color= thiscolor, thick=1.0, symsize= 1.0
oplot, totmass, alog10(sfrps_b), psym=thispsym, color= thiscolor, thick=1.0, symsize= 1.0


ni= n_elements(totmass)
sfrps_avg= fltarr(ni)
sfrps_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [sfrps_d[i], sfrps_b[i]]
   x= alog10(x)
   sfrps_avg[i]= mean(x)
   sfrps_err[i]= sqrt(variance(x))
endfor

oplot, totmass, sfrps_avg, psym=-thispsym, color= thiscolor, thick=3.0, symsize= 1.5
oploterror, totmass, sfrps_avg, sfrps_err, errcolor=thiscolor, color=thiscolor, thick= 3.0, errthick= 3.0

ylbl= -1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0, symsize= 1.1
xyouts, xpt3, ylbl-0.01, 'sb13', size=1.2, color=thiscolor, /data, charthick=4.0





; ---------------------------------------------






;  extras
; ---------
x=[xmin,xmax]
y=[1.0,1.0]
oplot, x, y, psym=-3, linestyle=2, color=0


;  done
; ---------
device, /close




end







;==========================================================================




