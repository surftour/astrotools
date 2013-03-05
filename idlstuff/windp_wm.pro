
; =================================================
;
;  Wind mass versus Galaxy Mass
;
;
;
;
; ==================================================



function get_wind_mass, fruns, subdir=subdir

nfruns= n_elements(fruns)

windmass= fltarr(nfruns)

for i= 0, nfruns-1 do begin

        frun= '/raid4/tcox/'+subdir+'/'+fruns[i]

        ;  get wind info
        ;------------------
        read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
        unboundmass= mass_egt[n_elements(time)-1]

        windmass[i]= 1.0e+10 * unboundmass / 0.7
endfor

return, windmass

end




;==========================================================================


function get_windfraction, fruns, subdir=subdir

nfruns= n_elements(fruns)

windf= fltarr(nfruns)

if strmid(fruns[0],0,1) eq 'b' then default= [0.267872, 0.780969, 2.31984, 6.24775, 17.3745, 49.9820, 190.666]
if strmid(fruns[0],0,1) eq 'd' then default= [0.133936, 0.390484, 1.15992, 3.12387, 8.12085, 24.9910, 95.3331]
;if subdir eq 'es' then default= [0.0669681, 0.195242, 0.579959, 1.56194, 4.06043, 12.4955, 47.6665]
;if subdir eq 'zs' then default= [0.0167420, 0.0488105, 0.144990, 0.390484, 1.08590, 3.12387, 11.9166]

for i= 0, nfruns-1 do begin

        frun= '/raid4/tcox/'+subdir+'/'+fruns[i]

        ;  get wind info
        ;------------------
        read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
        unboundmass= mass_egt[n_elements(time)-1]


        ;  get sfr info 
        ;------------------
        open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        originalmass= sfrgasmass[0]

        ; trap for zero original gas mass (i.e., if run isn't finished)
        if originalmass eq 10.0 then originalmass= default[i]

        print, fruns[i], unboundmass, originalmass

        windf[i]= unboundmass/originalmass
endfor

return, windf

end




;==========================================================================



function get_windmetal_mass, fruns, subdir=subdir

nfruns= n_elements(fruns)

windzmass= fltarr(nfruns)

for i= 0, nfruns-1 do begin

        frun= '/raid4/tcox/'+subdir+'/'+fruns[i]

        ;  get wind info
        ;------------------
        read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
        unboundmetals= mass_egt[n_elements(time)-1] * 0.02 * wind_z[n_elements(time)-1] * 1.0d+10 / 0.7 ; z is in solar

        ;  get sfr info
        ;------------------
	; already in msolar
	read_gasinfo_file, frun, Ttime, totgmass, totgzmass, totnsmass, totnszmass, totwindmass, totwindz
	unboundmetals2= totwindz
	;totalmetals= totgzmass
	;totalmetals= totgzmass
	;totalmetals= totgzmass + totnszmass


        windzmass[i]= unboundmetals2
endfor

return, windzmass

end







;==========================================================================



function get_windmetalf, fruns, subdir=subdir

nfruns= n_elements(fruns)

windf= fltarr(nfruns)

if strmid(fruns[0],0,1) eq 'b' then default= [0.267872, 0.780969, 2.31984, 6.24775, 17.3745, 49.9820, 190.666]
if strmid(fruns[0],0,1) eq 'd' then default= [0.133936, 0.390484, 1.15992, 3.12387, 8.12085, 24.9910, 95.3331]
;if subdir eq 'es' then default= [0.0669681, 0.195242, 0.579959, 1.56194, 4.06043, 12.4955, 47.6665]
;if subdir eq 'zs' then default= [0.0167420, 0.0488105, 0.144990, 0.390484, 1.08590, 3.12387, 11.9166]

for i= 0, nfruns-1 do begin

        frun= '/raid4/tcox/'+subdir+'/'+fruns[i]

        ;  get wind info
        ;------------------
        read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
        unboundmetals= mass_egt[n_elements(time)-1] * 0.02 * wind_z[n_elements(time)-1] * 1.0d+10 / 0.7 ; z is in solar


        ;  get sfr info 
        ;------------------
        ; already in msolar
        read_gasinfo_file, frun, Ttime, totgmass, totgzmass, totnsmass, totnszmass, totwindmass, totwindz
        unboundmetals2= totwindz
        ;totalmetals= totgzmass
        ;totalmetals= totgzmass
        totalmetals= totgzmass + totnszmass


        print, fruns[i], unboundmetals, unboundmetals2, totalmetals

        windf[i]= unboundmetals2/totalmetals
endfor

return, windf

end







;==========================================================================



;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro wind_eff_v1, junk


if not keyword_set(junk) then begin
	print, " "
	print, " wind_eff, junk"
	print, " "
	print, " requires:   .run time_wind"
	print, "             .run sfr_multi"
	print, " "
	return
endif

filename='wm.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



yaxistitle = "!6Wind Fraction"
ymax = 0.25
ymin = 0.0


xaxistitle = "!6Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10



;---------------------------

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









; 40% gas, e orbit
;--------------------
fruns=["d0e2_q", "d1e2_q", "d2e2_q", "d3e7", "d4e2_q", "d5e2_q", "d6e2_q"]
stmass_40e= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
wf_40e= get_windfraction(fruns, subdir='ds')
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_40e, wf_40e, psym=-thispsym, color=thiscolor, thick=3.0
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl, '40e_q', size=1.2, color=thiscolor, /data, charthick=4.0


; 40% gas, h orbit
;--------------------
fruns= ['d0h2_q', 'd1h2_q', 'd2h2_q', 'd3h7', 'd4h2_q', 'd5h2_q', 'd6h2_q']
stmass_40h= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
wf_40h= get_windfraction(fruns, subdir='ds')
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_40h, wf_40h, psym=-thispsym, color=thiscolor, thick=3.0


; 80% gas, e orbit
;--------------------
fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
stmass_80e= [8.1668, 23.81, 70.727, 190.48, 529.71, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
wf_80e= get_windfraction(fruns, subdir='bs')
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_80e, wf_80e, psym=-thispsym, color=thiscolor, thick=3.0


; 80% gas, h orbit
;--------------------
fruns=['b0h', 'b1h', 'b2h', 'b3h', 'b4h', 'b5h', 'b6h']
stmass_80h= [8.1668, 23.81, 70.727, 190.48, 529.71, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
wf_80h= get_windfraction(fruns, subdir='bs')
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_80h, wf_80h, psym=-thispsym, color=thiscolor, thick=3.0


; 20% gas, e orbit
;--------------------
fruns= ['e0e', 'e1e', 'e2e', 'e3e', 'e4e', 'e5e', 'e6e']
stmass_20e= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
wf_20e= get_windfraction(fruns, subdir='es')
pointt= 6
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_20e, wf_20e, psym=-thispsym, color=thiscolor, thick=3.0


; 5% gas, e orbit
;--------------------
fruns= ['z0e', 'z1e', 'z2e', 'z3e', 'z4e', 'z5e', 'z6e']
stmass_5e= [8.1668, 23.81, 70.727, 190.48, 529.709, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
wf_5e= get_windfraction(fruns, subdir='zs')
pointt= 7
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_5e, wf_5e, psym=-thispsym, color=thiscolor, thick=3.0






;  done
; ------
device, /close




end





;==========================================================================





;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro wind_eff_v2, junk


if not keyword_set(junk) then begin
	print, " "
	print, " wind_eff, junk"
	print, " "
	print, " requires:   .run time_wind"
	print, "             .run sfr_multi"
	print, " "
	return
endif

filename='wm.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



yaxistitle = "!6Wind Mass Fraction"
ymax = 0.88
ymin = 0.0


xaxistitle = "!6Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10



;---------------------------

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



xpt1= 4.0d+12
xpt2= 6.5d+12
xpt3= 9.0d+12





; BH + no sb winds
;--------------------
totmass= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass

fruns=["d0e2_q", "d1e2_q", "d2e2_q", "d3e7", "d4e2_q", "d5e2_q", "d6e2_q"]
wf_40e= get_windfraction(fruns, subdir='ds')

fruns= ['d0h2_q', 'd1h2_q', 'd2h2_q', 'd3h7', 'd4h2_q', 'd5h2_q', 'd6h2_q']
wf_40h= get_windfraction(fruns, subdir='ds')

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
wf_80e= get_windfraction(fruns, subdir='bs')

fruns=['b0h', 'b1h', 'b2h', 'b3h', 'b4h', 'b5h', 'b6h']
wf_80h= get_windfraction(fruns, subdir='bs')

fruns= ['e0e', 'e1e', 'e2e', 'e3e', 'e4e', 'e5e', 'e6e']
wf_20e= get_windfraction(fruns, subdir='es')

fruns= ['z0e', 'z1e', 'z2e', 'z3e', 'z4e', 'z5e', 'z6e']
wf_5e= get_windfraction(fruns, subdir='zs')

ni= n_elements(totmass)
wf_avg= fltarr(ni)
wf_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_40h[i], wf_80e[i], wf_80h[i], wf_20e[i], wf_5e[i]]
   wf_avg[i]= mean(x)
   wf_err[i]= sqrt(variance(x))
endfor

f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill 

oplot, totmass, wf_avg, psym=-8, color= 0, thick=3.0, symsize= 1.5
oploterror, totmass, wf_avg, wf_err, errcolor=0, color=0, thick= 3.0, errthick= 3.0

ylbl= 0.82
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-8, color=0, thick=3.0, symsize= 1.1
oplot, [xpt2], [ylbl], psym=8, color=0, thick=3.0, symsize= 1.7
xyouts, xpt3, ylbl-0.01, 'no SB winds', size=1.2, color=thiscolor, /data, charthick=4.0







; BH + sb10 winds
;--------------------

fruns=["d0e", "d1e", "d2e", "d3e", "d4e", "d5e", "d6e"]
wf_40e= get_windfraction(fruns, subdir='sb10_mass')

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
wf_80e= get_windfraction(fruns, subdir='sb10_mass')

ni= n_elements(totmass)
wf_avg= fltarr(ni)
wf_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg[i]= mean(x)
   wf_err[i]= sqrt(variance(x))
endfor 

oplot, totmass, wf_avg, psym=-2, color= 150, thick=3.0, symsize= 1.5
oploterror, totmass, wf_avg, wf_err, errcolor=150, color=150, thick= 3.0, errthick= 3.0

ylbl= 0.77
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-2, color=150, thick=3.0, symsize= 1.1
oplot, [xpt2], [ylbl], psym=2, color=150, thick=3.0, symsize= 1.7
xyouts, xpt3, ylbl-0.01, '!7g!6=0.5, !8v!6!DW!N=837', size=1.2, color=150, /data, charthick=4.0



; no BH + sb10 winds
;--------------------

fruns=["d0e_no", "d1e_no", "d2e_no", "d3e_no", "d4e_no", "d5e_no", "d6e_no"]
wf_40e= get_windfraction(fruns, subdir='sb10_mass')

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
wf_80e= get_windfraction(fruns, subdir='sb10_mass')

ni= n_elements(totmass)
wf_avg= fltarr(ni)
wf_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg[i]= mean(x)
   wf_err[i]= sqrt(variance(x))
endfor

oplot, totmass, wf_avg, psym=-2, color= 150, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wf_avg, wf_err, errcolor=150, color=150, thick= 3.0, errthick= 3.0, linestyle= 1







; BH + sb8 winds
;--------------------

fruns=["d0e", "d1e", "d2e", "d3e", "d4e", "d5e", "d6e"]
wf_40e= get_windfraction(fruns, subdir='sb8_mass')

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
wf_80e= get_windfraction(fruns, subdir='sb8_mass')

ni= n_elements(totmass)
wf_avg= fltarr(ni)
wf_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg[i]= mean(x)
   wf_err[i]= sqrt(variance(x))
endfor

oplot, totmass, wf_avg, psym=-6, color= 100, thick=3.0, symsize= 1.5 
oploterror, totmass, wf_avg, wf_err, errcolor=100, color=100, thick= 3.0, errthick= 3.0

ylbl= 0.72
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-6, color=100, thick=3.0, symsize= 1.1
oplot, [xpt2], [ylbl], psym=6, color=100, thick=3.0, symsize= 1.7
xyouts, xpt3, ylbl-0.01, '!7g!6=2.0, !8v!6!DW!N=837', size=1.2, color=100, /data, charthick=4.0



; no BH + sb8 winds
;--------------------

fruns=["d0e_no", "d1e_no", "d2e_no", "d3e_no", "d4e_no", "d5e_no", "d6e_no"]
wf_40e= get_windfraction(fruns, subdir='sb8_mass')

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
wf_80e= get_windfraction(fruns, subdir='sb8_mass')

ni= n_elements(totmass)
wf_avg= fltarr(ni)
wf_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg[i]= mean(x)
   wf_err[i]= sqrt(variance(x))
endfor

oplot, totmass, wf_avg, psym=-6, color= 100, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wf_avg, wf_err, errcolor=100, color=100, thick= 3.0, errthick= 3.0, linestyle= 1





; BH + sb13 winds
;--------------------
;fruns=["d0e", "d1e", "d2e", "d3e", "d4e", "d5e", "d6e"]
;wf_40e= get_windfraction(fruns, subdir='sb13_mass')

;fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
;wf_80e= get_windfraction(fruns, subdir='sb13_mass')
;wf_80e= wf_40e

;ni= n_elements(totmass)
;wf_avg= fltarr(ni)
;wf_err= fltarr(ni)
;for i=0, ni-1 do begin
;   x= [wf_40e[i], wf_80e[i]]
;   wf_avg[i]= mean(x)
;   wf_err[i]= sqrt(variance(x))
;endfor
;
;oplot, totmass, wf_avg, psym=-5, color= 50, thick=3.0, symsize= 1.5
;oploterror, totmass, wf_avg, wf_err, errcolor=50, color=50, thick= 3.0, errthick= 3.0

ylbl= 0.67
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-5, color=50, thick=3.0, symsize= 1.1
oplot, [xpt2], [ylbl], psym=5, color=50, thick=3.0, symsize= 1.7
xyouts, xpt3, ylbl-0.01, '!7g!6=2.0, !8v!6!DW!N=209', size=1.2, color=50, /data, charthick=4.0



; no BH + sb13 winds
;--------------------

fruns=["d0e_no", "d1e_no", "d2e_no", "d3e_no", "d4e_no", "d5e_no", "d6e_no"]
wf_40e= get_windfraction(fruns, subdir='sb13_mass')

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
wf_80e= get_windfraction(fruns, subdir='sb13_mass')
;wf_80e= wf_40e

ni= n_elements(totmass)
wf_avg= fltarr(ni)
wf_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg[i]= mean(x)
   wf_err[i]= sqrt(variance(x))
endfor

oplot, totmass, wf_avg, psym=-5, color= 50, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wf_avg, wf_err, errcolor=50, color=50, thick= 3.0, errthick= 3.0, linestyle= 1



;-----------------------


xyouts, 1.6e+13, 0.56, 'BH', size=1.2, color=0, /data, charthick=4.0
oplot, [xpt2,1.2d+13], [0.57,0.57], psym=-3, color=0, thick=6.0, linestyle= 0

xyouts, 1.6e+13, 0.51, 'without BH', size=1.2, color=0, /data, charthick=4.0
oplot, [xpt2,1.2d+13], [0.52,0.52], psym=-3, color=0, thick=6.0, linestyle= 1


;-----------------------


;xyouts, 0.23, 0.90, 'w/ BHs', size=1.9, color=0, /normal, charthick=4.0



;  done
; ------
device, /close


print, "wrote to file: ", filename


end





;==========================================================================



;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro wind_eff_v3, junk


if not keyword_set(junk) then begin
	print, " "
	print, " wind_eff_v3, junk"
	print, " "
	print, " requires:   .run time_wind"
	print, "             .run sfr_multi"
	print, " "
	return
endif




; BH + no sb winds
;--------------------
totmass= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass

fruns=["d0e2_q", "d1e2_q", "d2e2_q", "d3e7", "d4e2_q", "d5e2_q", "d6e2_q"]
wf_40e= get_windfraction(fruns, subdir='ds')
wm_40e= get_wind_mass(fruns, subdir='ds')

fruns= ['d0h2_q', 'd1h2_q', 'd2h2_q', 'd3h7', 'd4h2_q', 'd5h2_q', 'd6h2_q']
wf_40h= get_windfraction(fruns, subdir='ds')
wm_40h= get_wind_mass(fruns, subdir='ds')

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
wf_80e= get_windfraction(fruns, subdir='bs')
wm_80e= get_wind_mass(fruns, subdir='bs')

fruns=['b0h', 'b1h', 'b2h', 'b3h', 'b4h', 'b5h', 'b6h']
wf_80h= get_windfraction(fruns, subdir='bs')
wm_80h= get_wind_mass(fruns, subdir='bs')

fruns= ['e0e', 'e1e', 'e2e', 'e3e', 'e4e', 'e5e', 'e6e']
wf_20e= get_windfraction(fruns, subdir='es')
wm_20e= get_wind_mass(fruns, subdir='es')

fruns= ['z0e', 'z1e', 'z2e', 'z3e', 'z4e', 'z5e', 'z6e']
wf_5e= get_windfraction(fruns, subdir='zs')
wm_5e= get_wind_mass(fruns, subdir='zs')

ni= n_elements(totmass)
wf_avg_std= fltarr(ni)
wf_err_std= fltarr(ni)
wm_avg_std= fltarr(ni)
wm_err_std= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_40h[i], wf_80e[i], wf_80h[i], wf_20e[i], wf_5e[i]]
   wf_avg_std[i]= mean(x)
   wf_err_std[i]= sqrt(variance(x))
   x= [wm_40e[i], wm_40h[i], wm_80e[i], wm_80h[i], wm_20e[i], wm_5e[i]]
   wm_avg_std[i]= mean(x)
   wm_err_std[i]= sqrt(variance(x))
endfor

; trap for large errors - we don't have large enough sample
idx= where(wm_err_std gt wm_avg_std)
if idx(0) ne -1 then wm_err_std(idx)= wm_avg_std(idx) * 0.75




; BH + sb10 winds
;--------------------

fruns=["d0e", "d1e", "d2e", "d3e", "d4e", "d5e", "d6e"]
wf_40e= get_windfraction(fruns, subdir='sb10_mass')
wm_40e= get_wind_mass(fruns, subdir='sb10_mass')

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
wf_80e= get_windfraction(fruns, subdir='sb10_mass')
wm_80e= get_wind_mass(fruns, subdir='sb10_mass')

ni= n_elements(totmass)
wf_avg_sb10= fltarr(ni)
wf_err_sb10= fltarr(ni)
wm_avg_sb10= fltarr(ni)
wm_err_sb10= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg_sb10[i]= mean(x)
   wf_err_sb10[i]= sqrt(variance(x))
   x= [wm_40e[i], wm_80e[i]]
   wm_avg_sb10[i]= mean(x)
   wm_err_sb10[i]= sqrt(variance(x))
endfor 

; trap for large errors - we don't have large enough sample
idx= where(wm_err_sb10 gt wm_avg_sb10)
if idx(0) ne -1 then wm_err_sb10(idx)= wm_avg_sb10(idx) * 0.75



; no BH + sb10 winds
;--------------------

fruns=["d0e_no", "d1e_no", "d2e_no", "d3e_no", "d4e_no", "d5e_no", "d6e_no"]
wf_40e= get_windfraction(fruns, subdir='sb10_mass')
wm_40e= get_wind_mass(fruns, subdir='sb10_mass')

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
wf_80e= get_windfraction(fruns, subdir='sb10_mass')
wm_80e= get_wind_mass(fruns, subdir='sb10_mass')

ni= n_elements(totmass)
wf_avg_sb10no= fltarr(ni)
wf_err_sb10no= fltarr(ni)
wm_avg_sb10no= fltarr(ni)
wm_err_sb10no= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg_sb10no[i]= mean(x)
   wf_err_sb10no[i]= sqrt(variance(x))
   x= [wm_40e[i], wm_80e[i]]
   wm_avg_sb10no[i]= mean(x)
   wm_err_sb10no[i]= sqrt(variance(x))
endfor






; BH + sb8 winds
;--------------------

fruns=["d0e", "d1e", "d2e", "d3e", "d4e", "d5e", "d6e"]
wf_40e= get_windfraction(fruns, subdir='sb8_mass')
wm_40e= get_wind_mass(fruns, subdir='sb8_mass')

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
wf_80e= get_windfraction(fruns, subdir='sb8_mass')
wm_80e= get_wind_mass(fruns, subdir='sb8_mass')

ni= n_elements(totmass)
wf_avg_sb8= fltarr(ni)
wf_err_sb8= fltarr(ni)
wm_avg_sb8= fltarr(ni)
wm_err_sb8= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg_sb8[i]= mean(x)
   wf_err_sb8[i]= sqrt(variance(x))
   x= [wm_40e[i], wm_80e[i]]
   wm_avg_sb8[i]= mean(x)
   wm_err_sb8[i]= sqrt(variance(x))
endfor




; no BH + sb8 winds
;--------------------

fruns=["d0e_no", "d1e_no", "d2e_no", "d3e_no", "d4e_no", "d5e_no", "d6e_no"]
wf_40e= get_windfraction(fruns, subdir='sb8_mass')
wm_40e= get_wind_mass(fruns, subdir='sb8_mass')

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
wf_80e= get_windfraction(fruns, subdir='sb8_mass')
wm_80e= get_wind_mass(fruns, subdir='sb8_mass')

ni= n_elements(totmass)
wf_avg_sb8no= fltarr(ni)
wf_err_sb8no= fltarr(ni)
wm_avg_sb8no= fltarr(ni)
wm_err_sb8no= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg_sb8no[i]= mean(x)
   wf_err_sb8no[i]= sqrt(variance(x))
   x= [wm_40e[i], wm_80e[i]]
   wm_avg_sb8no[i]= mean(x)
   wm_err_sb8no[i]= sqrt(variance(x))
endfor




; BH + sb13 winds
;--------------------
;fruns=["d0e", "d1e", "d2e", "d3e", "d4e", "d5e", "d6e"]
;wf_40e= get_windfraction(fruns, subdir='sb13_mass')

;fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
;wf_80e= get_windfraction(fruns, subdir='sb13_mass')
;wf_80e= wf_40e

;ni= n_elements(totmass)
;wf_avg= fltarr(ni)
;wf_err= fltarr(ni)
;for i=0, ni-1 do begin
;   x= [wf_40e[i], wf_80e[i]]
;   wf_avg[i]= mean(x)
;   wf_err[i]= sqrt(variance(x))
;endfor
;
;oplot, totmass, wf_avg, psym=-5, color= 50, thick=3.0, symsize= 1.5
;oploterror, totmass, wf_avg, wf_err, errcolor=50, color=50, thick= 3.0, errthick= 3.0


; no BH + sb13 winds
;--------------------

fruns=["d0e_no", "d1e_no", "d2e_no", "d3e_no", "d4e_no", "d5e_no", "d6e_no"]
wf_40e= get_windfraction(fruns, subdir='sb13_mass')
wm_40e= get_wind_mass(fruns, subdir='sb13_mass')

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
wf_80e= get_windfraction(fruns, subdir='sb13_mass')
wm_80e= get_wind_mass(fruns, subdir='sb13_mass')
;wf_80e= wf_40e

ni= n_elements(totmass)
wf_avg_sb13no= fltarr(ni)
wf_err_sb13no= fltarr(ni)
wm_avg_sb13no= fltarr(ni)
wm_err_sb13no= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg_sb13no[i]= mean(x)
   wf_err_sb13no[i]= sqrt(variance(x))
   x= [wm_40e[i], wm_80e[i]]
   wm_avg_sb13no[i]= mean(x)
   wm_err_sb13no[i]= sqrt(variance(x))
endfor




;-----------------------
;  Print the stuff
;-----------------------


;
;
;      |-----------------------|
;      |                       |
;      |                       |
;      |                       |
;      |     mass              |
;      |                       |
;      |                       |
;      |                       |
;      |-----------------------|
;      |                       |
;      |                       |
;      |                       |
;      |     fraction          |
;      |                       |
;      |                       |
;      |                       |
;      |-----------------------|
;
;
;
;

filename='wm3.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newysize= 20, newxsize= 12


;---------------------------

yaxistitle = "!6Wind Mass (M!D!9n!6!N)"
;ymax = 10.0
;ymin = 0.001
ymax = 1.0e+11
ymin = 3.0e+7

xaxistitle = "!6Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10

x0= 0.18
x1= 0.98

y0= 0.08
y1= 0.50
y2= 0.92

!p.position= [x0, y1, x1, y2]

;plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
;	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0,  charthick=3.0, /xlog, /ylog, $
;        xtickformat='(a1)', ytitle=yaxistitle, ytickformat='exp_label', /nodata
plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=9, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0,  charthick=3.0, /xlog, /ylog, $
        xtickformat='(a1)', ytitle=yaxistitle, ytickformat='exp_label', /nodata

;xaxistitle = "!6Stellar Mass (M!D!9n!6!N)"
;xmax = 1.6e+14 * 0.04
;xmin = 6.0e+10 * 0.04
;axis, xaxis=1, /ylog, xtickformat='exp_label', /noerase, color= 0, $
;		xrange=[xmin,xmax], xtitle=xaxistitle, xstyle= 1, $
;		xcharsize=1.50, charthick=3.0, xthick=4.0


;xpt1= 4.0d+12
;xpt2= 6.5d+12
;xyouts, 1.6e+13, 0.56, 'BH', size=1.2, color=0, /data, charthick=4.0
;oplot, [xpt2,1.2d+13], [0.57,0.57], psym=-3, color=0, thick=6.0, linestyle= 0

;xyouts, 1.6e+13, 0.51, 'without BH', size=1.2, color=0, /data, charthick=4.0
;oplot, [xpt2,1.2d+13], [0.52,0.52], psym=-3, color=0, thick=6.0, linestyle= 1


;-----------------------


f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill 

oplot, totmass, wm_avg_std, psym=-8, color= 0, thick=3.0, symsize= 1.5
oploterror, totmass, wm_avg_std, wm_err_std, errcolor=0, color=0, thick= 3.0, errthick= 3.0

;ylbl= 0.82
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-8, color=0, thick=3.0, symsize= 1.1
;oplot, [xpt2], [ylbl], psym=8, color=0, thick=3.0, symsize= 1.7
;xyouts, xpt3, ylbl-0.01, 'no SB winds', size=1.2, color=thiscolor, /data, charthick=4.0


;-----------------------
;  std

;oplot, totmass, wm_avg_std, psym=-2, color= 150, thick=3.0, symsize= 1.5
;oploterror, totmass, wm_avg_std, wm_err_std, errcolor=150, color=150, thick= 3.0, errthick= 3.0

;ylbl= 0.77
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-2, color=150, thick=3.0, symsize= 1.1
;oplot, [xpt2], [ylbl], psym=2, color=150, thick=3.0, symsize= 1.7
;xyouts, xpt3, ylbl-0.01, '!7g!6=0.5, !8v!6!DW!N=837', size=1.2, color=150, /data, charthick=4.0


;-----------------------
;  sb10

oplot, totmass, wm_avg_sb10, psym=-2, color= 150, thick=3.0, symsize= 1.5
oploterror, totmass, wm_avg_sb10, wm_err_sb10, errcolor=150, color=150, thick= 3.0, errthick= 3.0

;ylbl= 0.77
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-2, color=150, thick=3.0, symsize= 1.1
;oplot, [xpt2], [ylbl], psym=2, color=150, thick=3.0, symsize= 1.7
;xyouts, xpt3, ylbl-0.01, '!7g!6=0.5, !8v!6!DW!N=837', size=1.2, color=150, /data, charthick=4.0


; no bh
oplot, totmass, wm_avg_sb10no, psym=-2, color= 150, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wm_avg_sb10no, wm_err_sb10no, errcolor=150, color=150, thick= 3.0, errthick= 3.0, linestyle= 1



;-----------------------
;  sb8

oplot, totmass, wm_avg_sb8, psym=-6, color= 100, thick=3.0, symsize= 1.5 
oploterror, totmass, wm_avg_sb8, wm_err_sb8, errcolor=100, color=100, thick= 3.0, errthick= 3.0

;ylbl= 0.72
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-6, color=100, thick=3.0, symsize= 1.1
;oplot, [xpt2], [ylbl], psym=6, color=100, thick=3.0, symsize= 1.7
;xyouts, xpt3, ylbl-0.01, '!7g!6=2.0, !8v!6!DW!N=837', size=1.2, color=100, /data, charthick=4.0

oplot, totmass, wm_avg_sb8no, psym=-6, color= 100, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wm_avg_sb8no, wm_err_sb8no, errcolor=100, color=100, thick= 3.0, errthick= 3.0, linestyle= 1



;-----------------------
;  sb13

;ylbl= 0.67
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-5, color=50, thick=3.0, symsize= 1.1
;oplot, [xpt2], [ylbl], psym=5, color=50, thick=3.0, symsize= 1.7
;xyouts, xpt3, ylbl-0.01, '!7g!6=2.0, !8v!6!DW!N=209', size=1.2, color=50, /data, charthick=4.0

oplot, totmass, wm_avg_sb13no, psym=-5, color= 50, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wm_avg_sb13no, wm_err_sb13no, errcolor=50, color=50, thick= 3.0, errthick= 3.0, linestyle= 1



xaxistitle = "!6Stellar Mass (M!D!9n!6!N)"
xmax = 1.6e+14 * 0.04
xmin = 6.0e+10 * 0.04
axis, xaxis=1, /ylog, xtickformat='exp_label', /noerase, color= 0, $
		xrange=[xmin,xmax], xtitle=xaxistitle, xstyle= 1, $
		xcharsize=1.50, charthick=3.0, xthick=4.0


;---------------------------
;---------------------------
;---------------------------

yaxistitle = "!6Wind Mass Fraction"
ymax = 0.88
ymin = 0.0

xaxistitle = "!6Total Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10

totmass= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0,  charthick=3.0, /xlog, $
        xtickformat='exp_label', xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase


xpt1= 4.0d+12
xpt2= 6.5d+12
xpt3= 9.3d+12

xyouts, 1.6e+13, 0.56, 'BH', size=1.2, color=0, /data, charthick=4.0
oplot, [xpt2,1.2d+13], [0.57,0.57], psym=-3, color=0, thick=6.0, linestyle= 0

xyouts, 1.6e+13, 0.51, 'without BH', size=1.2, color=0, /data, charthick=4.0
oplot, [xpt2,1.2d+13], [0.52,0.52], psym=-3, color=0, thick=6.0, linestyle= 1


;-----------------------


f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill 

oplot, totmass, wf_avg_std, psym=-8, color= 0, thick=3.0, symsize= 1.5
oploterror, totmass, wf_avg_std, wf_err_std, errcolor=0, color=0, thick= 3.0, errthick= 3.0

ylbl= 0.82
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-8, color=0, thick=3.0, symsize= 1.1
oplot, [xpt2], [ylbl], psym=8, color=0, thick=3.0, symsize= 1.7
xyouts, xpt3, ylbl-0.01, 'no SB winds', size=1.2, color=thiscolor, /data, charthick=4.0


;-----------------------
;  std

oplot, totmass, wf_avg_std, psym=-2, color= 150, thick=3.0, symsize= 1.5
oploterror, totmass, wf_avg_std, wf_err_std, errcolor=150, color=150, thick= 3.0, errthick= 3.0

ylbl= 0.77
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-2, color=150, thick=3.0, symsize= 1.1
oplot, [xpt2], [ylbl], psym=2, color=150, thick=3.0, symsize= 1.7
xyouts, xpt3, ylbl-0.01, '!7g!6=0.5, !8v!6!DW!N=837', size=1.2, color=150, /data, charthick=4.0


;-----------------------
;  sb10

oplot, totmass, wf_avg_sb10, psym=-2, color= 150, thick=3.0, symsize= 1.5
oploterror, totmass, wf_avg_sb10, wf_err_sb10, errcolor=150, color=150, thick= 3.0, errthick= 3.0

ylbl= 0.77
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-2, color=150, thick=3.0, symsize= 1.1
oplot, [xpt2], [ylbl], psym=2, color=150, thick=3.0, symsize= 1.7
xyouts, xpt3, ylbl-0.01, '!7g!6=0.5, !8v!6!DW!N=837', size=1.2, color=150, /data, charthick=4.0


; no bh
oplot, totmass, wf_avg_sb10no, psym=-2, color= 150, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wf_avg_sb10no, wf_err_sb10no, errcolor=150, color=150, thick= 3.0, errthick= 3.0, linestyle= 1



;-----------------------
;  sb8

oplot, totmass, wf_avg_sb8, psym=-6, color= 100, thick=3.0, symsize= 1.5 
oploterror, totmass, wf_avg_sb8, wf_err_sb8, errcolor=100, color=100, thick= 3.0, errthick= 3.0

ylbl= 0.72
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-6, color=100, thick=3.0, symsize= 1.1
oplot, [xpt2], [ylbl], psym=6, color=100, thick=3.0, symsize= 1.7
xyouts, xpt3, ylbl-0.01, '!7g!6=2.0, !8v!6!DW!N=837', size=1.2, color=100, /data, charthick=4.0

oplot, totmass, wf_avg_sb8no, psym=-6, color= 100, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wf_avg_sb8no, wf_err_sb8no, errcolor=100, color=100, thick= 3.0, errthick= 3.0, linestyle= 1



;-----------------------
;  sb13

ylbl= 0.67
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-5, color=50, thick=3.0, symsize= 1.1
oplot, [xpt2], [ylbl], psym=5, color=50, thick=3.0, symsize= 1.7
xyouts, xpt3, ylbl-0.01, '!7g!6=2.0, !8v!6!DW!N=209', size=1.2, color=50, /data, charthick=4.0

oplot, totmass, wf_avg_sb13no, psym=-5, color= 50, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wf_avg_sb13no, wf_err_sb13no, errcolor=50, color=50, thick= 3.0, errthick= 3.0, linestyle= 1





;---------------------------



;  done
; ------
device, /close


print, "wrote to file: ", filename


end







;=================================================================








pro wind_z, junk


if not keyword_set(junk) then begin
	print, " "
	print, " wind_z, junk"
	print, " "
	print, " requires:   .run time_wind"
	print, " "
	return
endif

filename='wmz.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



yaxistitle = "!6Ejected Metal Fraction"
ymax = 0.88
ymin = 0.0


xaxistitle = "!6Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10



;---------------------------

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



xpt1= 4.0d+12
xpt2= 6.5d+12
xpt3= 9.0d+12





; BH + no sb winds
;--------------------
totmass= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass

fruns=["d0e2_q", "d1e2_q", "d2e2_q", "d3e7", "d4e2_q", "d5e2_q", "d6e2_q"]
wf_40e= get_windmetalf(fruns, subdir='ds')

;fruns= ['d0h2_q', 'd1h2_q', 'd2h2_q', 'd3h7', 'd4h2_q', 'd5h2_q', 'd6h2_q']
;wf_40h= get_windmetalf(fruns, subdir='ds')

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
wf_80e= get_windmetalf(fruns, subdir='bs')

;fruns=['b0h', 'b1h', 'b2h', 'b3h', 'b4h', 'b5h', 'b6h']
;wf_80h= get_windmetalf(fruns, subdir='bs')

;fruns= ['e0e', 'e1e', 'e2e', 'e3e', 'e4e', 'e5e', 'e6e']
;wf_20e= get_windmetalf(fruns, subdir='es')

;fruns= ['z0e', 'z1e', 'z2e', 'z3e', 'z4e', 'z5e', 'z6e']
;wf_5e= get_windmetalf(fruns, subdir='zs')

ni= n_elements(totmass)
wf_avg= fltarr(ni)
wf_err= fltarr(ni)
for i=0, ni-1 do begin
   ;x= [wf_40e[i], wf_40h[i], wf_80e[i], wf_80h[i], wf_20e[i], wf_5e[i]]
   x= [wf_40e[i], wf_80e[i]]
   wf_avg[i]= mean(x)
   wf_err[i]= sqrt(variance(x))
endfor

f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill 

oplot, totmass, wf_avg, psym=-8, color= 0, thick=3.0, symsize= 1.5
oploterror, totmass, wf_avg, wf_err, errcolor=0, color=0, thick= 3.0, errthick= 3.0

;ylbl= 0.82
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-8, color=0, thick=3.0, symsize= 1.1
;xyouts, xpt3, ylbl-0.01, 'no SB winds', size=1.2, color=thiscolor, /data, charthick=4.0







; BH + sb10 winds
;--------------------
fruns=["d0e", "d1e", "d2e", "d3e", "d4e", "d5e", "d6e"]
wf_40e= get_windmetalf(fruns, subdir='sb10_mass')

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
wf_80e= get_windmetalf(fruns, subdir='sb10_mass')

ni= n_elements(totmass)
wf_avg= fltarr(ni)
wf_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg[i]= mean(x)
   wf_err[i]= sqrt(variance(x))
endfor 

oplot, totmass, wf_avg, psym=-2, color= 150, thick=3.0, symsize= 1.5
oploterror, totmass, wf_avg, wf_err, errcolor=150, color=150, thick= 3.0, errthick= 3.0

;ylbl= 0.77
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-2, color=150, thick=3.0, symsize= 1.1
;xyouts, xpt3, ylbl-0.01, '!7g!6=0.5, !8v!6!DW!N=837', size=1.2, color=150, /data, charthick=4.0


; no BH + sb10 winds
;--------------------
fruns=["d0e_no", "d1e_no", "d2e_no", "d3e_no", "d4e_no", "d5e_no", "d6e_no"]
wf_40e= get_windmetalf(fruns, subdir='sb10_mass')

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
wf_80e= get_windmetalf(fruns, subdir='sb10_mass')

ni= n_elements(totmass)
wf_avg= fltarr(ni)
wf_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg[i]= mean(x)
   wf_err[i]= sqrt(variance(x))
endfor 

oplot, totmass, wf_avg, psym=-2, color= 150, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wf_avg, wf_err, errcolor=150, color=150, thick= 3.0, errthick= 3.0, linestyle= 1






; BH + sb8 winds
;--------------------

fruns=["d0e", "d1e", "d2e", "d3e", "d4e", "d5e", "d6e"]
wf_40e= get_windmetalf(fruns, subdir='sb8_mass')

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
wf_80e= get_windmetalf(fruns, subdir='sb8_mass')

ni= n_elements(totmass)
wf_avg= fltarr(ni)
wf_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg[i]= mean(x)
   wf_err[i]= sqrt(variance(x))
endfor

oplot, totmass, wf_avg, psym=-6, color= 100, thick=3.0, symsize= 1.5 
oploterror, totmass, wf_avg, wf_err, errcolor=100, color=100, thick= 3.0, errthick= 3.0

;ylbl= 0.72
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-6, color=100, thick=3.0, symsize= 1.1
;xyouts, xpt3, ylbl-0.01, '!7g!6=2.0, !8v!6!DW!N=837', size=1.2, color=100, /data, charthick=4.0



; no BH + sb8 winds
;--------------------
fruns=["d0e_no", "d1e_no", "d2e_no", "d3e_no", "d4e_no", "d5e_no", "d6e_no"]
wf_40e= get_windmetalf(fruns, subdir='sb8_mass')

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
wf_80e= get_windmetalf(fruns, subdir='sb8_mass')

ni= n_elements(totmass)
wf_avg= fltarr(ni)
wf_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg[i]= mean(x)
   wf_err[i]= sqrt(variance(x))
endfor

oplot, totmass, wf_avg, psym=-6, color= 100, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wf_avg, wf_err, errcolor=100, color=100, thick= 3.0, errthick= 3.0, linestyle= 1






; BH + sb13 winds
;--------------------

;fruns=["d0e", "d1e", "d2e", "d3e", "d4e", "d5e", "d6e"]
fruns=["d0e_no", "d1e_no", "d2e_no", "d3e_no", "d4e_no", "d5e_no", "d6e_no"]
wf_40e= get_windmetalf(fruns, subdir='sb13_mass')

;fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
wf_80e= get_windmetalf(fruns, subdir='sb13_mass')
;wf_80e= wf_40e

ni= n_elements(totmass)
wf_avg= fltarr(ni)
wf_err= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg[i]= mean(x)
   wf_err[i]= sqrt(variance(x))
endfor

oplot, totmass, wf_avg, psym=-5, color= 50, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wf_avg, wf_err, errcolor=50, color=50, thick= 3.0, errthick= 3.0, linestyle= 1

;ylbl= 0.67
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-5, color=50, thick=3.0, symsize= 1.1
;xyouts, xpt3, ylbl-0.01, '!7g!6=2.0, !8v!6!DW!N=209', size=1.2, color=50, /data, charthick=4.0




;-----------------------


filename='wmz.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


yaxistitle = "!6Ejected Metal Fraction"
ymax = 0.88
ymin = 0.0


xaxistitle = "!6Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10



;---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, color= 0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /xlog, $
        xtickformat='exp_label', xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase


xpt1= 4.0d+12
xpt2= 6.5d+12
xpt3= 9.0d+12



;-----------------------


f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill 

oplot, totmass, wf_avg, psym=-8, color= 0, thick=3.0, symsize= 1.5
oploterror, totmass, wf_avg, wf_err, errcolor=0, color=0, thick= 3.0, errthick= 3.0

; ----



;-----------------------
;xyouts, 0.23, 0.90, 'w/ BHs', size=1.9, color=0, /normal, charthick=4.0



;  done
; ------
device, /close



print, "wrote to file: ", filename



end


;=================================================================








pro wz3, junk


if not keyword_set(junk) then begin
	print, " "
	print, " wz3, junk"
	print, " "
	print, " requires:   .run time_wind"
	print, " "
	return
endif



; BH + no sb winds
;--------------------
totmass= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass

fruns=["d0e2_q", "d1e2_q", "d2e2_q", "d3e7", "d4e2_q", "d5e2_q", "d6e2_q"]
wf_40e= get_windmetalf(fruns, subdir='ds')
wm_40e= get_windmetal_mass(fruns, subdir='ds')

;fruns= ['d0h2_q', 'd1h2_q', 'd2h2_q', 'd3h7', 'd4h2_q', 'd5h2_q', 'd6h2_q']
;wf_40h= get_windmetalf(fruns, subdir='ds')

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
wf_80e= get_windmetalf(fruns, subdir='bs')
wm_80e= get_windmetal_mass(fruns, subdir='bs')

;fruns=['b0h', 'b1h', 'b2h', 'b3h', 'b4h', 'b5h', 'b6h']
;wf_80h= get_windmetalf(fruns, subdir='bs')

;fruns= ['e0e', 'e1e', 'e2e', 'e3e', 'e4e', 'e5e', 'e6e']
;wf_20e= get_windmetalf(fruns, subdir='es')

;fruns= ['z0e', 'z1e', 'z2e', 'z3e', 'z4e', 'z5e', 'z6e']
;wf_5e= get_windmetalf(fruns, subdir='zs')

ni= n_elements(totmass)
wf_avg_std= fltarr(ni)
wf_err_std= fltarr(ni)
wm_avg_std= fltarr(ni)
wm_err_std= fltarr(ni)
for i=0, ni-1 do begin
   ;x= [wf_40e[i], wf_40h[i], wf_80e[i], wf_80h[i], wf_20e[i], wf_5e[i]]
   x= [wf_40e[i], wf_80e[i]]
   wf_avg_std[i]= mean(x)
   wf_err_std[i]= sqrt(variance(x))
   x= [wm_40e[i], wm_80e[i]]
   wm_avg_std[i]= mean(x)
   wm_err_std[i]= sqrt(variance(x))
endfor

; trap for large errors - we don't have large enough sample
idx= where(wm_err_std gt wm_avg_std)
if idx(0) ne -1 then wm_err_std(idx)= wm_avg_std(idx) * 0.75





; BH + sb10 winds
;--------------------
fruns=["d0e", "d1e", "d2e", "d3e", "d4e", "d5e", "d6e"]
wf_40e= get_windmetalf(fruns, subdir='sb10_mass')
wm_40e= get_windmetal_mass(fruns, subdir='sb10_mass')

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
wf_80e= get_windmetalf(fruns, subdir='sb10_mass')
wm_80e= get_windmetal_mass(fruns, subdir='sb10_mass')

ni= n_elements(totmass)
wf_avg_sb10= fltarr(ni)
wf_err_sb10= fltarr(ni)
wm_avg_sb10= fltarr(ni)
wm_err_sb10= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg_sb10[i]= mean(x)
   wf_err_sb10[i]= sqrt(variance(x))
   x= [wm_40e[i], wm_80e[i]]
   wm_avg_sb10[i]= mean(x)
   wm_err_sb10[i]= sqrt(variance(x))
endfor 


; trap for large errors - we don't have large enough sample
idx= where(wm_err_sb10 gt wm_avg_sb10)
if idx(0) ne -1 then wm_err_sb10(idx)= wm_avg_sb10(idx) * 0.75







; no BH + sb10 winds
;--------------------
fruns=["d0e_no", "d1e_no", "d2e_no", "d3e_no", "d4e_no", "d5e_no", "d6e_no"]
wf_40e= get_windmetalf(fruns, subdir='sb10_mass')
wm_40e= get_windmetal_mass(fruns, subdir='sb10_mass')

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
wf_80e= get_windmetalf(fruns, subdir='sb10_mass')
wm_80e= get_windmetal_mass(fruns, subdir='sb10_mass')

ni= n_elements(totmass)
wf_avg_sb10no= fltarr(ni)
wf_err_sb10no= fltarr(ni)
wm_avg_sb10no= fltarr(ni)
wm_err_sb10no= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg_sb10no[i]= mean(x)
   wf_err_sb10no[i]= sqrt(variance(x))
   x= [wm_40e[i], wm_80e[i]]
   wm_avg_sb10no[i]= mean(x)
   wm_err_sb10no[i]= sqrt(variance(x))
endfor 





; BH + sb8 winds
;--------------------

fruns=["d0e", "d1e", "d2e", "d3e", "d4e", "d5e", "d6e"]
wf_40e= get_windmetalf(fruns, subdir='sb8_mass')
wm_40e= get_windmetal_mass(fruns, subdir='sb8_mass')

fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
wf_80e= get_windmetalf(fruns, subdir='sb8_mass')
wm_80e= get_windmetal_mass(fruns, subdir='sb8_mass')

ni= n_elements(totmass)
wf_avg_sb8= fltarr(ni)
wf_err_sb8= fltarr(ni)
wm_avg_sb8= fltarr(ni)
wm_err_sb8= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg_sb8[i]= mean(x)
   wf_err_sb8[i]= sqrt(variance(x))
   x= [wm_40e[i], wm_80e[i]]
   wm_avg_sb8[i]= mean(x)
   wm_err_sb8[i]= sqrt(variance(x))
endfor




; no BH + sb8 winds
;--------------------
fruns=["d0e_no", "d1e_no", "d2e_no", "d3e_no", "d4e_no", "d5e_no", "d6e_no"]
wf_40e= get_windmetalf(fruns, subdir='sb8_mass')
wm_40e= get_windmetal_mass(fruns, subdir='sb8_mass')

fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
wf_80e= get_windmetalf(fruns, subdir='sb8_mass')
wm_80e= get_windmetal_mass(fruns, subdir='sb8_mass')

ni= n_elements(totmass)
wf_avg_sb8no= fltarr(ni)
wf_err_sb8no= fltarr(ni)
wm_avg_sb8no= fltarr(ni)
wm_err_sb8no= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg_sb8no[i]= mean(x)
   wf_err_sb8no[i]= sqrt(variance(x))
   x= [wm_40e[i], wm_80e[i]]
   wm_avg_sb8no[i]= mean(x)
   wm_err_sb8no[i]= sqrt(variance(x))
endfor





; BH + sb13 winds
;--------------------

;fruns=["d0e", "d1e", "d2e", "d3e", "d4e", "d5e", "d6e"]
fruns=["d0e_no", "d1e_no", "d2e_no", "d3e_no", "d4e_no", "d5e_no", "d6e_no"]
wf_40e= get_windmetalf(fruns, subdir='sb13_mass')
wm_40e= get_windmetal_mass(fruns, subdir='sb13_mass')

;fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
fruns= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
wf_80e= get_windmetalf(fruns, subdir='sb13_mass')
wm_80e= get_windmetal_mass(fruns, subdir='sb13_mass')
;wf_80e= wf_40e

ni= n_elements(totmass)
wf_avg_sb13no= fltarr(ni)
wf_err_sb13no= fltarr(ni)
wm_avg_sb13no= fltarr(ni)
wm_err_sb13no= fltarr(ni)
for i=0, ni-1 do begin
   x= [wf_40e[i], wf_80e[i]]
   wf_avg_sb13no[i]= mean(x)
   wf_err_sb13no[i]= sqrt(variance(x))
   x= [wm_40e[i], wm_80e[i]]
   wm_avg_sb13no[i]= mean(x)
   wm_err_sb13no[i]= sqrt(variance(x))
endfor

; trap for large errors - we don't have large enough sample
idx= where(wm_err_sb13no gt wm_avg_sb13no)
if idx(0) ne -1 then wm_err_sb13no(idx)= wm_avg_sb13no(idx) * 0.75



;-----------------------


filename='wmz.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newysize= 20, newxsize= 12

;---------------------------

yaxistitle = "!6Ejected Metal Mass (M!D!9n!6!N)"
;ymax = 10.0
;ymin = 0.001
ymax = 1.0e+9
ymin = 3.0e+4

xaxistitle = "!6Total Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10

x0= 0.18
x1= 0.98

y0= 0.08
y1= 0.50
y2= 0.92

!p.position= [x0, y1, x1, y2]

;plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
;       xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0,  charthick=3.0, /xlog, /ylog, $
;        xtickformat='(a1)', ytitle=yaxistitle, ytickformat='exp_label', /nodata
plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=9, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0,  charthick=3.0, /xlog, /ylog, $
        xtickformat='(a1)', ytitle=yaxistitle, ytickformat='exp_label', /nodata


;-----------------------

f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill
oplot, totmass, wm_avg_std, psym=-8, color= 0, thick=3.0, symsize= 1.5
oploterror, totmass, wm_avg_std, wm_err_std, errcolor=0, color=0, thick= 3.0, errthick= 3.0

; ----

oplot, totmass, wm_avg_sb10, psym=-2, color= 150, thick=3.0, symsize= 1.5
oploterror, totmass, wm_avg_sb10, wm_err_sb10, errcolor=150, color=150, thick= 3.0, errthick= 3.0

; ----

oplot, totmass, wm_avg_sb10no, psym=-2, color= 150, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wm_avg_sb10no, wm_err_sb10no, errcolor=150, color=150, thick= 3.0, errthick= 3.0, linestyle= 1

; ----

oplot, totmass, wm_avg_sb8, psym=-6, color= 100, thick=3.0, symsize= 1.5
oploterror, totmass, wm_avg_sb8, wm_err_sb8, errcolor=100, color=100, thick= 3.0, errthick= 3.0

; ----

oplot, totmass, wm_avg_sb8no, psym=-6, color= 100, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wm_avg_sb8no, wm_err_sb8no, errcolor=100, color=100, thick= 3.0, errthick= 3.0, linestyle= 1

; ----

oplot, totmass, wm_avg_sb13no, psym=-5, color= 50, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wm_avg_sb13no, wm_err_sb13no, errcolor=50, color=50, thick= 3.0, errthick= 3.0, linestyle= 1





xaxistitle = "!6Stellar Mass (M!D!9n!6!N)"
xmax = 1.6e+14 * 0.04
xmin = 6.0e+10 * 0.04
axis, xaxis=1, /ylog, xtickformat='exp_label', /noerase, color= 0, $
                xrange=[xmin,xmax], xtitle=xaxistitle, xstyle= 1, $
                xcharsize=1.50, charthick=3.0, xthick=4.0


;---------------------------
;---------------------------
;---------------------------



; ---------------------------------------------------------------------------------------------


yaxistitle = "!6Ejected Metal Fraction"
ymax = 0.88
ymin = 0.0


xaxistitle = "!6Total Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10



;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, color= 0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /xlog, $
        xtickformat='exp_label', xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase


xpt1= 4.0d+12
xpt2= 6.5d+12
xpt3= 9.0d+12



;-----------------------

f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill 
oplot, totmass, wf_avg_std, psym=-8, color= 0, thick=3.0, symsize= 1.5
oploterror, totmass, wf_avg_std, wf_err_std, errcolor=0, color=0, thick= 3.0, errthick= 3.0

; ----

oplot, totmass, wf_avg_sb10, psym=-2, color= 150, thick=3.0, symsize= 1.5
oploterror, totmass, wf_avg_sb10, wf_err_sb10, errcolor=150, color=150, thick= 3.0, errthick= 3.0

; ----

oplot, totmass, wf_avg_sb10no, psym=-2, color= 150, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wf_avg_sb10no, wf_err_sb10no, errcolor=150, color=150, thick= 3.0, errthick= 3.0, linestyle= 1

; ----

oplot, totmass, wf_avg_sb8, psym=-6, color= 100, thick=3.0, symsize= 1.5 
oploterror, totmass, wf_avg_sb8, wf_err_sb8, errcolor=100, color=100, thick= 3.0, errthick= 3.0

; ----

oplot, totmass, wf_avg_sb8no, psym=-6, color= 100, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wf_avg_sb8no, wf_err_sb8no, errcolor=100, color=100, thick= 3.0, errthick= 3.0, linestyle= 1

; ----

oplot, totmass, wf_avg_sb13no, psym=-5, color= 50, thick=3.0, symsize= 1.5, linestyle= 1
oploterror, totmass, wf_avg_sb13no, wf_err_sb13no, errcolor=50, color=50, thick= 3.0, errthick= 3.0, linestyle= 1




;-----------------------
;xyouts, 0.23, 0.90, 'w/ BHs', size=1.9, color=0, /normal, charthick=4.0



;  done
; ------
device, /close



print, "wrote to file: ", filename



end



;==========================================================================





pro gas_info, frun

if not keyword_set(frun) then begin
	print, " "
	print, " gas_info, frun "
	print, " "
	return
endif

spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])
lastsnap= nsnaps-1

ok=fload_snapshot_bh(frun,lastsnap)

; what time is it?
; ---------------------
Ttime= float(fload_time(1))

; gas masses
; ---------------------
gmass= fload_gas_mass(1)
zmass= fload_gas_metallicity(1)
totgmass= total(gmass) * 1.0d+10 / 0.7
totgzmass= total(gmass * zmass) * 1.0d+10 / 0.7

; stellar masses
; ---------------------
nsmass= fload_newstars_mass(1)
nsz= fload_newstars_z(1)
totnsmass= total(nsmass) * 1.0d+10 / 0.7
totnszmass= total(nsmass * nsz) * 1.0d+10 / 0.7

; wind masses
; ---------------------
ke= fload_gas_energy(1,/kinetic)
pe= fload_gas_energy(1,/potential)
the= fload_gas_energy(1,/thermal)
e= ke + pe + the

; positive energy?
egt0_idx= where(e gt 0.0)
if egt0_idx(0) ne -1 then begin
	totwindmass= total(gmass(egt0_idx)) * 1.0d+10 / 0.7
	totwindz= total(gmass(egt0_idx) * zmass(egt0_idx)) * 1.0d+10 / 0.7
endif else begin
	totwindmass= 0.0
	totwindz= 0.0
endelse




; ---------------------
openw, 1, frun+'/gasinfo.txt', ERROR=err

printf, 1, "#   gasinfo.txt"
printf, 1, "# "
printf, 1, "#            <--------------      total masses in M_solar      --------------> "
printf, 1, "#   time                  Metals                Metals in       Wind      Metals"
printf, 1, "# (Gyr/h)        Gas      in Gas    New Stars   New Stars        Gas     in Wind"
printf, 1, FORMAT= '(F8.3,"   ",6(E10.3,"  "))', $
                Ttime, totgmass, totgzmass, totnsmass, totnszmass, totwindmass, totwindz
close, 1




end



;----------------------------------------------------------------



;  Read gasinfo.txt file
; ----------------------------
pro read_gasinfo_file, frun, Ttime, totgmass, totgzmass, totnsmass, totnszmass, totwindmass, totwindz

hgasfile= frun+'/gasinfo.txt'

;spawn, "wc "+hgasfile,result
;lines=long(result)
;lines=lines(0)-5
;if lines GT 0 then hgas_data= fltarr(9,lines)

openr, 1, hgasfile
junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, Ttime, totgmass, totgzmass, totnsmass, totnszmass, totwindmass,totwindz
close, 1


end




