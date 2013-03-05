
;==========================================
;
;    Feedback Energy
;
;
;
;
;==========================================






; --------------------------------
;  Read FB Energy File
; ----------------------------------
pro read_fbenergy_file, frun, time, fb_std, fb_wind, fb_bh, $
				fbfile=fbfile

if keyword_set(fbfile) then fbfile= frun+'/'+fbfile else fbfile= frun+'/fb_energy.txt'

spawn, "wc "+fbfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then re_data= fltarr(4,lines)

openr, 1, fbfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, re_data
close, 1


time= re_data[0,*]
fb_std= re_data[1,*]
fb_wind= re_data[2,*]
fb_bh= re_data[3,*]


end




;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------





;-------------------------------------------------
;-------------------------------------------------
pro plot_fbem, junk, filename=filename, t_merg=t_merg


if not keyword_set(junk) then begin
        print, " "
        print, " plot_fbem, junk, filename=filename"
        print, " "
        print, " "
        return
endif

if not keyword_set(filename) then filename='fbenergymass.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=24.0, newysize= 20.0

;---------------------------


	frun= junk
        ;frun= "/raid4/tcox/sbw/sb10BH"
        ;frun= "/raid4/tcox/ds/d0e2_q"
        ;frun= "/raid4/tcox/ds/d1e2_q"
        ;frun= "/raid4/tcox/ds/d2e2_q"
        ;frun= "/raid4/tcox/ds/d3e7"
        ;frun= "/raid4/tcox/ds/d4e2_q"
        ;frun= "/raid4/tcox/ds/d5e2_q"
        ;frun= "/raid4/tcox/ds/d6e2_q"
        read_fbenergy_file, frun, time, fb_std, fb_wind, fb_bh


time=time/0.7


;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"
xmax = 4.25
xmin = 0.0

yaxistitle = "!6Log FB Energy Rate (erg s!E-1!N)"
ymax = 44.8
ymin = 39.2


;---------------------------

x0= 0.08
x1= 0.44

y0= 0.08
y1= 0.53
y2= 0.98

;---------------------------

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtickformat='(a1)', $
        ;xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------

if not keyword_set(t_merg) then begin
	t_merg= fload_bh_mergertime(frun)
	t_merg= t_merg/0.7
endif
btmin= t_merg-0.3
btmax= btmin+1.0
;btmin= 1.4
;btmax= 2.4

timerange= [btmin, btmin, btmax, btmax, btmin]
erange= [ymin,ymax,ymax,ymin,ymin]
polyfill, timerange, erange, /data, color= 250, /fill

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        charthick=3.0, xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase


;---------------------------


        oplot, time, fb_wind, thick=4.0, color= 50, psym=-2, linestyle=0
        ;oplot, time, fb_std, thick=4.0, color= 50, psym=-3
        oplot, time, fb_bh, thick=12.0, color= 150, psym=-3, linestyle=1

	;xyouts, 0.68, 0.85, "Std", /normal, color= 10, charthick=3.0, charsize=1.3
        oplot, [2.8,3.0,3.2], [44.0,44.0,44.0], thick=4.0, color= 50, psym=-2, linestyle=0
	xyouts, 0.37, 0.91, "SF", /normal, color= 50, charthick=3.0, charsize=1.5

        oplot, [2.8,3.2], [43.4,43.4], thick=12.0, color= 150, psym=-3, linestyle=1
	xyouts, 0.37, 0.86, "BH", /normal, color= 150, charthick=3.0, charsize=1.5
	

;---------------------------




yaxistitle = "!6Log Integrated FB Energy (erg)"
ymax = 60.4
ymin = 56.8


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, /noerase, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------

;t_merg= fload_bh_mergertime(frun)
;t_merg= t_merg/0.7
;btmin= t_merg-0.3
;btmax= btmin+1.0
;btmin= 1.4
;btmax= 2.4

timerange= [btmin, btmin, btmax, btmax, btmin]
erange= [ymin,ymax,ymax,ymin,ymin]
polyfill, timerange, erange, /data, color= 250, /fill


plot, [1.0],[1.0], psym=-3, /noerase, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        charthick=3.0, xtitle=xaxistitle, ytitle=yaxistitle, /nodata

;---------------------------

	; add in a zero point for feedback energy
	; to start the simulation

	ok= fload_snapshot_bh(frun,0)
	
	; starburst (vc3, 40% gas sim)
	;seed_stellarmass= 4.68581d+10

	seed_stellarmass= 1.0d+10 * total(fload_allstars_mass(1)) / 0.7
	zero_level_sb= alog10(0.5 * 1.4d+49 * seed_stellarmass)
	zero_level_sb= zero_level_sb - 40.0

	; black hole
	;seed_bhmass= 2.0d+5 / 0.7

	;see_bhmass= 1.0d+10 * total(fload_blackhole_mass(1)) / 0.07

	open_blackhole_txt, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, $
                                bh_totalmass, bh_mdot_edd
        seed_bhmass= bh_mass[0] / 0.7    ; already has the 1d10 factor  

	cc= 2.9979d+10
	zero_level_bh= alog10(0.05 * 0.1 * cc * cc * seed_bhmass * 1.989d+33)
	zero_level_bh= zero_level_bh - 40.0


;---------------------------


	gadunits_in_sec= 3.08568d+16

	ntime= n_elements(time)
	int_fb_wind= fltarr(ntime)
	int_fb_std= fltarr(ntime)
	int_fb_bh= fltarr(ntime)

	; put in more managable units
	unitsof40= 40.0
	fb_wind= fb_wind - unitsof40
	fb_std= fb_std - unitsof40
	fb_bh= fb_bh - unitsof40

	tot_wind_e= 0.0
	tot_std_e= 0.0
	tot_bh_e= 0.0

	int_fb_wind[0]= zero_level_sb + 40.0
	int_fb_std[0]= zero_level_sb + 40.0
	int_fb_bh[0]= zero_level_bh + 40.0

	for i=1,ntime-1 do begin

		;dt= (time[i]-time[i-1]) * gadunits_in_sec
		dt= (time[i]-time[i-1])

		avg_wind_e= dt * 0.5 * ((10^fb_wind[i-1]) + (10^fb_wind[i]))
		avg_std_e= dt * 0.5 * ((10^fb_std[i-1]) + (10^fb_std[i]))
		avg_bh_e= dt * 0.5 * ((10^fb_bh[i-1]) + (10^fb_bh[i]))

		if avg_wind_e gt 0.0 then tot_wind_e= tot_wind_e + avg_wind_e
		if avg_std_e gt 0.0 then tot_std_e= tot_std_e + avg_std_e
		if avg_bh_e gt 0.0 then tot_bh_e= tot_bh_e + avg_bh_e

		;int_fb_wind[i]= unitsof40 + alog10(gadunits_in_sec) + alog10(tot_wind_e)
		;int_fb_std[i]= unitsof40 + alog10(gadunits_in_sec) + alog10(tot_std_e)
		;int_fb_bh[i]= unitsof40 + alog10(gadunits_in_sec) + alog10(tot_bh_e)

		tmp_int_fb_wind= alog10(gadunits_in_sec) + alog10(tot_wind_e)
		tmp_int_fb_std= alog10(gadunits_in_sec) + alog10(tot_std_e)
		tmp_int_fb_bh= alog10(gadunits_in_sec) + alog10(tot_bh_e)

		tmp_int_fb_wind= 10^(tmp_int_fb_wind) + 10^(zero_level_sb)
		tmp_int_fb_std= 10^(tmp_int_fb_std) + 10^(zero_level_sb)
		tmp_int_fb_bh= 10^(tmp_int_fb_bh) + 10^(zero_level_bh)

		int_fb_wind[i]= unitsof40 + alog10(tmp_int_fb_wind)
		int_fb_std[i]= unitsof40 + alog10(tmp_int_fb_std)
		int_fb_bh[i]= unitsof40 + alog10(tmp_int_fb_bh)

	endfor

        oplot, time, int_fb_wind, thick=4.0, color= 50, psym=-2, linestyle=0
        ;oplot, time, int_fb_std, thick=4.0, color= 50, psym=-3
        oplot, time, int_fb_bh, thick=12.0, color= 150, psym=-3, linestyle=1

	lstidx= n_elements(time)-1
	print, "integrated totals:"
	print, "sf= ", int_fb_std(lstidx)
	print, "bh= ", int_fb_bh(lstidx)
	print, "sf/bh= ", 10^(int_fb_std(lstidx) - int_fb_bh(lstidx))

	;xyouts, 0.68, 0.85, "Std", /normal, color= 10, charthick=3.0, charsize=1.3
	;xyouts, 0.68, 0.80, "Wind", /normal, color= 50, charthick=3.0, charsize=1.3
	;xyouts, 0.68, 0.75, "BH", /normal, color= 150, charthick=3.0, charsize=1.0
	


;---------------------------

;t_merg= fload_bh_mergertime(frun)
;t_merg= t_merg/0.7

arrow, t_merg, 58.4, t_merg, 58.9, /data, thick=6.0, color= 0, hsize=-0.4
xyouts, 2.0, 58.6, "merger", /data, charthick= 2.0, color= 0, size= 1.1

oplot, [btmin, btmin], [57.6,57.8], thick=5.0, color= 0
oplot, [btmin, btmax], [57.7,57.7], thick=5.0, color= 0
oplot, [btmax, btmax], [57.6,57.8], thick=5.0, color= 0
;xyouts, 1.42, 57.4, "'starburst", /data, charthick= 2.0, color= 0, size= 1.1
xyouts, 1.42, 57.4, "'active", /data, charthick= 2.0, color= 0, size= 1.1
xyouts, 1.7, 57.25, "phase'", /data, charthick= 2.0, color= 0, size= 1.1


;---------------------------
;---------------------------

;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"
xmax = btmax
xmin = btmin


;yaxistitle = "!6FB Energy during the Burst (10!E60!N erg)"
yaxistitle = ' '
ymax = 3.25
ymin = 0.0


;---------------------------

x0= 0.55
x1= 0.95

y0= 0.32
y1= 0.76

;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, /noerase, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        ;xtickformat='(a1)', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



idx= where(time ge xmin)
btime= time(idx)
;bsfe= int_fb_std(idx) - 59.0
bsfe= int_fb_wind(idx) - 59.0
bbhe= int_fb_bh(idx) - 59.0

bsfe= 10^(bsfe)
bbhe= 10^(bbhe)

bsfe= bsfe - bsfe(0)
bbhe= bbhe - bbhe(0)

idx=where(btime ge btmax)
print, "xxxxxxxx"
print, "starburst sf energy input= ", bsfe(idx(0))
print, "starburst bh energy input= ", bbhe(idx(0))
print, "starburst bh/sf energy input= ", bbhe(idx(0))/bsfe(idx(0))
print, "xxxxxxxx"
sblbl= 'starburst bh/sf energy input= '+string(bbhe(idx(0))/bsfe(idx(0)))
;xyouts, 0.52, 0.16, 'frun= '+frun, /normal, color= 0, charthick=3.0, charsize=1.2
;xyouts, 0.52, 0.12, sblbl, /normal, color= 0, charthick=3.0, charsize=1.2

;oplot, time, int_fb_wind, thick=4.0, color= 50, psym=-3, linestyle=0
oplot, btime, bsfe, thick=4.0, color= 50, psym=-2
oplot, btime, bbhe, thick=12.0, color= 150, psym=-3, linestyle=1

xyouts, 0.59, 0.68, "!3x!610!E59!N ergs", /normal, color= 0, charthick=3.0, charsize=1.9

xyouts, 0.74, 0.42, "Star formation", /normal, color= 50, charthick=3.0, charsize=1.8
xyouts, 0.87, 0.71, "BH", /normal, color= 150, charthick=3.0, charsize=1.8

xyouts, 0.59, 0.82, "Integrated energy input during", /normal, color= 0, charthick=3.0, charsize=1.3
xyouts, 0.64, 0.79, "the 'active phase'", /normal, color= 0, charthick=3.0, charsize=1.3


; done
; ------
device, /close




end








;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------





;-------------------------------------------------
;-------------------------------------------------
pro fbescaling, junk


if not keyword_set(junk) then begin
        print, " "
        print, " fbescaling, junk"
        print, " "
        print, " "
        return
endif

filename='fbescaling.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

;---------------------------



;xaxistitle = "!6Log Halo Mass (M!D!9n!6!N)"
;xmax = 13.25
;xmin = 10.0
xaxistitle = "!6Log Stellar Mass (M!D!9n!6!N)"
xmax = 7.0e+12
xmin = 2.0e+9
;xaxistitle = "!7r!6 (km s!E-1!N)"
;xmax = 800.0
;xmin = 30.0

yaxistitle = "!6FB Energy Ratio (BH/Starburst) "
ymax = 350.0
ymin = 0.05


;---------------------------

x0= 0.18
x1= 0.98
y0= 0.15
y1= 0.98


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
	/xlog, $
	/ylog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtickformat='exp_label', $
        ;xtickformat='(a1)', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------

;xpt1= 120.0
;xpt2= 200.0
;xpt3= 230.0
xpt1= 3.0e+9
xpt2= 5.0e+9
xpt3= 7.0e+9

; both halos in 10^10 msolar/h
total_halo_masses= [8.2, 23.8, 70.7, 190.5, 529.7, 1523.8, 5813.0]
vcs= [56, 80, 115, 160, 225, 320, 500]

xvar= vcs

; d0e2_q, d1e2_q, d2e2_q, d3e7, d4e2_q, d5e2_q, d6e2_q
bhsbratio= [4.4, 6.1, 2.0, 3.1, 28.3, 30.6, 48.9]
stmass= [0.310, 0.899, 2.66, 7.22, 18.9, 57.8, 227.4]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, bhsbratio, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 50.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40e_q', size=1.2, color=thiscolor, /data, charthick=4.0



; d0h2_q, d1h2_q, d2h2_q, d3h7, d4h2_q, d5h2_q, d6h2_q
bhsbratio= [0.6, 8.8, 2.2, 3.8, 1.7, 21.3, 80.7]
stmass= [0.314, 0.871, 2.68, 7.37, 19.5, 57.8, 226.9]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, bhsbratio, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 80.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40h_q', size=1.2, color=thiscolor, /data, charthick=4.0


; b0e, b1e, b2e, b3e, b4e, b5e, b6e
bhsbratio= [0.1, 0.6, 0.2, 1.5, 3.1, 35.7, 16.2]
stmass= [0.291, 0.876, 2.53, 6.53, 18.0, 47.7, 212.9]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, bhsbratio, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 200.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80e', size=1.2, color=thiscolor, /data, charthick=4.0


; b0h, b1h, b2h, b3h, b4h, b5h, b6h
bhsbratio= [0.1, 0.3, 1.9, 1.4, 1.5, 5.4, 11.1]
stmass= [0.302, 0.814, 2.43, 6.63, 18.7, 53.8, 219.1]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, bhsbratio, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 100.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80h', size=1.2, color=thiscolor, /data, charthick=4.0


; e0e, e1e, e2e, e3e, e4e, e5e, e6e
bhsbratio= [0.8, 2.5, 3.9, 7.2, 8.4, 24.0, 96.8]
stmass= [0.321, 0.919, 2.69, 7.38, 19.1, 59.4, 228.4]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 6
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, bhsbratio, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 30.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '20e', size=1.2, color=thiscolor, /data, charthick=4.0


; z0e, z1e, z2e, z3e, z4e, z5e, z6e
bhsbratio= [0.5, 1.3, 6.9, 9.0, 40.0, 51.8, 220.0]
stmass= [0.328, 0.963, 2.85, 7.64, 21.3, 61.3, 234.0]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 7
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, bhsbratio, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 20.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '5e', size=1.2, color=thiscolor, /data, charthick=4.0



;---------------------------



; BLACK HOLE Feedback
; ---------------------
;sample_sigmas= [xmin, vcs, xmax]   ; x is sigma
;bhmasses= (1.2e+8)*((sample_sigmas/200.0)^(3.75))

sample_stmasses= [xmin, xmax]  ; x is stellar mass
bhmasses= sample_stmasses * 0.002   ; magorian

cc= 2.9979d+10               ; speed of light in cm/sec
solarmass_in_g = 1.989d+33
BHAccEfficiency= 0.1
FBCouplingEfficiency= 0.05            ; std value
bhfb_energy= FBCouplingEfficiency * BHAccEfficiency * cc * cc * bhmasses * solarmass_in_g

; assume 90% comes from the starburst phase
bhfb_energy= 0.9 * bhfb_energy



; Starburst Feedback
; ----------------------
; x axis is sigma
; ------
;G= 43007.1
;xmin_mass= xmin * xmin * xmin / G
;xmax_mass= xmax * xmax * xmax / G
;sample_masses= [xmin_mass, total_halo_masses, xmax_mass]
;sample_stmasses= sample_masses * 0.05
;sb_mass= sample_masses * 0.001

; x axis is mass
; ------
sample_masses= [xmin, xmax]
sample_stmasses= sample_masses
sb_mass= sample_stmasses * 0.001

FBEnergy_per_SN= 1.3955d+49
sbfb_energy= FBEnergy_per_SN * sb_mass



; ----------------------

bhsbratio= bhfb_energy / sbfb_energy

;bhsbratio= (sample_sigmas/100.0)^(2.5)   ; this is the approximate scaling
; 100 km/s = 23.3
;oplot, sample_sigmas, bhsbratio, psym=-3, color=0, thick=10.0, linestyle= 2

bhsbratio= (sample_stmasses/2.0d+10)^(0.9)   ; this is the approximate scaling
oplot, sample_stmasses, bhsbratio, psym=-3, color=0, thick=10.0, linestyle= 2


; done
; ------
device, /close




end






;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------





pro fbe_m, junk


if not keyword_set(junk) then begin
        print, " "
        print, " fbe_m, junk"
        print, " "
        print, " "
        return
endif

filename='fbe_m.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

;---------------------------



;xaxistitle = "!6Log Halo Mass (M!D!9n!6!N)"
;xmax = 13.25
;xmin = 10.0
xaxistitle = "!6Log Stellar Mass (M!D!9n!6!N)"
xmax = 7.0e+12
xmin = 2.0e+9
;xaxistitle = "!7r!6 (km s!E-1!N)"
;xmax = 800.0
;xmin = 30.0

;yaxistitle = "!6FB Energy BH (10!E59!N ergs)"
;yaxistitle = "!6FB Energy Starburst (10!E59!N ergs)"
;ymax = 350.0
;ymin = 0.005
;yaxistitle = "!6SB FB Energy (10!E59!N ergs) !3x !8f!6!Dg!N!E-1!N"
;ymax = 350.0
;ymin = 0.05
yaxistitle = "!6BH FB Energy (10!E59!N ergs) !3x !8f!6!Dg!N!E-1!N"
ymax = 1250.0
ymin = 0.005


;---------------------------

x0= 0.18
x1= 0.98
y0= 0.15
y1= 0.98


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
	/xlog, $
	/ylog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtickformat='exp_label', $
        ;xtickformat='(a1)', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------

;xpt1= 120.0
;xpt2= 200.0
;xpt3= 230.0
xpt1= 3.0e+9
xpt2= 5.0e+9
xpt3= 7.0e+9

; both halos in 10^10 msolar/h
total_halo_masses= [8.2, 23.8, 70.7, 190.5, 529.7, 1523.8, 5813.0]
vcs= [56, 80, 115, 160, 225, 320, 500]

xvar= vcs

; d0e2_q, d1e2_q, d2e2_q, d3e7, d4e2_q, d5e2_q, d6e2_q
;fbe_sf= [0.0742, 0.0681, 2.1789, 0.92098, 0.03654, 0.99134, 0.4665] & yvar= fbe_sf
fbe_bh= [0.3225, 0.414, 4.314, 2.855, 1.032, 30.3088, 22.8342] & yvar= fbe_bh
yvar = yvar / 0.4
;bhsbratio= [4.4, 6.1, 2.0, 3.1, 28.3, 30.6, 48.9]
stmass= [0.310, 0.899, 2.66, 7.22, 18.9, 57.8, 227.4]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 50.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40e_q', size=1.2, color=thiscolor, /data, charthick=4.0



; d0h2_q, d1h2_q, d2h2_q, d3h7, d4h2_q, d5h2_q, d6h2_q
;fbe_sf= [0.100685, 0.217414, 0.78156, 0.80296, 4.30762, 1.19811, 0.628265] & yvar= fbe_sf
fbe_bh= [0.0644018, 1.91121, 1.71065, 3.08271, 7.45865, 25.4779, 50.7168] & yvar= fbe_bh
yvar = yvar / 0.4
;bhsbratio= [0.6, 8.8, 2.2, 3.8, 1.7, 21.3, 80.7]
stmass= [0.314, 0.871, 2.68, 7.37, 19.5, 57.8, 226.9]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 80.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40h_q', size=1.2, color=thiscolor, /data, charthick=4.0


; b0e, b1e, b2e, b3e, b4e, b5e, b6e
;fbe_sf= [0.186388, 0.680715, 8.60526, 2.60144, 4.73157, 2.76085, 5.91937] & yvar= fbe_sf
fbe_bh= [0.0186825, 0.429821, 1.81813, 3.99579, 14.7627, 98.5307, 95.8707] & yvar= fbe_bh
yvar = yvar / 0.8
;bhsbratio= [0.1, 0.6, 0.2, 1.5, 3.1, 35.7, 16.2]
stmass= [0.291, 0.876, 2.53, 6.53, 18.0, 47.7, 212.9]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 200.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80e', size=1.2, color=thiscolor, /data, charthick=4.0


; b0h, b1h, b2h, b3h, b4h, b5h, b6h
;fbe_sf= [0.177083, 0.619332, 0.93397, 2.1366, 3.55698, 8.75327, 3.79956] & yvar= fbe_sf
fbe_bh= [0.00962603, 0.164973, 1.79306, 2.91232, 5.16586, 47.5779, 42.1003] & yvar= fbe_bh
yvar = yvar / 0.8
;bhsbratio= [0.1, 0.3, 1.9, 1.4, 1.5, 5.4, 11.1]
stmass= [0.302, 0.814, 2.43, 6.63, 18.7, 53.8, 219.1]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 100.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80h', size=1.2, color=thiscolor, /data, charthick=4.0


; e0e, e1e, e2e, e3e, e4e, e5e, e6e
;fbe_sf= [0.0529996, 0.0914196, 0.104064, 0.208446, 2.90191, 0.986931, 1.76358] & yvar= fbe_sf
fbe_bh= [0.0446973, 0.230271, 0.401799, 1.50744, 24.3720, 23.7173, 170.754] & yvar= fbe_bh
yvar = yvar / 0.2
;bhsbratio= [0.8, 2.5, 3.9, 7.2, 8.4, 24.0, 96.8]
stmass= [0.321, 0.919, 2.69, 7.38, 19.1, 59.4, 228.4]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 6
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 30.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '20e', size=1.2, color=thiscolor, /data, charthick=4.0


; z0e, z1e, z2e, z3e, z4e, z5e, z6e
;fbe_sf= [0.00892448, 0.023310, 0.0404499, 0.0480618, 0.0562077, 0.0772667, 0.209244] & yvar= fbe_sf
fbe_bh= [0.00449062, 0.0302372, 0.280130, 0.432570, 2.24829, 4.00141, 46.0352] & yvar= fbe_bh
yvar = yvar / 0.05
;bhsbratio= [0.5, 1.3, 6.9, 9.0, 40.0, 51.8, 220.0]
stmass= [0.328, 0.963, 2.85, 7.64, 21.3, 61.3, 234.0]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 7
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 20.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '5e', size=1.2, color=thiscolor, /data, charthick=4.0



;---------------------------



; BLACK HOLE Feedback
; ---------------------
sample_stmasses= [xmin, xmax]  ; x is stellar mass
bhmasses= sample_stmasses * 0.002   ; magorian

cc= 2.9979d+10               ; speed of light in cm/sec
solarmass_in_g = 1.989d+33
BHAccEfficiency= 0.1
FBCouplingEfficiency= 0.05            ; std value
bhfb_energy= FBCouplingEfficiency * BHAccEfficiency * cc * cc * bhmasses * solarmass_in_g

; assume 90% comes from the starburst phase
;bhfb_energy= 0.9 * bhfb_energy  / 1.0d+59
bhfb_energy= 0.5 * bhfb_energy  / 1.0d+59  * (sample_stmasses / 7.0d+10)^(0.1)

oplot, sample_stmasses, bhfb_energy, psym=-3, color=0, thick=10.0, linestyle= 2

xyouts, 0.6, 0.35, '!9?!6 M!D*!N!E1.1!N', /normal, color= 0, size=2.2, charthick=2.0
;xyouts, 0.7, 0.35, '!9?!6 M!D*!N', /normal, color= 0, size=2.2, charthick=2.0


; Starburst Feedback
; ----------------------
;sample_stmasses= [xmin, xmax]
;sb_mass= sample_stmasses * 0.1 * (sample_stmasses / 1.0d+11)^(-0.5)   ; this makes the actual power lower than 1

;FBEnergy_per_SN= 1.3955d+49
;sbfb_energy= FBEnergy_per_SN * sb_mass

;sbfb_energy= sbfb_energy / 1.0d+59

;oplot, sample_stmasses, sbfb_energy, psym=-3, color=0, thick=10.0, linestyle= 2

;xyouts, 0.6, 0.85, '!9?!6 M!D*!N!E1/2!N', /normal, color= 0, size=2.2, charthick=2.0
;xyouts, 0.7, 0.35, '!9?!6 M!D*!N', /normal, color= 0, size=2.2, charthick=2.0

; ----------------------


; done
; ------
device, /close




end







;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------





pro fbe_3, junk


if not keyword_set(junk) then begin
        print, " "
        print, " fbe_m, junk"
        print, " "
        print, " "
        return
endif

filename='fbe_3.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize= 12.0, newysize= 24.0

;---------------------------




;   |-----------|
;   |           |
;   |  BH       |
;   | FB Energy |
;   |           |
;   |           |
;   |-----------|
;   |           |
;   |  SF       |
;   | FB Energy |
;   |           |
;   |           |
;   |-----------|
;   |           |
;   |           |
;   |  ratio    |
;   |           |
;   |           |
;   |-----------|
;     M_stellar
;
;

x0= 0.19
x1= 0.98

y0= 0.09
y1= 0.39
y2= 0.69
y3= 0.99



;---------------------------



;xaxistitle = "!6Log Stellar Mass (M!D!9n!6!N)"
;xaxistitle = "!6Stellar Mass (M!D!9n!6!N)"
;xmax = 7.0e+12
;xmin = 2.0e+9

xaxistitle = "!6Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10


;---------------------------

; note, I went back and fixed the following by hand:
;    d4e2_q

;fruns=["d0e2_q", "d1e2_q", "d2e2_q", "d3e7", "d4e2_q", "d5e2_q", "d6e2_q"]:w
;stmass_40e= get_stellarmass(fruns, subdir='ds') * 1.0d+10 / 0.7
stmass_40e= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
;stmass_40e= [0.310, 0.899, 2.66, 7.22, 18.9, 57.8, 227.4] * 1.0d+10 / 0.7      ; stellar mass
fbe_sf_40e= [0.0742, 0.0681, 2.1789, 0.92098, 3.57952, 0.99134, 0.4665] / 0.4
fbe_bh_40e= [0.3225, 0.414, 4.314, 2.855, 34.5150, 30.3088, 22.8342] / 0.4


;fruns= ['d0h2_q, d1h2_q, d2h2_q, d3h7, d4h2_q, d5h2_q, d6h2_q
stmass_40h= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
;stmass_40h= [0.314, 0.871, 2.68, 7.37, 19.5, 57.8, 226.9] * 1.0d+10 / 0.7      ; stellar mass
fbe_sf_40h= [0.100685, 0.217414, 0.78156, 0.80296, 4.30762, 1.19811, 0.628265] / 0.4
fbe_bh_40h= [0.0644018, 1.91121, 1.71065, 3.08271, 7.45865, 25.4779, 50.7168] / 0.4


; b0e, b1e, b2e, b3e, b4e, b5e, b6e
stmass_80e= [8.1668, 23.81, 70.727, 190.48, 529.71, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
;stmass_80e= [0.291, 0.876, 2.53, 6.53, 18.0, 47.7, 212.9] * 1.0d+10 / 0.7      ; stellar mass
fbe_sf_80e= [0.186388, 0.680715, 8.60526, 2.60144, 4.73157, 2.76085, 5.91937] / 0.80
fbe_bh_80e= [0.0186825, 0.429821, 1.81813, 3.99579, 14.7627, 98.5307, 95.8707] / 0.80


; b0h, b1h, b2h, b3h, b4h, b5h, b6h
stmass_80h= [8.1668, 23.81, 70.727, 190.48, 529.71, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
;stmass_80h= [0.302, 0.814, 2.43, 6.63, 18.7, 53.8, 219.1] * 1.0d+10 / 0.7      ; stellar mass
fbe_sf_80h= [0.177083, 0.619332, 0.93397, 2.1366, 3.55698, 8.75327, 3.79956] / 0.80
fbe_bh_80h= [0.00962603, 0.164973, 1.79306, 2.91232, 5.16586, 47.5779, 42.1003] / 0.80


; e0e, e1e, e2e, e3e, e4e, e5e, e6e
stmass_20e= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
;stmass_20e= [0.321, 0.919, 2.69, 7.38, 19.1, 59.4, 228.4] * 1.0d+10 / 0.7      ; stellar mass
fbe_sf_20e= [0.0529996, 0.0914196, 0.104064, 0.208446, 2.90191, 0.986931, 1.76358] / 0.20
fbe_bh_20e= [0.0446973, 0.230271, 0.401799, 1.50744, 24.3720, 23.7173, 170.754] / 0.20


; z0e, z1e, z2e, z3e, z4e, z5e, z6e
stmass_5e= [8.1668, 23.81, 70.727, 190.48, 529.709, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass
;stmass_5e= [0.328, 0.963, 2.85, 7.64, 21.3, 61.3, 234.0] * 1.0d+10 / 0.7      ; stellar mass
fbe_sf_5e= [0.00892448, 0.023310, 0.0404499, 0.0480618, 0.0562077, 0.0772667, 0.209244] / 0.05
fbe_bh_5e= [0.00449062, 0.0302372, 0.280130, 0.432570, 2.24829, 4.00141, 46.0352] / 0.05




;---------------------------

;yaxistitle = "!6BH FB Energy (10!E59!N ergs) !3x !8f!6!Dg!N!E-1!N"
yaxistitle = "!13E!6!DBH!N (10!E59!N ergs) !3x !8f!6!Dg!N!E-1!N"
ymax = 1250.0
ymin = 0.005

!p.position= [x0, y2, x1, y3]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
	/xlog, $
	/ylog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        ;xtickformat='exp_label', $
        xtickformat='(a1)', $
        ;xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------

;xpt1= 3.0e+9
;xpt2= 5.0e+9
;xpt3= 7.0e+9
;----
;xpt1= 1.4e+11
;xpt2= 2.7e+11
;xpt3= 3.5e+11
;----
xpt1= 2.0e+12
xpt2= 3.5e+12
xpt3= 5.0e+12
;xyouts, xpt3, 1.2, '!8f!6!Dg!N, orientation', size=1.2, color=0, /data, charthick=4.0

;--------
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_40e, fbe_bh_40e, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 0.12
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.02, '0.40, tilted ', size=1.2, color=thiscolor, /data, charthick=4.0

;--------

pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_40h, fbe_bh_40h, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 0.06
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.01, '0.40, co-planar', size=1.2, color=thiscolor, /data, charthick=4.0

;--------

pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_80e, fbe_bh_80e, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 0.6
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.08, '0.80, tilted', size=1.2, color=thiscolor, /data, charthick=4.0


;--------
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_80h, fbe_bh_80h, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 0.26
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.04, '0.80, co-planar', size=1.2, color=thiscolor, /data, charthick=4.0


;--------
pointt= 6
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_20e, fbe_bh_20e, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 0.025
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.004, '0.20, tilted', size=1.2, color=thiscolor, /data, charthick=4.0


;--------
pointt= 7
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_5e, fbe_bh_5e, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 0.01
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.001, '0.05, tilted', size=1.2, color=thiscolor, /data, charthick=4.0



;---------------------------

do_second_fit= 1
;do_second_fit= 0
if do_second_fit eq 1 then begin
	x_tofit= [stmass_40e, stmass_40h, stmass_80e, stmass_80h, stmass_20e, stmass_5e]
	y_tofit= [fbe_bh_40e, fbe_bh_40h, fbe_bh_80e, fbe_bh_80h, fbe_bh_20e, fbe_bh_5e]
	x2_tofit= alog10(x_tofit)
	y2_tofit= alog10(y_tofit)
	;idx= where(x2_tofit gt 11.0)
	;idx= where(x2_tofit gt 11.3)
	;x2_tofit= x2_tofit(idx)
	;y2_tofit= y_tofit(idx)

	weight2_tofit= x2_tofit
	weight2_tofit(*)= 1.0

	guess= [1.0d+2,1.0]

	linear_result = MPFITFUN('func_linear', x2_tofit, y2_tofit, weight2_tofit, guess, $
                              BESTNORM=bestnorm, DOF=dof)

	redchi2= bestnorm/dof
	print, "Reduced Chi^2 = ", redchi2

	; plot
	x= [xmin, xmax]
	;x= [1.0d+11, xmax]
	y= linear_result[0] + alog10(x)*linear_result[1]
	y= 10^(y)
	oplot, x, y, psym= -3, color=0, thick= 8.0, linestyle= 2


	xyouts, x0+0.05, y3-0.06, '!9?!6M!D*!N!E+1.13!N', /normal, color= 0, size=2.2, charthick=4.0
endif



;---------------------------

;---------------------------



;yaxistitle = "!6SB FB Energy (10!E59!N ergs) !3x !8f!6!Dg!N!E-1!N"
yaxistitle = "!13E!6!DStarburst!N (10!E59!N ergs) !3x !8f!6!Dg!N!E-1!N"
;ymax = 350.0
;ymin = 0.05
ymax = 1250.0
ymin = 0.005

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
	/xlog, $
	/ylog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        ;xtickformat='exp_label', $
        xtickformat='(a1)', $
        ;xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata, /noerase



;---------------------------

pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_40e, fbe_sf_40e, psym=-thispsym, color=thiscolor, thick=3.0

;--------

pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_40h, fbe_sf_40h, psym=-thispsym, color=thiscolor, thick=3.0

;--------

pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_80e, fbe_sf_80e, psym=-thispsym, color=thiscolor, thick=3.0

;--------
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_80h, fbe_sf_80h, psym=-thispsym, color=thiscolor, thick=3.0

;--------
pointt= 6
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_20e, fbe_sf_20e, psym=-thispsym, color=thiscolor, thick=3.0

;--------
pointt= 7
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_5e, fbe_sf_5e, psym=-thispsym, color=thiscolor, thick=3.0


;---------------------------

do_second_fit= 1
;do_second_fit= 0
if do_second_fit eq 1 then begin
	x_tofit= [stmass_40e, stmass_40h, stmass_80e, stmass_80h, stmass_20e, stmass_5e]
	y_tofit= [fbe_sf_40e, fbe_sf_40h, fbe_sf_80e, fbe_sf_80h, fbe_sf_20e, fbe_sf_5e]
	x2_tofit= alog10(x_tofit)
	y2_tofit= alog10(y_tofit)
	;idx= where(x2_tofit gt 11.0)
	;idx= where(x2_tofit gt 11.3)
	;x2_tofit= x2_tofit(idx)
	;y2_tofit= y_tofit(idx)

	weight2_tofit= x2_tofit
	weight2_tofit(*)= 1.0

	guess= [1.0d+2,1.0]

	linear_result = MPFITFUN('func_linear', x2_tofit, y2_tofit, weight2_tofit, guess, $
                              BESTNORM=bestnorm, DOF=dof)

	redchi2= bestnorm/dof
	print, "Reduced Chi^2 = ", redchi2

	; plot
	x= [xmin, xmax]
	;x= [1.0d+11, xmax]
	y= linear_result[0] + alog10(x)*linear_result[1]
	y= 10^(y)
	oplot, x, y, psym= -3, color=0, thick= 8.0, linestyle= 2


	xyouts, x0+0.05, y2-0.06, '!9?!6M!D*!N!E+0.44!N', /normal, color= 0, size=2.2, charthick=4.0
endif



;---------------------------


;yaxistitle = "!6 Ratio (BH/Starburst) "
yaxistitle = "!13E!6!DBH!N / !13E!6!DStarburst!N"
ymax = 350.0
ymin = 0.05

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
	/xlog, $
	/ylog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtickformat='exp_label', $
        ;xtickformat='(a1)', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata, /noerase



;---------------------------

pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_40e, fbe_bh_40e/fbe_sf_40e, psym=-thispsym, color=thiscolor, thick=3.0

;--------

pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_40h, fbe_bh_40h/fbe_sf_40h, psym=-thispsym, color=thiscolor, thick=3.0

;--------

pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_80e, fbe_bh_80e/fbe_sf_80e, psym=-thispsym, color=thiscolor, thick=3.0

;--------
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_80h, fbe_bh_80h/fbe_sf_80h, psym=-thispsym, color=thiscolor, thick=3.0

;--------
pointt= 6
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_20e, fbe_bh_20e/fbe_sf_20e, psym=-thispsym, color=thiscolor, thick=3.0

;--------
pointt= 7
select_thispoint, pointt, thispsym, thiscolor
oplot, stmass_5e, fbe_bh_5e/fbe_sf_5e, psym=-thispsym, color=thiscolor, thick=3.0



;---------------------------

do_second_fit= 1
;do_second_fit= 0
if do_second_fit eq 1 then begin
	x_tofit= [stmass_40e, stmass_40h, stmass_80e, stmass_80h, stmass_20e, stmass_5e]
	allbhs= [fbe_bh_40e, fbe_bh_40h, fbe_bh_80e, fbe_bh_80h, fbe_bh_20e, fbe_bh_5e]
	allsfs= [fbe_sf_40e, fbe_sf_40h, fbe_sf_80e, fbe_sf_80h, fbe_sf_20e, fbe_sf_5e]
	y_tofit= allbhs/allsfs
	x2_tofit= alog10(x_tofit)
	y2_tofit= alog10(y_tofit)
	;idx= where(x2_tofit gt 11.0)
	;idx= where(x2_tofit gt 11.3)
	;x2_tofit= x2_tofit(idx)
	;y2_tofit= y_tofit(idx)

	weight2_tofit= x2_tofit
	weight2_tofit(*)= 1.0

	guess= [1.0d+2,1.0]

	linear_result = MPFITFUN('func_linear', x2_tofit, y2_tofit, weight2_tofit, guess, $
                              BESTNORM=bestnorm, DOF=dof)

	redchi2= bestnorm/dof
	print, "Reduced Chi^2 = ", redchi2

	; plot
	x= [xmin, xmax]
	;x= [1.0d+11, xmax]
	y= linear_result[0] + alog10(x)*linear_result[1]
	y= 10^(y)
	oplot, x, y, psym= -3, color=0, thick= 8.0, linestyle= 2


	xyouts, x0+0.05, y1-0.06, '!9?!6M!D*!N!E+0.69!N', /normal, color= 0, size=2.2, charthick=4.0
endif



;---------------------------

; ----------------------


; done
; ------
device, /close




end







;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------




pro fbe_gmass, junk


if not keyword_set(junk) then begin
        print, " "
        print, " fbe_gmass, junk"
        print, " "
        print, " "
        return
endif

filename='fbe_gmass.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

;---------------------------



xaxistitle = "!6Stellar Mass (M!D!9n!6!N)"
xmax = 7.0e+12
xmin = 2.0e+9

;yaxistitle = "!6Gas Mass prior to Burst (M!D!9n!6!N)"
;ymax = 8.0e+11
;ymin = 8.0e+7
;yaxistitle = "!6Gas Mass prior to Burst (M!D!9n!6!N) !3x !8f!6!Dg!N!E-1!N"
;ymax = 2.0e+12
;ymin = 1.0e+9
;yaxistitle = "!6Gas Mass prior to Burst / M!Dinitial!N"
;yaxistitle = "!6Cold Gas Mass pre-Burst / M!Dinitial!N"
;yaxistitle = "!6Gas Mass post-Burst / M!Dinitial!N"
;yaxistitle = "!6Hot Gas Mass pre-Burst / M!Dinitial!N"
;yaxistitle = "!6pre-Burst Hot Gas Fraction "
;ymax = 1.0
;ymax = 0.9
;ymin = 0.0
;yaxistitle = "!6Mass formed during Starburst / M!Dinitial!N"
;yaxistitle = "!6Starburst M!D*!N / M!Dinitial!N"
yaxistitle = "!6Starburst M!D*!N / Pre-Burst Gas Mass"
;ymax = 2.7
ymax = 1.0
ymin = 0.0
;yaxistitle = "!6Starburst M!D*!N"
;yaxistitle = "!6Starburst M!D*!N !3x !8f!6!Dg!N!E-1!N"
;ymax = 8.0e+11
;ymin = 7.0e+8


;---------------------------

x0= 0.18
x1= 0.98
y0= 0.15
y1= 0.98


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
	/xlog, $
	;/ylog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtickformat='exp_label', $
        ;xtickformat='(a1)', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------

;xpt1= 120.0
;xpt2= 200.0
;xpt3= 230.0
xpt1= 3.0e+9
xpt2= 5.0e+9
xpt3= 7.0e+9

; both halos in 10^10 msolar/h
total_halo_masses= [8.2, 23.8, 70.7, 190.5, 529.7, 1523.8, 5813.0]
vcs= [56, 80, 115, 160, 225, 320, 500]

xvar= vcs

; d0e2_q, d1e2_q, d2e2_q, d3e7, d4e2_q, d5e2_q, d6e2_q
;    fixed d4e2_q mass
fruns= ['d0e2_q', 'd1e2_q', 'd2e2_q', 'd3e7', 'd4e2_q', 'd5e2_q', 'd6e2_q']
initgmass= [0.133936, 0.390484, 1.15992, 3.12387, 8.12085, 24.9910, 95.3331]
gmass_hot_pre= get_hotgas_preburst(fruns, subdir='ds')
;gmass_hot_pre= get_hotgas_preburst(fruns, subdir='ds')
;gmass_presb= [0.08285, 0.196422, 0.601692, 1.59755, 3.62941, 5.45510, 11.5818] & yvar= gmass_presb * 1.0d+10 / 0.7
;gmass_postsb= [0.024476, 0.0767620, 0.242225, 0.605352, 1.41633, 4.61152, 10.8050] & yvar= gmass_postsb * 1.0d+10 / 0.7
gmass_presb= get_gasmass(fruns, subdir='ds', /pre)
gmass_postsb= get_gasmass(fruns, subdir='ds', /post)
;yvar= yvar / 0.4
;yvar= gmass_presb/initgmass
;yvar= gmass_postsb/initgmass
;yvar= gmass_hot_pre / initgmass
;yvar= gmass_hot_pre / gmass_presb
yvar= (gmass_presb - gmass_hot_pre) / gmass_presb      ; M_pre_cold / M_pre
;yvar= (gmass_presb - gmass_hot_pre) / initgmass        ; M_pre_cold / M_0
;yvar= (gmass_presb - gmass_postsb) / initgmass         ; M_sb / M_0
;yvar= (gmass_presb - gmass_postsb) * 1.0d+10 / 0.7
;yvar= yvar / 0.4
stmass= [0.310, 0.899, 2.66, 7.22, 18.9, 57.8, 227.4]  & xvar= stmass * 1.0d+10 / 0.7
x_tofit= xvar
y_tofit= yvar
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 1.4e+11
;ylbl= 0.2
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40e_q', size=1.2, color=thiscolor, /data, charthick=4.0



; d0h2_q, d1h2_q, d2h2_q, d3h7, d4h2_q, d5h2_q, d6h2_q
fruns= ['d0h2_q', 'd1h2_q', 'd2h2_q', 'd3h7', 'd4h2_q', 'd5h2_q', 'd6h2_q']
initgmass= [0.133936, 0.390484, 1.15992, 3.12387, 8.12085, 24.9910, 95.3331]
gmass_hot_pre= get_hotgas_preburst(fruns, subdir='ds')
;gmass_presb= [0.0817636, 0.188613, 0.629228, 1.66048, 2.86356, 5.47834, 11.9524] & yvar= gmass_presb * 1.0d+10 / 0.7
;gmass_postsb= [0.0231762, 0.109651, 0.222358, 0.478460, 0.848990, 4.63962, 11.2564] & yvar= gmass_postsb * 1.0d+10 / 0.7
gmass_presb= get_gasmass(fruns, subdir='ds', /pre)
gmass_postsb= get_gasmass(fruns, subdir='ds', /post)
;yvar= yvar / 0.4
;yvar= gmass_presb/initgmass
;yvar= gmass_postsb/initgmass
;yvar= gmass_hot_pre / initgmass
;yvar= gmass_hot_pre / gmass_presb
yvar= (gmass_presb - gmass_hot_pre) / gmass_presb
;yvar= (gmass_presb - gmass_hot_pre) / initgmass
;yvar= (gmass_presb - gmass_hot_pre - gmass_postsb)/initgmass
;yvar= (gmass_presb - gmass_postsb)/initgmass
;yvar= (gmass_presb-gmass_postsb) * 1.0d+10 / 0.7
;yvar= yvar / 0.4
stmass= [0.314, 0.871, 2.68, 7.37, 19.5, 57.8, 226.9]  & xvar= stmass * 1.0d+10 / 0.7
x_tofit= [xvar, x_tofit]
y_tofit= [yvar, y_tofit]
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 1.0e+11
;ylbl= 0.25
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40h_q', size=1.2, color=thiscolor, /data, charthick=4.0


; b0e, b1e, b2e, b3e, b4e, b5e, b6e
fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
initgmass= [0.267872, 0.780969, 2.31984, 6.24775, 17.3745, 49.9820, 190.666]
gmass_hot_pre= get_hotgas_preburst(fruns, subdir='bs')
;gmass_presb= [0.184691, 0.521868, 1.42230, 3.60632, 9.56094, 17.3935, 30.8439] & yvar= gmass_presb * 1.0d+10 / 0.7
;gmass_postsb= [0.0658505, 0.102141, 0.395794, 1.34846, 3.83704, 14.9035, 25.3432] & yvar= gmass_postsb * 1.0d+10 / 0.7
gmass_presb= get_gasmass(fruns, subdir='bs', /pre)
gmass_postsb= get_gasmass(fruns, subdir='bs', /post)
;yvar= yvar / 0.8
;yvar= gmass_presb/initgmass
;yvar= gmass_postsb/initgmass
;yvar= gmass_hot_pre / initgmass
;yvar= gmass_hot_pre / gmass_presb
yvar= (gmass_presb - gmass_hot_pre) / gmass_presb
;yvar= (gmass_presb - gmass_hot_pre) / initgmass
;yvar= (gmass_presb-gmass_postsb)/initgmass
;yvar= (gmass_presb-gmass_postsb) * 1.0d+10 / 0.7
;yvar= yvar / 0.8
stmass= [0.291, 0.876, 2.53, 6.53, 18.0, 47.7, 212.9]  & xvar= stmass * 1.0d+10 / 0.7
x_tofit= [xvar, x_tofit]
y_tofit= [yvar, y_tofit]
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 2.5e+11
;ylbl= 0.35
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80e', size=1.2, color=thiscolor, /data, charthick=4.0


; b0h, b1h, b2h, b3h, b4h, b5h, b6h
fruns=['b0h', 'b1h', 'b2h', 'b3h', 'b4h', 'b5h', 'b6h']
initgmass= [0.267872, 0.780969, 2.31984, 6.24775, 17.3745, 49.9820, 190.666]
gmass_hot_pre= get_hotgas_preburst(fruns, subdir='bs')
;gmass_presb= [0.198527, 0.525663, 1.47413, 3.53832, 8.78339, 15.3393, 23.3505] & yvar= gmass_presb * 1.0d+10 / 0.7
;gmass_postsb= [0.0934337, 0.187313, 0.575929, 1.40659, 3.57107, 9.32880, 19.2150] & yvar= gmass_postsb * 1.0d+10 / 0.7
gmass_presb= get_gasmass(fruns, subdir='bs', /pre)
gmass_postsb= get_gasmass(fruns, subdir='bs', /post)
;yvar= yvar / 0.8
;yvar= gmass_presb/initgmass
;yvar= gmass_postsb/initgmass
;yvar= gmass_hot_pre / initgmass
;yvar= gmass_hot_pre / gmass_presb
yvar= (gmass_presb - gmass_hot_pre) / gmass_presb
;yvar= (gmass_presb - gmass_hot_pre) / initgmass
;yvar= (gmass_presb-gmass_postsb)/initgmass
;yvar= (gmass_presb-gmass_postsb) * 1.0d+10 / 0.7
;yvar= yvar / 0.8
stmass= [0.302, 0.814, 2.43, 6.63, 18.7, 53.8, 219.1]  & xvar= stmass * 1.0d+10 / 0.7
x_tofit= [xvar, x_tofit]
y_tofit= [yvar, y_tofit]
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 5.0e+11
;ylbl= 0.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80h', size=1.2, color=thiscolor, /data, charthick=4.0


; e0e, e1e, e2e, e3e, e4e, e5e, e6e
fruns= ['e0e', 'e1e', 'e2e', 'e3e', 'e4e', 'e5e', 'e6e']
initgmass= [0.0669681, 0.195242, 0.579959, 1.56194, 4.06043, 12.4955, 47.6665]
gmass_hot_pre= get_hotgas_preburst(fruns, subdir='es')
;gmass_presb= [0.0440105, 0.0911420, 0.293024, 0.713833, 1.77501, 4.90624, 11.5713] & yvar= gmass_presb * 1.0d+10 / 0.7
;gmass_postsb= [0.0145264, 0.05666, 0.204390, 0.425940, 1.17387, 3.03745, 9.75377] & yvar= gmass_postsb * 1.0d+10 / 0.7
gmass_presb= get_gasmass(fruns, subdir='es', /pre)
gmass_postsb= get_gasmass(fruns, subdir='es', /post)
;yvar= yvar / 0.2
;yvar= gmass_presb/initgmass
;yvar= gmass_postsb/initgmass
;yvar= gmass_hot_pre / initgmass
;yvar= gmass_hot_pre / gmass_presb
yvar= (gmass_presb - gmass_hot_pre) / gmass_presb
;yvar= (gmass_presb - gmass_hot_pre) / initgmass
;yvar= (gmass_presb-gmass_postsb)/initgmass
;yvar= (gmass_presb-gmass_postsb) * 1.0d+10 / 0.7
;yvar= yvar / 0.2
stmass= [0.321, 0.919, 2.69, 7.38, 19.1, 59.4, 228.4]  & xvar= stmass * 1.0d+10 / 0.7
x_tofit= [xvar, x_tofit]
y_tofit= [yvar, y_tofit]
pointt= 6
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 6.0e+10
;ylbl= 0.15
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '20e', size=1.2, color=thiscolor, /data, charthick=4.0


; z0e, z1e, z2e, z3e, z4e, z5e, z6e
fruns= ['z0e', 'z1e', 'z2e', 'z3e', 'z4e', 'z5e', 'z6e']
initgmass= [0.0167420, 0.0488105, 0.144990, 0.390484, 1.08590, 3.12387, 11.9166]
gmass_hot_pre= get_hotgas_preburst(fruns, subdir='zs')
;gmass_presb= [0.0119414, 0.0259925, 0.0678239, 0.185226, 0.445891, 1.28673, 4.47340] & yvar= gmass_presb * 1.0d+10 / 0.7
;gmass_postsb= [0.00741065, 0.0132928, 0.0476191, 0.167281, 0.418283, 1.23108, 4.30079] & yvar= gmass_postsb * 1.0d+10 / 0.7
gmass_presb= get_gasmass(fruns, subdir='zs', /pre)
gmass_postsb= get_gasmass(fruns, subdir='zs', /post)
;yvar= yvar / 0.05
;yvar= gmass_presb/initgmass
;yvar= gmass_postsb/initgmass
;yvar= gmass_hot_pre / initgmass
;yvar= gmass_hot_pre / gmass_presb
yvar= (gmass_presb - gmass_hot_pre) / gmass_presb
;yvar= (gmass_presb - gmass_hot_pre) / initgmass
;yvar= (gmass_presb-gmass_postsb)/initgmass
;yvar= (gmass_presb-gmass_postsb) * 1.0d+10 / 0.7
;yvar= yvar / 0.05
stmass= [0.328, 0.963, 2.85, 7.64, 21.3, 61.3, 234.0]  & xvar= stmass * 1.0d+10 / 0.7
x_tofit= [xvar, x_tofit]
y_tofit= [yvar, y_tofit]
pointt= 7
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 3.0e+10
;ylbl= 0.1
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '5e', size=1.2, color=thiscolor, /data, charthick=4.0



;---------------------------

do_a_fit= 1
;do_a_fit= 0
if do_a_fit eq 1 then begin
	; normalize mass to 1.0d+10
	;x_tofit= x_tofit * 1.0d-10
	x1_tofit= alog10(x_tofit)
	;y1_tofit= alog10(y_tofit)
	y1_tofit= y_tofit
	;idx= where(x1_tofit lt 11.0)
	idx= where(x1_tofit lt 11.3)
	x1_tofit= x1_tofit(idx)
	y1_tofit= y1_tofit(idx)

	; now we fit to
	;
	;   gm = N * (m_*)^m     ; nope, don't normalize anymore (m_*/10^10)^m
	;
	;   log(gm) = log(N) + m*log(m_*)
	;
	;    y = a + m*x

	weight1_tofit= x1_tofit
	weight1_tofit(*)= 1.0
	; initial guess
	;
	;   [0] - intercept
	;   [1] - slope
	;
	guess= [1.0d+2,1.0]

	; markwardt mpfit procedure
	linear_result = MPFITFUN('func_linear', x1_tofit, y1_tofit, weight1_tofit, guess, $
                              BESTNORM=bestnorm, DOF=dof)

	redchi2= bestnorm/dof
	print, "Reduced Chi^2 = ", redchi2

	; plot
	;x= [xmin, xmax]
	x= [xmin, 1.2d+11]
	y= linear_result[0] + alog10(x)*linear_result[1]
	;y= 10^(y)
	oplot, x, y, psym= -3, color=0, thick= 8.0, linestyle= 2


	;xyouts, 0.7, 0.85, '!9?!6 M!D*!N!E-0.13!N', /normal, color= 0, size=2.2, charthick=2.0
	xyouts, 0.28, 0.53, '!9?!6 M!D*!N!E-0.02!N', /normal, color= 0, size=2.2, charthick=2.0
endif


;---------------------------

do_second_fit= 1
;do_second_fit= 0
if do_second_fit eq 1 then begin
	x2_tofit= alog10(x_tofit)
	;y2_tofit= alog10(y_tofit)
	;idx= where(x2_tofit gt 11.0)
	idx= where(x2_tofit gt 11.3)
	x2_tofit= x2_tofit(idx)
	y2_tofit= y_tofit(idx)

	weight2_tofit= x2_tofit
	weight2_tofit(*)= 1.0

	guess= [1.0d+2,1.0]

	linear_result = MPFITFUN('func_linear', x2_tofit, y2_tofit, weight2_tofit, guess, $
                              BESTNORM=bestnorm, DOF=dof)

	redchi2= bestnorm/dof
	print, "Reduced Chi^2 = ", redchi2

	; plot
	;x= [xmin, xmax]
	x= [1.2d+11, xmax]
	y= linear_result[0] + alog10(x)*linear_result[1]
	oplot, x, y, psym= -3, color=0, thick= 8.0, linestyle= 2


	xyouts, 0.70, 0.32, '!9?!6 M!D*!N!E-0.26!N', /normal, color= 0, size=2.2, charthick=2.0
endif



;---------------------------

;sample_stmasses= [xmin,xmax]
;bhsbratio= (sample_stmasses/1.0d+9)^(-0.5)   ; this is the approximate scaling
;bhsbratio= 6.0d+9 * (sample_stmasses/1.0d+10)^(1.0)   ; this is the approximate scaling
;bhsbratio= 7.0d+9 * (sample_stmasses/1.0d+10)^(0.9)   ; this is the approximate scaling
;bhsbratio= sample_stmasses
;oplot, sample_stmasses, bhsbratio, psym=-3, color=0, thick=10.0, linestyle= 2

;xyouts, 0.7, 0.85, '!9?!6 M!D*!N!E-0.2!N', /normal, color= 0, size=2.2, charthick=2.0
;xyouts, 0.7, 0.35, '!9?!6 M!D*!N!E0.9!N', /normal, color= 0, size=2.2, charthick=2.0
;xyouts, 0.7, 0.35, '!9?!6 M!D*!N', /normal, color= 0, size=2.2, charthick=2.0


; done
; ------
device, /close




end








;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------




pro fbe_bhmass, junk


if not keyword_set(junk) then begin
        print, " "
        print, " fbe_bhmass, junk"
        print, " "
        print, " "
        return
endif

filename='fbe_bhmass.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

;---------------------------



xaxistitle = "!6Stellar Mass (M!D!9n!6!N)"
xmax = 7.0e+12
xmin = 2.0e+9

yaxistitle = "!6Fraction of BH Mass Acc. during Active"
ymax = 0.9
;ymin = 0.01
ymin = 0.0


;---------------------------

x0= 0.18
x1= 0.98
y0= 0.15
y1= 0.98


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
	/xlog, $
	;/ylog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtickformat='exp_label', $
        ;xtickformat='(a1)', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------

;xpt1= 120.0
;xpt2= 200.0
;xpt3= 230.0
xpt1= 3.0e+9
xpt2= 5.0e+9
xpt3= 7.0e+9

; both halos in 10^10 msolar/h
total_halo_masses= [8.2, 23.8, 70.7, 190.5, 529.7, 1523.8, 5813.0]
vcs= [56, 80, 115, 160, 225, 320, 500]

xvar= vcs

; d0e2_q, d1e2_q, d2e2_q, d3e7, d4e2_q, d5e2_q, d6e2_q
bh_pre= [452660., 3.54192e+06, 9.55430e+06, 1.75119e+07, 2.25149e+08, 4.08617e+08, 2.05534e+09]
bh_post= [ 2.24288e+06, 7.34337e+06, 2.78140e+07, 5.05170e+07, 2.06005e+08, 6.28631e+08, 1.74839e+09]
bh_final= [2.25676e+06, 7.35489e+06, 2.78530e+07, 5.09988e+07, 2.06020e+08, 6.29633e+08, 1.75629e+09]
yvar= (bh_post-bh_pre)/bh_final
stmass= [0.310, 0.899, 2.66, 7.22, 18.9, 57.8, 227.4]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
;ylbl= 1.4e+11
ylbl= 0.2
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40e_q', size=1.2, color=thiscolor, /data, charthick=4.0



; d0h2_q, d1h2_q, d2h2_q, d3h7, d4h2_q, d5h2_q, d6h2_q
bh_pre= [ 484364., 4.55856e+06, 7.80890e+06, 2.12808e+07, 1.69754e+07, 4.68910e+08, 2.29284e+09]
bh_post= [ 1.37406e+06, 7.15526e+06, 2.50843e+07, 5.80834e+07, 9.90585e+07, 6.97035e+08, 2.21257e+09]
bh_final= [ 1.38734e+06, 7.18735e+06, 2.51885e+07, 5.82702e+07, 1.02955e+08, 6.97372e+08, 2.21568e+09]

yvar= (bh_post-bh_pre)/bh_final
stmass= [0.314, 0.871, 2.68, 7.37, 19.5, 57.8, 226.9]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
;ylbl= 1.0e+11
ylbl= 0.25
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40h_q', size=1.2, color=thiscolor, /data, charthick=4.0


; b0e, b1e, b2e, b3e, b4e, b5e, b6e
bh_pre= [ 254164., 398596., 2.96277e+06, 6.04874e+07, 6.45693e+07, 2.94191e+08, 1.18148e+09]
bh_post= [ 434605., 4.23624e+06, 1.53070e+07, 8.97064e+07, 2.57521e+08, 8.94738e+08, 1.86945e+09]
bh_final= [ 726210., 4.29924e+06, 1.54294e+07, 9.06945e+07, 2.62990e+08, 9.05095e+08, 1.87154e+09]
yvar= (bh_post-bh_pre)/bh_final
stmass= [0.291, 0.876, 2.53, 6.53, 18.0, 47.7, 212.9]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
;ylbl= 2.5e+11
ylbl= 0.35
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80e', size=1.2, color=thiscolor, /data, charthick=4.0


; b0h, b1h, b2h, b3h, b4h, b5h, b6h
bh_pre= [ 247384., 385756., 2.14952e+06, 7.15200e+07, 6.94630e+07, 2.86089e+08, 1.41094e+09]
bh_post= [ 336600., 1.78909e+06, 9.35678e+06, 7.53870e+07, 2.06921e+08, 9.67482e+08, 1.66638e+09]
bh_final= [ 618250., 1.91196e+06, 9.90712e+06, 7.70136e+07, 2.13563e+08, 9.70239e+08, 1.68025e+09]
yvar= (bh_post-bh_pre)/bh_final
stmass= [0.302, 0.814, 2.43, 6.63, 18.7, 53.8, 219.1]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
;ylbl= 5.0e+11
ylbl= 0.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80h', size=1.2, color=thiscolor, /data, charthick=4.0


; e0e, e1e, e2e, e3e, e4e, e5e, e6e
bh_pre= [ 367304., 2.55490e+06, 7.46752e+06, 2.00963e+07, 9.65755e+07, 3.80805e+08, 1.70910e+09]
bh_post= [ 849580., 3.58125e+06, 1.69626e+07, 6.39176e+07, 2.18664e+08, 7.82972e+08, 2.77722e+09]
bh_final= [ 924160., 3.58748e+06, 1.69787e+07, 6.39330e+07, 2.18743e+08, 7.83225e+08, 2.78397e+09]
yvar= (bh_post-bh_pre)/bh_final
stmass= [0.321, 0.919, 2.69, 7.38, 19.1, 59.4, 228.4]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 6
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
;ylbl= 6.0e+10
ylbl= 0.15
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '20e', size=1.2, color=thiscolor, /data, charthick=4.0


; z0e, z1e, z2e, z3e, z4e, z5e, z6e
bh_pre= [ 251956., 461268., 2.30592e+06, 1.13233e+07, 2.45681e+07, 1.17247e+08, 4.93031e+08]
bh_post= [ 281393., 655990., 2.74041e+06, 1.27224e+07, 3.74759e+07, 1.32837e+08, 6.87876e+08]
bh_final= [ 291950., 715450., 3.08483e+06, 1.54410e+07, 5.10591e+07, 1.78277e+08, 9.33218e+08]
yvar= (bh_post-bh_pre)/bh_final
stmass= [0.328, 0.963, 2.85, 7.64, 21.3, 61.3, 234.0]  & xvar= stmass * 1.0d+10 / 0.7
pointt= 7
select_thispoint, pointt, thispsym, thiscolor
oplot, xvar, yvar, psym=-thispsym, color=thiscolor, thick=3.0
;ylbl= 3.0e+10
ylbl= 0.1
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '5e', size=1.2, color=thiscolor, /data, charthick=4.0



;---------------------------



sample_stmasses= [xmin,xmax]
bhsbratio= (sample_stmasses/1.0d+9)^(-0.5)   ; this is the approximate scaling
;bhsbratio= 6.0d+9 * (sample_stmasses/1.0d+10)^(1.0)   ; this is the approximate scaling
;bhsbratio= 7.0d+9 * (sample_stmasses/1.0d+10)^(0.9)   ; this is the approximate scaling
;bhsbratio= sample_stmasses
oplot, sample_stmasses, bhsbratio, psym=-3, color=0, thick=10.0, linestyle= 2

xyouts, 0.7, 0.85, '!9?!6 M!D*!N!E-0.2!N', /normal, color= 0, size=2.2, charthick=2.0
;xyouts, 0.7, 0.35, '!9?!6 M!D*!N!E0.9!N', /normal, color= 0, size=2.2, charthick=2.0
;xyouts, 0.7, 0.35, '!9?!6 M!D*!N', /normal, color= 0, size=2.2, charthick=2.0


; done
; ------
device, /close




end





;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------


function get_stellarmass, fruns, subdir=subdir

nfruns= n_elements(fruns)

final_stellar_mass= fltarr(nfruns)

for i= 0, nfruns-1 do begin

	frun= '/raid4/tcox/'+subdir+'/'+fruns[i]

	; determine the number of
	; snapshots in frun directory
	;spawn, "ls "+frun+"/snapshot* | wc ",result
	;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
	spawn, "/bin/ls "+frun+"/snap* | wc ",result
	nsnaps=long(result[0])
	lastsnap= nsnaps-1

	ok= fload_snapshot_bh(frun, lastsnap)

	final_stellar_mass[i]= total(fload_allstars_mass(1))
endfor

return, final_stellar_mass

end







pro compute_stmass, junk


print, "stella mass= ", do_stmass("/raid4/tcox/ds/d0e2_q")
print, "stella mass= ", do_stmass("/raid4/tcox/ds/d1e2_q")
print, "stella mass= ", do_stmass("/raid4/tcox/ds/d2e2_q")
print, "stella mass= ", do_stmass("/raid4/tcox/ds/d3e7")
print, "stella mass= ", do_stmass("/raid4/tcox/ds/d4e2_q")
print, "stella mass= ", do_stmass("/raid4/tcox/ds/d5e2_q")
print, "stella mass= ", do_stmass("/raid4/tcox/ds/d6e2_q")

print, "stella mass= ", do_stmass("/raid4/tcox/ds/d0h2_q")
print, "stella mass= ", do_stmass("/raid4/tcox/ds/d1h2_q")
print, "stella mass= ", do_stmass("/raid4/tcox/ds/d2h2_q")
print, "stella mass= ", do_stmass("/raid4/tcox/ds/d3h7")
print, "stella mass= ", do_stmass("/raid4/tcox/ds/d4h2_q")
print, "stella mass= ", do_stmass("/raid4/tcox/ds/d5h2_q")
print, "stella mass= ", do_stmass("/raid4/tcox/ds/d6h2_q")

print, "stella mass= ", do_stmass("/raid4/tcox/bs/b0e")
print, "stella mass= ", do_stmass("/raid4/tcox/bs/b1e")
print, "stella mass= ", do_stmass("/raid4/tcox/bs/b2e")
print, "stella mass= ", do_stmass("/raid4/tcox/bs/b3e")
print, "stella mass= ", do_stmass("/raid4/tcox/bs/b4e")
print, "stella mass= ", do_stmass("/raid4/tcox/bs/b5e")
print, "stella mass= ", do_stmass("/raid4/tcox/bs/b6e")

print, "stella mass= ", do_stmass("/raid4/tcox/bs/b0h")
print, "stella mass= ", do_stmass("/raid4/tcox/bs/b1h")
print, "stella mass= ", do_stmass("/raid4/tcox/bs/b2h")
print, "stella mass= ", do_stmass("/raid4/tcox/bs/b3h")
print, "stella mass= ", do_stmass("/raid4/tcox/bs/b4h")
print, "stella mass= ", do_stmass("/raid4/tcox/bs/b5h")
print, "stella mass= ", do_stmass("/raid4/tcox/bs/b6h")

print, "stella mass= ", do_stmass("/raid4/tcox/es/e0e")
print, "stella mass= ", do_stmass("/raid4/tcox/es/e1e")
print, "stella mass= ", do_stmass("/raid4/tcox/es/e2e")
print, "stella mass= ", do_stmass("/raid4/tcox/es/e3e")
print, "stella mass= ", do_stmass("/raid4/tcox/es/e4e")
print, "stella mass= ", do_stmass("/raid4/tcox/es/e5e")
print, "stella mass= ", do_stmass("/raid4/tcox/es/e6e")

print, "stella mass= ", do_stmass("/raid4/tcox/zs/z0e")
print, "stella mass= ", do_stmass("/raid4/tcox/zs/z1e")
print, "stella mass= ", do_stmass("/raid4/tcox/zs/z2e")
print, "stella mass= ", do_stmass("/raid4/tcox/zs/z3e")
print, "stella mass= ", do_stmass("/raid4/tcox/zs/z4e")
print, "stella mass= ", do_stmass("/raid4/tcox/zs/z5e")
print, "stella mass= ", do_stmass("/raid4/tcox/zs/z6e")


end







;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------



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




;------------------------------------------------------------------





function get_hotgas_preburst, fruns, subdir=subdir

nfruns= n_elements(fruns)

total_gas_pre= fltarr(nfruns)
hot_gas_pre= fltarr(nfruns)
cold_gas_pre= fltarr(nfruns)

for i= 0, nfruns-1 do begin

	frun= '/raid4/tcox/'+subdir+'/'+fruns[i]
	t_merg= fload_bh_mergertime(frun)
	;t_merg= t_merg/0.7
	print, "t_merg= ", t_merg
	t_sb= t_merg-0.3
	print, "beginning of SB time= ", t_sb

	read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf

	idx=where(time ge t_sb)

	f= (t_sb - time(idx(0)-1)) / (time(idx(0)) - time(idx(0)-1))

	total_gas_pre[i]= (1-f)*gas_tot(idx(0)-1) + f*gas_tot(idx(0))
	hot_gas_pre[i]= (1-f)*gas_hot(idx(0)-1) + f*gas_hot(idx(0))
	cold_gas_pre[i]= (1-f)*gas_cold(idx(0)-1) + f*gas_cold(idx(0))
endfor

return, hot_gas_pre


end



;------------------------------------------------------------------



function get_hotgas_postburst, fruns, subdir=subdir

nfruns= n_elements(fruns)

total_gas_post= fltarr(nfruns)
hot_gas_post= fltarr(nfruns)
cold_gas_post= fltarr(nfruns)

for i= 0, nfruns-1 do begin

        frun= '/raid4/tcox/'+subdir+'/'+fruns[i]
        t_merg= fload_bh_mergertime(frun)
        ;t_merg= t_merg/0.7
        ;print, "t_merg= ", t_merg
        t_sb= t_merg-0.3
        t_postsb= t_merg-0.3+1.0
        ;print, "beginning of SB time= ", t_sb

        read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf

        idx=where(time ge t_postsb)

        f= (t_postsb - time(idx(0)-1)) / (time(idx(0)) - time(idx(0)-1))

        total_gas_post[i]= (1-f)*gas_tot(idx(0)-1) + f*gas_tot(idx(0))
        hot_gas_post[i]= (1-f)*gas_hot(idx(0)-1) + f*gas_hot(idx(0))
        cold_gas_post[i]= (1-f)*gas_cold(idx(0)-1) + f*gas_cold(idx(0))
endfor

return, hot_gas_post


end



;------------------------------------------------------------------






function get_gasmass, fruns, subdir=subdir, pre=pre, post=post

; default is pre-burst

nfruns= n_elements(fruns)

total_gas_post= fltarr(nfruns)
total_gas_pre= fltarr(nfruns)

for i= 0, nfruns-1 do begin

        frun= '/raid4/tcox/'+subdir+'/'+fruns[i]
        t_merg= fload_bh_mergertime(frun)
        ;t_merg= t_merg/0.7
        ;print, "t_merg= ", t_merg
        t_sb= t_merg-0.3
        t_postsb= t_merg-0.3+1.0
        ;t_postsb= t_merg-0.3+0.7
        ;print, "beginning of SB time= ", t_sb

	read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf

	; pre
	idx=where(time ge t_sb)
	f= (t_sb - time(idx(0)-1)) / (time(idx(0)) - time(idx(0)-1))
	total_gas_pre[i]= (1-f)*gas_tot(idx(0)-1) + f*gas_tot(idx(0))

	; post
	idx=where(time ge t_postsb)
	f= (t_postsb - time(idx(0)-1)) / (time(idx(0)) - time(idx(0)-1))
	total_gas_post[i]= (1-f)*gas_tot(idx(0)-1) + f*gas_tot(idx(0))

endfor

if keyword_set(post) then return, total_gas_post
return, total_gas_pre


end





;-------------------------------------------------------------------------





pro compute_gmass, junk


print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d0e2_q")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d1e2_q")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d2e2_q")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d3e7")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d4e2_q")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d5e2_q")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d6e2_q")

print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d0h2_q")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d1h2_q")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d2h2_q")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d3h7")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d4h2_q")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d5h2_q")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/ds/d6h2_q")

print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b0e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b1e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b2e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b3e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b4e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b5e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b6e")

print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b0h")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b1h")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b2h")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b3h")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b4h")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b5h")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/bs/b6h")

print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/es/e0e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/es/e1e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/es/e2e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/es/e3e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/es/e4e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/es/e5e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/es/e6e")

print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/zs/z0e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/zs/z1e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/zs/z2e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/zs/z3e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/zs/z4e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/zs/z5e")
print, "gas mass prior to burst= ", do_gmass("/raid4/tcox/zs/z6e")


end





;-----------------------------------------------------------------





function do_gmass_post, frun

t_merg= fload_bh_mergertime(frun)
;t_merg= t_merg/0.7
print, "t_merg= ", t_merg
t_sb= t_merg-0.3
t_sb_post= t_merg-0.3+1.0
print, "beginning of SB time= ", t_sb
print, "end of SB time= ", t_sb_post

read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf

idx=where(time ge t_sb_post)

f= (t_sb_post - time(idx(0)-1)) / (time(idx(0)) - time(idx(0)-1))

gasmass= (1-f)*gas_tot(idx(0)-1) + f*gas_tot(idx(0))

return, gasmass


end







pro compute_gmass, junk


print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d0e2_q")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d1e2_q")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d2e2_q")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d3e7")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d4e2_q")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d5e2_q")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d6e2_q")

print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d0h2_q")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d1h2_q")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d2h2_q")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d3h7")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d4h2_q")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d5h2_q")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/ds/d6h2_q")

print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b0e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b1e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b2e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b3e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b4e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b5e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b6e")

print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b0h")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b1h")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b2h")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b3h")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b4h")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b5h")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/bs/b6h")

print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/es/e0e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/es/e1e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/es/e2e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/es/e3e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/es/e4e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/es/e5e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/es/e6e")

print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/zs/z0e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/zs/z1e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/zs/z2e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/zs/z3e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/zs/z4e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/zs/z5e")
print, "gas mass after to burst= ", do_gmass_post("/raid4/tcox/zs/z6e")


end






;-----------------------------------------------------------------




function do_0gmass, frun

ok= fload_snapshot_bh(frun,0)
return, total(fload_gas_mass(1))

end







pro compute_0gmass, junk


print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d0e2_q")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d1e2_q")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d2e2_q")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d3e7")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d4e2_q")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d5e2_q")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d6e2_q")

;print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d0h2_q")
;print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d1h2_q")
;print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d2h2_q")
;print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d3h7")
;print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d4h2_q")
;print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d5h2_q")
;print, "initial gas mass= ", do_0gmass("/raid4/tcox/ds/d6h2_q")

print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b0e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b1e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b2e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b3e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b4e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b5e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b6e")

;print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b0h")
;print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b1h")
;print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b2h")
;print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b3h")
;print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b4h")
;print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b5h")
;print, "initial gas mass= ", do_0gmass("/raid4/tcox/bs/b6h")

print, "initial gas mass= ", do_0gmass("/raid4/tcox/es/e0e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/es/e1e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/es/e2e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/es/e3e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/es/e4e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/es/e5e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/es/e6e")

print, "initial gas mass= ", do_0gmass("/raid4/tcox/zs/z0e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/zs/z1e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/zs/z2e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/zs/z3e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/zs/z4e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/zs/z5e")
print, "initial gas mass= ", do_0gmass("/raid4/tcox/zs/z6e")


end








;-----------------------------------------------------------------





function do_bhmass_prepostSB, frun

t_merg= fload_bh_mergertime(frun)
;t_merg= t_merg/0.7
print, "t_merg= ", t_merg
t_sb= t_merg-0.3
t_sb_post= t_merg-0.3+1.0
print, "beginning of SB time= ", t_sb
print, "end of SB time= ", t_sb_post

; used blackholes.txt, if possible
;
; need to .run bh_multi
open_blackhole_txt, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, $
			bh_totalmass, bh_mdot_edd


idx=where(bhtime ge t_sb)
f= (t_sb_post - bhtime(idx(0)-1)) / (bhtime(idx(0)) - bhtime(idx(0)-1))
bhmass_preburst= (1-f)*bh_mass(idx(0)-1) + f*bh_mass(idx(0))
print, "preburst BH mass= ", bhmass_preburst

idx=where(bhtime ge t_sb_post)
f= (t_sb_post - bhtime(idx(0)-1)) / (bhtime(idx(0)) - bhtime(idx(0)-1))
bhmass_postburst= (1-f)*bh_mass(idx(0)-1) + f*bh_mass(idx(0))
print, "postburst BH mass= ", bhmass_postburst


print, "final BH mass= ", max(bh_mass)

return, bhmass_postburst-bhmass_preburst


end







pro compute_bhmass, junk


print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d0e2_q")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d1e2_q")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d2e2_q")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d3e7")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d4e2_q")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d5e2_q")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d6e2_q")

print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d0h2_q")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d1h2_q")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d2h2_q")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d3h7")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d4h2_q")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d5h2_q")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/ds/d6h2_q")

print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b0e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b1e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b2e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b3e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b4e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b5e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b6e")

print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b0h")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b1h")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b2h")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b3h")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b4h")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b5h")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/bs/b6h")

print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/es/e0e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/es/e1e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/es/e2e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/es/e3e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/es/e4e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/es/e5e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/es/e6e")

print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/zs/z0e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/zs/z1e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/zs/z2e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/zs/z3e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/zs/z4e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/zs/z5e")
print, "BH Mass formed during burst= ", do_bhmass_prepostSB("/raid4/tcox/zs/z6e")


end






;-----------------------------------------------------------------



