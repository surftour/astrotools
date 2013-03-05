
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






; ================================================================================







;-------------------------------------------------
;-------------------------------------------------
pro plot_fbenergy, junk


if not keyword_set(junk) then begin
        print, " "
        print, " plot_fbenergy, junk"
        print, " "
        print, " "
        return
endif

filename='fbenergy.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=12.0, newysize= 20.0

;---------------------------


        ;frun= "/raid4/tcox/vc3vc3e_2"
        ;frun= "/raid4/tcox/vc3vc3e_sb8"
        ;frun= "/raid4/tcox/vc3vc3e_sb1"
        ;frun= "/raid4/tcox/vc3vc3e_sb10"
        frun= "/raid4/tcox/sbw/sb10BH"
        ;frun= "/raid4/tcox/vc3vc3e_sb13"
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

x0= 0.15
x1= 0.98

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

        ;oplot, time, fb_wind, thick=4.0, color= 50, psym=-3, linestyle=0
        oplot, time, fb_std, thick=4.0, color= 50, psym=-3
        oplot, time, fb_bh, thick=12.0, color= 150, psym=-3, linestyle=1

	;xyouts, 0.68, 0.85, "Std", /normal, color= 10, charthick=3.0, charsize=1.3
	xyouts, 0.56, 0.90, "Starformation", /normal, color= 50, charthick=3.0, charsize=1.8
	xyouts, 0.66, 0.85, "BH", /normal, color= 150, charthick=3.0, charsize=1.8
	

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

	; add in a zero point for feedback energy
	; to start the simulation
	
	; starburst (vc3, 40% gas sim)
	zero_level_sb= alog10(0.5 * 1.4d+49 * 4.68581d+10)
	zero_level_sb= zero_level_sb - 40.0

	; black hole
	cc= 2.9979d+10
	zero_level_bh= alog10(0.05 * 0.1 * cc * cc * 2.0d+5 * 1.989d+33)
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

        ;oplot, time, int_fb_wind, thick=4.0, color= 50, psym=-3, linestyle=0
        oplot, time, int_fb_std, thick=4.0, color= 50, psym=-3
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




; done
; ------
device, /close




end








;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------





;-------------------------------------------------
;-------------------------------------------------
pro plot_fbenergy2, junk


if not keyword_set(junk) then begin
        print, " "
        print, " plot_fbenergy2, junk"
        print, " "
        print, " "
        return
endif

filename='fbenergy2.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=24.0, newysize= 20.0

;---------------------------


        ;frun= "/raid4/tcox/vc3vc3e_2"
        ;frun= "/raid4/tcox/vc3vc3e_sb8"
        ;frun= "/raid4/tcox/vc3vc3e_sb1"
        ;frun= "/raid4/tcox/vc3vc3e_sb10"
        frun= "/raid4/tcox/sbw/sb10BH"
        ;frun= "/raid4/tcox/vc3vc3e_sb13"
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
;
;   upper-left figure
;

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
;  put in yellow shaded region
;  and redraw axes
;

btmin= 1.4
btmax= 2.4

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
;
;    lower-left figure
;

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

btmin= 1.4
btmax= 2.4

timerange= [btmin, btmin, btmax, btmax, btmin]
erange= [ymin,ymax,ymax,ymin,ymin]
polyfill, timerange, erange, /data, color= 250, /fill


plot, [1.0],[1.0], psym=-3, /noerase, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        charthick=3.0, xtitle=xaxistitle, ytitle=yaxistitle, /nodata

;---------------------------

	; add in a zero point for feedback energy
	; to start the simulation
	
	; starburst (vc3, 40% gas sim)
	zero_level_sb= alog10(0.5 * 1.4d+49 * 4.68581d+10)
	zero_level_sb= zero_level_sb - 40.0

	; black hole
	cc= 2.9979d+10
	zero_level_bh= alog10(0.05 * 0.1 * cc * cc * 2.0d+5 * 1.989d+33)
	zero_level_bh= zero_level_bh - 40.0


;---------------------------

	; determine the integrated quantities

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
	


;---------------------------------------
;  draw arrow and denote "active phase"
;

t_merg= fload_bh_mergertime(frun)
t_merg= t_merg/0.7

arrow, t_merg, 58.4, t_merg, 58.9, /data, thick=6.0, color= 0, hsize=-0.4
xyouts, 2.0, 58.6, "merger", /data, charthick= 2.0, color= 0, size= 1.1

oplot, [btmin, btmin], [57.6,57.8], thick=5.0, color= 0
oplot, [btmin, btmax], [57.7,57.7], thick=5.0, color= 0
oplot, [btmax, btmax], [57.6,57.8], thick=5.0, color= 0
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


