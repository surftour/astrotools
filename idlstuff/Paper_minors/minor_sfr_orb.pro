;======================================================================
;
;
;
;
;
;       this script compares the sfr (and other variants)
;       between one merger and it's isolated constituents
;
;
;
;
;
;
;
;======================================================================
pro resample_sfrhist, oldsfr, oldtime, newsfr, newtime, Nelements=Nelements

	n_old= n_elements(oldsfr)

	newsfr= fltarr(Nelements)
	newtime= fltarr(Nelements)

	timemax= max(oldtime)
	timemin= min(oldtime)
	dt= (timemax - timemin) / (1.0*Nelements)

	; now resample to determine the new sfrhist
	for i=0, Nelements-1 do begin

		t= dt * (i + 0.5)

		jidx= where(oldtime ge t)
		j= jidx[0]
		sfr_0= oldsfr[j]

		if j ge 2 then sfr_m2= oldsfr[j-2] else sfr_m2= sfr_0
		if j ge 1 then sfr_m1= oldsfr[j-1] else sfr_m1= sfr_0
		if j le (n_old-2) then sfr_p1= oldsfr[j+1] else sfr_p1= sfr_0
		if j le (n_old-3) then sfr_p2= oldsfr[j+2] else sfr_p2= sfr_0

		newsfr[i]= (1.0 * sfr_m2 + $
			    2.0 * sfr_m1 + $
			    5.0 * sfr_0  + $
			    2.0 * sfr_p1 + $
			    1.0 * sfr_p2)/11.0

		newtime[i]= t
	endfor

end









;================================================================================================================
;
;
;      |----------------|----------------|----------------|----------------|----------------|----------------|
;      |                |                |                |                |                |                |
;      |                |                |                |                |                |                |
;      |      SFR       |      SFR       |      SFR       |      SFR       |      SFR       |      SFR       |
;      |                |                |                |                |                |                |
;      |                |                |                |                |                |                |
;      |----------------|----------------|----------------|----------------|----------------|----------------|
;
;
;
;================================================================================================================


pro msfr, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr, junk"
   print, "  "
   print, "  "
   return
endif


filename='minorsfr_row.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=22, newysize=10
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------


Rtitle= "Pericentric Distance"
R_peris= ['R!Dperi!N= 0.01', '6.8', '13.6', '27.2', '64.4']  ; R_peri (same for both)

msg= 'G3G2'
G3fruns= ["G3G2h-u1", "G3G2e-u1", "G3G2-u3", "G3G2f-u1", "G3G2g-u1"]
iso1="G3il-u1a"
iso2="G2im-u1a"

;msg= 'Mass ratio 1:5.8'
;G3fruns= ["G3G1h-u1", "G3G1e-u1", "G3G1-u3", "G3G1f-u1", "G3G1g-u1"]



;--------------------------------------
;--------------------------------------



y0= 0.18
y1= 0.98

x0= 0.08
x99= 0.99
xcols= 5
xs= (x99-x0)/xcols
x1= x0+xs
x2= x0+xs+xs
x3= x0+xs+xs+xs
x4= x0+xs+xs+xs+xs
x5= x0+xs+xs+xs+xs+xs




;--------------------------------------
; Left-most Plot
;--------------------------------------

xaxistitle= "Time (Gyr)"
xmax = 6.6
xmin = 0.0

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
ymax = 8.0
ymin = 0

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtitle=xaxistitle, ytitle=yaxistitle

	; 0
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[0], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100

	oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0 
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 0.3, 6.46, R_peris[0], /data, charthick=4, size=1.33, color=0


        xyouts, 0.3, 7.26, msg, /data, charthick=4, size=1.33, color=0



; Plot #2
;--------------------------------------

!p.position= [x1, y0, x2, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=xaxistitle, /noerase


        ; 1
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[1], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


	oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0 
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 0.3, 6.46, R_peris[1], /data, charthick=4, size=1.33, color=0



; Plot #3
;--------------------------------------

!p.position= [x2, y0, x3, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=xaxistitle, /noerase


        ; 1
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[2], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 0.3, 6.46, R_peris[2], /data, charthick=4, size=1.33, color=0




; Plot #4
;--------------------------------------

!p.position= [x3, y0, x4, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=xaxistitle, /noerase


        ; 1
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[3], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 0.3, 6.46, R_peris[3], /data, charthick=4, size=1.33, color=0




; Plot #5
;--------------------------------------

xmax = 9.6

!p.position= [x4, y0, x5, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=xaxistitle, /noerase


        ; 1
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[4], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 0.3, 6.46, R_peris[4], /data, charthick=4, size=1.33, color=0





;--------------------------------------
;--------------------------------------




; done
; ----------
device, /close

end










;================================================================================================================
;
;
;      |----------------|----------------|----------------|----------------|----------------|
;      |                |                |                |                |                |
;      |                |                |                |                |                |
;      |      SFR       |      SFR       |      SFR       |      SFR       |      SFR       |
;      |                |                |                |                |                |
;      |                |                |                |                |                |
;      |----------------|----------------|----------------|----------------|----------------|
;      |                |                |                |                |                |
;      |                |                |                |                |                |
;      |  GasConsump    |   GasConsump   |   GasConsump   |   GasConsump   |   GasConsump   |
;      |                |                |                |                |                |
;      |                |                |                |                |                |
;      |----------------|----------------|----------------|----------------|----------------|
;
;
;
;
;================================================================================================================


pro msfr2, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr2, junk"
   print, "  "
   print, "  "
   return
endif


filename='minorsfr_row2.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=22, newysize=10



;-------------------------------------------
;  set variables
;-------------------------------------------


Rtitle= "Pericentric Distance"
;R_peris= ['circular', 'R!Dperi!N= 6.8', '13.6', '27.2', '64.4']  ; R_peri (same for both)
R_peris= ['circular', '        6.8', '13.6', '27.2', '64.4']  ; R_peri (same for both)

;msg= 'G3G2'
;G3fruns= ["G3G2h-u1", "G3G2e-u1", "G3G2-u3", "G3G2f-u1", "G3G2g-u1"]
;iso1="G3il-u1a"
;iso2="G2im-u1a"

;msg= 'Mass ratio 1:5.8'
msg= 'G3G1'
G3fruns= ["G3G1h-u1", "G3G1e-u1", "G3G1-u3", "G3G1f-u1", "G3G1g-u1"]
iso1="G3il-u1a"
iso2="G1i-u1a"



;--------------------------------------
;--------------------------------------



y0= 0.14
y99= 0.98
ycols= 2
ys= (y99-y0)/ycols
y1= y0+ys
y2= y0+ys+ys

x0= 0.08
x99= 0.99
xcols= 5
xs= (x99-x0)/xcols
x1= x0+xs
x2= x0+xs+xs
x3= x0+xs+xs+xs
x4= x0+xs+xs+xs+xs
x5= x0+xs+xs+xs+xs+xs




;--------------------------------------
; Top Row, Left-most Plot
;--------------------------------------

xaxistitle= "Time (Gyr)"
xmax = 6.6
xmin = 0.0

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
ymax = 8.0
ymin = 0

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle

	; 0
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[0], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100

	;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0 


        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100

        sfri= sfr1+sfr2
        lastn= n_elements(sfri)

	;-------------------------------------------

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)] 
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
	oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0

	;-------------------------------------------

        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 0.3, 5.76, R_peris[0], /data, charthick=4, size=1.33, color=0


        xyouts, 0.3, 6.86, msg, /data, charthick=4, size=1.33, color=0



; Plot #2
;--------------------------------------

!p.position= [x1, y1, x2, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtickformat='(a1)', /noerase


        ; 1
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[1], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


	;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0 
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 1.3, 5.76, R_peris[1], /data, charthick=4, size=1.33, color=0
	xyouts, 3.8, 6.86, 'R!Dperi!N=', /data, charthick=4, size=1.33, color=0



; Plot #3
;--------------------------------------

!p.position= [x2, y1, x3, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtickformat='(a1)', /noerase


        ; 1
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[2], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 1.3, 5.76, R_peris[2], /data, charthick=4, size=1.33, color=0




; Plot #4
;--------------------------------------

!p.position= [x3, y1, x4, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtickformat='(a1)', /noerase


        ; 1
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[3], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 1.3, 5.76, R_peris[3], /data, charthick=4, size=1.33, color=0




; Plot #5
;--------------------------------------

xmax = 10.6

!p.position= [x4, y1, x5, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtickformat='(a1)', /noerase


        ; 1
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[4], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 1.3, 5.76, R_peris[4], /data, charthick=4, size=1.33, color=0





;--------------------------------------------------------------------------------------------




;--------------------------------------
; Bottom Row, Left-most Plot
;--------------------------------------

;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

xmax = 6.6
xmin = 0.0

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtitle=' ', ytitle=yaxistitle, /noerase

        minor_open_sfr_file, G3fruns[0], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0

	gasc= 1.0-min(sfrm)
	gasc= strcompress(string(100.0*gasc),/remove_all)
	gasc= strmid(gasc,0,4)+'%'
        xyouts, 2.3, 0.26, gasc, /data, charthick=4, size=1.33, color=0



; Plot #2
;--------------------------------------

!p.position= [x1, y0, x2, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=' ', /noerase

        minor_open_sfr_file, G3fruns[1], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0

	gasc= 1.0-min(sfrm)
	gasc= strcompress(string(100.0*gasc),/remove_all)
	gasc= strmid(gasc,0,4)+'%'
        xyouts, 2.3, 0.26, gasc, /data, charthick=4, size=1.33, color=0



; Plot #3
;--------------------------------------

!p.position= [x2, y0, x3, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=xaxistitle, /noerase

        minor_open_sfr_file, G3fruns[2], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0

	gasc= 1.0-min(sfrm)
	gasc= strcompress(string(100.0*gasc),/remove_all)
	gasc= strmid(gasc,0,4)+'%'
        xyouts, 2.3, 0.26, gasc, /data, charthick=4, size=1.33, color=0




; Plot #4
;--------------------------------------

!p.position= [x3, y0, x4, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=' ', /noerase

        minor_open_sfr_file, G3fruns[3], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0

	gasc= 1.0-min(sfrm)
	gasc= strcompress(string(100.0*gasc),/remove_all)
	gasc= strmid(gasc,0,4)+'%'
        xyouts, 2.3, 0.26, gasc, /data, charthick=4, size=1.33, color=0



; Plot #5
;--------------------------------------

xmax = 10.6

!p.position= [x4, y0, x5, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=' ', /noerase

        minor_open_sfr_file, G3fruns[4], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0

	gasc= 1.0-min(sfrm)
	gasc= strcompress(string(100.0*gasc),/remove_all)
	gasc= strmid(gasc,0,4)+'%'
        xyouts, 2.3, 0.26, gasc, /data, charthick=4, size=1.33, color=0



;--------------------------------------
;--------------------------------------




; done
; ----------
device, /close

end














;================================================================================================================
;
;
;      |----------------|----------------|----------------|----------------|----------------|----------------|
;      |                |                |                |                |                |                |
;      |                |                |                |                |                |                |
;      |      SFR       |      SFR       |      SFR       |      SFR       |      SFR       |      SFR       |
;      |                |                |                |                |                |                |
;      |                |                |                |                |                |                |
;      |----------------|----------------|----------------|----------------|----------------|----------------|
;      |                |                |                |                |                |                |
;      |                |                |                |                |                |                |
;      |  GasConsump    |   GasConsump   |   GasConsump   |   GasConsump   |   GasConsump   |   GasConsump   |
;      |                |                |                |                |                |                |
;      |                |                |                |                |                |                |
;      |----------------|----------------|----------------|----------------|----------------|----------------|
;
;
;
;
;================================================================================================================


pro msfr3, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr3, junk"
   print, "  "
   print, "  "
   return
endif


filename='minorsfr_row3.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=24, newysize=10



;-------------------------------------------
;  set variables
;-------------------------------------------


Rtitle= "Satellite Inclination (degrees)"

lastpanelxmax = 0.0

inclinations= ['0!Eo!N', '30!Eo!N', '60!Eo!N', '90!Eo!N', '150!Eo!N', '180!Eo!N']  ; theta

;msg1= 'Mass ratio 1:2.3'
msg= 'G3G2'
G3fruns= ["G3G2a-u1", "G3G2-u3", "G3G2c-u1", "G3G2b-u1", "G3G2r-u3", "G3G2d-u1"]
fruniso1="G3il-u1a"
fruniso2="G2im-u1a"


;msg2= 'Mass ratio 1:5.8'
;msg= 'G3G1'
;G3fruns= ["G3G1a-u1", "G3G1-u3", "G3G1c-u1", "G3G1b-u1", "G3G1r-u3", "G3G1d-u1"]
;fruniso1="G3il-u1a"
;fruniso2="G1i-u1a"



;--------------------------------------
;--------------------------------------


;Rtitle= "Pericentric Distance"

;lastpanelxmax = 10.6

;inclinations= ['circular', 'R!Dperi!N= 6.8', '13.6', '27.2', '64.4']  ; R_peri (same for both)
;inclinations= ['1.7', '3.4', '6.8', '13.6', '27.2', '64.4']  ; R_peri (same for both)

;msg= 'G3G2'
;G3fruns= ["G3G2j-u1", "G3G2i-u1", "G3G2e-u1", "G3G2-u3", "G3G2f-u1", "G3G2g-u1"]
;oG3fruns= ["G3G2j-u2", "G3G2i-u2"]
;fruniso1="G3il-u1a"
;fruniso2="G2im-u1a"


;msg= 'Mass ratio 1:5.8'
;msg= 'G3G1'
;G3fruns= ["G3G1j-u1", "G3G1i-u1", "G3G1e-u1", "G3G1-u3", "G3G1f-u1", "G3G1g-u1"]
;oG3fruns= ["G3G1j-u2", "G3G1i-u2"]
;fruniso1="G3il-u1a"
;fruniso2="G1i-u1a"


; note that the last plot needs to be manually fixed
; to xmax=10.0



;--------------------------------------
;--------------------------------------



y0= 0.14
y99= 0.98
ycols= 2
ys= (y99-y0)/ycols
y1= y0+ys
y2= y0+ys+ys

x0= 0.08
x99= 0.99
xcols= 6
xs= (x99-x0)/xcols
x1= x0+xs
x2= x0+xs+xs
x3= x0+xs+xs+xs
x4= x0+xs+xs+xs+xs
x5= x0+xs+xs+xs+xs+xs
x6= x0+xs+xs+xs+xs+xs+xs




;--------------------------------------
; Top Row, Left-most Plot
;--------------------------------------

xaxistitle= "Time (Gyr)"
xmax = 6.6
xmin = 0.0

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
ymax = 8.0
ymin = 0

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle

	; 0
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[0], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100

	;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0 

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100

        sfri= sfr1+sfr2
        lastn= n_elements(sfri)

        ;-------------------------------------------

	tmax= max(timem)
        x= [xmin, timem, tmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)] 
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0

        ;-------------------------------------------


        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 0.3, 5.76, inclinations[0], /data, charthick=4, size=1.33, color=0


        xyouts, 0.3, 6.86, msg, /data, charthick=4, size=1.33, color=0


	if n_elements(oG3fruns) ge 1 then begin
        	minor_open_sfr_file, oG3fruns[0], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        	oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        	resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        	oplot, timem, sfrm, psym=-3, color= 50, linestyle= 2, thick= 6.0
	endif



; Plot #2
;--------------------------------------

!p.position= [x1, y1, x2, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtickformat='(a1)', /noerase


        ; 1
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[1], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100

        ;-------------------------------------------

	tmax= max(timem)
        x= [xmin, timem, tmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)] 
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0

        ;-------------------------------------------

	;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0 
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150

	xyouts, 0.9, 6.86, 'satellite', /data, charthick=4, size=1.33, color=0
	xyouts, 0.9, 5.76, 'inclination', /data, charthick=4, size=1.33, color=0
        xyouts, 1.3, 4.66, inclinations[1], /data, charthick=4, size=1.33, color=0

	;xyouts, 0.3, 6.86, 'R!Dperi!N=', /data, charthick=4, size=1.33, color=0
        ;xyouts, 1.3, 5.76, inclinations[1], /data, charthick=4, size=1.33, color=0

	if n_elements(oG3fruns) ge 2 then begin
        	minor_open_sfr_file, oG3fruns[1], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        	oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        	resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        	oplot, timem, sfrm, psym=-3, color= 50, linestyle= 2, thick= 6.0
	endif



; Plot #3
;--------------------------------------

!p.position= [x2, y1, x3, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtickformat='(a1)', /noerase


        ; 1
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[2], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100

        ;-------------------------------------------

	tmax= max(timem)
        x= [xmin, timem, tmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)] 
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0

        ;-------------------------------------------

        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 0.7, 5.76, inclinations[2], /data, charthick=4, size=1.33, color=0




; Plot #4
;--------------------------------------

!p.position= [x3, y1, x4, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtickformat='(a1)', /noerase


        ; 1
        ;-------------------------------------------
        minor_open_sfr_file, G3fruns[3], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100

        ;-------------------------------------------

	tmax= max(timem)
        x= [xmin, timem, tmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)] 
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0

        ;-------------------------------------------

        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 1.3, 5.76, inclinations[3], /data, charthick=4, size=1.33, color=0




; Plot #5
;--------------------------------------

!p.position= [x4, y1, x5, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtickformat='(a1)', /noerase

        minor_open_sfr_file, G3fruns[4], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100

        ;-------------------------------------------

	tmax= max(timem)
        x= [xmin, timem, tmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)] 
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0

        ;-------------------------------------------
        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        xyouts, 1.3, 5.76, inclinations[4], /data, charthick=4, size=1.33, color=0




; Plot #6
;--------------------------------------

if lastpanelxmax gt 0 then xmax= 10.6

!p.position= [x5, y1, x6, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtickformat='(a1)', /noerase

        minor_open_sfr_file, G3fruns[5], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100

        ;-------------------------------------------

	tmax= max(timem)
        x= [xmin, timem, tmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)] 
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
	time1= [time1, tmax]
	sfri= [sfri, sfri[lastn-1]]
        oplot, time1, sfri, psym=-3, color= 100, linestyle= 0, thick= 9.0

        ;-------------------------------------------
        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        xyouts, 0.9, 5.76, inclinations[5], /data, charthick=4, size=1.33, color=0






;--------------------------------------------------------------------------------------------




;--------------------------------------
; Bottom Row, Left-most Plot
;--------------------------------------

yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

xmax = 6.6
xmin = 0.0

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtitle=' ', ytitle=yaxistitle, /noerase

        minor_open_sfr_file, G3fruns[0], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100 

        sfri= sfr1+sfr2
        sfri= sfri / sfri[0]
        lastn= n_elements(sfri)

        ;-------------------------------------------

	tmax= max(timem)
        x= [xmin, timem, tmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= '!8e!6='+strmid(belbl,0,4)
        xyouts, 2.5, 0.14, belbl, /data, charthick=4, size=1.33, color=0

	;-------------------------------------------

	;gasc= 1.0-min(sfrm)
	;gasc= strcompress(string(100.0*gasc),/remove_all)
	;gasc= strmid(gasc,0,4)+'%'
        ;xyouts, 2.3, 0.26, gasc, /data, charthick=4, size=1.33, color=150

	if n_elements(oG3fruns) ge 1 then begin
        	minor_open_sfr_file, oG3fruns[0], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        	oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        	resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        	oplot, timem, sfrm, psym=-3, color= 50, linestyle= 0, thick= 6.0

        	gasc= 1.0-min(sfrm)
        	gasc= strcompress(string(100.0*gasc),/remove_all)
        	gasc= strmid(gasc,0,4)+'%'
        	xyouts, 2.3, 0.12, gasc, /data, charthick=4, size=1.33, color=50
	endif



; Plot #2
;--------------------------------------

!p.position= [x1, y0, x2, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=' ', /noerase

        minor_open_sfr_file, G3fruns[1], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0

        ;-------------------------------------------

	tmax= max(timem)
        x= [xmin, timem, tmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= '!8e!6='+strmid(belbl,0,4)
        xyouts, 2.5, 0.14, belbl, /data, charthick=4, size=1.33, color=0

        ;-------------------------------------------

	;gasc= 1.0-min(sfrm)
	;gasc= strcompress(string(100.0*gasc),/remove_all)
	;gasc= strmid(gasc,0,4)+'%'
        ;xyouts, 2.3, 0.26, gasc, /data, charthick=4, size=1.33, color=150

	if n_elements(oG3fruns) ge 2 then begin
        	minor_open_sfr_file, oG3fruns[1], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        	oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        	resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        	oplot, timem, sfrm, psym=-3, color= 50, linestyle= 0, thick= 6.0

        	gasc= 1.0-min(sfrm)
        	gasc= strcompress(string(100.0*gasc),/remove_all)
        	gasc= strmid(gasc,0,4)+'%'
        	xyouts, 2.3, 0.12, gasc, /data, charthick=4, size=1.33, color=50
	endif



; Plot #3
;--------------------------------------

!p.position= [x2, y0, x3, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=xaxistitle, /noerase

        minor_open_sfr_file, G3fruns[2], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0

        ;-------------------------------------------

	tmax= max(timem)
        x= [xmin, timem, tmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= '!8e!6='+strmid(belbl,0,4)
        xyouts, 2.5, 0.14, belbl, /data, charthick=4, size=1.33, color=0

        ;-------------------------------------------

	;gasc= 1.0-min(sfrm)
	;gasc= strcompress(string(100.0*gasc),/remove_all)
	;gasc= strmid(gasc,0,4)+'%'
        ;xyouts, 2.3, 0.26, gasc, /data, charthick=4, size=1.33, color=0




; Plot #4
;--------------------------------------

!p.position= [x3, y0, x4, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=' ', /noerase

        minor_open_sfr_file, G3fruns[3], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0

        ;-------------------------------------------

	tmax= max(timem)
        x= [xmin, timem, tmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= '!8e!6='+strmid(belbl,0,4)
        xyouts, 2.5, 0.14, belbl, /data, charthick=4, size=1.33, color=0

        ;-------------------------------------------

	;gasc= 1.0-min(sfrm)
	;gasc= strcompress(string(100.0*gasc),/remove_all)
	;gasc= strmid(gasc,0,4)+'%'
        ;xyouts, 2.3, 0.26, gasc, /data, charthick=4, size=1.33, color=0



; Plot #5
;--------------------------------------

!p.position= [x4, y0, x5, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=' ', /noerase

        minor_open_sfr_file, G3fruns[4], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0

        ;-------------------------------------------

	tmax= max(timem)
        x= [xmin, timem, tmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= '!8e!6='+strmid(belbl,0,4)
        xyouts, 2.5, 0.14, belbl, /data, charthick=4, size=1.33, color=0

        ;-------------------------------------------


	;gasc= 1.0-min(sfrm)
	;gasc= strcompress(string(100.0*gasc),/remove_all)
	;gasc= strmid(gasc,0,4)+'%'
        ;xyouts, 2.3, 0.26, gasc, /data, charthick=4, size=1.33, color=0




; Plot #6
;--------------------------------------

if lastpanelxmax gt 0 then xmax= 10.6

!p.position= [x5, y0, x6, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ytickformat='(a1)', xtitle=' ', /noerase

        minor_open_sfr_file, G3fruns[5], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
        ;oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0

        ;-------------------------------------------

	tmax= max(timem)
        x= [xmin, timem, tmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0


        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        ;oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100

	time1= [time1, tmax]
	sfri= [sfri, sfri[lastn-1]]
        oplot, time1, sfri, psym=-3, color= 100, linestyle= 0, thick= 6.0


        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= '!8e!6='+strmid(belbl,0,4)
        xyouts, 2.5, 0.14, belbl, /data, charthick=4, size=1.33, color=0

        ;-------------------------------------------

        ;gasc= 1.0-min(sfrm)
        ;gasc= strcompress(string(100.0*gasc),/remove_all)
        ;gasc= strmid(gasc,0,4)+'%'
        ;xyouts, 2.3, 0.26, gasc, /data, charthick=4, size=1.33, color=0





;--------------------------------------
;--------------------------------------




; done
; ----------
device, /close

end









