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








;==================================================================
;
;
;      |----------------------|
;      |                      |
;      |                      |
;      |                      |
;      |      SFR             |
;      |                      |
;      |                      |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;      |                      |
;      |                      |
;      |      Gas Mass        |
;      |                      |
;      |                      |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;
;
;
;
;======================================================================


pro msfr, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr, junk"
   print, "  "
   print, "  "
   return
endif


filename='minorsfr4.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=10, newysize=16
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------

xaxistitle = "!6Time (Gyr)"
xmax = 6.0
xmin = 0


x0= 0.18
x1= 0.97

y0= 0.10
y1= 0.545
y2= 0.99

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 80.0
;ymax = 20.0
ymax = 13.0
;ymax = 4.0
ymin = 0

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle


        ;-------------------------------------------
	;process_one_sfr, "G3G2-u3", 150, '!8n2med!6', 17.6, 0.91, lstyle= 0
	;process_one_sfr, "G3G2-u4", 50, '!8n0med!6', 15.7, 0.89, lstyle= 1
	;xyouts, 0.28, 0.17, "G3G2, 1:2.3 merger", /normal, charthick=4, size=1.33, color=0


        ;-------------------------------------------
	;process_one_sfr, "G3G1-u3", 150, '!8n2med!6', 7.1, 0.93, lstyle= 0
	;process_one_sfr, "G3G1-u4", 50, '!8n0med!6', 6.3, 0.89, lstyle= 1
	;xyouts, 0.28, 0.17, "G3G1, 1:5.8 merger", /normal, charthick=4, size=1.33, color=0

        ;-------------------------------------------
	;process_one_sfr, "G3G0e-u3", 150, '!8n2med!6', 3.5, 0.93, lstyle= 0
	;process_one_sfr, "G3G0e-u4", 50, '!8n0med!6', 3.15, 0.89, lstyle= 1
	;xyouts, 0.28, 0.17, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=0

        ;-------------------------------------------
	;process_one_sfr, "G3blv5G3blv5-u1", 0, 'B/D=0', 18.0, 0.94, lstyle= 0
	;process_one_sfr, "G3G3b-u1", 50, 'B/D=0.2', 16.7, 0.91, lstyle= 1
	;process_one_sfr, "G3BT1G3BT1-u1", 100, 'B/D=0.25', 15.3, 0.88, lstyle= 2
	;process_one_sfr, "G3BT2G3BT2-u1", 150, 'B/D=0.5', 14.0, 0.85, lstyle= 3
	;process_one_sfr, "G3BT3G3BT3-u1", 200, 'B/D=1.0', 12.7, 0.82, lstyle= 4
	;xyouts, 0.28, 0.17, "G3G3, 1:1 merger", /normal, charthick=4, size=1.33, color=0

        ;-------------------------------------------
	;process_one_sfr, "G3blv5G1-u3", 0, 'B/D=0', 18.0, 0.94, lstyle= 0
	;process_one_sfr, "G3G1-u3", 50, 'B/D=0.2', 16.7, 0.91, lstyle= 1
	;process_one_sfr, "G3BT1G1-u1", 100, 'B/D=0.25', 15.3, 0.88, lstyle= 2
	;process_one_sfr, "G3BT2G1-u1", 150, 'B/D=0.5', 14.0, 0.85, lstyle= 3
	;process_one_sfr, "G3BT3G1-u1", 200, 'B/D=1.0', 12.7, 0.82, lstyle= 4
	;xyouts, 0.28, 0.17, "G3G1, 5.8:1 merger", /normal, charthick=4, size=1.33, color=0

        ;-------------------------------------------
	;process_one_sfr, "G3G3b-u1", 0, '!8f!6= 0.2', 72.0, 0.94, lstyle= 0
	;process_one_sfr, "G3gf1G3gf1b-u1", 50, '!8f!6= 0.5', 67.0, 0.91, lstyle= 1
	;process_one_sfr, "G3gf2G3gf2b-u1", 150, '!8f!6= 0.78', 61.0, 0.88, lstyle= 2
	;xyouts, 0.28, 0.17, "G3G3, 1:1 merger", /normal, charthick=4, size=1.33, color=0
	;process_one_sfr, "G3G1-u1", 0, '!8f!6= 0.2', 72.0, 0.94, lstyle= 0
	;process_one_sfr, "G3gf1G1-u1", 50, '!8f!6= 0.5', 67.0, 0.91, lstyle= 1
	;process_one_sfr, "G3gf2G2-u1", 150, '!8f!6= 0.78', 61.0, 0.88, lstyle= 2
	;xyouts, 0.28, 0.17, "G3G1, 5.8:1 merger", /normal, charthick=4, size=1.33, color=0

        ;-------------------------------------------
        ;process_one_sfr, "G3Rd4eG1", 0, 'B/D=0', 18.0, 0.94, lstyle= 0
        ;process_one_sfr, "G3bd3G1", 50, 'B/D=0.4', 16.7, 0.91, lstyle= 1
	;process_one_sfr, "G3Rd4eG1_n0", 100, 'B/D=0 (n0)', 15.3, 0.88, lstyle= 2
        ;process_one_sfr, "G3bd3G1_n0", 150, 'B/D=0.4 (n0)', 14.0, 0.85, lstyle= 3

        ;-------------------------------------------
        ;process_one_sfr, "G3Rd4eG2", 0, 'no R!Dcutoff!N', 18.0, 0.94, lstyle= 0
        ;process_one_sfr, "G3Rd4fG2", 50, 'R!Dcutoff!N=40 kpc', 16.7, 0.91, lstyle= 1
	;process_one_sfr, "G3Rd4gG2", 100, 'R!Dcutoff!N=30 kpc', 15.3, 0.88, lstyle= 2
        ;process_one_sfr, "G3Rd4hG2", 150, 'R!Dcutoff!N=20 kpc', 14.0, 0.85, lstyle= 3

	;process_one_sfr, "G3Rd4G2", 50, '!9a!6= 1', 16.7, 0.91, lstyle= 0

        ;-------------------------------------------
	;xyouts, 0.29, 'B/D', 
        process_one_sfr, "SbIm_o2_n0", 0, 'B/D= 0', 11.8, 0.94, lstyle= 0
        process_one_sfr, "Sbbd1Im_o2_n0", 50, 'B/D= 0.2', 11.0, 0.91, lstyle= 1
	process_one_sfr, "Sbbd2Im_o2_n0", 150, 'B/D= 0.5', 10.1, 0.88, lstyle= 2
        ;xyouts, 0.34, 0.14, '8:1 merger', /normal, charthick=4, size=1.33, color=0

        ;process_one_sfr, "SbG1_o2_n0", 50, 'B/D= 0', 11.0, 0.91, lstyle= 1



;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle


        ;-------------------------------------------
	;process_one_gasmass, "G3G2-u3", 150, lstyle=0
	;process_one_gasmass, "G3G2-u4", 50, lstyle=1

        ;-------------------------------------------
	;process_one_gasmass, "G3G1-u3", 150, lstyle=0
	;process_one_gasmass, "G3G1-u4", 50, lstyle=1

        ;-------------------------------------------
	;process_one_gasmass, "G3G0e-u3", 150, lstyle=0
	;process_one_gasmass, "G3G0e-u4", 50, lstyle=1

        ;-------------------------------------------
	;process_one_gasmass, "G3blv5G3blv5-u1", 0, lstyle=0
	;process_one_gasmass, "G3G3b-u1", 50, lstyle=1
	;process_one_gasmass, "G3BT1G3BT1-u1", 100, lstyle=2
	;process_one_gasmass, "G3BT2G3BT2-u1", 150, lstyle=3
	;process_one_gasmass, "G3BT3G3BT3-u1", 200, lstyle=4

        ;-------------------------------------------
        ;process_one_gasmass, "G3blv5G1-u3", 0, lstyle=0
        ;process_one_gasmass, "G3G1-u1", 50, lstyle=1
        ;process_one_gasmass, "G3BT1G1-u1", 100, lstyle=2
        ;process_one_gasmass, "G3BT2G1-u1", 150, lstyle=3
        ;process_one_gasmass, "G3BT3G1-u1", 200, lstyle=4

        ;-------------------------------------------
	;process_one_gasmass, "G3G3b-u1", 0, lstyle=0
	;process_one_gasmass, "G3gf1G3gf1b-u1", 50, lstyle=1
	;process_one_gasmass, "G3gf2G3gf2b-u1", 150, lstyle=2
	;process_one_gasmass, "G3G1-u1", 0, lstyle=0
	;process_one_gasmass, "G3gf1G1-u1", 50, lstyle=1
	;process_one_gasmass, "G3gf2G2-u1", 150, lstyle=2

        ;-------------------------------------------
        ;process_one_gasmass, "G3Rd4eG1", 0, lstyle=0
        ;process_one_gasmass, "G3bd3G1", 50, lstyle=1
	;process_one_gasmass, "G3Rd4eG1_n0", 100, lstyle=2
        ;process_one_gasmass, "G3bd3G1_n0", 150, lstyle=3

        ;-------------------------------------------
        ;process_one_gasmass, "G3Rd4eG2", 0, lstyle=0
        ;process_one_gasmass, "G3Rd4fG2", 50, lstyle=1
	;process_one_gasmass, "G3Rd4gG2", 100, lstyle=2
        ;process_one_gasmass, "G3Rd4hG2", 150, lstyle=3

        ;process_one_gasmass, "G3Rd4G2", 50, lstyle=1

        ;-------------------------------------------
        process_one_gasmass, "SbIm_o2_n0", 0, lstyle=0
        process_one_gasmass, "Sbbd1Im_o2_n0", 50, lstyle=1
	process_one_gasmass, "Sbbd2Im_o2_n0", 150, lstyle=2

        ;process_one_gasmass, "SbG1_o2_n0", 50, lstyle=1


;--------------------------------------
;--------------------------------------




; done
; ----------
device, /close

end











;-------------------------------------------------------



pro process_one_sfr, frun, thiscolor, msg, ypt, ymsg, lstyle=lstyle

        ;-------------------------------------------
        ;
        ;  #1
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        ;resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=30

        oplot, time1, sfr1, psym=-3, color= thiscolor, linestyle= lstyle, thick= 6.0
        oplot, [0.3, 0.7], [ypt, ypt], psym=-3, color= thiscolor, linestyle= lstyle, thick= 6.0
        xyouts, 0.28, ymsg, msg, /normal, charthick=4, size=1.33, color=thiscolor

end





pro process_one_gasmass, frun, thiscolor, lstyle=lstyle

        ;-------------------------------------------
        ;
        ;  Merger - frunmerg
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        ;resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        ;resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=50
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=30

        oplot, time1, sfr1, psym=-3, color= thiscolor, linestyle= lstyle, thick= 6.0


end




