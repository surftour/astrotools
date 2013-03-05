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

print, "check sfr resampling"
print, "total int of sfr (the following two should be equal)"
print, timemax/n_old * total(oldsfr)
print, dt * total(newsfr)

end










;======================================================================
;
;
;
;
;
;======================================================================
pro process_one_minor_sfr, frun, lbl, lcolor, ylbl, ylbl_key, old=old, xmax=xmax, xlbl=xlbl


nsmths= 100

	if not keyword_set(old) then begin
        	minor_open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
	endif else begin
		open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
	endelse
        oldtime= transpose(sfrtime)
        oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=nsmths

	thisthick= 6.0
	if lcolor eq 0   then thisthick= 14.0
	if lcolor eq 50  then thisthick= 10.0
	if lcolor eq 150 then thisthick= 6.0
	if lcolor eq 200 then thisthick= 2.0

	if lbl ne 'none' then begin
        	;oplot, timem, sfrm, psym=-3, color= lcolor, linestyle= 0, thick= 6.0
        	oplot, timem, sfrm, psym=-3, color= lcolor, linestyle= 0, thick= thisthick
        	;oplot, [xmax*0.548,xmax*0.613], [ylbl_key,ylbl_key], psym=-3, color= lcolor, linestyle= 0, thick=6.0

		; x0 is 0.7 for most plots, however, we're moving this
		; for the new figure 7
        	;xyouts, 0.7, ylbl, lbl, /normal, charthick=4, size=1.33, color=lcolor
        	xyouts, xlbl, ylbl, lbl, /normal, charthick=4, size=1.33, color=lcolor
	endif else begin
        	oplot, timem, sfrm, psym=-3, color= lcolor, linestyle= 1, thick= 6.0
	endelse

end






pro process_one_sfr_atsnaptimes, frun, lbl, lcolor, ylbl, ylbl_key, old=old, xmax=xmax


nsmths= 100

	read_sfratsnapt_file, frun, sfrtime, sfrtot, sfr1, sfr2
        oldtime= transpose(sfrtime)
        oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=nsmths

        if lbl ne 'none' then begin
                oplot, timem, sfrm, psym=-3, color= lcolor, linestyle= 0, thick= 6.0
                ;oplot, [xmax*0.548,xmax*0.613], [ylbl_key,ylbl_key], psym=-3, color= lcolor, linestyle= 0, thick=6.0
                xyouts, 0.7, ylbl, lbl, /normal, charthick=4, size=1.33, color=lcolor
        endif else begin
                oplot, timem, sfrm, psym=-3, color= lcolor, linestyle= 1, thick= 6.0
        endelse

end







pro process_one_sfr, frun, thiscolor, msg, ypt, ymsg, lstyle=lstyle


nsmths= 100

        ;-------------------------------------------
        ;
        ;  #1
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=nsmths

        oplot, time1, sfr1, psym=-3, color= thiscolor, linestyle= lstyle, thick= 6.0
        oplot, [0.3, 0.7], [ypt, ypt], psym=-3, color= thiscolor, linestyle= lstyle, thick= 6.0
        xyouts, 0.28, ymsg, msg, /normal, charthick=4, size=1.33, color=thiscolor

end





pro process_one_gasmass, frun, thiscolor, lstyle=lstyle


nsmths= 100

        ;-------------------------------------------
        ;
        ;  Merger - frunmerg
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=nsmths

        oplot, time1, sfr1, psym=-3, color= thiscolor, linestyle= lstyle, thick= 6.0


end














;======================================================================
;      PLOT FORMAT
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
;      |----------------------|
;
;
;
;
;
;     Fig. 7 in minors paper
;
;
;
;======================================================================

pro msfr_f7, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr_f7, junk"
   print, "  "
   print, "  "
   return
endif



;-------------------------------------------
;  set variables
;-------------------------------------------

frun1= "-1"
frun2= "-1"
frun3= "-1"
frun4= "-1"

frunold1= "-1"
frunold2= "-1"
frunold3= "-1"
frunold4= "-1"

; ------------------------------------------
;
;   SFR figure in Paper

;filename='minorsfr_g3.eps' & ymax= 15.0 & xmax= 6.2
;frun1= "G3G3b-u1" & lbl1= "G3G3, 1:1"
;frun2= "G3G2-u3"  & lbl2= "G3G2, 2.3:1"
;frun3= "G3G1-u3"  & lbl3= "G3G1, 5.8:1"
;frun4= "G3G0e-u3" & lbl4= "G3G0, 22.7:1"
;filename='minorsfr_g3.eps' & ymax= 28.0 & xmax= 6.2
filename='minorsfr_f7a.eps' & ymax= 28.0 & xmax= 6.2
frun1= "G3G3b-u2" & lbl1= "G3G3, 1:1"
frun2= "G3G2-u4"  & lbl2= "G3G2, 2.3:1"
frun3= "G3G1-u4"  & lbl3= "G3G1, 5.8:1"
frun4= "G3G0e-u4" & lbl4= "G3G0, 22.7:1"

;filename='minorsfr_g2.eps' & ymax= 12.0 & xmax= 4.2
;frun1= "G2G2-u1" & lbl1= "G2G2, 1:1"
;frun2= "G2G1-u3" & lbl2= "G2G1, 2.6:1"
;frun3= "G2G0-u3" & lbl3= "G2G0, 10.0:1"
;filename='minorsfr_g2.eps' & ymax= 18.0 & xmax= 4.2
;frun1= "G2G2-u2" & lbl1= "G2G2, 1:1"
;frun2= "G2G1-u4" & lbl2= "G2G1, 2.6:1"
;frun3= "G2G0-u4" & lbl3= "G2G0, 10.0:1"

;filename='minorsfr_g1.eps' & ymax= 5.0 & xmax= 4.2
;frun1= "G1G1a-u1" & lbl1= "G1G1, 1:1"
;frun2= "G1G0-u3"  & lbl2= "G1G0, 3.9:1"
;filename='minorsfr_g1.eps' & ymax= 13.0 & xmax= 4.2
;frun1= "G1G1a-u2" & lbl1= "G1G1, 1:1"
;frun2= "G1G0-u4"  & lbl2= "G1G0, 3.9:1"

;filename='minorsfr_g0.eps' & ymax= 1.2 & xmax= 4.2
;frun1= "G0G0a-u1" & lbl1= "G0G0, 1:1"
;filename='minorsfr_g0.eps' & ymax= 3.9 & xmax= 4.2
;filename='msfr_f7d.eps' & ymax= 3.9 & xmax= 4.2
;frun1= "G0G0a-u2" & lbl1= "G0G0, 1:1"

; ------------------------------------------
;
;     Various Test isolated galaxies

;filename='minorsfr_test.eps' & ymax= 2.0 & xmax= 6.0
;frun1= "G3il-u1a" & lbl1= "G3 (upsand)"
;frunold2= "isolated/G3"  & lblold= "G3 (sauron)"
;frun1= "G3gf1i-u1" & lbl1= "G3gf1 (upsand)"
;frunold2= "isolated/G3gf1"  & lblold= "G3gf1 (sauron)"
;
;filename='minorsfr_test.eps' & ymax= 3.0 & xmax= 6.0
;frunold1= "isolated/G3"  & lbl1= "G3 (n2med)"
;frunold2= "isolated/G3_n0"  & lbl2= "G3_n0 (n0med)"
;frunold3= "isolated/G3gf1"  & lbl3= "G3gf1 (n2med)"
;frunold4= "isolated/G3gf1_n0"  & lbl4= "G3gf1_n0 (n0med)"
;
;filename='minorsfr_test.eps' & ymax= 3.0 & xmax= 6.0
;frunold1= "isolated/G3"  & lbl1= "G3 (n2med)"
;frunold2= "isolated/G3_n0"  & lbl2= "G3_n0 (n0med)"
;frunold3= "isolated/G3bl"  & lbl3= "G3bl (n2med)"
;frunold4= "isolated/G3bl_n0"  & lbl4= "G3bl_n0 (n0med)"



; -------------------------------------------------


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



; ---------------------------------------------------
; ---------------------------------------------------
;
;  Print the Shit
;
; ---------------------------------------------------
; ---------------------------------------------------


yaxistitle="!6SFR (M!D!9n!6!N Yr!E-1!N)"
xaxistitle = "!6Time (Gyr)"

xmin = 0.0
ymin = 0.0

;nsmths= 70 , is now set in the processing


; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]
;!p.font= 0

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata



;-------------------------------------------
;
;-------------------------------------------
if frun1 ne "-1" then process_one_minor_sfr, frun1, lbl1, 0, 0.90, ymax*0.92, xmax=xmax
if frun2 ne "-1" then process_one_minor_sfr, frun2, lbl2, 50, 0.86, ymax*0.87, xmax=xmax
if frun3 ne "-1" then process_one_minor_sfr, frun3, lbl3, 150, 0.82, ymax*0.82, xmax=xmax
if frun4 ne "-1" then process_one_minor_sfr, frun4, lbl4, 200, 0.78, ymax*0.77, xmax=xmax


;-------------------------------------------
;
;-------------------------------------------
if frunold1 ne "-1" then process_one_minor_sfr, frunold1, lbl1, 0, 0.90, ymax*0.92, /old, xmax=xmax
if frunold2 ne "-1" then process_one_minor_sfr, frunold2, lbl2, 66, 0.86, ymax*0.87, /old, xmax=xmax
if frunold3 ne "-1" then process_one_minor_sfr, frunold3, lbl3, 133, 0.82, ymax*0.82, /old, xmax=xmax
if frunold4 ne "-1" then process_one_minor_sfr, frunold4, lbl4, 200, 0.78, ymax*0.77, /old, xmax=xmax





;--------------------------------------
;--------------------------------------

device, /close


end







pro msfr_f7_log, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr_f7_log, junk"
   print, "  "
   print, "  "
   return
endif



;-------------------------------------------
;  set variables
;-------------------------------------------

frun1= "-1"
frun2= "-1"
frun3= "-1"
frun4= "-1"

frunold1= "-1"
frunold2= "-1"
frunold3= "-1"
frunold4= "-1"

; ------------------------------------------
;
;   SFR figure in Paper

;filename='minorsfr_g3.eps' & ymax= 15.0 & xmax= 6.2
;frun1= "G3G3b-u1" & lbl1= "G3G3, 1:1"
;frun2= "G3G2-u3"  & lbl2= "G3G2, 2.3:1"
;frun3= "G3G1-u3"  & lbl3= "G3G1, 5.8:1"
;frun4= "G3G0e-u3" & lbl4= "G3G0, 22.7:1"
;filename='minorsfr_g3.eps' & ymax= 28.0 & xmax= 6.2
filename='minorsfr_f7a_log.eps' & ymax= 50.0 & xmax= 6.2
frun1= "G3G3b-u2" & lbl1= "G3G3, 1:1"
frun2= "G3G2-u4"  & lbl2= "G3G2, 2.3:1"
frun3= "G3G1-u4"  & lbl3= "G3G1, 5.8:1"
frun4= "G3G0e-u4" & lbl4= "G3G0, 22.7:1"

;filename='minorsfr_g2.eps' & ymax= 12.0 & xmax= 4.2
;frun1= "G2G2-u1" & lbl1= "G2G2, 1:1"
;frun2= "G2G1-u3" & lbl2= "G2G1, 2.6:1"
;frun3= "G2G0-u3" & lbl3= "G2G0, 10.0:1"
;filename='minorsfr_g2.eps' & ymax= 18.0 & xmax= 4.2
;frun1= "G2G2-u2" & lbl1= "G2G2, 1:1"
;frun2= "G2G1-u4" & lbl2= "G2G1, 2.6:1"
;frun3= "G2G0-u4" & lbl3= "G2G0, 10.0:1"

;filename='minorsfr_g1.eps' & ymax= 5.0 & xmax= 4.2
;frun1= "G1G1a-u1" & lbl1= "G1G1, 1:1"
;frun2= "G1G0-u3"  & lbl2= "G1G0, 3.9:1"
;filename='minorsfr_g1.eps' & ymax= 13.0 & xmax= 4.2
;frun1= "G1G1a-u2" & lbl1= "G1G1, 1:1"
;frun2= "G1G0-u4"  & lbl2= "G1G0, 3.9:1"

;filename='minorsfr_g0.eps' & ymax= 1.2 & xmax= 4.2
;frun1= "G0G0a-u1" & lbl1= "G0G0, 1:1"
;filename='minorsfr_g0.eps' & ymax= 3.9 & xmax= 4.2
;filename='msfr_f7d.eps' & ymax= 3.9 & xmax= 4.2
;frun1= "G0G0a-u2" & lbl1= "G0G0, 1:1"

; ------------------------------------------
;
;     Various Test isolated galaxies

;filename='minorsfr_test.eps' & ymax= 2.0 & xmax= 6.0
;frun1= "G3il-u1a" & lbl1= "G3 (upsand)"
;frunold2= "isolated/G3"  & lblold= "G3 (sauron)"
;frun1= "G3gf1i-u1" & lbl1= "G3gf1 (upsand)"
;frunold2= "isolated/G3gf1"  & lblold= "G3gf1 (sauron)"
;
;filename='minorsfr_test.eps' & ymax= 3.0 & xmax= 6.0
;frunold1= "isolated/G3"  & lbl1= "G3 (n2med)"
;frunold2= "isolated/G3_n0"  & lbl2= "G3_n0 (n0med)"
;frunold3= "isolated/G3gf1"  & lbl3= "G3gf1 (n2med)"
;frunold4= "isolated/G3gf1_n0"  & lbl4= "G3gf1_n0 (n0med)"
;
;filename='minorsfr_test.eps' & ymax= 3.0 & xmax= 6.0
;frunold1= "isolated/G3"  & lbl1= "G3 (n2med)"
;frunold2= "isolated/G3_n0"  & lbl2= "G3_n0 (n0med)"
;frunold3= "isolated/G3bl"  & lbl3= "G3bl (n2med)"
;frunold4= "isolated/G3bl_n0"  & lbl4= "G3bl_n0 (n0med)"



; -------------------------------------------------


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



; ---------------------------------------------------
; ---------------------------------------------------
;
;  Print the Shit
;
; ---------------------------------------------------
; ---------------------------------------------------


yaxistitle="!6SFR (M!D!9n!6!N Yr!E-1!N)"
xaxistitle = "!6Time (Gyr)"

xmin = 0.0
ymin = 0.03

;nsmths= 70 , is now set in the processing


; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]
;!p.font= 0

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, ytickformat='exp_label', $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata



;-------------------------------------------
;
;-------------------------------------------
if frun1 ne "-1" then process_one_minor_sfr, frun1, lbl1, 0, 0.90, ymax*0.92, xmax=xmax
if frun2 ne "-1" then process_one_minor_sfr, frun2, lbl2, 50, 0.86, ymax*0.87, xmax=xmax
if frun3 ne "-1" then process_one_minor_sfr, frun3, lbl3, 150, 0.82, ymax*0.82, xmax=xmax
if frun4 ne "-1" then process_one_minor_sfr, frun4, lbl4, 200, 0.78, ymax*0.77, xmax=xmax


;-------------------------------------------
;
;-------------------------------------------
if frunold1 ne "-1" then process_one_minor_sfr, frunold1, lbl1, 0, 0.90, ymax*0.92, /old, xmax=xmax
if frunold2 ne "-1" then process_one_minor_sfr, frunold2, lbl2, 66, 0.86, ymax*0.87, /old, xmax=xmax
if frunold3 ne "-1" then process_one_minor_sfr, frunold3, lbl3, 133, 0.82, ymax*0.82, /old, xmax=xmax
if frunold4 ne "-1" then process_one_minor_sfr, frunold4, lbl4, 200, 0.78, ymax*0.77, /old, xmax=xmax





;--------------------------------------
;--------------------------------------

device, /close


end















;======================================================================
;      PLOT FORMAT
;
;
;      |----------------------|
;      |                      |
;      |                      |
;      |    Isolated          |
;      |      SFR             |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;
;
;
;
;======================================================================

pro isosfr, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "isosfr, junk"
   print, "  "
   print, "  "
   return
endif



;-------------------------------------------
;  set variables
;-------------------------------------------

;
;   Isolated SFR figure in Paper

filename='minorsfr_iso.eps'

ymax= 3.0
ymin= 6.0e-5

xmax= 1.0
xmin= 0.0


; -------------------------------------------------


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



; ---------------------------------------------------
; ---------------------------------------------------
;
;  Print the Shit
;
; ---------------------------------------------------
; ---------------------------------------------------


yaxistitle="!6SFR (M!D!9n!6!N Yr!E-1!N)"
xaxistitle = "!6Time (Gyr)"


;nsmths= 70 , is now set in the processing


; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;---------------------------

!p.position= [0.18, 0.15, 0.96, 0.98]
;!p.font= 0

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	ytickformat='exp_label', $
	/nodata



;-------------------------------------------
;   n2med
;-------------------------------------------
process_one_minor_sfr, 'G3il-u1a', 'G3', 200, 0.92, 1.2, xmax=xmax
process_one_minor_sfr, 'G2im-u1a', 'G2', 100, 0.84, 0.5, xmax=xmax
process_one_minor_sfr, 'G1i-u1a', 'G1', 150, 0.72, 0.02, xmax=xmax
process_one_minor_sfr, 'G0i-u1a', 'G0', 50, 0.48, 0.0006, xmax=xmax


;-------------------------------------------
;   n0med
;-------------------------------------------
process_one_minor_sfr, "isolated/G3_n0", 'none', 200, 0.0, 0.0, /old, xmax=xmax
process_one_minor_sfr, 'G2im-u2', 'none', 100, 0.0, 0.0, xmax=xmax
process_one_minor_sfr, 'G1i-u2', 'none', 150, 0.0, 0.0, xmax=xmax
process_one_minor_sfr, 'G0i-u2', 'none', 50, 0.0, 0.0, xmax=xmax





;--------------------------------------
;--------------------------------------

device, /close


end









;======================================================================
;      PLOT FORMAT
;
;
;      |----------------------|
;      |                      |
;      |     Merger SFR       |
;      |         &            |
;      |                      |
;      |     gal 1 SFR        |
;      |                      |
;      |----------------------|
;      |                      |
;      |                      |
;      |     gal 2 SFR        |
;      |                      |
;      |                      |
;      |----------------------|
;
;
;
;
;======================================================================

pro sfrg12, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "sfrg12, junk"
   print, "  "
   print, "  "
   return
endif



;-------------------------------------------
;  set variables
;-------------------------------------------

;
;   Isolated SFR figure in Paper

filename='minorsfrg12.eps'


xmax= 6.0
xmin= 0.0


; -------------------------------------------------


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newysize= 14.0, newxsize= 12.0



; ---------------------------------------------------
; ---------------------------------------------------
;
;  Print the Shit
;
; ---------------------------------------------------
; ---------------------------------------------------


yaxistitle="!6SFR (M!D!9n!6!N Yr!E-1!N)"
xaxistitle = "!6Time (Gyr)"


;nsmths= 40
nsmths= 100

x0= 0.16
x1= 0.98

y0= 0.12
y1= 0.55
y2= 0.98


;---------------------------

ymax= 12.0
ymin= 0.0

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
	;xtickformat='(a1)', ytitle=yaxistitle
	xtickformat='(a1)', ytitle=' '


	; ------------------------------------------------------------

	;frun= "G3G1-u4"

        ;minor_open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime)
        ;oldsfr= transpose(sfrsfr)
        ;resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=nsmths

        ;oplot, timem, sfrm, psym=-3, color= 0, linestyle= 0, thick= 6.0

        ; ------------------------------------------------------------

        ;frun= "G3_n0"

        ;minor_open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime)
        ;oldsfr= transpose(sfrsfr)
        ;resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=nsmths

        ;oplot, timem, sfrm, psym=-3, color= 0, linestyle= 0, thick= 6.0


	; ------------------------------------------------------------

        ;xyouts, 0.19, 0.93, "G3G1, 1:5.8 merger", /normal, charthick=4, size=1.5, color=0
        ;xyouts, 0.19, 0.93, "SbG1, 4.8:1 merger", /normal, charthick=4, size=1.5, color=0
        xyouts, 0.19, 0.93, "SbIm, 8:1 merger", /normal, charthick=4, size=1.5, color=0

        oplot, [0.3, 0.7], [3.2, 3.2], psym=-3, color= 50, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G2, 1:2.3 merger", /normal, charthick=4, size=1.33, color=50
        ;xyouts, 0.29, 0.88, "Total (both G3 and G1)", /normal, charthick=4, size=1.33, color=50
        ;xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=50
        ;xyouts, 0.29, 0.88, "Total (both Sb and G1)", /normal, charthick=4, size=1.33, color=50
        xyouts, 0.29, 0.88, "Total (both Sb and Im)", /normal, charthick=4, size=1.33, color=50

        oplot, [0.3, 0.7], [2.7, 2.7], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.84, "contribution from G3", /normal, charthick=4, size=1.33, color=100
        xyouts, 0.29, 0.84, "contribution from Sb", /normal, charthick=4, size=1.33, color=100

	; ------------------------------------------------------------




	;frun= "/raid12/tcox/Gs/G3G1-u4"
	;frun= "/raid2/tcox/G3/SbG1_o2_n0"
	frun= "/raid2/tcox/G3/SbIm_o2_n0"

        read_sfratsnapt_file, frun, sfrtime, sfrtot, sfr1, sfr2
        oldtime= transpose(sfrtime)
        oldsfr= transpose(sfrtot)
        oldsfr1= transpose(sfr1)
        oldsfr2= transpose(sfr2)

        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=nsmths
	oplot, timem, sfrm, psym=-3, color= 50, linestyle= 0, thick= 6.0

        resample_sfrhist, oldsfr1, oldtime, sfrm, timem, Nelements=nsmths
	oplot, timem, sfrm, psym=-3, color= 100, linestyle= 1, thick= 6.0




xyouts, 0.04, 0.43, yaxistitle, /normal, charthick= 3.0, size= 1.7, color= 0, orientation= 90.0


;---------------------------

ymax= 1.74
ymin= 0.0

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ;xtitle=xaxistitle, ytitle=yaxistitle, /noerase
        xtitle=xaxistitle, ytitle=' ', /noerase


        ; ------------------------------------------------------------

        resample_sfrhist, oldsfr2, oldtime, sfrm, timem, Nelements=nsmths
        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 2, thick= 6.0


        ; ------------------------------------------------------------

        ;frun= "G1i-u2"

        ;minor_open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime)
        ;oldsfr= transpose(sfrsfr)
        ;resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=nsmths

        ;oplot, timem, sfrm, psym=-3, color= 0, linestyle= 0, thick= 6.0


        ; ------------------------------------------------------------

        oplot, [0.3, 0.7], [0.6, 0.6], psym=-3, color= 150, linestyle= 2, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.46, "contribution from G1", /normal, charthick=4, size=1.33, color=150
        xyouts, 0.29, 0.46, "contribution from Im", /normal, charthick=4, size=1.33, color=150

        ; ------------------------------------------------------------


;--------------------------------------
;--------------------------------------

device, /close


end



















;======================================================================
;      PLOT FORMAT
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
;      |----------------------|
;
;
;
;
;======================================================================

pro msfr3, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr3, junk"
   print, "  "
   print, "  "
   return
endif


filename='minorsfr3.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------

frunmerg= "G3G2-u3"
fruniso1= "G3il-u1a"
fruniso2= "G2im-u1a"



; ---------------------------------------------------
; ---------------------------------------------------
;
;  Print the Shit
;
; ---------------------------------------------------
; ---------------------------------------------------


yaxistitle="!6SFR (M!D!9n!6!N Yr!E-1!N)"
;yaxistitle="!6SFR (M!D!9n!6!N Yr!E-1!N) / Mass (10!e10!n M!D!9n!6!N)"
;if keyword_set(normalize) then yaxistitle="!6Gas Mass / M!D0!N "
;xaxistitle = "Time (Gyr/h)"
xaxistitle = "!6Time (Gyr)"

xmax = 6.0
;xmax = 5.0
;xmax = 4.0
;xmax= 2.0
xmin = 0

;ymax = 60.0
;ymax = 20.0
;ymax = 15
;ymax = 12
ymax = 8.0
;ymax = 5.0
;ymax = 2.0
;ymax = 1.2
;ymax = 1.0
ymin = 0 


; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]
;!p.font= 0

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata



	;-------------------------------------------
	;
	;  Merger
	;
	;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass

	oldtime= transpose(sfrtime)
	oldsfr= transpose(sfrsfr)

	resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100
	
        ;oplot, time, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
	;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150


        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime)
        oldsfr= transpose(sfrsfr)
	resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        ;oplot, time, sfr, psym=-3, color= 50, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.75, "gal 1", /normal, charthick=4, size=1.33, color=50


        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime)
        oldsfr= transpose(sfrsfr)
	resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        ;oplot, time, sfr, psym=-3, color= 100, linestyle= 0, thick= 6.0
	;xyouts, 0.7, 0.70, "gal 2", /normal, charthick=4, size=1.33, color=100

	sfri= sfr1+sfr2
	lastn= n_elements(sfri)

	;-------------------------------------------

	x= [xmin, timem, xmax, reverse(time1)]
	y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
	polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
	xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 0, thick= 6.0
	xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100




	;-------------------------------------------
	;
	;  Labels
	;
	;-------------------------------------------

	; cooling
	; ----------
	;xyouts, 0.7, 0.86, "Cooling", /normal, charthick=1, size=1.33, color=0
	;oplot, [1.65,2.2],[13.0,13.0], psym=-3, color= 0, linestyle= 0
	;xyouts, 0.7, 0.80, "Explicit", /normal, charthick=1, size=1.33, color=50
	;xyouts, 0.7, 0.75, "mzero.cie", /normal, charthick=1, size=1.33, color=100
	;xyouts, 0.7, 0.70, "m-00.cie", /normal, charthick=1, size=1.33, color=150








;--------------------------------------
;--------------------------------------

device, /close


end






;======================================================================
;      PLOT FORMAT
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
;      |----------------------|
;      |                      |
;      |                      |
;      |                      |
;      |      Gas Mass        |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;
;
;
;
;======================================================================


pro msfr3_2, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr3_2, junk"
   print, "  "
   print, "  "
   return
endif


filename='minorsfr3_2.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=10, newysize=16
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------

frunmerg= "G3G2-u3"
fruniso1= "G3il-u1a"
fruniso2= "G2im-u1a"

;frunmerg= "G3G1-u3"
;fruniso1= "G3il-u1a"
;fruniso2= "G1i-u1a"

;frunmerg= "G3G0e-u3"
;fruniso1= "G3il-u1a"
;fruniso2= "G0i-u1a"


;--------------------------------------
;--------------------------------------

;xaxistitle = "Time (Gyr h!E-1!D)"
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
ymax = 30.0
;ymax = 20.0
;ymax = 13.0
;ymax = 8.0
;ymax = 4.0
ymin = 0

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        oplot, [0.2, 0.6], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        xyouts, 0.28, 0.93, "G3G2, 1:2.3 merger", /normal, charthick=4, size=1.33, color=150
        ;xyouts, 0.28, 0.93, "G3G1, 1:5.8 merger", /normal, charthick=4, size=1.33, color=150
        ;xyouts, 0.28, 0.93, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        oplot, [0.2, 0.6], [6.4, 6.4], psym=-3, color= 100, linestyle= 1, thick= 6.0
        xyouts, 0.28, 0.89, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.28, 0.89, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.28, 0.89, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100





;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
yaxistitle="!6Gas Mass / Initial Gas Mass "
ymax = 1.1
ymin = 0

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



;--------------------------------------
;--------------------------------------




; done
; ----------
device, /close

end









;======================================================================
;      PLOT FORMAT
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
;      |----------------------|
;      |                      |
;      |                      |
;      |                      |
;      |      Gas Mass        |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;
;
;
;   currently fig 10  and fig 20
;
;
;   * very similar to msfr3_2, only this compares other
;     things.  note that this used to be in the file
;     minor_sfr_comp.pro, but was eventually combined.
;
;======================================================================


pro msfr_f10_f20, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "msfr_f10_f20, junk"
   print, "  "
   print, "  "
   return
endif


filename='minor_f10_f20.eps'

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
ymax = 20.0
;ymax = 13.0
;ymax = 4.0
ymin = 0

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle


        ;-------------------------------------------
	;
	;  f10a
	;
        ;process_one_sfr, "G3G2-u3", 150, '!8n2med!6', 17.6, 0.93, lstyle= 0
        ;process_one_sfr, "G3G2-u4", 50, '!8n0med!6', 15.7, 0.89, lstyle= 1
        ;xyouts, 0.28, 0.17, "G3G2, 1:2.3 merger", /normal, charthick=4, size=1.33, color=0

        ;-------------------------------------------
	;
	;  f10b
	;
        ;process_one_sfr, "G3G1-u3", 150, '!8n2med!6', 7.1, 0.93, lstyle= 0
        ;process_one_sfr, "G3G1-u4", 50, '!8n0med!6', 6.3, 0.89, lstyle= 1
        ;process_one_sfr, "G3G1-u3", 150, '!8n2med!6', 17.6, 0.93, lstyle= 0
        ;process_one_sfr, "G3G1-u4", 50, '!8n0med!6', 15.7, 0.89, lstyle= 1
        ;xyouts, 0.28, 0.17, "G3G1, 1:5.8 merger", /normal, charthick=4, size=1.33, color=0

        ;-------------------------------------------
	;
	;  f10c
	;
        ;process_one_sfr, "G3G0e-u3", 150, '!8n2med!6', 3.5, 0.93, lstyle= 0
        ;process_one_sfr, "G3G0e-u4", 50, '!8n0med!6', 3.15, 0.89, lstyle= 1
        process_one_sfr, "G3G0e-u3", 150, '!8n2med!6', 17.6, 0.93, lstyle= 0
        process_one_sfr, "G3G0e-u4", 50, '!8n0med!6', 15.7, 0.89, lstyle= 1
        xyouts, 0.28, 0.17, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=0

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
	;
	;  f20
	;
        ;;xyouts, 0.29, 'B/D', 
        ;process_one_sfr, "SbIm_o2_n0", 0, 'B/D= 0', 11.8, 0.94, lstyle= 0
        ;process_one_sfr, "Sbbd1Im_o2_n0", 50, 'B/D= 0.2', 11.0, 0.91, lstyle= 1
        ;process_one_sfr, "Sbbd2Im_o2_n0", 150, 'B/D= 0.5', 10.1, 0.88, lstyle= 2
        ;;xyouts, 0.34, 0.14, '8:1 merger', /normal, charthick=4, size=1.33, color=0

        ;;process_one_sfr, "SbG1_o2_n0", 50, 'B/D= 0', 11.0, 0.91, lstyle= 1


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
	;
	;  f10a
	;
        ;process_one_gasmass, "G3G2-u3", 150, lstyle=0
        ;process_one_gasmass, "G3G2-u4", 50, lstyle=1

        ; burst efficiency
        ;be= 0.32
        ;belbl= strcompress(string(be),/remove_all)
        ;belbl= '!8e!6='+strmid(belbl,0,4)
        ;xyouts, 0.5, 0.50, belbl, /data, charthick=4, size=1.5, color=150
        ;be= 0.38
        ;belbl= strcompress(string(be),/remove_all)
        ;belbl= '!8e!6='+strmid(belbl,0,4)
        ;xyouts, 0.5, 0.35, belbl, /data, charthick=4, size=1.5, color=50


        ;-------------------------------------------
	;
	;  f10b
	;
        ;process_one_gasmass, "G3G1-u3", 150, lstyle=0
        ;process_one_gasmass, "G3G1-u4", 50, lstyle=1

        ; burst efficiency
        ;be= 0.09
        ;belbl= strcompress(string(be),/remove_all)
        ;belbl= '!8e!6='+strmid(belbl,0,4)
        ;xyouts, 0.5, 0.50, belbl, /data, charthick=4, size=1.5, color=150
        ;be= 0.11
        ;belbl= strcompress(string(be),/remove_all)
        ;belbl= '!8e!6='+strmid(belbl,0,4)
        ;xyouts, 0.5, 0.35, belbl, /data, charthick=4, size=1.5, color=50

        ;-------------------------------------------
	;
	;  f10c
	;
        process_one_gasmass, "G3G0e-u3", 150, lstyle=0
        process_one_gasmass, "G3G0e-u4", 50, lstyle=1

        ; burst efficiency
        be= 0.02
        belbl= strcompress(string(be),/remove_all)
        belbl= '!8e!6='+strmid(belbl,0,4)
        xyouts, 0.5, 0.50, belbl, /data, charthick=4, size=1.5, color=150
        be= 0.01
        belbl= strcompress(string(be),/remove_all)
        belbl= '!8e!6='+strmid(belbl,0,4)
        xyouts, 0.5, 0.35, belbl, /data, charthick=4, size=1.5, color=50

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
	;
	;   f20
	;
        ;process_one_gasmass, "SbIm_o2_n0", 0, lstyle=0
        ;process_one_gasmass, "Sbbd1Im_o2_n0", 50, lstyle=1
        ;process_one_gasmass, "Sbbd2Im_o2_n0", 150, lstyle=2

        ;process_one_gasmass, "SbG1_o2_n0", 50, lstyle=1

	; burst efficiency
        ;be= 0.19
        ;belbl= strcompress(string(be),/remove_all)
        ;belbl= '!8e!6='+strmid(belbl,0,4)
        ;xyouts, 0.5, 0.35, belbl, /data, charthick=4, size=1.5, color=0
        ;be= 0.13
        ;belbl= strcompress(string(be),/remove_all)
        ;belbl= '!8e!6='+strmid(belbl,0,4)
        ;xyouts, 0.5, 0.25, belbl, /data, charthick=4, size=1.5, color=50
        ;be= 0.08
        ;belbl= strcompress(string(be),/remove_all)
        ;belbl= '!8e!6='+strmid(belbl,0,4)
        ;xyouts, 0.5, 0.15, belbl, /data, charthick=4, size=1.5, color=150


;--------------------------------------
;--------------------------------------




; done
; ----------
device, /close

end
















;======================================================================
;      PLOT FORMAT
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
;      |----------------------|
;      |                      |
;      |                      |
;      |                      |
;      |      SFR / iso       |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;      |                      |
;      |                      |
;      |                      |
;      |      SFR / mass      |
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
;      |----------------------|
;
;
;
;
;======================================================================


pro msfr3_3, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr3_3, junk"
   print, "  "
   print, "  "
   return
endif


filename='minorsfr3_3.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=10, newysize=22
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------

;frunmerg= "G3G2-u3"      & mmerg= 6.5
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

;frunmerg= "G3G1-u3"      & mmerg= 5.5
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G1i-u1a"      & miso2= 0.5

frunmerg= "G3G0e-u3"     & mmerg= 5.1
fruniso1= "G3il-u1a"     & miso1= 5.0
fruniso2= "G0i-u1a"      & miso2= 0.1


;--------------------------------------
;--------------------------------------

;xaxistitle = "Time (Gyr h!E-1!D)"
xaxistitle = "!6Time (Gyr)"
xmax = 6.0
xmin = 0


x0= 0.18
x1= 0.97

y0= 0.08
y1= 0.3275
y2= 0.55
y3= 0.7725
y4= 0.99

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
ymax = 8.0
ymin = 0

!p.position= [x0, y3, x1, y4]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.29, 0.96, "G3G2, 1:2.3 merger", /normal, charthick=4, size=1.33, color=150
        ;xyouts, 0.29, 0.96, "G3G1, 1:5.8 merger", /normal, charthick=4, size=1.33, color=150
        xyouts, 0.29, 0.96, "G3G0, 1:22.7 merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; Top - Middle  Plot
;--------------------------------------

yaxistitle="!6SFR / Isolated "
ymax = 13.4
;ymax = 9.4
ymin = 0.0

!p.position= [x0, y2, x1, y3]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle, /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

	sfrmsfri= sfrm / sfri

	oplot, time1, sfrmsfri, psym=-3, color= 0, linestyle= 0, thick= 9.0

	; draw line at 1
	oplot, [xmin,xmax], [1.0,1.0], psym=-3, color= 0, linestyle= 1, thick= 2.0


;--------------------------------------
; Middle-Bottom Plot
;--------------------------------------

;yaxistitle="!6SFR / Stellar Mass (0.1 Gyr!E-1!N)"
;yaxistitle="!6SFR / M!D*!N (0.1 Gyr!E-1!N)"
yaxistitle="!6SFR / M!D*!N (Gyr!E-1!N)"
ymax = 0.18
ymin = 0.0006

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ;xtickformat='(a1)', ytitle=yaxistitle, /noerase
        xtickformat='(a1)', ytitle=yaxistitle, /noerase, ytickformat='exp_label', /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
	sfrmfs= sfrmfs + mmerg
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrmfs)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

	sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
	mfsi= mfs1+mfs2 + miso1 + miso2
	mfssfri= sfri/mfsi
	mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100






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
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



;--------------------------------------
;--------------------------------------




; done
; ----------
device, /close

end














;======================================================================
;      PLOT FORMAT
;
;
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      |    SFR      |             |             |             |
;      |             |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      |  SFR / iso  |             |             |             |
;      |             |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      | SFR / mass  |             |             |             |
;      |             |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      |   Gas Mass  |             |             |             |
;      |             |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;
;
;
;
;
;
;
;
;======================================================================


pro msfr4, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr3_3, junk"
   print, "  "
   print, "  "
   return
endif


filename='minorsfr4.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=28, newysize=22
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------

;xaxistitle = "Time (Gyr h!E-1!D)"
xaxistitle = "!6Time (Gyr)"
xmax = 6.0
xmin = 0

x0= 0.08
x1= 0.3275
x2= 0.55
x3= 0.7725
x4= 0.99

y0= 0.08
y1= 0.3275
y2= 0.55
y3= 0.7725
y4= 0.99


;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  111111111111111111
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G3b-u1"      & mmerg= 10.0     &   ymax = 12.0
;fruniso1= "G3il-u1a"      & miso1= 5.0
;fruniso2= "G3il-u1a"      & miso2= 5.0   

frunmerg= "G3G3b-u2"      & mmerg= 10.0     &   ymax = 28.0
fruniso1= "G3_n0"         & miso1= 5.0
fruniso2= "G3_n0"         & miso2= 5.0   

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x0, y3, x1, y4]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        xyouts, x0+0.12, 0.96, "G3G3, 1:1", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; Top-Middle Plot
;--------------------------------------

yaxistitle="!6SFR / Isolated "
ymax = 90.0
;ymax = 13.4
;ymax = 9.4
ymin = 0.2
;ymin = 0.0

!p.position= [x0, y2, x1, y3]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle, /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

	sfrmsfri= sfrm / sfri

	oplot, time1, sfrmsfri, psym=-3, color= 0, linestyle= 0, thick= 9.0

	; draw line at 1
	oplot, [xmin,xmax], [1.0,1.0], psym=-3, color= 0, linestyle= 1, thick= 2.0


;--------------------------------------
; Middle-Bottom Plot
;--------------------------------------

;yaxistitle="!6SFR / Stellar Mass (0.1 Gyr!E-1!N)"
;yaxistitle="!6SFR / M!D*!N (0.1 Gyr!E-1!N)"
yaxistitle="!6SFR / M!D*!N (Gyr!E-1!N)"
ymax = 0.28
;ymax = 0.18
ymin = 0.0006

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ;xtickformat='(a1)', ytitle=yaxistitle, /noerase
        xtickformat='(a1)', ytitle=yaxistitle, /noerase, ytickformat='exp_label', /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
	sfrmfs= sfrmfs + mmerg
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrmfs)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

	sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
	mfsi= mfs1+mfs2 + miso1 + miso2
	mfssfri= sfri/mfsi
	mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100






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
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  22222222222222
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G2-u3"      & mmerg= 6.5    &  ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

frunmerg= "G3G2-u4"      & mmerg= 6.5    &  ymax = 28.0
fruniso1= "G3_n0"        & miso1= 5.0
fruniso2= "G2im-u2"      & miso2= 1.5   

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x1, y3, x2, y4]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        xyouts, x1+0.12, 0.96, "G3G2, 1:2.3", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; Top-Middle Plot
;--------------------------------------

yaxistitle="!6SFR / Isolated "
ymax = 90.0
;ymax = 13.4
;ymax = 9.4
ymin = 0.2
;ymin = 0.0

!p.position= [x1, y2, x2, y3]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

	sfrmsfri= sfrm / sfri

	oplot, time1, sfrmsfri, psym=-3, color= 0, linestyle= 0, thick= 9.0

	; draw line at 1
	oplot, [xmin,xmax], [1.0,1.0], psym=-3, color= 0, linestyle= 1, thick= 2.0


;--------------------------------------
; Middle-Bottom Plot
;--------------------------------------

;yaxistitle="!6SFR / Stellar Mass (0.1 Gyr!E-1!N)"
;yaxistitle="!6SFR / M!D*!N (0.1 Gyr!E-1!N)"
yaxistitle="!6SFR / M!D*!N (Gyr!E-1!N)"
ymax = 0.28
;ymax = 0.18
ymin = 0.0006

!p.position= [x1, y1, x2, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ;xtickformat='(a1)', ytickformat='(a1)', /noerase
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
	sfrmfs= sfrmfs + mmerg
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrmfs)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

	sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
	mfsi= mfs1+mfs2 + miso1 + miso2
	mfssfri= sfri/mfsi
	mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100






;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x1, y0, x2, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  33333333333333
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G1-u3"      & mmerg= 5.5     &  ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G1i-u1a"      & miso2= 0.5

frunmerg= "G3G1-u4"      & mmerg= 5.5     &  ymax = 28.0
fruniso1= "G3_n0"        & miso1= 5.0
fruniso2= "G1i-u2"       & miso2= 0.5

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x2, y3, x3, y4]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        xyouts, x2+0.11, 0.96, "G3G1, 1:5.8", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; Top-Middle Plot
;--------------------------------------

yaxistitle="!6SFR / Isolated "
ymax = 90.0
;ymax = 13.4
;ymax = 9.4
ymin = 0.2
;ymin = 0.0

!p.position= [x2, y2, x3, y3]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

	sfrmsfri= sfrm / sfri

	oplot, time1, sfrmsfri, psym=-3, color= 0, linestyle= 0, thick= 9.0

	; draw line at 1
	oplot, [xmin,xmax], [1.0,1.0], psym=-3, color= 0, linestyle= 1, thick= 2.0


;--------------------------------------
; Middle-Bottom Plot
;--------------------------------------

;yaxistitle="!6SFR / Stellar Mass (0.1 Gyr!E-1!N)"
;yaxistitle="!6SFR / M!D*!N (0.1 Gyr!E-1!N)"
yaxistitle="!6SFR / M!D*!N (Gyr!E-1!N)"
ymax = 0.28
;ymax = 0.18
ymin = 0.0006

!p.position= [x2, y1, x3, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ;xtickformat='(a1)', ytickformat='(a1)', /noerase
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
	sfrmfs= sfrmfs + mmerg
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrmfs)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

	sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
	mfsi= mfs1+mfs2 + miso1 + miso2
	mfssfri= sfri/mfsi
	mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100






;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x2, y0, x3, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  4444444444444444
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G0e-u3"     & mmerg= 5.1    & ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G0i-u1a"      & miso2= 0.1

frunmerg= "G3G0e-u4"     & mmerg= 5.1    & ymax = 28.0
fruniso1= "G3_n0"        & miso1= 5.0
fruniso2= "G0i-u2"       & miso2= 0.1

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x3, y3, x4, y4]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        xyouts, x3+0.10, 0.96, "G3G0, 1:22.7", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; Top-Middle Plot
;--------------------------------------

yaxistitle="!6SFR / Isolated "
ymax = 90.0
;ymax = 9.4
ymin = 0.2
;ymin = 0.0

!p.position= [x3, y2, x4, y3]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

	sfrmsfri= sfrm / sfri

	oplot, time1, sfrmsfri, psym=-3, color= 0, linestyle= 0, thick= 9.0

	; draw line at 1
	oplot, [xmin,xmax], [1.0,1.0], psym=-3, color= 0, linestyle= 1, thick= 2.0


;--------------------------------------
; Middle-Bottom Plot
;--------------------------------------

;yaxistitle="!6SFR / Stellar Mass (0.1 Gyr!E-1!N)"
;yaxistitle="!6SFR / M!D*!N (0.1 Gyr!E-1!N)"
yaxistitle="!6SFR / M!D*!N (Gyr!E-1!N)"
ymax = 0.28
;ymax = 0.18
ymin = 0.0006

!p.position= [x3, y1, x4, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ;xtickformat='(a1)', ytickformat='(a1)', /noerase
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
	sfrmfs= sfrmfs + mmerg
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrmfs)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

	sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
	mfsi= mfs1+mfs2 + miso1 + miso2
	mfssfri= sfri/mfsi
	mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100






;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x3, y0, x4, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



;--------------------------------------
;--------------------------------------



;--------------------------------------
;--------------------------------------
; done
; ----------
device, /close

end








;======================================================================
;      PLOT FORMAT
;
;
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      |    SFR      |             |             |             |
;      |             |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      |  SFR / iso  |             |             |             |
;      |             |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      | SFR / st.   |             |             |             |
;      |      mass   |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      | SFR / gas   |             |             |             |
;      |      mass   |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      |   Gas Mass  |             |             |             |
;      |             |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;
;
;
;
;  this is now figure 14 in minors paper
;
;
;
;======================================================================


pro msfr_f14, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr_f14, junk"
   print, "  "
   print, "  "
   return
endif


filename='minor_f14.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=28, newysize=27
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------

;xaxistitle = "Time (Gyr h!E-1!D)"
xaxistitle = "!6Time (Gyr)"
xmax = 6.0
xmin = 0

x0= 0.08
x1= 0.3275
x2= 0.55
x3= 0.7725
x4= 0.99

y0= 0.08
y1= 0.262
y2= 0.444
y3= 0.626
y4= 0.808
y5= 0.99


;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  111111111111111111
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G3b-u1"      & mmerg= 10.0     &   ymax = 12.0
;fruniso1= "G3il-u1a"      & miso1= 5.0
;fruniso2= "G3il-u1a"      & miso2= 5.0   

frunmerg= "G3G3b-u2"      & mmerg= 10.0     &   ymax = 28.0
fruniso1= "G3_n0"         & miso1= 5.0
fruniso2= "G3_n0"         & miso2= 5.0   

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x0, y4, x1, y5]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        xyouts, x0+0.12, 0.96, "G3G3, 1:1", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; Top-Middle Plot
;--------------------------------------

yaxistitle="!6SFR / Isolated "
ymax = 90.0
;ymax = 13.4
;ymax = 9.4
ymin = 0.2
;ymin = 0.0

!p.position= [x0, y3, x1, y4]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle, /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

	sfrmsfri= sfrm / sfri

	oplot, time1, sfrmsfri, psym=-3, color= 0, linestyle= 0, thick= 9.0

	; draw line at 1
	oplot, [xmin,xmax], [1.0,1.0], psym=-3, color= 0, linestyle= 1, thick= 2.0




;--------------------------------------
; Middle-Bottom Plot
;--------------------------------------

;yaxistitle="!6SFR / Stellar Mass (0.1 Gyr!E-1!N)"
;yaxistitle="!6SFR / M!D*!N (0.1 Gyr!E-1!N)"
yaxistitle="!6SFR / M!D*!N (Gyr!E-1!N)"
ymax = 2.8
;ymax = 0.18
ymin = 0.0006

!p.position= [x0, y2, x1, y3]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ;xtickformat='(a1)', ytitle=yaxistitle, /noerase
        xtickformat='(a1)', ytitle=yaxistitle, /noerase, ytickformat='exp_label', /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
	sfrmfs= sfrmfs + mmerg
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrmfs)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

	sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
	mfsi= mfs1+mfs2 + miso1 + miso2
	mfssfri= sfri/mfsi
	mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; New Plot
;--------------------------------------

yaxistitle="!6SFR / M!Dgas!N (Gyr!E-1!N)"
ymax = 2.8
ymin = 0.0006

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle, /noerase, ytickformat='exp_label', /ylog


        ;-------------------------------------------
        ;
        ;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

        sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr) 
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr) 
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
        mfsi= mfs1+mfs2
        mfssfri= sfri/mfsi
        mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)  

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100








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
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
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
        xyouts, 0.5, 0.21, belbl, /data, charthick=4, size=1.5, color=0



;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  22222222222222
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G2-u3"      & mmerg= 6.5    &  ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

frunmerg= "G3G2-u4"      & mmerg= 6.5    &  ymax = 28.0
fruniso1= "G3_n0"        & miso1= 5.0
fruniso2= "G2im-u2"      & miso2= 1.5   

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x1, y4, x2, y5]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        xyouts, x1+0.12, 0.96, "G3G2, 1:2.3", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; Top-Middle Plot
;--------------------------------------

yaxistitle="!6SFR / Isolated "
ymax = 90.0
;ymax = 13.4
;ymax = 9.4
ymin = 0.2
;ymin = 0.0

!p.position= [x1, y3, x2, y4]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

	sfrmsfri= sfrm / sfri

	oplot, time1, sfrmsfri, psym=-3, color= 0, linestyle= 0, thick= 9.0

	; draw line at 1
	oplot, [xmin,xmax], [1.0,1.0], psym=-3, color= 0, linestyle= 1, thick= 2.0


;--------------------------------------
; Middle-Bottom Plot
;--------------------------------------

;yaxistitle="!6SFR / Stellar Mass (0.1 Gyr!E-1!N)"
;yaxistitle="!6SFR / M!D*!N (0.1 Gyr!E-1!N)"
yaxistitle="!6SFR / M!D*!N (Gyr!E-1!N)"
ymax = 2.8
;ymax = 0.18
ymin = 0.0006

!p.position= [x1, y2, x2, y3]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ;xtickformat='(a1)', ytickformat='(a1)', /noerase
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
	sfrmfs= sfrmfs + mmerg
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrmfs)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

	sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
	mfsi= mfs1+mfs2 + miso1 + miso2
	mfssfri= sfri/mfsi
	mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; New Plot
;--------------------------------------

yaxistitle="!6SFR / M!Dgas!N (Gyr!E-1!N)"
ymax = 2.8
ymin = 0.0006

!p.position= [x1, y1, x2, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog 


        ;-------------------------------------------
        ;
        ;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        sfrmfs= sfrmfs + mmerg
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

        sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr) 
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
        mfsi= mfs1+mfs2
        mfssfri= sfri/mfsi
        mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)  

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100






;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x1, y0, x2, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
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
        xyouts, 0.5, 0.21, belbl, /data, charthick=4, size=1.5, color=0



;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  33333333333333
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G1-u3"      & mmerg= 5.5     &  ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G1i-u1a"      & miso2= 0.5

frunmerg= "G3G1-u4"      & mmerg= 5.5     &  ymax = 28.0
fruniso1= "G3_n0"        & miso1= 5.0
fruniso2= "G1i-u2"       & miso2= 0.5

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x2, y4, x3, y5]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        xyouts, x2+0.11, 0.96, "G3G1, 1:5.8", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; Top-Middle Plot
;--------------------------------------

yaxistitle="!6SFR / Isolated "
ymax = 90.0
;ymax = 13.4
;ymax = 9.4
ymin = 0.2
;ymin = 0.0

!p.position= [x2, y3, x3, y4]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

	sfrmsfri= sfrm / sfri

	oplot, time1, sfrmsfri, psym=-3, color= 0, linestyle= 0, thick= 9.0

	; draw line at 1
	oplot, [xmin,xmax], [1.0,1.0], psym=-3, color= 0, linestyle= 1, thick= 2.0


;--------------------------------------
; Middle-Bottom Plot
;--------------------------------------

;yaxistitle="!6SFR / Stellar Mass (0.1 Gyr!E-1!N)"
;yaxistitle="!6SFR / M!D*!N (0.1 Gyr!E-1!N)"
yaxistitle="!6SFR / M!D*!N (Gyr!E-1!N)"
ymax = 2.8
;ymax = 0.18
ymin = 0.0006

!p.position= [x2, y2, x3, y3]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ;xtickformat='(a1)', ytickformat='(a1)', /noerase
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
	sfrmfs= sfrmfs + mmerg
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrmfs)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

	sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
	mfsi= mfs1+mfs2 + miso1 + miso2
	mfssfri= sfri/mfsi
	mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



;--------------------------------------
; New Plot
;--------------------------------------

yaxistitle="!6SFR / M!Dgas!N (Gyr!E-1!N)"
ymax = 2.8
ymin = 0.0006

!p.position= [x2, y1, x3, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog 


        ;-------------------------------------------
        ;
        ;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        sfrmfs= sfrmfs + mmerg
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

        sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr) 
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
        mfsi= mfs1+mfs2
        mfssfri= sfri/mfsi
        mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)  

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100







;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x2, y0, x3, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
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
        xyouts, 0.5, 0.21, belbl, /data, charthick=4, size=1.5, color=0



;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  4444444444444444
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G0e-u3"     & mmerg= 5.1    & ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G0i-u1a"      & miso2= 0.1

frunmerg= "G3G0e-u4"     & mmerg= 5.1    & ymax = 28.0
fruniso1= "G3_n0"        & miso1= 5.0
fruniso2= "G0i-u2"       & miso2= 0.1

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x3, y4, x4, y5]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        xyouts, x3+0.10, 0.96, "G3G0, 1:22.7", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; Top-Middle Plot
;--------------------------------------

yaxistitle="!6SFR / Isolated "
ymax = 90.0
;ymax = 9.4
ymin = 0.2
;ymin = 0.0

!p.position= [x3, y3, x4, y4]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

	sfrmsfri= sfrm / sfri

	oplot, time1, sfrmsfri, psym=-3, color= 0, linestyle= 0, thick= 9.0

	; draw line at 1
	oplot, [xmin,xmax], [1.0,1.0], psym=-3, color= 0, linestyle= 1, thick= 2.0


;--------------------------------------
; Middle-Bottom Plot
;--------------------------------------

;yaxistitle="!6SFR / Stellar Mass (0.1 Gyr!E-1!N)"
;yaxistitle="!6SFR / M!D*!N (0.1 Gyr!E-1!N)"
yaxistitle="!6SFR / M!D*!N (Gyr!E-1!N)"
ymax = 2.8
;ymax = 0.18
ymin = 0.0006

!p.position= [x3, y2, x4, y3]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        ;xtickformat='(a1)', ytickformat='(a1)', /noerase
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
	sfrmfs= sfrmfs + mmerg
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrmfs)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

	sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
	mfsi= mfs1+mfs2 + miso1 + miso2
	mfssfri= sfri/mfsi
	mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100


;--------------------------------------
; New Plot
;--------------------------------------

yaxistitle="!6SFR / M!Dgas!N (Gyr!E-1!N)"
ymax = 2.8
ymin = 0.0006

!p.position= [x3, y1, x4, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.30, ycharsize=1.30, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase, /ylog 


        ;-------------------------------------------
        ;
        ;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        sfrmfs= sfrmfs + mmerg
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr/sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, sfrmfsm, timem, Nelements=100

        sfrmfsm= sfrmfsm / 10.0   ; convert to  Gyr^-1

        ;-------------------------------------------
        ;
        ;  2 x Isolated
        ;
        ;-------------------------------------------
        minor_open_sfr_file, fruniso1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfr1, time1, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, mfs1, time1, Nelements=100

        minor_open_sfr_file, fruniso2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr) 
        resample_sfrhist, oldsfr, oldtime, sfr2, time2, Nelements=100
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        resample_sfrhist, oldsfr, oldtime, mfs2, time2, Nelements=100

        sfri= sfr1+sfr2
        mfsi= mfs1+mfs2
        mfssfri= sfri/mfsi
        mfssfri= mfssfri / 10.0  ; convert to Gyr^-1
        lastn= n_elements(sfri)  

        ;-------------------------------------------

        ;x= [xmin, timem, xmax, reverse(time1)]
        ;y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        ;polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrmfsm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, mfssfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100







;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x3, y0, x4, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
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
        xyouts, 0.5, 0.21, belbl, /data, charthick=4, size=1.5, color=0



;--------------------------------------
;--------------------------------------



;--------------------------------------
;--------------------------------------
; done
; ----------
device, /close

end













;======================================================================
;      PLOT FORMAT
;
;
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      |    SFR      |             |             |             |
;      |             |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      |   Gas Mass  |             |             |             |
;      |             |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;
;
;
;
;
;
;======================================================================


pro msfr5, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr5, junk"
   print, "  "
   print, "  "
   return
endif


filename='minorsfr5.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=28, newysize=14
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------

;xaxistitle = "Time (Gyr h!E-1!D)"
xaxistitle = "!6Time (Gyr)"
xmax = 6.0
xmin = 0

x0= 0.08
x1= 0.3275
x2= 0.55
x3= 0.7725
x4= 0.99

y0= 0.12
y2= 0.99
ys= (y2-y0)/2.
y1= y0+ys


;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  111111111111111111
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G3b-u1"      & mmerg= 10.0     &   ymax = 12.0
;fruniso1= "G3il-u1a"      & miso1= 5.0
;fruniso2= "G3il-u1a"      & miso2= 5.0   

;frunmerg= "G3G3b-u2"      & mmerg= 10.0     &   ymax = 28.0
;fruniso1= "G3_n0"         & miso1= 5.0
;fruniso2= "G3_n0"         & miso2= 5.0   

frunmerg= "G3Rd4eG2_n0"      & mmerg= 5.5     &  ymax = 32.0
fruniso1= "G3Rd4e_n0"        & miso1= 5.0
fruniso2= "G2im-u2"          & miso2= 0.5


;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


print, "-------------------------------------------------------------"
print, "WARNING: manually fixing this so the range looks good."
print, "it doesn't change the results that we present, and also saves"
print, "us from having to change the exes."
print, "-------------------------------------------------------------"
idx=where(sfrm gt ymax)
if idx(0) ne -1 then sfrm(idx)= ymax

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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x0+0.12, 0.96, "G3G3, 1:1", /normal, charthick=4, size=1.33, color=0
        xyouts, x0+0.12, 0.96, "G3Rd4eG2_n0", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100





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
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100


        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= strmid(belbl,0,4)
        xyouts, 4.0, 0.21, belbl, /data, charthick=4, size=1.5, color=0




;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  22222222222222
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G2-u3"      & mmerg= 6.5    &  ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

;frunmerg= "G3G2-u4"      & mmerg= 6.5    &  ymax = 28.0
;fruniso1= "G3_n0"        & miso1= 5.0
;fruniso2= "G2im-u2"      & miso2= 1.5   

frunmerg= "G3Rd4hG2_n0"      & mmerg= 6.5    &  ymax = 32.0
fruniso1= "G3Rd4h_n0"        & miso1= 5.0
fruniso2= "G2im-u2"          & miso2= 1.5   

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x1, y1, x2, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x1+0.12, 0.96, "G3G2, 1:2.3", /normal, charthick=4, size=1.33, color=0
        xyouts, x1+0.12, 0.96, "G3Rd4hG2_n0", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100






;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x1, y0, x2, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100


        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= strmid(belbl,0,4)
        xyouts, 4.0, 0.21, belbl, /data, charthick=4, size=1.5, color=0



;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  33333333333333
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G1-u3"      & mmerg= 5.5     &  ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G1i-u1a"      & miso2= 0.5

;frunmerg= "G3G1-u4"      & mmerg= 5.5     &  ymax = 28.0
;fruniso1= "G3_n0"        & miso1= 5.0
;fruniso2= "G1i-u2"       & miso2= 0.5

frunmerg= "G3Rd4G2_n0"      & mmerg= 6.5    &  ymax = 32.0
fruniso1= "G3Rd4"           & miso1= 5.0
fruniso2= "G2im-u2"         & miso2= 1.5   

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x2, y1, x3, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x2+0.11, 0.96, "G3G1, 1:5.8", /normal, charthick=4, size=1.33, color=0
        xyouts, x2+0.11, 0.96, "G3Rd4G2_n0", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100








;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x2, y0, x3, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= strmid(belbl,0,4)
        xyouts, 4.0, 0.21, belbl, /data, charthick=4, size=1.5, color=0



;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  4444444444444444
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G0e-u3"     & mmerg= 5.1    & ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G0i-u1a"      & miso2= 0.1

;frunmerg= "G3G0e-u4"     & mmerg= 5.1    & ymax = 28.0
;fruniso1= "G3_n0"        & miso1= 5.0
;fruniso2= "G0i-u2"       & miso2= 0.1

frunmerg= "G3Rd4iG2_n0"      & mmerg= 6.5    &  ymax = 32.0
fruniso1= "G3Rd4i_n0"        & miso1= 5.0
fruniso2= "G2im-u2"          & miso2= 1.5   

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x3, y1, x4, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x3+0.10, 0.96, "G3G0, 1:22.7", /normal, charthick=4, size=1.33, color=0
        xyouts, x3+0.10, 0.96, "G3Rd4iG2_n0", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100





;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x3, y0, x4, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= strmid(belbl,0,4)
        xyouts, 4.0, 0.21, belbl, /data, charthick=4, size=1.5, color=0


;--------------------------------------
;--------------------------------------



;--------------------------------------
;--------------------------------------
; done
; ----------
device, /close

end












;======================================================================
;      PLOT FORMAT
;
;
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      |    SFR      |             |             |             |
;      |             |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;      |             |             |             |             |
;      |             |             |             |             |
;      |   Gas Mass  |             |             |             |
;      |             |             |             |             |
;      |             |             |             |             |
;      |-------------|-------------|-------------|-------------|
;
;
;
;
;
;   NOT currently used - couldn't figure out a good way to 
;                        deal with the varying SFR y-axis
;
;
;
;======================================================================


pro msfr_nf10, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr_nf10, junk"
   print, "  "
   print, "  "
   return
endif


filename='minor_nf10.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=28, newysize=14
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------

;xaxistitle = "Time (Gyr h!E-1!D)"
xaxistitle = "!6Time (Gyr)"
xmax = 6.0
xmin = 0

x0= 0.08
x1= 0.3275
x2= 0.55
x3= 0.7725
x4= 0.99

y0= 0.12
y2= 0.99
ys= (y2-y0)/2.
y1= y0+ys


;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  111111111111111111
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G3b-u1"      & mmerg= 10.0     &   ymax = 12.0
;fruniso1= "G3il-u1a"      & miso1= 5.0
;fruniso2= "G3il-u1a"      & miso2= 5.0   

;frunmerg= "G3G3b-u2"      & mmerg= 10.0     &   ymax = 28.0
;fruniso1= "G3_n0"         & miso1= 5.0
;fruniso2= "G3_n0"         & miso2= 5.0   

frunmerg= "G3Rd4eG2_n0"      & mmerg= 5.5     &  ymax = 32.0
fruniso1= "G3Rd4e_n0"        & miso1= 5.0
fruniso2= "G2im-u2"          & miso2= 0.5


;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle

        ;-------------------------------------------
        process_one_sfr, "G3G2-u3", 150, '!8n2med!6', 17.6, 0.93, lstyle= 0
        process_one_sfr, "G3G2-u4", 50, '!8n0med!6', 15.7, 0.89, lstyle= 1
        xyouts, 0.28, 0.17, "G3G2, 1:2.3 merger", /normal, charthick=4, size=1.33, color=0





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
        process_one_gasmass, "G3G2-u3", 150, lstyle=0, be=be
        process_one_gasmass, "G3G2-u4", 50, lstyle=1, be=be

        ;-------------------------------------------
        ;process_one_gasmass, "G3G1-u3", 150, lstyle=0, be=be
        ;process_one_gasmass, "G3G1-u4", 50, lstyle=1, be=be

        ;-------------------------------------------
        ;process_one_gasmass, "G3G0e-u3", 150, lstyle=0
        ;process_one_gasmass, "G3G0e-u4", 50, lstyle=1



        ; calculate burst efficiency
	; will need to manually do this
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= strmid(belbl,0,4)
        xyouts, 4.0, 0.21, belbl, /data, charthick=4, size=1.5, color=0




;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  22222222222222
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G2-u3"      & mmerg= 6.5    &  ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

;frunmerg= "G3G2-u4"      & mmerg= 6.5    &  ymax = 28.0
;fruniso1= "G3_n0"        & miso1= 5.0
;fruniso2= "G2im-u2"      & miso2= 1.5   

frunmerg= "G3Rd4hG2_n0"      & mmerg= 6.5    &  ymax = 32.0
fruniso1= "G3Rd4h_n0"        & miso1= 5.0
fruniso2= "G2im-u2"          & miso2= 1.5   

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x1, y1, x2, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x1+0.12, 0.96, "G3G2, 1:2.3", /normal, charthick=4, size=1.33, color=0
        xyouts, x1+0.12, 0.96, "G3Rd4hG2_n0", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100






;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x1, y0, x2, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100


        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= strmid(belbl,0,4)
        xyouts, 4.0, 0.21, belbl, /data, charthick=4, size=1.5, color=0



;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  33333333333333
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G1-u3"      & mmerg= 5.5     &  ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G1i-u1a"      & miso2= 0.5

;frunmerg= "G3G1-u4"      & mmerg= 5.5     &  ymax = 28.0
;fruniso1= "G3_n0"        & miso1= 5.0
;fruniso2= "G1i-u2"       & miso2= 0.5

frunmerg= "G3Rd4G2_n0"      & mmerg= 6.5    &  ymax = 32.0
fruniso1= "G3Rd4"           & miso1= 5.0
fruniso2= "G2im-u2"         & miso2= 1.5   

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x2, y1, x3, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x2+0.11, 0.96, "G3G1, 1:5.8", /normal, charthick=4, size=1.33, color=0
        xyouts, x2+0.11, 0.96, "G3Rd4G2_n0", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100








;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x2, y0, x3, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= strmid(belbl,0,4)
        xyouts, 4.0, 0.21, belbl, /data, charthick=4, size=1.5, color=0



;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  4444444444444444
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G0e-u3"     & mmerg= 5.1    & ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G0i-u1a"      & miso2= 0.1

;frunmerg= "G3G0e-u4"     & mmerg= 5.1    & ymax = 28.0
;fruniso1= "G3_n0"        & miso1= 5.0
;fruniso2= "G0i-u2"       & miso2= 0.1

frunmerg= "G3Rd4iG2_n0"      & mmerg= 6.5    &  ymax = 32.0
fruniso1= "G3Rd4i_n0"        & miso1= 5.0
fruniso2= "G2im-u2"          & miso2= 1.5   

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x3, y1, x4, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x3+0.10, 0.96, "G3G0, 1:22.7", /normal, charthick=4, size=1.33, color=0
        xyouts, x3+0.10, 0.96, "G3Rd4iG2_n0", /normal, charthick=4, size=1.33, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100





;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x3, y0, x4, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



        ; calculate burst efficiency
        be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= strmid(belbl,0,4)
        xyouts, 4.0, 0.21, belbl, /data, charthick=4, size=1.5, color=0


;--------------------------------------
;--------------------------------------



;--------------------------------------
;--------------------------------------
; done
; ----------
device, /close

end

















;======================================================================
;      PLOT FORMAT
;
;
;      |-------------|-------------|-------------|
;      |             |             |             |
;      |             |             |             |
;      |    SFR      |             |             |
;      |             |             |             |
;      |             |             |             |
;      |-------------|-------------|-------------|
;      |             |             |             |
;      |             |             |             |
;      |   Gas Mass  |             |             |
;      |             |             |             |
;      |             |             |             |
;      |-------------|-------------|-------------|
;
;
;
;
;   currently is figure 21 in minors paper
;
;
;======================================================================


pro msfr_f21, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr_f21, junk"
   print, "  "
   print, "  "
   return
endif


filename='minor_f21.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=28, newysize=18
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------

;xaxistitle = "Time (Gyr h!E-1!D)"
xaxistitle = "!6Time (Gyr)"
xmax = 6.0
xmin = 0

x0= 0.10
x3= 0.99
xs= (x3-x0)/3.
x1= x0+xs
x2= x0+xs+xs

y0= 0.12
y1= 0.555
y2= 0.99


;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  111111111111111111
;
;
;--------------------------------------
;--------------------------------------

frunmerg= "G3G2-u3"      & mmerg= 6.5    &  ymax = 20.0
fruniso1= "G3il-u1a"     & miso1= 5.0
fruniso2= "G2im-u1a"     & miso2= 1.5   

;frunmerg= "G3G2-u4"      & mmerg= 6.5    &  ymax = 28.0
;fruniso1= "G3_n0"        & miso1= 5.0
;fruniso2= "G2im-u2"      & miso2= 1.5   

;frunmerg= "G3Rd4eG2"      & mmerg= 5.5     &  ymax = 12.0
;fruniso1= "G3Rd4e"        & miso1= 5.0
;fruniso2= "G2im-u1a"        & miso2= 0.5

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=2.10, ycharsize=2.10, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x0+0.12, 0.96, "G3G2, 2.3:1", /normal, charthick=4, size=1.33, color=0
        ;xyouts, x0+0.12, 0.96, frunmerg+', 2.3:1', /normal, charthick=4, size=1.33, color=0

        xyouts, x0+0.12, 0.90, '!8f!6= 0.20', /normal, charthick=4, size=2.6, color=0
        xyouts, x0+0.12, 0.82, '(fiducial)', /normal, charthick=4, size=2.6, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100





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
	xcharsize=2.10, ycharsize=2.10, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
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
        xyouts, 2.5, 0.21, belbl, /data, charthick=4, size=2.6, color=0




;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  22222222222222
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3Rd4G2"      & mmerg= 6.5    &  ymax = 12.0
;fruniso1= "G3Rd4"        & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

;frunmerg= "G3Rd4hG2"      & mmerg= 6.5    &  ymax = 12.0
;fruniso1= "G3Rd4h"        & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

frunmerg= "G3gf1G2-u1"      & mmerg= 6.5    &  ymax = 20.0
fruniso1= "G3gf1i-u1"        & miso1= 5.0
fruniso2= "G2im-u1a"     & miso2= 1.5   

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x1, y1, x2, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x1+0.12, 0.96, "G3Rd4G2, 1:2.3", /normal, charthick=4, size=1.33, color=0
        ;xyouts, x1+0.12, 0.96, frunmerg+', 2.3:1', /normal, charthick=4, size=1.33, color=0

        xyouts, x1+0.12, 0.90, '!8f!6= 0.50', /normal, charthick=4, size=2.6, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100






;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x1, y0, x2, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=2.10, ycharsize=2.10, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
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
        xyouts, 2.5, 0.21, belbl, /data, charthick=4, size=2.6, color=0



;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  33333333333333
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3Rd4eG2"      & mmerg= 5.5     &  ymax = 28.0
;fruniso1= "G3Rd4e"        & miso1= 5.0
;fruniso2= "G1i-u2"        & miso2= 0.5

;frunmerg= "G3Rd4G2"      & mmerg= 6.5    &  ymax = 12.0
;fruniso1= "G3Rd4"        & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

frunmerg= "G3gf2G2-u1"      & mmerg= 6.5    &  ymax = 20.0
fruniso1= "G3gf2i-u1"        & miso1= 5.0
fruniso2= "G2im-u1a"     & miso2= 1.5   

;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x2, y1, x3, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x2+0.11, 0.96, "G3Rd4eG2, 1:5.8", /normal, charthick=4, size=1.33, color=0
        ;xyouts, x2+0.12, 0.96, frunmerg+', 2.3:1', /normal, charthick=4, size=1.33, color=0

        xyouts, x2+0.12, 0.90, '!8f!6= 0.78', /normal, charthick=4, size=2.6, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x2, y0, x3, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=2.10, ycharsize=2.10, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
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
        xyouts, 2.5, 0.21, belbl, /data, charthick=4, size=2.6, color=0



;--------------------------------------
;--------------------------------------



;--------------------------------------
;--------------------------------------
; done
; ----------
device, /close

end












;======================================================================
;      PLOT FORMAT
;
;
;      |-------------|-------------|
;      |             |             |
;      |             |             |
;      |    SFR      |             |
;      |             |             |
;      |             |             |
;      |-------------|-------------|
;      |             |             |
;      |             |             |
;      |   Gas Mass  |             |
;      |             |             |
;      |             |             |
;      |-------------|-------------|
;
;
;
;   this is fig 18 currently, and now its not
;
;
;======================================================================


pro msfr_f18, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr_f18, junk"
   print, "  "
   print, "  "
   return
endif


filename='minorsfr_f18.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=20, newysize=18
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------

;xaxistitle = "Time (Gyr h!E-1!D)"
xaxistitle = "!6Time (Gyr)"
xmax = 6.0
xmin = 0

x0= 0.11
x2= 0.99
xs= (x2-x0)/2.
x1= x0+xs

y0= 0.09
y2= 0.99
ys= (y2-y0)/2.
y1= y0+ys


;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  111111111111111111
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G2-u3"      & mmerg= 6.5    &  ymax = 12.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

;frunmerg= "G3G2-u4"      & mmerg= 6.5    &  ymax = 28.0
;fruniso1= "G3_n0"        & miso1= 5.0
;fruniso2= "G2im-u2"      & miso2= 1.5   

;frunmerg= "G3Rd4eG2"      & mmerg= 5.5     &  ymax = 12.0
;fruniso1= "G3Rd4e"        & miso1= 5.0
;fruniso2= "G2im-u1a"        & miso2= 0.5

;frunmerg= "G3Rd4eG2_n0"      & mmerg= 5.5     &  ymax = 32.0
;fruniso1= "G3Rd4e_n0"        & miso1= 5.0
;fruniso2= "G2im-u2"        & miso2= 0.5

;frunmerg= "G3Rd4hG2_n0"      & mmerg= 6.5    &  ymax = 32.0
;fruniso1= "G3Rd4h_n0"        & miso1= 5.0
;fruniso2= "G2im-u2"          & miso2= 1.5   

; use this one
frunmerg= "G3Rd4iG2_n0"      & mmerg= 6.5    &  ymax = 32.0 ; & ymax = 28.0
fruniso1= "G3Rd4i_n0"        & miso1= 5.0
fruniso2= "G2im-u2"     & miso2= 1.5   


;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


print, "-------------------------------------------------------------"
print, "WARNING: manually fixing this so the range looks good."
print, "it doesn't change the results that we present, and also saves"
print, "us from having to change the exes."
print, "-------------------------------------------------------------"
idx=where(sfrm gt ymax)
if idx(0) ne -1 then sfrm(idx)= ymax

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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x0+0.12, 0.96, "G3G2, 2.3:1", /normal, charthick=4, size=1.33, color=0
        ;xyouts, x0+0.12, 0.96, frunmerg+', 2.3:1', /normal, charthick=4, size=1.33, color=0

	xyouts,  x0+0.05, 0.90, '!7a!6= 1', /normal, charthick=4, size=2.2, color=0
	;xyouts,  x0+0.05, 0.90, '!7a!6= 3', /normal, charthick=4, size=2.2, color=0
	;xyouts,  x0+0.05, 0.85, 'with', /normal, charthick=4, size=2.2, color=0
	;xyouts,  x0+0.03, 0.80, 'R!Dcutoff!N=', /normal, charthick=4, size=2.2, color=0
	;xyouts,  x0+0.06, 0.75, '20 kpc', /normal, charthick=4, size=2.2, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100





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
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
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
        xyouts, 2.5, 0.14, belbl, /data, charthick=4, size=2.2, color=0




;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  22222222222222
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3Rd4G2"      & mmerg= 6.5    &  ymax = 12.0
;fruniso1= "G3Rd4"        & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

;frunmerg= "G3Rd4iG2_n0"      & mmerg= 6.5    &  ymax = 32.0 ; & ymax = 28.0
;fruniso1= "G3Rd4i_n0"        & miso1= 5.0
;fruniso2= "G2im-u2"     & miso2= 1.5   

; use this one
frunmerg= "G3Rd4hG2_n0"      & mmerg= 6.5    &  ymax = 32.0
fruniso1= "G3Rd4h_n0"        & miso1= 5.0
fruniso2= "G2im-u2"          & miso2= 1.5   



;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x1, y1, x2, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x1+0.12, 0.96, "G3Rd4G2, 1:2.3", /normal, charthick=4, size=1.33, color=0
        ;xyouts, x1+0.12, 0.96, frunmerg+', 2.3:1', /normal, charthick=4, size=1.33, color=0

	;xyouts,  x1+0.05, 0.90, '!7a!6= 1', /normal, charthick=4, size=2.2, color=0
	xyouts,  x1+0.05, 0.90, '!7a!6= 3', /normal, charthick=4, size=2.2, color=0
	xyouts,  x1+0.05, 0.85, 'with', /normal, charthick=4, size=2.2, color=0
	xyouts,  x1+0.03, 0.80, 'R!Dcutoff!N=', /normal, charthick=4, size=2.2, color=0
	xyouts,  x1+0.06, 0.75, '20 kpc', /normal, charthick=4, size=2.2, color=0


        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100






;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x1, y0, x2, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

	; manual manipulation
	sfrm= sfrm - (1.0-sfri)*0.5
	sfri=sfri - (1.0-sfri)*0.5


        lastn= n_elements(sfri)

        ;-------------------------------------------

        x= [xmin, timem, xmax, reverse(time1)]
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
        xyouts, 2.5, 0.14, belbl, /data, charthick=4, size=2.2, color=0




;--------------------------------------
;--------------------------------------



;--------------------------------------
;--------------------------------------
; done
; ----------
device, /close

end
















;======================================================================
;      PLOT FORMAT
;
;
;      |-------------|-------------|-------------|
;      |             |             |             |
;      |             |             |             |
;      |    SFR      |             |             |
;      |             |             |             |
;      |             |             |             |
;      |-------------|-------------|-------------|
;      |             |             |             |
;      |             |             |             |
;      |   Gas Mass  |             |             |
;      |             |             |             |
;      |             |             |             |
;      |-------------|-------------|-------------|
;
;
;
;
;   currently is not figure 18 in minors paper
;
;
;======================================================================


pro msfr_notf18, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr_notf18, junk"
   print, "  "
   print, "  "
   return
endif


filename='minor_notf18.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=28, newysize=18
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------

;xaxistitle = "Time (Gyr h!E-1!D)"
xaxistitle = "!6Time (Gyr)"
xmax = 6.0
xmin = 0

x0= 0.10
x3= 0.99
xs= (x3-x0)/3.
x1= x0+xs
x2= x0+xs+xs

y0= 0.12
y1= 0.555
y2= 0.99


;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  111111111111111111
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G2-u3"      & mmerg= 6.5    &  ymax = 20.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

;frunmerg= "G3G2-u4"      & mmerg= 6.5    &  ymax = 28.0
;fruniso1= "G3_n0"        & miso1= 5.0
;fruniso2= "G2im-u2"      & miso2= 1.5   

;frunmerg= "G3Rd4eG2"      & mmerg= 5.5     &  ymax = 12.0
;fruniso1= "G3Rd4e"        & miso1= 5.0
;fruniso2= "G2im-u1a"        & miso2= 0.5

;frunmerg= "G3Rd4eG2_n0"      & mmerg= 5.5     &  ymax = 32.0
;fruniso1= "G3Rd4e_n0"        & miso1= 5.0
;fruniso2= "G2im-u2"        & miso2= 0.5

frunmerg= "G3Rd4hG2_n0"      & mmerg= 6.5    &  ymax = 32.0 
fruniso1= "G3Rd4h_n0"        & miso1= 5.0
fruniso2= "G2im-u2"          & miso2= 1.5   



;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=2.10, ycharsize=2.10, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x0+0.12, 0.96, "G3G2, 2.3:1", /normal, charthick=4, size=1.33, color=0
        ;xyouts, x0+0.12, 0.96, frunmerg+', 2.3:1', /normal, charthick=4, size=1.33, color=0


        xyouts,  x1+0.05, 0.90, '!7a!6= 3', /normal, charthick=4, size=2.2, color=0


        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100





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
	xcharsize=2.10, ycharsize=2.10, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100


	; calculate burst efficiency
	be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= strmid(belbl,0,4)
        xyouts, 4.0, 0.21, belbl, /data, charthick=4, size=2.6, color=0




;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  22222222222222
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3G2-u3"      & mmerg= 6.5    &  ymax = 20.0
;fruniso1= "G3il-u1a"     & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

;frunmerg= "G3G2-u4"      & mmerg= 6.5    &  ymax = 28.0
;fruniso1= "G3_n0"        & miso1= 5.0
;fruniso2= "G2im-u2"      & miso2= 1.5   

;frunmerg= "G3Rd4eG2"      & mmerg= 5.5     &  ymax = 12.0
;fruniso1= "G3Rd4e"        & miso1= 5.0
;fruniso2= "G2im-u1a"        & miso2= 0.5

;frunmerg= "G3Rd4eG2_n0"      & mmerg= 5.5     &  ymax = 32.0
;fruniso1= "G3Rd4e_n0"        & miso1= 5.0
;fruniso2= "G2im-u2"        & miso2= 0.5

frunmerg= "G3Rd4hG2_n0"      & mmerg= 6.5    &  ymax = 32.0
fruniso1= "G3Rd4h_n0"        & miso1= 5.0
fruniso2= "G2im-u2"          & miso2= 1.5



;--------------------------------------


;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x1, y1, x2, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x1+0.12, 0.96, "G3Rd4G2, 1:2.3", /normal, charthick=4, size=1.33, color=0
        ;xyouts, x1+0.12, 0.96, frunmerg+', 2.3:1', /normal, charthick=4, size=1.33, color=0

        xyouts,  x1+0.05, 0.90, '!7a!6= 3', /normal, charthick=4, size=2.2, color=0

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100






;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x1, y0, x2, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=2.10, ycharsize=2.10, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



	; calculate burst efficiency
	be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= strmid(belbl,0,4)
        xyouts, 4.0, 0.21, belbl, /data, charthick=4, size=2.6, color=0



;--------------------------------------
;--------------------------------------
;
;
;     COLUMN  33333333333333
;
;
;--------------------------------------
;--------------------------------------

;frunmerg= "G3Rd4eG2"      & mmerg= 5.5     &  ymax = 28.0
;fruniso1= "G3Rd4e"        & miso1= 5.0
;fruniso2= "G1i-u2"        & miso2= 0.5

;frunmerg= "G3Rd4G2"      & mmerg= 6.5    &  ymax = 12.0
;fruniso1= "G3Rd4"        & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   

;frunmerg= "G3gf2G2-u1"      & mmerg= 6.5    &  ymax = 20.0
;fruniso1= "G3gf2i-u1"        & miso1= 5.0
;fruniso2= "G2im-u1a"     & miso2= 1.5   


frunmerg= "G3Rd4iG2_n0"      & mmerg= 6.5    &  ymax = 32.0 & ymax = 28.0
fruniso1= "G3Rd4i_n0"        & miso1= 5.0
fruniso2= "G2im-u2"     & miso2= 1.5


;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 15.0
ymin = 0

!p.position= [x2, y1, x3, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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
        ;oplot, [0.3, 0.7], [7.1, 7.1], psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, x2+0.11, 0.96, "G3Rd4eG2, 1:5.8", /normal, charthick=4, size=1.33, color=0
        ;xyouts, x2+0.12, 0.96, frunmerg+', 2.3:1', /normal, charthick=4, size=1.33, color=0


        xyouts,  x1+0.05, 0.90, '!7a!6= 1', /normal, charthick=4, size=2.2, color=0


        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 9.0
        ;oplot, [0.3, 0.7], [6.0, 6.0], psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.29, 0.93, "G3 + G2 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G1 isolated", /normal, charthick=4, size=1.33, color=100
        ;xyouts, 0.29, 0.93, "G3 + G0 isolated", /normal, charthick=4, size=1.33, color=100




;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6New Star Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
;yaxistitle="!6Gas Mass / Initial Gas Mass "
yaxistitle="!6M!Dgas!N / M!Dgas!N(0)"
ymax = 1.1
ymin = 0

!p.position= [x2, y0, x3, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=2.10, ycharsize=2.10, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


        ;-------------------------------------------
        ;
	;  Merger
        ;
        ;-------------------------------------------
        minor_open_sfr_file, frunmerg, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrmfs)
        ;oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass)
        oldtime= transpose(sfrtime) & oldsfr= transpose(sfrgasmass) & oldsfr= oldsfr/oldsfr[0]
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=100


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

        x= [xmin, timem, xmax, reverse(time1)]
        y= [sfri[0], sfrm, sfri[lastn-1], reverse(sfri)]
        polyfill, x, y, /data, color= 200, /fill, thick= 3.0

        oplot, timem, sfrm, psym=-3, color= 150, linestyle= 0, thick= 6.0
        ;xyouts, 0.7, 0.80, "merger", /normal, charthick=4, size=1.33, color=150

        oplot, time1, sfri, psym=-3, color= 100, linestyle= 1, thick= 6.0
        ;xyouts, 0.7, 0.70, "combined isolated", /normal, charthick=4, size=1.33, color=100



	; calculate burst efficiency
	be= sfri[99] - sfrm[99]      ; order is reversed because sfrx is gas remaining rather than consumed
        belbl= strcompress(string(be),/remove_all)
        belbl= strmid(belbl,0,4)
        xyouts, 4.0, 0.21, belbl, /data, charthick=4, size=2.6, color=0



;--------------------------------------
;--------------------------------------



;--------------------------------------
;--------------------------------------
; done
; ----------
device, /close

end










;======================================================================
;      PLOT FORMAT
;
;
;      |----------------------|
;      |                      |
;      |                      |
;      |                      |
;      |      SFR (G3Gx)      |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;      |                      |
;      |                      |
;      |                      |
;      |      SFR (G2Gx)      |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;      |                      |
;      |                      |
;      |                      |
;      |      SFR (G1Gx)      |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;      |                      |
;      |                      |
;      |                      |
;      |      SFR (G0Gx)      |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;
;
;
;   new fig. 7  (Fig. 7 in minors)
;
;
;======================================================================


pro msfr_f7_new, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "msfr_f7_new, junk"
   print, "  "
   print, "  "
   return
endif


filename='minorsfr_f7new.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=10, newysize=22



;-------------------------------------------
;  set variables
;-------------------------------------------

xaxistitle = "!6Time (Gyr)"
;xmax = 6.0
;xmax = 5.0
xmax = 4.6
xmin = 0

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
ymax= 20.0
ymin= 0.0


x0= 0.18
x1= 0.97

y0= 0.08
y1= 0.3275
y2= 0.55
y3= 0.7725
y4= 0.99



; makeitlog= 0
makeitlog= 1




;--------------------------------------
; Top Plot
;--------------------------------------

ymax= 28.0
frun1= "G3G3b-u2" & lbl1= "G3G3, 1:1"
frun2= "G3G2-u4"  & lbl2= "G3G2, 2.3:1"
frun3= "G3G1-u4"  & lbl3= "G3G1, 5.8:1"
frun4= "G3G0e-u4" & lbl4= "G3G0, 22.7:1"


!p.position= [x0, y3, x1, y4]

if not makeitlog eq 1 then begin
	plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle

	process_one_minor_sfr, frun1, lbl1, 0, 0.95, ymax*0.92, xmax=xmax, xlbl= 0.25
	process_one_minor_sfr, frun2, lbl2, 50, 0.92, ymax*0.87, xmax=xmax, xlbl= 0.25
	process_one_minor_sfr, frun3, lbl3, 150, 0.89, ymax*0.82, xmax=xmax, xlbl= 0.25
	process_one_minor_sfr, frun4, lbl4, 200, 0.86, ymax*0.77, xmax=xmax, xlbl= 0.25

endif else begin
	ymax= 30.0
	ymin= 1.5e-2
	plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle, /noerase, ytickformat='exp_label', /ylog

	process_one_minor_sfr, frun1, lbl1, 0, 0.86, ymax*0.92, xmax=xmax, xlbl= 0.25
	process_one_minor_sfr, frun2, lbl2, 50, 0.83333, ymax*0.87, xmax=xmax, xlbl= 0.25
	process_one_minor_sfr, frun3, lbl3, 150, 0.80667, ymax*0.82, xmax=xmax, xlbl= 0.25
	process_one_minor_sfr, frun4, lbl4, 200, 0.78, ymax*0.77, xmax=xmax, xlbl= 0.25

endelse




;--------------------------------------
; Top - Middle  Plot
;--------------------------------------

ymax= 18.0
frun1= "G2G2-u2" & lbl1= "G2G2, 1:1"
frun2= "G2G1-u4" & lbl2= "G2G1, 2.6:1"
frun3= "G2G0-u4" & lbl3= "G2G0, 10.0:1"

!p.position= [x0, y2, x1, y3]

if not makeitlog eq 1 then begin
	plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle, /noerase
endif else begin
	ymax= 40.0
	ymin= 2.0e-2
	plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle, /noerase, ytickformat='exp_label', /ylog
endelse

process_one_minor_sfr, frun1, lbl1, 0, 0.74, ymax*0.92, xmax=xmax, xlbl= 0.6
process_one_minor_sfr, frun2, lbl2, 50, 0.71, ymax*0.87, xmax=xmax, xlbl= 0.6
process_one_minor_sfr, frun3, lbl3, 150, 0.68, ymax*0.82, xmax=xmax, xlbl= 0.6



;--------------------------------------
; Middle-Bottom Plot
;--------------------------------------

ymax= 13.0
frun1= "G1G1a-u2" & lbl1= "G1G1, 1:1"
frun2= "G1G0-u4"  & lbl2= "G1G0, 3.9:1"

!p.position= [x0, y1, x1, y2]

if not makeitlog eq 1 then begin
	plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle, /noerase
        ;xtickformat='(a1)', ytitle=yaxistitle, /noerase, ytickformat='exp_label', /ylog
endif else begin
	ymax= 40.0
	ymin= 1.0e-2
	plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle, /noerase, ytickformat='exp_label', /ylog
endelse


process_one_minor_sfr, frun1, lbl1, 0, 0.52, ymax*0.92, xmax=xmax, xlbl= 0.6
process_one_minor_sfr, frun2, lbl2, 50, 0.49, ymax*0.87, xmax=xmax, xlbl= 0.6


;--------------------------------------
; Bottom Plot
;--------------------------------------

ymax= 3.9
frun1= "G0G0a-u2" & lbl1= "G0G0, 1:1"

!p.position= [x0, y0, x1, y1]

if not makeitlog eq 1 then begin
	plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle
endif else begin
	ymax= 8.0
	ymin= 4.0e-4
	plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtitle=xaxistitle, ytitle=yaxistitle, /noerase, ytickformat='exp_label', /ylog
endelse


process_one_minor_sfr, frun1, lbl1, 0, 0.29, ymax*0.92, xmax=xmax, xlbl= 0.6


;--------------------------------------
;--------------------------------------




; done
; ----------
device, /close

end









