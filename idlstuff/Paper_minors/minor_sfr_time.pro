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






;
pro doit, junk


sfr, "G3G3b-u1", "G3G3, 1:1",  2.7, 15.0, 6.0, 'minorsfr_g3g3.eps'
sfr, "G3G2-u3", "G3G2, 2.3:1", 3.0, 15.0, 6.0, 'minorsfr_g3g2.eps'
sfr, "G3G1-u3", "G3G1, 5.8:1", 4.0, 15.0, 6.0, 'minorsfr_g3g1.eps'
sfr, "G3G0e-u3", "G3G0, 22.7:1", 4.0, 15.0, 6.0, 'minorsfr_g3g0.eps'

sfr, "G2G2-u1", "G2G2, 1:1",  1.4, 12.0, 6.0, 'minorsfr_g2g2.eps'
sfr, "G2G1-u3", "G2G1, 2.6:1", 1.8, 12.0, 6.0, 'minorsfr_g2g1.eps'
sfr, "G2G0-u3", "G2G0, 10.0:1", 2.5, 12.0, 6.0, 'minorsfr_g2g0.eps'

sfr, "G1G1a-u1", "G1G1, 1:1", 1.2, 5.0, 6.0, 'minorsfr_g1g1.eps'
sfr, "G1G0-u3", "G1G0, 3.9:1", 2.1, 5.0, 6.0, 'minorsfr_g1g0.eps'

sfr, "G0G0a-u1", "G0G0, 1:1", 1.8, 1.2, 6.0, 'minorsfr_g0g0.eps'




;-------------------


;sfr, "G3G3b-u2", "G3G3, 1:1",  2.6, 40.0, 6.0, 'minorsfr_g3g3.eps'
;sfr, "G3G2-u4", "G3G2, 2.3:1", 3.0, 40.0, 6.0, 'minorsfr_g3g2.eps'
;sfr, "G3G1-u4", "G3G1, 5.8:1", 4.0, 40.0, 6.0, 'minorsfr_g3g1.eps'
;sfr, "G3G0e-u4", "G3G0, 22.7:1", 4.0, 40.0, 6.0, 'minorsfr_g3g0.eps'

;sfr, "G2G2-u2", "G2G2, 1:1",  1.4, 38.0, 6.0, 'minorsfr_g2g2.eps'
;sfr, "G2G1-u4", "G2G1, 2.6:1", 1.6, 38.0, 6.0, 'minorsfr_g2g1.eps'
;sfr, "G2G0-u4", "G2G0, 10.0:1", 2.5, 38.0, 6.0, 'minorsfr_g2g0.eps'

;sfr, "G1G1a-u2", "G1G1, 1:1", 1.2, 10.0, 6.0, 'minorsfr_g1g1.eps'
;sfr, "G1G0-u4", "G1G0, 3.9:1", 2.0, 10.0, 6.0, 'minorsfr_g1g0.eps'

;sfr, "G0G0a-u2", "G0G0, 1:1", 1.5, 3.4, 6.0, 'minorsfr_g0g0.eps'




end







;======================================================================
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

pro sfr, frun, lbl, t_merg, ymax, xmax, filename, quiet=quiet
;pro junk

if not keyword_set(frun) then begin
;if not keyword_set(junk) then begin
   print, "  "
   print, "sfr, junk"
   print, "  "
   print, "  "
   return
endif



;-------------------------------------------
;  set variables
;-------------------------------------------


;filename='minorsfr_time.eps'


;frun= "G3G3b-u1" & lbl= "G3G3, 1:1" & t_merg= 2.6  & ymax= 15.0  & xmax= 6.0 & filename='minorsfr_g3g3.eps'
;frun= "G3G2-u3"  & lbl= "G3G2, 2.3:1" & t_merg= 3.0  & ymax= 15.0  & xmax= 6.0 & filename='minorsfr_g3g2.eps'
;frun= "G3G1-u3"  & lbl= "G3G1, 5.8:1" & t_merg= 4.0  & ymax= 15.0  & xmax= 6.0 & filename='minorsfr_g3g1.eps'
;frun= "G3G0e-u3" & lbl= "G3G0, 22.7:1" & t_merg= 4.0  & ymax= 15.0  & xmax= 6.0 & filename='minorsfr_g3g0.eps'


;frun1= "G3G3b-u2" & lbl1= "G3G3, 1:1"
;frun= "G3G2-u4"  & lbl= "G3G2, 2.3:1"  & t_merg= 2.8  & ymax= 40.0  & xmax= 6.0
;frun3= "G3G1-u4"  & lbl3= "G3G1, 5.8:1"
;frun4= "G3G0e-u4" & lbl4= "G3G0, 22.7:1"

;filename='minorsfr_g2.eps' & ymax= 12.0 & xmax= 4.2
;frun1= "G2G2-u1" & lbl1= "G2G2, 1:1"
;frun2= "G2G1-u3" & lbl2= "G2G1, 2.6:1"
;frun3= "G2G0-u3" & lbl3= "G2G0, 10.0:1"
;filename='minorsfr_g2.eps' & ymax= 38.0 & xmax= 4.2
;frun1= "G2G2-u2" & lbl1= "G2G2, 1:1"
;frun2= "G2G1-u4" & lbl2= "G2G1, 2.6:1"
;frun3= "G2G0-u4" & lbl3= "G2G0, 10.0:1"

;filename='minorsfr_g1.eps' & ymax= 5.0 & xmax= 4.2
;frun1= "G1G1a-u1" & lbl1= "G1G1, 1:1"
;frun2= "G1G0-u3"  & lbl2= "G1G0, 3.9:1"
;filename='minorsfr_g1.eps' & ymax= 10.0 & xmax= 4.2
;frun1= "G1G1a-u2" & lbl1= "G1G1, 1:1"
;frun2= "G1G0-u4"  & lbl2= "G1G0, 3.9:1"

;filename='minorsfr_g0.eps' & ymax= 1.2 & xmax= 4.2
;frun1= "G0G0a-u1" & lbl1= "G0G0, 1:1"
;filename='minorsfr_g0.eps' & ymax= 3.4 & xmax= 4.2
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

;nsmths= 200
nsmths= 100

        if not keyword_set(old) then begin
                minor_open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        endif else begin
                open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        endelse
;oplot, sfrtime, sfrsfr, psym=-3, color=100, linestyle=0, thick=3.0
        oldtime= transpose(sfrtime)
        oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=nsmths

        if lbl ne 'none' then begin
                oplot, timem, sfrm, psym=-3, color= 0, linestyle= 0, thick= 12.0
                ;oplot, [xmax*0.548,xmax*0.613], [ylbl_key,ylbl_key], psym=-3, color= lcolor, linestyle= 0, thick=6.0
                ;xyouts, 0.7, 0.85, lbl, /normal, charthick=4, size=1.33, color=lcolor
                xyouts, 0.22, 0.85, lbl, /normal, charthick=4, size=2.0, color=0

                xyouts, 0.22, 0.78, 'Nuclear', /normal, charthick=2, size=1.8, color=0
                xyouts, 0.22, 0.73, 'Coalescence', /normal, charthick=2, size=1.8, color=0
        endif else begin
                oplot, timem, sfrm, psym=-3, color= 0, linestyle= 1, thick= 6.0
        endelse






;-------------------------------------------
;
;-------------------------------------------
;do_sb_fit= 0
do_sb_fit= 1
if do_sb_fit eq 1 then begin
	; ----------------
	; use full data
        ;t_merg= fload_bh_mergertime('/raid4/tcox/'+frun)
        ;idx= where((sfrtime gt t_merg-0.8) and (sfrtime lt t_merg+0.8))
        ;time= sfrtime(idx)
        ;sfr= sfrsfr(idx)
	; ----------------
	; use smoothed data
	idx= where((timem gt t_merg-0.8) and (timem lt t_merg+0.8))
        time= timem(idx)
        sfr= sfrm(idx)
        ; ----------------

        fitandoverplot_gaussian, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime, quiet=quiet
endif

do_sb_fit= 0
;do_sb_fit= 1
if do_sb_fit eq 1 then begin
        ;t_merg= fload_bh_mergertime('/raid4/tcox/'+frun)
        idx= where((sfrtime gt t_merg-0.5) and (sfrtime lt t_merg+0.5))
        time= sfrtime(idx)
        sfr= sfrsfr(idx)
        fitandoverplot_gaussian_pluszeropt, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime
endif


;do_expdecay_fit= 0
do_expdecay_fit= 1
if do_expdecay_fit eq 1 then begin
        ;idx= where(sfrtime gt SBtime)
        idx= where((sfrtime gt SBtime) and (sfrtime lt SBtime+1.3))
        time= sfrtime(idx)
        time= time- time(0)
        sfr= sfrsfr(idx)
        fitandoverplot_expdecay, time, sfr, SBtime=SBtime, tau=tau
endif

;do_fwhm_fit= 0
do_fwhm_fit= 1
if do_fwhm_fit eq 1 then begin
	; ----------------
	; use full data
        ;t_merg= fload_bh_mergertime('/raid4/tcox/'+frun)
        ;idx= where((sfrtime gt t_merg-0.5) and (sfrtime lt t_merg+0.5))
        ;time= sfrtime
        ;sfr= sfrsfr
	; ----------------
	; use smoothed data
        time= timem
        sfr= sfrm
	; ----------------
        findandoverplot_fwhm, time, sfr, $
                        ;SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, quiet=quiet
                        ;SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, /quiet
                        SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime
endif



;--------------------------------------
;--------------------------------------

device, /close


end









;======================================================================
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

pro sfr_firstp, frun, lbl, t_merg, ymax, xmax, filename, quiet=quiet
;pro junk

if not keyword_set(frun) then begin
;if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_firstp, junk"
   print, "  "
   print, "  "
   return
endif



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

;nsmths= 200
nsmths= 100

        if not keyword_set(old) then begin
                minor_open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        endif else begin
                open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
        endelse
;oplot, sfrtime, sfrsfr, psym=-3, color=100, linestyle=0, thick=3.0
        oldtime= transpose(sfrtime)
        oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=nsmths

        if lbl ne 'none' then begin
                oplot, timem, sfrm, psym=-3, color= 0, linestyle= 0, thick= 12.0
                ;oplot, [xmax*0.548,xmax*0.613], [ylbl_key,ylbl_key], psym=-3, color= lcolor, linestyle= 0, thick=6.0
                ;xyouts, 0.7, 0.85, lbl, /normal, charthick=4, size=1.33, color=lcolor
                xyouts, 0.22, 0.85, lbl, /normal, charthick=4, size=2.0, color=0

                xyouts, 0.22, 0.78, 'First', /normal, charthick=2, size=1.8, color=0
                xyouts, 0.22, 0.73, 'Passage', /normal, charthick=2, size=1.8, color=0
        endif else begin
                oplot, timem, sfrm, psym=-3, color= 0, linestyle= 1, thick= 6.0
        endelse






;-------------------------------------------
;
;-------------------------------------------
;do_sb_fit= 0
do_sb_fit= 1
if do_sb_fit eq 1 then begin
	; ----------------
	; use smoothed data, and do for first pass
	idx= where(timem lt t_merg-0.5)
        time= timem(idx)
        sfr= sfrm(idx)
        ; ----------------

        fitandoverplot_gaussian, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime, quiet=quiet
endif

do_sb_fit= 0
;do_sb_fit= 1
if do_sb_fit eq 1 then begin
        ;t_merg= fload_bh_mergertime('/raid4/tcox/'+frun)
        idx= where((sfrtime gt t_merg-0.5) and (sfrtime lt t_merg+0.5))
        time= sfrtime(idx)
        sfr= sfrsfr(idx)
        fitandoverplot_gaussian_pluszeropt, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime
endif


;do_expdecay_fit= 0
do_expdecay_fit= 1
if do_expdecay_fit eq 1 then begin
        ;idx= where(sfrtime gt SBtime)
        idx= where((sfrtime gt SBtime) and (sfrtime lt SBtime+1.0))
        time= sfrtime(idx)
        time= time- time(0)
        sfr= sfrsfr(idx)
        fitandoverplot_expdecay, time, sfr, SBtime=SBtime, tau=tau
endif

;do_fwhm_fit= 0
do_fwhm_fit= 1
if do_fwhm_fit eq 1 then begin
	; ----------------
	; use smoothed data
	idx= where(timem lt t_merg-0.5)
        time= timem(idx)
        sfr= sfrm(idx)
	; ----------------
        findandoverplot_fwhm, time, sfr, $
                        ;SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, quiet=quiet
                        ;SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, /quiet
                        SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime
endif



;--------------------------------------
;--------------------------------------

device, /close


end













;======================================================================
;
;
;  TEST   TEST    TEST   TEST   TEST
;  TEST   TEST    TEST   TEST   TEST
;  TEST   TEST    TEST   TEST   TEST
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

pro test, junk


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename='test.eps', colortable= 4



; --------------------------

yaxistitle="  "
xaxistitle = "  "

xmax = 7.0
xmin = 0.0

ymax = 1.5
ymin = 0.0



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





;
;  define gaussian
; --------------------------

x= findgen(1000)/999. * (xmax - xmin) + xmin

sigma= 0.3
x0= 3.5
y= exp( - (x-x0) * (x-x0) / 2. / sigma / sigma)


;y(where(y lt 1.0e-3))= 0.0


oplot, x, y, psym=-3, color= 0, linestyle= 0, thick= 2.0


; --------------------------




t_merg= 3.5
timem= x
sfrtime= x
sfrm= y
sfrsfr= y


;-------------------------------------------
;
;-------------------------------------------
;do_sb_fit= 0
do_sb_fit= 1
if do_sb_fit eq 1 then begin
	; ----------------
	; use full data
        ;t_merg= fload_bh_mergertime('/raid4/tcox/'+frun)
        ;idx= where((sfrtime gt t_merg-0.8) and (sfrtime lt t_merg+0.8))
        ;time= sfrtime(idx)
        ;sfr= sfrsfr(idx)
	; ----------------
	; use smoothed data
	idx= where((timem gt t_merg-3.8) and (timem lt t_merg+3.8))
        time= timem(idx)
        sfr= sfrm(idx)
        ; ----------------

        fitandoverplot_gaussian, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime, quiet=quiet
endif

do_sb_fit= 0
;do_sb_fit= 1
if do_sb_fit eq 1 then begin
        ;t_merg= fload_bh_mergertime('/raid4/tcox/'+frun)
        idx= where((sfrtime gt t_merg-0.5) and (sfrtime lt t_merg+0.5))
        time= sfrtime(idx)
        sfr= sfrsfr(idx)
        fitandoverplot_gaussian_pluszeropt, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime
endif


;do_expdecay_fit= 0
do_expdecay_fit= 1
if do_expdecay_fit eq 1 then begin
        ;idx= where(sfrtime gt SBtime)
        idx= where((sfrtime gt SBtime) and (sfrtime lt SBtime+3.3))
        time= sfrtime(idx)
        time= time- time(0)
        sfr= sfrsfr(idx)
        fitandoverplot_expdecay, time, sfr, SBtime=SBtime
endif

;do_fwhm_fit= 0
do_fwhm_fit= 1
if do_fwhm_fit eq 1 then begin
	; ----------------
	; use full data
        ;t_merg= fload_bh_mergertime('/raid4/tcox/'+frun)
        ;idx= where((sfrtime gt t_merg-0.5) and (sfrtime lt t_merg+0.5))
        ;time= sfrtime
        ;sfr= sfrsfr
	; ----------------
	; use smoothed data
        time= timem
        sfr= sfrm
	; ----------------
        findandoverplot_fwhm, time, sfr, $
                        ;SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, quiet=quiet
                        SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, /quiet
endif



;--------------------------------------
;--------------------------------------

device, /close


end








;======================================================================
;
;
;
;
;
;======================================================================




; Fit Gaussian to Starburst
;-----------------------------
pro fitandoverplot_gaussian, time, sfr, $
			SFR_max= SFR_max, width=width, SBtime=SBtime, $
			quiet=quiet, noplot=noplot

        x_tofit= time
        y_tofit= sfr
        weight_tofit= 1.0 + 0.0*sfr


        ;  parameters
        ;---------------- 
	; p[0] = Normalization, i.e., the SFR_max
        ; p[1] = Width of the dist.; sigma
        ; p[2] = mean of the dist., i.e., time of SFR_max

        ; initial guess
        ;guess= [10.0,0.2,1.5]
	normz= max(sfr)
	idd= where(sfr eq normz)
	meanz= time(idd[0])
	wid= max(time)/3.0
	guess= [normz, wid, meanz]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
        ; Normalization is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ; width is greater than 0
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.001
        ; mean is greater than 0
        pi[2].limited(0) = 1
        pi[2].limits(0) = 0.001


        ; markwardt mpfit procedure
        sb_result = MPFITFUN('func_gaussian', x_tofit, y_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

	print, "--------------------------"
        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        chi2= strcompress(string(redchi2),/remove_all)
        chi2= strmid(chi2,0,4)

        SFR_max= sb_result[0]
	print, "SFR_max= ", SFR_max
        sfrmaxlbl= strcompress(string(SFR_max),/remove_all)
        sfrmaxlbl= strmid(sfrmaxlbl,0,5)

        width= sb_result[1]
	print, "Starburst Width= ", width," Gyr"
        wlbl= strcompress(string(width),/remove_all)
        wlbl= strmid(wlbl,0,4)

        SBtime= sb_result[2]
	print, "Time at SFR_max= ", SBtime," Gyr"
	sbtm= strcompress(string(SBtime),/remove_all)
	sbtm= strmid(sbtm,0,4)


        ; overplot fit
        ; ----------
        ;x=min(time) - 0.2 + (0.2+max(time)-min(time))*findgen(100)/100.0
        ;x= SBtime + 1.5*(findgen(100)/100.0 - 0.5)
        x= SBtime + 4.0*(findgen(100)/100.0 - 0.5)
        y= func_gaussian(x,sb_result)
        if not keyword_set(quiet) then begin
		oplot, x, y, psym=-3, linestyle= 0, color= 150, thick=3.0
	endif


	if not keyword_set(quiet) then begin
		;xyouts, 0.22, 0.75, 'Gaussian ('+chi2+')', size=1.2, color=0, /normal
		;xyouts, 0.23, 0.71, 'SFRmax= '+sfrmaxlbl, size=1.0, color=0, /normal
		;xyouts, 0.23, 0.67, 'sigma=  '+wlbl, size=1.0, color=0, /normal
		;xyouts, 0.23, 0.63, 'sbtime= '+sbtm, size=1.0, color=0, /normal
		xyouts, 0.65, 0.84, 'Gaussian', size=1.6, color=0, /normal
		xyouts, 0.65, 0.80, '!7r!6= '+wlbl+' Gyr', size=1.5, color=0, /normal
	endif


        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end




;=====================================================================




; Fit Gaussian to Starburst
;-----------------------------
pro fitandoverplot_gaussian_pluszeropt, time, sfr, $
			SFR_max= SFR_max, width=width, SBtime=SBtime

        x_tofit= time
        y_tofit= sfr
        weight_tofit= 1.0 + 0.0*sfr


        ;  parameters
        ;---------------- 
	; p[0] = Normalization, i.e., the SFR_max
        ; p[1] = Width of the dist.; sigma
        ; p[2] = mean of the dist., i.e., time of SFR_max

        ; initial guess
        guess= [100.0,0.2,1.5,1.0]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
        ; Normalization is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ; width is greater than 0
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.01
        ; mean is greater than 0
        pi[2].limited(0) = 1
        pi[2].limits(0) = 0.001
        ; zero point is greater than 0
        pi[3].limited(0) = 1
        pi[3].limits(0) = 0.1



        ; markwardt mpfit procedure
        sb_result = MPFITFUN('func_gaussian_pluszeropt', x_tofit, y_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        ;chi2= strcompress(string(redchi2),/remove_all)
        ;chi2= strmid(chi2,0,4)

        SFR_max= sb_result[0]
	print, "SFR_max= ", SFR_max
        ;rho0lbl= strcompress(string(rho0lbl),/remove_all)
        ;rho0lbl= strmid(rho0lbl,0,5)

        width= sb_result[1]
	print, "Starburst Width= ", width," Gyr"
        ;rslbl= strcompress(string(rslbl),/remove_all)
        ;rslbl= strmid(rslbl,0,5)

        SBtime= sb_result[2]
	print, "Time at SFR_max= ", SBtime," Gyr"

        zeropt= sb_result[3]
	print, "Zero Pt= ", zeropt," M_solar/Yr"


        ; overplot fit
        ; ----------
        x=min(time) - 0.2 + (0.2+max(time)-min(time))*findgen(100)/100.0
        y= func_gaussian_pluszeropt(x,sb_result)
        oplot, x, y, psym=-3, linestyle= 0, color= 200, thick=2.0

        ;xyouts, 0.65, 0.75, 'NFW', /normal, size= 1.0, color=150
        ;xyouts, 0.65, 0.71, 'r!Ds!N='+rslbl, /normal, size= 1.0, color=150
        ;xyouts, 0.77, 0.71, ', c='+clbl, /normal, size= 1.0, color=150


        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end





;====================================================================




; Fit Exponential to Fading Starburst
;--------------------------------------
pro fitandoverplot_expdecay, time, sfr, $
			SFR_SB= SFR_SB, tau=tau, SBtime=SBtime, $
			quiet=quiet

        x_tofit= time
        y_tofit= sfr
        ;weight_tofit= 1.0 + 0.0*sfr
        weight_tofit= 0.8*sfr + 1.0
stop


        ;  parameters
        ;---------------- 
	; p[0] = Normalization, i.e., the SFR_max ( - p[2])
        ; p[1] = Exp. decay, i.e., tau
        ; p[2] = constant, so it decays to this, rather than 0

        ; initial guess
        guess= [100.0,1.0,0.1]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
        ; SFR_SB is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ; tau is greater than 0
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.001
        ; const is greater than 0
        pi[2].limited(0) = 1
        pi[2].limits(0) = 0.00001


        ; markwardt mpfit procedure
        tau_result = MPFITFUN('func_exp_sf', x_tofit, y_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        ;chi2= strcompress(string(redchi2),/remove_all)
        ;chi2= strmid(chi2,0,4)

        SFR_SB= tau_result[0]
	print, "SFR_SB= ", SFR_SB, " M_solar/Yr"
        ;rho0lbl= strcompress(string(rho0lbl),/remove_all)
        ;rho0lbl= strmid(rho0lbl,0,5)

        tau= tau_result[1]
	print, "SB tau= ", tau," Gyr"
        taulbl= strcompress(string(tau),/remove_all)
        taulbl= strmid(taulbl,0,4)

        decay_const= tau_result[2]
	print, "decay const= ", decay_const," M_solar/Yr"
        ;rslbl= strcompress(string(rslbl),/remove_all)
        ;rslbl= strmid(rslbl,0,5)


        ; overplot fit
        ; ----------
        x=time
        y= func_exp_sf(x,tau_result)
        x=x+SBtime
	if not keyword_set(quiet) then begin
        	oplot, x, y, psym=-3, linestyle= 0, color= 100
	endif

        ;xyouts, 0.65, 0.75, 'NFW', /normal, size= 1.0, color=150
        ;xyouts, 0.65, 0.71, 'r!Ds!N='+rslbl, /normal, size= 1.0, color=150
        ;xyouts, 0.77, 0.71, ', c='+clbl, /normal, size= 1.0, color=150

	if not keyword_set(quiet) then begin
		xyouts, 0.65, 0.51, 'Tau Model', color=0, size=1.5, /normal
		xyouts, 0.65, 0.46, '!7s!6!DSF!N= '+taulbl+' Gyr', color=0, size=1.5, /normal
	endif

        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end








;====================================================================





; Find FWHM to Starburst
;-----------------------------
pro findandoverplot_fwhm, time, sfr, $
			SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, $
			quiet=quiet

	print, "xxxxxxxxxxx from file xxxxxxxxxxx"
	Ntime= n_elements(time)
	print, "Ntime= ", Ntime

	SFR_max= max(sfr)
	print, "SFR_max= ", SFR_max
	maxsfr= strcompress(string(SFR_max),/remove_all)
	maxsfr= strmid(maxsfr,0,5)

	idx= where(sfr eq SFR_max)
	maxidx= idx(0)
	print, "maxidx= ", maxidx

	SBtime= time(maxidx)
	print, "SBtime= ", SBtime
	sbtm= strcompress(string(SBtime),/remove_all)
	sbtm= strmid(sbtm,0,4)

	if not keyword_set(quiet) then begin
		;xyouts, 0.22, 0.92, 'from sfr file', size=1.2, color=0, /normal
		;xyouts, 0.22, 0.88, 'max(SFR)= '+maxsfr, size=1.0, color=0, /normal
		;xyouts, 0.22, 0.84, 'time@max= '+sbtm, size=1.0, color=0, /normal
	endif

	if maxidx lt 2 or maxidx gt Ntime-3 then begin
		print, " "
		print, "Problem: SFR_max is at boundary of SB region"
		print, " "
		fwhm= 1.0e-2
		return
	endif

	SFR_half= SFR_max * 0.5
	print, "SFR_half= ", SFR_half

	print, "----------------------------"
	; sfr_half on the pre-SB side
	;pretime= time(0:idx(0))
	;presfr= sfr(0:idx(0))
	;halfidx= where(presfr gt SFR_half)
	;print, halfidx(0)
	;if halfidx(0) eq -1 then begin & fwhm= 1.0e-2 & return & endif
	;time1= time(halfidx(0))
	;print, "time1= ", time1
	for halfidx=maxidx-1, 0, -1 do if sfr(halfidx) lt SFR_half then break
	;if halfidx le 0 then begin & fwhm= 1.0e-2 & return & endif
	if halfidx le 0 then halfidx= 1

	x1= time(halfidx)
	x2= time(halfidx+1)
	y1= sfr(halfidx)
	y2= sfr(halfidx+1)
	m= (y2-y1)/(x2-x1)
	b= y2 - m*x2 ;& print, b
	;b= y1 - m*x1 & print, b
	time1= (SFR_half - b)/m
	print, "time1= ", time1
	;print, "other meth= ", time1

	; sfr_half on the post-SB side
	;posttime= time(idx(0):Ntime-1)
	;postsfr= sfr(idx(0):Ntime-1)
	;halfidx= where(postsfr lt SFR_half)
	;print, halfidx(0)
	;if halfidx(0) eq -1 then begin & fwhm= 1.0e-2 & return & endif
	;time2= time(halfidx(0))
	;print, "time2= ", time2
	for halfidx=maxidx+1, Ntime-1 do if sfr(halfidx) lt SFR_half then break
	;if halfidx le 0 then begin & fwhm= 1.0e-2 & return & endif
	if halfidx ge Ntime then begin
		halfidx= Ntime-1
		time2= time(halfidx)
	endif else begin
		x1= time(halfidx-1)
        	x2= time(halfidx)
        	y1= sfr(halfidx-1)
        	y2= sfr(halfidx)
        	m= (y2-y1)/(x2-x1)
        	b= y2 - m*x2 ;& print, b
        	;b= y1 - m*x1 & print, b
        	time2= (SFR_half - b)/m
	endelse
	print, "time2= ", time2
        ;print, "other meth= ", time2

        ; overplot fit
        ; ----------
        x= [time1,time1]
        y= [0.0,40.0]
        ;y= [0.0,24.0]
	if not keyword_set(quiet) then begin
        	oplot, x, y, psym=-3, linestyle= 1, color= 180, thick=6.0
	endif

        x= [time2,time2]
        y= [0.0,40.0]
        ;y= [0.0,24.0]
	if not keyword_set(quiet) then begin
        	oplot, x, y, psym=-3, linestyle= 1, color= 180, thick=6.0
	endif


	fwhm= time2-time1
	fwhmlbl= strcompress(string(fwhm),/remove_all)
	fwhmlbl= strmid(fwhmlbl,0,4)
	print, "fwhm= ", fwhm
	if not keyword_set(quiet) then begin
		;xyouts, 0.22, 0.50, 'FWHM= '+fwhmlbl+' Gyr', color=0, size=1.2, /normal
		xyouts, 0.65, 0.67, 'FWHM ', color=0, size=1.5, /normal
		xyouts, 0.65, 0.62, fwhmlbl+' Gyr', color=0, size=1.5, /normal
	endif
	print, "----------------------------"


end




;===============================================================================





; Fit to fwhm
;-----------------------------
pro fitandoverplot_sbtime, massratio, fwhm, fwhmwt, $
			fNorm= fNorm, alpha=alpha, $
			lcolor=lcolor

        x_tofit= massratio
        y_tofit= fwhm
        ;weight_tofit= 1.0 + 0.0*massratio
        if keyword_set(fwhmwt) then weight_tofit=fwhmwt else weight_tofit= 0.5*fwhm


        ;  parameters
        ;---------------- 
	; p[0] = Normalization, i.e., the major merger FWHM
        ; p[1] = mass ratio dependence

        ; initial guess
        guess= [0.5,1.0]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},2)
        ; Normalization is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 1.0e-10
        ; mass dependence is greater than 0
        ;pi[1].limited(0) = 1
        ;pi[1].limits(0) = 0.001



        ; markwardt mpfit procedure
        sb_result = MPFITFUN('minor_func_sbtime', x_tofit, y_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

	print, "--------------------------"
        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        chi2= strcompress(string(redchi2),/remove_all)
        chi2= strmid(chi2,0,4)

        fNorm= sb_result[0]
	print, "fNorm= ", fNorm
        fnormlbl= strcompress(string(fNorm),/remove_all)
        fnormlbl= strmid(fnormlbl,0,5)

        alpha= sb_result[1]
	print, "Alpha= ", alpha
        albl= strcompress(string(alpha),/remove_all)
        albl= strmid(albl,0,4)


        ; overplot fit
        ; ----------
        ;x=min(time) - 0.2 + (0.2+max(time)-min(time))*findgen(100)/100.0
        x= 10^(4.*(findgen(100)/100.0) - 3.0)
        y= minor_func_sbtime(x,sb_result)
        oplot, x, y, psym=-3, linestyle= 0, color= lcolor, thick=2.0

	;xyouts, 0.22, 0.75, 'MassDep ('+chi2+')', size=1.2, color=0, /normal
	;xyouts, 0.23, 0.71, 'fNorm= '+fnormlbl, size=1.0, color=0, /normal
	;xyouts, 0.23, 0.67, 'alpha=  '+albl, size=1.0, color=0, /normal

        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end




;=====================================================================




;=====================================================================




function grab_fwhm, fruns, sbts=sbts, sfrms=sfrms

nfs= n_elements(fruns)
fwhms= fltarr(nfs)
sbts= fltarr(nfs)
sfrms= fltarr(nfs)

for i=0, nfs-1 do begin
	minor_open_sfr_file, fruns[i], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass

        time= sfrtime
        sfr= sfrsfr
        findandoverplot_fwhm, time, sfr, $
                        SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime

	fwhms[i]= fwhm
	sbts[i]= SBtime
	sfrms[i]= SFR_max
endfor

return, fwhms

end





;=====================================================================




function grab_sbtimescale, fruns, sbtserr=sbtserr, firstpass=firstpass, T_cut= T_cut

nfs= n_elements(fruns)
sbtimescale= fltarr(nfs)
sbtserr= fltarr(nfs)

for i=0, nfs-1 do begin
	print, "processing: ", fruns[i]
        minor_open_sfr_file, fruns[i], sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass


	; fwhm
	; -----
	if not keyword_set(firstpass) then begin
           time= sfrtime
           sfr= sfrsfr
           findandoverplot_fwhm, time, sfr, $
                        SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, /quiet
	endif else begin
	   idx= where(sfrtime lt T_cut[i])
	   time= sfrtime(idx)
	   sfr= sfrsfr(idx)
           findandoverplot_fwhm, time, sfr, $
                        SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, /quiet
	endelse

	if fwhm gt 5.0 then fwhm= 5.0


	; gaussian
        ; ----------------
        ; use smoothed data
	nsmths= 100
	t_merg= SBtime
        oldtime= transpose(sfrtime)
        oldsfr= transpose(sfrsfr)
        resample_sfrhist, oldsfr, oldtime, sfrm, timem, Nelements=nsmths
	if not keyword_set(firstpass) then begin
           idx= where((timem gt t_merg-0.8) and (timem lt t_merg+0.8))
           time= timem(idx)
           sfr= sfrm(idx)
	endif else begin
           idx= where(timem lt T_cut[i])
           time= timem
           sfr= sfrm
	endelse
        ; ----------------

        fitandoverplot_gaussian, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime, /quiet

	if width gt 5.0 then fwhm= 5.0


	; exponential decay
	; --------------------
        ;idx= where(sfrtime gt SBtime)
        idx= where((sfrtime gt SBtime) and (sfrtime lt SBtime+1.8))
        time= sfrtime(idx)
        sfr= sfrsfr(idx)
        time= time- time(0)
        fitandoverplot_expdecay, time, sfr, SBtime=SBtime, tau=tau, /quiet

	if tau gt 5.0 then tau= 5.0

	;
	; tau = 1.1 * sig
	; fwhm = 2.35 * sig
	;	
	;
	; tau = 1.1 * sig = 1.1 * ( fwhm / 2.35 ) = 0.47 * fwhm
	;
	x= [0.47 * fwhm, 1.1 * width, tau]

	sbtimescale[i]= mean(x)
	sbtserr[i]= sqrt(variance(x))

endfor

return, sbtimescale

end










;======================================================================
;
;
;
;    ---------------------------
;    |                         |
;    |                         |
;    |                         |
;    |                         |
;    |                         |
;    |   Tsfr v. mass ratio    |
;    |                         |
;    |                         |
;    |                         |
;    |                         |
;    |                         |
;    ---------------------------
;
;
;
;======================================================================
pro fwhm, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "fwhm, junk"
   return
endif


;============================================
;============================================

;Isofruns= ["G3il-u1a", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
;G3fruns= ["G3G3b-u1", "G3G2-u3", "G3G1-u3", "G3G0e-u3"]

;Isofruns= ["G2im-u1a", "G1i-u1a", "G0i-u1a"]
;G2fruns= ["G2G2-u1", "G2G1-u3", "G2G0-u3"] 

;Isofruns= ["G1i-u1a", "G0i-u1a"]
;G1fruns= ["G1G1a-u1", "G1G0-u3"]

;Isofruns= ["G0i-u1a"]
;G0fruns= ["G0G0a-u1"]


; all of them
Mratios= 1.0/[1.0, 2.3, 5.8, 22.7, $
	1.0, 2.6, 10.0, $
	1.0, 3.9, $
	1.0]


Gstds= ["G3G3b-u1", "G3G2-u3", "G3G1-u3", "G3G0e-u3", $
	"G2G2-u1", "G2G1-u3", "G2G0-u3", $
	"G1G1a-u1", "G1G0-u3", $
	"G0G0a-u1"]


; now take out some of the non-burst models
Mratios= 1.0/[1.0, 2.3, 5.8, $
        1.0, 2.6, 10.0, $
        1.0, 3.9, $
        1.0]


Gstds= ["G3G3b-u1", "G3G2-u3", "G3G1-u3", $
        "G2G2-u1", "G2G1-u3", "G2G0-u3", $
        "G1G1a-u1", "G1G0-u3", $
        "G0G0a-u1"]


;fwhms_n2= grab_fwhm(Gstds, sbts=sbts, sfrms=sfrms)
FWHM_n2= grab_sbtimescale(Gstds, sbtserr=sbtserr)
FWHM_wt_n2= sbtserr


;print, "Mratios= ", Mratios
;print, "Gstds= ", Gstds
;print, "fwhms= ", fwhms_n2
;print, "sbts= ", sbts
;print, "sfrms= ", sfrms


; 
; get these by hand from the above sfr procedure.  the
; columns are
;
;	sfr_max from file
;	time of sfr_max from file
;	gaussian fit - sfr_max      (set to 1.00 if fit fails)
;	gaussian fit - sigma, width      (set to 1.00 if fit fails)
;	gaussian fit - time at sfr_max      (set to 1.00 if fit fails)
;	fwhm - 
;
;
;
;RunParameters= [13.29, 2.75, 9.80, 0.403, 2.822, 0.521, $       ; g3g3
;		5.00,  2.99, 4.08, 0.465, 3.075, 0.806, $       ; g3g2
;		1.54,  1.47, 1.00, 1.000, 1.000, 2.000, $       ; g3g1 - no real burst so fits suck
;		1.97,  0.00, 1.00, 1.000, 1.000, 2.000, $       ; g3g0 - no real burst so fits sucks
;		11.14, 1.40, 5.40, 0.437, 1.518, 0.588, $       ; g2g2 - too noisy for fwhm
;		4.02,  1.70, 2.95, 0.466, 1.735, 0.594, $       ; g2g1
;		0.74,  1.75, 1.00, 1.000, 1.000, 0.670, $       ; g2g0 - no real burst so fit sucks
;		4.15,  1.59, 2.12, 0.384, 1.547, 0.325, $       ; g1g1 - too noisy for fwhm
;		0.62,  1.87, 0.51, 0.527, 2.137, 0.758, $       ; g1g0
;		0.95,  1.81, 0.51, 0.394, 1.768, 0.553]         ; g0g0 - too noisy for fwhm
;

;
; XXXXX Original Version of Figure XXXXXX
; -----------------------------------------
;
; from minor_n2_[raw,smoothed].ps
;
;SBmax_n2= [12.5, 4.8, -1.0, -1.0, 8.0, 3.5, -1.0, 3.7, 0.6, 0.8]
;SBtime_n2= [2.7, 2.9, -1.0, -1.0, 1.5, 1.7, -1.0, 1.5, 2.0, 1.8]
;
;SBwidth_n2= [0.37, 0.46, 10.0, 10.0, 0.43, 0.46, 10.0, 0.36, 0.54, 0.40]
;FWHM_n2= [0.58, 0.85, 10.0, 10.0, 0.81, 0.85, 10.0, 0.75, 0.85, 0.80]
;
; rather than uniform weighting, or something that is based upon
; the fwhm, we use the difference between the fwhm and the width
;FWHM_wt_n2= abs(FWHM_n2 - SBwidth_n2) * 0.5
;;FWHM_wt_n2= abs(FWHM_n2 - SBwidth_n2) * 1.0

;
; conclusions?
; --------------
; maybe slight trend for larger sbwidth with smaller mass
; maybe slight trend for larger sbwidth with larger mass ratio
; tail of SF in G1's and below
;
;


;----------------------------------------------------------------

; all of them
Gstds= ["G3G3b-u2", "G3G2-u4", "G3G1-u4", "G3G0e-u4", $
	"G2G2-u2", "G2G1-u4", "G2G0-u4", $
	"G1G1a-u2", "G1G0-u4", $
	"G0G0a-u2"]

; on those with SB
Gstds= ["G3G3b-u2", "G3G2-u4", "G3G1-u4", $
	"G2G2-u2", "G2G1-u4", "G2G0-u4", $
	"G1G1a-u2", "G1G0-u4", $
	"G0G0a-u2"]

;fwhms_n0= grab_fwhm(Gstds, sbts=sbts, sfrms=sfrms)
FWHM_n0= grab_sbtimescale(Gstds, sbtserr=sbtserr)
FWHM_wt_n0= sbtserr
;print, "Mratios= ", Mratios
;print, "Gstds= ", Gstds
;print, "fwhms= ", fwhms_n0
;print, "sbts= ", sbts
;print, "sfrms= ", sfrms


; g3g2 is kind of strange here, and i think a bit
; underestimated, the errors too
FWHM_n0[1]= 0.1
FWHM_wt_n0[1]= 0.06
FWHM_wt_n0[5]= 0.2

;
; XXXXX Original Version of Figure XXXXXX
; -----------------------------------------
;
; from minor_n0_[raw,smoothed_50,100,200].ps
; i ended up liking the smoothed 100 best
;
;SBmax_n0= [26.0, 27.0, -1.0, -1.0, 14.1, 11.5, -1.0, 12.0, 0.7, 3.8]  ; final
;SBtime_n0= [2.5, 2.8, -1.0, -1.0, 1.25, 1.45, -1.0, 1.25, 1.85, 1.5]  ; final

; v0
;SBwidth_n0= [0.28, 0.03, 10.0, 10.0, 0.37, 0.41, 10.0, 0.09, 0.53, 0.08]
;FWHM_n0= [0.11, 0.08, 10.0, 10.0, 0.13, 0.06, 10.0, 0.10, 0.93, 0.14]
; v1
;SBwidth_n0= [0.28, 0.03, 10.0, 10.0, 0.37, 0.41, 10.0, 0.09, 0.53, 0.08]
;FWHM_n0= [0.11, 0.08, 10.0, 10.0, 0.13, 0.06, 10.0, 0.10, 0.25, 0.14]
; v2
;SBwidth_n0= [0.28, 0.01, 10.0, 10.0, 0.37, 0.41, 10.0, 0.09, 0.53, 0.08]  ; final
;FWHM_n0=    [0.11, 0.08, 10.0, 10.0, 0.13, 0.06, 10.0, 0.10, 0.25, 0.14]  ; final

; rather than uniform weighting, or something that is based upon
; the fwhm, we use the difference between the fwhm and the width
;FWHM_wt_n0= abs(FWHM_n0 - SBwidth_n0) * 0.5
;FWHM_wt_n0= abs(FWHM_n0 - SBwidth_n0) * 1.0

;
; conclusions?
; -------------
; not well fit by gaussian, often there is one large peak with
; a broad base of SF, or multiple peaks that are closely spaced, 
; or other assymetries.  this is also reflected in the fact that
; the FWHM is not 2x sigma, as in the n2's that are well approximated
; by a gaussian.
;
;


;============================================
;============================================


xaxistitle= "!6Mass Ratio (M!Dsat!N/M!Dprimary!N)"
xmax= 2.0
xmin= 0.01

;yaxistitle= "!6Starburst FWHM (Gyr)"
;yaxistitle= "!6Starburst Timescale (Gyr)"
yaxistitle= "!7s!6!DSF!N (Gyr)"
ymax= 5.0
ymin= 0.005

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename='minorfwhm.eps', colortable=4

!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, charthick=3.0, xthick=4.0, ythick=4.0, color= 0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /xlog, /ylog, ytickformat='exp_label'


;----------------------------------------


oplot, Mratios, FWHM_n2, psym=5, thick=3.0, color=150, symsize=1.5
idx= where(FWHM_n2 lt 10.0)
Mr= Mratios(idx)
fw= FWHM_n2(idx)
fwwt= FWHM_wt_n2(idx)
oploterror, Mr, fw, fwwt, errcolor=150, color=150, thick= 3.0, errthick= 3.0, psym=3
fitandoverplot_sbtime, mr, fw, fwwt, fNorm= fNorm, alpha=alpha, lcolor=150

        ;x= [0.001,0.01,0.1,1.0,10.0]
        ;y= 0.0 * x + 0.85
        ;oplot, x, y, psym=-3, linestyle= 2, color= 150, thick=2.0

oplot, [0.02], [0.14], psym=5, color=150, symsize=1.5, thick=4.0
xyouts, [0.027], [0.12], '!8n2med!6', color=150, charthick=3.0, /data

;xyouts, [0.023], [0.08], '!8F!6!D0!N=0.36', color=0, charthick=3.0, /data, size=1.2
xyouts, [0.023], [0.08], '!7s!6!D0!N=0.36', color=0, charthick=3.0, /data, size=1.2
xyouts, [0.023], [0.06], '!7a!6=-0.5', color=0, charthick=3.0, /data, size=1.2



oplot, Mratios, FWHM_n0, psym=7, thick=3.0, color=50, symsize=1.5
idx= where(FWHM_n0 lt 10.0)
Mr= Mratios(idx)
fw= FWHM_n0(idx)
fwwt= FWHM_wt_n0(idx)
oploterror, Mr, fw, fwwt, errcolor=50, color=50, thick= 3.0, errthick= 3.0, psym=3
fitandoverplot_sbtime, mr, fw, fwwt, fNorm= fNorm, alpha=alpha, lcolor=50

        ;x= 10^(4.*(findgen(10)/10.0) - 3.0)
        ;y= minor_func_sbtime(x,[0.12,0.7])
        ;oplot, x, y, psym=-3, linestyle= 2, color= 50, thick=2.0
        ;x= [0.001,0.01,0.1,1.0,10.0]
        ;y= 0.0 * x + 0.15
        ;oplot, x, y, psym=-3, linestyle= 2, color= 50, thick=2.0

oplot, [0.02], [0.027], psym=7, color=50, symsize=1.5, thick=4.0
xyouts, [0.027], [0.025], '!8n0med!6', color=50, charthick=3.0, /data

;xyouts, [0.023], [0.017], '!8F!6!D0!N=0.09', color=0, charthick=3.0, /data, size=1.2
xyouts, [0.023], [0.017], '!7s!6!D0!N=0.09', color=0, charthick=3.0, /data, size=1.2
xyouts, [0.023], [0.0125], '!7a!6=-0.7', color=0, charthick=3.0, /data, size=1.2


;oplot, [0.015,0.022], [0.011,0.011], psym=-3, color=0, linestyle=0, thick=2.0
;xyouts, [0.03], [0.0017], 'best fit', color=0, charthick=3.0, /data

;
; print the formula on there
;
;xyouts, 0.22, 0.25, '!6FWHM = !8F!6!D0!N', color=0, charthick=3.0, /normal, size=1.3
;xyouts, 0.44, 0.24, '(', color=0, charthick=3.0, /normal, size=2.6
;xyouts, 0.50, 0.275, 'M!Dsat!N', color=0, charthick=3.0, /normal, size=1.3
;oplot, [0.07,0.18], [0.013,0.013], psym=-3, linestyle=0, thick=4.0
;xyouts, 0.47, 0.23, 'M!Dprimary!N', color=0, charthick=3.0, /normal, size=1.3
;xyouts, 0.60, 0.24, ')', color=0, charthick=3.0, /normal, size=2.6
;xyouts, 0.62, 0.29, '!7a!6', color=0, charthick=3.0, /normal, size=1.2


xyouts, 0.33, 0.88, 'Nuclear Coalescence', color=0, charthick=3.0, size=1.2, /normal


;============================================

device, /close




end









;======================================================================
;
;
;
;    ---------------------------
;    |                         |
;    |                         |
;    |                         |
;    |                         |
;    |                         |
;    |   Tsfr v. mass ratio    |
;    |                         |
;    |                         |
;    |                         |
;    |                         |
;    |                         |
;    ---------------------------
;
;
;
;======================================================================
pro fwhm_firstpass, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "fwhm_firstpass, junk"
   return
endif


;============================================
;============================================



; all of them
Mratios= 1.0/[1.0, 2.3, 5.8, 22.7, $
	1.0, 2.6, 10.0, $
	1.0, 3.9, $
	1.0]


Gstds= ["G3G3b-u1", "G3G2-u3", "G3G1-u3", "G3G0e-u3", $
	"G2G2-u1", "G2G1-u3", "G2G0-u3", $
	"G1G1a-u1", "G1G0-u3", $
	"G0G0a-u1"]


; now take out some of the non-burst models
Mratios= 1.0/[1.0, 2.3, 5.8, $
        1.0, 2.6, 10.0, $
        1.0, 3.9, $
        1.0]


Gstds= ["G3G3b-u1", "G3G2-u3", "G3G1-u3", $
        "G2G2-u1", "G2G1-u3", "G2G0-u3", $
        "G1G1a-u1", "G1G0-u3", $
        "G0G0a-u1"]

T_cut= [2.0, 2.3, 3.0, $
	1.0, 1.1, 2.2, $
	0.6, 1.5, $
	1.3]

;fwhms_n2= grab_fwhm(Gstds, sbts=sbts, sfrms=sfrms)
FWHM_n2= grab_sbtimescale(Gstds, sbtserr=sbtserr, /firstpass, T_cut= T_cut)
FWHM_wt_n2= sbtserr


;----------------------------------------------------------------

; all of them
Gstds= ["G3G3b-u2", "G3G2-u4", "G3G1-u4", "G3G0e-u4", $
	"G2G2-u2", "G2G1-u4", "G2G0-u4", $
	"G1G1a-u2", "G1G0-u4", $
	"G0G0a-u2"]

; on those with SB
Gstds= ["G3G3b-u2", "G3G2-u4", "G3G1-u4", $
	"G2G2-u2", "G2G1-u4", "G2G0-u4", $
	"G1G1a-u2", "G1G0-u4", $
	"G0G0a-u2"]

;fwhms_n0= grab_fwhm(Gstds, sbts=sbts, sfrms=sfrms)
FWHM_n0= grab_sbtimescale(Gstds, sbtserr=sbtserr, /firstpass, T_cut= T_cut)
FWHM_wt_n0= sbtserr


; g3g2 is kind of strange here, and i think a bit
; underestimated in term of the errors
FWHM_n0[1]= 0.1
FWHM_wt_n0[1]= 0.06
FWHM_wt_n0[5]= 0.2



;============================================
;============================================


xaxistitle= "!6Mass Ratio (M!Dsat!N/M!Dprimary!N)"
xmax= 2.0
xmin= 0.01

;yaxistitle= "!6Starburst FWHM (Gyr)"
;yaxistitle= "!6Starburst Timescale (Gyr)"
yaxistitle= "!7s!6!DSF!N (Gyr)"
ymax= 5.0
ymin= 0.005

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename='minorfwhm_fp.eps', colortable=4

!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, charthick=3.0, xthick=4.0, ythick=4.0, color= 0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /xlog, /ylog, ytickformat='exp_label'


;----------------------------------------


oplot, Mratios, FWHM_n2, psym=5, thick=3.0, color=150, symsize=1.5
idx= where(FWHM_n2 lt 10.0)
Mr= Mratios(idx)
fw= FWHM_n2(idx)
fwwt= FWHM_wt_n2(idx)
oploterror, Mr, fw, fwwt, errcolor=150, color=150, thick= 3.0, errthick= 3.0, psym=3
fitandoverplot_sbtime, mr, fw, fwwt, fNorm= fNorm, alpha=alpha, lcolor=150

        ;x= [0.001,0.01,0.1,1.0,10.0]
        ;y= 0.0 * x + 0.85
        ;oplot, x, y, psym=-3, linestyle= 2, color= 150, thick=2.0

oplot, [0.02], [0.14], psym=5, color=150, symsize=1.5, thick=4.0
xyouts, [0.027], [0.12], '!8n2med!6', color=150, charthick=3.0, /data

;xyouts, [0.023], [0.08], '!8F!6!D0!N=0.36', color=0, charthick=3.0, /data, size=1.2
xyouts, [0.023], [0.08], '!7s!6!D0!N=0.57', color=0, charthick=3.0, /data, size=1.2
xyouts, [0.023], [0.06], '!7a!6=-0.3', color=0, charthick=3.0, /data, size=1.2



oplot, Mratios, FWHM_n0, psym=7, thick=3.0, color=50, symsize=1.5
idx= where(FWHM_n0 lt 10.0)
Mr= Mratios(idx)
fw= FWHM_n0(idx)
fwwt= FWHM_wt_n0(idx)
oploterror, Mr, fw, fwwt, errcolor=50, color=50, thick= 3.0, errthick= 3.0, psym=3
fitandoverplot_sbtime, mr, fw, fwwt, fNorm= fNorm, alpha=alpha, lcolor=50

        ;x= 10^(4.*(findgen(10)/10.0) - 3.0)
        ;y= minor_func_sbtime(x,[0.12,0.7])
        ;oplot, x, y, psym=-3, linestyle= 2, color= 50, thick=2.0
        ;x= [0.001,0.01,0.1,1.0,10.0]
        ;y= 0.0 * x + 0.15
        ;oplot, x, y, psym=-3, linestyle= 2, color= 50, thick=2.0

oplot, [0.02], [0.027], psym=7, color=50, symsize=1.5, thick=4.0
xyouts, [0.027], [0.025], '!8n0med!6', color=50, charthick=3.0, /data

;xyouts, [0.023], [0.017], '!8F!6!D0!N=0.09', color=0, charthick=3.0, /data, size=1.2
xyouts, [0.023], [0.017], '!7s!6!D0!N=0.10', color=0, charthick=3.0, /data, size=1.2
xyouts, [0.023], [0.0125], '!7a!6=-0.6', color=0, charthick=3.0, /data, size=1.2


;oplot, [0.015,0.022], [0.011,0.011], psym=-3, color=0, linestyle=0, thick=2.0
;xyouts, [0.03], [0.0017], 'best fit', color=0, charthick=3.0, /data

;
; print the formula on there
;
;xyouts, 0.22, 0.25, '!6FWHM = !8F!6!D0!N', color=0, charthick=3.0, /normal, size=1.3
;xyouts, 0.44, 0.24, '(', color=0, charthick=3.0, /normal, size=2.6
;xyouts, 0.50, 0.275, 'M!Dsat!N', color=0, charthick=3.0, /normal, size=1.3
;oplot, [0.07,0.18], [0.013,0.013], psym=-3, linestyle=0, thick=4.0
;xyouts, 0.47, 0.23, 'M!Dprimary!N', color=0, charthick=3.0, /normal, size=1.3
;xyouts, 0.60, 0.24, ')', color=0, charthick=3.0, /normal, size=2.6
;xyouts, 0.62, 0.29, '!7a!6', color=0, charthick=3.0, /normal, size=1.2


xyouts, 0.33, 0.88, 'First Passage', color=0, charthick=3.0, size=1.2, /normal


;============================================

device, /close




end



