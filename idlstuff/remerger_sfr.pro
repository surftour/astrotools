pro sfr_multi, junk, $
		filename=filename, $
		cumulative=cumulative, $
		gasmass=gasmass, $
		ynotlog=ynotlog, $
		h=h


if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_multi, junk, filename=filename, /h, "
   print, "           /cumulative, /gasmass, /ynotlog"
   print, "  "
   print, "  "
   print, "  WARNING: cumulative and gasmass do not currently work!  "
   print, "  "
   print, "  default filename: sfr.eps"
   print, "  "
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfr.eps'
if keyword_set(gasmass) then filename='gmass.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
if keyword_set(cumulative) then yaxistitle="!6Mass (10!e10!n h!E-1!NM!D!9n!6!N)"
if keyword_set(gasmass) then yaxistitle="!6Gas Mass (10!e10!n h!E-1!NM!D!9n!6!N)"

;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

;xmax = 14.0
;xmax = 12.0
;xmax = 10.0
;xmax = 7.5
;xmax = 6.0
xmax = 5.0
;xmax = 4.25
;xmax = 4.0
;xmax = 3.5
;xmax = 3.0
;xmax = 2.8
;xmax = 2.5
;xmax = 2.4
;xmax = 2.0
;xmax = 1.5
;xmax = 1.4
;xmax = 1.3
;xmax = 1.0
;xmax = 0.5
xmin = 0

;ymax= 8e+5
;ymax= 8e+4
;ymax= 2e+4
;ymax = 8000
;ymax = 4000
;ymax = 2500
;ymax = 1750
;ymax = 1500
;ymax = 1000
;ymax = 800
;ymax = 600
ymax = 400
;ymax = 300
;ymax = 250.0
;ymax = 202
;ymax = 180
;ymax = 120
;ymax = 100
;ymax = 90
;ymax = 80
;ymax = 60
;ymax = 50.0
;ymax = 40.0
;ymax = 30.0
;ymax = 20.0
;ymax = 16.0
;ymax = 15
;ymax = 13
;ymax = 12
;ymax = 10.0
;ymax = 8.0
;ymax = 6.0
;ymax = 5.0
;ymax = 3.5
;ymax = 2.0
;ymax = 1.5
;ymax = 1.3
;ymax = 1.2
;ymax = 0.75
;ymax = 0.6
;ymax = 0.5
;ymin= 8.0
;ymin= 1.0
;ymin = 0  & ynotlog= 1
;ymin = 0.8
;ymin = 0.4
;ymin = 0.1
;ymin= 0.07
ymin= 0.01
;ymin= 0.001
;ymin= 0.0001
;ymin= 0.00001

if keyword_set(ynotlog) then ymin = 0 

; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	;h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;---------------------------

!p.position= [0.16, 0.14, 0.98, 0.98]

if keyword_set(ynotlog) then begin
    plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	;ytickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata
endif else begin
    plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	ytickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata
endelse




;-------------------------------------------
;   Load the runs to display
;-------------------------------------------


process_one_sfr, "remergers/b3eb3e2", lcolor=50, lthick= 2.0, msg='b3eb3e2', x0= 0.22, y0= 0.88
process_one_sfr, "remergers/b3ev2", lcolor=150, lthick= 2.0, msg='b3ev2', x0= 0.22, y0= 0.84
process_one_sfr, "remergers/b3ev3", lcolor=100, lthick= 2.0, msg='b3ev3', x0= 0.22, y0= 0.80
process_one_sfr, "remergers/b3ev4", lcolor=200, lthick= 2.0, msg='b3ev4', x0= 0.22, y0= 0.76
process_one_sfr, "remergers/b3esph", lcolor=0, lthick= 2.0, msg='b3esph', x0= 0.22, y0= 0.76





;--------------------------------------
;--------------------------------------

device, /close


end





;
;
;
;
;=====================================================================
;
;
;
;
;
;
;
;
;
;
;


pro sfr_one, frun, $
		xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, $
		lcolor=lcolor, lstyle=lstyle, lthick=lthick, $
		ctab=ctab, msg=msg, y0=y0, x0=x0, lpsym=lpsym, $
		filename=filename, $
		cumulative=cumulative, $
		gasmass=gasmass, $
		ynotlog=ynotlog, $
		show_specific_times=show_specific_times, $
		h=h


if not keyword_set(frun) then begin
   print, "  "
   print, "sfr_one, frun, filename=filename, /h, "
   print, "           /cumulative, /gasmass, /ynotlog"
   print, "  "
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
if keyword_set(cumulative) then yaxistitle="!6Mass (10!e10!n h!E-1!NM!D!9n!6!N)"
if keyword_set(gasmass) then yaxistitle="!6Gas Mass (10!e10!n h!E-1!NM!D!9n!6!N)"

;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7
;if keyword_set(h) then xaxistitle = "Time (Gyr)"

if not keyword_set(xmax) then xmax = 2.4
if not keyword_set(xmin) then xmin = 0

if not keyword_set(ymax) then ymax = 250.0
if not keyword_set(ymin) then ymin = 0.1

if keyword_set(ynotlog) then ymin = 0 


if not keyword_set(lcolor) then lcolor= 150
if not keyword_set(lthick) then lthick= 2.0
if not keyword_set(lstyle) then lstyle= 0

;if not keyword_set(msg) then msg=' '

if not keyword_set(x0) then x0= 0.70
if not keyword_set(y0) then y0= 0.85

;---------------------------

!p.position= [0.16, 0.14, 0.98, 0.98]

if keyword_set(ynotlog) then begin
    plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	;ytickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata
endif else begin
    plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	ytickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata
endelse




;--------------------------------------
;--------------------------------------



process_one_sfr, frun, lcolor=lcolor, lthick= lthick, ctab=ctab, $
			 msg=msg, x0= x0, y0= y0, lstyle= lstyle, h=h, $
			 show_specific_times=show_specific_times





;--------------------------------------
;--------------------------------------

device, /close


end







    




    


;===============================================================================
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

pro sfr, frun, filename=filename, quiet=quiet, $
		show_specific_time=show_specific_time
;pro junk

if not keyword_set(frun) then begin
;if not keyword_set(junk) then begin
   print, "  "
   print, "sfr, frun, filename=filename, quiet=quiet, $"
   print, "          show_specific_time=show_specific_time"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename=frun+"/sfr.eps"

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
xaxistitle = "!6Time (Gyr/h)"


open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac

xmax= max(sfrtime)
xmin= min(sfrtime)

ymax= max(sfrsfr) * 1.2
ymin= min(sfrsfr(where(sfrsfr gt 0.0))) * 0.2
;ymin= 0.002

; physical units
if keyword_set(h) then begin
        xaxistitle = "Time (Gyr)"
        h = fload_cosmology('h')

	sfrtime = sfrtime / h
        xmax = xmax / h
        xmin = xmin / h
endif


;---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]
;!p.font= 0

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, ytickforma='exp_label', $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;-------------------------------------------\

; plot sfr history

oplot, sfrtime, sfrsfr, psym=-3, color= 0, linestyle= 0, thick= 5.0



;-------------------------------------------

t_merg= fload_bh_mergertime(frun)
sfr_max_raw= max(sfrsfr)
t_sfr_max= sfrtime(where(sfrsfr eq sfr_max_raw))


;do_meansf= 0
do_meansf= 1
if do_meansf eq 1 then begin

	; initial period
	idx= where(sfrtime lt t_merg-0.2)
	meansfr_initial= mean(sfrsfr(idx))
	x= [min(sfrtime), t_merg-0.2]
	y= [meansfr_initial, meansfr_initial]
	oplot, x, y, psym=-3, linestyle=1, color= 150, thick=8.0

	; starburst
	idx= where((sfrtime ge t_merg-0.2) and (sfrtime le t_merg+0.2))
	meansfr_sb= mean(sfrsfr(idx))
	x= [t_merg-0.2, t_merg+0.2]
	y= [meansfr_sb, meansfr_sb]
	oplot, x, y, psym=-3, linestyle=1, color= 150, thick=8.0

	; final period
	idx= where(sfrtime gt t_merg+0.2)
	if idx(0) ne -1 then begin
		meansfr_final= mean(sfrsfr(idx))
		x= [t_merg+0.2,max(sfrtime)]
		y= [meansfr_final, meansfr_final]
		oplot, x, y, psym=-3, linestyle=1, color= 150, thick=8.0
	endif else begin
		meansfr_final= 1e-10
	endelse

endif


;do_sb_fit= 0
do_sb_fit= 1        ; first passage starburst
if do_sb_fit eq 1 then begin
        idx= where((sfrtime gt 0.1) and (sfrtime lt t_merg-0.2))
        time= sfrtime(idx)
        sfr= sfrsfr(idx)

	sfr_max_first_peak= max(sfr)
	t_first_peak= time(where(sfr eq sfr_max_first_peak))

        fitandoverplot_gaussian, time, sfr, $
                        SFR_max= SFR_max_fp, width=width_fp, SBtime=SBtime_fp, quiet=quiet
endif


;do_sb_fit= 0
do_sb_fit= 1        ; final merger starburst
if do_sb_fit eq 1 then begin
        ; ----------------
        ; use full data
        idx= where((sfrtime gt t_merg-0.2) and (sfrtime lt t_merg+0.2))
        time= sfrtime(idx)
        sfr= sfrsfr(idx)
        ; ----------------
        ; use smoothed data
        ;idx= where((timem gt t_merg-0.8) and (timem lt t_merg+0.8))
        ;time= timem(idx)
        ;sfr= sfrm(idx)
        ; ----------------

	sfr_max_sb= max(sfr)
	t_max_sb= time(where(sfr eq sfr_max_sb))

        fitandoverplot_gaussian, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime, quiet=quiet
endif

do_sb_fit= 0
;do_sb_fit= 1
if do_sb_fit eq 1 then begin
        idx= where((sfrtime gt t_merg-0.2) and (sfrtime lt t_merg+0.2))
        time= sfrtime(idx)
        sfr= sfrsfr(idx)
        fitandoverplot_gaussian_pluszeropt, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime
endif


;do_expdecay_fit= 0
do_expdecay_fit= 1
if do_expdecay_fit eq 1 then begin
        idx= where(sfrtime gt t_merg)
        ;idx= where((sfrtime gt SBtime) and (sfrtime lt SBtime+1.3))
        time= sfrtime(idx)
        time= time- time(0)
        sfr= sfrsfr(idx)
        fitandoverplot_expdecay, time, sfr, SBtime=t_merg, tau=tau, decay_const=decay_const
endif

do_fwhm_fit= 0
;do_fwhm_fit= 1
if do_fwhm_fit eq 1 then begin
        ; ----------------
        ; use full data
        idx= where((sfrtime gt t_merg-0.5) and (sfrtime lt t_merg+0.5))
        time= sfrtime
        sfr= sfrsfr
        ; ----------------
        ; use smoothed data
        ;time= timem
        ;sfr= sfrm
        ; ----------------
        findandoverplot_fwhm, time, sfr, $
                        ;SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, quiet=quiet
                        ;SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, /quiet
                        SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime
endif



;--------------------------------------
;--------------------------------------



openw, 1, frun+'/sfrfits.txt', ERROR=err

printf, 1, "#   "
printf, 1, "#    Fits to the Star Formation History  "
printf, 1, "#   (all quantities are in Gadget units) "
printf, 1, "#   "
printf, 1, "#   "

printf, 1, "t_merg        ", t_merg
printf, 1, "sfr_max_raw   ", sfr_max_raw
printf, 1, "t_sfr_max     ", t_sfr_max
printf, 1, "sfr_max_fp    ", sfr_max_first_peak
printf, 1, "t_fp          ", t_first_peak
printf, 1, "sfr_max_sb    ", sfr_max_sb
printf, 1, "t_sb          ", t_max_sb
printf, 1, " "

printf, 1, "Mean SFRs (0 -> t_merg - 0.2 Myr/h -> t_merg + 0.2 Myr/h -> end)"
printf, 1, "meansfr_init  ", meansfr_initial
printf, 1, "meansfr_sb    ", meansfr_sb
printf, 1, "meansfr_final ", meansfr_final
printf, 1, " "

printf, 1, "Gaussian Fit to First Passage Starburst"
printf, 1, "SFR_max_fp    ", SFR_max_fp
printf, 1, "SBtime_fp     ", SBtime_fp
printf, 1, "width_fp      ", width_fp
printf, 1, " "

printf, 1, "Gaussian Fit to Starburst at t_merg"
printf, 1, "SFR_max       ", SFR_max
printf, 1, "SBtime        ", SBtime
printf, 1, "width         ", width
printf, 1, " "

printf, 1, "Exp. Fit to Post-SB SF"
printf, 1, "tau           ", tau
printf, 1, "decay_const   ", decay_const
printf, 1, " "

close, 1



;--------------------------------------
;--------------------------------------

device, /close


end















;===============================================================================


pro process_one_sfr, frun, lcolor=lcolor, lstyle=lstyle, lthick=lthick, $
				ctab=ctab, msg=msg, y0=y0, x0=x0, lpsym=lpsym, $
				cumulative=cumulative, $
				gasmass=gasmass, $
				gasfraction=gasfraction, $
				normalize=normalize, h=h, $
				multiby=multiby, $
				sfrtime=sfrtime, sfrsfr=sfrsfr, $
				show_specific_times=show_specific_times, $
				raw=raw


	if not keyword_set(lthick) then lthick= 1.0
	if not keyword_set(lstyle) then lstyle= 0
	if not keyword_set(lcolor) then lcolor= 0
	if not keyword_set(lpsym) then lpsym= -3


	if keyword_set(ctab) then begin
                        loadct, ctab
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0
	endif


	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	if keyword_set(raw) then begin
		open_raw_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac
	endif else begin
		open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac
	endelse


	; physical units
	;if i eq 1 then sfrtime = sfrtime/(0.7)    ; did this for comparing h and noh
	if keyword_set(h) then begin
		print, "correcting time by h= ", h
		sfrtime = sfrtime / h
	endif

	if keyword_set(normalize) then sfrsfr=sfrsfr/sfrsfr(0)

	if keyword_set(multiby) then sfrsfr=sfrsfr * multiby


	    lthick= 4.0 * lthick

	if keyword_set(cumulative) then begin
		oplot, sfrtime, sfrmfs, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick
		goto, moveon
	endif


	if keyword_set(gasmass) then begin
		oplot, sfrtime, sfrgasmass, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick
		goto, moveon
	endif


	if keyword_set(gasfraction) then begin
		gasf= gasfraction * gasfrac
		oplot, sfrtime, gasf, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick
		goto, moveon
	endif


	; default option
	oplot, sfrtime, sfrsfr, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick


moveon:


if not keyword_set(y0) then y0= 0.92
if not keyword_set(x0) then x0= 0.75

if keyword_set(msg) then begin
    xyouts, x0, y0, msg, /normal, charthick=3.0, size=1.33, color= lcolor
endif else begin
    xyouts, x0, y0, frun, /normal, charthick=3.0, size=1.33, color= lcolor
endelse



do_sb_fit= 0
;do_sb_fit= 1
if do_sb_fit eq 1 then begin
	t_merg= fload_bh_mergertime('/odyssey/tcox/'+frun)
	idx= where((sfrtime gt t_merg-0.3) and (sfrtime lt t_merg+0.3))
	time= sfrtime(idx)
	sfr= sfrsfr(idx)
	fitandoverplot_gaussian, time, sfr, $
			SFR_max= SFR_max, width=width, SBtime=SBtime
endif

do_expdecay_fit= 0
;do_expdecay_fit= 1
if do_expdecay_fit eq 1 then begin
	;idx= where(sfrtime gt SBtime)
	idx= where((sfrtime gt SBtime) and (sfrtime lt SBtime+1.3))
	time= sfrtime(idx)
	time= time- time(0)
	sfr= sfrsfr(idx)
	fitandoverplot_expdecay, time, sfr, SBtime=SBtime
endif


if keyword_set(show_specific_times) eq 1 then begin

	for ii= 0, n_elements(show_specific_times)-1 do begin
		;specific_time= 1.205
		specific_time= show_specific_times[ii]
		idx=where(sfrtime ge specific_time)
		oplot, [sfrtime[idx[0]]], [sfrsfr[idx[0]]], psym= 4, color= 150, thick=7.0, symsize=2.0
	endfor
endif

end







;======================================================================






pro open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac

	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	    ; get sfr data (from file)
	    ; -------------------------
	    loopnum= 0
	    repeat begin
	    ;sfrfile= '/data/tcox/sfr/'+fload_getid(fruns[i])+'.sfr'        ; twopiter
	    ;sfrfile= '/home/tcox/data/sfr/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/raid4/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/raid2/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/data6/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard

                ;if strmid(frun,1,4) eq 'raid' then begin
		;print, strmid(frun,0,7)
                if strmid(frun,0,7) eq '/n/data' then begin
                        spawn, '/bin/ls '+frun+'/*.sfr ',result
                endif else begin
			case loopnum of
			  0: datadir='/n/home/tcox/data'
			  1: datadir='/n/circelfs/hernquist_lab'
			  2: datadir='/n/home/tcox'
			  3: datadir='/n/home'
			  4: datadir='/n/scratch/hernquist_lab/tcox'
			  5: datadir='/data7'
			  6: datadir=''
			  else: break
			endcase

			;sfrfile= datadir+'/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'

			nfrun= fload_getid(frun)
			;spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/*.sfr ',result
			spawn, '/bin/ls '+datadir+"/"+nfrun+'/*.sfr ',result
		endelse

		sfrfile=strcompress(result[0],/remove_all)



		get_lun, unit
		openr, unit, sfrfile, ERROR=err
		close, unit
		free_lun, unit

		if (err NE 0) then begin
		   ERR= 0
		   foundit= 0
		endif else begin
		   foundit= 1
		endelse

		loopnum= loopnum+1

	    endrep until (foundit eq 1)

	    if foundit eq 0 then begin
		msg= 'Problem: sfr file no found  file=[data]/tcox/'+fload_getid(frun)+'/'+fload_getid(frun)+'.sfr'
	        print, "  "
	        print, msg
	        print, "  "
		return
	    endif

	    ; read to open
	    print, "opening: ",sfrfile
	    sfrdata= read_ascii(sfrfile)
	    sfrtime= sfrdata.field1[0,*]
	    sfrsfr= sfrdata.field1[1,*]
	    sfrmfs= sfrdata.field1[3,*]    ; amount of new stars
	    gasfrac= sfrdata.field1[4,*]
	    n_cols= n_elements(sfrdata.field1[*,0])
	    n_rows= n_elements(sfrdata.field1[0,*])
	    finalnewstarmass= sfrmfs[n_rows-1]
	    ; gas mass
	    if n_cols ge 6 then sfrgasmass= sfrdata.field1[6,*] else sfrgasmass=[20.0,20.0-finalnewstarmass]
	    sfrmsg= ''


            t1per= 0
            idx=where(sfrtime ge 1.0)
            if idx(0) ne -1 and n_cols gt 5 then begin
                gm= sfrgasmass[0]-sfrgasmass(idx(0))
                t1per= 100.0*gm/sfrgasmass[0]
            endif

            t4per= 0
            idx=where(sfrtime ge 4.0)
            if idx(0) ne -1 and n_cols gt 5 then begin
                gm= sfrgasmass[0]-sfrgasmass(idx(0))
                t4per= 100.0*gm/sfrgasmass[0]
            endif

            maxsfr= max(sfrsfr)
            idx= where(sfrsfr eq maxsfr)
            maxsfrtime= sfrtime(idx)

            sfraftertmerg= 0.0
            smaftertmerg= 0.0
            taftertmerg= 0.0
            tmerger= 0.0
            if n_elements(tmerg) gt 1 then begin
                idx= where(sfrtime ge tmerg[i]+0.2)
                if idx(0) eq -1 then begin
                        sfraftertmerg= sfrsfr[n_rows-1]
                        smaftertmerg= sfrmfs[n_rows-1]
                        tmerger= sfrtime[n_rows-1]
                        taftertmerg= tmerger
                endif else begin
                        sfraftertmerg= sfrsfr[idx(0)]
                        smaftertmerg= sfrmfs[idx(0)]
                        tmerger= tmerg[i]
                        taftertmerg= sfrtime[idx(0)]
                endelse
            endif

            print, "-------------------------------------"
            print, "maximum sfr= ", maxsfr, " occurs at ", maxsfrtime
            print, "average sfr= ", total(sfrsfr)/n_rows
            print, "  "
            print, "      tmerg= ",tmerger,'  and 200 Myr later= ',taftertmerg
            print, "        sfr= ", sfraftertmerg,'  200 Myr after t_merger'
            print, "         sm= ", smaftertmerg,'  200 Myr after t_merger'
            n_mfs= n_elements(sfrgasmass)
            print, "-------------------------------------"
            print, "original gas mass   =", sfrgasmass[0]
            print, "remnant gas mass    =", sfrgasmass[n_mfs-1]
            print, "gas consumed        =", sfrgasmass[0]-sfrgasmass[n_mfs-1]
            print,"                      ", 100.0*(sfrgasmass[0]-sfrgasmass[n_mfs-1])/sfrgasmass[0],'  %'
            print,"                      ", t1per,' % after 1 Gyr'
            print,"                      ", t4per,' % after 4 Gyr'
            print, "-------------------------------------"


end





;======================================================================






pro open_raw_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac

	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	    ; get sfr data (from file)
	    ; -------------------------
	    loopnum= 0
	    repeat begin
	    ;sfrfile= '/data/tcox/sfr/'+fload_getid(fruns[i])+'.sfr'        ; twopiter
	    ;sfrfile= '/home/tcox/data/sfr/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/raid4/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/raid2/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/data6/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard

                if strmid(frun,1,4) eq 'raid' then begin
                        spawn, '/bin/ls '+frun+'/sfr.txt ',result
                endif else begin
			case loopnum of
                          0: datadir='/n/home/tcox/data'
                          1: datadir='/n/circelfs/hernquist_lab'
                          2: datadir='/n/home/tcox'
                          3: datadir='/n/home'
                          4: datadir='/n/scratch/hernquist_lab/tcox'
                          5: datadir='/data7'
                          6: datadir=''
                          else: break
			endcase

			;sfrfile= datadir+'/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'

			nfrun= fload_getid(frun)
			spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/sfr.txt ',result
		endelse

		sfrfile=strcompress(result[0],/remove_all)



		get_lun, unit
		openr, unit, sfrfile, ERROR=err
		close, unit
		free_lun, unit

		if (err NE 0) then begin
		   ERR= 0
		   foundit= 0
		endif else begin
		   foundit= 1
		endelse

		loopnum= loopnum+1

	    endrep until (foundit eq 1)

	    if foundit eq 0 then begin
		msg= 'Problem: sfr file no found  file=[data]/tcox/'+fload_getid(frun)+'/sfr.txt'
	        print, "  "
	        print, msg
	        print, "  "
		return
	    endif

	    ; read to open
	    print, "opening: ",sfrfile
	    sfrdata= read_ascii(sfrfile)
	    sfrtime= sfrdata.field1[0,*]
	    sfrmfs= sfrdata.field1[1,*]
	    sfrsfr_gadu= sfrdata.field1[2,*]    ; amount of new stars
	    sfrsfr= sfrdata.field1[3,*]
	    sfrmfs_spawned= sfrdata.field1[4,*]
	    n_cols= n_elements(sfrdata.field1[*,0])
	    n_rows= n_elements(sfrdata.field1[0,*])

            maxsfr= max(sfrsfr)
            idx= where(sfrsfr eq maxsfr)
            maxsfrtime= sfrtime(idx)

            print, "-------------------------------------"
            print, "maximum sfr= ", maxsfr, " occurs at ", maxsfrtime
            print, "average sfr= ", total(sfrsfr)/n_rows
            print, "-------------------------------------"

	    N_sample= 5
            print, "-------------------------------------"
            print, "-------------------------------------"
	    print, "   RESAMPLING N_sample= ", N_sample
	    sfrsfr= transpose(sfrsfr)
	    sfrtime= transpose(sfrtime)
	    n= n_elements(sfrsfr)
	    n_sub= n/N_sample
	    n_extra= n-N_sample*n_sub
	    for i=0,(N_sample-n_extra-1) do begin
		sfrsfr= [sfrsfr, sfrsfr[n-1]]
		sfrtime= [sfrtime, sfrtime[n-1]]
	    endfor
	    sfrsfr= rebin(sfrsfr, n_sub+1)
	    sfrtime= rebin(sfrtime, n_sub+1)
            print, "-------------------------------------"
            print, "-------------------------------------"

end





;--------------------------------------------------------------------



pro read_sfrfits_file, frun, $
		t_merg, sfr_max, t_sfr_max, $
		sfr_max_fp, t_sfr_max_fp, $
		sfr_max_sb, t_sfr_max_sb, $
		meansfr_initial, meansfr_sb, meansfr_final, $
		gSFR_max_fp, gSBtime_fp, gwidth_fp, $
		gSFR_max, gSBtime, gwidth, $
		tau, decay_const


;sfrfitsfile= '/raid4/tcox/'+frun+'/sfrfits.txt'
sfrfitsfile= frun+'/sfrfits.txt'

openr, 1, sfrfitsfile
junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_merg= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & sfr_max= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_sfr_max= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & sfr_max_fp= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_sfr_max_fp= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & sfr_max_sb= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_sfr_max_sb= float(tempjunk(1))
readf, 1, junk

readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & meansfr_initital= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & meansfr_sb= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & meansfr_final= float(tempjunk(1))
readf, 1, junk

readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gSFR_max_fp= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gSBtime_fp= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gwidth_fp= float(tempjunk(1))
readf, 1, junk

readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gSFR_max= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gSBtime= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gwidth= float(tempjunk(1))
readf, 1, junk

readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & tau= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & decay_const= float(tempjunk(1))
readf, 1, junk

close, 1


end










;==============================================================================
;==============================================================================


pro sfr_test, junk, $
		filename=filename, $
		cumulative=cumulative, $
		gasmass=gasmass, $
		ynotlog=ynotlog, $
		h=h


if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_test, junk, filename=filename, /h, "
   print, "           /cumulative, /gasmass, /ynotlog"
   print, "  "
   print, "  "
   print, "  WARNING: cumulative and gasmass do not currently work!  "
   print, "  "
   print, "  default filename: sfrtest.eps"
   print, "  "
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfrtest.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
if keyword_set(cumulative) then yaxistitle="!6Mass (10!e10!n h!E-1!NM!D!9n!6!N)"
if keyword_set(gasmass) then yaxistitle="!6Gas Mass (10!e10!n h!E-1!NM!D!9n!6!N)"

xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
;xaxistitle = "!6Time (Gyr)"  & h= 0.7

;xmax = 14.0
;xmax = 12.0
;xmax = 10.0
;xmax = 7.5
;xmax = 6.0
;xmax = 5.0
;xmax = 4.25
;xmax = 4.0
;xmax = 3.5
;xmax = 3.0
;xmax = 2.8
;xmax = 2.4
;xmax = 2.0
xmax= 1.8
;xmax = 1.5
;xmax = 1.3
;xmax = 1.0
;xmax= 0.8
;xmax= 0.75
;xmax = 0.5
xmin = 0

;ymax= 8e+5
;ymax= 8e+4
;ymax = 8000
;ymax = 2500
;ymax = 1750
;ymax = 1500
;ymax = 1000
;ymax = 600
;ymax = 400
;ymax = 300
;ymax = 250.0
;ymax = 180
ymax = 120
;ymax = 100
;ymax = 90
;ymax = 80
;ymax = 60
;ymax = 50.0
;ymax = 40.0
;ymax = 20.0
;ymax = 16.0
;ymax = 15
;ymax = 13
;ymax = 12
;ymax = 10.0
;ymax = 8.0
;ymax = 6.0
;ymax = 5.0
;ymax = 3.5
;ymax = 2.0
;ymax = 1.5
;ymax = 1.3
;ymax = 1.2
;ymax = 0.75
;ymax = 0.6
;ymax = 0.5
;ymin= 8.0
;ymin= 1.0
ymin = 0  & ynotlog= 1
;ymin = 0.1
;ymin= 0.07
;ymin= 0.01
;ymin= 0.001
;ymin= 0.0001
;ymin= 0.00001

if keyword_set(ynotlog) then ymin = 0 

; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	;h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;---------------------------

!p.position= [0.16, 0.14, 0.98, 0.98]

if keyword_set(ynotlog) then begin
    plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	;ytickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata
endif else begin
    plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	ytickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata
endelse



;frun_orig= "ssft/iso2a"
frun_orig= junk
;snapnum= 30



;-------------------------------------------
;   Load the sfr.txt file
;-------------------------------------------
frun= frun_orig
;process_one_sfr, frun, lcolor=50, lthick= 2.0, msg='from sfr.txt file', x0= 0.60, y0= 0.83, /raw




;-------------------------------------------
;   Load the *.sfr file
;-------------------------------------------
frun= frun_orig
process_one_sfr, frun, lcolor=150, lthick= 2.0, msg='average SFR - from gas density', x0= 0.25, y0= 0.88, $
			sfrtime=sfrtime, sfrsfr=sfrsfr




;-------------------------------------------
;   Infer SFR from final snapshot 
;-------------------------------------------


spawn, "/bin/ls "+frun+"/snapshot_* | wc ",result
snapnum=long(result[0])-1

ok=fload_snapshot_bh(frun,snapnum)
id=fload_newstars_id(1)
time=fload_time(1)
formationtime=fload_newstars_age(1)
m=fload_newstars_mass(1)
print, min(m), max(m), mean(m)

timebinsize= 0.01     ; in Gadget units = 10 Myr
oplotit= 3
histlevels= (xmax-xmin)/timebinsize

; normalization ??
scaling= mean(m) * 1.0d10 / (timebinsize * 1.0d9)     ; conver particle mass to m_solar, and then we're uring 10 Myr bins
normalization= 1.0/scaling

idx= where(formationtime le xmax)
formationtime= formationtime(idx)

;stop
; now the histogram
temp= process_histogram(formationtime, xmax=xmax, xmin=xmin, levels=histlevels, $
				 oplotit=oplotit, mannorm=normalization, bins=bins)
;print, min(temp), max(temp)
xyouts, 0.25, 0.75, "new star formation times", size=1.5, color=0, /normal, charthick=3.0


;--------------------------------------
;--------------------------------------

device, /close






;=========================================
;
; now, display differences
;

	diffsfr=temp*0.0
	for i=0,long(n_elements(diffsfr)-2) do begin
	   idx=where(sfrtime gt (bins[i]-1.0e-6))
	   if idx(0) ne -1 then diffsfr[i]= temp[i]-sfrsfr(idx(0)) else print, "bad i= ", i
	endfor

	diffrange=1.0*long(max(abs(diffsfr))*1.2)
	momentdiff= moment(diffsfr)
	p= [1.0, sqrt(momentdiff(1)), momentdiff(0)]
	

	wheredot= strlen(filename)-4
	filename= strmid(filename,0,wheredot)+'_diff.eps'
	initialize_plotinfo, 1
	setup_plot_stuff, 'ps', filename=filename, colortable= 4

	plot, [1.0],[1.0], psym=-3, xrange=[-diffrange,diffrange], yrange=[0.0,1.0], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        charthick=3.0, xtitle='!7D!3 (M!D!9n!6!N yr!E-1!N)', ytitle=' ', /nodata

	temp=process_histogram(diffsfr, xmax=diffrange, xmin=-diffrange, levels=20, oplotit=5)

	x= diffrange * (2.0 * findgen(101)/100.0 - 1.0)
	oplot, x, func_gaussian(x,p), color= 150, thick=2.0

	device, /close

;=========================================


end

    









;==============================================================================
;==============================================================================






pro sfr_4, junk, $
		filename=filename, $
		cumulative=cumulative, $
		gasmass=gasmass, $
		h=h


if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_4, junk, filename=filename, /h"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=14.0, newysize=14.0
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
if keyword_set(cumulative) then yaxistitle="!6Mass (10!e10!n h!E-1!NM!D!9n!6!N)"
if keyword_set(gasmass) then yaxistitle="!6Gas Mass (10!e10!n h!E-1!NM!D!9n!6!N)"

xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
;xaxistitle = "!6Time (Gyr)"  & h= 0.7

;xmax = 14.0
;xmax = 12.0
;xmax = 10.0
;xmax = 7.5
;xmax = 6.0
;xmax = 5.0
;xmax = 4.25
;xmax = 4.0
xmax = 3.0
;xmax = 2.8
;xmax = 2.0
;xmax = 1.5
;xmax = 1.3
;xmax = 1.0
;xmax = 0.5
xmin = 0

;ymax = 8000
;ymax = 2500
;ymax = 1750
;ymax = 1500
;ymax = 400
ymax = 300
;ymax = 250.0
;ymax = 180
;ymax = 120
;ymax = 100
;ymax = 90
;ymax = 80
;ymax = 60
;ymax = 50.0
;ymax = 20.0
;ymax = 16.0
;ymax = 15
;ymax = 13
;ymax = 12
;ymax = 10.0
;ymax = 8.0
;ymax = 6.0
;ymax = 5.0
;ymax = 3.5
;ymax = 2.0
;ymax = 1.5
;ymax = 1.3
;ymax = 1.2
;ymax = 0.75
;ymax = 0.6
;ymax = 0.5
;ymin = 0 
;ymin= 0.01
ymin= 0.002
;ymin= 0.001
;ymin= 0.0001
;ymin= 0.00001


; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;-----------------------------------------
;
;   Load the panels
;
;-----------------------------------------
;
;     |---------------------|
;     |          |          |
;     |          |          |
;     |    1     |    2     |
;     |          |          |
;     |          |          |
;     |---------------------|
;     |          |          |
;     |          |          |
;     |    3     |    4     |
;     |          |          |
;     |          |          |
;     |---------------------|
;
;

x0= 0.14
xs= 0.5*(0.98-x0)
x1= x0+xs
x2= x0+xs+xs

y0= 0.10
ys= 0.5*(0.98-y0)
y1= y0+ys
y2= y0+ys+ys



;  1
; ---
!p.position= [x0,y1,x1,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, $
	xtickformat='(a1)', ytitle=yaxistitle, ytickformat='exp_label'

fload_newcolortable, 4
process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0
process_one_sfr, "sbw/sb6", lcolor=220, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb2", lcolor=180, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb1", lcolor=140, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb3", lcolor=100, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb4", lcolor=60, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb5", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0

xyouts, x1-0.15, y2-0.06, '!7g!6=0.005', /normal, charthick=3, size=1.33, color=0



;  2
; ---
!p.position= [x1,y1,x2,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

fload_newcolortable, 4
process_one_sfr, "vc3vc3e_no", lcolor= 200, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb10", lcolor=200, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0

xyouts, x2-0.15, y2-0.06, '!7g!6=0.5', /normal, charthick=3, size=1.33, color=0



;  3
; ---
!p.position= [x0,y0,x1,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytitle=yaxistitle, ytickformat='exp_label'

fload_newcolortable, 4
process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0
process_one_sfr, "sbw/sb8", lcolor=200, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb12", lcolor=140, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb13", lcolor=80, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb14", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0

xyouts, x1-0.15, y1-0.06, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0




;  4
; ---
!p.position= [x1,y0,x2,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytickformat='(a1)'

fload_newcolortable, 4
process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0
process_one_sfr, "sbw/sb11", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0

xyouts, x2-0.15, y1-0.06, '!7g!6=5.0', /normal, charthick=3, size=1.33, color=0




;--------------------------------------
device, /close


end













;======================================================================








pro sfr_at_snaptimes, frun, h=h


if not keyword_set(frun) then begin
   print, "  "
   print, "sfr_at_snaptimes, frun, /h"
   print, "  "
   print, "  "
   return
endif


;-------------------------------------------
;-------------------------------------------


; determine the number of
; snapshots in frun directory
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snapshot* | wc ",result
nsnaps=long(result[0])




; ---------------------------------------------------


time= fltarr(nsnaps)

sfr_inst= fltarr(nsnaps)
sfr_avg= fltarr(nsnaps)
sfr_snap= fltarr(nsnaps)


ymin = 0 


; physical units
if keyword_set(h) then begin
	h = fload_cosmology('h')
endif


; ---------------------------------------------------


	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac



	    ; physical units
	    if keyword_set(h) then sfrtime = sfrtime / h
	    ;if i eq 1 then sfrtime = sfrtime/(0.7)    ; did this for comparing h and noh


; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"


        ; open snapshot
        ;ok=fload_snapshot(frun,snapnum)
        ;ok=fload_snapshot_bh(frun,i,/nopot_in_snap)
        ok=fload_snapshot_bh(frun,i,/skip_center)


	sfr_snap[i]= total(fload_gas_sfr(1))

        ; what time is it?
        time[i]= fload_time(1)

	idx_gtr_snaptime= where(sfrtime ge (time[i]-0.0001))
	idx_curr_snaptime= idx_gtr_snaptime[0]


	; instantaneous sfr
	sfr_inst[i]= -1
	if idx_curr_snaptime(0) ne -1 then sfr_inst[i]= sfrsfr(idx_curr_snaptime)


	; set time window
	dt= 0.01
	idx_window= where((sfrtime ge (time[i]-dt)) and (sfrtime le (time[i]+dt)))
	sfr_avg[i]= -1
	if idx_window(0) ne -1 then sfr_avg[i]= mean(sfrsfr(idx_window))

	print, "T= ", time[i], sfr_inst[i], sfr_avg[i]

endfor



;--------------------------------------
;--------------------------------------


; SFR(snaptimes) FILE
; --------------------
openw, 1, frun+'/sfr_snaptimes.txt', ERROR=err

printf, 1, "#   sfr_snaptimes.txt"
printf, 1, "# "
printf, 1, "#           inst.       avg     from snap   "
printf, 1, "# time       sfr        sfr        sfr      "
printf, 1, "# (Gyr)    (Mo /Yr)   (Mo /Yr)   (Mo /Yr)   "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"    ",3(F8.3,"   "))', $
                time[i], sfr_inst[i], sfr_avg[i], sfr_snap[i]
endfor
close, 1




; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            "
print, "---------------------------------"
print, "---------------------------------"





end








;====================================================================



; Fit Gaussian to Starburst
;-----------------------------
pro fitandoverplot_gaussian, time, sfr, $
			SFR_max= SFR_max, width=width, SBtime=SBtime, $
			quiet=quiet

        x_tofit= time
        y_tofit= sfr
        weight_tofit= 1.0 + 0.0*sfr
        ;weight_tofit= 1.0 + 0.1*sfr


        ;  parameters
        ;---------------- 
	; p[0] = Normalization, i.e., the SFR_max
        ; p[1] = Width of the dist.; sigma
        ; p[2] = mean of the dist., i.e., time of SFR_max

        ; initial guess
        ;guess= [100.0,0.2,1.5]
	normz= max(sfr)
        idd= where(sfr eq normz)
        meanz= time(idd[0])
        wid= max(time)/3.0
        guess= [normz, wid, meanz]
print, guess

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

	idx= where(y gt 0.3*min(sfr))
	y= y(idx)
	x= x(idx)

	if not keyword_set(quiet) then begin
        	oplot, x, y, psym=-3, linestyle= 2, color= 50, thick= 6.0
	endif


        if not keyword_set(quiet) then begin
                ;xyouts, 0.22, 0.75, 'Gaussian ('+chi2+')', size=1.2, color=0, /normal
                ;xyouts, 0.23, 0.71, 'SFRmax= '+sfrmaxlbl, size=1.0, color=0, /normal
                ;xyouts, 0.23, 0.67, 'sigma=  '+wlbl, size=1.0, color=0, /normal
                ;xyouts, 0.23, 0.63, 'sbtime= '+sbtm, size=1.0, color=0, /normal
                ;xyouts, 0.65, 0.84, 'Gaussian', size=1.6, color=0, /normal
                ;xyouts, 0.65, 0.80, '!7r!6= '+wlbl+' Gyr', size=1.5, color=0, /normal
        endif

        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end




;====================================================================


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
pro fitandoverplot_expdecay_0, time, sfr, $
			SFR_SB= SFR_SB, tau=tau, SBtime=SBtime

        x_tofit= time
        y_tofit= sfr
        weight_tofit= 1.0 + 0.0*sfr


        ;  parameters
        ;---------------- 
	; p[0] = Normalization, i.e., the SFR_max
        ; p[1] = Width of the dist.; sigma
        ; p[2] = mean of the dist., i.e., time of SFR_max

        ; initial guess
        guess= [100.0,1.0]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},2)
        ; SFR_SB is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ; tau is greater than 0
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.001


        ; markwardt mpfit procedure
        tau_result = MPFITFUN('func_exp', x_tofit, y_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        ;chi2= strcompress(string(redchi2),/remove_all)
        ;chi2= strmid(chi2,0,4)

        SFR_SB= tau_result[0]
	print, "SFR_SB= ", SFR_SB
        ;rho0lbl= strcompress(string(rho0lbl),/remove_all)
        ;rho0lbl= strmid(rho0lbl,0,5)

        tau= tau_result[1]
	print, "SB tau= ", tau," Gyr"
        ;rslbl= strcompress(string(rslbl),/remove_all)
        ;rslbl= strmid(rslbl,0,5)


        ; overplot fit
        ; ----------
        x=time
        y= func_exp(x,tau_result)
        x=x+SBtime
        oplot, x, y, psym=-3, linestyle= 0, color= 100

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
                        SFR_SB= SFR_SB, tau=tau, SBtime=SBtime, decay_const=decay_const, $
                        quiet=quiet

        x_tofit= time
        y_tofit= sfr
        ;weight_tofit= 1.0 + 0.0*sfr
        weight_tofit= 0.1*sfr + 1.0


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
                oplot, x, y, psym=-3, linestyle= 0, color= 200, thick= 6.0
        endif

        ;xyouts, 0.65, 0.75, 'NFW', /normal, size= 1.0, color=150
        ;xyouts, 0.65, 0.71, 'r!Ds!N='+rslbl, /normal, size= 1.0, color=150
        ;xyouts, 0.77, 0.71, ', c='+clbl, /normal, size= 1.0, color=150


        if not keyword_set(quiet) then begin
                ;xyouts, 0.65, 0.71, 'Tau Model', color=0, size=1.5, /normal
                ;xyouts, 0.65, 0.66, '!7s!6!DSF!N= '+taulbl+' Gyr', color=0, size=1.5, /normal
        endif

        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end





;====================================================================




