;===============================================================================
;===============================================================================
;===============================================================================
;===============================================================================
;===============================================================================

pro process_and_plot_sfr_history, frun, lcolor=lcolor, lstyle=lstyle, lthick=lthick, $
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
		read_file_sfr_raw, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac
	endif else begin
		read_file_sfr_processed, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac
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



;---------------------------------------------

moveon:


if not keyword_set(y0) then y0= 0.92
if not keyword_set(x0) then x0= 0.75

if keyword_set(msg) then begin
    xyouts, x0, y0, msg, /normal, charthick=3.0, size=1.33, color= lcolor
endif else begin
    xyouts, x0, y0, frun, /normal, charthick=3.0, size=1.33, color= lcolor
endelse




;---------------------------------------------



if keyword_set(show_specific_times) eq 1 then begin

        for ii= 0, n_elements(show_specific_times)-1 do begin
                ;specific_time= 1.205
                specific_time= show_specific_times[ii]
                idx=where(sfrtime ge specific_time)
                oplot, [sfrtime[idx[0]]], [sfrsfr[idx[0]]], psym= 4, color= 150, thick=7.0, symsize=2.0
        endfor
endif



;---------------------------------------------

;
;  all of the below is fitting the history,
;  plotting it, and potentially writing it to a file
;


t_merg= fload_bh_mergertime(frun)
sfr_max_raw= max(sfrsfr)
t_sfr_max= sfrtime(where(sfrsfr eq sfr_max_raw))



do_meansf= 0
;do_meansf= 1
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



do_sb_fit= 0
;do_sb_fit= 1        ; first passage starburst
if do_sb_fit eq 1 then begin
        idx= where((sfrtime gt 0.1) and (sfrtime lt t_merg-0.2))
        time= sfrtime(idx)
        sfr= sfrsfr(idx)

        sfr_max_first_peak= max(sfr)
        t_first_peak= time(where(sfr eq sfr_max_first_peak))

        fitandoverplot_gaussian, time, sfr, $
                        SFR_max= SFR_max_fp, width=width_fp, SBtime=SBtime_fp, quiet=quiet
endif


do_sb_fit= 0
;do_sb_fit= 1        ; final merger starburst
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

        ;fitandoverplot_gaussian_pluszeropt, time, sfr, $
        fitandoverplot_gaussian, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime, quiet=quiet
endif



do_expdecay_fit= 0
;do_expdecay_fit= 1
if do_expdecay_fit eq 1 then begin
	idx= where(sfrtime gt t_merg)
	;idx= where(sfrtime gt SBtime)
	;idx= where((sfrtime gt SBtime) and (sfrtime lt SBtime+1.3))
	time= sfrtime(idx)
	time= time- time(0)
	sfr= sfrsfr(idx)
	fitandoverplot_expdecay, time, sfr, SBtime=SBtime, tau=tau, decay_const=decay_const
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


have_and_print_fits=0
;have_and_print_fits=1
if have_and_print_fits eq 1 then begin

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

end

;--------------------------------------
;--------------------------------------



;-------------------------------------------
;   Infer SFR from final snapshot 
;-------------------------------------------


;spawn, "/bin/ls "+frun+"/snapshot_* | wc ",result
;snapnum=long(result[0])-1
;
;ok=fload_snapshot_bh(frun,snapnum)
;id=fload_newstars_id(1)
;time=fload_time(1)
;formationtime=fload_newstars_age(1)
;m=fload_newstars_mass(1)
;print, min(m), max(m), mean(m)

;timebinsize= 0.01     ; in Gadget units = 10 Myr
;oplotit= 3
;histlevels= (xmax-xmin)/timebinsize

; normalization ??
;scaling= mean(m) * 1.0d10 / (timebinsize * 1.0d9)     ; conver particle mass to m_solar, and then we're uring 10 Myr bins
;normalization= 1.0/scaling

;idx= where(formationtime le xmax)
;formationtime= formationtime(idx)

;stop
; now the histogram
;temp= process_histogram(formationtime, xmax=xmax, xmin=xmin, levels=histlevels, $
;                                 oplotit=oplotit, mannorm=normalization, bins=bins)
;print, min(temp), max(temp)

;xyouts, 0.25, 0.75, "new star formation times", size=1.5, color=0, /normal, charthick=3.0


;=========================================
;
; now, display differences
;

;        diffsfr=temp*0.0
;        for i=0,long(n_elements(diffsfr)-2) do begin
;           idx=where(sfrtime gt (bins[i]-1.0e-6))
;           if idx(0) ne -1 then diffsfr[i]= temp[i]-sfrsfr(idx(0)) else print, "bad i= ", i
;        endfor

;        diffrange=1.0*long(max(abs(diffsfr))*1.2)
;        momentdiff= moment(diffsfr)
;        p= [1.0, sqrt(momentdiff(1)), momentdiff(0)]


;        wheredot= strlen(filename)-4
;        filename= strmid(filename,0,wheredot)+'_diff.eps'
;        initialize_plotinfo, 1
;        setup_plot_stuff, 'ps', filename=filename, colortable= 4

;        plot, [1.0],[1.0], psym=-3, xrange=[-diffrange,diffrange], yrange=[0.0,1.0], color= 0, $
;        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
;        charthick=3.0, xtitle='!7D!3 (M!D!9n!6!N yr!E-1!N)', ytitle=' ', /nodata

;        temp=process_histogram(diffsfr, xmax=diffrange, xmin=-diffrange, levels=20, oplotit=5)

;        x= diffrange * (2.0 * findgen(101)/100.0 - 1.0)
;        oplot, x, func_gaussian(x,p), color= 150, thick=2.0


;=========================================



;--------------------------------------
;--------------------------------------



end







;======================================================================



