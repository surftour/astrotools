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
   print, "  WARNING: cumulative and gasmass do not current work!  "
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfr.eps'

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
xaxistitle = "!6Time (Gyr)"  & h= 0.7

;xmax = 3.0
;xmax = 2.4
;xmax = 1.5
xmax = 1.0
xmin = 0

;ymax = 1000
ymax = 300
;ymax = 80
;ymax = 30
ymin= 1.0
;ymin= 0.7
;ymin= 0.01

ynotlog= 1

if keyword_set(ynotlog) then ymin = 0 

; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	h = fload_cosmology('h')
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


;process_one_sfr, "vc3vc3e_2", lcolor= 50, ctab= 4, lthick=2.0, msg='tilted', x0= 0.20, y0= 0.36, h=0.7
;process_one_sfr, "vc3vc3h_2", lcolor= 150, ctab= 4, lthick=2.0, msg='pro-pro', x0= 0.20, y0= 0.32, h=0.7


; ---------------------------------------------------

; isolated
;process_one_sfr, "isolated/test0", lcolor=   0, ctab= 4, lthick=2.0, msg='std', x0= 0.20, y0= 0.44, h=0.7
;process_one_sfr, "isolated/test1", lcolor= 150, ctab= 4, lthick=2.0, msg='m=0.5, n=1.0', x0= 0.20, y0= 0.40, h=0.7
;process_one_sfr, "isolated/test2", lcolor=  50, ctab= 4, lthick=2.0, msg='m=0.0, n=1.5', x0= 0.20, y0= 0.36, h=0.7
;process_one_sfr, "isolated/test3", lcolor= 100, ctab= 4, lthick=2.0, msg='m=1.0, n=1.0', x0= 0.20, y0= 0.32, h=0.7
;process_one_sfr, "isolated/test4", lcolor= 200, ctab= 4, lthick=2.0, msg='m=0.5, n=1.5', x0= 0.20, y0= 0.28, h=0.7
;process_one_sfr, "isolated/test5", lcolor= 135, ctab= 4, lthick=2.0, msg='m=1.0, n=1.0 (higher t!dsf!n)', x0= 0.20, y0= 0.24, h=0.7
;process_one_sfr, "isolated/test6", lcolor=  70, ctab= 4, lthick=2.0, msg='m=0.5, n=1.5 (higher t!dsf!n)', x0= 0.20, y0= 0.20, h=0.7


; ---------------------------------------------------

; Josh's z=5 mergers
;process_one_sfr, "z5/jyA", lcolor=  50, ctab= 4, lthick=2.0, msg='direct', x0= 0.20, y0= 0.34, h=0.7
;process_one_sfr, "z5/jyB", lcolor= 150, ctab= 4, lthick=2.0, msg='higher J', x0= 0.20, y0= 0.30, h=0.7


;xx0= 0.2 & yy0= 0.38   ; good for log
xx0= 0.65 & yy0=0.90   ; good for linear
process_one_sfr, "z5/jyDfg0.4_1", lcolor= 150, ctab= 4, lthick=2.0, msg='shock, m=0.5', x0= xx0, y0= yy0-0.04, h=0.7
process_one_sfr, "z5/jyDfg0.4_2", lcolor= 100, ctab= 4, lthick=2.0, msg='shock, m=1.0', x0= xx0, y0= yy0-0.08, h=0.7
process_one_sfr, "z5/jyDfg0.4_3", lcolor= 200, ctab= 4, lthick=2.0, msg='shock + rho', x0= xx0, y0= yy0-0.12, h=0.7
process_one_sfr, "z5/jyDfg0.4",   lcolor=  50, ctab= 4, lthick=2.0, msg='std', x0= xx0, y0= yy0, h=0.7

;process_one_sfr, "z5/jyCfg0.4_1", lcolor= 150, ctab= 4, lthick=2.0, msg='shock, m=0.5', x0= xx0, y0= yy0-0.04, h=0.7
;process_one_sfr, "z5/jyCfg0.4_2", lcolor= 100, ctab= 4, lthick=2.0, msg='shock, m=1.0', x0= xx0, y0= yy0-0.08, h=0.7
;process_one_sfr, "z5/jyCfg0.4_3", lcolor= 200, ctab= 4, lthick=2.0, msg='shock + rho', x0= xx0, y0= yy0-0.12, h=0.7
;process_one_sfr, "z5/jyCfg0.4",   lcolor=  50, ctab= 4, lthick=2.0, msg='std', x0= xx0, y0= yy0, h=0.7


; ---------------------------------------------------

; direct d's
;process_one_sfr, "ds/d3e1", lcolor=   0, ctab= 4, lthick=2.0, msg='d3e1', x0= 0.20, y0= 0.44, h=0.7
;process_one_sfr, "ds/d3f1", lcolor= 150, ctab= 4, lthick=2.0, msg='d3f1', x0= 0.20, y0= 0.40, h=0.7
;process_one_sfr, "ds/d3h1", lcolor=  50, ctab= 4, lthick=2.0, msg='d3h1', x0= 0.20, y0= 0.36, h=0.7
;process_one_sfr, "ds/d3k1", lcolor= 100, ctab= 4, lthick=2.0, msg='d3k1', x0= 0.20, y0= 0.32, h=0.7


; ---------------------------------------------------




;--------------------------------------
; done
;--------------------------------------
device, /close


end

    

















pro process_one_sfr, frun, lcolor=lcolor, lstyle=lstyle, lthick=lthick, $
				ctab=ctab, msg=msg, y0=y0, x0=x0, lpsym=lpsym, $
				cumulative=cumulative, $
				gasmass=gasmass, $
				normalize=normalize, h=h


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
	open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass


	; physical units
	if keyword_set(h) then sfrtime = sfrtime / h
	;if i eq 1 then sfrtime = sfrtime/(0.7)    ; did this for comparing h and noh

	if keyword_set(normalize) then sfrsfr=sfrsfr/sfrsfr(0)


	    lthick= 4.0 * lthick

	    if keyword_set(cumulative) then begin
		oplot, sfrtime, sfrmfs, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick
	    endif else begin
		if keyword_set(gasmass) then begin
			oplot, sfrtime, sfrgasmass, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick
		endif else begin
			oplot, sfrtime, sfrsfr, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick
		endelse
	    endelse



if not keyword_set(y0) then y0= 0.92
if not keyword_set(x0) then x0= 0.75

if keyword_set(msg) then begin
    xyouts, x0, y0, msg, /normal, charthick=3.0, size=1.33, color= lcolor
endif else begin
    xyouts, x0, y0, frun, /normal, charthick=3.0, size=1.33, color= lcolor
endelse


tpsym= 5
if lcolor eq 50 then tpsym= 4
if lcolor eq 150 then tpsym= 6


if ((frun eq 'vc3vc3e_2') or (frun eq 'vc3vc3h_2')) then begin
;if ((frun eq 'z5/jyA') or (frun eq 'z5/jyB')) then begin
;if ((frun eq 'ds/d3e1') or (frun eq 'ds/d3f1') or (frun eq 'ds/d3h1') or (frun eq 'ds/d3k1')) then begin
;if ((strmid(frun,0,6) eq 'z5/jyC') or (strmid(frun,0,6) eq 'z5/jyD')) then begin
	read_centerpositions, cp_time, cp_cen1, cp_cen2, filename="/raid4/tcox/"+frun+"/centers.txt"

	if keyword_set(h) then cp_time = cp_time / h

	x1= transpose(cp_cen1(0,*))
	y1= transpose(cp_cen1(1,*))
	z1= transpose(cp_cen1(2,*))

	x2= transpose(cp_cen2(0,*))
	y2= transpose(cp_cen2(1,*))
	z2= transpose(cp_cen2(2,*))

	rdiff= sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))

	rd= 1.0

	; loop until we find close ones
	repeat begin
	    rd= rd+1.0
	    idx= where((rdiff gt 7.0 - rd) and (rdiff lt 7.0 + rd))
	endrep until idx(0) ne -1

	print, idx
	print, rdiff(idx)

	for i= 0, n_elements(idx)-1 do begin
	    tdx= where(sfrtime ge cp_time(idx(i)))
	    oplot, [sfrtime(tdx[0])], [sfrsfr(tdx[0])], psym=tpsym, color= lcolor, symsize=2.0, thick= 4.0
	endfor

endif


end






;======================================================================






pro open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass

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
                        spawn, '/bin/ls '+frun+'/*.sfr ',result
                endif else begin
			case loopnum of
			  0: datadir='/raid4'
			  1: datadir='/home'
			  2: datadir='/raid2'
			  3: datadir='/data'
			  4: datadir='/data6'
			  5: datadir='/data7'
			  else: break
			endcase

			;sfrfile= datadir+'/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'

			nfrun= fload_getid(frun)
			spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/*.sfr ',result
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






;==============================================================================






;-----------------------------
pro fitandoverplot_gaussian, time, sfr, $
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
        guess= [100.0,0.2,1.5]

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


        ; overplot fit
        ; ----------
        x=min(time) - 0.2 + (0.2+max(time)-min(time))*findgen(100)/100.0
        y= func_gaussian(x,sb_result)
        oplot, x, y, psym=-3, linestyle= 0, color= 50

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




