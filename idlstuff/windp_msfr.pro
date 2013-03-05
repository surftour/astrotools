pro sfr1, junk


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

filename='sfr1.eps'
filename='sfr2.eps'
filename='sfr3.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

xmax = 4.25
xmin = 0

ymax = 1000.0
ymin = 0.01



;---------------------------

!p.position= [0.16, 0.14, 0.98, 0.98]

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




;-------------------------------------------
;   Load the runs to display
;-------------------------------------------

;process_one_sfr, "ds/d1e2_q", lcolor=0, lthick= 2.0, msg='nosb', x0= 0.75, y0= 0.88, h=h
;;process_one_sfr, "sb8_mass/d1e", lcolor=200, lthick= 2.0, msg='sb8', x0= 0.75, y0= 0.84, h=h
;process_one_sfr, "sb8_mass/d1e_no", lcolor=190, lthick= 2.0, msg='sb8no', x0= 0.75, y0= 0.80, h=h
;;process_one_sfr, "sb10_mass/d1e", lcolor=50, lthick= 2.0, msg='sb10', x0= 0.75, y0= 0.76, h=h
;process_one_sfr, "sb10_mass/d1e_no", lcolor=40, lthick= 2.0, msg='sb10no', x0= 0.75, y0= 0.72, h=h
;;process_one_sfr, "sb13_mass/d1e", lcolor=100, lthick= 2.0, msg='sb13', x0= 0.75, y0= 0.68, h=h
;process_one_sfr, "sb13_mass/d1e_no", lcolor=90, lthick= 2.0, msg='sb13no', x0= 0.75, y0= 0.66, h=h

;process_one_sfr, "ds/d2e2_q", lcolor=0, lthick= 2.0, msg='nosb', x0= 0.75, y0= 0.88, h=h
;process_one_sfr, "sb8_mass/d2e", lcolor=200, lthick= 2.0, msg='sb8', x0= 0.75, y0= 0.84, h=h
;process_one_sfr, "sb10_mass/d2e", lcolor=50, lthick= 2.0, msg='sb10', x0= 0.75, y0= 0.80, h=h
;process_one_sfr, "sb13_mass/d2e", lcolor=100, lthick= 2.0, msg='sb13', x0= 0.75, y0= 0.76, h=h

;process_one_sfr, "ds/d4e2_q", lcolor=0, lthick= 2.0, msg='nosb', x0= 0.75, y0= 0.88, h=h
;process_one_sfr, "sb8_mass/d4e_no", lcolor=200, lthick= 2.0, msg='sb8', x0= 0.75, y0= 0.84, h=h
;process_one_sfr, "sb10_mass/d4e_no", lcolor=50, lthick= 2.0, msg='sb10', x0= 0.75, y0= 0.80, h=h
;process_one_sfr, "sb13_mass/d4e_no", lcolor=100, lthick= 2.0, msg='sb13', x0= 0.75, y0= 0.76, h=h

process_one_sfr, "ds/d5e2_q", lcolor=0, lthick= 2.0, msg='nosb', x0= 0.75, y0= 0.88, h=h
process_one_sfr, "sb8_mass/d5e_no", lcolor=200, lthick= 2.0, msg='sb8', x0= 0.75, y0= 0.84, h=h
process_one_sfr, "sb10_mass/d5e_no", lcolor=50, lthick= 2.0, msg='sb10', x0= 0.75, y0= 0.80, h=h
process_one_sfr, "sb13_mass/d5e_no", lcolor=100, lthick= 2.0, msg='sb13', x0= 0.75, y0= 0.76, h=h

;process_one_sfr, "ds/d6e2_q", lcolor=0, lthick= 2.0, msg='nosb', x0= 0.75, y0= 0.88
;process_one_sfr, "sb8_mass/d6e", lcolor=200, lthick= 2.0, msg='sb8', x0= 0.75, y0= 0.84
;process_one_sfr, "sb10_mass/d6e", lcolor=50, lthick= 2.0, msg='sb10', x0= 0.75, y0= 0.80
;process_one_sfr, "sb13_mass/d6e", lcolor=100, lthick= 2.0, msg='sb13', x0= 0.75, y0= 0.76


;process_one_sfr, "ds/d6e", lcolor=100, lthick= 2.0, msg='d6e', x0= 0.75, y0= 0.88
;process_one_sfr, "ds/d6e2", lcolor=50, lthick= 2.0, msg='d6e2', x0= 0.75, y0= 0.84
;process_one_sfr, "ds/d6e2_q", lcolor=150, lthick= 2.0, msg='d6e2_q', x0= 0.75, y0= 0.80

;process_one_sfr, "bs/b5e", lcolor=0, lthick= 2.0, msg='b5e', x0= 0.75, y0= 0.92
;process_one_sfr, "ds/d5e", lcolor=100, lthick= 2.0, msg='d5e', x0= 0.75, y0= 0.88
;process_one_sfr, "ds/d5e2", lcolor=50, lthick= 2.0, msg='d5e2', x0= 0.75, y0= 0.84
;process_one_sfr, "ds/d5e2_q", lcolor=150, lthick= 2.0, msg='d5e2_q', x0= 0.75, y0= 0.80

;process_one_sfr, "bs/b4e", lcolor=0, lthick= 2.0, msg='b4e', x0= 0.75, y0= 0.92
;process_one_sfr, "ds/d4e", lcolor=100, lthick= 2.0, msg='d4e', x0= 0.75, y0= 0.88
;process_one_sfr, "ds/d4e2", lcolor=50, lthick= 2.0, msg='d4e2', x0= 0.75, y0= 0.84
;process_one_sfr, "ds/d4e2_q", lcolor=150, lthick= 2.0, msg='d4e2_q', x0= 0.75, y0= 0.80


;
;process_one_sfr, "bs/b3e_no", lcolor=0, lthick= 2.0, msg='b3e_no', x0= 0.75, y0= 0.88
;process_one_sfr, "bs/b3e", lcolor=200, lthick= 2.0, msg='b3e', x0= 0.75, y0= 0.84




;--- xxxxx  55555  xxxxx ---
;process_one_sfr, "bs/b5e", lcolor=0, lthick= 2.0, msg='b5e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb10_mass/b5e", lcolor=100, lthick= 2.0, msg='sb10/b5e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb10_mass/b5e_no", lcolor=150, lthick= 2.0, msg='sb10/b5e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d5e2_q", lcolor=0, lthick= 2.0, msg='d5e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb10_mass/d5e", lcolor=100, lthick= 2.0, msg='sb10/d5e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb10_mass/d5e_no", lcolor=150, lthick= 2.0, msg='sb10/d5e_no', x0= 0.70, y0= 0.80
;---
;process_one_sfr, "bs/b5e", lcolor=0, lthick= 2.0, msg='b5e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb13_mass/b5e", lcolor=100, lthick= 2.0, msg='sb13/b5e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb13_mass/b5e_no", lcolor=150, lthick= 2.0, msg='sb13/b5e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d5e2_q", lcolor=0, lthick= 2.0, msg='d5e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb13_mass/d5e", lcolor=100, lthick= 2.0, msg='sb13/d5e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb13_mass/d5e_no", lcolor=150, lthick= 2.0, msg='sb13/d5e_no', x0= 0.70, y0= 0.80
;---
;process_one_sfr, "bs/b5e", lcolor=0, lthick= 2.0, msg='b5e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb8_mass/b5e", lcolor=100, lthick= 2.0, msg='sb8/b5e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb8_mass/b5e_no", lcolor=150, lthick= 2.0, msg='sb8/b5e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d5e2_q", lcolor=0, lthick= 2.0, msg='d5e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb8_mass/d5e", lcolor=100, lthick= 2.0, msg='sb8/d5e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb8_mass/d5e_no", lcolor=150, lthick= 2.0, msg='sb8/d5e_no', x0= 0.70, y0= 0.80


;process_one_sfr, "bs/b5e", lcolor=0, lthick= 2.0, msg='b5e', x0= 0.70, y0= 0.88
;process_one_sfr, "bs/b5e_no", lcolor=50, lthick= 2.0, msg='b5e_no', x0= 0.70, y0= 0.84
;process_one_sfr, "bs/b5e_no2", lcolor=100, lthick= 2.0, msg='b5e_no2', x0= 0.70, y0= 0.80
;process_one_sfr, "bs/b5e_no3", lcolor=150, lthick= 2.0, msg='b5e_no3', x0= 0.70, y0= 0.76




;--- xxxxx  33333  xxxxx ---
;process_one_sfr, "bs/b3e", lcolor=0, lthick= 2.0, msg='b3e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb10_mass/b3e", lcolor=100, lthick= 2.0, msg='sb10/b3e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb10_mass/b3e_no", lcolor=150, lthick= 2.0, msg='sb10/b3e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d3e7", lcolor=0, lthick= 2.0, msg='d3e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb10_mass/d3e", lcolor=100, lthick= 2.0, msg='sb10/d3e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb10_mass/d3e_no", lcolor=150, lthick= 2.0, msg='sb10/d3e_no', x0= 0.70, y0= 0.80
;---
;process_one_sfr, "bs/b3e", lcolor=0, lthick= 2.0, msg='b3e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb13_mass/b3e", lcolor=100, lthick= 2.0, msg='sb13/b3e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb13_mass/b3e_no", lcolor=150, lthick= 2.0, msg='sb13/b3e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d3e7", lcolor=0, lthick= 2.0, msg='d3e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb13_mass/d3e", lcolor=100, lthick= 2.0, msg='sb13/d3e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb13_mass/d3e_no", lcolor=150, lthick= 2.0, msg='sb13/d3e_no', x0= 0.70, y0= 0.80
;---
;process_one_sfr, "bs/b3e", lcolor=0, lthick= 2.0, msg='b3e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb8_mass/b3e", lcolor=100, lthick= 2.0, msg='sb8/b3e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb8_mass/b3e_no", lcolor=150, lthick= 2.0, msg='sb8/b3e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d3e7", lcolor=0, lthick= 2.0, msg='d3e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb8_mass/d3e", lcolor=100, lthick= 2.0, msg='sb8/d3e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb8_mass/d3e_no", lcolor=150, lthick= 2.0, msg='sb8/d3e_no', x0= 0.70, y0= 0.80
;--- xxxxx  11111  xxxxx ---
;process_one_sfr, "bs/b1e", lcolor=0, lthick= 2.0, msg='b1e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb10_mass/b1e", lcolor=100, lthick= 2.0, msg='sb10/b1e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb10_mass/b1e_no", lcolor=150, lthick= 2.0, msg='sb10/b1e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d1e2_q", lcolor=0, lthick= 2.0, msg='d1e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb10_mass/d1e", lcolor=100, lthick= 2.0, msg='sb10/d1e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb10_mass/d1e_no", lcolor=150, lthick= 2.0, msg='sb10/d1e_no', x0= 0.70, y0= 0.80
;---
;process_one_sfr, "bs/b1e", lcolor=0, lthick= 2.0, msg='b1e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb13_mass/b1e", lcolor=100, lthick= 2.0, msg='sb13/b1e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb13_mass/b1e_no", lcolor=150, lthick= 2.0, msg='sb13/b1e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d1e2_q", lcolor=0, lthick= 2.0, msg='d1e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb13_mass/d1e", lcolor=100, lthick= 2.0, msg='sb13/d1e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb13_mass/d1e_no", lcolor=150, lthick= 2.0, msg='sb13/d1e_no', x0= 0.70, y0= 0.80
;---
;process_one_sfr, "bs/b1e", lcolor=0, lthick= 2.0, msg='b1e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb8_mass/b1e", lcolor=100, lthick= 2.0, msg='sb8/b1e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb8_mass/b1e_no", lcolor=150, lthick= 2.0, msg='sb8/b1e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d1e2_q", lcolor=0, lthick= 2.0, msg='d1e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb8_mass/d1e", lcolor=100, lthick= 2.0, msg='sb8/d1e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb8_mass/d1e_no", lcolor=150, lthick= 2.0, msg='sb8/d1e_no', x0= 0.70, y0= 0.80




; -------------------------------------
;  Volker's Runs

;fruns= ["aA1","aA1w"] & lbls=['coll_1/vc1/nbh/','coll_1/vc1/bh']
;fruns= ["aA2","aA2w"] & lbls=['coll_1/vc2/nbh/','coll_1/vc2/bh']
;fruns= ["aA3","aA3w"] & lbls=['coll_1/vc3/nbh/','coll_1/vc3/bh']
;fruns= ["aA4","aA4w"] & lbls=['coll_1/vc4/nbh/','coll_1/vc4/bh']

;fruns= ["bA1","bA1w"] & lbls=['coll_2/vc1/nbh/','coll_2/vc1/bh']
;fruns= ["bA2","bA2w"] & lbls=['coll_2/vc2/nbh/','coll_2/vc2/bh']
;fruns= ["bA3","bA3w"] & lbls=['coll_2/vc3/nbh/','coll_2/vc3/bh']
;fruns= ["bA4","bA4w"] & lbls=['coll_2/vc4/nbh/','coll_2/vc4/bh']

;msg=''


;fruns= ["As/A3", "As/A3_10x"]
;fruns= ["As/A3_1x", "As/A3_10x", "As/A3_10xa"]
;fruns= ["As/A3", "As/A3_1x", "As/A3_10x"]
;fruns= ["As/A3", "As/A3nopot"]
;fruns= ["As/A3", "As/A3duprr1"]
msg= ' '




; -------------------------------------

; remnant mergers

; vc3 remnants
;fruns= ["vc3vc3_wBH", "vc3rem_vc3","vc3rem_vc3rem"]
;	tmreg=[1.1,1.1,1.1]
;	msg=''


;fruns= ["vc3vc3", "vc3rem_vc3","vc3rem_vc2"]
;fruns= ["vc3vc3", "vc3rem_vc3","vc3rem_vc2","vc3rem_vc1","vc3rem_vc1a","vc3rem_vc1b","vc3rem_vc1c"]
;	tmreg=[1.1,1.1,1.1]
;	msg=''



; -------------------------------------

; gaseous balls

; group A
;fruns= ["grpA_1","grpA_2","grpA_3","grpA_4","grpA_5", "grpA_no","grpA_no_1"]
;lbls= ["1e5","1e2","1e1","1e9","1e8", "grpA_no","grpA_no_1"]
;fruns= ["grpA_no_1","grpA_3","grpA_2","grpA_1","grpA_5","grpA_4"]
;lbls= ['none','1e1','1e2','1e5','1e8','1e9']
;msg= ' '




; -------------------------------------

; five to one mergers

;fruns= ["fivetoone/vc3V"]
;msg= ' '



; -------------------------------------

; local group

;fruns= ["localgroup/v2", "localgroup/v2_noigm"]
;fruns= ["localgroup/v3", "localgroup/v3_noigm"]
;lbls= ["!6The Local Group", "Isolated MW + M31"]
msg= ' '


; -------------------------------------

; Avishai's direct hit

;process_one_sfr, "crazy", lcolor=150, lthick= 2.0, msg='direct hit sim', x0= 0.20, y0= 0.88
msg= ' '



; ---------------------------------------------------




; std
std= 0
;std= 1
if std eq 1 then begin
    x0= 0.70
    ;x0= 0.65
    ;x0= 0.55
    ;y0= 0.40
    ;x0= 0.35
    y0= 0.92
    zcolor= zcolor_orig
    for i=0, n_elements(fruns)-1 do begin
        y0= y0-0.04
        if n_elements(lbls) gt 0 then begin
                xyouts, x0, y0, lbls[i], /normal, charthick=3, size=1.33, color=zcolor
                    zcolor= zcolor+deltacolor
        endif else begin
                xyouts, x0, y0, fruns[i], /normal, charthick=3, size=1.33, color=zcolor
                    zcolor= zcolor+deltacolor
        endelse
    endfor
endif





; ---------------------------------------------------








;--------------------------------------
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



; cooling
;xyouts, 0.7, 0.86, "Cooling", /normal, charthick=1, size=1.33, color=0
;oplot, [1.65,2.2],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.7, 0.80, "Explicit", /normal, charthick=1, size=1.33, color=50
;xyouts, 0.7, 0.75, "mzero.cie", /normal, charthick=1, size=1.33, color=100
;xyouts, 0.7, 0.70, "m-00.cie", /normal, charthick=1, size=1.33, color=150

; sfr law
;;xyouts, 0.7, 0.86, "SF Law", /normal, charthick=1, size=1.33, color=0
;xyouts, 0.74, 0.86, '!7s!3!N!D!20E!3!N', /normal, charthick=1, size=1.33, color=0
;oplot, [1.65,2.2],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.7, 0.80, "constant", /normal, charthick=1, size=1.33, color=50
;xyouts, 0.7, 0.75, "dynamical", /normal, charthick=1, size=1.33, color=100
;xyouts, 0.7, 0.70, '1/!7q!3!N', /normal, charthick=1, size=1.33, color=150

; resolution
;xyouts, 0.7, 0.86, "Resolution", /normal, charthick=1, size=1.33, color=0
;oplot, [1.65,2.2],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.7, 0.80, '4', /normal, charthick=1, size=1.33, color=50
;xyouts, 0.7, 0.76, '2', /normal, charthick=1, size=1.33, color=75
;xyouts, 0.7, 0.72, '1', /normal, charthick=1, size=1.33, color=100
;xyouts, 0.7, 0.68, '1/2', /normal, charthick=1, size=1.33, color=125
;xyouts, 0.7, 0.64, '1/4', /normal, charthick=1, size=1.33, color=150
;xyouts, 0.7, 0.60, '1/8', /normal, charthick=1, size=1.33, color=175

; feedback - within on n
;xyouts, 0.65, 0.86, 'n=2 Feedback', /normal, charthick=1, size=1.33, color=0
;;xyouts, 0.65, 0.86, '!7b!3 for n=2', /normal, charthick=1, size=1.33, color=0
;oplot, [1.50,2.35],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.65, 0.80, 'high', /normal, charthick=1, size=1.33, color=50
;xyouts, 0.65, 0.75, 'medium high', /normal, charthick=1, size=1.33, color=87
;xyouts, 0.65, 0.70, 'medium low', /normal, charthick=1, size=1.33, color=124
;xyouts, 0.65, 0.65, 'low', /normal, charthick=1, size=1.33, color=161

; alternate feedback n's
;xyouts, 0.55, 0.86, 'Alternate Feedback n', /normal, charthick=1, size=1.33, color=0
;oplot, [1.20,2.35],[13.0,13.0], psym=-3, color= 0, linestyle= 0
;xyouts, 0.6, 0.80, 'n=1 high fb', /normal, charthick=1, size=1.33, color=50
;xyouts, 0.6, 0.76, 'n=1 medium fb', /normal, charthick=1, size=1.33, color=75
;xyouts, 0.6, 0.72, 'n=1 low fb', /normal, charthick=1, size=1.33, color=100
;xyouts, 0.6, 0.68, 'n=0 high fb', /normal, charthick=1, size=1.33, color=125
;xyouts, 0.6, 0.64, 'n=0 medium fb', /normal, charthick=1, size=1.33, color=150
;xyouts, 0.6, 0.60, 'n=0 low fb', /normal, charthick=1, size=1.33, color=175

; mh comparison
;xyouts, 0.6, 0.85, "Isothermal Gas", /normal, charthick=1, size=1.33, color=50
;xyouts, 0.6, 0.80, "Full Gas Treatment", /normal, charthick=1, size=1.33, color=150


; springel comparison
;xyouts, 0.3, 0.80, "Conservative Entropy", /normal, charthick=1, size=1.33, color=50
;xyouts, 0.3, 0.75, "Arithmetic", /normal, charthick=1, size=1.33, color=100
;xyouts, 0.3, 0.70, "Geometric", /normal, charthick=1, size=1.33, color=150


; minor paper - iso figure
; ------------------------------
;xyouts, 0.80, 0.84, 'G3', /normal, charthick=3, size=1.33, color=zcolor_orig+deltacolor+deltacolor+deltacolor
;xyouts, 0.80, 0.38, 'G2', /normal, charthick=3, size=1.33, color=zcolor_orig+deltacolor+deltacolor
;xyouts, 0.80, 0.22, 'G1', /normal, charthick=3, size=1.33, color=zcolor_orig+deltacolor
;xyouts, 0.61, 0.36, 'G0', /normal, charthick=3, size=1.33, color=zcolor_orig
;toppos= 0.32
;botpos= 0.01
;xarrow= 0.60
;yarrow= 0.71
;arrow, xarrow, toppos, yarrow, botpos, /data, COLOR=zcolor_orig, THICK=2.0, hthick=3.0



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
	t_merg= fload_bh_mergertime('/raid4/tcox/'+frun)
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








pro sfr_at_snaptimes, frun, $
		filename=filename, $
		h=h


if not keyword_set(frun) then begin
   print, "  "
   print, "sfr_at_snaptimes, frun, filename=filename, /h"
   print, "  "
   print, "  "
   return
endif


if not keyword_set(filename) then filename="sfr_snaptimes.txt"

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;-------------------------------------------


; determine the number of
; snapshots in frun directory
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])




; ---------------------------------------------------


time= fltarr(nsnaps)

sfr_inst= fltarr(nsnaps)
sfr_avg= fltarr(nsnaps)


ymin = 0 


; physical units
if keyword_set(h) then begin
	h = fload_cosmology('h')
endif


; ---------------------------------------------------


	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass



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
        ok=fload_snapshot_bh(frun,i)


        ; what time is it?
        time[i]= fload_time(1)

	idx_gtr_snaptime= where(sfrtime ge time[i])
	idx_curr_snaptime= idx_gtr_snaptime[0]


	; instantaneous sfr
	sfr_inst[i]= sfrsfr(idx_curr_snaptime)


	; set time window
	dt= 0.01
	idx_window= where((sfrtime ge (time[i]-dt)) and (sfrtime le (time[i]+dt)))
	sfr_dt= sfrsfr(idx_window)
	sfr_avg[i]= mean(sfr_dt)

	print, "T= ", time[i], sfr_inst[i], sfr_avg[i]

endfor



;--------------------------------------
;--------------------------------------


; SFR(snaptimes) FILE
; --------------------
openw, 1, frun+'/sfr_snaptimes.txt', ERROR=err

printf, 1, "#   sfr_snaptimes.txt"
printf, 1, "# "
printf, 1, "#           inst.       avg      "
printf, 1, "# time       sfr        sfr      "
printf, 1, "# (Gyr)    (Mo /Yr)   (Mo /Yr)   "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"    ",2(F8.3,"   "))', $
                time[i], sfr_inst[i], sfr_avg[i]
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




