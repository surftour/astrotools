pro sfr_multi, sendto, filename=filename, $
		cumulative=cumulative, $
		h=h


if not keyword_set(sendto) then begin
   print, "  "
   print, "sfr_multi, sendto, filename=filename, /h"
   print, "  "
   print, "  "
   return
endif


initialize_plotinfo, 1


setup_plot_stuff, sendto, filename=filename, colortable= 4
;setup_plot_stuff, sendto, filename=filename, colortable= 0



;-------------------------------------------
;   Load the runs to display
;-------------------------------------------


; ------------------
; Sbc major mergers
; ------------------


	tmerg= [1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8]
	msg= ' '
	lbltitleon= 0


; fig 6
; -------------
;fruns= ["Sbc201a-u4","Sbc201a-u43"] & lbls=["n2low","n2high"]   ; paper version
;fruns= ["Sbc201a-u53","Sbc201a-u52"] & lbls=["n1low","n1high"]   ; paper version
;fruns= ["Sbc201a-u50","Sbc201a-u51"] & lbls=["n0low","n0high"]   ; paper version
;  fruns= ["Sbc201a-u4","Sbc201a-u43","Sbc201a-u53","Sbc201a-u52","Sbc201a-u50","Sbc201a-u51"]  - all of them

;fruns= ["Sbc201a-u57","Sbc201a-u4","Sbc201a-u43"] & lbls=["!8n2low!6","!8n2med!6","!8n2high!6"]   ; paper version
;fruns= ["Sbc201a-u56","Sbc201a-u53","Sbc201a-u52"] & lbls=["!8n1low!6","!8n1med!6","!8n1high!6"]   ; paper version
;fruns= ["Sbc201a-u55","Sbc201a-u50","Sbc201a-u51"] & lbls=["!8n0low!6","!8n0med!6","!8n0high!6"]   ; paper version



; fig 12
; ---------
; resolution
;fruns= ["Sbc201a-u4", "Sbc201a-n4", "Sbc201a2x-n4", "Sbc201a4x-n4", "Sbc201a10x-n4"]
;fruns= ["Sbc201a-u4", "Sbc201a2x-n4", "Sbc201a4x-n4", "Sbc201a10x-n4", "Sbc201a10x-u4"]   ; paper version
;lbls= ['1x (n2low)', '2x', '4x', '10x', '10x, h=0.047']    ; paper version
;fruns= ["Sbc201a-u4", "Sbc201a10x-n4", "Sbc201a10x-u4"]
;lbls= ['n2low', '10x','10x-h']
;fruns= ["Sbc201a10x-n4", "Sbc201a10x-u4"]
;lbls= ['h=0.100','h=0.047']




; fig 13
; --------
; sph version
;fruns= ["Sbc201a-u4", "Sbc201a-u4d", "Sbc201a-u4c"]
;lbls= ["entropy (n2low)","energy, standard","energy, geometric"]
;fruns= ["Sbc201a-u4", "Sbc201a-u4d", "Sbc201a-u4c","Sbc201a-u4r"]
;lbls= ["entropy (n2low)","energy, standard","energy, geometric","energy, asymmetric"]
;fruns= ["Sbc201a-u4", "Sbc201a-u4d", "Sbc201a-u4c","Sbc201a-u4u","Sbc201a-u4v"]
;lbls= ["entropy (n2low)","energy, standard","energy, geometric","energy, asymmetric","entropy, standard"]

; test processor number
;fruns= ["Sbc201a-u4", "Sbc201a-u4s", "Sbc201a-u4r"]
;lbls= ["std","std - 4 nodes", "e,a - 4 nodes"]



; fig 14
; ---------
; cstar variations
;fruns= ["Sbc201a-u74","Sbc201a-u64","Sbc201a-u4","Sbc201a-u44","Sbc201a-u14","Sbc201a-u54"]
;       lbls= ["c!I*!N=0.3","c!I*!N=0.06","n2low","c!I*!N=0.015","c!I*!N=0.01","c!I*!N=0.004"]
;fruns= ["Sbc201a-u74","Sbc201a-u64","Sbc201a-u4","Sbc201a-u44","Sbc201a-u54"]
;       lbls= ["c!I*!N=0.3","c!I*!N=0.06","c!I*!N=0.03 (n2low)","c!I*!N=0.015","c!I*!N=0.004"]
;        lbls= ["c!I* !N=0.3","  !I !N 0.06","  !I !N 0.03 (n2low)","  !I !N 0.015","  !I !N 0.004"]




; fig 16
; -----------
;fruns= ["Sbc201a-u4","Sbc201a-u4e"]
;fruns= ["Sbc201a-u4","Sbc201a-u4q"]
;lbls= ["n2low","isothermal"]




; fig 17
; --------
; different cooling
;fruns= ["Sbc201a-u4", "Sbc201a-u4f", "Sbc201a-u4g", "Sbc201a-u4h", "Sbc201a-u4i"]
;	lbls= ['KWH','mzero.cie','m-30.cie','m-05.cir','m-00.cie']
;fruns= ["Sbc201a-u4", "Sbc201a-u4f", "Sbc201a-u4g", "Sbc201a-u4k", "Sbc201a-u4j", "Sbc201a-u4i"]
;       lbls= ['KWH','mzero.cie','m-30.cie','m-20.cie','m-10.cie','m-00.cie']
;fruns= ["Sbc201a-u4", "Sbc201a-u4f", "Sbc201a-u4g", "Sbc201a-u4k", "Sbc201a-u4l", "Sbc201a-u4j", "Sbc201a-u4i"]
       ;lbls= ['KWH','zero','10!E-3!N','10!E-2!N','10!E-1.5!N','10!E-1!N','solar']
;fruns= ["Sbc201a-u4", "Sbc201a-u4g", "Sbc201a-u4k", "Sbc201a-u4l", "Sbc201a-u4j", "Sbc201a-u4i"]  ; - paper version
;lbls= ['0 (n2low)','10!E-3!N','10!E-2!N','10!E-1.5!N','10!E-1!N','1']  ; - paper version

;lbltitleon= 1
;lbltitle='Metallicity (Z!D!9n!3!N)'
       ;lbls= ['zero (n2low)','-3','-2','-1.5','-1','solar']
;fruns= ["Sbc201a-u4f", "Sbc201a-u4i"]
;       lbls= ['mzero.cie','m-00.cie']






; fig A1
; --------
;fruns= ["mhmaj-u6","mhmaj-u7","mhmaj-u5"]
;lbls= ["isothermal","energy: geometric","entropy"]



; fig A2
; --------
; springel comparison
;fruns= ["A1m-u1","A1m-u2","A1m-u3"]
;lbls= ["entropy", "energy: standard (S00)", "energy: geometric"]
fruns= ["A1m-u2", "A1m-u1","A1m-u7","A1m-u8"]
lbls= ["S00", "S00, new SPH", "!8n2low!6, new SPH", "!8n2med!6, new SPH"]



; fig A3
; --------
;fruns= ["A1m-u1","A1m-u4","A1m-u7"]
;lbls= ["c!I*!N=0.004 (S00)","c!I*!N=0.012", "c!I*!N=0.03"]
;fruns= ["A1m-u1","A1m-u7"]
;lbls= ["c!I*!N=0.004 (S00)","  !I !N 0.03"]






; -------------------------------------
; nothing
; -------
; different bulge
;fruns= ["Sbc201a-u4", "Sbc201a-u4hbl"]
;       lbls= ['Exponential','Hernquist']
;
; -------------------------------------



zcolor= 50
deltacolor= (250-zcolor)/n_elements(fruns)
;deltacolor= 50
lstyle= 0
zcolor_orig= zcolor

;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6 SFR (M!D!9n!6!N yr!E-1!N)"
if keyword_set(cumulative) then yaxistitle="Mass (10!e10!n h!E-1!NM!D!9n!6!N)"

xaxistitle = "!6Time (!8h!6!e-1!N Gyr)"
;xaxistitle = "Time (Gyr)"

;xmax = 14.0
;xmax = 10.0
;xmax = 7.5
;xmax = 6.0
;xmax = 5.0
;xmax = 4.25
;xmax = 4.0
;xmax = 3.0
xmax = 2.8
;xmax = 1.5
;xmax = 1.3
;xmax = 1.0
xmin = 0

;ymax = 300
;ymax = 180
;ymax = 120
;ymax = 90
;ymax = 80
;ymax = 60
ymax = 15
;ymax = 13
;ymax = 12
;ymax = 10.0
;ymax = 8.0
;ymax = 6.0
;ymax = 5.0
;ymax = 2.0
;ymax = 1.5
;ymax = 1.3
;ymax = 1.2
;ymax = 0.75
;ymax = 0.6
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

;stop


if sendto EQ 'x' then begin
;	oplot, sfrtime, sfrsfr, psym=-3, color= getcolor('blue')
endif else begin
	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	ifinal= n_elements(fruns)-1
	for i=0,ifinal do begin
	    ; get sfr data
	    sfrfile= '/data/tcox/sfr/'+fload_getid(fruns[i])+'.sfr'

	    get_lun, unit
	    openr, unit, sfrfile, ERROR=err

	    if (err NE 0) then begin
	        print, "  "
	        print, "Problem: ",!ERR_STRING
	        print, "  "
	        close, unit
	        ERR= 0
	        sfrtime= [0]
	        sfrsfr= [0]
	        sfrmsg= "No Star Formation"
	    endif else begin
	        close, unit
	        print, "opening: ",sfrfile
	        sfrdata= read_ascii(sfrfile)
	        sfrtime= sfrdata.field1[0,*]
	        sfrsfr= sfrdata.field1[1,*]
		sfrmfs= sfrdata.field1[3,*]
		n_cols= n_elements(sfrdata.field1[*,0])
		n_rows= n_elements(sfrdata.field1[0,*])
		finalnewstarmass= sfrmfs[n_rows-1]
		if n_cols ge 6 then sfrgasmass= sfrdata.field1[6,*] else sfrgasmass=[20.0,20.0-finalnewstarmass]
	        sfrmsg= ''
	    endelse

	    ; physical units
	    if keyword_set(h) then sfrtime = sfrtime / h
	    ;if i eq 1 then sfrtime = sfrtime/(0.7)    ; did this for comparing h and noh

; this was for mhmaj-u* information, don't really need 
;  on a regular basis.
;idx= where(sfrtime le 1.3)
;if idx(0) ne -1 then begin
;	sfrtime= sfrtime(idx)
;	sfrsfr= sfrsfr(idx)
;	sfrmfs= sfrmfs(idx)
;	if n_elements(sfrgasmass) gt 2 then sfrgasmass= sfrgasmass(idx)
;endif

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

            if keyword_set(cumulative) then begin
                oplot, sfrtime, sfrmfs, psym=-3, color= zcolor, linestyle= lstyle, thick= 4.0
            endif else begin
                oplot, sfrtime, sfrsfr, psym=-3, color= zcolor, linestyle= lstyle, thick= 4.0
            endelse

	    lstyle= lstyle+1
	    zcolor= zcolor+deltacolor
	endfor
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




xyouts, 0.25, 0.85, msg, /normal, charthick=3.0, size=1.5, color= 0



; std
;std= 0
std= 1
if std eq 1 then begin
    ;x0= 0.75
    ;x0= 0.23
    x0= 0.21
    ;y0= 0.92
    ;y0= 0.90
    y0= 0.93
    zcolor= zcolor_orig
    if lbltitleon eq  1 then begin
	xyouts, x0, y0, lbltitle, /normal, charthick=3, size=1.5, color= 0
    endif
    for i=0, n_elements(fruns)-1 do begin
	y0= y0-0.05
	if n_elements(lbls) gt 0 then begin
                xyouts, x0, y0, lbls[i], /normal, charthick=3, size=1.50, color=zcolor
                    zcolor= zcolor+deltacolor
	endif else begin
		xyouts, x0, y0, fruns[i], /normal, charthick=3, size=1.33, color=zcolor
	            zcolor= zcolor+deltacolor
	endelse
    endfor
endif










;--------------------------------------
;--------------------------------------

if (sendto EQ 'ps') then device, /close


end


