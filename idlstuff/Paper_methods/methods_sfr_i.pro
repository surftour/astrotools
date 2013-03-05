pro methods_sfr_i, sendto, filename=filename, $
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


; methods paper
; -------------
;fruns= ["Sbc11i4-u40","Sbc11i4-u46"] & lbls=["n1low","n1high"]
;  fruns= ["Sbc11i4-u53","Sbc11i4-u52"] & lbls=["n1low","n1high"]
;  fruns= ["Sbc11i4-u56","Sbc11i4-u53","Sbc11i4-u52"] & lbls=["n1low","n1med","n1high"]

;fruns= ["Sbc11i4-u8","Sbc11i4-u45"] & lbls=["n0low","n0high"]
;  fruns= ["Sbc11i4-u50","Sbc11i4-u51"] & lbls=["n0low","n0high"]
;fruns= ["Sbc11i4-u67","Sbc11i4-u50"] & lbls=["n0toolow","n0low"]
;  fruns= ["Sbc11i4-u55","Sbc11i4-u50","Sbc11i4-u51"] & lbls=["n0low","n0med","n0high"]

;fruns= ["Sbc11i4-u4","Sbc11i4-u43"] & lbls=["n2low","n2high"]
;fruns= ["Sbc11i4-u57","Sbc11i4-u4","Sbc11i4-u43"] & lbls=["n2low","n2med","n2high"]

;fruns= ["Sbc11i4-u4","Sbc11i4-u4a","Sbc11i4-u4i","Sbc11i4-u43"] & lbls=["zero (n2low)","solar","n2high"]
; -------------------------------------


;fruns= ["I1i-u16","I1i-u3","I1i-u17"] & lbls=["n1low","n1med","n1high"]
;fruns= ["I1i-u8","I1i-u2","I1i-u9"] & lbls=["n0low","n0med","n0high"]
;fruns= ["I1i-u6","I1i-u1","I1i-u7"] & lbls=["n2low","n2med","n2high"]

;fruns= ["I1gf1i-u16","I1gf1i-u3","I1gf1i-u17"] & lbls=["n1low","n1med","n1high"]
;fruns= ["I1gf1i-u8","I1gf1i-u2","I1gf1i-u9"] & lbls=["n0low","n0med","n0high"]
;fruns= ["I1gf1i-u6","I1gf1i-u1","I1gf1i-u7"] & lbls=["n2low","n2med","n2high"]

;fruns= ["I1gf2i-u16","I1gf2i-u3","I1gf2i-u17"] & lbls=["n1low","n1med","n1high"]
;fruns= ["I1gf2i-u8","I1gf2i-u2","I1gf2i-u9"] & lbls=["n0low","n0med","n0high"]
;fruns= ["I1gf2i-u6","I1gf2i-u1","I1gf2i-u7"] & lbls=["n2low","n2med","n2high"]

fruns= ["I1i-u4a","I1i-u4","I1i-u6","I1i-u1"]
lbls=["S00","S00, new SPH","!8n2low!6, new SPH","!8n2med!6, new SPH"]
;fruns= ["I1i-u4a","I1i-u4","I1i-u6","I1i-u1","I1i-u7"]
;lbls=["S00","S00, new SPH","!8n2low!6, new SPH","!8n2med!6, new SPH","!8n2high!6, new SPH"]


; -------------------------------------


zcolor= 50
deltacolor= (250-zcolor)/n_elements(fruns)
;deltacolor= 50
lstyle= 0
zcolor_orig= zcolor
thickstyle= 6.0

;--------------------------------------
;  Print the Shit
;--------------------------------------


xaxistitle = "!6Time (!8h!6!e-1!N Gyr)"
;xaxistitle = "!6Time (Gyr)"

xmax = 1.0
xmin = 0

;ymax = 5.0
ymax = 7.25
;ymax = 50.0
ymin = 0 
;ymin = 0.2


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
	ytitle="!6SFR (M!D!9n!6!N Yr!U-1!N)", $
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

	    oplot, sfrtime, sfrsfr, psym=-3, color= zcolor, linestyle= lstyle, thick= thickstyle
	    lstyle= lstyle+1
	    zcolor= zcolor+deltacolor
	    ;thickstyle= thickstyle-4.0
	endfor
endelse


xyouts, 0.55, 0.88, lbls[0], /normal, charthick=6.0, size=1.50, color= zcolor_orig
;xyouts, 0.7, 0.71, "n2low", /normal, charthick=6.0, size=1.50, color= zcolor_orig
;xyouts, 0.7, 0.76, "n1low", /normal, charthick=6.0, size=1.50, color= zcolor_orig
;xyouts, 0.7, 0.76, "n0low", /normal, charthick=6.0, size=1.50, color= zcolor_orig
;xyouts, 0.7, 0.55, "n0toolow", /normal, charthick=6.0, size=1.50, color= zcolor_orig

xyouts, 0.55, 0.83, lbls[1], /normal, charthick=6.0, size=1.50, color= zcolor_orig+deltacolor
;xyouts, 0.7, 0.47, "n2med", /normal, charthick=6.0, size=1.50, color= zcolor_orig+deltacolor
;xyouts, 0.7, 0.47, "n1med", /normal, charthick=6.0, size=1.50, color= zcolor_orig+deltacolor
;xyouts, 0.7, 0.47, "n0med", /normal, charthick=6.0, size=1.50, color= zcolor_orig+deltacolor

xyouts, 0.55, 0.78, lbls[2], /normal, charthick=5.0, size=1.50, color= zcolor_orig+deltacolor+deltacolor
;xyouts, 0.7, 0.27, "n2high", /normal, charthick=2.0, size=1.50, color= zcolor_orig+deltacolor+deltacolor
;xyouts, 0.7, 0.29, "n1high", /normal, charthick=2.0, size=1.50, color= zcolor_orig+deltacolor+deltacolor
;xyouts, 0.7, 0.29, "n0high", /normal, charthick=4.0, size=1.50, color= zcolor_orig+deltacolor+deltacolor
;xyouts, 0.7, 0.25, "n0low", /normal, charthick=6.0, size=1.50, color= zcolor_orig+deltacolor

xyouts, 0.55, 0.73, lbls[3], /normal, charthick=5.0, size=1.50, color= zcolor_orig+deltacolor+deltacolor+deltacolor

xyouts, 0.55, 0.68, lbls[4], /normal, charthick=5.0, size=1.50, color= zcolor_orig+deltacolor+deltacolor+deltacolor+deltacolor




;--------------------------------------
;--------------------------------------

if (sendto EQ 'ps') then device, /close


end


