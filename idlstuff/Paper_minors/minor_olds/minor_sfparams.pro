pro minor_sfp, junk, $
		filename=filename, $
		gasmass= gasmass, $
		starmass=starmass, $
		normalize=normalize, $
		h=h

if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_multi, junk, filename=filename, h=0.7, /starmass, /gasmass, /normalize"
   print, "  "
   print, "  "
   return
endif


if not keyword_set(filename) then filename='minorsfp1.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0



;-------------------------------------------
;  set variables
;-------------------------------------------

frun1= "G3G3b-u1"
frun2= "G3il-u1a"


;  G Models

;fruns= ["G2i-u1", "G2in-u1", "G2im-u1"]
;fruns= ["G0i-u1", "G1i-u1", "G2im-u1", "G3il-u1"] & lbls= ['G0','G1','G2','G3']
;fruns= ["G0i-u1a", "G1i-u1a", "G2im-u1a", "G3il-u1a"] & lbls= ['G0','G1','G2','G3']
;fruns= ["G3G0-u1", "G3G0a-u1", "G3G0b-u1", "G3G0c-u3"]
;fruns= ["G3G3a-u1", "G3G3b-u1"]

; prograde
;fruns= ["G3G3b-u1", "G3G2-u3", "G3G1-u3", "G3G0e-u3"] & lbls= ['G3G3 (2.5)','G3G2 (2.9)','G3G1 (4.2)','G3G0 (6.0)']
;fruns= ["G2G2-u1", "G2G1-u3", "G2G0-u3"] & lbls= ['G2G2 (1.3)','G2G1 (1.5)','G2G0 (2.7)']
;fruns= ["G1G1a-u1", "G1G0-u3"] & lbls= ['G1G1 (1.3)','G1G0 (2.1)']
;fruns=  ["G0G0a-u1"] & lbls= ['G0G0 (1.4)']

; smaller majors with lower smoothing length
;fruns=  ["G0G0a-u1","G0G0a-u3"] & lbls= ['G0G0 (1.4)','G0G0 (h=0.035)']
;fruns=  ["G1G1a-u1","G1G1a-u3"] & lbls= ['G1G1 (1.3)','G1G1 (h=0.055)']

; retrograte
;fruns= ["G3G3r-u1","G3G2r-u3","G3G1r-u3"] & lbls= ['G3G3r','G3G2r','G3G1r']
;fruns= ["G2G2r-u1","G2G1r-u3","G2G0r-u3"] & lbls= ['G2G2r','G2G1r','G2G0r']
;fruns= ["G1G1r-u1","G1G0r-u3"] & lbls= ['G1G1r','G1G0r']
;fruns=  ["G0G0r-u1b"] & lbls= ['G0G0r']


;fruns= ["G3il-u1","G3bli-u1","G3bliv2-u1","G3bliv3-u1","G3bliv4-u1","G3bliv5-u1","G3bliv6-u1"]
;fruns= ["G3il-u1","G3bliv5-u1"]
;fruns= ["G3G3b-u1","G3blG3bl-u1"]
;fruns= ["G3G3b-u1","G3blv5G3blv5-u1"]


; just majors
;fruns= ["G3G3a-u1", "G3G3b-u1","G3G3r-u1"] & tmerg=[3.4,2.5,2.6] & lbls=['G3a','G3b','G3r']
;fruns= ["G2G2-u1","G2G2r-u1"] & tmerg=[1.3,1.3] & lbls=['G2','G2r']
;fruns= ["G1G1a-u1","G1G1r-u1"] & tmerg=[1.3,1.3] & lbls=['G1','G1r']
;fruns= ["G0G0a-u1","G0G0r-u1b"] & tmerg=[1.4,1.4] & lbls=['G0','G0r']
;fruns= ["G3blv5G3blv5-u1"] & tmerg= [2.4,2.4] & lblbs=['G3bl']

; ------------------------------------------------------

; vary gas fraction
;fruns= ["G3G3b-u1","G3gf1G3gf1b-u1","G3gf2G3gf2b-u1"]
;lbls= ['20%','50%','70%']
;fruns= ["G3gf4G3gf4","G3G3b-u1","G3gf1G3gf1b-u1","G3gf2G3gf2b-u1","G3gf3G3gf3"]
;lbls= ['11%','20%','42%','58%','75%']

; vary bulge fraction
;fruns= ["G3G3b-u1","G3BT1G3BT1-u1","G3BT2G3BT2-u1","G3BT3G3BT3-u1"]
;fruns= ["G3G3b-u1","G3blv5G3blv5-u1","G3BT1G3BT1-u1","G3BT2G3BT2-u1","G3BT3G3BT3-u1"]
;fruns= ["G3BT1G3BT1-u1","G3BT1G2-u1","G3BT1G1-u1","G3BT1G0-u1"]
;fruns= ["G3BT2G3BT2-u1","G3BT2G2-u1","G3BT2G1-u1","G3BT2G0-u1"]
;fruns= ["G3BT3G3BT3-u1","G3BT3G2-u1","G3BT3G1-u1","G3BT3G0-u1"]

;fruns= ["G3il-u1a","G3bliv5-u1","G3BT1i-u1","G3BT2i-u1","G3BT3i-u1"]

; minors with different n (feedback n)
;fruns= ["G3G3b-u1","G3G1-u3","G3G1-u4"]
;fruns= ["G3G3b-u1","G3G2-u3","G3G2-u4"]
;fruns= ["G3G3b-u1","G3G3b-u2","G3G0e-u3","G3G0e-u4"]

;fruns= ["G3G3b-u2", "G3G2-u4", "G3G1-u4", "G3G0e-u4"]
;fruns= ["G2G2-u2", "G2G1-u4", "G2G0-u4"]
;fruns= ["G1G1a-u2", "G1G0-u4"]
;fruns=  ["G0G0a-u2"]

;fruns= ["G2G2-u1", "G2G0-u3", "G2G0-u4"]

;fruns= ["G3gf1G3gf1b-u1","G3gf1G2-u1","G3gf1G1-u1","G3gf1G0-u1"]
;fruns= ["G3gf2G3gf2b-u1","G3gf2G2-u1","G3gf2G1-u1","G3gf2G0-u1"]

;fruns= ["G3G2a-u1","G3G2-u3","G3G2c-u1","G3G2b-u1","G3G2r-u3","G3G2d-u1"]
;fruns= ["G3G2e-u1","G3G2-u3","G3G2f-u1","G3G2g-u1","G3G2h-u1"]

;fruns= ["G3G1a-u1","G3G1-u3","G3G1c-u1","G3G1b-u1","G3G1r-u3","G3G1d-u1"]


; Mako & UpsAnd comparison
;fruns= ["G3G3b-u1", "G3blv5G3blv5-u1","G3bG3b"]
;fruns= ["G3BT1G3BT1-u1", "G3b1G3b1"]
;fruns= ["G3BT2G3BT2-u1", "G3b2G3b2"]
;fruns= ["G3BT3G3BT3-u1", "G3b3G3b3"]

; Bulge-to-disk ratios
;    -- exp. bulge (upsand) ---
;fruns= ["G3blv5G3blv5-u1", "G3BT1G3BT1-u1","G3BT2G3BT2-u1","G3BT3G3BT3-u1"]
;    -- exp. bulge (mako) ---
;fruns= ["G3bG3b", "G3b1G3b1","G3b2G3b2","G3b3G3b3"]
;    -- hernquist bulge (mako) ---
;fruns= ["G3bG3b", "G3b6G3b6","G3b5G3b5","G3b4G3b4"]


; Bulge-to-disk ratios  & gas fraction
;fruns= ["G3b5_gf3_mm","G3b5_gf2_mm","G3b5_gf1_mm","G3b5G3b5", "G3b5_gf4_mm"]
;fruns= ["G3b5_gf3_mm","G3b5_gf2_mm","G3b5G3b5", "G3b5_gf4_mm"]





; ---------------------------------------------------
; ---------------------------------------------------
;
;  Print the Shit
;
; ---------------------------------------------------
; ---------------------------------------------------


yaxistitle="!6SFR (M!D!9n!6!N Yr!E-1!N)"
if keyword_set(starmass) then yaxistitle="!6Mass (10!e10!n M!D!9n!6!N)"
if keyword_set(gasmass) then yaxistitle="!6Mass (10!e10!n M!D!9n!6!N)"
if keyword_set(normalize) then yaxistitle="!6Gas Mass / M!D0!N "
;xaxistitle = "Time (Gyr/h)"
xaxistitle = "!6Time (Gyr)"

;xmax = 6.0
xmax = 5.0
;xmax = 4.0
;xmax= 2.0
xmin = 0

;ymax = 60.0
;ymax = 15
;ymax = 12
;ymax = 5.0
;ymax = 2.0
;ymax = 1.2
ymax = 1.0
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
	;  Merger - frun1  (should be larger)
	;
	;-------------------------------------------
        open_sfr_file, frun1, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass

	xvar1= transpose(sfrtime)
	yvar1= transpose(sfrsfr)
        if keyword_set(starmass) then yvar1= transpose(sfrmfs)
	if keyword_set(gasmass) then yvar1= transpose(sfrgasmass)
	if keyword_set(normalize) then yvar1=yvar1/yvar1(0)
 
        oplot, xvar1, yvar1, psym=-3, color= 150, linestyle= 0, thick= 6.0




        ;-------------------------------------------
        ;
        ;  Isolated - frun2  (should be smaller)
        ;
        ;-------------------------------------------
        open_sfr_file, frun2, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass

        xvar2= transpose(sfrtime)
        yvar2= transpose(sfrsfr)
        if keyword_set(starmass) then yvar2= transpose(sfrmfs)
        if keyword_set(gasmass) then yvar2= transpose(sfrgasmass)
	yvar2= 2.0 * yvar2
        if keyword_set(normalize) then yvar2=yvar2/yvar2(0)

        oplot, xvar2, yvar2, psym=-3, color= 0, linestyle= 0, thick= 2.0

        oplot, xvar2, yvar2, psym=-3, color= 0, linestyle= 0, thick= 2.0




        ;-------------------------------------------
        ;
        ;  Shade in region
        ;
        ;-------------------------------------------
	idx= where(xvar1 le xmax)
	xvar1= xvar1(idx)
	yvar1= yvar1(idx)

	idx= where(xvar2 le xmax)
	xvar2= xvar2(idx)
	yvar2= yvar2(idx)

	y20= yvar2[0]
	sort2= bsort(xvar2, /reverse)
	xvar2_rev= xvar2(sort2)
	yvar2_rev= yvar2(sort2)
	
	xrng= [ xmin, xvar1, xvar2_rev, xmin]
	yrng= [  y20, yvar1, yvar2_rev,  y20] 

	polyfill, xrng, yrng, /data, color= 250, /fill


	plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, charthick=3.0, xthick=4.0, ythick=4.0, $
	xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase



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

                ;if strmid(frun,1,4) eq 'raid' then begin
                ;        spawn, '/bin/ls '+frun+'/*.sfr ',result
                ;endif else begin
                ;        case loopnum of
                ;          0: datadir='/data/tcox/sfr'
                ;          1: datadir='/home'
                ;          2: datadir='/raid2'
                ;          3: datadir='/data'
                ;          4: datadir='/data6'
                ;          5: datadir='/data7'
                ;          else: break
                ;        endcase

                        ;sfrfile= datadir+'/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'

                ;        nfrun= fload_getid(frun)
                ;        spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/*.sfr ',result
                ;endelse
		;
                ;sfrfile=strcompress(result[0],/remove_all)
		;

		sfrfile= '/data/tcox/sfr/'+fload_getid(frun)+'.sfr'



                openr, 1, sfrfile, ERROR=err
                close, 1

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






