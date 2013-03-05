pro sfcomp, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "sfcomp, junk"
   return
endif


;============================================
;============================================
;
;     Star Formation Comparison
;
;============================================
;============================================

;--------------------------------------
;  Print the Shit
;--------------------------------------

ymax = 1.0
ymin = 0.0
xmax = 1.0
xmin = 0.0 

;---------------------------

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename='min_comp1.eps', colortable=4

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, charthick=3.0, xthick=4.0, ythick=4.0, color= 0, $
        ytitle="!6gas consumption (!8n0med!6)", $
        xtitle="!6gas consumption (!8n2med!6)", $
        /nodata

	; Do G3, G2, G1, and G0
	;G3's
	G3fruns= ['G3G3b-u1', 'G3G2-u3', 'G3G1-u3', 'G3G0e-u3']
	grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, j7, GCFlistn2
	G3fruns= ['G3G3b-u2', 'G3G2-u4', 'G3G1-u4', 'G3G0e-u4']
	grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, j7, GCFlistn0
        oplot, GCFlistn2, GCFlistn0, psym=2, symsize=1.5, thick=3.0, color= 150

	G2fruns= ['G2G2-u1', 'G2G1-u3', 'G2G0-u3']
	grab_sf_info, G2fruns, j1, j2, j3, j4, j5, j6, j7, GCFlistn2
	G2fruns= ['G2G2-u2', 'G2G1-u4', 'G2G0-u4']
	grab_sf_info, G2fruns, j1, j2, j3, j4, j5, j6, j7, GCFlistn0
        oplot, GCFlistn2, GCFlistn0, psym=5, symsize=1.5, thick=3.0, color= 100

	G1fruns= ['G1G1a-u1', 'G1G0-u3']
	grab_sf_info, G1fruns, j1, j2, j3, j4, j5, j6, j7, GCFlistn2
	G1fruns= ['G1G1a-u2', 'G1G0-u4']
	grab_sf_info, G1fruns, j1, j2, j3, j4, j5, j6, j7, GCFlistn0
        oplot, GCFlistn2, GCFlistn0, psym=7, symsize=1.5, thick=3.0, color= 50

	G0fruns= ['G0G0a-u1']
	grab_sf_info, G0fruns, j1, j2, j3, j4, j5, j6, j7, GCFlistn2
	G0fruns= ['G0G0a-u2']
	grab_sf_info, G0fruns, j1, j2, j3, j4, j5, j6, j7, GCFlistn0
        oplot, GCFlistn2, GCFlistn0, psym=6, symsize=1.5, thick=3.0, color= 0




	; Key for G3, G2, G1, and G0
	; --------------------------------
        x0= 0.29
        y0= 0.85
        xyouts, x0, y0, 'G3Gx', /normal, charthick=3.0, size=1.3, color= 150
        xyouts, x0, y0-0.05, 'G2Gx', /normal, charthick=3.0, size=1.3, color= 100
        xyouts, x0, y0-0.10, 'G1Gx', /normal, charthick=3.0, size=1.3, color= 50
        xyouts, x0, y0-0.15, 'G0G0', /normal, charthick=3.0, size=1.3, color= 0
        oplot, [0.08], [0.89], psym=2, thick=3.0, symsize=1.5, color=150
        oplot, [0.08], [0.83], psym=5, thick=3.0, symsize=1.5, color=100
        oplot, [0.08], [0.765], psym=7, thick=3.0, symsize=1.5, color=50
        oplot, [0.08], [0.70], psym=6, thick=3.0, symsize=1.5, color= 0

        ;x0= 0.03 & xs= 0.25
        ;y0= 0.92 & ys= 0.25
        ;oplot, [x0,x0+xs,x0+xs,x0,x0],[y0,y0,y0-ys,y0-ys,y0], thick=3.0, color= 0, psym=-3


	; draw equivalency line
	; ------------------------
	oplot, [0,1], [0,1], psym=-3, thick=2.0, color= 0, linestyle= 1




	; isolated comparisons
	; -----------------------
	Ifruns= ['G3il-u1a','G2im-u1a','G1i-u1a','G0i-u1a']
        grab_sf_info, Ifruns, j1, j2, j3, j4, j5, j6, j7, GCFlistn2
        Ifruns= ['G3_n0','G2im-u2','G1i-u2','G0i-u2']
        grab_sf_info, Ifruns, j1, j2, j3, j4, j5, j6, j7, GCFlistn0
        oplot, GCFlistn2, GCFlistn0, psym=4, symsize=1.5, thick=3.0, color= 185

	xyouts, x0, y0-0.20, 'isolated', /normal, charthick=3.0, size=1.3, color= 185
        oplot, [0.08], [0.64], psym=4, thick=3.0, symsize=1.5, color=185




device, /close


;============================================



end





;--------------------------------------------------------------------------



pro calculate_burst_params, frun, MaxSFRTime= MaxSFRTime, AvgSFR= AvgSFR, $
	    MassFormedStars= MassFormedStars, MaximumSFR= MaximumSFR, $
            OrigGasMass= OrigGasMass, FinalGasMass= FinalGasMass, $
	    GasConsumption= GasConsumption, GasCFrac= GasCFrac


	; get sfr data
            ;sfrfile= '/data/tcox/sfr/'+fload_getid(frun)+'.sfr'
            sfrfile= '/home/tcox/upsand/sfr/'+fload_getid(frun)+'.sfr'

            openr, 1, sfrfile, ERROR=err

            if (err NE 0) then begin
                print, "  "
                print, "Problem: ",!ERR_STRING
                print, "  " 
                close, 1 
                ERR= 0
                sfrtime= [0]
                sfrsfr= [0] 
                sfrmsg= "No Star Formation"
            endif else begin
                close, 1
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


            sfraftertmerg= 0.0
            smaftertmerg= 0.0
            taftertmerg= 0.0
            tmerger= 0.0
            ;if n_elements(tmerg) gt 1 then begin
            ;    idx= where(sfrtime ge tmerg[i]+0.2)
            ;    if idx(0) eq -1 then begin
            ;            sfraftertmerg= sfrsfr[n_rows-1]
            ;            smaftertmerg= sfrmfs[n_rows-1]
            ;            tmerger= sfrtime[n_rows-1]
            ;            taftertmerg= tmerger
            ;    endif else begin
            ;            sfraftertmerg= sfrsfr[idx(0)]
            ;            smaftertmerg= sfrmfs[idx(0)]
            ;            tmerger= tmerg[i]
            ;            taftertmerg= sfrtime[idx(0)]
            ;    endelse
            ;endif


            MaximumSFR= max(sfrsfr)
            idx= where(sfrsfr eq MaximumSFR)
	    MaxSFRTime= sfrtime(idx)
	    AvgSFR= total(sfrsfr)/n_rows
            n_mfs= n_elements(sfrgasmass)
	    MassFormedStars= sfrmfs[n_mfs-1]
            print, "-------------------------------------"
            print, "maximum sfr= ", MaximumSFR, " occurs at ", MaxSFRTime
            print, "average sfr= ", AvgSFR
            print, "mass of new stars formed= ", MassFormedStars

            OrigGasMass= sfrgasmass[0]
            FinalGasMass= sfrgasmass[n_mfs-1]
	    GasConsumption= OrigGasMass - FinalGasMass
	    GasCFrac= GasConsumption / OrigGasMass
            print, "-------------------------------------"
            print, "original gas mass   =", OrigGasMass
            print, "remnant gas mass    =", FinalGasMass
            print, "gas consumed        =", GasConsumption
            print, "gas consump f       =", GasCFrac
            print, "-------------------------------------"




end



;--------------------------------------------------------------------------




;--------------------------------------------------------------------------


pro grab_sf_info, fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        CGlist, GCFlist

   numfrun= n_elements(fruns)

   MaxSFRlist= fltarr(numfrun)
   MaxSFRTimelist= fltarr(numfrun)
   AvgSFRlist= fltarr(numfrun)
   MassFSlist= fltarr(numfrun)
   OrigGMlist= fltarr(numfrun)
   FinalGMlist= fltarr(numfrun)
   CGlist= fltarr(numfrun)
   GCFlist= fltarr(numfrun)

   for i=0, numfrun-1 do begin

        calculate_burst_params, fruns[i], MaxSFRTime= MaxSFRTime, AvgSFR= AvgSFR, $
            MassFormedStars= MassFormedStars, $
            OrigGasMass= OrigGasMass, FinalGasMass= FinalGasMass, $
            MaximumSFR= MaximumSFR, $
            GasConsumption= GasConsumption, GasCFrac= GasCFrac

        MaxSFRlist[i]= MaximumSFR
        MaxSFRTimelist[i]= MaxSFRTime
        AvgSFRlist[i]= AvgSFR
        MassFSlist[i]= MassFormedStars
        OrigGMlist[i]= OrigGasMass
        FinalGMlist[i]= FinalGasMass
        CGlist[i]= GasConsumption
        GCFlist[i]= GasCFrac
   endfor

end







;--------------------------------------------------------------------------



