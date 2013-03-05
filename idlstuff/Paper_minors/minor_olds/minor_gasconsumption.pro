
;--------------------------------------------------------------------------



pro calculate_burst_params, frun, MaxSFRTime= MaxSFRTime, AvgSFR= AvgSFR, $
	    MassFormedStars= MassFormedStars, MaximumSFR= MaximumSFR, $
            OrigGasMass= OrigGasMass, FinalGasMass= FinalGasMass, $
	    GasConsumption= GasConsumption, GasCFrac= GasCFrac


	; get sfr data
            sfrfile= '/data/tcox/sfr/'+fload_getid(frun)+'.sfr'

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




pro minor_sfinfo, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "minor_sfinfo, junk"
   return
endif



l1= 0
l2= 0
l3= 0
l4= 0
ymax= -1

; --------------------------------------------------
;
;     Bulge Size
;do_BtoD_ratio= 1
do_BtoD_ratio= 0
if do_BtoD_ratio eq 1 then begin
	xaxistitle= "Bulge-to-Disk Ratio"
	xmax = 1.1
	xmin =-0.05


	l1= 1
	msg1= 'Exponential Bulge'
	G3fruns= ["G3G3b-u1", "G3BT1G3BT1-u1", "G3BT2G3BT2-u1", "G3BT3G3BT3-u1"]
	xvar= [0.20, 0.25, 0.5, 1.0]
	grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
			MassFSlist, OrigGMlist, FinalGMlist, $
			GClist, GCFlist

	sfmax_1= MaxSFRlist
	gc_1= GCFlist


	l2= 1
	msg2= 'Hernquist Bulge'
	G3fruns= ["G3bG3b", "G3b6G3b6", "G3b5G3b5", "G3b4G3b4"]
	xvar= [0.20, 0.25, 0.5, 1.0]
	grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
			MassFSlist, OrigGMlist, FinalGMlist, $
			GClist, GCFlist

	sfmax_2= MaxSFRlist
	gc_2= GCFlist
endif




; --------------------------------------------------
;
;     Oribits
;
;do_orbit= 1
do_orbit= 0
if do_orbit eq 1 then begin
        xaxistitle= "Pericentric Distance"
        xmax = 68.0
        xmin = -2.0

        xvar= [0.01, 6.8, 13.6, 27.2, 64.4]  ; R_peri (same for both)

        l1= 1
        msg1= 'Mass ratio 1:2.3'
        G3fruns= ["G3G2h-u1", "G3G2e-u1", "G3G2-u3", "G3G2f-u1", "G3G2g-u1"]
        grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        GClist, GCFlist

        sfmax_1= MaxSFRlist
        gc_1= GCFlist


        l2= 1
        msg2= 'Mass ratio 1:5.8'
        G3fruns= ["G3G1h-u1", "G3G1e-u1", "G3G1-u3", "G3G1f-u1", "G3G1g-u1"]
        grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        GClist, GCFlist

        sfmax_2= MaxSFRlist
        gc_2= GCFlist
endif




; --------------------------------------------------
;
;     Orientations
;
;do_ori= 1
do_ori= 0
if do_ori eq 1 then begin
        xaxistitle= "Satellite Inclination (degrees)"
        xmax = 185.0
        xmin = -5.0

        xvar= [0.0, 30.0, 60.0, 90.0, 150.0, 180.0]  ; theta

        l1= 1 
        msg1= 'Mass ratio 1:2.3'
        G3fruns= ["G3G2a-u1", "G3G2-u3", "G3G2c-u1", "G3G2b-u1", "G3G2r-u3", "G3G2d-u1"]
        grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        GClist, GCFlist

        sfmax_1= MaxSFRlist
        gc_1= GCFlist


        l2= 1 
        msg2= 'Mass ratio 1:5.8'
        G3fruns= ["G3G1a-u1", "G3G1-u3", "G3G1c-u1", "G3G1b-u1", "G3G1r-u3", "G3G1d-u1"]
        grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        GClist, GCFlist

        sfmax_2= MaxSFRlist
        gc_2= GCFlist
endif




; --------------------------------------------------
;
;     Gas Fraction
;
do_gf= 1
;do_gf= 0
if do_gf eq 1 then begin
        xaxistitle= "Progenitor Gas Fraction"
        xmax = 1.05
        xmin = -0.05

	ymax= 100.0

        ;xvar= [.11, 0.2, 0.42, 0.58, 0.75]  ; theta
        xvar= [.75, 0.58, 0.42, 0.2, 0.11]  ; theta

        l1= 1 
        msg1= 'B/T= 0.2'
        ;G3fruns= ["G3gf4G3gf4", "G3G3b-u1", "G3gf1G3gf1b-u1", "G3gf2G3gf2b-u1", "G3gf3G3gf3"]
        G3fruns= ["G3gf3G3gf3", "G3gf2G3gf2b-u1", "G3gf1G3gf1b-u1", "G3G3b-u1", "G3gf4G3gf4"]
        grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        GClist, GCFlist

        sfmax_1= MaxSFRlist
        gc_1= GCFlist


        l2= 1 
        msg2= 'B/T= 0.5 (Hb)'
        G3fruns= ["G3b5_gf3_mm", "G3b5_gf2_mm", "G3b5_gf2_mm", "G3b5G3b5", "G3b5_gf4_mm"]
        grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        GClist, GCFlist

        sfmax_2= MaxSFRlist
        gc_2= GCFlist
endif





;============================================
;============================================
;
;     Now, we plot!!
;
;============================================
;============================================


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename='sfinfo.eps', colortable=4, newxsize= 11.0, newysize= 18.0


;
;    |-----------|
;    |           |
;    |  SFR_max  |
;    |           |
;    |           |
;    |-----------|
;    |           |
;    |  Gas      |
;    |  Consump. |
;    |           |
;    |-----------|
;
;
;
;
x0= 0.18
x1= 0.98

y0= 0.09
y1= 0.54
y2= 0.99


;--------------------------------------

;  x axis information is set above



;============================================
;
;     Panel 1
;
;============================================

yaxistitle= "SFR!Dmax!N (M!D!9n!6!N Yr!E-1!N)"
if ymax lt 0 then ymax = 26.0
ymin = 1.0

!p.position= [x0, y1, x1, y2]

;---------------------------

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, charthick=3.0, xthick=4.0, ythick=4.0, color= 0, $
	xtickformat='(a1)', ytitle=yaxistitle, /nodata

        if l1 ne 0 then oplot, xvar, sfmax_1, psym=-2, thick=3.0, color= 150, linestyle= 0
        if l2 ne 0 then oplot, xvar, sfmax_2, psym=-5, thick=3.0, color= 50, linestyle= 0
        if l3 ne 0 then oplot, xvar, sfmax_3, psym=-7, thick=3.0, color= 100, linestyle= 0
        if l4 ne 0 then oplot, xvar, sfmax_4, psym=-6, thick=3.0, color= 0, linestyle= 0



        xl0= 0.45
        yl0= 0.72
        if l1 ne 0 then xyouts, xl0, yl0, msg1, /normal, charthick=3.0, size=1.3, color= 150
        if l2 ne 0 then xyouts, xl0, yl0-0.05, msg2, /normal, charthick=3.0, size=1.3, color= 50
        if l3 ne 0 then xyouts, xl0, yl0-0.10, msg3, /normal, charthick=3.0, size=1.3, color= 100
        if l4 ne 0 then xyouts, xl0, yl0-0.15, msg4, /normal, charthick=3.0, size=1.3, color= 0
        if l1 ne 0 then oplot, [0.05], [0.88], psym=2, thick=3.0, symsize=1.5, color=150
        if l2 ne 0 then oplot, [0.05], [0.82], psym=5, thick=3.0, symsize=1.5, color=50
        if l3 ne 0 then oplot, [0.05], [0.76], psym=7, thick=3.0, symsize=1.5, color=100
        if l4 ne 0 then oplot, [0.05], [0.70], psym=6, thick=3.0, symsize=1.5, color= 0

        ;x0= 0.01 & xs= 0.21
        ;y0= 0.92 & ys= 0.25
        ;oplot, [x0,x0+xs,x0+xs,x0,x0],[y0,y0,y0-ys,y0-ys,y0], thick=3.0, color= 0, psym=-3





;============================================
;
;     Panel 1
;
;============================================

yaxistitle= "Gas Consumption Fraction"
ymax = 1.0
ymin = 0.0

!p.position= [x0, y0, x1, y1]

;---------------------------

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, charthick=3.0, xthick=4.0, ythick=4.0, color= 0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase


        if l1 ne 0 then oplot, xvar, gc_1, psym=-2, thick=3.0, color= 150, symsize=1.5
        if l2 ne 0 then oplot, xvar, gc_2, psym=-5, thick=3.0, color= 50, symsize=1.5
        if l3 ne 0 then oplot, xvar, gc_3, psym=-7, thick=3.0, color= 100, symsize=1.5
        if l4 ne 0 then oplot, xvar, gc_4, psym=-6, thick=3.0, color= 0, symsize=1.5




;============================================
device, /close



end








;====================================================================================
;
;====================================================================================







pro be, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "be, junk"
   return
endif



l1= 0
l2= 0
l3= 0
l4= 0
ymax= -1

; --------------------------------------------------
;
;     Bulge Size
;do_BtoD_ratio= 1
do_BtoD_ratio= 0
if do_BtoD_ratio eq 1 then begin
	xaxistitle= "Bulge-to-Disk Ratio"
	xmax = 1.1
	xmin =-0.05


	l1= 1
	msg1= 'Exponential Bulge'
	G3fruns= ["G3G3b-u1", "G3BT1G3BT1-u1", "G3BT2G3BT2-u1", "G3BT3G3BT3-u1"]
	xvar= [0.20, 0.25, 0.5, 1.0]
	grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
			MassFSlist, OrigGMlist, FinalGMlist, $
			GClist, GCFlist

	sfmax_1= MaxSFRlist
	gc_1= GCFlist


	l2= 1
	msg2= 'Hernquist Bulge'
	G3fruns= ["G3bG3b", "G3b6G3b6", "G3b5G3b5", "G3b4G3b4"]
	xvar= [0.20, 0.25, 0.5, 1.0]
	grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
			MassFSlist, OrigGMlist, FinalGMlist, $
			GClist, GCFlist

	sfmax_2= MaxSFRlist
	gc_2= GCFlist
endif




; --------------------------------------------------
;
;     Oribits
;
do_orbit= 1
;do_orbit= 0
if do_orbit eq 1 then begin
        xaxistitle= "Pericentric Distance"
        xmax = 68.0
        xmin = -2.0

        xvar= [0.01, 6.8, 13.6, 27.2, 64.4]  ; R_peri (same for both)

        l1= 1
        ;msg1= 'Mass ratio 1:2.3'
        msg1= 'G3G2, mass ratio 2.3:1'
        G3fruns= ["G3G2h-u1", "G3G2e-u1", "G3G2-u3", "G3G2f-u1", "G3G2g-u1"]
	isolated_gas_consumption= fltarr(n_elements(G3fruns))
	isolated_gas_consumption(*)= 0.232424
        grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        GClist, GCFlist

        sfmax_1= MaxSFRlist
        be_1= GCFlist - isolated_gas_consumption


        l2= 1
        ;msg2= 'Mass ratio 1:5.8'
        msg2= 'G3G1, mass ratio 5.8:1'
        G3fruns= ["G3G1h-u1", "G3G1e-u1", "G3G1-u3", "G3G1f-u1", "G3G1g-u1"]
	isolated_gas_consumption= fltarr(n_elements(G3fruns))
	isolated_gas_consumption(*)= 0.219662
        grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        GClist, GCFlist

        sfmax_2= MaxSFRlist
        be_2= GCFlist - isolated_gas_consumption
endif




; --------------------------------------------------
;
;     Orientations
;
;do_ori= 1
do_ori= 0
if do_ori eq 1 then begin
        xaxistitle= "Satellite Inclination (degrees)"
        xmax = 185.0
        xmin = -5.0

        xvar= [0.0, 30.0, 60.0, 90.0, 150.0, 180.0]  ; theta

        l1= 1 
        ;msg1= 'Mass ratio 1:2.3'
        msg1= 'G3G2, mass ratio 2.3:1'
        G3fruns= ["G3G2a-u1", "G3G2-u3", "G3G2c-u1", "G3G2b-u1", "G3G2r-u3", "G3G2d-u1"]
	isolated_gas_consumption= fltarr(n_elements(G3fruns))
	isolated_gas_consumption(*)= 0.232424
        grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        GClist, GCFlist

        sfmax_1= MaxSFRlist
        be_1= GCFlist - isolated_gas_consumption


        l2= 1 
        ;msg2= 'Mass ratio 1:5.8'
        msg2= 'G3G1, mass ratio 5.8:1'
        G3fruns= ["G3G1a-u1", "G3G1-u3", "G3G1c-u1", "G3G1b-u1", "G3G1r-u3", "G3G1d-u1"]
	isolated_gas_consumption= fltarr(n_elements(G3fruns))
	isolated_gas_consumption(*)= 0.219662
        grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        GClist, GCFlist

        sfmax_2= MaxSFRlist
        be_2= GCFlist - isolated_gas_consumption
endif




; --------------------------------------------------
;
;     Gas Fraction
;
;do_gf= 1
do_gf= 0
if do_gf eq 1 then begin
        xaxistitle= "Progenitor Gas Fraction"
        xmax = 1.05
        xmin = -0.05

	ymax= 100.0

        ;xvar= [.11, 0.2, 0.42, 0.58, 0.75]  ; theta
        xvar= [.75, 0.58, 0.42, 0.2, 0.11]  ; theta

        l1= 1 
        msg1= 'B/T= 0.2'
        ;G3fruns= ["G3gf4G3gf4", "G3G3b-u1", "G3gf1G3gf1b-u1", "G3gf2G3gf2b-u1", "G3gf3G3gf3"]
        G3fruns= ["G3gf3G3gf3", "G3gf2G3gf2b-u1", "G3gf1G3gf1b-u1", "G3G3b-u1", "G3gf4G3gf4"]
	isolated_gas_consumption= fltarr(n_elements(G3fruns))
	isolated_gas_consumption(*)= 0.232424 & print, "PROBLEM  PROBLEM  PROBLEM  PROBLEM "
        grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        GClist, GCFlist

        sfmax_1= MaxSFRlist
        be_1= GCFlist - isolated_gas_consumption


        l2= 1 
        msg2= 'B/T= 0.5 (Hb)'
        G3fruns= ["G3b5_gf3_mm", "G3b5_gf2_mm", "G3b5_gf2_mm", "G3b5G3b5", "G3b5_gf4_mm"]
	isolated_gas_consumption= fltarr(n_elements(G3fruns))
	isolated_gas_consumption(*)= 0.232424 & print, "PROBLEM  PROBLEM  PROBLEM  PROBLEM "
        grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        GClist, GCFlist

        sfmax_2= MaxSFRlist
        be_2= GCFlist - isolated_gas_consumption
endif





;============================================
;============================================
;
;     Now, we plot!!
;
;============================================
;============================================


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename='beinfo.eps', colortable=4


;
;    |-----------|
;    |           |
;    |  Gas      |
;    |  Consump. |
;    |           |
;    |-----------|
;
;
;
;
x0= 0.18
x1= 0.98

y0= 0.12
y1= 0.98


;--------------------------------------

;  x axis information is set above




;============================================
;
;
;============================================

yaxistitle= "burst efficiency"
ymax = 1.0
ymin = 0.0

!p.position= [x0, y0, x1, y1]

;---------------------------

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, charthick=3.0, xthick=4.0, ythick=4.0, color= 0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase


        if l1 ne 0 then oplot, xvar, be_1, psym=-2, thick=3.0, color= 150, symsize=1.5
        if l2 ne 0 then oplot, xvar, be_2, psym=-5, thick=3.0, color= 50, symsize=1.5
        if l3 ne 0 then oplot, xvar, be_3, psym=-7, thick=3.0, color= 100, symsize=1.5
        if l4 ne 0 then oplot, xvar, be_4, psym=-6, thick=3.0, color= 0, symsize=1.5


        xl0= 0.28
        yl0= 0.87
        if l1 ne 0 then xyouts, xl0, yl0, msg1, /normal, charthick=3.0, size=1.3, color= 150
        if l2 ne 0 then xyouts, xl0, yl0-0.05, msg2, /normal, charthick=3.0, size=1.3, color= 50
        if l3 ne 0 then xyouts, xl0, yl0-0.10, msg3, /normal, charthick=3.0, size=1.3, color= 100
        if l4 ne 0 then xyouts, xl0, yl0-0.15, msg4, /normal, charthick=3.0, size=1.3, color= 0
	; orientation
        ;if l1 ne 0 then oplot, [8.0], [0.885], psym=2, thick=3.0, symsize=1.5, color=150
        ;if l2 ne 0 then oplot, [8.0], [0.825], psym=5, thick=3.0, symsize=1.5, color=50
	; r_peri
        if l1 ne 0 then oplot, [3.0], [0.885], psym=2, thick=3.0, symsize=1.5, color=150
        if l2 ne 0 then oplot, [3.0], [0.825], psym=5, thick=3.0, symsize=1.5, color=50
        if l3 ne 0 then oplot, [8.0], [0.76], psym=7, thick=3.0, symsize=1.5, color=100
        if l4 ne 0 then oplot, [8.0], [0.70], psym=6, thick=3.0, symsize=1.5, color= 0


;============================================
device, /close



end










; ----------------------------------------------------------------
;
;
;
;   random information about the isolated disks
;
;
;-------------------------------------
;opening: /data/tcox/sfr/G3il-u1a.sfr
;-------------------------------------
;maximum sfr=       1.95668 occurs at    0.00514800
;average sfr=      0.485828
;  
;      tmerg=       1.80000  and 200 Myr later=       2.00062
;        sfr=      0.844958  200 Myr after t_merger
;         sm=      0.190481  200 Myr after t_merger
;-------------------------------------
;original gas mass   =      1.22000
;remnant gas mass    =     0.936298 at T=      6.00000
;gas consumed        =     0.283702
;                            23.2543  %
;                            7.54778 % after 1 Gyr
;                            22.2912 % after 4 Gyr
;-------------------------------------
;
;
;
;
;-------------------------------------
;opening: /data/tcox/sfr/G2im-u1a.sfr
;-------------------------------------
;maximum sfr=      0.500792 occurs at    0.00516200
;average sfr=      0.191099
;  
;      tmerg=       1.80000  and 200 Myr later=       2.00342
;        sfr=      0.297535  200 Myr after t_merger
;         sm=     0.0599360  200 Myr after t_merger
;-------------------------------------
;original gas mass   =     0.480000
;remnant gas mass    =     0.368579 at T=      6.00000
;gas consumed        =     0.111421
;                            23.2127  %
;                            5.10542 % after 1 Gyr
;                            21.0963 % after 4 Gyr
;-------------------------------------
;
;
;
;
;
;-------------------------------------
;opening: /data/tcox/sfr/G1i-u1a.sfr
;-------------------------------------
;maximum sfr=      0.152132 occurs at    0.00523900
;average sfr=     0.0487608
;  
;      tmerg=       1.80000  and 200 Myr later=       2.00043
;        sfr=     0.0955120  200 Myr after t_merger
;         sm=     0.0122640  200 Myr after t_merger
;-------------------------------------
;original gas mass   =     0.200000
;remnant gas mass    =     0.171754 at T=      6.00000
;gas consumed        =    0.0282460
;                            14.1230  %
;                            3.01100 % after 1 Gyr
;                            11.3940 % after 4 Gyr
;-------------------------------------
;
;
;
;
;
;-------------------------------------
;opening: /data/tcox/sfr/G0i-u1a.sfr
;-------------------------------------
;maximum sfr=     0.0246900 occurs at       4.90792
;average sfr=    0.00648689
;  
;      tmerg=       1.80000  and 200 Myr later=       2.00315
;        sfr=    0.00349800  200 Myr after t_merger
;         sm=   0.000416000  200 Myr after t_merger
;-------------------------------------
;original gas mass   =    0.0600000
;remnant gas mass    =    0.0562320 at T=      6.00000
;gas consumed        =   0.00376800
;                            6.27999  %
;                           0.231663 % after 1 Gyr
;                            2.91833 % after 4 Gyr
;-------------------------------------
;
;
;
;
;
