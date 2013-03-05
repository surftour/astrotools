;==================================================================================
;
;
;    Compute the gas consumption and fit to it
;
;
;
;==================================================================================
pro calculate_burst_params, frun, MaxSFRTime= MaxSFRTime, AvgSFR= AvgSFR, $
	    MassFormedStars= MassFormedStars, MaximumSFR= MaximumSFR, $
            OrigGasMass= OrigGasMass, FinalGasMass= FinalGasMass, $
	    GasConsumption= GasConsumption, GasCFrac= GasCFrac


	; get sfr data
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


pro grab_sf_info, fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
			MassFSlist, OrigGMlist, FinalGMlist, $
			GClist, GCFlist

   numfrun= n_elements(fruns)

   MaxSFRlist= fltarr(numfrun)
   MaxSFRTimelist= fltarr(numfrun)
   AvgSFRlist= fltarr(numfrun)
   MassFSlist= fltarr(numfrun)
   OrigGMlist= fltarr(numfrun)
   FinalGMlist= fltarr(numfrun)
   GClist= fltarr(numfrun)
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
	GClist[i]= GasConsumption
	GCFlist[i]= GasCFrac
   endfor

end




;--------------------------------------------------------------------------


function perform_sbefficiency_fit, mratios, burstes


mratios_tofit= mratios
burstes_tofit= burstes
weight= mratios_tofit*0.0 + 1.0
p1= [1.0,1.0]
;spec_sf_eff= MPFITFUN('minor_func_ssflaw', mratios_tofit, burstes_tofit, weight, p1, BESTNORM=bestnorm, DOF=dof)
burst_eff_params= MPFITFUN('minor_func_samsbl', mratios_tofit, burstes_tofit, weight, p1, BESTNORM=bestnorm, DOF=dof)

redchi2= bestnorm/dof
print, "Reduced Chi^2 = ", redchi2
print, burst_eff_params

return, burst_eff_params


end



;--------------------------------------------------------------------------



pro minor_sf, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "minor_sfr, junk"
   return
endif


; this is more max sfr's than average
;isodata=    [[0.08, 2.0], $    ; G3
;		[0.05, 0.5],  $    ; G2
;		[0.02, 0.2],  $    ; G1
;		[0.005, 0.01]]     ; G0
; this is average sfr, but the efficiency is for 1 gyr
;isodata=    [[0.08, 0.95], $    ; G3
;               [0.05, 0.25],  $    ; G2
;               [0.03, 0.06],  $    ; G1
;               [0.002, 0.001]]     ; G0
; this is efficiency at 6 (G3) and 4 (the rest) Gyr,
;   and an average sfr
isodata=    [[0.23, 0.95], $    ; G3
               [0.21, 0.25],  $    ; G2
               [0.11, 0.06],  $    ; G1
               [0.029, 0.001]]     ; G0

Isofruns= ["G3il-u1", "G2im-u1", "G1i-u1", "G0i-u1"]




; total mass ratio,
;  stel mass ratio,
;  bar  mass ratio,
;  burst efficiency e,
;  sfr 200 Myr after merger
;  stellar mass at same time


; G3-G3 (1.1x10^12) Major Mergers
G3data=[[1.0, 1.0, 1.0, 0.65, 12.0, 10.85],  $      ; G3G3b-u1
        [2.3, 3.3, 3.1, 0.55, 4.5, 7.098],  $       ; G3G2-u3
        [5.8, 10.0, 8.9,  0.31, 0.65, 5.877],  $     ; G3G1-u3
	[22.7, 50.0, 38.9, 0.23, 0.06, 5.306]]       ; G3G0c-u3

G3fruns= ["G3G3b-u1", "G3G2-u3", "G3G1-u3", "G3G0c-u3"]
grab_sf_info, G3fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
			MassFSlist, OrigGMlist, FinalGMlist, $
			CGlist, GCFlist


G3rdata=[[1.0, 1.0, 1.0, 0.65, 9.0, 10.6604],  $      ; G3G3r-u1
        [2.3, 3.3, 3.1, 0.47, 4.4, 6.98783],  $       ; G3G2r-u3
        [5.8, 10.0, 8.9,  0.29, 2.1, 5.66511],  $     ; G3G1r-u3
        [22.7, 50.0, 38.9, 0.23, 0.05, 5.3]]       ; G3G0r-u3

G3bldata=[[1.0, 1.0, 1.0, 0.69, 14.7, 9.23256],  $      ; G3blG3bl-u1
        [2.3, 3.3, 3.1, 0.55, 6.9, 6.2006],  $       ; G3blG2-u3
        [5.8, 10.0, 8.9,  0.38, 2.5, 4.70153],  $     ; G3blG1-u3
        [22.7, 50.0, 38.9, 0.00, 0.0, 0.00000]]       ; G3blG0c-u3



; G2-G2 (5x10^11) Major Mergers
G2data=[[1.0, 1.0, 1.0, 0.73, 7.3, 3.32],  $      ; G2G2-u1
	[2.6, 3.0, 2.8, 0.63, 3.8, 2.191],   $      ; G2G1-u3
	[10.0, 15.0, 12.4, 0.28, 0.5, 1.693]]       ; G2G0-u3

G2rdata=[[1.0, 1.0, 1.0, 0.66, 5.8, 3.14579],  $      ; G2G2r-u1
        [2.6, 3.0, 2.8, 0.49, 2.6, 2.17141],   $      ; G2G1r-u3
        [10.0, 15.0, 12.4, 0.28, 0.7, 1.67446]]       ; G2G0r-u3



; G1-G1 (2x10^11) Major Mergers
G1data=[[1.0, 1.0, 1.0, 0.74, 2.60, 1.101],  $     ; G1G1a-u1
        [3.9, 5.0, 4.4, 0.34, 0.43, 0.633091]]        ; G1G0-u3

G1rdata=[[1.0, 1.0, 1.0, 0.66, 2.1, 1.10084],  $     ; G1G1r-u1
        [3.9, 5.0, 4.4, 0.30, 0.5, 0.642274]]        ; G1G0r-u3




; G1-G1 (2x10^11) Major Mergers
G0data=[[1.0, 1.0, 1.0, 0.67, 0.57, 0.2157]]        ; G0G0a-u1

G0rdata=[[1.0, 1.0, 1.0, 0.53, 0.7, 0.205655]]        ; G0G0r-u1b







;============================================
;============================================
;
;     Star Formation Efficiency
;
;============================================
;============================================

;--------------------------------------
;  Print the Shit
;--------------------------------------

ymax = 1.0
ymin = 0.0
xmax = 1.05
xmin =-0.05

;---------------------------

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename='min_sfr_gaseat.eps', colortable=4

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, charthick=3.0, xthick=4.0, ythick=4.0, color= 0, $
        xtitle="Mass Ratio (M!Dsat!N/M!Dprimary!N)", $
        ytitle="gas consumption", $
        /nodata

        oplot, 1.0/G3data[0,*], G3data[3,*], psym=-2, thick=3.0, color= 150, linestyle= 0
	oplot, 1.0/G2data[0,*], G2data[3,*], psym=-5, thick=3.0, color= 100, linestyle= 1
	oplot, 1.0/G1data[0,*], G1data[3,*], psym=-7, thick=3.0, color= 50, linestyle= 2
	oplot, 1.0/G0data[0,*], G0data[3,*], psym=-6, thick=3.0, color= 0



        x0= 0.30
        y0= 0.85
        xyouts, x0, y0, 'G3', /normal, charthick=3.0, size=1.3, color= 150
        xyouts, x0, y0-0.05, 'G2', /normal, charthick=3.0, size=1.3, color= 100
        xyouts, x0, y0-0.10, 'G1', /normal, charthick=3.0, size=1.3, color= 50
        xyouts, x0, y0-0.15, 'G0', /normal, charthick=3.0, size=1.3, color= 0
        oplot, [0.05], [0.88], psym=2, thick=3.0, symsize=1.5, color=150
        oplot, [0.05], [0.82], psym=5, thick=3.0, symsize=1.5, color=100
        oplot, [0.05], [0.76], psym=7, thick=3.0, symsize=1.5, color=50
        oplot, [0.05], [0.70], psym=6, thick=3.0, symsize=1.5, color= 0

        x0= 0.01 & xs= 0.21
        y0= 0.92 & ys= 0.25
        oplot, [x0,x0+xs,x0+xs,x0,x0],[y0,y0,y0-ys,y0-ys,y0], thick=3.0, color= 0, psym=-3



device, /close


;============================================

G3mr= 1.0/G3data[0,*]
G2mr= 1.0/G2data[0,*]
G1mr= 1.0/G1data[0,*]
G0mr= 1.0/G0data[0,*]

G3be= G3data[3,*]-isodata[0,0]
G2be= G2data[3,*]-isodata[0,1]
G1be= G1data[3,*]-isodata[0,2]
G0be= G0data[3,*]-isodata[0,3]

Gmasses= [G3data[5,0], G2data[5,0], G1data[5,0],G0data[5,0]]
Gbes= [max(G3be),max(G2be),max(G1be),max(G0be)]


; let's do some fitting



; first we'll fit major merger be
print, "first we'll fit the major merger burst efficiencies"
; ---------------------------------------------
weight= Gmasses*0.0 + 1.0
p1= [1.0,1.0]
;spec_sf_eff= MPFITFUN('minor_func_ssflaw', Gmasses, Gbes, weight, p1, PARINFO=pi, BESTNORM=bestnorm, DOF=dof)
spec_sf_eff= MPFITFUN('minor_func_ssflaw', Gmasses, Gbes, weight, p1, BESTNORM=bestnorm, DOF=dof)

redchi2= bestnorm/dof
print, "Reduced Chi^2 = ", redchi2
print, spec_sf_eff



; next we'll fit the burst efficiency
; ------------------------------------


; do them all at once
print, "fit all 10 simulations burst efficiency vs. merer ratio"
Gmr= [transpose(G3mr[0,*]),transpose(G2mr[0,*]),transpose(G1mr[0,*]),transpose(G0mr[0,*])]
Gbe= [transpose(G3be[0,*]),transpose(G2be[0,*]),transpose(G1be[0,*]),transpose(G0be[0,*])]
weight= Gmr*0.0 + 1.0
p1= [1.0,1.0]
G_result= MPFITFUN('minor_func_samsbl', Gmr, Gbe, weight, p1, PARINFO=pi, BESTNORM=bestnorm, DOF=dof)
redchi2= bestnorm/dof
print, "Reduced Chi^2 = ", redchi2
print, G_result



stop


; do them one at a time
; fixes the p[1] 
;pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},2)
;pi[1].fixed = 1


p1= [1.0,Gbes[0]]
weight= G3mr*0.0 + 1.0
G3be_result = MPFITFUN('minor_func_samsbl', G3mr, G3be, weight, p1, PARINFO=pi, BESTNORM=bestnorm, DOF=dof)
redchi2= bestnorm/dof
print, "Reduced Chi^2 = ", redchi2
print, G3be_result


stop


;p1= [1.0,Gbes[1]]
;weight= G2mr*0.0 + 1.0
;G2be_result = MPFITFUN('minor_func_samsbl', G2mr, G2be, weight, p1, PARINFO=pi, BESTNORM=bestnorm, DOF=dof)
;redchi2= bestnorm/dof
;print, "Reduced Chi^2 = ", redchi2
;print, G2be_result
;

;p1= [1.0,Gbes[2]]
;weight= G1mr*0.0 + 1.0
;G1be_result = MPFITFUN('minor_func_samsbl', G1mr, G1be, weight, p1, PARINFO=pi, BESTNORM=bestnorm, DOF=dof)
;redchi2= bestnorm/dof
;print, "Reduced Chi^2 = ", redchi2
;print, G1be_result











initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename='sfr_minne.eps', colortable=4

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, charthick=3.0, xthick=4.0, ythick=4.0, color= 0, $
        xtitle="Mass Ratio (M!Dsat!N/M!Dprimary!N)", $
        ytitle="burst efficiency", $
        /nodata

        ;oplot, G3mr, G3be, psym=-2, thick=3.0, color= 150, linestyle= 0
        ;oplot, G2mr, G2be, psym=-5, thick=3.0, color= 100, linestyle= 1
        ;oplot, G1mr, G1be, psym=-7, thick=3.0, color= 50, linestyle= 2
        oplot, G3mr, G3be, psym=2, thick=3.0, color= 150, symsize=1.5
        oplot, G2mr, G2be, psym=5, thick=3.0, color= 100, symsize=1.5
        oplot, G1mr, G1be, psym=7, thick=3.0, color= 50, symsize=1.5
        oplot, G0mr, G0be, psym=6, thick=3.0, color= 0, symsize=1.5


        x0= 0.30
        y0= 0.85
        xyouts, x0, y0, 'G3', /normal, charthick=3.0, size=1.3, color= 150
        xyouts, x0, y0-0.05, 'G2', /normal, charthick=3.0, size=1.3, color= 100
        xyouts, x0, y0-0.10, 'G1', /normal, charthick=3.0, size=1.3, color= 50
	xyouts, x0, y0-0.15, 'G0', /normal, charthick=3.0, size=1.3, color= 0
	oplot, [0.05], [0.89], psym=2, thick=3.0, symsize=1.5, color=150
	oplot, [0.05], [0.83], psym=5, thick=3.0, symsize=1.5, color=100
	oplot, [0.05], [0.77], psym=7, thick=3.0, symsize=1.5, color=50
	oplot, [0.05], [0.71], psym=6, thick=3.0, symsize=1.5, color= 0

	x0= 0.01 & xs= 0.21
	y0= 0.92 & ys= 0.25
	oplot, [x0,x0+xs,x0+xs,x0,x0],[y0,y0,y0-ys,y0-ys,y0], thick=3.0, color= 0, psym=-3


	; plot fits
	;------------
	x= findgen(21)/20.0
	;y= minor_func_samsbl(x,G3be_result)
	;oplot, x, y, psym=-3, color= 0, linestyle= 0

        ;y= minor_func_samsbl(x,G2be_result)
        ;oplot, x, y, psym=-3, color= 0, linestyle= 0

        ;y= minor_func_samsbl(x,G1be_result)
        ;oplot, x, y, psym=-3, color= 0, linestyle= 0

	y= minor_func_samsbl(x,G_result)
	oplot, x, y, psym=-3, color= 0, thick=3.0, linestyle= 2


device, /close


;============================================


G3mr= 1.0/G3data[0,*]
G2mr= 1.0/G2data[0,*]
G1mr= 1.0/G1data[0,*]
G0mr= 1.0/G0data[0,*]

G3rbe= G3rdata[3,*]-isodata[0,0]
G2rbe= G2rdata[3,*]-isodata[0,1]
G1rbe= G1rdata[3,*]-isodata[0,2]
G0rbe= G0rdata[3,*]-isodata[0,3]

Gmasses= [G3rdata[5,0], G2rdata[5,0], G1rdata[5,0],G0rdata[5,0]]
Gbes= [max(G3rbe),max(G2rbe),max(G1rbe),max(G0rbe)]


; let's do some fitting



; first we'll fit major merger be
print, "retrograde: first we'll fit the major merger burst efficiencies"
; ---------------------------------------------
weight= Gmasses*0.0 + 1.0
p1= [1.0,1.0]
;spec_sf_eff= MPFITFUN('minor_func_ssflaw', Gmasses, Gbes, weight, p1, PARINFO=pi, BESTNORM=bestnorm, DOF=dof)
spec_sf_eff= MPFITFUN('minor_func_ssflaw', Gmasses, Gbes, weight, p1, BESTNORM=bestnorm, DOF=dof)

redchi2= bestnorm/dof
print, "Reduced Chi^2 = ", redchi2
print, spec_sf_eff



; next we'll fit the burst efficiency
; ------------------------------------


; do them all at once
print, "retrograde: fit all 10 simulations burst efficiency vs. merer ratio"
Gmr= [transpose(G3mr[0,*]),transpose(G2mr[0,*]),transpose(G1mr[0,*]),transpose(G0mr[0,*])]
Gbe= [transpose(G3rbe[0,*]),transpose(G2rbe[0,*]),transpose(G1rbe[0,*]),transpose(G0rbe[0,*])]
weight= Gmr*0.0 + 1.0
p1= [1.0,1.0]
G_result= MPFITFUN('minor_func_samsbl', Gmr, Gbe, weight, p1, PARINFO=pi, BESTNORM=bestnorm, DOF=dof)
redchi2= bestnorm/dof
print, "Reduced Chi^2 = ", redchi2
print, G_result




xmax= 1.0 & xmin= 0.0 & ymax= 1.0 & ymin= 0.0
initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename='sfr_pere.eps', colortable=4

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, charthick=2.0, xthick=4.0, ythick=4.0, color= 0, $
        xtitle="prograde e (%)", $
        ytitle="retrograde e (%)", $
        /nodata

        oplot, G3data[3,*], G3rdata[3,*], psym=-2, thick=3.0, color= 150
        oplot, G2data[3,*], G2rdata[3,*], psym=-5, thick=3.0, color= 100
        oplot, G1data[3,*], G1rdata[3,*], psym=-7, thick=3.0, color= 50
        ;oplot, G0data[3,*], G0rdata[3,*], psym=-6, thick=3.0, color= 150

        x0= 0.25
        y0= 0.85
        xyouts, x0, y0, 'G3', /normal, charthick=3.0, size=1.3, color= 150
        xyouts, x0, y0-0.05, 'G2', /normal, charthick=3.0, size=1.3, color= 100
        xyouts, x0, y0-0.10, 'G1', /normal, charthick=3.0, size=1.3, color= 50

	x=[xmin,xmax] & y=x
	oplot, x, y, psym= -3, thick= 1.0, linestyle= 1, color= 0

device, /close



;============================================

xmax= 1.0 & xmin= 0.0 & ymax= 1.0 & ymin= 0.0
initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename='sfr_eble.eps', colortable=4

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, charthick=2.0, xthick=4.0, ythick=4.0, color= 0, $
        xtitle="G3 e (%)", $
        ytitle="bulge-less G3 e (%)", $
        /nodata

        oplot, G3data[3,*], G3bldata[3,*], psym=-2, thick=3.0, color= 150
        ;oplot, G2data[3,*], G2rdata[3,*], psym=-5, thick=3.0, color= 100
        ;oplot, G1data[3,*], G1rdata[3,*], psym=-7, thick=3.0, color= 50
        ;oplot, G0data[3,*], G0rdata[3,*], psym=-6, thick=3.0, color= 150

        x0= 0.25
        y0= 0.85
        ;xyouts, x0, y0, 'G3', /normal, charthick=3.0, size=1.3, color= 150
        ;xyouts, x0, y0-0.05, 'G2', /normal, charthick=3.0, size=1.3, color= 100
        ;xyouts, x0, y0-0.10, 'G1', /normal, charthick=3.0, size=1.3, color= 50

        x=[xmin,xmax] & y=x
        oplot, x, y, psym= -3, thick= 1.0, linestyle= 1, color= 0

device, /close



;--------------------------------------------







;============================================
;============================================
;
;    Max Star Formation Rate
;
;============================================
;============================================

;--------------------------------------
;  Print the Shit
;--------------------------------------

ymax = 15.0
ymin = 0.0
xmax = 1.05
xmin =-0.05

;---------------------------

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename='sfr_minsfr.eps', colortable=4

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, color=0, charthick=2.0, xthick=4.0, ythick=4.0, $
        xtitle="Mass Ratio (M!Dsat!N/M!Dprimary!N)", $
        ytitle="SFR!Imax!N (M!D!9n!3!N/yr)", $
        /nodata

        oplot, 1.0/G3data[0,*], G3data[4,*], psym=-2, thick=3.0, color= 150
        oplot, 1.0/G2data[0,*], G2data[4,*], psym=-5, thick=3.0, color= 100
        oplot, 1.0/G1data[0,*], G1data[4,*], psym=-7, thick=3.0, color= 50
        ;oplot, 1.0/G0[*,0], G0[*,4], psym=-6, thick=3.0, color= 150

        x0= 0.25
        y0= 0.85
        xyouts, x0, y0, 'G3', /normal, charthick=3.0, size=1.3, color= 150
        xyouts, x0, y0-0.05, 'G2', /normal, charthick=3.0, size=1.3, color= 100
        xyouts, x0, y0-0.10, 'G1', /normal, charthick=3.0, size=1.3, color= 50

device, /close

; ----------------------------------------

ymax= 0.3

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename='sfr_spsfr.eps', colortable=4

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, color=0, charthick=3.0, xthick=4.0, ythick=4.0, $
        xtitle="Mass Ratio (M!Dsat!N/M!Dprimary!N)", $
        ;ytitle="SFR (M!D!9n!3!N/yr) / Stellar Mass (10!E10!N M!D!9n!3!N)", $
	ytitle="specific SFR (Gyr!E-1!N)", $
        /nodata

	; factor of 10 is because data are in units of msun/yr divided by 10^10 msun,
	;     or 1/10Gyr, so divide data by 10 to get 1/Gyr
        oplot, 1.0/G3data[0,*], G3data[4,*]/G3data[5,*]/10.0, psym=-2, thick=3.0, color= 150, linestyle= 0
        oplot, 1.0/G2data[0,*], G2data[4,*]/G2data[5,*]/10.0, psym=-5, thick=3.0, color= 100, linestyle= 1
        oplot, 1.0/G1data[0,*], G1data[4,*]/G1data[5,*]/10.0, psym=-7, thick=3.0, color= 50, linestyle= 2
        oplot, 1.0/G0data[0,*], G0data[4,*]/G0data[5,*]/10.0, psym=-6, thick=3.0, color= 0

        x0= 0.30
        y0= 0.85
        xyouts, x0, y0, 'G3', /normal, charthick=3.0, size=1.3, color= 150
        xyouts, x0, y0-0.05, 'G2', /normal, charthick=3.0, size=1.3, color= 100
        xyouts, x0, y0-0.10, 'G1', /normal, charthick=3.0, size=1.3, color= 50
        xyouts, x0, y0-0.15, 'G0', /normal, charthick=3.0, size=1.3, color= 0
        oplot, [0.05], [0.266], psym=2, thick=3.0, symsize=1.5, color=150
        oplot, [0.05], [0.247], psym=5, thick=3.0, symsize=1.5, color=100
        oplot, [0.05], [0.229], psym=7, thick=3.0, symsize=1.5, color=50
        oplot, [0.05], [0.211], psym=6, thick=3.0, symsize=1.5, color= 0

        x0= 0.01 & xs= 0.21
        y0= 0.280 & ys= 0.080
        oplot, [x0,x0+xs,x0+xs,x0,x0],[y0,y0,y0-ys,y0-ys,y0], thick=3.0, color= 0, psym=-3


device, /close



; ---------------------------------------






end




