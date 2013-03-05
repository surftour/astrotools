;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   Krumholz & McKee SF versus Ours
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



; =========================
;
;
; =========================

;   SPH Particle Information



;  calculate KM05 SF rate
; ------------------------------------
pro KF_SFR, frun, snapnum, $
		extra=extra, $
		TurbulentFactor=TurbulentFactor, $
		grid_radius=grid_radius


if not keyword_set(frun) then begin
	print, "  "
	print, " KF_SFR, frun, snapnum"
	print, "  "
	return
endif

if not keyword_set(snapnum) then snapnum= 12

if not keyword_set(extra) then extra=''



;-------------------------------------------
;
;   Open Snapshot
;

;spawn, "/bin/ls "+frun+"/snap* | wc ",result                                    
;lastsnapshot=long(result[0])-1                                                  
;ok=fload_snapshot_bh(frun,lastsnapshot)         
ok=fload_snapshot_bh(frun,snapnum)         




; ----------------------------------------------------------------------
;
;   Center
;
; ----------------------------------------------------------------------

;print, " ---------------------"
;print, "  DETERMINE: center   "
;print, " ---------------------"

; determine_center
; -----------------
;center= [0,0,0]
;center_bh= fload_center_alreadycomp(1)

; two black holes
;if fload_npart(5) gt 1 then begin
;	bhid= fload_blackhole_id(1)
;	bhid1= bhid[0]
;	bhid2= bhid[1]
;	print, "Blackhole ID's: ", bhid1, bhid2
;	center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
;	center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
;	center_bh= center1
	;center_bh= center2
;endif

; one black hole
;if fload_npart(5) eq 1 then begin
;	bhid= fload_blackhole_id(1)
;	center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
;endif




; ----------------------------------------------------------------------
;
;  grab variables
;
; ----------------------------------------------------------------------

;print, "using center= ", center_bh
;x= fload_gas_xyz('x',center=center_bh)
;y= fload_gas_xyz('y',center=center_bh)
;z= fload_gas_xyz('z',center=center_bh)
m= fload_gas_mass(1)
;ids= fload_gas_id(1)
hsml= fload_gas_hsml(1)


sfr= fload_gas_sfr(1)


; ----------------------------------------------------------------------
;
;  determine cold gas
;
; ----------------------------------------------------------------------


print, " "
print, " "
print, " ================================ "
print, " "
print, " "
print, " N_Gas= ", n_elements(m)
print, " Total Gas Mass= ", total(m)
print, " "
print, " Total SFR= ", total(sfr)
print, " "

; -------------------
; density information
rho= fload_gas_rho(1)
print, " rho: min/max= ", min(rho), max(rho)

u= fload_gas_u(1)
print, " u: min/max= ", min(u), max(u)


Pressure= rho * u

; -------------------
; -------------------

  ; define the cut-off density
  ; for multi-phase model

  rho_crit= 0.000854924     ; std GADGET value

; -------------------
; -------------------


idx= where(rho gt rho_crit)
print, " "
print, "   Star-forming Gas "
print, "      N_sf= ", n_elements(idx)
print, "      Mass= ", total(m(idx))
print, " "


; --------------------
;  This is
;
;     TurbulentFactor = 1 + (P_turb / P_thermal)
; 
if not keyword_set(TurbulentFactor) then begin
	TurbulentFactor= 1.0
	;TurbulentFactor= 6.0
	;TurbulentFactor= 11.0
	;TurbulentFactor= 101.0
endif
; --------------------
h= 0.7
; --------------------


        ;  Contants
        ; ---------------------------
        HYDROGEN_MASSFRAC= 0.76
        GAMMA_MINUS1= (5./3.) - 1.0
        BOLTZMANN=   1.3806d-16
        PROTONMASS=  1.6726d-24
        UnitMass_in_g= 1.989d+43
        UnitDensity_in_cgs = 6.76991d-22
        UnitEnergy_in_cgs = 1.989d+53

        ;  Cold Phase Temperature
        ; ---------------------------
        TempClouds    = 1000.0           ; in K
        meanweight = 4. / (1. + 3. * HYDROGEN_MASSFRAC)
        EgySpecCold = 1. / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * TempClouds
        u_cold = EgySpecCold * UnitMass_in_g / UnitEnergy_in_cgs   ;  now in gadget units 


rho_cold= fload_gas_rho(1,/cold,TurbulentFactor=TurbulentFactor)
idx= where(rho_cold gt 0.0)
print, "   Cold Gas (should be same as Star-forming) "
print, "      N_cold (rho_cold > 0)= ", n_elements(idx)
;print, "      Mass= ", total(m(idx))
print, "      rho_cold max/min: ", max(rho_cold), min(rho_cold(idx)), min(rho_cold)

coldf= fload_gas_coldfraction(1)
print, "      cold gas fraction max/min: ", max(coldf), min(coldf)

ColdGasMass= m*coldf
print, "      Mass= ", total(ColdGasMass)

; volume_cold / volume_hot
coldgasff= fload_gas_coldgasfillingfactor(1,TurbulentFactor=TurbulentFactor)
print, "      coldgasff max/min: ", max(coldgasff), min(coldgasff)

; volume_cold / volume_total
coldgasvf= fload_gas_coldgasvolumefraction(1,TurbulentFactor=TurbulentFactor)
print, "      coldgasvf max/min: ", max(coldgasvf), min(coldgasvf)
print, " "


; BAD clumpsize calculation.  The hsml is not a good 
; indicator of the size of the SPH particle
;
;clumpsize= hsml(idx) * (coldgasvf(idx)^(1./3.))
;print, "      clumpsize max/min: ", max(clumpsize), min(clumpsize)
;print, "                   mean: ", mean(clumpsize)
;

; cold clump size
newclumpsize= ( (4.*!PI/ 3.) * rho_cold(idx) / ColdGasMass(idx))^(-1./3.)
clumpsize= newclumpsize
print, "(new) clumpsize max/min: ", max(clumpsize), min(clumpsize)
print, "                   mean: ", mean(clumpsize)

print, " "
print, " "




; -----------------------------------------------
;
;    Just the cold/SF gas information

	N_ColdGas= n_elements(idx)
	ColdGas_Mass= ColdGasMass(idx) / h
	ColdGas_Density= rho_cold(idx) * h * h
	ColdGas_Size= clumpsize / h
	ColdGas_Pressure= Pressure(idx) * h * h
	SphGas_Hsml= hsml(idx) / h

	Snap_Sfr= sfr(idx)

; -----------------------------------------------

print, "xxxxxxxxxxxxxxxxxx"
print, "  Check SPH info  "
print, "  "
cgmass= 4.*!PI/3. * ColdGas_Size * ColdGas_Size * ColdGas_Size * ColdGas_Density
print, "total mass= ", total(cgmass)
print, "  "
print, "xxxxxxxxxxxxxxxxxx"

; -----------------------------------------------
;


;  Now, calculate the Krumholz and McKee SFR

coldg_ff_time = sqrt(3. * !PI / 32. / 43007.1 / ColdGas_Density)
;coldg_ff_time = 1 / sqrt(4. * !PI * 43007.1 * ColdGas_Density)
;coldg_ff_time= 0.0047 * (ColdGas_Mass / 1.0e+6)^(0.25)
print, "ff time min/max= ", min(coldg_ff_time),max(coldg_ff_time)

radius= (3. * 43007.1 * ColdGas_Mass * ColdGas_Mass / 40. / !PI / ColdGas_Pressure)^(0.25)
print, "radius min/max= ", min(radius), max(radius)

sigma= sqrt(43007.1 * ColdGas_Mass / 5. / radius)
print, "sigma min/max= ", min(sigma), max(sigma)

;sigma= sqrt(43007.1 * ColdGas_Mass / 5.0 / ColdGas_Radius)
;c_s= 5.0   ; in km/sec  (a free variable)
;c_s= 0.5   ; in km/sec  (a free variable)
c_s= 0.3   ; in km/sec  (a free variable)
Mach= sigma / c_s

print, "Mach min/max= ", min(Mach), max(Mach)

KM_sfr= 0.014 * (0.76)^(-0.68) * (Mach/100)^(-0.32) * ColdGas_Mass / coldg_ff_time


print, "  "
print, "  "
print, "  KM Sfr = ", total(KM_sfr)
print, "       (reg SFR factor ", total(KM_sfr)/total(sfr), " )"
print, "  "
print, "  snap Sfr= ", total(Snap_sfr)
print, "  frun: ", frun
print, "  snapnum: ", snapnum
print, "  simulation time: ", fload_timelbl(1,4,/noteq)
print, "  "


; -----------------------------------------------

; generate histogram of mach numbers
mach_hist, Mach


; compare the two SFR's
sfr_v_sfr, Snap_sfr, KM_sfr



; stop

end





; =============================================================================
; =============================================================================
; =============================================================================
; =============================================================================



;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Histogram of Mach Numbers
;     ------------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro mach_hist, machs

if not keyword_set(machs) then begin
   print, "  "
   print, "mach_hist, machs"
   print, "  "
   return
endif

filename='machhist.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; radius, linear scale
xaxistitle= "!8Mach Number!6 (!7r!6/c!Ds!N)"
xmax = 100.0
xmin = 0.0

; number (histogram)
;yaxistitle= "!6V!Dr!N (km s!E-1!N)"
yaxistitle= ' '
ymax = 1.2
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.05
y0= 0.15
x1= 0.95
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
	;/ylog, $
	;/xlog, $
	color= 0, $
	xrange=[xmin,xmax], $
	yrange=[ymin,ymax], $
	xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
	xtitle=xaxistitle, $
	ytickformat='(a1)', $
	;ytitle=yaxistitle, $
	/nodata, /noerase


; -----------------------------------------------


; ------------------
; compute histogram
; ------------------



temp= process_histogram(machs, xmax=xmax, xmin=xmin, levels=100, oplotit=50)
print, min(temp), max(temp)

; -----------------------------------------------

avg_machs= mean(machs)
print, "avg= ", avg_machs
avglbl= strcompress(string(avg_machs),/remove_all)
xyouts, 0.10, 0.88, "<!8M!6> = "+avglbl, size=1.2, color=0, /normal, charthick=3.0



;xyouts, 0.75, 0.90, fload_timelbl(1,2), size=1.2, color=0, /normal, charthick=4.0
;xyouts, 0.55, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0
;xyouts, 0.75, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0

;xyouts, 0.10, 0.85, "!6V!D200!N = 50 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0
;xyouts, 0.10, 0.85, "!6V!D200!N = 500 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0
;xyouts, 0.10, 0.85, "!6V!D200!N = 160 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0



; print extras
; -------------


; done
; -----
device, /close


end












;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Sfr (KM) vs. Sfr (snap)
;     ------------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro sfr_v_sfr, sfr, KM_sfr

if not keyword_set(sfr) then begin
   print, "  "
   print, "sfr_v_sfr, sfr, KM_sfr"
   print, "  "
   return
endif

filename='sfrvsfr.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


xaxistitle= "!6Log SFR from Snapshot (M!D!9n!6!N Yr!E-1!N)"
xmax = 1.0
xmin = -5.0

yaxistitle= "!6Log SFR from KM05 (M!D!9n!6!N Yr!E-1!N)"
ymax = 1.0
ymin = -5.0



; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
	;/ylog, $
	;/xlog, $
	color= 0, $
	xrange=[xmin,xmax], $
	yrange=[ymin,ymax], $
	xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
	xtitle=xaxistitle, $
	;ytickformat='(a1)', $
	ytitle=yaxistitle, $
	/nodata, /noerase


; -----------------------------------------------

;oplot, sfr, KM_sfr, psym=2, color=150

KM_Factor= total(KM_sfr)/total(sfr)

sfr_log= alog10(sfr)
KM_sfr_log= alog10(KM_sfr)
oplot, sfr_log, KM_sfr_log, psym=2, color=150



; print extras
; -------------

x= [xmin,xmax]
y= x
oplot, x, y, psym=-3, color=0, linestyle= 2, thick= 2.0


; SFR factor
; -----------
KMlbl= strcompress(string(KM_Factor),/remove_all)
xyouts, 0.52, 0.25, "KM/Snap = "+KMlbl, size=1.2, color=0, /normal, charthick=3.0



; done
; -----
device, /close


end











