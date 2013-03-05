;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;  Figure motivated by disucssions with Todd Thomson.
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro toddt, frun, snapnum, filename=filename
        

if not keyword_set(frun) then begin
   print, "  "
   print, " toddt, frun, snapnum, filename=filename, "
   print, "  "
   print, "  "
   print, "  "
   return
endif   
        
if not keyword_set(snapnum) then snapnum=12
if not keyword_set(filename) then filename="toddt.eps"


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4




; ----------------------------- 
; set up constants      
; -----------------------------
                        
bins = 50
                

xaxistitle= '!6R (kpc)'
;xmax = 50.0
xmax = 10.0
xmin = 0.1
;xmin = 0.04
;xmin= 0.0


yaxistitle= '!6L!Dtot!N/L!DEdd!N'
ymax = 50
ymin = 0.00005


;---------------------------
;  Print it
;---------------------------
   
!p.position= [0.18, 0.15, 0.98, 0.98]
   
plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
        /xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
	ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; ------------------------
;  determine gas profile
; ------------------------

;plot_one_ltotledd, frun, snapnum, xmin, xmax, bins
;plot_one_ltotledd2, frun, snapnum, xmin, xmax, bins


;xyouts, 0.25, 0.38, "!7j!6= 10.0 cm!E2!N/g", /normal, color= 0, charthick=2.0, size=1.2

xyouts, 0.25, 0.40, "!7j!6 = !7R!6!E-1!N ", /normal, color= 0, charthick=2.0, size=1.5


;plot_one_ltotledd, "/raid4/tcox/ds/vc3vc3e_2", 34, xmin, xmax, bins, lcolor= 50
plot_one_ltotledd_new, "/raid4/tcox/ds/vc3vc3e_2", 34, xmin, xmax, bins, lcolor= 50
xyouts, 0.25, 0.33, "pre-sb", /normal, color= 50, charthick=2.0, size=1.2

;plot_one_ltotledd, "/raid4/tcox/ds/vc3vc3e_2", 40, xmin, xmax, bins, lcolor= 150
plot_one_ltotledd_new, "/raid4/tcox/ds/vc3vc3e_2", 40, xmin, xmax, bins, lcolor= 150
xyouts, 0.25, 0.28, "peak-sb", /normal, color= 150, charthick=2.0, size=1.2
plot_one_ltotledd_new, "/raid4/tcox/ds/vc3vc3e_2", 40, xmin, xmax, bins, lcolor= 0

;plot_one_ltotledd, "/raid4/tcox/ds/vc3vc3e_2", 52, xmin, xmax, bins, lcolor= 100
plot_one_ltotledd_new, "/raid4/tcox/ds/vc3vc3e_2", 52, xmin, xmax, bins, lcolor= 100
xyouts, 0.25, 0.23, "post-sb", /normal, color= 100, charthick=2.0, size=1.2



; -----------------
;  Plot Extras
; -----------------

; draw 1
x=[xmin,xmax]
y=[1.0,1.0]
oplot, x, y, linestyle=1, color= 0



;--------------------------------------
;--------------------------------------

device, /close



end





;==================================================================================






; actually do the work
; --------------------------
pro plot_one_ltotledd, frun, snapnum, xmin, xmax, bins, lcolor=lcolor


	ok=fload_snapshot_bh(frun,snapnum)

	center=fload_center_alreadycomp(1)
	print, "using center= ", center


	; stellar info
	m_stars= fload_allstars_mass(1)
	print, "total stellar mass= ", total(m_stars)
	m_stars= m_stars / 0.7
	r_stars=fload_allstars_xyz('r',center=center) / 0.7


	; luminosities
	bololum=fload_allstars_bololum(1)
	TTime= float(fload_time(1))
	bhlum= fload_blackhole_lum(frun,Ttime,/bolometric)


	; gaseous info
	m_gas= fload_gas_mass(1)
	print, "total gas mass= ", total(m_gas)
	m_gas= m_gas / 0.7
	r_gas= fload_gas_xyz('r', center=center) / 0.7

	; dark mass info
	m_dark= fload_halo_mass(1)
	print, "total dark mass= ", total(m_dark)
	m_dark= m_dark / 0.7
	r_dark= fload_halo_xyz('r', center=center) / 0.7



	; now, calculate profiles
	; --------------------------
	this_xmax= alog10(xmax)
	this_xmin= alog10(xmin)

	; cumulative stellar mass
	r_stars= alog10(r_stars)
	;process_prof_tot, r_stars, m_stars, bins, this_xmax, this_xmin, rs_stars, ms_stars
	process_cummass_profile, r_stars, m_stars, this_xmin, this_xmax, bins, rs_stars, ms_stars
	print, "stars max/min= ", max(ms_stars), min(ms_stars)

	; cumulative gas mass
	r_gas= alog10(r_gas)
	;process_prof_tot, r_gas, m_gas, bins, this_xmax, this_xmin, rs_gas, ms_gas
	process_cummass_profile, r_gas, m_gas, this_xmin, this_xmax, bins, rs_gas, ms_gas
	print, "gas max/min= ", max(ms_gas), min(ms_gas)

	; cumulative dark mass
	r_dark= alog10(r_dark)
	;process_prof_tot, r_dark, m_dark, bins, this_xmax, this_xmin, rs_dark, ms_dark
	process_cummass_profile, r_dark, m_dark, this_xmin, this_xmax, bins, rs_dark, ms_dark
	print, "dark mass max/min= ", max(ms_dark), min(ms_dark)

	; cumulative stellar luminosity
	;process_prof_tot, r_stars, bololum, bins, this_xmax, this_xmin, rs_stars, lum_stars
	process_cummass_profile, r_stars, bololum, this_xmin, this_xmax, bins, rs_stars, lum_stars
	print, "stellar luminosity (solar masses) max/min= ", max(lum_stars), min(lum_stars)



	G= 43007.1  ; (gadget units)
	G_cgs= 6.67259d-8 	; cm3 g-1 s-2
	cc= 3.0e+5   ; km/s  (gadget units)
	cc_cgs= 3.0e+10 ; cm/s
	UnitLength_in_cm= 3.08568d+21  ; cm 
	UnitSolarLuminosity_in_ergs= 3.827d+33   ; erg/s.
	UnitMass_in_g = 1.989d+43 
	UnitTime_in_s = 3.08568d+16 
	UnitVelocity_in_cm_per_s = 100000 
	UnitDensity_in_cgs = 6.7699d-22 
	UnitEnergy_in_cgs = 1.989d+53 


	; eddington luminosity (gadget units)
	mass_total = ms_stars + ms_gas + ms_dark
	kappa= 1.0    ; cm^2/g   
	kappa= kappa * UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm
	l_edd= 4 * !PI * G * mass_total * cc / kappa

	; eddington luminosity (cgs units)
	kappa= 1.0    ; cm^2/g   
	l_edd_cgs= 4 * !PI * G_cgs * mass_total * UnitMass_in_g *  cc_cgs / kappa

	; total luminosity
	lum_total= lum_stars
	for i= 0, bins-1 do lum_total[i]= lum_total[i] + total(bhlum)

	lum_total_cgs= lum_total * UnitSolarLuminosity_in_ergs

	solarlum_in_gadget_units= UnitSolarLuminosity_in_ergs * UnitTime_in_s / UnitEnergy_in_cgs
	lum_total= lum_total * solarlum_in_gadget_units ; gadget units



	;
	; l_tot / l_edd
	;
	ltotledd= lum_total/l_edd
	ltotledd_cgs= lum_total_cgs/l_edd_cgs

	rs= 10^(rs_stars)
	;oplot, rs, ltotledd, psym=-3, linestyle=0, color= 150, thick=4.0
	;xyouts, 0.65, 0.88, "!7j!6= 1.0 cm!E2!N/g", /normal, color= 150, charthick=2.0, size=1.2


	; now use kappa= 10 cm2/g
	ltotledd_10= ltotledd * 10.0
	oplot, rs, ltotledd_10, psym=-3, linestyle=0, color= lcolor, thick=4.0
	;xyouts, 0.65, 0.83, "!7j!6= 10.0 cm!E2!N/g", /normal, color= 50, charthick=2.0, size=1.2


	; now use kappa= 100 cm2/g
	;ltotledd_100= ltotledd * 100.0
	;oplot, rs, ltotledd_100, psym=-3, linestyle=0, color= 100, thick=4.0
	;xyouts, 0.65, 0.78, "!7j!6= 100.0 cm!E2!N/g", /normal, color= 100, charthick=2.0, size=1.2


	;
	; mass of gas, divided by mass of stars
	;
	;gas_to_stars= ms_gas/ms_stars
	;rs_gas= 10^(rs_gas)
	;oplot, rs_gas, gas_to_stars, psym=-3, linestyle=0, color= 150, thick=4.0
	;xyouts, 0.65, 0.88, 'stars', /normal, color= 150, charthick=2.0, size=1.2


	


end





;==================================================================================




; actually do the work
; --------------------------
pro plot_one_ltotledd2, frun, snapnum, xmin, xmax, bins


	ok=fload_snapshot_bh(frun,snapnum)

	center=fload_center_alreadycomp(1)
	print, "using center= ", center


	; stellar info
	m_stars= fload_allstars_mass(1)
	print, "total stellar mass= ", total(m_stars)
	m_stars= m_stars / 0.7
	r_stars=fload_allstars_xyz('r',center=center) / 0.7


	; luminosities
	bololum=fload_allstars_bololum(1)
	TTime= float(fload_time(1))
	bhlum= fload_blackhole_lum(frun,Ttime,/bolometric)


	; gaseous info
	m_gas= fload_gas_mass(1)
	print, "total gas mass= ", total(m_gas)
	m_gas= m_gas / 0.7
	r_gas= fload_gas_xyz('r', center=center) / 0.7

	; dark mass info
	m_dark= fload_halo_mass(1)
	print, "total dark mass= ", total(m_dark)
	m_dark= m_dark / 0.7
	r_dark= fload_halo_xyz('r', center=center) / 0.7



	; now, calculate profiles
	; --------------------------
	this_xmax= alog10(xmax)
	this_xmin= alog10(xmin)

	; cumulative stellar mass
	r_stars= alog10(r_stars)
	;process_prof_tot, r_stars, m_stars, bins, this_xmax, this_xmin, rs_stars, ms_stars
	process_cummass_profile, r_stars, m_stars, this_xmin, this_xmax, bins, rs_stars, ms_stars
	print, "stars max/min= ", max(ms_stars), min(ms_stars)

	; cumulative gas mass
	r_gas= alog10(r_gas)
	;process_prof_tot, r_gas, m_gas, bins, this_xmax, this_xmin, rs_gas, ms_gas
	process_cummass_profile, r_gas, m_gas, this_xmin, this_xmax, bins, rs_gas, ms_gas
	print, "gas max/min= ", max(ms_gas), min(ms_gas)

	; cumulative dark mass
	r_dark= alog10(r_dark)
	;process_prof_tot, r_dark, m_dark, bins, this_xmax, this_xmin, rs_dark, ms_dark
	process_cummass_profile, r_dark, m_dark, this_xmin, this_xmax, bins, rs_dark, ms_dark
	print, "dark mass max/min= ", max(ms_dark), min(ms_dark)

	; cumulative stellar luminosity
	;process_prof_tot, r_stars, bololum, bins, this_xmax, this_xmin, rs_stars, lum_stars
	process_cummass_profile, r_stars, bololum, this_xmin, this_xmax, bins, rs_stars, lum_stars
	print, "stellar luminosity (solar masses) max/min= ", max(lum_stars), min(lum_stars)



	G= 43007.1  ; (gadget units)
	G_cgs= 6.67259d-8 	; cm3 g-1 s-2
	cc= 3.0e+5   ; km/s  (gadget units)
	cc_cgs= 3.0e+10 ; cm/s
	UnitLength_in_cm= 3.08568d+21  ; cm 
	UnitSolarLuminosity_in_ergs= 3.827d+33   ; erg/s.
	UnitMass_in_g = 1.989d+43 
	UnitTime_in_s = 3.08568d+16 
	UnitVelocity_in_cm_per_s = 100000 
	UnitDensity_in_cgs = 6.7699d-22 
	UnitEnergy_in_cgs = 1.989d+53 


	; eddington luminosity (gadget units)
	mass_total = ms_stars + ms_gas + ms_dark
	kappa_cgs= 1.0    ; cm^2/g   
	kappa= kappa_cgs * UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm
	l_edd= 4 * !PI * G * mass_total * cc / kappa

	; eddington luminosity (cgs units)
	l_edd_cgs= 4 * !PI * G_cgs * mass_total * UnitMass_in_g *  cc_cgs / kappa_cgs

	; total luminosity
	lum_total= lum_stars
	for i= 0, bins-1 do lum_total[i]= lum_total[i] + bhlum

	lum_total_cgs= lum_total * UnitSolarLuminosity_in_ergs

	solarlum_in_gadget_units= UnitSolarLuminosity_in_ergs * UnitTime_in_s / UnitEnergy_in_cgs
	lum_total= lum_total * solarlum_in_gadget_units ; gadget units



	;
	; l_tot / l_edd
	;
	ltotledd= lum_total/l_edd
	ltotledd_cgs= lum_total_cgs/l_edd_cgs

	; now use kappa= 10 cm2/g
	rs= 10^(rs_stars)
	ltotledd_10= ltotledd * 10.0
	oplot, rs, ltotledd_10, psym=-3, linestyle=0, color= 50, thick=4.0
	xyouts, 0.55, 0.38, "!7j!6= 10.0 cm!E2!N/g", /normal, color= 0, charthick=2.0, size=1.2
	xyouts, 0.55, 0.33, "all masses", /normal, color= 50, charthick=2.0, size=1.2


	; now use kappa= 10 cm2/g, no dark matter
	mass_total = ms_gas
	l_edd= 4 * !PI * G * mass_total * cc / kappa
	ltotledd_nodm= lum_total/l_edd
	ltotledd_nodm= ltotledd_nodm * 10.0
	oplot, rs, ltotledd_nodm, psym=-3, linestyle=0, color= 100, thick=4.0
	xyouts, 0.55, 0.28, "gas mass only", /normal, color= 100, charthick=2.0, size=1.2


	;
	; mass of gas, divided by mass of stars
	;
	;gas_to_stars= ms_gas/ms_stars
	;rs_gas= 10^(rs_gas)
	;oplot, rs_gas, gas_to_stars, psym=-3, linestyle=0, color= 150, thick=4.0
	;xyouts, 0.65, 0.88, 'stars', /normal, color= 150, charthick=2.0, size=1.2


	


end




;==================================================================================



; actually do the work
; --------------------------
pro plot_one_ltotledd_new, frun, snapnum, xmin, xmax, bins, lcolor= lcolor


	ok=fload_snapshot_bh(frun,snapnum)

	center=fload_center_alreadycomp(1)
	print, "using center= ", center


	; stellar info
	m_stars= fload_allstars_mass(1)
	print, "total stellar mass= ", total(m_stars)
	m_stars= m_stars / 0.7
	r_stars=fload_allstars_xyz('r',center=center) / 0.7


	; luminosities
	bololum=fload_allstars_bololum(1)
	TTime= float(fload_time(1))
	bhlum= fload_blackhole_lum(frun,Ttime,/bolometric)


	; gaseous info
	m_gas= fload_gas_mass(1)
	print, "total gas mass= ", total(m_gas)
	m_gas= m_gas / 0.7
	r_gas= fload_gas_xyz('r', center=center) / 0.7


	; diffuse gas info
	cf= fload_gas_coldfraction(1)
	m_gas_diffuse= m_gas * (1.0 - cf)
	r_gas_diffuse= r_gas


	; dark mass info
	m_dark= fload_halo_mass(1)
	print, "total dark mass= ", total(m_dark)
	m_dark= m_dark / 0.7
	r_dark= fload_halo_xyz('r', center=center) / 0.7



	; now, calculate profiles
	; --------------------------
	this_xmax= alog10(xmax)
	this_xmin= alog10(xmin)

	; cumulative stellar mass
	r_stars= alog10(r_stars)
	;process_prof_tot, r_stars, m_stars, bins, this_xmax, this_xmin, rs_stars, ms_stars
	process_cummass_profile, r_stars, m_stars, this_xmin, this_xmax, bins, rs_stars, ms_stars
	print, "stars max/min= ", max(ms_stars), min(ms_stars)

	; cumulative gas mass
	r_gas= alog10(r_gas)
	;process_prof_tot, r_gas, m_gas, bins, this_xmax, this_xmin, rs_gas, ms_gas
	process_cummass_profile, r_gas, m_gas, this_xmin, this_xmax, bins, rs_gas, ms_gas
	print, "gas max/min= ", max(ms_gas), min(ms_gas)


	; gas density
	radius= r_gas
	mass= m_gas
	rho_gas= fltarr(bins)
	radius_gas= fltarr(bins)
	weight= fltarr(bins)
	process_prof_rho, radius, mass, bins, this_xmax, this_xmin,  radius_gas, rho_gas, weight


	; cumulative diffuse gas mass
	r_gas_diffuse= alog10(r_gas_diffuse)
	process_cummass_profile, r_gas_diffuse, m_gas_diffuse, this_xmin, this_xmax, bins, rs_gas_diffuse, ms_gas_diffuse
	print, "diffuse gas max/min= ", max(ms_gas_diffuse), min(ms_gas_diffuse)


	; diffuse gas density
	radius_diffuse= r_gas_diffuse
	mass_diffuse= m_gas_diffuse
	rho_gas_diffuse= fltarr(bins)
	radius_gas_diffuse= fltarr(bins)
	weight_diffuse= fltarr(bins)
	process_prof_rho, radius_diffuse, mass_diffuse, bins, this_xmax, this_xmin,  radius_gas_diffuse, rho_gas_diffuse, weight_diffuse


	; cumulative dark mass
	r_dark= alog10(r_dark)
	;process_prof_tot, r_dark, m_dark, bins, this_xmax, this_xmin, rs_dark, ms_dark
	process_cummass_profile, r_dark, m_dark, this_xmin, this_xmax, bins, rs_dark, ms_dark
	print, "dark mass max/min= ", max(ms_dark), min(ms_dark)

	; cumulative stellar luminosity
	;process_prof_tot, r_stars, bololum, bins, this_xmax, this_xmin, rs_stars, lum_stars
	process_cummass_profile, r_stars, bololum, this_xmin, this_xmax, bins, rs_stars, lum_stars
	print, "stellar luminosity (solar masses) max/min= ", max(lum_stars), min(lum_stars)



	G= 43007.1  ; (gadget units)
	G_cgs= 6.67259d-8 	; cm3 g-1 s-2
	cc= 3.0e+5   ; km/s  (gadget units)
	cc_cgs= 3.0e+10 ; cm/s
	UnitLength_in_cm= 3.08568d+21  ; cm 
	UnitSolarLuminosity_in_ergs= 3.827d+33   ; erg/s.
	UnitMass_in_g = 1.989d+43 
	UnitTime_in_s = 3.08568d+16 
	UnitVelocity_in_cm_per_s = 100000 
	UnitDensity_in_cgs = 6.7699d-22 
	UnitEnergy_in_cgs = 1.989d+53 
	solarlum_in_gadget_units= UnitSolarLuminosity_in_ergs * UnitTime_in_s / UnitEnergy_in_cgs


	; eddington luminosity 
	; ---------------------

	; total mass
	mass_total = ms_stars + ms_gas + ms_dark

	; opacity
	;kappa_cgs= 1.0    ; cm^2/g   
	;kappa= kappa_cgs * UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm
	sigma_gas= rs_stars*0.0
	sigma_gas_diffuse= rs_stars*0.0
	r0= 10^(rs_stars[0])
	sigma_gas[0]= rho_gas[0]*r0
	sigma_gas_diffuse[0]= rho_gas_diffuse[0]*r0
	for i= 1, bins-1 do begin
		;rmid= 10^(0.5*(rs_stars[i]-rs_stars[i-1]))
		rmin= 10^(rs_stars[i-1])
		rmax= 10^(rs_stars[i])
		dr= rmax-rmin
		sigma_gas[i]= sigma_gas[i-1] + rho_gas[i]*dr
		sigma_gas_diffuse[i]= sigma_gas_diffuse[i-1] + rho_gas_diffuse[i]*dr
	endfor

	kappa= 1/sigma_gas
	kappa_diffuse= 1/sigma_gas_diffuse

	sigma_gas_cgs = sigma_gas * UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm
	kappa_cgs= 1/sigma_gas_cgs
	sigma_gas_diffuse_cgs = sigma_gas_diffuse * UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm
	kappa_diffuse_cgs= 1/sigma_gas_diffuse_cgs


	; eddington luminosity
	l_edd= 4 * !PI * G * mass_total * cc / kappa
	l_edd_cgs= 4 * !PI * G_cgs * mass_total * UnitMass_in_g *  cc_cgs / kappa_cgs

	l_edd_diffuse= 4 * !PI * G * mass_total * cc / kappa_diffuse
	l_edd_diffuse_cgs= 4 * !PI * G_cgs * mass_total * UnitMass_in_g *  cc_cgs / kappa_diffuse_cgs


	; total luminosity
	lum_total= lum_stars
	for i= 0, bins-1 do begin
	   lum_total[i]= total(lum_stars[i] + bhlum)
	endfor


	lum_total= lum_total * solarlum_in_gadget_units ; gadget units
	lum_total_cgs= lum_total * UnitSolarLuminosity_in_ergs



	;
	; l_tot / l_edd
	;
	ltotledd= lum_total/l_edd
	print, "*** ltotledd ***"
	print, ltotledd
	ltotledd_cgs= lum_total_cgs/l_edd_cgs
	print, "*** ltotledd_cgs ***"
	print, ltotledd_cgs

	ltotledd_diffuse= lum_total/l_edd_diffuse


; the above (cgs vs. non) are off by exactly 1684.4713 ??

	; now plot
	rs= 10^(rs_stars)


	if lcolor gt 0 then begin
		oplot, rs, ltotledd, psym=-3, linestyle=0, color= lcolor, thick=4.0
	endif else begin
		oplot, rs, ltotledd_diffuse, psym=-3, linestyle=1, color= lcolor, thick=2.0
;stop
	endelse

	;xyouts, 0.55, 0.38, "!7j!6 = !7R!6!E-1!N ", /normal, color= 0, charthick=2.0, size=1.2
	;xyouts, 0.55, 0.33, "std", /normal, color= 50, charthick=2.0, size=1.2



end




;==================================================================================



; actually do the work
; --------------------------
pro plot_one_sigma, frun, snapnum, xmin, xmax, bins, $
			lcolor=lcolor, $
			loadedsnap=loadedsnap


	if not keyword_set(loadedsnap) then ok=fload_snapshot_bh(frun,snapnum)

	center=fload_center_alreadycomp(1)
	print, "using center= ", center


	; gaseous info
	m_gas= fload_gas_mass(1)
	print, "total gas mass= ", total(m_gas)
	m_gas= m_gas / 0.7
	r_gas= fload_gas_xyz('r', center=center) / 0.7



	; now, calculate profiles
	; --------------------------
	this_xmax= alog10(xmax)
	this_xmin= alog10(xmin)


	; gas density
	radius= alog10(r_gas)
	mass= m_gas
	rho_gas= fltarr(bins)
	radius_gas= fltarr(bins)
	weight= fltarr(bins)
	process_prof_rho, radius, mass, bins, this_xmax, this_xmin,  radius_gas, rho_gas, weight



	G= 43007.1  ; (gadget units)
	G_cgs= 6.67259d-8 	; cm3 g-1 s-2
	cc= 3.0e+5   ; km/s  (gadget units)
	cc_cgs= 3.0e+10 ; cm/s
	PROTONMASS=  1.6726e-24   ; g
	UnitLength_in_cm= 3.08568d+21  ; cm 
	UnitSolarLuminosity_in_ergs= 3.827d+33   ; erg/s.
	UnitMass_in_g = 1.989d+43 
	UnitTime_in_s = 3.08568d+16 
	UnitVelocity_in_cm_per_s = 100000 
	UnitDensity_in_cgs = 6.7699d-22 
	UnitEnergy_in_cgs = 1.989d+53 
	solarlum_in_gadget_units= UnitSolarLuminosity_in_ergs * UnitTime_in_s / UnitEnergy_in_cgs


	; column density
	; ---------------------

	sigma_gas= fltarr(bins)

	r0= 10^(radius_gas[0])
	sigma_gas[0]= rho_gas[0]*r0
	for i= 1, bins-1 do begin
		;rmid= 10^(0.5*(rs_stars[i]-rs_stars[i-1]))
		rmin= 10^(radius_gas[i-1])
		rmax= 10^(radius_gas[i])
		dr= rmax-rmin
		sigma_gas[i]= sigma_gas[i-1] + rho_gas[i]*dr
	endfor
	sigma_gas = sigma_gas * UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm

	; actually in cm-2 (not g)
	sigma_gas = sigma_gas / MASSPROTON

	sigma_gas= alog10(sigma_gas)
	rs= 10^(radius_gas)
	oplot, rs, sigma_gas, psym=-3, linestyle=0, color= lcolor, thick=4.0


end





;==================================================================================


; --------------------------
; process the cumulative profile

pro process_cummass_profile, x, y, xmin, xmax, bins, $
                                xs, ys, $
                                xlog=xlog

xs = fltarr(bins)
ys= fltarr(bins)

tempxmax= xmax
tempxmin= xmin
if keyword_set(xlog) then begin
        ;tempxmin= alog10(xmin)
        ;tempxmax= alog10(xmax)
        x= alog10(x)
endif

binsize = float((tempxmax-tempxmin))/bins

for i=0, bins-1 do begin 
    currentr= tempxmin + (i+1)*binsize
    idx= where(x le currentr)
    xs[i]= currentr 
    if idx(0) ne -1 then ys[i]= total(y(idx)) else ys[i]= 0.0
endfor



end



;==================================================================================




pro testsig, frun, snapnum, filename=filename
        

if not keyword_set(frun) then begin
   print, "  "
   print, " testsig, frun, snapnum, filename=filename, "
   print, "  "
   print, "  "
   print, "  "
   return
endif   
        
if not keyword_set(snapnum) then snapnum=12
if not keyword_set(filename) then filename="testsig.eps"


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4




; ----------------------------- 
; set up constants      
; -----------------------------
                        
bins = 50
                

xaxistitle= '!6R (kpc)'
;xmax = 50.0
xmax = 10.0
xmin = 0.1
;xmin = 0.04
;xmin= 0.0


yaxistitle= 'Log !7R!6!Dgas!N'
ymax = 24.0
ymin = 18.0


;---------------------------
;  Print it
;---------------------------
   
!p.position= [0.18, 0.15, 0.98, 0.98]
   
plot, [0], [0], psym=-3, linestyle=0, $
        ;/ylog, $
        /xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
	;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; ------------------------
;  determine gas profile
; ------------------------


plot_one_sigma, "/raid4/tcox/ds/vc3vc3e_2", 52, xmin, xmax, bins, lcolor= 50
xyouts, 0.25, 0.23, "post-starburst", /normal, color= 50, charthick=2.0, size=1.2



; -----------------
;  Plot Extras
; -----------------

; draw 1
;x=[xmin,xmax]
;y=[1.0,1.0]
;oplot, x, y, linestyle=1, color= 0



;--------------------------------------
;--------------------------------------

device, /close



end





;==================================================================================


