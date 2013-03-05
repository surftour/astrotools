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



pro do_seq, frun

   starti= 0
   ;starti= 1
   ;starti= 35


   spawn, "/bin/ls "+frun+"/snap* | wc ",result
   endi=long(result[0])-1
   ;endi= 5

   print, "frun= ", frun
   print, "starti= ", starti
   print, "endi= ", endi

   for i=starti,endi do begin

	thisi= i

	exts='0000'+strcompress(string(thisi),/remove_all)
	;thisfile= frun+'/toddt/toddt_'+strmid(exts,strlen(exts)-3,3)+'.eps'
	thisfile= frun+'/toddt/toddt_'+strmid(exts,strlen(exts)-3,3)+'.png'
	tempfile= 'temp.eps'

	;print, "toddt_9, "+frun+", "+string(i)+", filename="+thisfile
	toddt_9, frun, thisi, filename=tempfile

	cmd="gs -sDEVICE=png16m -dSAFER -dGraphicsAlphaBits=4 -r150 -dBATCH -dNOPAUSE -dEPSCrop -sOutputFile="
	cmd=cmd+thisfile+"  temp.eps"
	spawn, cmd

   endfor

end













;
;
;   ----------------            ---------------
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   ----------------            ---------------
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   ----------------            ---------------
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   ----------------            ---------------
;
;

pro toddt_6, frun, snapnum, filename=filename
        

if not keyword_set(frun) then begin
   print, "  "
   print, " toddt_6, frun, snapnum, filename=filename, "
   print, "  "
   print, "  "
   print, "  "
   return
endif   
        
if not keyword_set(snapnum) then snapnum=12
if not keyword_set(filename) then filename="toddt6.eps"




;--------------------------------------------


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=22, newysize=18



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




;---------------------------
;  Calc it
;---------------------------
   

	ok=fload_snapshot_bh(frun,snapnum)

	center=fload_center_alreadycomp(1)
	print, "using center= ", center


	; stellar info
	m_stars= fload_allstars_mass(1)
	print, "total stellar mass= ", total(m_stars)
	m_stars= m_stars / 0.7
	r_stars=fload_allstars_xyz('r',center=center) / 0.7


	; luminosities (in solar units)
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
	process_cumulative_profile, r_stars, m_stars, this_xmin, this_xmax, bins, rs_stars, ms_stars
	print, "stars max/min= ", max(ms_stars), min(ms_stars)

	; cumulative gas mass
	r_gas= alog10(r_gas)
	;process_prof_tot, r_gas, m_gas, bins, this_xmax, this_xmin, rs_gas, ms_gas
	process_cumulative_profile, r_gas, m_gas, this_xmin, this_xmax, bins, rs_gas, ms_gas
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
	process_cumulative_profile, r_gas_diffuse, m_gas_diffuse, this_xmin, this_xmax, bins, rs_gas_diffuse, ms_gas_diffuse
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
	process_cumulative_profile, r_dark, m_dark, this_xmin, this_xmax, bins, rs_dark, ms_dark
	print, "dark mass max/min= ", max(ms_dark), min(ms_dark)

	; cumulative stellar luminosity
	;process_prof_tot, r_stars, bololum, bins, this_xmax, this_xmin, rs_stars, lum_stars
	process_cumulative_profile, r_stars, bololum, this_xmin, this_xmax, bins, rs_stars, lum_stars
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


	; total luminosity (in solar units)
	lum_total= lum_stars
	for i= 0, bins-1 do begin
	   lum_total[i]= total(lum_stars[i] + bhlum)
	endfor


	; trial 1
	;lum_total= lum_total * solarlum_in_gadget_units ; gadget units
	;lum_total_cgs= lum_total * UnitSolarLuminosity_in_ergs
	; trial 2
	lum_total_lsun= lum_total
	lum_total_cgs= lum_total * UnitSolarLuminosity_in_ergs
	lum_total= lum_total_cgs / UnitEnergy_in_cgs * UnitTime_in_s



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


; the above (cgs vs. non) are off by exactly 1684.4713 (exactly 1/solarlum_in_gadget_units)

	; now plot
	rs= 10^(rs_stars)




;---------------------------
;  Now, plot it
;---------------------------
   

x0= 0.10 & xsize= 0.36
x1= x0 + xsize

; put the same x0 inbetween the two columns
x2= 2.*x0 + xsize + 0.03
x3= x2 + xsize


y0= 0.09 & ysize=0.30
y1= y0+ysize
y2= y0+ysize+ysize
y3= y0+ysize+ysize+ysize



; left-hand panel
;-----------------

; top (cumulative mass)
;
;ymax = 50
;ymin = 0.07
mass_total= 10.0 + alog10(mass_total)
ymax= max(mass_total) + 0.5
ymin= min(mass_total) - 0.5
plot_panel, rs, mass_total, x0, y2, x1, y3, xmax, xmin, ymax, ymin, $
                yaxistitle='Log Mass (<r) [M!9n!6]', /xlog, pointt=-0

mass_dark= 10.0 + alog10(ms_dark)
oplot, rs, mass_dark, psym=-3, linestyle= 1, color= 50, thick=1.0
mass_gas= 10.0 + alog10(ms_gas)
oplot, rs, mass_gas, psym=-3, linestyle= 2, color= 100, thick=1.0
mass_stars= 10.0 + alog10(ms_stars)
oplot, rs, mass_stars, psym=-3, linestyle= 2, color= 150, thick=1.0

xyouts, x0+0.02, y3-0.03, "Total", /normal, color= 0, charthick=1.5, size=1.0
xyouts, x0+0.02, y3-0.06, "Dark Matter", /normal, color= 50, charthick=1.5, size=1.0
xyouts, x0+0.02, y3-0.09, "Stars", /normal, color= 150, charthick=1.5, size=1.0
xyouts, x0+0.02, y3-0.12, "Gas", /normal, color= 100, charthick=1.5, size=1.0


; middle (gas density)
;
ymax =  1.4
ymin = -3.2
plot_panel, rs, rho_gas, x0, y1, x1, y2, xmax, xmin, ymax, ymin, $
                ;yaxistitle='!7q!6 (cm!E-3!N)', /xlog, /ylog, pointt=-50
                yaxistitle='Log !7q!6 (cm!E-3!N)', /xlog, pointt=-50

; bottom (luminosity)
;
;ymax = 1.0e+12
;ymin = 1.0e+8
;ymax= 13.2
;ymin= 9.8
lum_total_lsun= alog10(lum_total_lsun)
ymax= max(lum_total_lsun) + 0.5
ymin= min(lum_total_lsun) - 0.5
plot_panel, rs, lum_total_lsun, x0, y0, x1, y1, xmax, xmin, ymax, ymin, $
                yaxistitle='Log L!Dtot!N (<r) [L!9n!6]', /xlog, pointt=-50, $
		xaxistitle=xaxistitle

bh_to_total= total(bhlum) / (total(bhlum) + total(bololum))
lblz = strcompress(string(bh_to_total),/remove_all)
lblz = strmid(lblz,0,5)        ; 0.xxx
xyouts, x0+0.02, y1-0.04, "L!DBH!N/L!Dtot!N = "+lblz, /normal, color= 0, charthick=1.5, size=1.2

lblz = strcompress(string(total(fload_gas_sfr(1))),/remove_all)
lblz = strmid(lblz,0,5)        ; 0.xxx
xyouts, x0+0.02, y1-0.08, "SFR = "+lblz+" M!D!9n!6!N Yr!E-1!N", /normal, color= 0, charthick=1.5, size=1.2





; right-hand panel
;-----------------

; top (L_tot/L_edd)
;
ymax = 50.0
ymin = 0.0005
plot_panel, rs, ltotledd, x2, y2, x3, y3, xmax, xmin, ymax, ymin, $
                yaxistitle='!6L!Dtot!N/L!DEdd!N', /xlog, /ylog, pointt=-50

; middle (sigma)
;
ymax = 5.00
ymin = 0.05
plot_panel, rs, sigma_gas, x2, y1, x3, y2, xmax, xmin, ymax, ymin, $
                yaxistitle='!7R!6 (g cm!E-2!N)', /xlog, /ylog, pointt=-50

; bottom (kappa)
;
ymax = 10
ymin = 0.1
plot_panel, rs, kappa, x2, y0, x3, y1, xmax, xmin, ymax, ymin, $
                yaxistitle='!7j!6 (cm!E2!N/g)', /xlog, /ylog, pointt=-50, $
		xaxistitle=xaxistitle




;xyouts, 0.55, 0.38, "!7j!6 = !7R!6!E-1!N ", /normal, color= 0, charthick=2.0, size=1.2
;xyouts, 0.55, 0.33, "std", /normal, color= 50, charthick=2.0, size=1.2



;------------------
device, /close

end




;==================================================================================





;
;
;   ----------------            ---------------       -----------------------------
;   |              |            |             |       |              |            |
;   |              |            |             |       |              |            |
;   |              |            |             |       |              |            |
;   |              |            |             |       |              |            |
;   |              |            |             |       |              |            |
;   |              |            |             |       |              |            |
;   ----------------            ---------------       -----------------------------
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |            ------------------
;   |              |            |             |            |                |
;   |              |            |             |            |                |
;   ----------------            ---------------            |                |
;   |              |            |             |            |                |
;   |              |            |             |            |                |
;   |              |            |             |            ------------------
;   |              |            |             |
;   |              |            |             |
;   |              |            |             |
;   ----------------            ---------------
;
;

pro toddt_9, frun, snapnum, filename=filename
        

if not keyword_set(frun) then begin
   print, "  "
   print, " toddt_9, frun, snapnum, filename=filename, "
   print, "  "
   print, "  "
   print, "  "
   return
endif   
        
if not keyword_set(snapnum) then snapnum=0
if not keyword_set(filename) then filename="toddt9.eps"




;--------------------------------------------


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=35, newysize=20



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




;---------------------------
;  Calc it
;---------------------------
   

	ok=fload_snapshot_bh(frun,snapnum,/nopot_in_snap)

	; option 1 - use computed center
	center=fload_center_alreadycomp(1)
	orig_center= center
	print, "alreadycomp center= ", center

	; option 2 - use BH center
        bhid= max(fload_blackhole_id(1))
        ;bhid= bhid[0]
        ;bhid= bhid[1]
        ;bhid= 200001L    ; used for ds/vc3vc3h_2
        ;bhid= 280002L   ; used for z3/b4e
        ;bhid= 400002L   ; used for ds/vc3vc3e_2
        center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)

        print, "Blackhole ID: ", bhid
        print, "Blackhole center: ", center_bh
	orig_center= center_bh


	; stellar info
	center= orig_center
	m_stars= fload_allstars_mass(1)
	print, "total stellar mass= ", total(m_stars)
	m_stars= m_stars / 0.7
	r_stars=fload_allstars_xyz('r',center=center) / 0.7


	; luminosities (in solar units)
	bololum=fload_allstars_bololum(1)
	TTime= float(fload_time(1))
	bhlum= fload_blackhole_lum(frun,Ttime,/bolometric)


	; gaseous info
	center=orig_center
	m_gas= fload_gas_mass(1)
	print, "total gas mass= ", total(m_gas)
	m_gas= m_gas / 0.7
	r_gas= fload_gas_xyz('r', center=center) / 0.7


	; diffuse gas info
	cf= fload_gas_coldfraction(1)
	m_gas_diffuse= m_gas * (1.0 - cf)
	r_gas_diffuse= r_gas


	; dark mass info
	center=orig_center
	m_dark= fload_halo_mass(1)
	print, "total dark mass= ", total(m_dark)
	m_dark= m_dark / 0.7
	r_dark= fload_halo_xyz('r', center=center) / 0.7



	; now, calculate profiles
	; --------------------------
	this_xmax= alog10(xmax)
	this_xmin= alog10(xmin)



	;---------------------------
	; cumulative stellar mass
	r_stars= alog10(r_stars)
	;process_prof_tot, r_stars, m_stars, bins, this_xmax, this_xmin, rs_stars, ms_stars
	process_cumulative_profile, r_stars, m_stars, this_xmin, this_xmax, bins, rs_stars, ms_stars
	print, "stars max/min= ", max(ms_stars), min(ms_stars)

	;---------------------------
	; cumulative gas mass
	r_gas= alog10(r_gas)
	;process_prof_tot, r_gas, m_gas, bins, this_xmax, this_xmin, rs_gas, ms_gas
	process_cumulative_profile, r_gas, m_gas, this_xmin, this_xmax, bins, rs_gas, ms_gas
	print, "gas max/min= ", max(ms_gas), min(ms_gas)


	;---------------------------
	; gas density (v1)
	radius= r_gas
	mass= m_gas
	rho_gas= fltarr(bins)
	radius_gas= fltarr(bins)
	weight= fltarr(bins)
	process_prof_rho, radius, mass, bins, this_xmax, this_xmin,  radius_gas, rho_gas, weight, /r_is_log
	print, "gas rho (v1) max/min= ", max(rho_gas), min(rho_gas)

	; this is the straight-up conversion to cm^-3
        rho_gas_cgs= rho_gas * 404.5


	;---------------------------
	; gas density (v2)
	radius= r_gas
	rho= fload_gas_rho(1)
	rho_v2_gas= fltarr(bins)
	radius_v2_gas= fltarr(bins)
	weight_v2= fltarr(bins)
	; returns values in log
	process_prof_avg, radius, rho, bins, this_xmax, this_xmin, radius_v2_gas, rho_v2_gas, weight_v2, /r_is_log, /log_avg
	print, "gas rho (v2) max/min= ", max(rho_v2_gas), min(rho_v2_gas)

	idx= where(rho_v2_gas ne 0.0)    ; hopefully no density is exactly 0.0 (remember that this is log)
	if idx(0) ne -1 then begin
		rho_v2_gas= rho_v2_gas(idx)
		weight_v2= weight_v2(idx)
		radius_v2_gas= radius_v2_gas(idx)
	endif

	; this is the straight-up conversion to cm^-3
        rho_v2_gas_cgs= rho_v2_gas + alog10(404.5)
	rho_v2_gas_cgs_plus= rho_v2_gas + weight_v2 + alog10(404.5)
	rho_v2_gas_cgs_minus= rho_v2_gas - weight_v2 + alog10(404.5)


	;---------------------------
	; gas density (v3)
	radius= r_gas
	rho= fload_gas_rho(1)
	rho_v3_gas= fltarr(bins)
	radius_v3_gas= fltarr(bins)
	weight_v3= fltarr(bins)
	; returns values in log
	process_prof_avg, radius, rho, bins, this_xmax, this_xmin, radius_v3_gas, rho_v3_gas, weight_v3, /r_is_log
	print, "gas rho (v3) max/min= ", max(rho_v3_gas), min(rho_v3_gas)

	idx= where(rho_v3_gas gt 0.0)    ; hopefully no density is exactly 0.0 (remember that this is log)
	if idx(0) ne -1 then begin
		rho_v3_gas= rho_v3_gas(idx)
		weight_v3= weight_v3(idx)
		radius_v3_gas= radius_v3_gas(idx)
	endif

	; this is the straight-up conversion to cm^-3
        rho_v3_gas_cgs= rho_v3_gas * 404.5
	rho_v3_gas_cgs_plus= (rho_v3_gas + weight_v3) * 404.5
	rho_v3_gas_cgs_minus= (rho_v3_gas - weight_v3) * 404.5
	idx= where(rho_v3_gas_cgs_minus le 0.0)
	if idx(0) ne -1 then begin
		rho_v3_gas_cgs_minus(idx)= rho_v3_gas_cgs(idx)/1000.0
	endif
	rho_v3_gas_cgs= alog10(rho_v3_gas_cgs)
	rho_v3_gas_cgs_plus= alog10(rho_v3_gas_cgs_plus)
	rho_v3_gas_cgs_minus= alog10(rho_v3_gas_cgs_minus)



	;---------------------------
	; cumulative diffuse gas mass
	r_gas_diffuse= alog10(r_gas_diffuse)
	process_cumulative_profile, r_gas_diffuse, m_gas_diffuse, this_xmin, this_xmax, bins, rs_gas_diffuse, ms_gas_diffuse
	print, "diffuse gas max/min= ", max(ms_gas_diffuse), min(ms_gas_diffuse)

	;---------------------------
	; diffuse gas density
	radius_diffuse= r_gas_diffuse
	mass_diffuse= m_gas_diffuse
	rho_gas_diffuse= fltarr(bins)
	radius_gas_diffuse= fltarr(bins)
	weight_diffuse= fltarr(bins)
	process_prof_rho, radius_diffuse, mass_diffuse, bins, this_xmax, this_xmin,  radius_gas_diffuse, rho_gas_diffuse, weight_diffuse, /r_is_log
	print, "diffuse gas rho max/min= ", max(rho_gas_diffuse), min(rho_gas_diffuse)

	; this is the straight-up conversion to cm^-3
        rho_gas_diffuse_cgs= rho_gas_diffuse * 404.5

	;---------------------------
	; cumulative dark mass
	r_dark= alog10(r_dark)
	;process_prof_tot, r_dark, m_dark, bins, this_xmax, this_xmin, rs_dark, ms_dark
	process_cumulative_profile, r_dark, m_dark, this_xmin, this_xmax, bins, rs_dark, ms_dark
	print, "dark mass max/min= ", max(ms_dark), min(ms_dark)

	;---------------------------
	; cumulative stellar luminosity
	;process_prof_tot, r_stars, bololum, bins, this_xmax, this_xmin, rs_stars, lum_stars
	process_cumulative_profile, r_stars, bololum, this_xmin, this_xmax, bins, rs_stars, lum_stars
	print, "stellar luminosity (solar masses) max/min= ", max(lum_stars), min(lum_stars)



	;
	;  units and other brew-ha-ha
	; ------------------------------------

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


	; total luminosity (in solar units)
	lum_total= lum_stars
	for i= 0, bins-1 do begin
	   lum_total[i]= total(lum_stars[i] + bhlum)
	endfor


	; trial 1
	;lum_total= lum_total * solarlum_in_gadget_units ; gadget units
	;lum_total_cgs= lum_total * UnitSolarLuminosity_in_ergs
	; trial 2
	lum_total_lsun= lum_total
	lum_total_cgs= lum_total * UnitSolarLuminosity_in_ergs
	lum_total= lum_total_cgs / UnitEnergy_in_cgs * UnitTime_in_s



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


; the above (cgs vs. non) are off by exactly 1684.4713 (exactly 1/solarlum_in_gadget_units)

	; now plot
	rs= 10^(rs_stars)




;---------------------------
;  Now, plot it
;---------------------------
   

x0= 0.08 & xsize= 0.18
x1= x0 + xsize

; put the same x0 inbetween the two columns
x2= 2.*x0 + xsize + 0.01
x3= x2 + xsize


y0= 0.09 & ysize=0.30
y1= y0+ysize
y2= y0+ysize+ysize
y3= y0+ysize+ysize+ysize



; left-hand panel
;-----------------

; top (cumulative mass)
;
;ymax = 50
;ymin = 0.07
ymax = 11.55
ymin = 7.07
mass_total= 10.0 + alog10(mass_total)
;ymax= max(mass_total) + 0.5
;ymin= min(mass_total) - 0.5
plot_panel, rs, mass_total, x0, y2, x1, y3, xmax, xmin, ymax, ymin, $
                yaxistitle='Log Mass (<r) [M!D!9n!6!N]', /xlog, pointt=-0

mass_dark= 10.0 + alog10(ms_dark)
oplot, rs, mass_dark, psym=-3, linestyle= 0, color= 50, thick=5.0
mass_gas= 10.0 + alog10(ms_gas)
oplot, rs, mass_gas, psym=-3, linestyle= 0, color= 100, thick=5.0
mass_stars= 10.0 + alog10(ms_stars)
oplot, rs, mass_stars, psym=-3, linestyle= 0, color= 150, thick=5.0

xyouts, x1-0.07, y2+0.08, "Total", /normal, color= 0, charthick=1.5, size=1.0
xyouts, x1-0.07, y2+0.06, "Dark Matter", /normal, color= 50, charthick=1.5, size=1.0
xyouts, x1-0.07, y2+0.04, "Stars", /normal, color= 150, charthick=1.5, size=1.0
xyouts, x1-0.07, y2+0.02, "Gas", /normal, color= 100, charthick=1.5, size=1.0


; middle (gas density)
;
;ymax = 50
;ymin = 0.0005
ymax =  3.4
ymin = -3.2
idx=where(rho_gas_cgs gt 0.0)
if idx(0) ne -1 then begin
	rho_gas_cgs= rho_gas_cgs(idx)
	radius_gas= radius_gas(idx)
endif
rho_gas_cgs= alog10(rho_gas_cgs)
radius_gas= 10^(radius_gas)
plot_panel, radius_gas, rho_gas_cgs, x0, y1, x1, y2, xmax, xmin, ymax, ymin, $
                ;yaxistitle='!7q!6 (cm!E-3!N)', /xlog, /ylog, pointt=-50
                yaxistitle='Log !7q!6!Dgas!N (cm!E-3!N)', /xlog, pointt=-50

; diffuse component
idx=where(rho_gas_diffuse_cgs gt 0.0)
if idx(0) ne -1 then begin
        rho_gas_diffuse_cgs= rho_gas_diffuse_cgs(idx)
        radius_gas_diffuse= radius_gas_diffuse(idx)
endif
rho_gas_diffuse_cgs= alog10(rho_gas_diffuse_cgs)
radius_gas_diffuse= 10^(radius_gas_diffuse)
oplot, radius_gas_diffuse, rho_gas_diffuse_cgs, psym=-3, color= 200, thick= 3.0
;oplot, radius_gas_diffuse, rho_gas_diffuse_cgs, psym=2, linestyle= 2, color= 200, thick= 5.0, symsize= 0.7

xyouts, x0+0.02, y1+0.05, "total", /normal, color= 50, charthick=1.5, size=1.0
xyouts, x0+0.02, y1+0.03, "`diffuse-phase' only", /normal, color= 200, charthick=1.5, size=1.0


; v2
;radius_v2_gas= 10^(radius_v2_gas)
;oplot, radius_v2_gas, rho_v2_gas_cgs, psym=-3, linestyle= 0, color= 100, thick= 5.0
;oplot, radius_v2_gas, rho_v2_gas_cgs_plus, psym=-3, linestyle= 1, color= 100, thick= 3.0
;oplot, radius_v2_gas, rho_v2_gas_cgs_minus, psym=-3, linestyle= 1, color= 100, thick= 3.0

; v3
;radius_v3_gas= 10^(radius_v3_gas)
;oplot, radius_v3_gas, rho_v3_gas_cgs, psym=-3, linestyle= 0, color= 150, thick= 5.0
;oplot, radius_v3_gas, rho_v3_gas_cgs_plus, psym=-3, linestyle= 1, color= 150, thick= 3.0
;oplot, radius_v3_gas, rho_v3_gas_cgs_minus, psym=-3, linestyle= 1, color= 150, thick= 3.0


; bottom (luminosity)
;
;ymax = 1.0e+12
;ymin = 1.0e+8
ymax= 12.8
ymin= 9.6
lum_total_lsun= alog10(lum_total_lsun)
;ymax= max(lum_total_lsun) + 0.5
;ymin= min(lum_total_lsun) - 0.5
plot_panel, rs, lum_total_lsun, x0, y0, x1, y1, xmax, xmin, ymax, ymin, $
                yaxistitle='Log L!Dtot!N (<r) [L!D!9n!6!N]', /xlog, pointt=-50, $
		xaxistitle=xaxistitle

bh_to_total= total(bhlum) / (total(bhlum) + total(bololum))
lblz = strcompress(string(bh_to_total),/remove_all)
lblz = strmid(lblz,0,5)        ; 0.xxx
xyouts, x0+0.02, y1-0.04, "L!DBH!N/L!Dtot!N = "+lblz, /normal, color= 0, charthick=1.5, size=1.2





; right-hand panel
;-----------------

; top (L_tot/L_edd)
;
ymax = 2.2
ymin = -3.85
ltotledd= alog10(ltotledd)
plot_panel, rs, ltotledd, x2, y2, x3, y3, xmax, xmin, ymax, ymin, $
                ;yaxistitle='!6L!Dtot!N/L!DEdd!N', /xlog, /ylog, pointt=-50
                yaxistitle='Log !6L!Dtot!N/L!DEdd!N', /xlog, pointt=-50

x= [xmin,xmax] & y= [0.0, 0.0]
oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0

ltotledd_diffuse= alog10(ltotledd_diffuse)
oplot, rs, ltotledd_diffuse, psym=-3, linestyle= 0, color= 200, thick= 3.0


; middle (sigma)
;
ymax = 0.3
ymin = -5.4
sigma_gas= alog10(sigma_gas)
plot_panel, rs, sigma_gas, x2, y1, x3, y2, xmax, xmin, ymax, ymin, $
                yaxistitle='Log !7R!6 (g cm!E-2!N)', /xlog, pointt=-50

sigma_gas_diffuse= alog10(sigma_gas_diffuse)
oplot, rs, sigma_gas_diffuse, psym=-3, linestyle= 0, color= 200, thick= 3.0


; bottom (kappa)
;
ymax = 4.4
ymin = -0.2
kappa= alog10(kappa)
plot_panel, rs, kappa, x2, y0, x3, y1, xmax, xmin, ymax, ymin, $
                yaxistitle='Log !7j!6 (cm!E2!N g!E-1!N)', /xlog, pointt=-50, $
		xaxistitle=xaxistitle

kappa_diffuse= alog10(kappa_diffuse)
oplot, rs, kappa_diffuse, psym=-3, linestyle= 0, color= 200, thick= 3.0




; now images
;-----------------

xlen= 10.0

x0= 0.56
x1= 0.77
x2= 0.98

y0= 0.58
y1= 0.96


;do_gas= 0
do_gas= 1

if do_gas eq 1 then begin
    npart= fload_npart(0)
    N= long(npart)

    center= orig_center
    x= fload_gas_xyz('x',center=center) / 0.7
    y= fload_gas_xyz('y',center=center) / 0.7
    z= fload_gas_xyz('z',center=center) / 0.7
    m= fload_gas_mass(1) / 0.7
    hsml= fload_gas_hsml(1) / 0.7

endif



; ----------------------------------------------
;  Call our generic contour plotting procedure
; ----------------------------------------------

;set_maxden= 10.0  
set_maxden= 1.0    ; 10^4 msolar/pc2
;set_dynrng= 1.0e+3
;set_dynrng= 1.0e+4
set_dynrng= 1.0e+5
;set_dynrng= 1.0e+6
;set_maxden= 0.5e-1
;set_dynrng= 1.0e+6


center= [0.0, 0.0, 0.0]

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
        hsml=hsml, $
        xthickness=xthickness, ythickness=ythickness, $
        pixels=pixels, zthickness=zthickness, $
        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        crude=crude, center=center, $
        set_maxden= set_maxden, $
        set_dynrng= set_dynrng, $
        NxNImage=NxNImage

        XYImage= NxNImage


xz= 1
center= [0.0, 0.0, 0.0]

contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
        hsml=hsml, $
        xthickness=xthickness, ythickness=ythickness, $
        pixels=pixels, zthickness=zthickness, $
        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        crude=crude, center=center, $
        set_maxden= set_maxden, $
        set_dynrng= set_dynrng, $
        NxNImage=NxNImage

        XZImage= NxNImage



tv, XYImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal
tv, XZImage, x1, y0, xsize=(x2-x1), ysize=(y1-y0), /normal


xmin= -xlen
xmax=  xlen
ymin= -xlen
ymax=  xlen

; 
; xy
;
;------
!p.position=[x0,y0,x1,y1]
!p.ticklen=0.03

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'

xyouts, x0+0.02, y1-0.04, '!6xy', size=1.2, /normal, color= 0, charthick=1.5

n_bh= fload_npart(5)
bhid= fload_blackhole_id(1)

for i=0, n_bh-1 do begin
	;usersym,2.0*cos(findgen(49)/49*2*!pi),2.0*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[1]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[1]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor


; 
; xy
;
;------
!p.position=[x1,y0,x2,y1]

plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
	xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
	xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'


xyouts, x1+0.02, y1-0.04, '!6xz', size=1.2, /normal, color= 0, charthick=1.5

for i=0, n_bh-1 do begin
	center= orig_center
	bhid1= bhid[i]
	bh_xyz= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid1) / 0.7
	;oplot, [bh_xyz[0]], [bh_xyz[2]], psym=8, thick=6.0, color= 0, symsize= 2.0
	oplot, [bh_xyz[0]], [bh_xyz[2]], psym=7, thick=5.0, color= 220, symsize= 2.0
endfor







; and SFR
;-----------------

; needs .run sfr_multi

x0= 0.65
x1= 0.95
y0= 0.09
y1= 0.50

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

xmin= 0.0
xmax= 3.0
ymax= 300.0
ymin= 0.05


!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, /ylog, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        ytickformat='exp_label', xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase, /normal



open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac

; physical units
sfrtime = sfrtime / 0.7

oplot, sfrtime, sfrsfr, psym=-3, color= 150, linestyle= 0, thick= 4.0

idx=where(sfrtime ge (fload_time(1)/0.7))
if idx(0) eq -1 then idx[0]= n_elements(sfrtime)-1
oplot, [sfrtime[idx[0]]], [sfrsfr[idx[0]]], psym= 4, color= 0, thick=5.0, symsize=1.6


sfrtot= total(fload_gas_sfr(1))
lblz = strcompress(string(sfrtot),/remove_all)
digz= 5
if sfrtot lt 100.0 then digz= 4
if sfrtot lt 10.0 then digz= 3
lblz = strmid(lblz,0,digz)        ; 0.xxx
xyouts, x0+0.04, y0+0.08, lblz+" M!D!9n!6!N Yr!E-1!N", /normal, color= 0, charthick=1.5, size=1.5




;------------------
device, /close

end




;==================================================================================





;==================================================================================


