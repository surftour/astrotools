;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Compute Surface BRIGHTNESS Profiles
;     ---------------------------------------------
;     This uses the above procedures and the new color
;    code we've developed (with the help of Brant and Paul).
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------





; =============================================================================







; procedure to do most of the work
; ---------------------------------
pro process_and_plot_brightnessprofile, N, Time, x, y, z, m, age, zmets, $
			xmin, xmax, bins, $
			linecolor=linecolor, includeerr=includeerr, $
			fitdevac=fitdevac, $
			fitsersic=fitsersic, $
			x_is_devac=x_is_devac, x_is_log=x_is_log, $
			oneproj=oneproj

	print, "processing N=",N," particles"

	if N le 0 then return

	m= m*1.0e+10    ; needs conversion to solar masses

	;
	; first get luminosities, will need
	; to load the determine_lums.pro
	; procedure
	;
	;load_stellar_luminosities, N, Time, m, age, zmets, $
	;			Lum_Bol, Lum_u, Lum_b, Lum_v, Lum_k, /notsolar
	;load_stellar_luminosities, N, Time, m, age, zmets, $
        ;                       Lum_Bol, Lum_u, Lum_b, Lum_v, Lum_k

        ; get the luminosities
        ;  - in units of solar luminosities
        print, "load luminosities"
        load_all_stellar_luminosities, N, Time, m, age, zmets, $
                                Lum_Bol, Lum_U, Lum_B, Lum_V, Lum_R, Lum_I, Lum_J, Lum_H, Lum_K, $
                                Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z   ;, $
                                ;/notsolar

        ; trap for any NaN's
        ;idx=where(finite(Lum_B) eq 0)
        idx=where(finite(Lum_K) eq 0)
        if idx(0) ne -1 then begin
                Lum_Bol(idx)= 100.0 & Lum_U(idx)= 100.0 & Lum_B(idx)= 100.0 & Lum_V(idx)= 100.0 & Lum_R(idx)= 100.0
                Lum_I(idx)= 100.0 & Lum_J(idx)= 100.0 & Lum_H(idx)= 100.0 & Lum_K(idx)= 100.0
                Lum_sdss_u(idx)= 100.0 & Lum_sdss_g(idx)= 100.0 & Lum_sdss_r(idx)= 100.0
                Lum_sdss_i(idx)= 100.0 & Lum_sdss_z(idx)= 100.0
        print, "Fixed ",n_elements(idx), " luminosities when trapping for NaN values."
        endif

        print, "Bolometric luminosity= ", total(Lum_Bol)
        print, "Total K band luminosity= ", total(Lum_K)
        print, "Total B band luminosity= ", total(Lum_B)




	; print B band
	;process_and_plot_onebandprofile, x, y, z, Lum_B, xmin, xmax, bins, $
	;				linecolor=linecolor, includeerr=includeerr, $
	;				x_is_devac=x_is_devac

	; print K band
	process_and_plot_onebandprofile, x, y, z, Lum_K, xmin, xmax, bins, $
					linecolor=linecolor, includeerr=includeerr, $
					fitdevac=fitdevac, $
					fitsersic=fitsersic, x_is_log=x_is_log, $
					x_is_devac=x_is_devac, oneproj=oneproj

end






; =============================================================================






; OK, not process the profile
; ----------------------------
pro process_and_plot_onebandprofile, x, y, z, lum_band, xmin, xmax, bins, $
				linecolor=linecolor, includeerr=includeerr, $
				fitdevac=fitdevac, $
				fitsersic=fitsersic, x_is_log=x_is_log, $
				x_is_devac=x_is_devac, oneproj=oneproj

        r_avg= fltarr(bins)
        mag_avg= fltarr(bins)
        mag_1sig= fltarr(bins)

	if keyword_set(x_is_devac) then begin
	     if keyword_set(oneproj) then begin
		a= sqrt(x*x + y*y)
                c= lum_band

		process_devac_profile, a, c, bins, r_avg, mag_avg, $
					xmin, xmax, mag_1sig
	     endif else begin
		tempxmin= xmin
		tempxmax= xmax
        	;process_averages_formanyprojections, x, y, z, lum_band, $
        	process_prof_frommanyprojections, x, y, z, lum_band, $
                                        tempxmin, tempxmax, bins, $
                                        mag_avg, mag_1sig, r_avg, /x_is_devac
	     endelse
	endif else begin
	     if keyword_set(oneproj) then begin
		r= sqrt(x*x + y*y)
print, "Radii  min/max= ", min(r), max(r)
		idx= where(r gt 0.0)
		a= alog10(r(idx))
                c= lum_band(idx)

		process_prof_sd, a, c, bins, xmax, xmin, r_avg, mag_avg, $
				mass_sd_nlog_1sig=mag_1sig, /x_is_log
stop

	     endif else begin
		if x_is_log eq 1 then begin
			tempxmin= xmin
			tempxmax= xmax
		endif else begin
			tempxmin= alog10(xmin)
			tempxmax= alog10(xmax)
		endelse
        	;process_averages_formanyprojections, x, y, z, lum_band, $
        	process_prof_frommanyprojections, x, y, z, lum_band, $
                                        tempxmin, tempxmax, bins, $
                                        mag_avg, mag_1sig, r_avg, /x_is_log
	     endelse
	endelse

        idx=where(mag_avg gt 0)
        if idx(0) ne -1 then begin
                mag_avg=mag_avg(idx)
                mag_1sig=mag_1sig(idx)
                r_avg= r_avg(idx)
        endif

	; this changes it to lum/pc2
        mag_1sig= mag_1sig * 1.0e-6
        mag_avg= mag_avg * 1.0e-6

	; now change to magnitudes
        mag_p1sig= -2.5 * alog10(mag_avg+mag_1sig)
        mag_m1sig= -2.5 * alog10(mag_avg-mag_1sig)
        mag_avg= -2.5 * alog10(mag_avg)

	; this changes it to mag/arcsec2
	mag_p1sig= mag_p1sig + 21.572
	mag_m1sig= mag_m1sig + 21.572
	mag_avg= mag_avg + 21.572

	;if not keyword_set(x_is_devac) and not keyword_set(x_is_devac) then r_avg= 10^(r_avg)

        oplot, r_avg, mag_avg, psym=-3, thick=4.0, color= linecolor

        if keyword_set(includeerr) then begin
                oplot, r_avg, mag_p1sig, thick=2.0, psym=-3, color=linecolor, linestyle=1
                oplot, r_avg, mag_m1sig, thick=2.0, psym=-3, color=linecolor, linestyle=1
        endif

	; try to fit de Vac. to it
	; -------------------------
	if keyword_set(fitdevac) then begin

                if keyword_set(x_is_devac) then r_tofit= r_avg^(4.0)
                if keyword_set(x_is_log) then r_tofit= 10^r_avg
		sd_tofit= 10^(mag_avg)
		weight= mag_1sig
		;weight= 1.0e4*abs(mass_sd_1sig/mass_sd)

		fitandoverplot_devacprofile, r_tofit, sd_tofit, weight, /ylogaxis

	endif



        ; try to fit Sersic Profile to it
        ; --------------------------------
        if keyword_set(fitsersic) then begin

                if keyword_set(x_is_devac) then r_tofit= r_avg^(4.0)
                if keyword_set(x_is_log) then r_tofit= 10^r_avg
                ;sd_tofit= 10^(mag_avg)
                sd_tofit= 10^(-1.0*mag_avg/2.5)
                weight= mag_1sig
                ;weight= 1.0e4*abs(mass_sd_1sig/mass_sd)

		;idx=where(r_tofit gt 1.0)
		;r_tofit= r_tofit(idx)
		;sd_tofit= sd_tofit(idx)
		;weight= weight(idx)
		minfitradius= 0.05

                fitandoverplot_sersicprofile, r_tofit, sd_tofit, weight, /ymagaxis, $
						x_is_devac=x_is_devac, x_is_log=x_is_log, $
						minfitradius=minfitradius

        endif




end





; =============================================================================





; should be (exactly) similar to the above
; only we allow this routine to determine the quantities
; and pass them back  (mainly for colors routine)
;------------------------------------------------------
pro process_and_return_onebandprofile, x, y, z, lum_band, xmin, xmax, bins, $
				r_avg, mag_avg

	; comment out these because we'll
	; define them outside and pass them back
        ;r_avg= fltarr(bins)
        ;mag_avg= fltarr(bins)
        mag_1sig= fltarr(bins)

	tempxmin= alog10(xmin)
	tempxmax= alog10(xmax)
        ;process_averages_formanyprojections, x, y, z, lum_band, $
        process_prof_frommanyprojections, x, y, z, lum_band, $
                                        tempxmin, tempxmax, bins, $
                                        mag_avg, mag_1sig, r_avg, /x_is_log

	; don't do all this stuff in case we 
	; run into mismatches in the radial 
	; direction when we B-V them

        ;idx=where(mag_avg gt 0)
        ;if idx(0) ne -1 then begin
        ;        mag_avg=mag_avg(idx)
        ;        mag_1sig=mag_1sig(idx)
        ;        r_avg= r_avg(idx)
        ;endif

	; this changes it to lum/pc2
        ;mag_1sig= mag_1sig * 1.0e-6
        ;mag_avg= mag_avg * 1.0e-6

	; now change to magnitudes
        ;mag_p1sig= -2.5 * alog10(mag_avg+mag_1sig)
        ;mag_m1sig= -2.5 * alog10(mag_avg-mag_1sig)
        ;mag_avg= -2.5 * alog10(mag_avg)

	; this changes it to mag/arcsec2
	;mag_p1sig= mag_p1sig + 21.572
	;mag_m1sig= mag_m1sig + 21.572
	;mag_avg= mag_avg + 21.572

	; converts back to real numbers (not log)
	r_avg= 10^(r_avg)


end





; ================================================================




; comparison calls this to process one galaxy
; at a time
; ----------------------------------------------

pro do_one_galaxy_brightnessprofile, frun, snapnum, $
			linecolor=linecolor, $
			do_newstars=do_newstars, $
			do_disk=do_disk, $
			do_bulge=do_bulge, $
			msg=msg, $
			fitdevac=fitdevac, $
			fitsersic=fitsersic, $
			oneproj=oneproj, $
			loadedsnap= loadedsnap


; make sure these are the same as the
; above calling procedure
;
bins = 50

;makeitdevac= 1
makeitdevac= 0

makeitlog= 1
;makeitlog= 0

; x is radius in kpc/h
;xmax = 20.0
xmax = 100.0
xmin = 0.05
if makeitdevac eq 1 then begin
	xmax = 2.4
	xmin = 0.3
endif
if makeitlog eq 1 then begin
	xmax = 2.0
	xmin = -1.3
endif

; y is surface brightness
;ymax = 1e+4
;ymin = 1e-1
ymax = 12.0
ymin = 25.0


        if not keyword_set(loadedsnap) then ok= fload_snapshot_bh(frun,snapnum)
	Time= fload_time(1)


       	x= fload_allstars_xyz('x')
       	y= fload_allstars_xyz('y')
       	z= fload_allstars_xyz('z')
       	m= fload_allstars_mass(1)
       	age=fload_allstars_age(1)
       	age=float(Time-age)
       	zmets=fload_allstars_z(1)
	N= fload_npart(2) + fload_npart(3) + fload_npart(4)


        ; newstar light profile
        ; --------------------------
	if keyword_set(do_newstars) then begin
           Nnewstars= fload_npart(4)
	   if Nnewstars gt 0 then begin
        	x= fload_newstars_xyz('x')
        	y= fload_newstars_xyz('y')
        	z= fload_newstars_xyz('z')
        	m= fload_newstars_mass(1)
        	age=fload_newstars_age(1)
        	age=float(Time-age)
        	zmets=fload_newstars_z(1)
		N=Nnewstars
	   endif
	endif


        ; disk light profile
        ; ---------------------
	if keyword_set(do_disk) then begin
           Ndisk= fload_npart(2)
	   if Ndisk gt 0 then begin
        	x=[fload_disk_xyz('x')]
        	y=[fload_disk_xyz('y')]
        	z=[fload_disk_xyz('z')]
        	m=[fload_disk_mass(1)]
        	dage=fload_disk_age(1)
        	age=[float(Time-dage)]
        	zmets= [fload_disk_z(1)]
        	N=Ndisk
	   endif
	endif


        ; bulge light profile
        ; ---------------------
	if keyword_set(do_bulge) then begin
           Nbulge= fload_npart(3)
	   if Nbulge gt 0 then begin
        	x=fload_bulge_xyz('x')
        	y=fload_bulge_xyz('y')
        	z=fload_bulge_xyz('z')
        	m=fload_bulge_mass(1)
        	age=0.0*m + float(Time) + 2.0
        	zmets= 0.0*m + 0.02
        	N=Nbulge
	   endif
	endif


	if makeitdevac eq 1 then begin
        	process_and_plot_brightnessprofile, N, Time, x, y, z, m, age, zmets, $
                        xmin, xmax, bins, /x_is_devac, $
			fitdevac=fitdevac, $
			fitsersic=fitsersic, $
			linecolor=linecolor, oneproj=oneproj   ;, /includeerr
	endif else begin
        	process_and_plot_brightnessprofile, N, Time, x, y, z, m, age, zmets, $
                        xmin, xmax, bins, /x_is_log, $
			fitdevac=fitdevac, $
			fitsersic=fitsersic, $
			linecolor=linecolor, oneproj=oneproj   ;, /includeerr
	endelse


	if keyword_set(msg) then begin
         ; xyouts, 0.70, 0.85, msg, /normal, size= 1.5, color= linecolor, charthick=3.0
	endif



end












;============================================================================







pro plot_surface_brightness, frun, filename=filename, snapnum=snapnum, smoothlen=smoothlen

if not keyword_set(frun) then begin
   print, "  "
   print, "plot_surface_brightness, frun, filename=filename, snapnum=snapnum, smoothlen=smoothlen"
   print, "  "
   return
endif


if not keyword_set(snapnum) then snapnum= 30

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 100

; x is radius in kpc/h
;xmax = 20.0
xmax = 100.0
xmin = 0.05

; y is surface brightness
;ymax = 1e+4
;ymin = 1e-1
ymax = 12.0
ymin = 25.0


;xtit="R (h!E-1!Nkpc)"
xtit="R (kpc)"
;ytit="U mag kpc!E-2!N"
;ytit="B mag arcsec!E-2!N"
ytit="K mag arcsec!E-2!N"


; ------------------
; Plot this up
; ------------------
!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, $
        ;/ylog, $
	/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, $
	xthick=4.0, ythick=4.0, $
	charthick=3.0, $
	xtitle=xtit, $
	ytitle=ytit, $
        /nodata


        ok= fload_snapshot_bh(frun,snapnum)

	x0= 0.75
	y0= 0.90
	xyouts, x0, y0-0.03, fload_fid(1), /normal, size= 1.5, color= 0, charthick=3.0
	xyouts, x0, y0-0.08, fload_timelbl(1,2), /normal, size= 1.5, color= 0, charthick=3.0



	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=50, msg='new stars', /do_newstars, /loadedsnap
	xyouts, x0, y0-0.13, fload_getid(frun), /normal, size= 1.5, color= 50, charthick= 2.0

	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=100, msg='disk', /do_disk, /loadedsnap
	xyouts, x0, y0-0.18, fload_getid(frun), /normal, size= 1.5, color= 100, charthick= 2.0

	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=150, msg='bulge', /do_bulge, /loadedsnap
	xyouts, x0, y0-0.23, fload_getid(frun), /normal, size= 1.5, color= 150, charthick= 2.0







; put time and label on there
; ----------------------------
if keyword_set(smoothlen) then begin
	x=[smoothlen,smoothlen]
	y=[ymin,ymax]
	oplot, x, y, linestyle=1, color= 0
endif



; done, close up shop
; --------------------
device, /close



end





; =============================================================================





; 
; do a comparison of surface brightness for
; two different simulations
;
; we're making this devac profile for the moment
;
pro plot_surface_brightness_comparison, junk, filename=filename, smoothlen=smoothlen, $
					snapnum=snapnum

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: plot_surface_brightness_comparison"
   print, "  "
   return
endif


if not keyword_set(filename) then filename= 'profile.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=1
;setup_plot_stuff, 'ps', filename=filename, colortable=4


;makeitdevac= 1
makeitdevac= 0

makeitlog= 1
;makeitlog= 0

; -----------------------------
; set up constants
; -----------------------------

bins = 100

; x is radius in kpc/h
xmax = 20.0
;xmax = 100.0
xmin = 0.05
if makeitdevac eq 1 then begin
	xmax = 2.4
	xmin = 0.3
endif
if makeitlog eq 1 then begin
	xmax = 2.0
	xmin = -1.3
endif

; y is surface brightness
;ymax = 1e+4
;ymin = 1e-1
ymax = 12.0
ymin = 25.0


;xtit="R (h!E-1!Nkpc)"
xtit="R (kpc)"
if makeitdevac eq 1 then xtit="!6[R (kpc)]!E1/4!N"
if makeitlog eq 1 then xtit="!6Log R (kpc)"
;ytit="U mag kpc!E-2!N"
;ytit="B mag arcsec!E-2!N"
ytit="K mag arcsec!E-2!N"


; ------------------
; Plot this up
; ------------------
!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, $
        ;/ylog, $
	;/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, $
	xthick=4.0, ythick=4.0, $
	charthick=3.0, $
	xtitle=xtit, $
	ytitle=ytit, $
        /nodata


	if not keyword_set(snapnum) then snapnum= 25


	; Galaxy 1
	;
	; ------------------------------------------
	;frun= "pool/vc3vc3"
	;frun= "pool/vc3bvc3b"
	;frun= "/raid4/tcox/vc3vc3b"
	;do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=40, msg=fload_getid(frun)


        ; Galaxy 2
        ;
        ; ------------------------------------------
        ;frun= "pool/vc3vc3_wBH"
        ;frun= "pool/vc3bvc3b_wBH"
	;frun= "/raid4/tcox/vc3vc3c"
	;do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=50, msg=fload_getid(frun)

  do_orbits_and_orient= 0
  if do_orbits_and_orient eq 1 then begin
	frun= "/raid4/tcox/vc3vc3d"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=60, msg=fload_getid(frun)

	frun= "/raid4/tcox/vc3vc3e"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=80, msg=fload_getid(frun)

	frun= "/raid4/tcox/vc3vc3f"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=100, msg=fload_getid(frun)

	frun= "/raid4/tcox/vc3vc3g"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=120, msg=fload_getid(frun)

	frun= "/raid4/tcox/vc3vc3h"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=140, msg=fload_getid(frun)

	frun= "/raid4/tcox/vc3vc3i"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=160, msg=fload_getid(frun)

	frun= "/raid4/tcox/vc3vc3j"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=180, msg=fload_getid(frun)

	frun= "/raid4/tcox/vc3vc3k"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=200, msg=fload_getid(frun)

	frun= "/raid4/tcox/vc3vc3l"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=220, msg=fload_getid(frun)

	frun= "/raid4/tcox/vc3vc3m"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=240, msg=fload_getid(frun)

	frun= "/raid4/tcox/vc3vc3n"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=0, msg=fload_getid(frun)

	frun= "/raid4/tcox/vc3vc3o"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=50, msg=fload_getid(frun)

	frun= "/raid4/tcox/vc3vc3p"
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor=10, msg=fload_getid(frun)
   endif

   do_larger_vcs= 0
   if do_larger_vcs eq 1 then begin
	do_one_galaxy_brightnessprofile, "/raid4/tcox/vc3vc3f", snapnum, linecolor= 150
	do_one_galaxy_brightnessprofile, "/raid4/tcox/vc3vc3h", snapnum, linecolor= 150
	do_one_galaxy_brightnessprofile, "/raid4/tcox/vc3vc3j", snapnum, linecolor= 150

	do_one_galaxy_brightnessprofile, "/raid4/tcox/vc3vc3f_no", snapnum, linecolor= 50
	do_one_galaxy_brightnessprofile, "/raid4/tcox/vc3vc3h_no", snapnum, linecolor= 50
	do_one_galaxy_brightnessprofile, "/raid4/tcox/vc3vc3j_no", snapnum, linecolor= 50
   endif


   do_nobh= 0
   if do_nobh eq 1 then begin
	do_one_galaxy_brightnessprofile, "/raid4/tcox/vc3vc3e_2", 107, linecolor= 150, /oneproj
	do_one_galaxy_brightnessprofile, "/raid4/tcox/vc3vc3e_no", 107, linecolor= 50, /oneproj
	xyouts, 0.55, 0.85, '!6black hole', size=1.5, color=150, /normal, charthick=4.0
	xyouts, 0.55, 0.80, 'no black hole ', size=1.5, color=50, /normal, charthick=4.0
   endif


   do_vc3_time= 0
   if do_vc3_time eq 1 then begin
	frun= "/raid4/tcox/vc3vc3h"
	snapnum= 15
	initlc = 10
	do_one_galaxy_brightnessprofile, frun, snapnum, linecolor= initlc, /oneproj
	;do_one_galaxy_brightnessprofile, frun, snapnum+1, linecolor= initlc+10, /oneproj
	do_one_galaxy_brightnessprofile, frun, snapnum+2, linecolor= initlc+20, /oneproj
	;do_one_galaxy_brightnessprofile, frun, snapnum+3, linecolor= initlc+30, /oneproj
	do_one_galaxy_brightnessprofile, frun, snapnum+4, linecolor= initlc+40, /oneproj
	;do_one_galaxy_brightnessprofile, frun, snapnum+5, linecolor= initlc+50, /oneproj
	do_one_galaxy_brightnessprofile, frun, snapnum+6, linecolor= initlc+60, /oneproj
	;do_one_galaxy_brightnessprofile, frun, snapnum+7, linecolor= initlc+70, /oneproj
	do_one_galaxy_brightnessprofile, frun, snapnum+8, linecolor= initlc+80, /oneproj
	;do_one_galaxy_brightnessprofile, frun, snapnum+9, linecolor= initlc+90, /oneproj
	do_one_galaxy_brightnessprofile, frun, snapnum+10, linecolor= initlc+100, /oneproj
	;do_one_galaxy_brightnessprofile, frun, snapnum+11, linecolor= initlc+110, /oneproj
	do_one_galaxy_brightnessprofile, frun, snapnum+12, linecolor= initlc+120, /oneproj
	;do_one_galaxy_brightnessprofile, frun, snapnum+13, linecolor= initlc+130, /oneproj
	;do_one_galaxy_brightnessprofile, frun, snapnum+14, linecolor= initlc+140, /oneproj
	do_one_galaxy_brightnessprofile, frun, snapnum+15, linecolor= initlc+150, /oneproj
	xyouts, 0.55, 0.85, '!6time course', size=1.5, color=0, /normal, charthick=4.0
	xyouts, 0.55, 0.80, 'dark -> ', size=1.5, color=20, /normal, charthick=4.0
	xyouts, 0.75, 0.80, ' light', size=1.5, color=160, /normal, charthick=4.0
   endif


   do_vc3_nobhtime= 0
   if do_vc3_nobhtime eq 1 then begin
	;frun= "/raid4/tcox/vc3vc3e_no"
	initlc = 10
	;do_one_galaxy_brightnessprofile, frun, 86, linecolor= initlc, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 95, linecolor= initlc+20, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 98, linecolor= initlc+40, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 99, linecolor= initlc+60, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 101, linecolor= initlc+80, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 103, linecolor= initlc+100, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 105, linecolor= initlc+120, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 107, linecolor= initlc+150, /oneproj, /fitsersic

	;frun= "/raid4/tcox/ds/d6e"
	;do_one_galaxy_brightnessprofile, frun, 10, linecolor= initlc, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 12, linecolor= initlc+20, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 15, linecolor= initlc+40, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 18, linecolor= initlc+60, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 20, linecolor= initlc+80, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 22, linecolor= initlc+100, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 28, linecolor= initlc+120, /oneproj, /fitsersic
	;do_one_galaxy_brightnessprofile, frun, 30, linecolor= initlc+150, /oneproj, /fitsersic

	frun= "/raid4/tcox/ds/d0e"
	do_one_galaxy_brightnessprofile, frun, 30, linecolor= initlc, /oneproj, /fitsersic
	do_one_galaxy_brightnessprofile, frun, 32, linecolor= initlc+20, /oneproj, /fitsersic
	do_one_galaxy_brightnessprofile, frun, 35, linecolor= initlc+40, /oneproj, /fitsersic
	do_one_galaxy_brightnessprofile, frun, 38, linecolor= initlc+60, /oneproj, /fitsersic
	do_one_galaxy_brightnessprofile, frun, 40, linecolor= initlc+80, /oneproj, /fitsersic
	do_one_galaxy_brightnessprofile, frun, 42, linecolor= initlc+100, /oneproj, /fitsersic
	do_one_galaxy_brightnessprofile, frun, 48, linecolor= initlc+120, /oneproj, /fitsersic
	do_one_galaxy_brightnessprofile, frun, 50, linecolor= initlc+150, /oneproj, /fitsersic

	;xyouts, 0.55, 0.90, '!6no black hole', size=1.5, color=0, /normal, charthick=4.0
	xyouts, 0.55, 0.85, '!6time course', size=1.5, color=0, /normal, charthick=4.0
	xyouts, 0.55, 0.80, 'dark -> ', size=1.5, color=20, /normal, charthick=4.0
	xyouts, 0.75, 0.80, ' light', size=1.5, color=160, /normal, charthick=4.0
   endif


   do_vc3_gf= 0
   if do_vc3_gf eq 1 then begin
	do_one_galaxy_brightnessprofile, "/raid4/tcox/vc3vc3", 19, linecolor= 50
	do_one_galaxy_brightnessprofile, "/raid4/tcox/vc3vc3h", 19, linecolor= 100
	do_one_galaxy_brightnessprofile, "/raid4/tcox/As/A3", 140, linecolor= 150
   endif


   do_brant= 1
   if do_brant eq 1 then begin
	do_one_galaxy_brightnessprofile, "brantprofs/A", 99, linecolor= 50, /oneproj, /fitsersic
   endif





; put time and label on there
; ----------------------------
if keyword_set(smoothlen) then begin
	x=[smoothlen,smoothlen]
	if makeitdevac eq 1 then x=x^(0.25)
	if makeitlog eq 1 then x=x^(0.25)
	y=[ymin,ymax]
	oplot, x, y, linestyle=1, color= 0
endif


;xyouts, 0.70, 0.85, fload_timelbl(1,2), /normal, size= 1.5, color= 0, charthick=3.0


; done, close up shop
; --------------------
device, /close



end







; =============================================================================







