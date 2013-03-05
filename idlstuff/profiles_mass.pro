
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Surface Density Plots
;   -------------------------
;
;  Two options, log and de Vaucouleurs (i.e., r^0.25)
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------


pro plot_sd_profile, frun, snapnum, $
			smoothlen=smoothlen, $
			filename=filename, $
			x_is_devac=x_is_devac, $
			x_is_linear=x_is_linear

if not keyword_set(frun) then begin
   print, "  "
   print, "plot_sd_profile, frun, snapnum, smoothlen=smoothlen, filename=filename, $"
   print, "                /x_is_devac, /x_is_linear  (default is log)"
   print, "  "
   return
endif

if not keyword_set(smoothlen) then smoothlen= 0.1
if not keyword_set(filename) then filename='devac.eps'
if not keyword_set(snapnum) then snapnum=30


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 40

yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"


; default is log
;
xaxistitle="!6Log R (kpc)"
xmax = 1.5
xmin = -1.0
ymax = 6.0
ymin = -1.0
x_is_log= 1

if keyword_set(x_is_devac) then begin
	xaxistitle="!6[R (kpc)]!E1/4!N"
	xmax = 2.4
	;xmax = 4.2
	xmin = 0.5
	ymax = 6.0
	ymin = -1.0
	x_is_log= 0
endif

if keyword_set(x_is_linear) then begin
        xaxistitle="!6 R (kpc)"
        xmax = 20.0
        xmin = 0.0
        ymax = 4.2
        ymin = -1.0
	x_is_log= 0
endif




; ------------------
; Plot this up
; ------------------
!p.position= [0.2, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, /nodata, $
        ;/ylog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle




; -----------------------------------------------


	process_and_plot_sd, frun, snapnum, "allstars", xmin, xmax, bins, linecolor=150, $
	;process_and_plot_sd, x, y, z, m, xmin, xmax, bins, linecolor=150
	;process_and_plot_sd, x, y, z, m, xmin, xmax, bins, linecolor=150, $
						x_is_devac=x_is_devac, $
	;					/includeerr, $
	;					x_is_log=x_is_log, $
						/fitsersic




        ;xyouts, 0.60, 0.86, '!6Merger Remnant', /normal, charthick=3.0, size= 1.2, color= 150
        ;xyouts, 0.65, 0.85, frun, /normal, charthick=3.0, size= 1.2, color= 0


	do_special= 0
	;do_special= 1
	if do_special eq 1 then begin


	        process_and_plot_sd, "/raid4/tcox/minor/iso_Sb", 0, "allstars", xmin, xmax, bins, linecolor=50
	        ;process_and_plot_sd, x, y, z, m, xmin, xmax, bins, linecolor=50
	        ;process_and_plot_sd, x, y, z, m, xmin, xmax, bins, linecolor=50, $
                ;                                x_is_devac=x_is_devac, $
                ;                                /includeerr, $
                ;                                x_is_log=x_is_log, $
                ;                                /fitsersic

        	xyouts, 0.60, 0.80, '!6Initial Disk Profile', /normal, charthick=3.0, size= 1.2, color= 50
	endif



; extras
; --------

if keyword_set(devac) then begin
    ok = oplot_smoothing_length(smoothlen, /devac)
endif else begin
    ok = oplot_smoothing_length(smoothlen, /plog)
endelse

;xyouts, 0.35, 0.85, fload_fid(1), size=1.7, color= 0, /normal, charthick=2.5
;xyouts, 0.35, 0.8, 'Standard', size=1.7, color= 0, /normal



;  done
; -------
device, /close



end






;============================================================================





pro plot_sd_profile_all, frun, snapnum, $
			smoothlen=smoothlen, $
			filename=filename, $
			x_is_devac=x_is_devac, $
			x_is_linear=x_is_linear

if not keyword_set(frun) then begin
   print, "  "
   print, "plot_sd_profile, frun, snapnum, smoothlen=smoothlen, filename=filename, $"
   print, "                /x_is_devac, /x_is_linear  (default is log)"
   print, "  "
   return
endif

if not keyword_set(smoothlen) then smoothlen= 0.1
if not keyword_set(filename) then filename='devac.eps'
if not keyword_set(snapnum) then snapnum=0


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 40

yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"


; default is log
;
xaxistitle="!6Log R (kpc)"
xmax = 2.0
xmin = -1.5
ymax = 6.0
ymin = -1.0
x_is_log= 1

if keyword_set(x_is_devac) then begin
	xaxistitle="!6[R (kpc)]!E1/4!N"
	xmax = 2.4
	;xmax = 4.2
	xmin = 0.5
	ymax = 6.0
	ymin = -1.0
	x_is_log= 0
endif

if keyword_set(x_is_linear) then begin
        xaxistitle="!6 R (kpc)"
        xmax = 20.0
        xmin = 0.0
        ymax = 4.2
        ymin = -1.0
	x_is_log= 0
endif




; ------------------
; Plot this up
; ------------------
!p.position= [0.2, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, /nodata, $
        ;/ylog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle




; -----------------------------------------------


	process_and_plot_sd, frun, snapnum, "allstars", xmin, xmax, bins, linecolor=0, $
						x_is_devac=x_is_devac, $
						x_is_log=x_is_log



	; gas profile
	; --------------
	process_and_plot_sd, "loaded", 0, "gas", xmin, xmax, bins, linecolor=50, $
                                                x_is_devac=x_is_devac, $
                                                x_is_log=x_is_log


	; disk profile
	; --------------
	process_and_plot_sd, "loaded", 0, "disk", xmin, xmax, bins, linecolor=100, $
                                                x_is_devac=x_is_devac, $
                                                x_is_log=x_is_log


	; bulge profile
	; --------------
	process_and_plot_sd, "loaded", 0, "bulge", xmin, xmax, bins, linecolor=200, $
                                                x_is_devac=x_is_devac, $
                                                x_is_log=x_is_log


	; new stars profile
	; -------------------
	process_and_plot_sd, "loaded", 0, "newstars", xmin, xmax, bins, linecolor=150, $
                                                x_is_devac=x_is_devac, $
                                                x_is_log=x_is_log


	do_special= 0
	;do_special= 1
	if do_special eq 1 then begin
		startid= 1600001
		numpart= 170001
                x=fload_1gal_allstars_xyz('x', startid, numpart)
                y=fload_1gal_allstars_xyz('y', startid, numpart)
                z=fload_1gal_allstars_xyz('z', startid, numpart)
                m=fload_1gal_allstars_mass(startid, numpart)

                process_and_plot_sd, "special", 0, "special", m, xmin, xmax, bins, linecolor=50, $
                                                x_is_devac=x_is_devac, $
                                                x_is_log=x_is_log, $
						special_x= x, $
						special_y= y, $
						special_z= z, $
						special_m= m

	endif




        x0= 0.65
        y0= 0.9
        xyouts, x0, y0, '!6Total', /normal, charthick=3.0, size= 1.2, color= 0
        if fload_npart(0) GT 0 then xyouts, x0, y0-0.04, 'Gas', /normal, charthick= 3.0, size= 1.2, color= 50
        if fload_npart(2) GT 0 then xyouts, x0, y0-0.08, 'Old Stars', /normal, charthick= 3.0, size= 1.2, color= 100
        if fload_npart(3) GT 0 then xyouts, x0, y0-0.12, 'Bulge', /normal, charthick= 3.0, size= 1.2, color= 200
        if fload_npart(4) GT 0 then xyouts, x0, y0-0.16, 'New Stars', /normal, charthick= 3.0, size= 1.2, color= 150
        if do_special GT 0 then xyouts, x0, y0-0.20, 'Sattelite', /normal, charthick= 3.0, size= 1.2, color= 50



; extras
; --------

if keyword_set(x_is_devac) then begin
    ok = oplot_smoothing_length(smoothlen, /devac)
endif else begin
    ok = oplot_smoothing_length(smoothlen, /plog)
endelse

;xyouts, 0.35, 0.85, fload_fid(1), size=1.7, color= 0, /normal, charthick=2.5
;xyouts, 0.35, 0.8, 'Standard', size=1.7, color= 0, /normal



;  done
; -------
device, /close



end








;===============================================================================








; --------------------------------
;  Compute Merger Remnant (Avg.)
; --------------------------------
pro plot_devacouleur_profile_avg, frun, snapnum, smoothlen=smoothlen, filename=filename

if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_devacouleur_profile_avg"
   print, "  "
   return
endif

if not keyword_set(smoothlen) then smoothlen= 0.1
if not keyword_set(snapnum) then snapnum= -99
if not keyword_set(filename) then filename='devac.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 100

xmax = 2.4
;xmax = 4.2
xmin = 0.3
ymax = 6.0
ymin = -1.0


xaxistitle= "!6[R (kpc)]!E1/4!N"
yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"


; ------------------
; Plot this up
; ------------------
!p.position= [0.2, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
	;/ylog, $
	color= 0, $
	xrange=[xmin,xmax], $
	yrange=[ymin,ymax], $
	xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata


; -----------------------------------------------


	;process_and_plot_sd, frun, snapnum, "allstars", xmin, xmax, bins, linecolor=0, /includeerr
	process_and_plot_sd, frun, snapnum, "allstars", xmin, xmax, bins, linecolor=0, /x_is_devac


	; old stellar disk profile
	; ------------------------
	process_and_plot_sd, "loaded", 0, "disk", xmin, xmax, bins, linecolor=50, /x_is_devac


	; old stellar bulge profile
	; -------------------------
	process_and_plot_sd, "loaded", 0, "bulge", xmin, xmax, bins, linecolor=100, /x_is_devac


	; new stars profile
	; -------------------
	process_and_plot_sd, "loaded", 0, "newstars", xmin, xmax, bins, linecolor=150, /x_is_devac




        ;if fload_npart(0) GT 0 then oplot, r_devac_gas, mass_sd_gas, psym=-3, linestyle=0, color= 50, thick= 2.5
        ;if fload_npart(2) GT 0 then oplot, r_devac_disk, mass_sd_disk, psym=-3, linestyle=0, color= 100, thick= 2.5
        ;if fload_npart(3) GT 0 then oplot, r_devac_bulge, mass_sd_bulge, psym=-3, linestyle=0, color= 200, thick= 2.5
        ;if fload_npart(4) GT 0 then oplot, r_devac_stars, mass_sd_stars, psym=-3, linestyle=0, color= 150, thick= 2.5

        x0= 0.65
        y0= 0.9
        ;xyouts, x0, y0, 'Total', /normal, charthick=3.0, size= 1.7, color= 0
        ;if fload_npart(0) GT 0 then xyouts, x0, y0-0.04, 'Gas', /normal, charthick= 3.0, size= 1.7, color= 50
        ;if fload_npart(2) GT 0 then xyouts, x0, y0-0.08, 'Disk', /normal, charthick= 3.0, size= 1.7, color= 100
        ;if fload_npart(3) GT 0 then xyouts, x0, y0-0.16, 'Bulge', /normal, charthick= 3.0, size= 1.7, color= 200
        ;if fload_npart(4) GT 0 then xyouts, x0, y0-0.12, 'New Stars', /normal, charthick= 3.0, size= 1.7, color= 150




; print extras
; -------------
x=[smoothlen,smoothlen]
x=x^(0.25)
y=[ymin,ymax]
oplot, x, y, linestyle=1, color= 0

;xyouts, 0.75, 0.90, fload_fid(1), size=1.5, color= 0, /normal, charthick=4.0
;xyouts, 0.75, 0.85, fload_timelbl(1,1), size=1.5, color= 0, /normal, charthick=4.0

xyouts, 0.60, 0.90, "Total", size=1.3, color= 0, /normal, charthick=4.0
xyouts, 0.60, 0.86, "Old stars (disk)", size=1.3, color= 50, /normal, charthick=4.0
xyouts, 0.60, 0.82, "Old stars (bulge)", size=1.3, color= 100, /normal, charthick=4.0
xyouts, 0.60, 0.78, "New stars", size=1.3, color= 150, /normal, charthick=4.0


; done
; -----
device, /close


end





;================================================================================
;================================================================================








; -------------------------------
;
;  Overplot devac profile for one
;  galaxy (and one snap).
;
; -------------------------------
pro do_one_galaxy_devac, frun, snapnum, linecolor=linecolor, msg=msg, $
					fitdevac=fitdevac, $
					fitsersic=fitsersic, $
					includeerr=includeerr, $
					frommanyproj=frommanyproj

; make sure these are same as 
; calling routine
;
bins = 40

xmax = 2.4
;xmax = 4.2
xmin = 0.4
;xmin = 0.5
;ymax = 6.0
;ymin = -1.0

	;
	; total (all star) profile
	; --------------------------

	;process_and_plot_sd, x, y, z, m, xmin, xmax, bins, linecolor=linecolor, /calcfit
	process_and_plot_sd, frun, snapnum, "allstars", xmin, xmax, bins, linecolor=linecolor, $
					fitdevac=fitdevac, /x_is_devac, $
					fitsersic=fitsersic, $
					includeerr=includeerr, $
					frommanyproj=frommanyproj, $
					h= 0.7

	if keyword_set(msg) then begin
	  ;xyouts, 0.75, 0.90, msg, size=1.5, color=linecolor, /normal, charthick=4.0
	endif

end





;================================================================================





; --------------------------------
;  Compute Merger Remnant (Avg.)
;
; --------------------------------
pro plot_devacouleur_profile_avg_comparison, junk, smoothlen=smoothlen, filename=filename, snapnum=snapnum

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: plot_devacouleur_profile_avg_comparison"
   print, "  "
   return
endif

if not keyword_set(smoothlen) then smoothlen= 0.1
if not keyword_set(filename) then filename= 'devac.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=1
;setup_plot_stuff, 'ps', filename=filename, colortable=2
;setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 100

xmax = 2.4
;xmax = 4.2
xmin = 0.4
ymax = 6.0
ymin = -1.0


xaxistitle= "!6[R (kpc)]!E1/4!N"
yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"


if not keyword_set(snapnum) then snapnum= 25

; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
	;/ylog, $
	color= 0, $
	xrange=[xmin,xmax], $
	yrange=[ymin,ymax], $
	xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata


; -----------------------------------------------

	;
	; Galaxy 1
	;
	; total (all star) profile
	; --------------------------
	;frun= "pool/vc3vc3"
	;frun= "pool/vc3bvc3b"
	;frun= "/raid4/tcox/vc3rem_vc3"
	;frun= "/raid4/tcox/vc3vc3"
	;frun= "/raid4/tcox/vc3vc3b"
	;do_one_galaxy_devac, frun, snapnum, linecolor=20, msg=fload_getid(frun)

        ;
        ; Galaxy 2
        ; --------------------------

        ;frun= "pool/vc3vc3_wBH"
        ;frun= "pool/vc3bvc3b_wBH"
	;frun= "/raid4/tcox/vc3rem_vc3"
	;frun= "/raid4/tcox/vc3rem_vc3rem"
	;frun= "/raid4/tcox/vc3vc3c"

   do_remergers = 0
   if do_remergers eq 1 then begin
	do_one_galaxy_devac, "/raid4/tcox/vcs/vc3vc3", 30, linecolor= 150, /fitsersic, /includeerr
	do_one_galaxy_devac, "/raid4/tcox/vc3rem_vc3", 40, linecolor= 100, /fitsersic, /includeerr
	do_one_galaxy_devac, "/raid4/tcox/vc3rem_vc3rem", 40, linecolor= 50, /fitsersic, /includeerr
   endif


   do_mhmergers = 0
   if do_mhmergers eq 1 then begin
	;;do_one_galaxy_devac, "/raid4/tcox/vc3vc3h_2", 107, linecolor= 150, /fitsersic, /includeerr
	;;do_one_galaxy_devac, "/raid4/tcox/vc3vc3h_mh_no", 107, linecolor= 50, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/vc3vc3h_2", 107, linecolor= 0, /oneproj, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/vc3vc3h_no", 107, linecolor= 100, /oneproj, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/vc3vc3h_mh_no", 107, linecolor= 150, /oneproj, /includeerr
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3h_2", 107, linecolor= 0, /includeerr
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3h_no", 107, linecolor= 100, /includeerr
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3h_mh_no", 107, linecolor= 150, /includeerr
	xyouts, 0.55, 0.85, 'std', size=1.5, color=0, /normal, charthick=4.0
	xyouts, 0.55, 0.80, 'no', size=1.5, color=100, /normal, charthick=4.0
	xyouts, 0.55, 0.75, 'MH, no', size=1.5, color=150, /normal, charthick=4.0
   endif



   do_all_vcs= 0
   if do_all_vcs eq 1 then begin
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3b", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3c", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3d", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3e", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3f", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3g", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3h", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3i", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3j", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3k", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3l", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3m", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3n", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3o", snapnum, linecolor=20
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3p", snapnum, linecolor=20
   ;endif

   ;do_all_vcs= 1
   ;if do_all_vcs eq 1 then begin
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3b", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3c", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3d", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3e", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3f", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3g", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3h", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3i", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3j", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3k", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3l", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3m", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3n", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3o", snapnum, linecolor=200
        do_one_galaxy_devac, "/raid4/tcox/collisionless/cvc3vc3p", snapnum, linecolor=200
   endif



   do_vcs_withandwithout= 0
   if do_vcs_withandwithout eq 1 then begin
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3f", snapnum, linecolor= 150
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3h", snapnum, linecolor= 150
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3j", snapnum, linecolor= 150

	do_one_galaxy_devac, "/raid4/tcox/vc3vc3f_no", snapnum, linecolor= 50
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3h_no", snapnum, linecolor= 50
	do_one_galaxy_devac, "/raid4/tcox/vc3vc3j_no", snapnum, linecolor= 50
   endif
	
   do_large_vcs= 0
   if do_large_vcs eq 1 then begin
	;do_one_galaxy_devac, "/raid4/tcox/vc4vc4a", snapnum, linecolor= 150
	;do_one_galaxy_devac, "/raid4/tcox/vc5vc5a", snapnum, linecolor= 150
	;do_one_galaxy_devac, "/raid4/tcox/vc6vc6a", snapnum, linecolor= 150

	;do_one_galaxy_devac, "/raid4/tcox/vc4vc4_no", snapnum, linecolor= 50
	;do_one_galaxy_devac, "/raid4/tcox/vc5vc5_no", snapnum, linecolor= 50
	;do_one_galaxy_devac, "/raid4/tcox/vc6vc6_no", snapnum, linecolor= 50

	do_one_galaxy_devac, "/raid4/tcox/ds/d6k", 30, linecolor= 180, /oneproj, /fitsersic
        do_one_galaxy_devac, "/raid4/tcox/ds/d5k", 30, linecolor= 150, /oneproj, /fitsersic
        do_one_galaxy_devac, "/raid4/tcox/ds/d4k", 30, linecolor= 120, /oneproj, /fitsersic
        ;do_one_galaxy_devac, "/raid4/tcox/vc3vc3k", 30, linecolor= 90, /oneproj, /fitsersic
        ;do_one_galaxy_devac, "/raid4/tcox/ds/d3k", 30, linecolor= 90, /oneproj, /fitsersic
        do_one_galaxy_devac, "/raid4/tcox/ds/d2k", 30, linecolor= 60, /oneproj, /fitsersic
        do_one_galaxy_devac, "/raid4/tcox/ds/d1k", 30, linecolor= 30, /oneproj, /fitsersic
        do_one_galaxy_devac, "/raid4/tcox/ds/d0k", 50, linecolor= 10, /oneproj, /fitsersic

   endif

   do_vc3_time= 1
   if do_vc3_time eq 1 then begin
	;frun= "/raid4/tcox/vc3vc3h"
	;frun= "/raid4/tcox/ds/d6h"
	;frun= "/raid4/tcox/ds/d1h"
	;frun= "/raid4/tcox/ds/d0h"
	;frun= "/raid4/tcox/vc3vc3h_no"
	frun="data/z3/b5e"
	;snapnum= 15
	snapnum= 60
	initlc = 10
	;do_one_galaxy_devac, frun, snapnum, linecolor= initlc
	;do_one_galaxy_devac, frun, snapnum+1, linecolor= initlc+10
	;do_one_galaxy_devac, frun, snapnum+2, linecolor= initlc+20
	;do_one_galaxy_devac, frun, snapnum+3, linecolor= initlc+30
	;do_one_galaxy_devac, frun, snapnum+4, linecolor= initlc+40
	;do_one_galaxy_devac, frun, snapnum+5, linecolor= initlc+50
	;do_one_galaxy_devac, frun, snapnum+6, linecolor= initlc+60
	;do_one_galaxy_devac, frun, snapnum+7, linecolor= initlc+70
	;do_one_galaxy_devac, frun, snapnum+8, linecolor= initlc+80
	;do_one_galaxy_devac, frun, snapnum+9, linecolor= initlc+90
	;do_one_galaxy_devac, frun, snapnum+10, linecolor= initlc+100
	;do_one_galaxy_devac, frun, snapnum+11, linecolor= initlc+110
	;do_one_galaxy_devac, frun, snapnum+12, linecolor= initlc+120
	;do_one_galaxy_devac, frun, snapnum+13, linecolor= initlc+130
	;do_one_galaxy_devac, frun, snapnum+14, linecolor= initlc+140
	do_one_galaxy_devac, frun, snapnum, linecolor= initlc
	do_one_galaxy_devac, frun, snapnum+10, linecolor= initlc+20
	do_one_galaxy_devac, frun, snapnum+20, linecolor= initlc+40
	do_one_galaxy_devac, frun, snapnum+30, linecolor= initlc+60
	do_one_galaxy_devac, frun, snapnum+40, linecolor= initlc+80
	do_one_galaxy_devac, frun, snapnum+50, linecolor= initlc+100
	do_one_galaxy_devac, frun, snapnum+60, linecolor= initlc+120
	do_one_galaxy_devac, frun, snapnum+70, linecolor= initlc+140

	;xyouts, 0.55, 0.90, 'no black hole', size=1.5, color=0, /normal, charthick=4.0
	xyouts, 0.55, 0.85, '!6time course', size=1.5, color=0, /normal, charthick=4.0
	xyouts, 0.55, 0.80, 'dark -> ', size=1.5, color=20, /normal, charthick=4.0
	xyouts, 0.75, 0.80, ' light', size=1.5, color=160, /normal, charthick=4.0
   endif

   do_vc3_gf= 0
   if do_vc3_gf eq 1 then begin
	;do_one_galaxy_devac, "/raid4/tcox/vc3vc3", 19, linecolor= 50
	;do_one_galaxy_devac, "/raid4/tcox/vc3vc3h", 19, linecolor= 100
	;do_one_galaxy_devac, "/raid4/tcox/As/A3", 140, linecolor= 150
	;
	do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3u_h", 30, linecolor= 180
	do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3v_h", 30, linecolor= 150
	do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3w_h", 30, linecolor= 120
	do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3x2_h", 30, linecolor= 90
	do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3y_h", 30, linecolor= 60
	do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3z_h", 30, linecolor= 30
	do_one_galaxy_devac, "/raid4/tcox/cvc3vc3h", 30, linecolor= 10
	;

	; ------
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3u_e", 30, linecolor= 180, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3v_e", 30, linecolor= 150, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3w_e", 30, linecolor= 120, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3x2_e", 30, linecolor= 90, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3y_e", 30, linecolor= 60, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3z_e", 30, linecolor= 30, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/cvc3vc3e", 30, linecolor= 10, /fitsersic, /includeerr
	; ------
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3u_h", 30, linecolor= 180, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3v_h", 30, linecolor= 150, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3w_h", 30, linecolor= 120, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3x2_h", 30, linecolor= 90, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3y_h", 30, linecolor= 60, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3z_h", 30, linecolor= 30, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/cvc3vc3h", 30, linecolor= 10, /fitsersic, /includeerr
	; ------
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3u_k", 30, linecolor= 180, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3v_k", 30, linecolor= 150, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3w_k", 30, linecolor= 120, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3x2_k", 30, linecolor= 90, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3y_k", 30, linecolor= 60, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/gfs/vc3vc3z_k", 30, linecolor= 30, /fitsersic, /includeerr
	;do_one_galaxy_devac, "/raid4/tcox/cvc3vc3k", 30, linecolor= 10, /fitsersic, /includeerr
	; ------
   endif


   do_vc3_hdmdiff= 0
   if do_vc3_hdmdiff eq 1 then begin
        do_one_galaxy_devac, "/raid4/tcox/vc3vc3", snapnum, linecolor= 50
        do_one_galaxy_devac, "/raid4/tcox/vc3vc3a", snapnum, linecolor= 150
	xyouts, 0.75, 0.85, 'vc3vc3', size=1.5, color=50, /normal, charthick=4.0
	xyouts, 0.75, 0.80, 'vc3vc3a', size=1.5, color=150, /normal, charthick=4.0
   endif


   do_vc3_cless= 0
   if do_vc3_cless eq 1 then begin
        do_one_galaxy_devac, "/raid4/tcox/cvc3vc3b", snapnum, linecolor= 50
        do_one_galaxy_devac, "/raid4/tcox/cvc3vc3b_2", snapnum, linecolor= 150
	xyouts, 0.55, 0.85, 'cvc3vc3b, h=0.2', size=1.5, color=50, /normal, charthick=4.0
	xyouts, 0.55, 0.80, 'cvc3vc3b_2, h=0.4', size=1.5, color=150, /normal, charthick=4.0
   endif


   do_vc3_clessvnot= 0
   if do_vc3_clessvnot eq 1 then begin
        do_one_galaxy_devac, "/raid4/tcox/vc3vc3e", snapnum, linecolor= 50
        do_one_galaxy_devac, "/raid4/tcox/cvc3vc3e", snapnum, linecolor= 150
        xyouts, 0.55, 0.85, 'vc3vc3e', size=1.5, color=50, /normal, charthick=4.0
        xyouts, 0.55, 0.80, 'cvc3vc3e', size=1.5, color=150, /normal, charthick=4.0
   endif


   do_sbwind= 0
   if do_sbwind eq 1 then begin
        do_one_galaxy_devac, "/raid4/tcox/vc3vc3e_no", 107, linecolor= 0
        do_one_galaxy_devac, "/raid4/tcox/vc3vc3e_2", 107, linecolor= 100
        do_one_galaxy_devac, "/raid4/tcox/sbw/sb10BH", 107, linecolor= 50
        do_one_galaxy_devac, "/raid4/tcox/sbw/sb10", 107, linecolor= 150
        xyouts, 0.65, 0.85, 'std', size=1.5, color=0, /normal, charthick=4.0
        xyouts, 0.65, 0.80, 'stdBH', size=1.5, color=100, /normal, charthick=4.0
        xyouts, 0.65, 0.75, 'sb10BH', size=1.5, color=50, /normal, charthick=4.0
        xyouts, 0.65, 0.70, 'sb10', size=1.5, color=150, /normal, charthick=4.0
   endif


   do_fb= 0
   if do_fb eq 1 then begin
        do_one_galaxy_devac, "/raid4/tcox/ds/d5e2", 30, linecolor= 50
        do_one_galaxy_devac, "/raid4/tcox/ds/d5e2_q", 30, linecolor= 150
        xyouts, 0.65, 0.85, 'high feedback', size=1.5, color=50, /normal, charthick=4.0
        xyouts, 0.65, 0.80, 'low feedback', size=1.5, color=150, /normal, charthick=4.0
   endif



   do_inta= 0
   if do_inta eq 1 then begin
        do_one_galaxy_devac, "/raid4/tcox/bs/b3e_i2", 30, linecolor= 0
        do_one_galaxy_devac, "/raid4/tcox/bs/b3e_i1", 30, linecolor= 50
        do_one_galaxy_devac, "/raid4/tcox/bs/b3e", 30, linecolor= 150
        do_one_galaxy_devac, "/raid4/tcox/bs/b3e_i3", 30, linecolor= 100
	xyouts, 0.60, 0.92, 'ErrTolIntAccuracy', /normal, charthick=3, size=1.33, color=0
        xyouts, 0.65, 0.85, '1.0', size=1.5, color=0, /normal, charthick=4.0
        xyouts, 0.65, 0.80, '0.1', size=1.5, color=50, /normal, charthick=4.0
        xyouts, 0.65, 0.75, '0.025', size=1.5, color=150, /normal, charthick=4.0
        xyouts, 0.65, 0.70, '0.0025', size=1.5, color=100, /normal, charthick=4.0
   endif


   do_inta= 0
   if do_inta eq 1 then begin
        ;do_one_galaxy_devac, "/raid4/tcox/bs/b3e_i3", 12, linecolor= 0
        ;do_one_galaxy_devac, "/raid4/tcox/bs/b3e_i3", 15, linecolor= 50
        ;do_one_galaxy_devac, "/raid4/tcox/bs/b3e_i3", 20, linecolor= 150
        ;do_one_galaxy_devac, "/raid4/tcox/bs/b3e_i3", 30, linecolor= 100
        ;do_one_galaxy_devac, "/raid4/tcox/bs/b3e", 12, linecolor= 0
        ;do_one_galaxy_devac, "/raid4/tcox/bs/b3e", 15, linecolor= 50
        ;do_one_galaxy_devac, "/raid4/tcox/bs/b3e", 20, linecolor= 150
        ;do_one_galaxy_devac, "/raid4/tcox/bs/b3e", 30, linecolor= 100
        ;xyouts, 0.65, 0.85, '12', size=1.5, color=0, /normal, charthick=4.0
        ;xyouts, 0.65, 0.80, '15', size=1.5, color=50, /normal, charthick=4.0
        ;xyouts, 0.65, 0.75, '20', size=1.5, color=150, /normal, charthick=4.0
        ;xyouts, 0.65, 0.70, '30', size=1.5, color=100, /normal, charthick=4.0
	;xyouts, 0.40, 0.92, 'ErrTolIntAccuracy= 0.025', /normal, charthick=3, size=1.33, color=0
	;xyouts, 0.40, 0.92, 'ErrTolIntAccuracy= 0.0025', /normal, charthick=3, size=1.33, color=0
        ;xyouts, 0.60, 0.85, 'T!Dmerg!N + 0.2 Gyr', size=1.5, color=0, /normal, charthick=4.0
        ;xyouts, 0.60, 0.80, 'T!Dmerg!N + 0.6 Gyr', size=1.5, color=50, /normal, charthick=4.0
        ;xyouts, 0.60, 0.75, 'T!Dmerg!N + 1.3 Gyr', size=1.5, color=150, /normal, charthick=4.0
        ;xyouts, 0.60, 0.70, 'T!Dmerg!N + 2.8 Gyr', size=1.5, color=100, /normal, charthick=4.0
        do_one_galaxy_devac, "/raid4/tcox/es/e3e", 13, linecolor= 0
        do_one_galaxy_devac, "/raid4/tcox/es/e3e", 16, linecolor= 50
        do_one_galaxy_devac, "/raid4/tcox/es/e3e", 20, linecolor= 150
        do_one_galaxy_devac, "/raid4/tcox/es/e3e", 30, linecolor= 100
        xyouts, 0.65, 0.85, '13', size=1.5, color=0, /normal, charthick=4.0
        xyouts, 0.65, 0.80, '16', size=1.5, color=50, /normal, charthick=4.0
        xyouts, 0.65, 0.75, '20', size=1.5, color=150, /normal, charthick=4.0
        xyouts, 0.65, 0.70, '30', size=1.5, color=100, /normal, charthick=4.0
   endif


   do_Sbc= 0
   if do_Sbc eq 1 then begin
        do_one_galaxy_devac, "/raid12/tcox/Sbc/Sbc201a", 20, linecolor= 0
        do_one_galaxy_devac, "/raid12/tcox/Sbc/Sbc201a-u4", 40, linecolor= 50
        xyouts, 0.65, 0.85, 'Gadget-2', size=1.5, color=0, /normal, charthick=4.0
        xyouts, 0.65, 0.80, 'Gadget-1', size=1.5, color=50, /normal, charthick=4.0
   endif



   ;xyouts, 0.75, 0.90, fload_timelbl(1,2), size=1.5, color=0, /normal, charthick=4.0


; print extras
; -------------
if not keyword_set(smoothlen) then smoothlen= 0.1/0.7 else smoothlen= smoothlen/0.7
x=[smoothlen,smoothlen]
x=x^(0.25)
y=[ymin,ymax]
oplot, x, y, linestyle=1, color= 0

;smoothlen= 0.4
;x=[smoothlen,smoothlen]
;x=x^(0.25)
;y=[ymin,ymax]
;oplot, x, y, linestyle=1, color= 0




; done
; -----
device, /close


end







;=================================================================================
;=================================================================================
;
;
;
;
;=================================================================================
;=================================================================================
;
;
;
;
;=================================================================================
;=================================================================================



; --------------------------------
;  6 panel devac profile figure
;
; --------------------------------
pro plot_devacouleur_6, junk, smoothlen=smoothlen, filename=filename, snapnum=snapnum

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: plot_devacouleur_6, junk"
   print, "  "
   return
endif

if not keyword_set(smoothlen) then smoothlen= 0.1
if not keyword_set(filename) then filename= 'devac.eps'


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1, newxsize=30, newysize=20
;setup_plot_stuff, 'ps', filename=filename, colortable=2, newxsize=30, newysize=20
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=30, newysize=20



; -----------------------------
; set up constants
; -----------------------------

bins = 50

xmax = 2.4
;xmax = 4.2
xmin = 0.5
ymax = 5.3
ymin = -1.0


xaxistitle= "!6[R (kpc)]!E1/4!N"
yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"


frun1= "/raid4/tcox/gfs/vc3vc3u_k"   & msg1='1.0'
frun2= "/raid4/tcox/gfs/vc3vc3v_k"   & msg2='0.8'
frun3= "/raid4/tcox/gfs/vc3vc3w_k"   & msg3='0.6'
frun4= "/raid4/tcox/gfs/vc3vc3x2_k"   & msg4='0.4'
frun5= "/raid4/tcox/gfs/vc3vc3y_k"   & msg5='0.2'
frun6= "/raid4/tcox/gfs/vc3vc3z_k"   & msg6='0.05'
cfrun= "/raid4/tcox/cvc3vc3k"

if not keyword_set(snapnum) then snapnum= 30

x0= 0.08
xs= 0.30   ; assumes 3 panels
x1= x0+xs
x2= x0+xs+xs
x3= x0+xs+xs+xs
    
y0= 0.08
ys= 0.45   ; assumes 2 panels
y1= y0+ys
y2= y0+ys+ys



; smoothing length
; -----------------
smoothlen= 0.2
x=[smoothlen,smoothlen]
xsmooth=x^(0.25)
ysmooth=[ymin,ymax]


; ------------------
;    1
; ------------------
!p.position= [x0, y1, x1, y2]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	xtickformat='(a1)', ytitle=yaxistitle, /noerase

;do_one_galaxy_devac, frun1, snapnum, linecolor=150, msg=fload_getid(frun1)
do_one_galaxy_devac, frun1, snapnum, linecolor=150, /oneproj
do_one_galaxy_devac, cfrun, snapnum, linecolor= 0, /oneproj

oplot, xsmooth, ysmooth, linestyle=1, color= 0

oplot, [1.2,1.4], [4.6,4.6], psym=-3, color= 0, thick=6.0, linestyle= 1
xyouts, 1.50, 4.5, 'collisionless', /data, size= 1.5, color= 0, charthick=5.0

oplot, [1.2,1.4], [4.1,4.1], psym=-3, color= 150, thick=8.0, linestyle= 0
xyouts, 1.50, 4.0, 'f= ', /data, size= 1.5, color= linecolor, charthick=5.0
xyouts, 1.80, 4.0, msg1, /data, size= 1.5, color= linecolor, charthick=5.0


; ------------------
;    2
; ------------------
!p.position= [x1, y1, x2, y2]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase

do_one_galaxy_devac, frun2, snapnum, linecolor=150, /oneproj
;do_one_galaxy_devac, frun2, snapnum, linecolor= 150, /fitsersic, /includeerr
do_one_galaxy_devac, cfrun, snapnum, linecolor= 0, /oneproj

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 1.80, 4.0, msg2, /data, size= 1.5, color= linecolor, charthick=5.0


; ------------------
;    3
; ------------------
!p.position= [x2, y1, x3, y2]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase

do_one_galaxy_devac, frun3, snapnum, linecolor=150, /oneproj
do_one_galaxy_devac, cfrun, snapnum, linecolor= 0, /oneproj

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 1.80, 4.0, msg3, /data, size= 1.5, color= linecolor, charthick=5.0


; ------------------
;    4
; ------------------
!p.position= [x0, y0, x1, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /noerase

do_one_galaxy_devac, frun4, snapnum, linecolor=150, /oneproj
do_one_galaxy_devac, cfrun, snapnum, linecolor= 0, /oneproj

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 1.80, 4.0, msg4, /data, size= 1.5, color= linecolor, charthick=5.0


; ------------------
;    5
; ------------------
!p.position= [x1, y0, x2, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytickformat='(a1)', /noerase

do_one_galaxy_devac, frun5, snapnum, linecolor=150, /oneproj
do_one_galaxy_devac, cfrun, snapnum, linecolor= 0, /oneproj

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 1.80, 4.0, msg5, /data, size= 1.5, color= linecolor, charthick=5.0


; ------------------
;    6
; ------------------
!p.position= [x2, y0, x3, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytickformat='(a1)', /noerase

do_one_galaxy_devac, frun6, snapnum, linecolor=150, /oneproj
do_one_galaxy_devac, cfrun, snapnum, linecolor= 0, /oneproj

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 1.80, 4.0, msg6, /data, size= 1.5, color= linecolor, charthick=5.0





; done
; -----
device, /close


end





;================================================================================
;================================================================================
;================================================================================





; --------------------------------
;  Profiles
;
; --------------------------------
pro profs, junk, smoothlen=smoothlen, filename=filename

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: profs"
   print, "  "
   return
endif

if not keyword_set(smoothlen) then smoothlen= 0.1
if not keyword_set(filename) then filename= 'profiles.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 40

xmax = 85.0
xmin = 0.0
ymax = 2.5
ymin = -1.0


xaxistitle= "!6R (kpc)"
yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"
;yaxistitle= "!4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"



; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
	;/ylog, $
	color= 0, $
	xrange=[xmin,xmax], $
	yrange=[ymin,ymax], $
	xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata


; -----------------------------------------------



process_and_plot_sd, "data1/Sbc/Sbc_10x_wBH", 0, "gas", xmin, xmax, bins, linecolor=0, thislinest= 2
xyouts, 0.2, 0.48, "Sbc_10x_wBH", /normal, charthick=3.0, size=1.2, color= 0

process_and_plot_sd, "data1/Sbc/Sbcm1", 0, "gas", xmin, xmax, bins, linecolor=50
xyouts, 0.2, 0.44, "Sbcm1", /normal, charthick=3.0, size=1.2, color= 50
process_and_plot_sd, "data1/Sbc/Sbcm2", 0, "gas", xmin, xmax, bins, linecolor=100
xyouts, 0.2, 0.40, "Sbcm2", /normal, charthick=3.0, size=1.2, color= 100
process_and_plot_sd, "data1/Sbc/Sbcm3", 0, "gas", xmin, xmax, bins, linecolor=150
xyouts, 0.2, 0.36, "Sbcm3", /normal, charthick=3.0, size=1.2, color= 150

process_and_plot_sd, "data1/Sbc/Sbccut", 0, "gas", xmin, xmax, bins, linecolor=200
xyouts, 0.2, 0.32, "Sbccut", /normal, charthick=3.0, size=1.2, color= 200



; -----------------------------------------------


; print extras
; -------------
if not keyword_set(smoothlen) then smoothlen= 0.1/0.7 else smoothlen= smoothlen/0.7
x=[smoothlen,smoothlen]
;x=x^(0.25)
y=[ymin,ymax]
oplot, x, y, linestyle=1, color= 0

;smoothlen= 0.4
;x=[smoothlen,smoothlen]
;x=x^(0.25)
;y=[ymin,ymax]
;oplot, x, y, linestyle=1, color= 0




; done
; -----
device, /close


end


