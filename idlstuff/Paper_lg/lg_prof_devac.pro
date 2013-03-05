
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


pro lg_prof_dev, junk

frun='/raid4/tcox/localgroup/v5'
snapnum=105
smoothlen=0.1
filename=frun+'/profdev.eps'
x_is_devac=1


if not keyword_set(frun) then begin
   print, "  "
   print, "lg_prof_dev, junk"
   print, "  "
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
xmax = 2.0
xmin = -1.5
ymax = 4.3
ymin = -1.0

if keyword_set(x_is_devac) then begin
	xaxistitle="!6[R (kpc)]!E1/4!N"
	xmax = 2.4
	;xmax = 4.2
	xmin = 0.5
	ymax = 4.3
	ymin = -1.0
	x_is_log= 0
endif

if keyword_set(x_is_linear) then begin
        xaxistitle="!6 R (kpc)"
        xmax = 20.0
        xmin = 0.0
        ymax = 6.0
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

   ok=fload_snapshot_bh(frun,snapnum)

; -----------------------------------------------


	; total (all star) profile
	; --------------------------
	x=fload_allstars_xyz('x')
	y=fload_allstars_xyz('y')
	z=fload_allstars_xyz('z')
        m=fload_allstars_mass(1)


	;process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=0, /includeerr
	process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=0, $
						x_is_devac=x_is_devac, $
						x_is_log=x_is_log   ;, /includeerr





	; new stars profile
	; -------------------
	if fload_npart(4) GT 0 then begin
		x=fload_newstars_xyz('x')
		y=fload_newstars_xyz('y')
		z=fload_newstars_xyz('z')
		m=fload_newstars_mass(1)

		process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=150, $
                                                x_is_devac=x_is_devac, $
                                                x_is_log=x_is_log  ;, /includeerr
	endif


        x0= 0.42
        y0= 0.40
        ;xyouts, x0, y0, '!6Total', /normal, charthick=3.0, size= 1.7, color= 0
        ;if fload_npart(0) GT 0 then xyouts, x0, y0-0.04, 'Gas', /normal, charthick= 3.0, size= 1.7, color= 50
        ;if fload_npart(2) GT 0 then xyouts, x0, y0-0.08, 'Disk', /normal, charthick= 3.0, size= 1.7, color= 100
        ;if fload_npart(3) GT 0 then xyouts, x0, y0-0.16, 'Bulge', /normal, charthick= 3.0, size= 1.7, color= 200
        ;if fload_npart(4) GT 0 then xyouts, x0, y0-0.05, 'New Stars', /normal, charthick= 3.0, size= 1.7, color= 0
	oplot, [0.7,0.9], [0.4,0.4], thick=4.0, color= 0, linestyle= 0
        xyouts, x0, y0-0.05, 'All Stars', /normal, charthick= 3.0, size= 1.7, color= 0
        ;if fload_npart(4) GT 0 then xyouts, x0, y0-0.10, 'All Stars', /normal, charthick= 3.0, size= 1.7, color= 150
	oplot, [0.7,0.9], [0.05,0.05], thick=4.0, color= 150, linestyle= 2
        xyouts, x0, y0-0.10, 'New Stars', /normal, charthick= 3.0, size= 1.7, color= 150



; extras
; --------

if keyword_set(devac) then begin
    ok = oplot_smoothing_length(smoothlen, /devac)
endif else begin
    ok = oplot_smoothing_length(smoothlen, /plog)
endelse

;xyouts, 0.35, 0.85, fload_fid(1), size=1.7, color= 0, /normal, charthick=2.5
;xyouts, 0.35, 0.8, 'Standard', size=1.7, color= 0, /normal



xyouts, 0.8, 0.9, '(b)', /normal, size=1.8, color= 0, charthick=3.0




;  done
; -------
device, /close



end








;===============================================================================








;----------------------------------------
;  Does some of the hard work
;----------------------------------------
pro process_and_plot_profile, x, y, z, m, $
			xmin, xmax, bins, $
			x_is_devac=x_is_devac, $
			x_is_log=x_is_log, $
			linecolor=linecolor, $
			includeerr=includeerr, $
			fitdevac=fitdevac, $
			fitsersic=fitsersic, $
			frommanyproj=frommanyproj


	; compute (appropriate) profile
	; ------------------------------
	r_s = fltarr(bins)
	mass_sd_avg = fltarr(bins)
	mass_sd_1sig = fltarr(bins)

	if not keyword_set(frommanyproj) then begin
		a= sqrt(x*x + y*y)
                c= m

		if keyword_set(x_is_devac) then a=a^(0.25)
		if keyword_set(x_is_log) then a=alog10(a)

		process_sd_profile, a, c, bins, xmax, xmin, r_s, mass_sd_avg, $
					sd_1sig=mass_sd_1sig, $
					x_is_devac=x_is_devac, $
					x_is_log=x_is_log
	endif else begin
		process_prof_frommanyprojections, x, y, z, m, $
					xmin, xmax, bins, $
					x_is_devac=x_is_devac, $
					x_is_log=x_is_log, $
					mass_sd_avg, mass_sd_1sig, r_s
	endelse



	; take out zeros
	; ----------------
	idx= where(mass_sd_avg gt 0)
	if idx(0) ne -1 then begin
		mass_sd_avg= mass_sd_avg(idx)
		mass_sd_1sig= mass_sd_1sig(idx)
		r_s= r_s(idx)
	endif


	; standard conversion to m_solar pc-2
	;-------------------------------------
	mass_sd_p1sig= alog10(mass_sd_avg+mass_sd_1sig) + 4.0
	mass_sd_m1sig= alog10(mass_sd_avg-mass_sd_1sig) + 4.0

	;mass_sd_p1sig= alog10(mass_sd_avg)+alog10(mass_sd_1sig) + 4.0
        ;mass_sd_m1sig= alog10(mass_sd_avg)-alog10(mass_sd_1sig) + 4.0
	mass_sd= alog10(mass_sd_avg) + 4.0

	
	; now printing
	; -------------
	if linecolor eq 200 then thislinest= 2 else thislinest= 0
	if linecolor eq 200 then thisthick= 1.0 else thisthick= 3.0
	if linecolor eq 200 then thispsym= 3 else thispsym= 3

	if linecolor eq 150 then thislinest= 2

	; plot the total
	oplot, r_s, mass_sd, psym=-thispsym, linestyle=thislinest, color= linecolor, thick=thisthick

	if keyword_set(includeerr) then begin
		oplot, r_s, mass_sd_p1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
		oplot, r_s, mass_sd_m1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
	endif


	; try to fit de Vac. to it
	; -------------------------
	if keyword_set(fitdevac) then begin

		if keyword_set(x_is_devac) then r_tofit= r_s^(4.0)
		if keyword_set(x_is_log) then r_tofit= 10^r_s

		sd_tofit= 10^(mass_sd)
		weight= mass_sd_1sig
		;weight= 1.0e4*abs(mass_sd_1sig/mass_sd)

		fitandoverplot_devacprofile, r_tofit, sd_tofit, $
					weight, /ylogaxis, $
					x_is_devac=x_is_devac, $
					x_is_log=x_is_log

	endif



        ; try to fit Sersic Profile to it
        ; --------------------------------
        if keyword_set(fitsersic) then begin

		r_tofit= r_s
                if keyword_set(x_is_devac) then r_tofit= r_s^(4.0)
		if keyword_set(x_is_log) then r_tofit= 10^r_s

                sd_tofit= 10^(mass_sd)
                weight= mass_sd_1sig
                ;weight= 1.0e4*abs(mass_sd_1sig/mass_sd)

		;idx=where(r_tofit gt 1.0)
		;r_tofit= r_tofit(idx)
		;sd_tofit= sd_tofit(idx)
		;weight= weight(idx)
		minfitradius= 0.5

                fitandoverplot_sersicprofile, r_tofit, sd_tofit, weight, $
						/ylogaxis, $
						minfitradius=minfitradius, $
						x_is_devac=x_is_devac, $
                                        	x_is_log=x_is_log


        endif




end






;================================================================================
;================================================================================




