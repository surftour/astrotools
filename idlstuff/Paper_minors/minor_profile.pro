
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------




pro plot_sd_profile_all, frun, snapnum, xmin=xmin, xmax=xmax, bins=bins


; -----------------------------------------------

   ok=fload_snapshot(frun,snapnum)

; -----------------------------------------------


	; total (all star) profile
	; --------------------------
	;x=fload_allstars_xyz('x')
	;y=fload_allstars_xyz('y')
	;z=fload_allstars_xyz('z')
        ;m=fload_allstars_mass(1)

	;process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=0, /includeerr
	;process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=0




	; gas profile
	; --------------
	;if fload_npart(0) GT 0 then begin
	;	x=fload_gas_xyz('x')
	;	y=fload_gas_xyz('y')
	;	z=fload_gas_xyz('z')
	;	m=fload_gas_mass(1)
	;	process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=50
	;endif


	; disk profile
	; --------------
        ;if fload_npart(2) GT 0 then begin
        ;        x=fload_disk_xyz('x')
        ;        y=fload_disk_xyz('y')
        ;        z=fload_disk_xyz('z')
        ;        m=fload_disk_mass(1)
	;	process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=100
        ;endif


	; bulge profile
	; --------------
        ;if fload_npart(3) GT 0 then begin
        ;        x=fload_bulge_xyz('x')
        ;        y=fload_bulge_xyz('y')
        ;        z=fload_bulge_xyz('z')
        ;        m=fload_bulge_mass(1)
	;	process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=200
        ;endif


	; new stars profile
	; -------------------
        if fload_npart(4) GT 0 then begin
                x=fload_newstars_xyz('x')
                y=fload_newstars_xyz('y')
                z=fload_newstars_xyz('z')
                m=fload_newstars_mass(1)
		process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=150
        endif



end






;===============================================================================





pro plot_sd_profile_newstars, frun, snapnum, xmin=xmin, xmax=xmax, bins=bins, thiscolor=thiscolor, $
				satnumpart= satnumpart


; -----------------------------------------------

   ok=fload_snapshot(frun,snapnum)

; -----------------------------------------------


        ; new stars profile
        ; -------------------
        x=fload_newstars_xyz('x')
        y=fload_newstars_xyz('y')
        z=fload_newstars_xyz('z')
        m=fload_newstars_mass(1)
        ;process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=1211
        process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=thiscolor
        ;process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=thiscolor, /use_dlogx


	;x=fload_1gal_newstars_xyz('x', 1, 240000, /useparentinfo)
	;y=fload_1gal_newstars_xyz('y', 1, 240000, /useparentinfo)
	;z=fload_1gal_newstars_xyz('z', 1, 240000, /useparentinfo)
	;m=fload_1gal_newstars_mass(1, 1, 240000, /useparentinfo)
        ;process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=1212


	;x=fload_1gal_newstars_xyz('x', 240001, satnumpart, /useparentinfo)
	;y=fload_1gal_newstars_xyz('y', 240001, satnumpart, /useparentinfo)
	;z=fload_1gal_newstars_xyz('z', 240001, satnumpart, /useparentinfo)
	;m=fload_1gal_newstars_mass(1, 240001, satnumpart, /useparentinfo)
        ;process_and_plot_profile, x, y, z, m, xmin, xmax, bins, linecolor=1213

end




;================================================================================





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
			frommanyproj=frommanyproj, $
			use_dlogx=use_dlogx


	; compute (appropriate) profile
	; ------------------------------
	r_s = fltarr(bins)
	mass_sd_avg = fltarr(bins)
	mass_sd_1sig = fltarr(bins)

	if not keyword_set(frommanyproj) then begin
		a= sqrt(x*x + y*y)
                c= m

		if keyword_set(x_is_devac) then a=a^(0.25)
		if keyword_set(use_dlogx) then begin
			xmax= alog10(xmax)
			xmin= alog10(xmin)
			x_is_log= 1
		endif
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


	; correct if we turned on log x spacing
	; --------------------------------------
	if keyword_set(use_dlogx) then begin
		xmax= 10^(xmax)
		xmin= 10^(xmin)
		x_is_log= 0
		r_s= 10^r_s
	endif

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
	thispsym= 3
	thislinest= 0
	thisthick= 3.0
	if linecolor eq 200 then thislinest= 1
	if linecolor eq 200 then thisthick= 12.0
	if linecolor eq 150 then thisthick= 2.0
	if linecolor eq 100 then thisthick= 10.0
	if linecolor eq 150 then thispsym= 2
	if linecolor eq  50 then thislinest= 2
	if linecolor eq  50 then thisthick= 8.0

	if linecolor eq 1211 then begin & linecolor= 0 & thisthick= 6.0 & thispsym= 3 & thislinest= 0 &  end
	if linecolor eq 1212 then begin & linecolor= 0 & thisthick= 2.0 & thispsym= 3 & thislinest= 1 &  end
	if linecolor eq 1213 then begin & linecolor= 0 & thisthick= 3.0 & thispsym= 3 & thislinest= 2 &  end

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

;================================================================================



pro minor_prof, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "minor_prof, junk"
   print, "  "
   return
endif

smoothlen= 0.1
filename= 'minorprof.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newysize=24, newxsize=12



; -----------------------------
; set up constants
; -----------------------------

bins = 50

xmax = 20.0
xmin = 0.0
ymax = 5.3
ymin = -3.3


xaxistitle= "!6R (kpc)"
yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"


frun1= "/raid12/tcox/Gs/G3G3b-u2"   & snap1= 120    & msg1='G3G3, 1:1'
frun2= "/raid12/tcox/Gs/G3G2-u4"    & snap2= 120    & msg2='G3G2, 2.3:1'
frun3= "/raid12/tcox/Gs/G3G1-u4"    & snap3= 120    & msg3='G3G1, 5.8:1'
frun4= "/raid12/tcox/Gs/G3G0e-u4"   & snap4= 120    & msg4='G3G0, 22.3:1'

;snap1 = 0
;snap2 = 0
;snap3 = 0
;snap4 = 0


x0= 0.18
x1= 0.96
    
y0= 0.08
ys= 0.225   ; assumes 4 panels
y1= y0+ys
y2= y0+ys+ys
y3= y0+ys+ys+ys
y4= y0+ys+ys+ys+ys







; ------------------
;    1
; ------------------
!p.position= [x0, y3, x1, y4]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	;xtickformat='(a1)', ytickformat='(a1)'
	;xtickformat='(a1)', ytitle=yaxistitle
	xtickformat='(a1)', ytitle=' '

;plot_sd_profile_all, frun1, snap1, xmin=xmin, xmax=xmax, bins=bins
plot_sd_profile_newstars, frun1, snap1, xmin=xmin, xmax=xmax, bins=bins, satnumpart=240000

xsmooth=[0.1,0.1]
ysmooth=[ymin,ymax]
oplot, xsmooth, ysmooth, linestyle=1, color= 0

;oplot, [1.2,1.4], [4.6,4.6], psym=-3, color= 0, thick=6.0, linestyle= 1
;xyouts, 1.50, 4.5, 'collisionless', /data, size= 1.5, color= 0, charthick=5.0

;oplot, [1.2,1.4], [4.1,4.1], psym=-3, color= 150, thick=8.0, linestyle= 0
;xyouts, 1.50, 4.0, 'f= ', /data, size= 1.5, color= linecolor, charthick=5.0
xyouts, 11.0, 4.0, msg1, /data, size= 1.5, color= linecolor, charthick=5.0

        ;xx0= 0.65
        ;yy0= 0.95
        ;xyouts, xx0, yy0, '!6Disk', /normal, charthick=3.0, size= 1.7, color= 100
        ;xyouts, xx0, yy0-0.03, 'Bulge', /normal, charthick= 3.0, size= 1.7, color= 200
        ;xyouts, xx0, yy0-0.06, 'New Stars', /normal, charthick= 3.0, size= 1.7, color= 150

xx0= 0.33
yy0= 0.81

oplot, [0.7,2.9], [-0.8,-0.8], psym=-3, color= 0, thick=6.0, linestyle= 0
xyouts, xx0, yy0, '!6Total', /normal, charthick=2.0, size= 1.3, color= 0

oplot, [0.7,2.9], [-1.6,-1.6], psym=-3, color= 0, thick=3.0, linestyle= 1
xyouts, xx0, yy0-0.02, 'Primary', /normal, charthick= 2.0, size= 1.3, color= 0

oplot, [0.7,2.9], [-2.4,-2.4], psym=-3, color= 0, thick=2.0, linestyle= 2
xyouts, xx0, yy0-0.04, 'Satellite', /normal, charthick= 2.0, size= 1.3, color= 0



xyouts, 0.04, 0.43, yaxistitle, /normal, charthick= 3.0, size= 1.7, color= 0, orientation= 90.0


; ------------------
;    2
; ------------------
!p.position= [x0, y2, x1, y3]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	;xtickformat='(a1)', ytickformat='(a1)', /noerase
        ;xtickformat='(a1)', ytitle=yaxistitle, /noerase
        xtickformat='(a1)', ytitle=' ', /noerase

;plot_sd_profile_all, frun2, snap2, xmin=xmin, xmax=xmax, bins=bins
plot_sd_profile_newstars, frun2, snap2, xmin=xmin, xmax=xmax, bins=bins, satnumpart=150000

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 11.0, 4.0, msg2, /data, size= 1.5, color= linecolor, charthick=5.0


; ------------------
;    3
; ------------------
!p.position= [x0, y1, x1, y2]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	;xtickformat='(a1)', ytickformat='(a1)', /noerase
        ;xtickformat='(a1)', ytitle=yaxistitle, /noerase
        xtickformat='(a1)', ytitle=' ', /noerase

;plot_sd_profile_all, frun3, snap3, xmin=xmin, xmax=xmax, bins=bins
plot_sd_profile_newstars, frun3, snap3, xmin=xmin, xmax=xmax, bins=bins, satnumpart=95000

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 11.0, 4.0, msg3, /data, size= 1.5, color= linecolor, charthick=5.0


; ------------------
;    4
; ------------------
!p.position= [x0, y0, x1, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	;xtickformat='(a1)', ytickformat='(a1)', /noerase
        ;xtitle=xaxistitle, ytitle=yaxistitle, /noerase
        xtitle=xaxistitle, ytitle=' ', /noerase

;plot_sd_profile_all, frun4, snap4, xmin=xmin, xmax=xmax, bins=bins
plot_sd_profile_newstars, frun4, snap4, xmin=xmin, xmax=xmax, bins=bins, satnumpart=51000

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 11.0, 4.0, msg4, /data, size= 1.5, color= linecolor, charthick=5.0





; done
; -----
device, /close


end





;================================================================================


;================================================================================















;================================================================================

;================================================================================

;================================================================================



pro minor_prof2, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "minor_prof2, junk"
   print, "  "
   return
endif

smoothlen= 0.1
filename= 'minorprof2.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 20

xmax = 15.0
xmin = 0.0
ymax = 5.3
ymin = -1.9


xaxistitle= "!6R (kpc)"
yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"


;frun= "/raid12/tcox/Gs/G3G3b-u2"   & snap= 120    & msg='G3G3, 1:1'
;frun= "/raid12/tcox/Gs/G3G2-u4"    & snap= 120    & msg='G3G2, 2.3:1'
frun= "/raid12/tcox/Gs/G3G1-u4"    & snap= 120    & msg='G3G1, 5.8:1'
;frun= "/raid12/tcox/Gs/G3G0e-u4"   & snap= 120    & msg='G3G0, 22.3:1'

;snap1 = 0
;snap2 = 0
;snap3 = 0
;snap4 = 0


x0= 0.16
x1= 0.98
    
y0= 0.12
y1= 0.98







; ------------------
; ------------------

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	;xtickformat='(a1)', ytickformat='(a1)'
	;xtickformat='(a1)', ytitle=yaxistitle
        xtitle=xaxistitle, ytitle=yaxistitle
	;xtickformat='(a1)', ytitle=' '

plot_sd_profile_all, frun, snap, xmin=xmin, xmax=xmax, bins=bins

;xsmooth=[0.1,0.1]
;ysmooth=[ymin,ymax]
;oplot, xsmooth, ysmooth, linestyle=1, color= 0

xx0= 0.65
yy0= 0.90

oplot, [7.2,8.2], [4.8,4.8], psym=-3, color= 100, thick=6.0, linestyle= 0
xyouts, xx0, yy0, '!6Disk', /normal, charthick=3.0, size= 1.7, color= 100

oplot, [7.2,8.2], [4.0,4.0], psym=-3, color= 200, thick=6.0, linestyle= 1
xyouts, xx0, yy0-0.09, 'Bulge', /normal, charthick= 3.0, size= 1.7, color= 200

oplot, [7.2,8.2], [3.2,3.2], psym=-2, color= 150, thick=3.0, linestyle= 0
xyouts, xx0, yy0-0.18, 'New Stars', /normal, charthick= 3.0, size= 1.7, color= 150


xyouts, 1.80, -1.0, msg, /data, size= 1.5, color= 0, charthick=5.0



; overplot isolated??
;plot_sd_profile_all, '/raid4/tcox/isolated/G3_n0', 120, xmin=xmin, xmax=xmax, bins=bins




; done
; -----
device, /close


end





;================================================================================


;================================================================================



pro minor_prof3, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "minor_prof3, junk"
   print, "  "
   return
endif

smoothlen= 0.1
filename= 'minorprof3.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 20

xmax = 15.0
xmin = 0.0
;xmin = 0.1
ymax = 4.8
ymin = -1.1


xaxistitle= "!6R (kpc)"
yaxistitle= "!6Log !7R!6(R) !N(M!D!9n!6!N pc!E-2!N)"


frun1= "/raid12/tcox/Gs/G3G3b-u2"   & snap1= 120    & msg1='G3G3, 1:1'
frun2= "/raid12/tcox/Gs/G3G2-u4"    & snap2= 120    & msg2='G3G2, 2.3:1'
frun3= "/raid12/tcox/Gs/G3G1-u4"    & snap3= 120    & msg3='G3G1, 5.8:1'
frun4= "/raid12/tcox/Gs/G3G0e-u4"   & snap4= 120    & msg4='G3G0, 22.3:1'


x0= 0.16
x1= 0.98
    
y0= 0.12
y1= 0.98



; ------------------
; ------------------

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	;xtickformat='(a1)', ytickformat='(a1)'
	;xtickformat='(a1)', ytitle=yaxistitle
        xtitle=xaxistitle, ytitle=yaxistitle
        ;xtitle=xaxistitle, ytitle=yaxistitle, /xlog
	;xtickformat='(a1)', ytitle=' '

plot_sd_profile_newstars, frun1, snap1, xmin=xmin, xmax=xmax, bins=bins, thiscolor= 100
plot_sd_profile_newstars, frun2, snap2, xmin=xmin, xmax=xmax, bins=bins, thiscolor= 200
plot_sd_profile_newstars, frun3, snap3, xmin=xmin, xmax=xmax, bins=bins, thiscolor= 150
plot_sd_profile_newstars, frun4, snap4, xmin=xmin, xmax=xmax, bins=bins, thiscolor= 50

;plot_sd_profile_newstars, '/raid4/tcox/isolated/G3_n0', 120, xmin=xmin, xmax=xmax, bins=bins, thiscolor= 0

;xsmooth=[0.1,0.1]
;ysmooth=[ymin,ymax]
;oplot, xsmooth, ysmooth, linestyle=1, color= 0

xx0= 0.65
yy0= 0.90

oplot, [7.0,8.4], [4.35,4.35], psym=-3, color= 100, thick=10.0, linestyle= 0
xyouts, xx0, yy0, msg1, /normal, charthick=3.0, size= 1.7, color= 100

oplot, [7.0,8.4], [4.03,4.03], psym=-3, color= 200, thick=12.0, linestyle= 1
xyouts, xx0, yy0-0.05, msg2, /normal, charthick= 3.0, size= 1.7, color= 200

oplot, [7.0,7.7,8.4], [3.66,3.66,3.66], psym=-2, color= 150, thick=3.0, linestyle= 0
xyouts, xx0, yy0-0.10, msg3, /normal, charthick= 3.0, size= 1.7, color= 150

oplot, [7.0,8.4], [3.33,3.33], psym=-3, color= 50, thick=8.0, linestyle= 2
xyouts, xx0, yy0-0.15, msg4, /normal, charthick= 3.0, size= 1.7, color= 50



; done
; -----
device, /close


end


















;================================================================================



pro star_loc, junk

do_star_loc, "/raid12/tcox/Gs/G3G3b-u2", 120
do_star_loc, "/raid12/tcox/Gs/G3G2-u4", 120
do_star_loc, "/raid12/tcox/Gs/G3G1-u4", 120
do_star_loc, "/raid12/tcox/Gs/G3G0e-u4", 120
do_star_loc, "/raid4/tcox/isolated/G3_n0", 120

end


; -----------------------------------------------


pro do_star_loc, frun, snapnum


ok=fload_snapshot(frun,snapnum)


r=fload_newstars_xyz('r')
m=fload_newstars_mass(1)
mg=fload_gas_mass(1)

rcheck= 0.5
;rcheck= 1.0
;rcheck= 2.0
;rcheck= 3.0
;rcheck= 4.0

in_idx= where(r lt rcheck)
out_idx= where(r ge rcheck)

print, "innerf (newstars)= ", total(m(in_idx))/total(m)
print, "innerf (gasmass)= ", total(m(in_idx))/(total(m)+total(mg))



end







