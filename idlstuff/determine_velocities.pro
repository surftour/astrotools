;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   Procedures to determine the velocity field
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------





pro do_all, junk




end




;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Beta Profile
;     -------------------------------------------
;     Two methods to do this, mine and felix's
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro beta_profile, junk, smoothlen=smoothlen, filename=filename, $
				snapnum=snapnum, $
				frun=frun

if not keyword_set(junk) then begin
   print, "  "
   print, "beta_profile, junk, smoothlen=smoothlen, filename=filename, snapnum=snapnum"
   print, "  "
   return
endif

if not keyword_set(smoothlen) then smoothlen= 0.1
if not keyword_set(filename) then filename='beta.eps'
if not keyword_set(snapnum) then snapnum=25


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 50

; radius, linear scale
xmax = 20.0
xmin =  0.0

ymax = 1.2
ymin = -2.0


xaxistitle= "r (h!E-1!Nkpc)"
yaxistitle= "!7b!3"


if not keyword_set(snapnum) then snapnum= 25

; ------------------
; Plot this up
; ------------------
!p.position= [0.2, 0.15, 0.95, 0.95]

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

	frun= "/raid4/tcox/vc3vc3e"
	;frun= "/raid4/tcox/vc3vc3h"
	process_one_beta_profile, frun, snapnum, xmin, xmax, bins, linecolor=50
	xyouts, 0.55, 0.35, '40% gas, black hole', size=1.33, color=50, /normal, charthick=4.0

	frun= "/raid4/tcox/cvc3vc3e"
	;frun= "/raid4/tcox/cvc3vc3h"
	process_one_beta_profile, frun, snapnum, xmin, xmax, bins, linecolor=150
	xyouts, 0.55, 0.30, 'collisionless', size=1.33, color=150, /normal, charthick=4.0


	;xyouts, 0.75, 0.90, fload_timelbl(1,2), size=1.5, color=0, /normal, charthick=4.0


; print extras
; -------------

; zero
x=[xmin,xmax]
y=x*0.0
oplot, x, y, linestyle= 1, color= 0

; softening length
;smoothlen= 0.2
;x=[smoothlen,smoothlen]
;x=x^(0.25)
;y=[ymin,ymax]
;oplot, x, y, linestyle=1, color= 0

;smoothlen= 0.4
;x=[smoothlen,smoothlen]
;x=x^(0.25)
;y=[ymin,ymax]
;oplot, x, y, linestyle=1, color= 0




; done
; -----
device, /close


end






; --------------------------
;  process beta profile
; --------------------------
pro process_one_beta_profile, frun, snapnum, xmin, xmax, bins, linecolor=linecolor, msg=msg


ok=fload_snapshot_bh(frun,snapnum)
radius= fload_allstars_xyz('r')
comvel= fload_all_comvel(1)
vel_radial= fload_allstars_v('r',comvel=comvel)
vel_theta= fload_allstars_v('theta',comvel=comvel)
vel_phi= fload_allstars_v('phi',comvel=comvel)
vel_tan= fload_allstars_v('tan',comvel=comvel)


; tj's method
; -------------
rs= fltarr(bins)
betas= fltarr(bins)

binsize = float((xmax-xmin))/bins

for i=1,bins do begin
	lg_r = i*binsize + xmin
	sm_r = (i-1)*binsize + xmin
	rs(i-1) = 0.5*(lg_r + sm_r)

	idx= where((radius GE sm_r) AND (radius LT lg_r))
	if idx(0) lt 0 then begin
		betas(i-1)= 5.0   ; flag for no data
	endif else begin
		;vtan= sqrt(variance(vel_tan(idx)))
		vrad= sqrt(variance(vel_radial(idx)))
		vtha= sqrt(variance(vel_theta(idx)))
		vphi= sqrt(variance(vel_phi(idx)))
		vtan= sqrt(vtha*vtha + vphi*vphi)      ; alternate calculation of v_tangential
		betas(i-1)= 1.0 - (vtan*vtan)/(2.0*vrad*vrad)
	endelse
;print, "i= ",i, "  vtan= ", vtan, "  vrad= ",vrad
;print, "i= ",i, " vtan= ", vtan, "  sqrt(vtha^2 + vphi^2)= ",sqrt(vtha*vtha + vphi*vphi)

endfor

; take out where there was no data
idx=where(betas lt 5.0)
if idx(0) ne -1 then begin
    rs= rs(idx)
    betas= betas(idx)
endif

oplot, rs, betas, psym=-3, color=linecolor, linestyle=0, thick=4.0




; felix's method
; ----------------
; I've checked a time or two, and this
; gives the same answer as my method.
;
one= radius
one(*)= 1.0

binsize= float((xmax-xmin))/bins

vrad= hist1d(radius, vel_radial*vel_radial, obin=x, binsize=binsize, max=xmax)
vrad= sqrt(vrad / hist1d(radius, one, obin=x, binsize=binsize, max=xmax))

vtan= hist1d(radius, vel_tan*vel_tan, obin=x, binsize=binsize, max=xmax)
vtan= sqrt(vtan / hist1d(radius, one, obin=x, binsize=binsize, max=xmax))

betas= 1.0 - (vtan*vtan) / (2.0*vrad*vrad)

;print, "v_rad"
;print, "-----"
;print, vrad
;print, "-----"
;print, "v_tan"
;print, vtan
;print, "-----"

;oplot, x, betas, psym=-3, color=linecolor, linestyle=1, thick=4.0



end







;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Velocity Dispersion Profiles
;     -------------------------------------------
;     Two methods to do this, mine and felix's
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro sig3D_profile, junk, smoothlen=smoothlen, filename=filename, $
				snapnum=snapnum, $
				frun=frun

if not keyword_set(junk) then begin
   print, "  "
   print, "sig3D_profile, junk, smoothlen=smoothlen, filename=filename, snapnum=snapnum"
   print, "  "
   return
endif

if not keyword_set(smoothlen) then smoothlen= 0.1
if not keyword_set(filename) then filename='sig3D.eps'
if not keyword_set(snapnum) then snapnum=25


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 50

; radius, linear scale
xmax = 20.0
xmin =  0.0

ymax = 300.0
ymin = 0.0


xaxistitle= "r (h!E-1!Nkpc)"
yaxistitle= "!7r!3"



; ------------------
; Plot this up
; ------------------
!p.position= [0.2, 0.15, 0.95, 0.95]

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

	frun= "/raid4/tcox/vc3vc3e"
	;frun= "/raid4/tcox/vc3vc3h"
	process_one_sig3D_profile, frun, snapnum, xmin, xmax, bins, linecolor=50
        ;xyouts, 0.65, 0.85, 'vc3vc3e', size=1.5, color=50, /normal, charthick=4.0
        xyouts, 0.55, 0.85, '40% gas, black hole', size=1.33, color=50, /normal, charthick=4.0

	;frun= "/raid4/tcox/cvc3vc3e"
	;process_one_sig3D_profile, frun, snapnum, xmin, xmax, bins, linecolor=50
        ;xyouts, 0.55, 0.80, 'cvc3vc3e', size=1.5, color=150, /normal, charthick=4.0
        ;xyouts, 0.55, 0.85, 'collisionless', size=1.33, color=50, /normal, charthick=4.0


	;xyouts, 0.75, 0.90, fload_timelbl(1,2), size=1.33, color=0, /normal, charthick=4.0


; print extras
; -------------



; done
; -----
device, /close


end






; --------------------------
;  process beta profile
; --------------------------
pro process_one_sig3D_profile, frun, snapnum, xmin, xmax, bins, linecolor=linecolor, msg=msg


ok=fload_snapshot_bh(frun,snapnum)

; allstars
;
radius= fload_allstars_xyz('r')
comvel= fload_all_comvel(1)
vel_radial= fload_allstars_v('r',comvel=comvel)
vel_theta= fload_allstars_v('theta',comvel=comvel)
vel_phi= fload_allstars_v('phi',comvel=comvel)
vel_tan= fload_allstars_v('tan',comvel=comvel)
vel_tot= fload_allstars_v('tot',comvel=comvel)

; dark matter
;
radius_dm= fload_halo_xyz('r')
vel_radial_dm= fload_halo_v('r',comvel=comvel)
vel_theta_dm= fload_halo_v('theta',comvel=comvel)
vel_phi_dm= fload_halo_v('phi',comvel=comvel)
vel_tan_dm= fload_halo_v('tan',comvel=comvel)
vel_tot_dm= fload_halo_v('tot',comvel=comvel)


; tj's method
; -------------
rs= fltarr(bins)
vtan= fltarr(bins)
vrad= fltarr(bins)
vtot= fltarr(bins)

rs_dm= fltarr(bins)
vtan_dm= fltarr(bins)
vrad_dm= fltarr(bins)
vtot_dm= fltarr(bins)

binsize = float((xmax-xmin))/bins

for i=1,bins do begin
	lg_r = i*binsize + xmin
	sm_r = (i-1)*binsize + xmin
	rs(i-1) = 0.5*(lg_r + sm_r)
	rs_dm(i-1) = 0.5*(lg_r + sm_r)

	; all stars
	idx= where((radius GE sm_r) AND (radius LT lg_r))
	if idx(0) le 2 then begin
		vtan(i-1)= -5.0   ; flag for no data
	endif else begin
		;vtan(i-1)= sqrt(variance(vel_tan(idx)))
		vtot(i-1)= sqrt(variance(vel_tot(idx)))
		vrad(i-1)= sqrt(variance(vel_radial(idx)))
		vtha= sqrt(variance(vel_theta(idx)))
		vphi= sqrt(variance(vel_phi(idx)))
		vtan(i-1)= sqrt(vtha*vtha + vphi*vphi)      ; alternate calculation of v_tangential
	endelse

	; dark matter
	idx= where((radius_dm GE sm_r) AND (radius_dm LT lg_r))
	if idx(0) le 2 then begin
		vtan_dm(i-1)= -5.0   ; flag for no data
	endif else begin
		;vtan_dm(i-1)= sqrt(variance(vel_tan_dm(idx)))
		vtot_dm(i-1)= sqrt(variance(vel_tot_dm(idx)))
		vrad_dm(i-1)= sqrt(variance(vel_radial_dm(idx)))
		vtha= sqrt(variance(vel_theta_dm(idx)))
		vphi= sqrt(variance(vel_phi_dm(idx)))
		vtan_dm(i-1)= sqrt(vtha*vtha + vphi*vphi)      ; alternate calculation of v_tangential
	endelse

endfor

; take out where there was no data
idx=where(vtan eq -5.0)
if idx(0) ne -1 then begin
    rs= rs(idx)
    vtan= vtan(idx)
    vrad= vrad(idx)
    vtot= vtot(idx)
endif

idx=where(vtan_dm eq -5.0)
if idx(0) ne -1 then begin
    rs_dm= rs_dm(idx)
    vtan_dm= vtan_dm(idx)
    vrad_dm= vrad_dm(idx)
    vtot_dm= vtot_dm(idx)
endif



; all the same color
oplot, rs, vtan, psym=-3, color=linecolor, linestyle=2, thick=4.0
oplot, rs, vrad, psym=-3, color=linecolor, linestyle=1, thick=4.0
;oplot, rs, vtot, psym=-5, color=linecolor, linestyle=0, thick=4.0
oplot, rs, vtot, psym=-3, color=linecolor, linestyle=0, thick=4.0

xyouts, 0.27, 0.25, 'stars', /normal, charthick= 3.0, size=1.2, color=linecolor

; manually set colors
;oplot, rs, vtan, psym=-3, color=100, linestyle=2, thick=4.0
;oplot, rs, vrad, psym=-3, color=150, linestyle=1, thick=4.0
;oplot, rs, vtot, psym=-5, color=0, linestyle=0, thick=4.0
;xyouts, 0.27, 0.31, 'tangential', /normal, charthick= 3.0, size=1.2, color=150
;xyouts, 0.27, 0.27, 'radial', /normal, charthick= 3.0, size=1.2, color=100
;xyouts, 0.27, 0.23, 'total', /normal, charthick= 3.0, size=1.2, color=0

oplot, rs_dm, vtan_dm, psym=-3, color=linecolor+100, linestyle=2, thick=4.0
oplot, rs_dm, vrad_dm, psym=-3, color=linecolor+100, linestyle=1, thick=4.0
oplot, rs_dm, vtot_dm, psym=-3, color=linecolor+100, linestyle=0, thick=4.0

xyouts, 0.27, 0.20, 'dark matter', /normal, charthick= 3.0, size=1.2, color=linecolor+100




; felix's method
; ----------------
one= radius
one(*)= 1.0

binsize= float((xmax-xmin))/bins

vrad= hist1d(radius, vel_radial*vel_radial, obin=x, binsize=binsize, max=xmax)
vrad= sqrt(vrad / hist1d(radius, one, obin=x, binsize=binsize, max=xmax))

vtan= hist1d(radius, vel_tan*vel_tan, obin=x, binsize=binsize, max=xmax)
vtan= sqrt(vtan / hist1d(radius, one, obin=x, binsize=binsize, max=xmax))

betas= 1.0 - (vtan*vtan) / (2.0*vrad*vrad)

;oplot, x, betas, psym=-3, color=linecolor, linestyle=1, thick=4.0



end









;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     2D Sigma Profile
;     -------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro sigma_profile, junk, smoothlen=smoothlen, filename=filename, $
				snapnum=snapnum, $
				frun=frun

if not keyword_set(junk) then begin
   print, "  "
   print, "sigma_profile, junk, smoothlen=smoothlen, filename=filename, snapnum=snapnum"
   print, "  "
   return
endif

if not keyword_set(smoothlen) then smoothlen= 0.1
if not keyword_set(filename) then filename='sigma.eps'
if not keyword_set(snapnum) then snapnum=25


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=1
;setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 50

; radius, linear scale
xmax = 40.0
xmin =  0.0

ymax = 200.0
ymin = 0.0


;xaxistitle= "R (h!E-1!Nkpc)"
xaxistitle= "R (kpc)"
yaxistitle= "!7r!3!DP!N (R)"


if not keyword_set(snapnum) then snapnum= 25

; ------------------
; Plot this up
; ------------------
!p.position= [0.2, 0.15, 0.95, 0.95]

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

	;frun= "/raid4/tcox/vc3vc3e"
	;frun= "/raid4/tcox/vc3vc3h"
	;process_one_2Dsig_profile, frun, snapnum, xmin, xmax, bins, linecolor=50, msg=fload_getid(frun)
        ;xyouts, 0.55, 0.85, 'vc3vc3e', size=1.5, color=50, /normal, charthick=4.0
        ;xyouts, 0.55, 0.85, '40% gas, black hole', size=1.33, color=50, /normal, charthick=4.0

	frun= "/raid4/tcox/cvc3vc3e"
	;frun= "/raid4/tcox/cvc3vc3h"
	process_one_2Dsig_profile, frun, snapnum, xmin, xmax, bins, linecolor=150, msg=fload_getid(frun)
        ;xyouts, 0.55, 0.80, 'cvc3vc3e', size=1.5, color=150, /normal, charthick=4.0
        xyouts, 0.55, 0.85, 'collisionless', size=1.5, color=150, /normal, charthick=4.0


	;xyouts, 0.75, 0.90, fload_timelbl(1,2), size=1.5, color=0, /normal, charthick=4.0


; print extras
; -------------

; zero
x=[xmin,xmax]
y=x*0.0
oplot, x, y, linestyle= 1, color= 0

; softening length
;smoothlen= 0.2
;x=[smoothlen,smoothlen]
;x=x^(0.25)
;y=[ymin,ymax]
;oplot, x, y, linestyle=1, color= 0

;smoothlen= 0.4
;x=[smoothlen,smoothlen]
;x=x^(0.25)
;y=[ymin,ymax]
;oplot, x, y, linestyle=1, color= 0




; done
; -----
device, /close


end







; --------------------------
;  process 2D sigma profile
; --------------------------
pro process_one_2Dsig_profile, frun, snapnum, xmin, xmax, bins, $
				linecolor=linecolor, msg=msg


ok=fload_snapshot_bh(frun,snapnum)
x= fload_allstars_xyz('x')
y= fload_allstars_xyz('y')
z= fload_allstars_xyz('z')
comvel= fload_all_comvel(1)
vel_x= fload_allstars_v('x',comvel=comvel)
vel_y= fload_allstars_v('y',comvel=comvel)
vel_z= fload_allstars_v('z',comvel=comvel)


; felix's method
; -------------
rs= fltarr(bins)
vx_sig= fltarr(bins)
vy_sig= fltarr(bins)
vz_sig= fltarr(bins)

binsize = float((xmax-xmin))/bins

for i=1,bins do begin
	lg_r = i*binsize + xmin
	sm_r = (i-1)*binsize + xmin
	rs(i-1) = 0.5*(lg_r + sm_r)

	; xy direction
	radius= sqrt(x*x + y*y)
	idx= where((radius GE sm_r) AND (radius LT lg_r))
	if n_elements(idx) gt 2 then vz_sig(i-1)= sqrt(variance(vel_z(idx)))

        ; xz direction
        radius= sqrt(x*x + z*z)
        idx= where((radius GE sm_r) AND (radius LT lg_r))
        if n_elements(idx) gt 2 then vy_sig(i-1)= sqrt(variance(vel_y(idx)))

        ; yz direction
        radius= sqrt(y*y + z*z)
        idx= where((radius GE sm_r) AND (radius LT lg_r))
        if n_elements(idx) gt 2 then vx_sig(i-1)= sqrt(variance(vel_x(idx)))

endfor

; take out where there was no data
;idx=where(betas lt 5.0)
;if idx(0) ne -1 then begin
;    rs= rs(idx)
;    betas= betas(idx)
;endif

oplot, rs, vx_sig, psym=-3, color=linecolor, linestyle=0, thick=4.0
oplot, rs, vy_sig, psym=-3, color=linecolor, linestyle=1, thick=4.0
oplot, rs, vz_sig, psym=-3, color=linecolor, linestyle=2, thick=4.0



xyouts, 0.27, 0.22, 'three projections', size=1.2, color=0, /normal, charthick=4.0


end










; ==========================================================================================

; ==========================================================================================






; We here test our procedure to get the nearest neighbors to any
; point in space, or equivalently for each particles.
; ---------------------------------------------------------------

;pro calc_velfield, junk
pro calc_velfield, frun


if not keyword_set(frun) then begin
	print, "  "
	print, " calc_velfield, frun"
	print, "  "
	return
endif

;frun= "/raid4/tcox/vc3vc3e"
;snapnum= 25
snapnum= 31

ok=fload_snapshot_bh(frun,snapnum)

x= fload_allstars_xyz('x')
y= fload_allstars_xyz('y')
z= fload_allstars_xyz('z')
m= fload_allstars_mass(1)
ids= fload_allstars_ids(1)
comvel= fload_all_comvel(1)
vx= fload_allstars_v('x',comvel=comvel)
vy= fload_allstars_v('y',comvel=comvel)
vz= fload_allstars_v('z',comvel=comvel)

N= long(n_elements(x))

;Coord=fltarr(3,N)
;Coord=fltarr(6,N)
Coord=fltarr(7,N)
Masses= fltarr(N)
;Id= lonarr(N)
Velx= fltarr(N)
Vely= fltarr(N)
Velz= fltarr(N)

Coord(0,*)= x
Coord(1,*)= y
Coord(2,*)= z
Masses(*)= m
;Id(*)= ids
;Velx(*)= vx
;Vely(*)= vy
;Velz(*)= vz
Coord(3,*)= vx
Coord(4,*)= vy
Coord(5,*)= vz
Coord(6,*)= ids


DesNgb=96L
Hmax= 100.0



if(N gt 0) then begin

        VelEllipsoid= fltarr(6,N)      ;  full velocity field

        S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/VelocityEllipsoid/velocityellipsoid.so', $
                'velocityellipsoid', $
		N, $
                Coord, $
                Masses, $
                ;Velx, $
                ;Vely, $
                ;Velz, $
                ;Id, $
		DesNgb, $
		Hmax, $
                VelEllipsoid)
endif else begin
        print,'No stars, no problem.'
        return
endelse

help, VelEllipsoid

AvgVx= VelEllipsoid[0,*]
AvgDispx= VelEllipsoid[1,*]
AvgVy= VelEllipsoid[2,*]
AvgDispy= VelEllipsoid[3,*]
AvgVz= VelEllipsoid[4,*]
AvgDispz= VelEllipsoid[5,*]

;print, "--------------------------"
;print, "PI_xx   ", total(AvgDispx)
;print, "PI_yy   ", total(AvgDispy)
;print, "PI_zz   ", total(AvgDispz)
;print, "--------------------------"

openw, 1, frun+'/tvt.txt', ERROR=err
printf, 1, "#   "
printf, 1, "# Tensor Virial Theorem   "
printf, 1, "#   "
printf, 1, "#   "
printf, 1, "#   "

sigx2= vx-AvgVx
printf, 1, "PI_xx    ", total(m*sigx2*sigx2)

sigy2= vy-AvgVy
printf, 1, "PI_yy    ", total(m*sigy2*sigy2)

sigz2= vz-AvgVz
printf, 1, "PI_zz    ", total(m*sigz2*sigz2)

printf, 1, "PI_xy    ", total(m*sigx2*sigy2)
printf, 1, "PI_xz    ", total(m*sigx2*sigz2)
printf, 1, "PI_yz    ", total(m*sigy2*sigz2)


printf, 1, "T_xx    ", total(m*AvgVx*AvgVx)
printf, 1, "T_yy    ", total(m*AvgVy*AvgVy)
printf, 1, "T_zz    ", total(m*AvgVz*AvgVz)
printf, 1, "T_xy    ", total(m*AvgVx*AvgVy)
printf, 1, "T_xz    ", total(m*AvgVx*AvgVz)
printf, 1, "T_yz    ", total(m*AvgVy*AvgVz)

close, 1

;print, "--------------------------"

end





; =============================================================================





pro read_tvt_file, frun, pi_xx, pi_yy, pi_zz, $
			pi_xy, pi_xz, pi_yz, $
			t_xx, t_yy, t_zz, $
			t_xy, t_xz, t_yz


tvtfile= '/raid4/tcox/'+frun+'/tvt.txt'

openr, 1, tvtfile
junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

pi_xx= 0.0 & pi_yy= 0.0 & pi_zz= 0.0
pi_xy= 0.0 & pi_xz= 0.0 & pi_yz= 0.0
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & pi_xx= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & pi_yy= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & pi_zz= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & pi_xy= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & pi_xz= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & pi_yz= float(tempjunk(1))

t_xx= 0.0 & t_yy= 0.0 & t_zz= 0.0
t_xy= 0.0 & t_xz= 0.0 & t_yz= 0.0
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_xx= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_yy= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_zz= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_xy= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_xz= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_yz= float(tempjunk(1))

close, 1


end




