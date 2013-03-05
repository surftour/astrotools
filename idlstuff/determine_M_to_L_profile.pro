;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   More Profiles
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------




;------------------------------------------------------------------------
;
;  The below procedure is **NOT** itentical to that in 
; determine_profiles.  This one calculates the total mass
; (i.e. dark and luminous and gaseous) and total L within
; a fixed aperture, and then divides the two.
;
;
;------------------------------------------------------------------------

pro process_M_to_L_formanyprojections, lx, ly, lz, llum, $
					mx, my, mz, mm, $
					xmin, xmax, bins, $
					result_avg, result_1sig, xaxisvals, $
					center=center, $
					x_is_log=x_is_log, x_is_devac=x_is_devac, $
					average=average, $
					y_weighting= y_weighting


;n_projections= 100    ; maybe do 1000 at some point
;n_projections= 40
;n_projections= 10
n_projections= 1
seed= 154L

; 2D array which stores projected quantity for
; each projection
temp_average= fltarr(n_projections+1,bins)

if keyword_set(center) then c=center else c=[0,0,0]


; for mass, we don't project
mass_r= sqrt((mx-c[0])*(mx-c[0]) + (my-c[1])*(my-c[1]) + (mz-c[2])*(mz-c[2]))
if keyword_set(x_is_log) then mass_r= alog10(mass_r)


binsize = float((xmax-xmin))/bins

; ---------------------------
;  begin projection loop
; ---------------------------
for i=0,n_projections do begin

        rdphi= randomu(seed)
        rdtheta= randomu(seed)

        ; in radians
	theta= rdtheta*!PI
	phi= rdphi*2*!PI


        ; rotate
        rot_x= (lx-c[0])*(cos(theta)*cos(phi)) + (ly-c[1])*(cos(theta)*sin(phi)) + (lz-c[2])*sin(theta)
        rot_y= -(lx-c[0])*sin(theta) + (ly-c[1])*cos(theta)
	rot_r= sqrt(rot_x*rot_x + rot_y*rot_y)

	if keyword_set(x_is_log) then rot_r=alog10(rot_r)

	; calculate the profile for this projection
	thisrad= fltarr(bins)
	thisquant= fltarr(bins)


	; go through the bins
	; --------------------------------------
        for ii=1,bins do begin
           lg_r = ii*binsize + xmin
           thisrad(ii-1) = lg_r

	   ; enclosed light
           idx= where((rot_r LE lg_r))
           if idx(0) lt 0 then begin
		enclosed_light= 1.0
	   endif else begin
		enclosed_light= total(llum(idx))
	   endelse

	   ; enclosed mass
           idx= where((mass_r LE lg_r))
           if idx(0) lt 0 then begin
		enclosed_mass= 0.0
	   endif else begin
		enclosed_mass= total(mm(idx))
	   endelse

;print, "i= ",ii, enclosed_mass, enclosed_light
	   thisquant(ii-1)= enclosed_mass/enclosed_light
        endfor
        ; ---------------------------------------


	temp_average[i,*]= thisquant

	xaxisvals= thisrad
	
        if (i mod 100) eq 0 then print, "i= ",i

endfor

for i=0,bins-1 do begin
	result_moment= moment(temp_average[*,i])
	result_avg[i]= result_moment[0]
	result_1sig[i]= sqrt(result_moment[1])
endfor


end








;================================================================================









;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Compute M/L Profiles
;     ---------------------------------------------
;     This uses the above procedure(s) and the new color
;    code we've developed (with the help of Brant and Paul).
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;  Main procedure to do a comparison of M/L ratios for multiple
;  simulations.
;
;
; =======================================================================================
pro plot_M_to_L_comparison, junk, $
			filename=filename, $
			smoothlen=smoothlen, $
			snapnum=snapnum

if not keyword_set(junk) then begin
   print, "  "
   print, "plot_M_to_L_comparison, junk"
   print, "  "
   print, "  "
   print, "  warning: needs determine_lums.pro loaded"
   print, "  "
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename="MtoL.eps"
if not keyword_set(snapnum) then snapnum= 25

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 50

; x is radius in kpc/h
;  (std is linear)
xmax = 10.0
xmin = 0.0

makeitlog= 0
if makeitlog eq 1 then begin
	xmin = 0.1
endif

; y is M/L ratio
ymax = 5.0
ymin = 0.0


xtit="R (h!E-1!Nkpc)"
;xtit="R (kpc)"

;ytit="M/L!DU!N"
;ytit="M/L!DB!N"
ytit="M/L!DK!N"


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



  do_regular= 0
  if do_regular eq 1 then begin
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3e", 30, linecolor=50
	xyouts, 0.25, 0.85, 'black hole', /normal, size= 1.33, color=50, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3e_no", 107, linecolor=150
	xyouts, 0.25, 0.80, 'no black hole', /normal, size= 1.33, color=150, charthick=3.0
  endif


  do_concentrations= 0
  if do_concentrations eq 1 then begin
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3e_c5", 30, linecolor=40
	xyouts, 0.25, 0.90, 'c=5', /normal, size= 1.33, color=40, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3e_c7", 30, linecolor=80
	xyouts, 0.25, 0.86, 'c=7', /normal, size= 1.33, color=80, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3e", 30, linecolor=120
	xyouts, 0.25, 0.82, 'c=9', /normal, size= 1.33, color=120, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3e_c11", 30, linecolor=160
	xyouts, 0.25, 0.78, 'c=11', /normal, size= 1.33, color=160, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3e_c13", 30, linecolor=200
	xyouts, 0.25, 0.74, 'c=13', /normal, size= 1.33, color=200, charthick=3.0
  endif



  do_orbits_and_orient= 0
  if do_orbits_and_orient eq 1 then begin
	snapnum= 30
	xyouts, 0.25, 0.75, 'orientations', /normal, size= 1.33, color=0, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3b", snapnum, linecolor=20
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3c", snapnum, linecolor=30
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3d", snapnum, linecolor=40
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3e", snapnum, linecolor=50
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3f", snapnum, linecolor=60
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3g", snapnum, linecolor=70
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3h", snapnum, linecolor=80
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3i", snapnum, linecolor=90
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3j", snapnum, linecolor=100
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3k", snapnum, linecolor=110
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3l", snapnum, linecolor=120
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3m", snapnum, linecolor=130
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3n", snapnum, linecolor=140
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3o", snapnum, linecolor=150
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3p", snapnum, linecolor=160
   endif

   do_larger_vcs= 1
   if do_larger_vcs eq 1 then begin
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3h", 30, linecolor= 50
	xyouts, 0.25, 0.90, 'vc3', /normal, size= 1.33, color=50, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc4vc4a", 30, linecolor= 100
	xyouts, 0.25, 0.86, 'vc4', /normal, size= 1.33, color=100, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc5vc5a", 30, linecolor= 150
	xyouts, 0.25, 0.82, 'vc5', /normal, size= 1.33, color=150, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc6vc6a", 30, linecolor= 200
	xyouts, 0.25, 0.78, 'vc6', /normal, size= 1.33, color=200, charthick=3.0
   endif


   do_bhseed_masses= 0
   if do_bhseed_masses eq 1 then begin
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3avc3a_3", 20, linecolor= 40
	xyouts, 0.25, 0.90, '1e3', /normal, size= 1.33, color=40, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3avc3a_2", 20, linecolor= 80
	xyouts, 0.25, 0.86, '1e4', /normal, size= 1.33, color=80, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3avc3a_1", 20, linecolor= 120
	xyouts, 0.25, 0.82, '1e5', /normal, size= 1.33, color=120, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3avc3a_4", 20, linecolor= 160
	xyouts, 0.25, 0.78, '1e6', /normal, size= 1.33, color=160, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3avc3a_6", 20, linecolor= 200
	xyouts, 0.25, 0.74, '1e7', /normal, size= 1.33, color=200, charthick=3.0
   endif


   do_vc3_time= 0
   if do_vc3_time eq 1 then begin
	frun= "/raid4/tcox/vc3vc3e"
	xyouts, 0.25, 0.75, 'time sequence', /normal, size= 1.33, color=0, charthick=3.0
	snapnum= 15
	initlc = 20
	do_one_galaxy_MtoLratio, frun, snapnum, linecolor= initlc
	do_one_galaxy_MtoLratio, frun, snapnum+1, linecolor= initlc+10
	do_one_galaxy_MtoLratio, frun, snapnum+2, linecolor= initlc+20
	do_one_galaxy_MtoLratio, frun, snapnum+3, linecolor= initlc+30
	do_one_galaxy_MtoLratio, frun, snapnum+4, linecolor= initlc+40
	do_one_galaxy_MtoLratio, frun, snapnum+5, linecolor= initlc+50
	do_one_galaxy_MtoLratio, frun, snapnum+6, linecolor= initlc+60
	do_one_galaxy_MtoLratio, frun, snapnum+7, linecolor= initlc+70
	do_one_galaxy_MtoLratio, frun, snapnum+8, linecolor= initlc+80
	do_one_galaxy_MtoLratio, frun, snapnum+9, linecolor= initlc+90
	do_one_galaxy_MtoLratio, frun, snapnum+10, linecolor= initlc+100
	do_one_galaxy_MtoLratio, frun, snapnum+11, linecolor= initlc+110
	do_one_galaxy_MtoLratio, frun, snapnum+12, linecolor= initlc+120
	do_one_galaxy_MtoLratio, frun, snapnum+13, linecolor= initlc+130
	do_one_galaxy_MtoLratio, frun, snapnum+14, linecolor= initlc+140
	do_one_galaxy_MtoLratio, frun, snapnum+15, linecolor= initlc+150
   endif


   do_vc3_gf= 0
   if do_vc3_gf eq 1 then begin
	; vc3vc3h's
	; ----------
	do_one_galaxy_MtoLratio, "/raid4/tcox/cvc3vc3h", 19, linecolor= 50
	xyouts, 0.25, 0.90, '0%', /normal, size= 1.33, color=50, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3", 19, linecolor= 100
	xyouts, 0.25, 0.86, '20%', /normal, size= 1.33, color=100, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3h", 19, linecolor= 150
	xyouts, 0.25, 0.82, '40%', /normal, size= 1.33, color=150, charthick=3.0
	do_one_galaxy_MtoLratio, "/raid4/tcox/As/A3", 140, linecolor= 200
	xyouts, 0.25, 0.78, '100%', /normal, size= 1.33, color=200, charthick=3.0

	;do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3e_c11", 30, linecolor=160
	;xyouts, 0.25, 0.78, 'c=11', /normal, size= 1.33, color=160, charthick=3.0
	;do_one_galaxy_MtoLratio, "/raid4/tcox/vc3vc3e_c13", 30, linecolor=200
	;xyouts, 0.25, 0.74, 'c=13', /normal, size= 1.33, color=200, charthick=3.0
   endif


; put time and label on there
; ----------------------------
;if keyword_set(smoothlen) then begin
;	x=[smoothlen,smoothlen]
;	x=x^(0.25)
;	y=[ymin,ymax]
;	oplot, x, y, linestyle=1, color= 0
;endif


;xyouts, 0.70, 0.85, fload_timelbl(1,2), /normal, size= 1.5, color= 0, charthick=3.0


; done, close up shop
; --------------------
device, /close



end






; comparison calls this to process one galaxy
; at a time
; ----------------------------------------------

pro do_one_galaxy_MtoLratio, frun, snapnum, linecolor=linecolor, msg=msg


; make sure these are the same as the
; above calling procedure
;
bins = 20

makeitdevac= 1

; x is radius in kpc/h
;xmax = 20.0
xmax = 10.0
xmin = 0.0

makeitlog= 0
if makeitlog eq 1 then begin
	xmin = 0.1
endif

; y is surface brightness
ymax = 5.0
ymin = 0.0


ok= fload_snapshot_bh(frun,snapnum)
TTime= float(fload_time(1))


N= fload_npart(2)+fload_npart(3)+fload_npart(4)
m= 1.0e+10*fload_allstars_mass(1)
age=fload_allstars_age(1)
age=float(TTime-age)
zmets=fload_allstars_z(1)


; get the luminosities
print, "load luminosities"
load_all_stellar_luminosities, N, TTime, m, age, zmets, $
	Lum_Bol, Lum_U, Lum_B, Lum_V, Lum_R, Lum_I, Lum_J, Lum_H, Lum_K, $
	Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z
	;/notsolar  - yes, we'll take it in solar now

; trap for any NaN's
idx=where(finite(Lum_B) eq 0)
if idx(0) ne -1 then begin
	Lum_Bol(idx)= 100.0 & Lum_U(idx)= 100.0 & Lum_B(idx)= 100.0 & Lum_V(idx)= 100.0 & Lum_R(idx)= 100.0
	Lum_I(idx)= 100.0 & Lum_J(idx)= 100.0 & Lum_H(idx)= 100.0 & Lum_K(idx)= 100.0
	Lum_sdss_u(idx)= 100.0 & Lum_sdss_g(idx)= 100.0 & Lum_sdss_r(idx)= 100.0
	Lum_sdss_i(idx)= 100.0 & Lum_sdss_z(idx)= 100.0
	print, "Fixed ",n_elements(idx), " luminosities when trapping for NaN values."
endif



; set variables to send

M_tosend= m
; luminous material
lx= fload_allstars_xyz('x')
ly= fload_allstars_xyz('y')
lz= fload_allstars_xyz('z')

L_tosend= Lum_K
;L_tosend= Lum_B

print, "Total L_K= ", total(Lum_K)
print, "Total L_B= ", total(Lum_B)



; all mass
mx= fload_all_xyz('x')
my= fload_all_xyz('y')
mz= fload_all_xyz('z')
mm= 1.0e+10*fload_all_mass(1)


; Now do the hard work
; ----------------------
rs= fltarr(bins)
ML_avg= fltarr(bins)
ML_1sig= fltarr(bins)

;tempxmin= alog10(xmin)
;tempxmax= alog10(xmax)
tempxmin= xmin
tempxmax= xmax


process_M_to_L_formanyprojections, lx, ly, lz, l_tosend, $
				mx, my, mz, mm, $
				tempxmin, tempxmax, bins, $
				ML_avg, ML_1sig, rs    ; , /x_is_log



idx=where(ML_avg gt 0)
if idx(0) ne -1 then begin
	ML_avg=ML_avg(idx)
	ML_1sig=ML_1sig(idx)
	rs= rs(idx)
endif


; print it up
oplot, rs, ML_avg, psym=-3, thick=4.0, color= linecolor


; print error if we want
if keyword_set(includeerr) then begin
	oplot, rs, ML_p1sig, thick=2.0, psym=-3, color=linecolor, linestyle=1
	oplot, rs, ML_m1sig, thick=2.0, psym=-3, color=linecolor, linestyle=1
endif


if keyword_set(msg) then begin
   ; xyouts, 0.70, 0.85, msg, /normal, size= 1.5, color= linecolor, charthick=3.0
endif



end






