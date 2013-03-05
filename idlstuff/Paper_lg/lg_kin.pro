;---------------------------------------------------------------
;
;
;
;---------------------------------------------------------------
pro lgkin, junk


frun= '/raid4/tcox/localgroup/v5'
snapnum= 105

xlen= 10.0


; slit dimensions
; -----------------
slit_width = 0.5
slit_bins = 26
slit_len = 5.0

; plot range
slitmax = 1.2*slit_len


set_maxden= 1.0e+0
set_dynrng= 1.0e+4


; open the file up
;   and load variables
; ---------------------
ok= fload_snapshot_bh(frun, snapnum)

x_orig= fload_allstars_xyz('x')
y_orig= fload_allstars_xyz('y')
z_orig= fload_allstars_xyz('z')
vx_orig= fload_allstars_v('x')
vy_orig= fload_allstars_v('y')
vz_orig= fload_allstars_v('z')
m_orig= fload_allstars_mass(1)



; --------------------------------------------------------

;  read the list of theta's and phi's

anglefile= '/home/tcox/unitsphereangles.txt'
;anglefile= '/home/tcox/unitsphereangles.txt.test'

spawn, "wc "+anglefile,result
lines=long(result)
Nangles=lines(0)-5
if Nangles GT 0 then angle_data= fltarr(2,Nangles)

openr, 1, anglefile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, angle_data
close, 1


theta= angle_data[0,*]
phi= angle_data[1,*]


; --------------------------------------------------------
;
;   Loop through the angles, project and calculate
;  the disky/boxy - ness
;

vmax= fltarr(Nangles)
sig= fltarr(Nangles)
pai= fltarr(Nangles)
ellipi= fltarr(Nangles)



for i= 0, Nangles-1 do begin

	rotate_phi= phi[i]
	rotate_theta= theta[i]

	if rotate_phi lt 1.0e-6 then rotate_phi= 1.0e-6
	if rotate_theta lt 1.0e-6 then rotate_theta= 1.0e-6

	; ------------------------------------------------------------------


	x= x_orig
	y= y_orig
	z= z_orig
	m= m_orig
        contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, $
                        filename=filename, fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, $
                        pixels=pixels, zthickness=zthickness, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        crude=crude, center=center, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage, HalfMassSB=HalfMassSB

        Pic= NxNImage
	hmsb= HalfMassSB


        
        istring= strcompress(string(i), /remove_all)
        filename= frun+'/eps/image'+istring+'.eps'    ; kind of funky becuase
                                                        ; contour_makepic assumes
                                                        ; the ending is eps


	initialize_plotinfo, 1
	setup_plot_stuff, 'ps', filename=filename, colortable= 0, newxsize=12, newysize=12

        x0= 0.01
        x1= 0.99
        y0= 0.01
        y1= 0.99

        tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

	!p.position=[x0,y0,x1,y1]
        !p.ticklen=0.03

	plot,[0],[0], psym=3, xstyle=1, ystyle=1, color= 0, /noerase, $
		xrange=[-xlen,xlen], yrange=[-xlen,xlen], $
		xcharsize=0.01, ycharsize=0.01, $
		xthick=4.0, ythick=4.0, /normal, /nodata, $
		xtickformat='(a1)', ytickformat='(a1)'

	xyouts, 0.07, 0.09, fload_timelbl(1,2,/noteq), /normal, size= 1.3, charthick=3.0, color= 0



	; ------------------------------------------------------------------
	;
	; this overplots contours and fits an ellipse to
	; the half-mass countour
	;
	oplot_contours, Pic, x0, y0, x1, y1, hmsb=hmsb, $
                                xlen=xlen, /fitellipse, $
                                pa=pa, ellip=ellip



	device, /close

	; ------------------------------------------------------------------


	; find semi-major axis & rotate
	; ----------------------------------------
	print, "passed pa= ", pa
	; there is a difference in definition, pa as above is
	; clockwise rotation, but below is opposite
	if pa gt 0.0 then pa= 180.0 - pa else pa= -pa
       	print, "new pa= ", pa
       	par = pa * !PI/180.0


	ellipi[i]= ellip
	pai[i]= pa

	x= x_orig
	y= y_orig
	z= z_orig
	vx= vx_orig
	vy= vy_orig
	vz= vz_orig
	mass=m_orig
        process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new
        process_rotation, vx, vy, vz, rotate_theta, rotate_phi, vx_new, vy_new, vz_new

       	x_newest= x_new*cos(par) + y_new*sin(par)      ; with this choice of signs, it automatically
       	y_newest= -x_new*sin(par) + y_new*cos(par)     ;  rotates clock-wise (hence don't need -par)
       	x_new= x_newest
       	y_new= y_newest



	; calculate slit velocities
	process_slit, x_new, y_new, vz_new, mass, slit, v_mean, v_disp, $
                        slit_width=slit_width, $
                        slit_len=slit_len, $
                        slit_bins=slit_bins, $
                        v_mean_unweighted=v_mean_unweighted, $
                        v_disp_unweighted=v_disp_unweighted, $
                        mags=mags, $
                        v_mean_magwt=v_mean_magwt, $
                        v_disp_magwt=v_disp_magwt



	; now, calculate v_max and sigma


	idx= where(slit le 1.0)
	sigggg= mean(v_disp(idx))
	sig[i]= sigggg

	v1= abs(max(v_mean))
	v2= abs(min(v_mean))
	vmaxxxx= 0.5 * (v1 + v2)
	vmax[i]= vmaxxxx

	; ------------------------------------------------------------------

endfor


; ------------------------------------------------------------------

openw, 1, frun+'/kinematics.txt', ERROR=err
        
printf, 1, "#   kinematics.txt"
printf, 1, "#   "
printf, 1, "#   "
printf, 1, "#  theta      phi     ellip       pa      vmax       sig  "
printf, 1, "#  (deg)     (deg)   (1-b/a)    (deg)    (km/s)    (km/s) "

for i=0,Nangles-1 do begin
	printf, 1, FORMAT= '(2(F8.3,"  "),F8.5,"  ",F8.3,"  ",2(F8.2,"  "))', $
                theta[i], phi[i], ellipi[i], pai[i], vmax[i], sig[i]
endfor
close, 1




; ------------------------------------------------------------------



end









; =================================================================
; =================================================================




pro read_kin_file, frun, theta, phi, ellipi, pai, vmax, sig

spawn, '/bin/ls '+frun+'/kinematics.txt ',result

kinfile=strcompress(result[0],/remove_all)

filename= kinfile

print, " --------------------------- "
print, 'opening: ',filename

spawn, "wc "+filename,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then kin_data= fltarr(6,lines)

openr, 1, filename

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, kin_data
close, 1


theta=  kin_data[0,*]
phi=    kin_data[1,*]
ellipi= kin_data[2,*]
pai=    kin_data[3,*]
vmax=   kin_data[4,*]
sig=    kin_data[5,*]

print, "<e>=",mean(ellipi), "  +/-", sqrt(variance(ellipi))
print, "<v>=",mean(vmax), "  +/-", sqrt(variance(vmax))
print, "<sigma>=",mean(sig), "  +/-", sqrt(variance(sig))


end




; =================================================================
; =================================================================





;
; returns v_mean, v_disp
; -------------------------
pro process_slit, x, y, v, mass, slit, v_mean, v_disp, $
                        slit_width=slit_width, $
                        slit_len=slit_len, $
                        slit_bins=slit_bins, $
                        v_mean_unweighted=v_mean_unweighted, $
                        v_disp_unweighted=v_disp_unweighted, $
                        mags=mags, $
                        v_mean_magwt=v_mean_magwt, $
                        v_disp_magwt=v_disp_magwt

slit  = fltarr(slit_bins)
slit(*) = (2*slit_len)*findgen(slit_bins)/float(slit_bins) - slit_len
dx = (2*slit_len)/float(slit_bins)

; trim particles
idx= where((y le slit_width/2.0) and (y ge -1.0*slit_width/2.0))
if idx(0) ne -1 then begin
        x_inslit= x(idx)
        mass_inslit= mass(idx)
        v_inslit= v(idx)
        if keyword_set(mags) then mags_inslit= mags(idx)
endif else begin
;       x_inslit= [14000]
endelse

v_mean  = fltarr(slit_bins)
v_mean_unweighted= fltarr(slit_bins)
v_mean_magwt= fltarr(slit_bins)
v_disp = fltarr(slit_bins)
v_disp_unweighted= fltarr(slit_bins)
v_disp_magwt= fltarr(slit_bins)



for i=0,slit_bins-1 do begin

        ;idx = where((x_inslit gt slit(i)-dx) and (x_inslit le slit(i)+dx))
        idx = where((x_inslit gt slit(i)-0.5*dx) and (x_inslit le slit(i)+0.5*dx))
        if idx(0) ne -1 then begin

                ; velocity
                v_mean(i)  = total(mass_inslit(idx)*v_inslit(idx))/total(mass_inslit(idx))
                v_mean_unweighted(i)= mean(v_inslit(idx))
                if keyword_set(mags) then v_mean_magwt(i)  = total(mags_inslit(idx)*v_inslit(idx))/total(mags_inslit(idx))

                ; dispersion
                v_disp(i) = sqrt(total(mass_inslit(idx)*(v_inslit(idx)-v_mean(i))^2)/total(mass_inslit(idx)))
                v_disp_unweighted(i)= sqrt(variance(v_inslit(idx)))
                if keyword_set(mags) then v_disp_magwt(i) = sqrt(total(mags_inslit(idx)*(v_inslit(idx)-v_mean_magwt(i))^2)/total(mags_inslit(idx)))

;print, "Slit x= ",slit(i), v_mean(i), mean(v_inslit(idx)), v_disp(i), sqrt(variance(v_inslit(idx)))
        endif else begin
                v_mean(i) = 0.0
                v_disp(i) = 0.0
        endelse

endfor

; average of central three points
; -------------------------------
c_i= slit_bins/2
central_avg= 0.333*(v_disp[c_i-1]+v_disp[c_i]+v_disp[c_i+1])
print, "Central Velocity Disperion: ", central_avg


; rotation velocity
; ------------------
;f = (2.0*!pi/16.0)*findgen(17)
;usersym,cos(f),sin(f),/fill
;oplot,slit,v_rot,psym=-8,color=140


; dispersion
; ------------
;oplot,slit,v_disp,psym=-2,color=40


end






; =================================================================
; =================================================================







;------------------------------------------------------------------------
;
;     Density Map of v/sigma vs. Ellipticity
;     ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro vsig_map, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "vsig_map, junk"
   print, "  "
   return 
endif

filename='lgvsig.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3


;--------------------------------------
;--------------------------------------


; ellipticity
;xaxistitle= 'ellipticity'
xaxistitle= '!7e!6'
xmax= 0.68
xmin= 0.0

; v/sigma
;yaxistitle= 'v/sigma'
;yaxistitle= '!7t!3/!7r!3'
yaxistitle= '!6V!Dmaj!N/!7r!6'
;yaxistitle= 'V!Dmin!N/!7r!3'
ymax = 1.7
ymin = 0.0


bins= 30



read_kin_file, "/raid4/tcox/localgroup/v5", theta, phi, ellipi, pai, vmax, sig

vsig= vmax/sig
ell= ellipi




; ------------------
; compute histogram
; ------------------

contour_makegeneralpic, ell, vsig, xmax, xmin, ymax, ymin, $
                                pixels= bins, $
                                NxNImage=NxNImage



;--------------------------------------- 
;  Print Tv image
;-----------------------------------------

x0= 0.20
y0= 0.15
x1= 0.98
y1= 0.98
x_size= (x1-x0)
y_size= (y1-y0)
           
                
tv, NxNImage, x0, y0, xsize=x_size, ysize=y_size, /normal



;----------------------------------------
; Generate plot 
;----------------------------------------

!p.position=[x0, y0, x0+x_size,y0+y_size]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, $
        xrange=[xmin,xmax],$
        yrange=[ymin,ymax],$
        color=0, $
        xcharsize=1.5, ycharsize=1.5, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;xticks= 14, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /noerase, $
        /nodata





; --------------------

; line of isotropic rotator
x=findgen(50)*(0.6/50.)
y=sqrt(x/(1-x))
oplot, x, y, psym=-3, color=0, thick=3.0


; --------------------

; these scripts are in:    Paper_RemI/visg_general
;
;
; add some real data!
;readandplot_bender_data, 1
;readandplot_davies_data, 1
;readandplot_deZeeuw_data, 1

readandplot_bender_data, 1, /addkey
readandplot_davies_data, 1, /addkey
readandplot_deZeeuw_data, 1, /addkey




; --------------------

xyouts, 0.8, 0.9, '(c)', /normal, size=1.8, color= 0, charthick=3.0

; --------------------

device, /close



end


