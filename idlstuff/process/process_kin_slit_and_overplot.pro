

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Calculate  Velocity & Sigma along Slit
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; returns v_mean, v_disp

pro process_slit, x, y, v, mass, slit, v_mean, v_disp, $
			slit_width=slit_width, $
			slit_len=slit_len, $
			slit_bins=slit_bins, $
			v_mean_unweighted=v_mean_unweighted, $
			v_disp_unweighted=v_disp_unweighted, $
			mags=mags, $
			v_mean_magwt=v_mean_magwt, $
			v_disp_magwt=v_disp_magwt, $
			slit_mass= slit_mass

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
;	x_inslit= [14000]
endelse

v_mean  = fltarr(slit_bins)
v_mean_unweighted= fltarr(slit_bins)
v_mean_magwt= fltarr(slit_bins)
v_disp = fltarr(slit_bins)
v_disp_unweighted= fltarr(slit_bins)
v_disp_magwt= fltarr(slit_bins)
slit_mass= fltarr(slit_bins)


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
		if n_elements(idx) gt 1 then v_disp_unweighted(i)= sqrt(variance(v_inslit(idx))) else v_disp_unweighted(i)= 0.8 * v_mean(i)
		if keyword_set(mags) then v_disp_magwt(i) = sqrt(total(mags_inslit(idx)*(v_inslit(idx)-v_mean_magwt(i))^2)/total(mags_inslit(idx)))

		; total mass in this slit element
		slit_mass(i)= total(mass_inslit(idx))

;print, "Slit x= ",slit(i), v_mean(i), mean(v_inslit(idx)), v_disp(i), sqrt(variance(v_inslit(idx)))
	endif else begin
		v_mean(i) = 0.0
		v_mean_unweighted(i) = 0.0
		v_mean_magwt(i) = 0.0
		v_disp(i) = 0.0
		v_disp_unweighted(i) = 0.0
		v_disp_magwt(i) = 0.0
		slit_mass(i)= 1.0e-10
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




;----------------------------------------------



pro write_slit_info, fitsdir, istring, rad, v, sig, $
				major=major, minor=minor

;
;    write information to file
;
if keyword_set(major) then filename= fitsdir+'/'+istring+'_lst_slit_major.txt'
if keyword_set(minor) then filename= fitsdir+'/'+istring+'_lst_slit_minor.txt'
openw, 1, filename, ERROR=err

printf, 1, "#   ###_lst_slit_[axis].txt, "+fitsdir
printf, 1, "#   "
printf, 1, "#                             "
printf, 1, "#  radius      v       sigma  "
printf, 1, "# (kpc/h)    (km/s)    (km/s) "


nbins= n_elements(rad)

for i=0,nbins-1 do begin
        printf, 1, FORMAT= '(" ",3(F8.3,"  "))', $
                rad[i], v[i], sig[i]
endfor
close, 1


end



;========================================================================
;========================================================================
;
;
;
;   Main Procedure
;
;
;========================================================================
;========================================================================


; process the slit, and overplot it
; -----------------------------------
pro process_kin_slit_and_overplot, x, y, z, vx, vy, vz, mass, $
				x0, y0, x1, y1, $
				center=center, labelit=labelit, fitsdir=fitsdir, istring=istring, $
				rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
				slit_len=slit_len, slit_width=slit_width, $
				slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
				show_slitlen=show_slitlen, $
				pa=pa, mags=mags, overplot=overplot, $
				Rin=Rin, Sig= Sig, R_e=R_e, $
				Vrot_maj=Vrot_maj, avgVrot_maj=avgVrot_maj, $
				Vrot_min=Vrot_min, avgVrot_min=avgVrot_min, $
				ditchsigma=ditchsigma, thiscolor= thiscolor, $
				outerkinematics=outerkinematics, $
				innerkinematics=innerkinematics, $
				masswtkinematics=masswtkinematics



print, "====================================="
print, "    process_kin_slit_and_overplot    "
print, "====================================="




if keyword_set(center) then begin
     print, "process slit, changing center to= ",center
     x=x-center[0]
     y=y-center[1]
     z=z-center[2]

endif


; Now we need to rotate for the slit
; ------------------------------------
;if (rotate_theta gt 0) or (rotate_phi gt 0) or (total(abs(center)) gt 0) then begin
if keyword_set(rotate_theta) or keyword_set(rotate_phi) then begin

        ; rotate things
        process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new
	process_rotation, vx, vy, vz, rotate_theta, rotate_phi, vx_new, vy_new, vz_new

endif else begin
	x_new= x
	y_new= y
	z_new= z
	vx_new= vx
	vy_new= vy
	vz_new= vz
endelse


; ------------------
f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill
if not keyword_set(thiscolor) then thiscolor= 150
if thiscolor eq 150 then thispsym= 8
if thiscolor eq 100 then thispsym= 2
if thiscolor eq 75 then thispsym= 7
if thiscolor eq 50 then thispsym= 5
if thiscolor eq 5 then thispsym= 3






; plot it up
; -----------
;xaxistitle='!6R (h!u-1!N kpc)'
xaxistitle='!6R (kpc)'   ;  &   slitr= slitr / 0.7  & slitmax = slitmax / 0.7
;yaxistitle='!6V(R) (km s!U-1!N)'
yaxistitle='!6Velocity (km s!U-1!N)'
;yaxistitle='!6(km s!U-1!N)'
;yaxistitle='!Ms!N(R) (km s!U-1!N)'
xmax = slitmax
xmin = -slitmax
ymax = vmax
ymin = -vmax

;xmax= 8.2  & xmin= -8.2   ; to match the 0.7 factor above

;---------------------------

if not keyword_set(overplot) then begin
	;!p.position= [0.18, 0.15, 0.95, 0.95]
	!p.position= [x0, y0, x1, y1]
	!p.ticklen=0.03

	plot, [1.0],[1.0], psym=-3, $
	        xrange=[xmin,xmax], yrange=[ymin,ymax], $
	        color= 0, $
	        ;/ylog, $
	        ;/xlog, $
	        xstyle=1, $
	        ystyle=1, $
	        ;ystyle=8, $     ; this will suppress the second y axis
	        ;xcharsize=1.10, $
		;ycharsize=1.10, $
	        xcharsize=1.5, $
		ycharsize=1.5, $
	        xthick=4.0, ythick=4.0, $
	        charthick=3.0, $
	        ;ytickformat='exp_label', $
		xtickformat='(a1)', $
		;ytickformat='(a1)', $
	        ;xtitle=xaxistitle, $
	        ytitle=yaxistitle, $
	        /nodata, /noerase

endif




;
;   Major Axis 
;------------------

; rotate just x and y so that it's
; along the major axis position angle (pa)
; ----------------------------------------
if keyword_set(pa) then begin
	print, "major axis"
        print, "actual slit pa= ", pa
        par = pa * !PI/180.0
        x_newest= x_new*cos(par) + y_new*sin(par)      ; with this choice of signs, it automatically
        y_newest= -x_new*sin(par) + y_new*cos(par)     ;  rotates clock-wise (hence don't need -par)
        x_new= x_newest
        y_new= y_newest
endif else begin
        print, "pa not set.  fixing to be 0"
        pa= 0.0
endelse


; calculate slit velocities
process_slit, x_new, y_new, vz_new, mass, slit, v_mean, v_disp, $
                        slit_width=slit_width, $
                        slit_len=slit_len, $
                        slit_bins=slit_bins, $
                        v_mean_unweighted=v_mean_unweighted, $
                        v_disp_unweighted=v_disp_unweighted, $
                        mags=mags, $
                        v_mean_magwt=v_mean_magwt, $
                        v_disp_magwt=v_disp_magwt, $
			slit_mass=slit_mass




; -----------------------------------------------
; now, process some of the slit information
;


;
;  full data
;
slitr= slit
slitv= v_mean
slitsigma= v_disp
if keyword_set(masswtkinematics) then begin
        slitv(*)= v_mean(*) * slit_mass(*) / total(slit_mass)                         ; mass weighted
	slitsigma(*)= v_disp(*) * slit_mass(*) / total(slit_mass)
;stop
endif

print, "slit velocity max/min= ",max(slitv),min(slitv)
print, "slit sigma max/min= ",max(slitsigma),min(slitsigma)
oplot,slitr,slitv,psym=-thispsym,color=thiscolor, thick=1.0, symsize=1.0


if keyword_set(fitsdir) then write_slit_info, fitsdir, istring, slitr, slitv, slitsigma, /major


; avg. maximum velocities (+/-) along slit
posidx= where(slitr gt 0)
negidx= where(slitr lt 0)
if keyword_set(outerkinematics) then begin
	posidx= where(slitr gt +2.5*Rin)     ; Rin = 0.5 R_e
	negidx= where(slitr lt -2.5*Rin)
endif
if keyword_set(innerkinematics) then begin
	posidx= where((slitr lt +2.5*Rin) and (slitr ge 0.0))     ; Rin = 0.5 R_e
	negidx= where((slitr gt -2.5*Rin) and (slitr lt 0.0))
endif
Vrot_maj = 0.5 * (max(abs(slitv(posidx))) + max(abs(slitv(negidx))))
print, "Vrot_maj= ", Vrot_maj

;idx= where(abs(slitv) gt 0.5*vrot)
idx= where(abs(slitr) gt Rin)
if keyword_set(outerkinematics) then idx= where(abs(slitr) gt 2.5*Rin)
if keyword_set(innerkinematics) then idx= where(abs(slitr) lt 2.5*Rin)
avgVrot_maj = mean(abs(slitv(idx)))
print, "avgVrot_maj= ", avgVrot_maj



;
;  smoothed data
;
print, " *** smoothed version *** "
slitr= smooth(slitr,2)
slitv= smooth(slitv,2)
slitsigma= smooth(slitsigma,2)
print, "slit velocity max/min= ",max(slitv),min(slitv)
print, "slit sigma max/min= ",max(slitsigma),min(slitsigma)
;oplot,slitr,slitv,psym=-thispsym,color=thiscolor, thick=5.0, symsize=1.5
oplot,slitr,slitv,psym=-3,color=thiscolor, thick=6.0
;oplot,slitr,slitv,psym=-thispsym,color=thiscolor, thick=2.0, symsize=1.00
;oplot,slitr,v_mean_unweighted,psym=-3,linestyle= 1, color=thiscolor

; process of velocity information
posidx= where(slitr gt 0)
negidx= where(slitr lt 0)
if keyword_set(outerkinematics) then begin
	posidx= where(slitr gt +2.5*Rin)     ; Rin = 0.5 R_e
	negidx= where(slitr lt -2.5*Rin)
endif
if keyword_set(innerkinematics) then begin
	posidx= where((slitr lt +2.5*Rin) and (slitr ge 0.0))     ; Rin = 0.5 R_e
	negidx= where((slitr gt -2.5*Rin) and (slitr lt 0.0))
endif
Vrot_maj = 0.5 * (max(abs(slitv(posidx))) + max(abs(slitv(negidx))))
print, "Vrot_maj= ", Vrot_maj

;idx= where(abs(slitv) gt 0.5*vrot)
idx= where(abs(slitr) gt Rin)
if keyword_set(outerkinematics) then idx= where(abs(slitr) gt 2.5*Rin)
if keyword_set(innerkinematics) then idx= where(abs(slitr) lt 2.5*Rin)
avgVrot_maj = mean(abs(slitv(idx)))
print, "avgVrot_maj= ", avgVrot_maj


idx= where(abs(slitr) lt Rin) 
Sig= mean(slitsigma(idx))


; add label
xyouts, x0+0.16, y0+0.12, 'major', color=thiscolor, charthick=3.0, size=1.33, /normal




; -----------------------------------------------


;
;   Minor Axis 
;------------------

; rotate just x and y so that it's
; along the major axis position angle (pa)
; ----------------------------------------
;if keyword_set(pa) then begin
;        print, "now, minor axis "
;	if pa lt 270.0 then par = pa + 90.0 else par = pa - 90.0
;	print, "minor axis pa = ", par
;        par = par * !PI/180.0
;        x_newest= x_new*cos(par) + y_new*sin(par)      ; with this choice of signs, it automatically
;        y_newest= -x_new*sin(par) + y_new*cos(par)     ;  rotates clock-wise (hence don't need -par)
;        x_new= x_newest
;        y_new= y_newest
;endif else begin
;        print, "pa not set.  fixing to be 90 for minor axis"
;        pa= 90.0
;endelse


; don't need to rotate since previous already did this,
; instead we'll just switch x and y

; calculate slit velocities
;process_slit, x_new, y_new, vz_new, mass, slit, v_mean, v_disp, $
process_slit, y_new, x_new, vz_new, mass, slit, v_mean, v_disp, $
                        slit_width=slit_width, $
                        slit_len=slit_len, $
                        slit_bins=slit_bins, $
                        v_mean_unweighted=v_mean_unweighted, $
                        v_disp_unweighted=v_disp_unweighted, $
                        mags=mags, $
                        v_mean_magwt=v_mean_magwt, $
                        v_disp_magwt=v_disp_magwt, $
			slit_mass=slit_mass


; -----------------------------------------------
; now, process some of the slit information
;

;
;  full data
;
slitr= slit
slitv= v_mean
slitsigma= v_disp
if keyword_set(masswtkinematics) then begin
        slitv(*)= v_mean(*) * slit_mass(*) / total(slit_mass)                         ; mass weighted
	slitsigma(*)= v_disp(*) * slit_mass(*) / total(slit_mass)
endif

print, "slit velocity max/min= ",max(slitv),min(slitv)
print, "slit sigma max/min= ",max(slitsigma),min(slitsigma)
thispsym= 2
thiscolor= thiscolor-100
oplot,slitr,slitv,psym=-thispsym,color=thiscolor, thick=1.0, symsize=1.0, linestyle= 1


if keyword_set(fitsdir) then write_slit_info, fitsdir, istring, slitr, slitv, slitsigma, /minor


; process of velocity information
posidx= where(slitr gt 0)
negidx= where(slitr lt 0)
if keyword_set(outerkinematics) then begin
	posidx= where(slitr gt +2.5*Rin)     ; Rin = 0.5 R_e
	negidx= where(slitr lt -2.5*Rin)
endif
if keyword_set(innerkinematics) then begin
	posidx= where((slitr lt +2.5*Rin) and (slitr ge 0.0))     ; Rin = 0.5 R_e
	negidx= where((slitr gt -2.5*Rin) and (slitr lt 0.0))
endif
Vrot_min = 0.5 * (max(abs(slitv(posidx))) + max(abs(slitv(negidx))))
print, "Vrot_min= ", Vrot_min

;idx= where(abs(slitv) gt 0.5*vrot)
idx= where(abs(slitr) gt Rin)
if keyword_set(outerkinematics) then idx= where(abs(slitr) gt 2.5*Rin)
if keyword_set(innerkinematics) then idx= where(abs(slitr) lt 2.5*Rin)
avgVrot_min = mean(abs(slitv(idx)))
print, "avgVrot_min= ", avgVrot_min



;
;  smoothed data
;
print, " *** smoothed version *** "
slitr= smooth(slitr,2)
slitv= smooth(slitv,2)
slitsigma= smooth(slitsigma,2)
print, "slit velocity max/min= ",max(slitv),min(slitv)
print, "slit sigma max/min= ",max(slitsigma),min(slitsigma)

;oplot,slitr,slitv,psym=-thispsym,color=thiscolor, thick=1.0, symsize=1.0, linestyle= 1
oplot,slitr,slitv,psym=-3,color=thiscolor, thick=6.0
;oplot,slitr,slitv,psym=-thispsym,color=thiscolor, thick=2.0, symsize=1.00
;oplot,slitr,v_mean_unweighted,psym=-3,linestyle= 1, color=thiscolor

; process of velocity information
posidx= where(slitr gt 0)
negidx= where(slitr lt 0)
if keyword_set(outerkinematics) then begin
	posidx= where(slitr gt +2.5*Rin)     ; Rin = 0.5 R_e
	negidx= where(slitr lt -2.5*Rin)
endif
if keyword_set(innerkinematics) then begin
	posidx= where((slitr lt +2.5*Rin) and (slitr ge 0.0))     ; Rin = 0.5 R_e
	negidx= where((slitr gt -2.5*Rin) and (slitr lt 0.0))
endif
Vrot_min = 0.5 * (max(abs(slitv(posidx))) + max(abs(slitv(negidx))))
print, "Vrot_min= ", Vrot_min

;idx= where(abs(slitv) gt 0.5*vrot)
idx= where(abs(slitr) gt Rin)
if keyword_set(outerkinematics) then idx= where(abs(slitr) gt 2.5*Rin)
if keyword_set(innerkinematics) then idx= where(abs(slitr) lt 2.5*Rin)
avgVrot_min = mean(abs(slitv(idx)))
print, "avgVrot_min= ", avgVrot_min


; add label
xyouts, x0+0.16, y0+0.04, 'minor', color=thiscolor, charthick=3.0, size=1.33, /normal





;-----------------------------------------------



if keyword_set(show_slitlen) then begin
	slitlenlbl= strcompress(string(2.0*slit_len),/remove_all)
	if (2.0*slit_len) ge 1.0 then digs= 1
        if (2.0*slit_len) ge 10.0 then digs= 2
        if (2.0*slit_len) ge 100.0 then digs= 3
        slitlenlbl = strmid(slitlenlbl,0,digs)        ; T=0.x (4+digits after decimal)
        slitlenlbl= 'slit= '+slitlenlbl+' kpc/h'
        xyouts, x1-0.10, y1-0.07, slitlenlbl, /normal, size= 1.2, charthick=3.0, color= 0
endif


if keyword_set(labelit) then begin
   ;yval= 0.65
   ;ylbl= '!6velocity'
   yval= 0.50
   ylbl= '!6all stars'
   if thiscolor eq 100 then yval= 0.59
   if thiscolor eq 100 then ylbl= 'young'
   if thiscolor eq 75 then yval= 0.68
   if thiscolor eq 75 then ylbl= 'intermediate'
   if thiscolor eq 50 then yval= 0.77
   if thiscolor eq 50 then ylbl= 'old'
   if thiscolor eq 5 then yval= 0.86
   if thiscolor eq 5 then ylbl= 'dissipationless'
   oplot, [slitr[0],slitr[2]], [ymin*yval,ymin*yval], psym=-thispsym,color=thiscolor, symsize= 1.5, thick= 5.0
   ;oplot, [slitr[0],slitr[2]], [ymin*yval,ymin*yval], psym=-thispsym,color=thiscolor, symsize= 1.00, thick= 2.0
   xyouts, slitr[4], ymin*(yval+0.02), ylbl, color=thiscolor, charthick=3.0, size=1.33
endif


if keyword_set(slitsigma) and not keyword_set(ditchsigma) then begin
	thispsym= 2
	thiscolor= 50
	;oplot, slitr, slitsigma, psym=-thispsym, color=thiscolor, symsize= 1.5, thick= 5.0
	oplot, slitr, slitsigma, psym=-thispsym, color=thiscolor, symsize= 1.33, thick= 2.0
	;oplot, slitr, v_disp_unweighted, psym=-3, linestyle=1, color=thiscolor

	if keyword_set(mags) then begin
	   oplot, slitr, v_disp_magwt, psym=-3, linestyle=2, color=thiscolor, thick=3.0
	endif

	print, "slit dispersion max/min= ", max(slitsigma), min(slitsigma)

	if keyword_set(labelit) then begin
	   oplot, [slitr[0],slitr[2]], [ymin*0.82,ymin*0.82], psym=-thispsym,color=thiscolor
	   xyouts, slitr[4], ymin*0.85, '!7r!6', color=thiscolor, charthick=3.0, size=1.33
	endif
endif


; ------------------
;  Print any extras
; ------------------
;xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0
;xyouts, 0.2, 0.90, fload_fid(frun), /normal, charthick=2.0, size=1.33, color=0


if not keyword_set(overplot) then begin
	; zero
	x=[xmin,xmax]
	y=[0,0]
	oplot, x, y, psym=-3, linestyle=2, color= 0, thick=2.0
endif



end







