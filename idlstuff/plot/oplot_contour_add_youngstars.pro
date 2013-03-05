pro oplot_contour_add_youngstars, junk, xz=xz, yz=yz, $
			center=center, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			agecutoff=agecutoff


;----------------------------------------------------------------
;  
; 
;
;   Add new stars on top of contour image.
;
;
;
;
;
;
;
;----------------------------------------------------------------



if keyword_set(center) then begin
	orig_center= center
endif else begin
	orig_center=[0,0,0]
	center=[0,0,0]
endelse


        x= fload_newstars_xyz('x',center=center)
        y= fload_newstars_xyz('y',center=center)
        z= fload_newstars_xyz('z',center=center)
        m= fload_newstars_mass(1)


	; what time is it?
        Ttime= float(fload_time(1))
        nsage= float(Ttime-fload_newstars_age(1))

    
	; rotate
	; ---------
	if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin
              process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new
              x= x_new
              y= y_new
              z= z_new
	endif


	; -------------------------------
	;
        ;   Plot the actual points!!!
	;
        ; -------------------------------
        if keyword_set(yz) then begin
              x= y
              y= z
        endif

        if keyword_set(xz) then y= z

	print, "Nstars= ", n_elements(x)

        ; next youngest
        idx= where((nsage lt 0.035) and (nsage ge 0.010))
        if idx(0) ne -1 then begin
               print, "Nstars (younger between 0.035 and 0.010)= ", n_elements(x(idx))
               select_thispoint, 24, thispsym, thiscolor
               oplot, x(idx), y(idx), psym=thispsym, color=thiscolor
        endif

	; youngest stars
	idx= where(nsage lt 0.010)
	if idx(0) ne -1 then begin
		print, "Nstars (younger than 0.010)= ", n_elements(x(idx))
        	select_thispoint, 23, thispsym, thiscolor
        	oplot, x(idx), y(idx), psym=thispsym, color=thiscolor
	endif



; -------------
;  Done
; -------------



end


