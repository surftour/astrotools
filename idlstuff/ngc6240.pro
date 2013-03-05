pro ngc6240, junk


frun= "/raid4/tcox/As/A3sauron"
snapnum= 90

;frun= "/raid4/tcox/As/A3"
;snapnum= 117

xlen= 10.0
;rotate_phi= 78.44
;rotate_theta= 180.0
;rotate_phi= 0.0
;rotate_theta= 0.0
;rotate_phi= 90.0
;rotate_theta= 90.0
rotate_phi= 45.0
rotate_theta= 0.0



	;ok=fload_snapshot(frun,snapnum)
	;ok=fload_snapshot_bh(frun,snapnum)
	ok=fload_snapshot_bh(frun,snapnum,/nopot_in_snap)


	; ----------------------------------------
	;  Do our thing
	; ----------------------------------------

	usethiscenter= [0,0,0]
	;usethiscenter= [13.0,56.0,-15.0]   ; /raid4/tcox/As/A3sauron, snap 90


	; gas quantities
	; --------------
	filename='gas_sd.eps'
	center= usethiscenter
	contour_gas, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, /loadedsnap, /colorbar
 
	filename='gas_vmap.eps'
	center= usethiscenter
	contour_gas_velocitymap, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, /loadedsnap, /colorbar

	filename='gas_smap.eps'
	center= usethiscenter
	contour_gas_dispersionmap, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, /loadedsnap, /colorbar


	; stellar quantities
	; ------------------
        filename='star_sd.eps'
	center= usethiscenter
        contour_allstars, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, /loadedsnap, /colorbar

        filename='star_vmap.eps'
	center= usethiscenter
        contour_allstars_velocitymap, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, /loadedsnap, /colorbar

        filename='star_smap.eps'
	center= usethiscenter
        contour_allstars_dispersionmap, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, /loadedsnap, /colorbar



print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end





