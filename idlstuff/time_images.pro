pro time_images, frun, snapnum


if not keyword_set(frun) then begin
	print, " "
	print, " time_images, frun, snapnum"
	print, " "
	print, " "
	return
endif


; setup parameters
; -----------------------
xlen= 50.0

barlen= 10.0

;rotate_phi= 78.44
;rotate_theta= 180.0
rotate_phi= 0.0
rotate_theta= 0.0


; determine number of snapshots
; ------------------------------
spawn, "/usr/bin/ls "+frun+"/snap* | wc", result
nsnaps=long(result[0])


; cycle through the snapshot
; -----------------------------
for i=0,nsnaps do begin

	;ok=fload_snapshot(frun,i)
	ok=fload_snapshot_bh(frun,i)

	time= fload_time(1)

	; ----------------------------------------
	;  Do our thing
	; ----------------------------------------

	thisnum= '000'+strcompress(string(i),/remove_all)
	thisnum= strmid(thisnum,strlen(thisnum)-2,2)
        
	;center= get_center_to_use(time)
	usethiscenter= [0,0,0]

	; gas quantities
	; --------------
	filename=frun+'/g'+thisnum+'_sd.eps'
	center= usethiscenter
	contour_gas, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, $
		/loadedsnap, /colorbar
 
	;filename=frun+'/g'+thisnum+'_vmap.eps'
	;center= usethiscenter
	;contour_gas_velocitymap, frun, snapnum, xlen, 'ps', filename=filename, $
	;	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	;	center= center, $
	;	/loadedsnap, /colorbar

	;filename=frun+'/g'+thisnum+'_smap.eps'
	;center= usethiscenter
	;contour_gas_dispersionmap, frun, snapnum, xlen, 'ps', filename=filename, $
	;	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	;	center= center, $
	;	/loadedsnap, /colorbar


	; stellar quantities
	; ------------------
        filename=frun+'/s'+thisnum+'_sd.eps'
	center= usethiscenter
        contour_allstars, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, $
                /loadedsnap, /colorbar

        ;filename=frun+'/s'+thisnum+'_vmap.eps'
	;center= usethiscenter
        ;contour_allstars_velocitymap, frun, snapnum, xlen, 'ps', filename=filename, $
	;	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	;	center= center, $
        ;        /loadedsnap, /colorbar

        ;filename=frun+'/s'+thisnum+'_smap.eps'
	;center= usethiscenter
        ;contour_allstars_dispersionmap, frun, snapnum, xlen, 'ps', filename=filename, $
	;	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	;	center= center, $
        ;        /loadedsnap, /colorbar

	;snapnum= snapnum+5


endfor


print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end










function get_center_touse, time



        if i eq 0 then usethiscenter= [89.0022,45.5882,0.00309830]
        if i eq 1 then usethiscenter= [54.0885,36.3291,0.312841]
        if i eq 2 then usethiscenter= [11.4516,20.4729,0.526275]
        if i eq 3 then usethiscenter= [-12.4595,-20.8376,0.621572]
        if i eq 4 then usethiscenter= [-13.2122,-36.5850,-0.0530290]
        if i eq 5 then usethiscenter= [-10.1949,-35.8209,-0.991644]
        if i eq 6 then usethiscenter= [-5.09338,-23.4534,-1.80362]
        if i eq 7 then usethiscenter= [0.849999,-0.0348494,-2.05953]
        if i eq 8 then usethiscenter= [1.80098,-0.558575,-2.63818]
        if i eq 9 then usethiscenter= [0,0,0]
        if i eq 10 then usethiscenter= [0,0,0]
        if i eq 11 then usethiscenter= [0,0,0]
        if i eq 12 then usethiscenter= [0,0,0]

	return, usethiscenter

end
