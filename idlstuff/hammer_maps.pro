pro hammer_maps, frun, snapeveryX, startid=startid, numpart=numpart, $
			rotatedangle=rotatedangle


;frun= "execute/G3gf1G3gf1b-u1"
;frun= "execute/G3gf2G3gf2b-u1"

;snapnum= 200
snapnum= 170
;snapnum= 30
;snapnum= 0
;xlen= 50.0

hdirg= '/hammer/top_gas'
hdirs= '/hammer/top_st'
if keyword_set(rotatedangle) then begin
   hdirg= '/hammer/tilt_gas'
   hdirs= '/hammer/tilt_st'
endif

rotate_phi= 0.0
rotate_theta= 0.0
if keyword_set(rotatedangle) then begin
   rotate_phi= 78.44
   rotate_theta= 180.0
   ;rotate_phi= 90.0
   ;rotate_theta= 90.0
endif


;for i=0,n_elements(frun)-1 do begin
;for i=0,12 do begin
;for i=11,12 do begin
;for i=1,20 do begin
for i=0,16 do begin
;for i=1,6 do begin

	print, "--------------------------------------"
	if fload_snapshot(frun,snapnum) then begin
	;if fload_snapshot_bh(frun,snapnum) then begin
		print, "PROBLEM: opening file",frun,snapnum
		;return
	endif else begin

	; ----------------------------------------
	;  Do our thing
	; ----------------------------------------

	print, "Time= ", fload_time(1)

	thisnum= '000'+strcompress(string(snapnum),/remove_all)
	thisnum= strmid(thisnum,strlen(thisnum)-3,3)
        
	usethiscenter= [0,0,0]
	usethiscenter= fload_1gal_center(startid,numpart)
	print, "center= ", usethiscenter

	;if i eq 0 then usethiscenter= [89.0022,45.5882,0.00309830]
	;if i eq 1 then usethiscenter= [54.0885,36.3291,0.312841]
	;if i eq 2 then usethiscenter= [11.4516,20.4729,0.526275]
	;if i eq 3 then usethiscenter= [-12.4595,-20.8376,0.621572]
	;if i eq 4 then usethiscenter= [-13.2122,-36.5850,-0.0530290]
	;if i eq 5 then usethiscenter= [-10.1949,-35.8209,-0.991644]
	;if i eq 6 then usethiscenter= [-5.09338,-23.4534,-1.80362]
	;if i eq 7 then usethiscenter= [0.849999,-0.0348494,-2.05953]
	;if i eq 8 then usethiscenter= [1.80098,-0.558575,-2.63818]
	;if i eq 9 then usethiscenter= [0,0,0]
	;if i eq 10 then usethiscenter= [0,0,0]
	;if i eq 11 then usethiscenter= [0,0,0]
	;if i eq 12 then usethiscenter= [0,0,0]

	;if i eq 11 then usethiscenter= [5.0,40.0,-23.0]   ; Sbc201b_wBH, snap 30
	;if i eq 11 then usethiscenter= [13.0,56.0,-15.0]   ; Sbc201b, snap 30

	; gas quantities
	; --------------
	filename=frun+hdirg+thisnum+'_sd.eps'
	center= usethiscenter
	contour_gas, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, $
		/loadedsnap, /fitstoo, /colorbar
 
	filename=frun+hdirg+thisnum+'_vmap.eps'
	center= usethiscenter
	contour_gas_velocitymap, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, $
		/loadedsnap, /fitstoo, /colorbar

	filename=frun+hdirg+thisnum+'_smap.eps'
	center= usethiscenter
	contour_gas_dispersionmap, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, $
		/loadedsnap, /fitstoo, /colorbar


	; stellar quantities
	; ------------------
        filename=frun+hdirs+thisnum+'_sd.eps'
	center= usethiscenter
        contour_allstars, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, $
                /loadedsnap, /fitstoo, /colorbar

        filename=frun+hdirs+thisnum+'_vmap.eps'
	center= usethiscenter
        contour_allstars_velocitymap, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, $
                /loadedsnap, /fitstoo, /colorbar

        filename=frun+hdirs+thisnum+'_smap.eps'
	center= usethiscenter
        contour_allstars_dispersionmap, frun, snapnum, xlen, 'ps', filename=filename, $
		rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		center= center, $
                /loadedsnap, /fitstoo, /colorbar

	endelse

	;snapnum= snapnum+5
	snapnum= snapnum+snapeveryX

endfor


print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end





