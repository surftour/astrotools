



;=================================================================================






pro compile_js, junk

;snapnum= long(30)
;determine_one_sim_js, "/raid4/tcox/vc3vc3b", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3c", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3d", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3e", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3f", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3g", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3h", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3i", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3j", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3k", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3l", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3m", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3n", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3o", snapnum
;determine_one_sim_js, "/raid4/tcox/vc3vc3p", snapnum

;snapnum= long(0)
;determine_one_sim_js, "/raid4/tcox/vc3vc3b", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3c", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3d", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3e", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3f", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3g", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3h", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3i", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3j", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3k", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3l", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3m", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3n", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3o", snapnum, /notmerged
;determine_one_sim_js, "/raid4/tcox/vc3vc3p", snapnum, /notmerged

;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3b", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3c", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3d", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3e", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3f", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3g", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3h", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3i", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3j", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3k", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3l", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3m", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3n", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3o", snapnum
;determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3p", snapnum

snapnum= long(0)
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3b", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3c", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3d", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3e", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3f", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3g", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3h", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3i", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3j", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3k", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3l", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3m", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3n", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3o", snapnum, /notmerged
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3p", snapnum, /notmerged

end








;============================================================
;
;
;
;    Calculate the angular momentum of
;	1. All Stars
;	2. Old Stars (total)
;	3.    "      (gal 1)
;	4.    "      (gal 2)
;	5. New Stars (total)
;	6.    "      (gal 1)
;	7.    "      (gal 2)
;	8. New Stars that are within 3 kpc
;		     (total)
;	9.    "      (gal 1)
;	10.   "      (gal 2)
;	11. New Stars that are within 3 kpc, and J is within 5 degrees of total
;		     (total)
;	12.   "      (gal 1)
;	13.   "      (gal 2)
;	14. Remnant Gas Disk (r<10 kpc)
;		     (total)
;	15.    "     (gal 1)
;	16.    "     (gal 2)
;
;
;
;============================================================
pro determine_one_sim_js, frun, snapnum, notmerged=notmerged

; ----------------------
;  Now Do the Work
; ----------------------

ok=fload_snapshot_bh(frun,snapnum)



; open angular momentum file
snaplbl= strcompress(string(snapnum),/remove_all)
openw, 1, frun+'/j_'+snaplbl+'.txt', ERROR=err
printf, 1, "#   "
printf, 1, "#  Angular Momentum  (Gadget Units) "
printf, 1, "#   "
printf, 1, "#              J_tot        J_x        J_y        J_z      theta       phi  "
printf, 1, "#               (GU)       (GU)       (GU)       (GU)     (deg.)     (deg.) "



; a)  All Stars
; ------------------
; ------------------
print, "---------"
print, "all stars"
print, " "

r_as= fload_allstars_xyz('r')

; total j's
jx= fload_allstars_j(1)   ; j_x
jy= fload_allstars_j(2)   ; j_y
jz= fload_allstars_j(3)   ; j_z
jtot= sqrt(jx*jx + jy*jy + jz*jz)
jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
jphi= atan(jy,jx) * 180.0 / !PI
print, "tot= ", jtot, "  (",jx,jy,jz," )"
print, "theta= ", jtheta
print, "phi= ", jphi

printf, 1, FORMAT= '("AllStars   ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi



; b)  New Stars
; ------------------
; ------------------
print, "---------"
print, "new stars"
print, " "

; total j's
if fload_npart(4) gt 0 then begin
	jx= fload_newstars_j(1)   ; j_x
	jy= fload_newstars_j(2)   ; j_y
	jz= fload_newstars_j(3)   ; j_z
	jtot= sqrt(jx*jx + jy*jy + jz*jz)
	jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
	jphi= atan(jy,jx) * 180.0 / !PI
	print, "tot= ", jtot, "  (",jx,jy,jz," )"
	print, "theta= ", jtheta
	print, "phi= ", jphi

	printf, 1, FORMAT= '("NewStars   ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

	; galaxy 1
	startid= 1L
	numpart= 200000L
	if keyword_set(notmerged) then begin
		center= fload_1gal_center(startid,numpart)
		comvel= fload_1gal_comvel(startid,numpart,center=center,rfact=0.01)
		jx= fload_1gal_newstars_j(1,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_x
        	jy= fload_1gal_newstars_j(2,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_y
        	jz= fload_1gal_newstars_j(3,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_z
	endif else begin
		jx= fload_1gal_newstars_j(1,startid=startid,numpart=numpart)   ; j_x
		jy= fload_1gal_newstars_j(2,startid=startid,numpart=numpart)   ; j_y
		jz= fload_1gal_newstars_j(3,startid=startid,numpart=numpart)   ; j_z
	endelse
	jtot= sqrt(jx*jx + jy*jy + jz*jz)
	jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
	jphi= atan(jy,jx) * 180.0 / !PI
	print, "tot= ", jtot, "  (",jx,jy,jz," )"
	;print, "theta= ", jtheta
	;print, "phi= ", jphi

	printf, 1, FORMAT= '("NewStars#1 ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

	; galaxy 2
	startid= 200001L
	numpart= 200000L
	if keyword_set(notmerged) then begin
		center= fload_1gal_center(startid,numpart)
		comvel= fload_1gal_comvel(startid,numpart,center=center,rfact=0.01)
		jx= fload_1gal_newstars_j(1,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_x
        	jy= fload_1gal_newstars_j(2,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_y
        	jz= fload_1gal_newstars_j(3,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_z
	endif else begin
		jx= fload_1gal_newstars_j(1,startid=startid,numpart=numpart)   ; j_x
		jy= fload_1gal_newstars_j(2,startid=startid,numpart=numpart)   ; j_y
		jz= fload_1gal_newstars_j(3,startid=startid,numpart=numpart)   ; j_z
	endelse
	jtot= sqrt(jx*jx + jy*jy + jz*jz)
	jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
	jphi= atan(jy,jx) * 180.0 / !PI
	print, "tot= ", jtot, "  (",jx,jy,jz," )"
	;print, "theta= ", jtheta
	;print, "phi= ", jphi

	printf, 1, FORMAT= '("NewStars#2 ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

endif else begin
	jtot= 0.0
	jx= 0.0
	jy= 0.0
	jz= 0.0
	jtheta= 0.0
	jphi= 0.0
	printf, 1, FORMAT= '("NewStars   ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi
	printf, 1, FORMAT= '("NewStars#1 ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi
	printf, 1, FORMAT= '("NewStars#2 ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi
endelse



; c)  Old Stars
; ------------------
; ------------------
print, "---------"
print, "old stars"
print, " "
r_old= fload_oldstars_xyz('r')

;specific j's
jx= fload_oldstars_j(1)
jy= fload_oldstars_j(2)
jz= fload_oldstars_j(3)
jtot= sqrt(jx*jx + jy*jy + jz*jz)
jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
jphi= atan(jy,jx) * 180.0 / !PI
print, "tot= ", jtot, "  (",jx,jy,jz," )"
print, "theta= ", jtheta
print, "phi= ", jphi

printf, 1, FORMAT= '("OldStars   ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

; galaxy 1
startid= 1L
numpart= 200000L
if keyword_set(notmerged) then begin
	center= fload_1gal_center(startid,numpart)
	comvel= fload_1gal_comvel(startid,numpart,center=center,rfact=0.01)
	jx= fload_1gal_oldstars_j(1,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_x
       	jy= fload_1gal_oldstars_j(2,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_y
       	jz= fload_1gal_oldstars_j(3,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_z
endif else begin
	jx= fload_1gal_oldstars_j(1,startid=startid,numpart=numpart)   ; j_x
	jy= fload_1gal_oldstars_j(2,startid=startid,numpart=numpart)   ; j_y
	jz= fload_1gal_oldstars_j(3,startid=startid,numpart=numpart)   ; j_z
endelse
jtot= sqrt(jx*jx + jy*jy + jz*jz)
jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
jphi= atan(jy,jx) * 180.0 / !PI
print, "tot= ", jtot, "  (",jx,jy,jz," )"
;print, "theta= ", jtheta
;print, "phi= ", jphi

printf, 1, FORMAT= '("OldStars#1 ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

; galaxy 2
startid= 200001L
numpart= 200000L
if keyword_set(notmerged) then begin
	center= fload_1gal_center(startid,numpart)
	comvel= fload_1gal_comvel(startid,numpart,center=center,rfact=0.01)
	jx= fload_1gal_oldstars_j(1,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_x
       	jy= fload_1gal_oldstars_j(2,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_y
       	jz= fload_1gal_oldstars_j(3,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_z
endif else begin
	jx= fload_1gal_oldstars_j(1,startid=startid,numpart=numpart)   ; j_x
	jy= fload_1gal_oldstars_j(2,startid=startid,numpart=numpart)   ; j_y
	jz= fload_1gal_oldstars_j(3,startid=startid,numpart=numpart)   ; j_z
endelse
jtot= sqrt(jx*jx + jy*jy + jz*jz)
jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
jphi= atan(jy,jx) * 180.0 / !PI
print, "tot= ", jtot, "  (",jx,jy,jz," )"
;print, "theta= ", jtheta
;print, "phi= ", jphi

printf, 1, FORMAT= '("OldStars#2 ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi



; d)  New Stars (r<5 kpc)
; -------------------------
; -------------------------
print, "-------------------"
print, "new stars (r<5 kpc)"
print, " "

xlen= 5.0

; total j's
if fload_npart(4) gt 0 then begin
	jx= fload_newstars_j(11)   ; j_x
	jy= fload_newstars_j(12)   ; j_y
	jz= fload_newstars_j(13)   ; j_z
	jr= fload_newstars_xyz('r')
	idx= where(jr lt xlen)
	jx= total(jx(idx))
	jy= total(jy(idx))
	jz= total(jz(idx))
	jtot= sqrt(jx*jx + jy*jy + jz*jz)
	jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
	jphi= atan(jy,jx) * 180.0 / !PI
	print, "tot= ", jtot, "  (",jx,jy,jz," )"
	print, "theta= ", jtheta
	print, "phi= ", jphi

	printf, 1, FORMAT= '("NewStars<5 ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

	; galaxy 1
	startid= 1L
	numpart= 200000L
	if keyword_set(notmerged) then begin
		center= fload_1gal_center(startid,numpart)
		comvel= fload_1gal_comvel(startid,numpart,center=center,rfact=0.01)
		jx= fload_1gal_newstars_j(11,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_x
        	jy= fload_1gal_newstars_j(12,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_y
        	jz= fload_1gal_newstars_j(13,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_z
		jr= fload_1gal_newstars_xyz('r',startid,numpart,center=center)
	endif else begin
		jx= fload_1gal_newstars_j(11,startid=startid,numpart=numpart)   ; j_x
		jy= fload_1gal_newstars_j(12,startid=startid,numpart=numpart)   ; j_y
		jz= fload_1gal_newstars_j(13,startid=startid,numpart=numpart)   ; j_z
		jr= fload_1gal_newstars_xyz('r',startid,numpart)
	endelse
	idx= where(jr lt xlen)
	jx= total(jx(idx))
	jy= total(jy(idx))
	jz= total(jz(idx))
	jtot= sqrt(jx*jx + jy*jy + jz*jz)
	jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
	jphi= atan(jy,jx) * 180.0 / !PI
	print, "tot= ", jtot, "  (",jx,jy,jz," )"
	;print, "theta= ", jtheta
	;print, "phi= ", jphi

	printf, 1, FORMAT= '("NewSt<5#1  ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

	; galaxy 2
	startid= 200001L
	numpart= 200000L
	if keyword_set(notmerged) then begin
		center= fload_1gal_center(startid,numpart)
		comvel= fload_1gal_comvel(startid,numpart,center=center,rfact=0.01)
		jx= fload_1gal_newstars_j(11,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_x
        	jy= fload_1gal_newstars_j(12,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_y
        	jz= fload_1gal_newstars_j(13,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_z
		jr= fload_1gal_newstars_xyz('r',startid,numpart,center=center)
	endif else begin
		jx= fload_1gal_newstars_j(11,startid=startid,numpart=numpart)   ; j_x
		jy= fload_1gal_newstars_j(12,startid=startid,numpart=numpart)   ; j_y
		jz= fload_1gal_newstars_j(13,startid=startid,numpart=numpart)   ; j_z
		jr= fload_1gal_newstars_xyz('r',startid,numpart)
	endelse
	idx= where(jr lt xlen)
	jx= total(jx(idx))
	jy= total(jy(idx))
	jz= total(jz(idx))
	jtot= sqrt(jx*jx + jy*jy + jz*jz)
	jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
	jphi= atan(jy,jx) * 180.0 / !PI
	print, "tot= ", jtot, "  (",jx,jy,jz," )"
	;print, "theta= ", jtheta
	;print, "phi= ", jphi

	printf, 1, FORMAT= '("NewSt<5#2  ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

endif else begin
        jtot= 0.0
        jx= 0.0
        jy= 0.0
        jz= 0.0
        jtheta= 0.0
        jphi= 0.0
	printf, 1, FORMAT= '("NewStars<5 ", 6(F9.2,"  "))', $
                jtot, jx, jy, jz, jtheta, jphi
	printf, 1, FORMAT= '("NewSt<5#1  ", 6(F9.2,"  "))', $
                jtot, jx, jy, jz, jtheta, jphi
	printf, 1, FORMAT= '("NewSt<5#2  ", 6(F9.2,"  "))', $
                jtot, jx, jy, jz, jtheta, jphi
endelse





; e)  New Stars (r< 5 kpc, aligned with total)
; ---------------------------------------------
; ---------------------------------------------
print, "---------"
print, "new stars (r< 5kpc, aligned with total)"
print, " "

xlen= 5.0


; total j's
; ----------
if fload_npart(4) gt 0 then begin
	ijx= fload_newstars_j(11)   ; j_x
	ijy= fload_newstars_j(12)   ; j_y
	ijz= fload_newstars_j(13)   ; j_z
	ijr= fload_newstars_xyz('r')
	ridx= where(ijr lt xlen)
	print, n_elements(ridx)," elements within 5 kpc"
	ijx= ijx(ridx)
	ijy= ijy(ridx)
	ijz= ijz(ridx)
	ijtot= sqrt(ijx*ijx + ijy*ijy + ijz*ijz)
	jx= total(ijx)
	jy= total(ijy)
	jz= total(ijz)
	jtot= sqrt(jx*jx + jy*jy + jz*jz)
	n_x= jx/jtot
	n_y= jy/jtot
	n_z= jz/jtot
	jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
	jphi= atan(jy,jx) * 180.0 / !PI
	print, "tot= ", jtot, "  (",jx,jy,jz," )"
	print, "theta= ", jtheta
	print, "phi= ", jphi
	;help, ridx

	ij_dot_n= ijx*n_x + ijy*n_y + ijz*n_z
	thetadiff= acos(ij_dot_n/ijtot)*180./!PI

	; grab aligned angular momenta??
	;idxx= where(thetadiff lt 10.0)
	idxx= where(thetadiff lt 5.0)
	;idxx= where(thetadiff lt 3.0)
	print, n_elements(idxx)," elements within 5 kpc, and aligned"

	ajx= total(ijx(idxx))
	ajy= total(ijy(idxx))
	ajz= total(ijz(idxx))
	ajtot= sqrt(ajx*ajx + ajy*ajy + ajz*ajz)
	ajtheta= atan(sqrt(ajx*ajx + ajy*ajy), ajz) * 180.0 / !PI
	ajphi= atan(ajy,ajx) * 180.0 / !PI

	printf, 1, FORMAT= '("NewStAlgned", 6(F9.2,"  "))', $
		ajtot, ajx, ajy, ajz, ajtheta, ajphi

	; galaxy 1
	startid= 1L
	numpart= 200000L
	if keyword_set(notmerged) then begin
		center= fload_1gal_center(startid,numpart)
		comvel= fload_1gal_comvel(startid,numpart,center=center,rfact=0.01)
		ijx= fload_1gal_newstars_j(11,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_x
        	ijy= fload_1gal_newstars_j(12,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_y
        	ijz= fload_1gal_newstars_j(13,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_z
		ijr= fload_1gal_newstars_xyz('r',startid,numpart,center=center)
	endif else begin
		ijx= fload_1gal_newstars_j(11,startid=startid,numpart=numpart)   ; j_x
		ijy= fload_1gal_newstars_j(12,startid=startid,numpart=numpart)   ; j_y
		ijz= fload_1gal_newstars_j(13,startid=startid,numpart=numpart)   ; j_z
		ijr= fload_1gal_newstars_xyz('r',startid,numpart)
	endelse
	ridx= where(ijr lt xlen)
	ijx= ijx(ridx)
	ijy= ijy(ridx)
	ijz= ijz(ridx)
	ijtot= sqrt(ijx*ijx + ijy*ijy + ijz*ijz)
	;help, ridx

	ij_dot_n= ijx*n_x + ijy*n_y + ijz*n_z
	thetadiff= acos(ij_dot_n/ijtot)*180./!PI

	; grab aligned angular momenta??
	;idxx= where(thetadiff lt 10.0)
	idxx= where(thetadiff lt 5.0)
	;idxx= where(thetadiff lt 3.0)
	print, n_elements(idxx)," elements within 5 kpc, and aligned"

	ajx= total(ijx(idxx))
	ajy= total(ijy(idxx))
	ajz= total(ijz(idxx))
	ajtot= sqrt(ajx*ajx + ajy*ajy + ajz*ajz)
	ajtheta= atan(sqrt(ajx*ajx + ajy*ajy), ajz) * 180.0 / !PI
	ajphi= atan(ajy,ajx) * 180.0 / !PI

	printf, 1, FORMAT= '("NewStAlgn#1", 6(F9.2,"  "))', $
		ajtot, ajx, ajy, ajz, ajtheta, ajphi

	; galaxy 2
	startid= 200001L
	numpart= 200000L
	if keyword_set(notmerged) then begin
		center= fload_1gal_center(startid,numpart)
		comvel= fload_1gal_comvel(startid,numpart,center=center,rfact=0.01)
		ijx= fload_1gal_newstars_j(11,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_x
        	ijy= fload_1gal_newstars_j(12,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_y
        	ijz= fload_1gal_newstars_j(13,startid=startid,numpart=numpart,center=center,vcom=comvel)   ; j_z
		ijr= fload_1gal_newstars_xyz('r',startid,numpart,center=center)
	endif else begin
		ijx= fload_1gal_newstars_j(11,startid=startid,numpart=numpart)   ; j_x
		ijy= fload_1gal_newstars_j(12,startid=startid,numpart=numpart)   ; j_y
		ijz= fload_1gal_newstars_j(13,startid=startid,numpart=numpart)   ; j_z
		ijr= fload_1gal_newstars_xyz('r',startid,numpart)
	endelse
	ridx= where(ijr lt xlen)
	ijx= ijx(ridx)
	ijy= ijy(ridx)
	ijz= ijz(ridx)
	ijtot= sqrt(ijx*ijx + ijy*ijy + ijz*ijz)
	;help, ridx

	ij_dot_n= ijx*n_x + ijy*n_y + ijz*n_z
	thetadiff= acos(ij_dot_n/ijtot)*180./!PI

	; grab aligned angular momenta??
	;idxx= where(thetadiff lt 10.0)
	idxx= where(thetadiff lt 5.0)
	;idxx= where(thetadiff lt 3.0)
	print, n_elements(idxx)," elements within 5 kpc, and aligned"

	ajx= total(ijx(idxx))
	ajy= total(ijy(idxx))
	ajz= total(ijz(idxx))
	ajtot= sqrt(ajx*ajx + ajy*ajy + ajz*ajz)
	ajtheta= atan(sqrt(ajx*ajx + ajy*ajy), ajz) * 180.0 / !PI
	ajphi= atan(ajy,ajx) * 180.0 / !PI

	printf, 1, FORMAT= '("NewStAlgn#2", 6(F9.2,"  "))', $
		ajtot, ajx, ajy, ajz, ajtheta, ajphi

endif else begin
        jtot= 0.0
        jx= 0.0
        jy= 0.0
        jz= 0.0
        jtheta= 0.0
        jphi= 0.0
	printf, 1, FORMAT= '("NewStAlgned", 6(F9.2,"  "))', $
                jtot, jx, jy, jz, jtheta, jphi
	printf, 1, FORMAT= '("NewStAlgn#1", 6(F9.2,"  "))', $
                jtot, jx, jy, jz, jtheta, jphi
	printf, 1, FORMAT= '("NewStAlgn#2", 6(F9.2,"  "))', $
                jtot, jx, jy, jz, jtheta, jphi
endelse





; f) Gas Disk (r<10 kpc) 
; ------------------
; ------------------
print, "---------"
print, "gas disk (r<10 kpc)"
print, " "

xlen= 10.0

;specific j's
if fload_npart(0) gt 0 then begin
	jx= fload_gas_j(11)
	jy= fload_gas_j(12)
	jz= fload_gas_j(13)
	jr= fload_gas_xyz('r')
	idx= where(jr lt xlen)
	if idx(0) eq -1 then begin
	   jtot= 0 & jx= 0 & jy= 0 & jz= 0 & jtheta= 0 & jphi= 0
	endif else begin
	   jx= total(jx(idx))
	   jy= total(jy(idx))
	   jz= total(jz(idx))
	   jtot= sqrt(jx*jx + jy*jy + jz*jz)
	   jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
	   jphi= atan(jy,jx) * 180.0 / !PI
	endelse
	print, "tot= ", jtot, "  (",jx,jy,jz," )"
	print, "theta= ", jtheta
	print, "phi= ", jphi

	printf, 1, FORMAT= '("GasDisk    ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

	; galaxy 1
	startid= 1L
	numpart= 200000L
	if keyword_set(notmerged) then begin
		center= fload_1gal_center(startid,numpart)
		comvel= fload_1gal_comvel(startid,numpart,center=center,rfact=0.01)
		jx= fload_1gal_gas_j(11,startid,numpart,center=center,vcom=comvel)   ; j_x
        	jy= fload_1gal_gas_j(12,startid,numpart,center=center,vcom=comvel)   ; j_y
        	jz= fload_1gal_gas_j(13,startid,numpart,center=center,vcom=comvel)   ; j_z
		jr= fload_1gal_gas_xyz('r',startid,numpart,center=center)
	endif else begin
		jx= fload_1gal_gas_j(11,startid,numpart)   ; j_x
		jy= fload_1gal_gas_j(12,startid,numpart)   ; j_y
		jz= fload_1gal_gas_j(13,startid,numpart)   ; j_z
		jr= fload_1gal_gas_xyz('r',startid,numpart)
	endelse
	idx= where(jr lt xlen)
	jx= total(jx(idx))
	jy= total(jy(idx))
	jz= total(jz(idx))
	jtot= sqrt(jx*jx + jy*jy + jz*jz)
	jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
	jphi= atan(jy,jx) * 180.0 / !PI
	print, "tot= ", jtot, "  (",jx,jy,jz," )"
	;print, "theta= ", jtheta
	;print, "phi= ", jphi

	printf, 1, FORMAT= '("GasDisk#1  ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

	; galaxy 2
	startid= 200001L
	numpart= 200000L
	if keyword_set(notmerged) then begin
		center= fload_1gal_center(startid,numpart)
		comvel= fload_1gal_comvel(startid,numpart,center=center,rfact=0.01)
		jx= fload_1gal_gas_j(11,startid,numpart,center=center,vcom=comvel)   ; j_x
        	jy= fload_1gal_gas_j(12,startid,numpart,center=center,vcom=comvel)   ; j_y
        	jz= fload_1gal_gas_j(13,startid,numpart,center=center,vcom=comvel)   ; j_z
		jr= fload_1gal_gas_xyz('r',startid,numpart,center=center)
	endif else begin
		jx= fload_1gal_gas_j(11,startid,numpart)   ; j_x
		jy= fload_1gal_gas_j(12,startid,numpart)   ; j_y
		jz= fload_1gal_gas_j(13,startid,numpart)   ; j_z
		jr= fload_1gal_gas_xyz('r',startid,numpart)
	endelse
	idx= where(jr lt xlen)
	jx= total(jx(idx))
	jy= total(jy(idx))
	jz= total(jz(idx))
	jtot= sqrt(jx*jx + jy*jy + jz*jz)
	jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
	jphi= atan(jy,jx) * 180.0 / !PI
	print, "tot= ", jtot, "  (",jx,jy,jz," )"
	;print, "theta= ", jtheta
	;print, "phi= ", jphi

	printf, 1, FORMAT= '("GasDisk#2  ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

endif else begin
        jtot= 0.0
        jx= 0.0
        jy= 0.0
        jz= 0.0
        jtheta= 0.0
        jphi= 0.0
	printf, 1, FORMAT= '("GasDisk    ", 6(F9.2,"  "))', $
                jtot, jx, jy, jz, jtheta, jphi
	printf, 1, FORMAT= '("GasDisk#1  ", 6(F9.2,"  "))', $
                jtot, jx, jy, jz, jtheta, jphi
	printf, 1, FORMAT= '("GasDisk#2  ", 6(F9.2,"  "))', $
                jtot, jx, jy, jz, jtheta, jphi
endelse





; done
;---------
close, 1




end







;===================================================================================
;
;
;===================================================================================





pro read_j_file, filename, jdata


spawn, "wc "+filename,result
lines=long(result)
datalines=lines(0)-5
;print, 'datalines= ', datalines, '   (there should be 16)'
jdata= fltarr(datalines,6)

openr, 1, filename
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
        readf, 1, junk
        tempjunk= strsplit(junk,/extract,count=count)
        whichj= tempjunk(0)
        jdata[i,0]= float(tempjunk(1))                 ; total j
        jdata[i,1]= float(tempjunk(2))                 ; jx
        jdata[i,2]= float(tempjunk(3))                 ; jy
        jdata[i,3]= float(tempjunk(4))                 ; jz
        jdata[i,4]= float(tempjunk(5))                 ; theta
        jdata[i,5]= float(tempjunk(6))                 ; phi

endfor

close, 1


end







;-----------------------------------------------------------------------------------






pro fload_alignments_foronerun, jfrun, $
		nslt5_x_os, $
		nslt5_x_gd, $
		os_x_zaxis, $
		t0spin1_x_osg1, $
		t0spin2_x_osg2, $
		t0spin1_x_nsg1, $
		t0spin2_x_nsg2, $
		t0spin1_x_gd, $
		t0spin2_x_gd, $
		t0spin1_x_nslt5g1, $
		t0spin2_x_nslt5g2, $
		t0spin1_x_nsal, $
		t0spin2_x_nsal



        jfile= jfrun+'/j_30.txt'
        read_j_file, jfile, jdata

	; all vectors are normalized
	as=           transpose(jdata[0,1:3] / jdata[0,0])
	ns=           transpose(jdata[1,1:3] / jdata[1,0])
	nsg1=         transpose(jdata[2,1:3] / jdata[2,0])
	nsg2=         transpose(jdata[3,1:3] / jdata[3,0])
	os=           transpose(jdata[4,1:3] / jdata[4,0])
	osg1=         transpose(jdata[5,1:3] / jdata[5,0])
	osg2=         transpose(jdata[6,1:3] / jdata[6,0])
	nslt5=        transpose(jdata[7,1:3] / jdata[7,0])
	nslt5g1=      transpose(jdata[8,1:3] / jdata[8,0])
	nslt5g2=      transpose(jdata[9,1:3] / jdata[9,0])
	nsal=         transpose(jdata[10,1:3] / jdata[10,0])
	nsalg1=       transpose(jdata[11,1:3] / jdata[11,0])
	nsalg2=       transpose(jdata[12,1:3] / jdata[12,0])
	gd=           transpose(jdata[13,1:3] / jdata[13,0])
	gdg1=         transpose(jdata[14,1:3] / jdata[14,0])
	gdg2=         transpose(jdata[15,1:3] / jdata[15,0])


	nslt5_x_os = 180. / !PI * acos(nslt5[0]*os[0] + nslt5[1]*os[1] + nslt5[2]*os[2])

	nslt5_x_gd = 180. / !PI * acos(nslt5[0]*gd[0] + nslt5[1]*gd[1] + nslt5[2]*gd[2])

	os_x_zaxis = 180. / !PI * acos(os[2])
	print, os_x_zaxis, '  vs. os theta= ', jdata[4,4]


	; original (t=0) j
	;---------------------- 
        jfile= jfrun+'/j_0.txt'
        read_j_file, jfile, jdata
   
        ; all vectors are normalized
        t0spin1=         transpose(jdata[5,1:3] / jdata[5,0])
        t0spin2=         transpose(jdata[6,1:3] / jdata[6,0])



	t0spin1_x_osg1 = 180. / !PI * acos(t0spin1[0]*osg1[0] + t0spin1[1]*osg1[1] + t0spin1[2]*osg1[2])
	t0spin2_x_osg2 = 180. / !PI * acos(t0spin2[0]*osg2[0] + t0spin2[1]*osg2[1] + t0spin2[2]*osg2[2])

	t0spin1_x_nsg1 = 180. / !PI * acos(t0spin1[0]*nsg1[0] + t0spin1[1]*nsg1[1] + t0spin1[2]*nsg1[2])
	t0spin2_x_nsg2 = 180. / !PI * acos(t0spin2[0]*nsg2[0] + t0spin2[1]*nsg2[1] + t0spin2[2]*nsg2[2])

	t0spin1_x_nslt5g1 = 180. / !PI * acos(t0spin1[0]*nslt5g1[0] + t0spin1[1]*nslt5g1[1] + t0spin1[2]*nslt5g1[2])
	t0spin2_x_nslt5g2 = 180. / !PI * acos(t0spin2[0]*nslt5g2[0] + t0spin2[1]*nslt5g2[1] + t0spin2[2]*nslt5g2[2])
	t0spin1_x_nsal = 180. / !PI * acos(t0spin1[0]*nsal[0] + t0spin1[1]*nsal[1] + t0spin1[2]*nsal[2])
	t0spin2_x_nsal = 180. / !PI * acos(t0spin2[0]*nsal[0] + t0spin2[1]*nsal[1] + t0spin2[2]*nsal[2])


	t0spin1_x_gd = 180. / !PI * acos(t0spin1[0]*gd[0] + t0spin1[1]*gd[1] + t0spin1[2]*gd[2])
	t0spin2_x_gd = 180. / !PI * acos(t0spin2[0]*gd[0] + t0spin2[1]*gd[1] + t0spin2[2]*gd[2])


end








;-----------------------------------------------------------------------------------






pro fload_all_alignments, jfruns, $
                a_nslt5_x_os, $
		a_nslt5_x_gd, $
                a_os_x_zaxis, $
                a_t0spin1_x_osg1, $
                a_t0spin2_x_osg2, $
                a_t0spin1_x_nsg1, $
                a_t0spin2_x_nsg2, $
                a_t0spin1_x_gd, $
		a_t0spin2_x_gd, $
		a_t0spin1_x_nslt5g1, $
		a_t0spin2_x_nslt5g2, $
		a_t0spin1_x_nsal, $
		a_t0spin2_x_nsal


	ns= n_elements(jfruns)

	a_nslt5_x_os= fltarr(ns)
	a_nslt5_x_gd= fltarr(ns)
	a_os_x_zaxis= fltarr(ns)
	a_t0spin1_x_osg1= fltarr(ns)
	a_t0spin2_x_osg2= fltarr(ns)
	a_t0spin1_x_nsg1= fltarr(ns)
	a_t0spin2_x_nsg2= fltarr(ns)
	a_t0spin1_x_gd= fltarr(ns)
	a_t0spin2_x_gd= fltarr(ns)
	a_t0spin1_x_nslt5g1= fltarr(ns)
	a_t0spin2_x_nslt5g2= fltarr(ns)
	a_t0spin1_x_nsal= fltarr(ns)
	a_t0spin2_x_nsal= fltarr(ns)

   	for i=0,ns-1 do begin

	    fload_alignments_foronerun, jfruns[i], $
		nslt5_x_os, $
		nslt5_x_gd, $
		os_x_zaxis, $
		t0spin1_x_osg1, $
		t0spin2_x_osg2, $
		t0spin1_x_nsg1, $
		t0spin2_x_nsg2, $
		t0spin1_x_gd, $
		t0spin2_x_gd, $
		t0spin1_x_nslt5g1, $
		t0spin2_x_nslt5g2, $
		t0spin1_x_nsal, $
		t0spin2_x_nsal


		a_nslt5_x_os[i]=          nslt5_x_os
		a_nslt5_x_gd[i]=          nslt5_x_gd
                a_os_x_zaxis[i]=          os_x_zaxis
                a_t0spin1_x_osg1[i]=      t0spin1_x_osg1
                a_t0spin2_x_osg2[i]=      t0spin2_x_osg2
                a_t0spin1_x_nsg1[i]=      t0spin1_x_nsg1
                a_t0spin2_x_nsg2[i]=      t0spin2_x_nsg2
                a_t0spin1_x_gd[i]=        t0spin1_x_gd
                a_t0spin2_x_gd[i]=        t0spin2_x_gd
		a_t0spin1_x_nslt5g1[i]=   t0spin1_x_nslt5g1
		a_t0spin2_x_nslt5g2[i]=   t0spin2_x_nslt5g2
		a_t0spin1_x_nsal[i]=      t0spin1_x_nsal
		a_t0spin2_x_nsal[i]=      t0spin2_x_nsal

	endfor

end









;-----------------------------------------------------------------------------------








pro alignment_hist, junk



if not keyword_set(junk) then begin
        print, " "
        print, " alignment_hist, junk"
        print, " "
        print, " "
        return 
endif




; -------------------
filename='ahist.eps'



initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


; -----------------------------
; set up constants
; -----------------------------


; alignement angle
xaxistitle= "!6Alignment (degrees)"
xmax = 180.0
xmin = 0.0

; number (histogram)
yaxistitle= ' '
ymax = 11.0
;ymax = 1.05
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.05
y0= 0.15
x1= 0.95
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
        ;/ylog, $
        ;/xlog, $ 
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
        xtitle=xaxistitle, $
        ytickformat='(a1)', $
        ;ytitle=yaxistitle, $
        /nodata   ;, /noerase






; =====================
; =====================


;  std's
; -------------
do_h1= 1
;do_h1= 0
if do_h1 eq 1 then begin
	fruns= ['vc3vc3h','vc3vc3b','vc3vc3c','vc3vc3d','vc3vc3e', $
        	'vc3vc3f','vc3vc3g','vc3vc3i','vc3vc3j','vc3vc3k', $
        	'vc3vc3l','vc3vc3m','vc3vc3n','vc3vc3o','vc3vc3p'] 


	jfruns= '/raid4/tcox/'+fruns(*)


	fload_all_alignments, jfruns, $
                a_nslt5_x_os, $
                a_nslt5_x_gd, $
                a_os_x_zaxis, $
                a_t0spin1_x_osg1, $
                a_t0spin2_x_osg2, $
                a_t0spin1_x_nsg1, $
                a_t0spin2_x_nsg2, $
                a_t0spin1_x_gd, $
                a_t0spin2_x_gd, $
		a_t0spin1_x_nslt5g1, $
		a_t0spin2_x_nslt5g2, $
		a_t0spin1_x_nsal, $
		a_t0spin2_x_nsal



	; print whatever histograms we want
	;-----------------------------------

	; are old stars aligned with z?
	;oplotit= 150
	;temp= process_histogram(a_os_x_zaxis, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
	;print, min(temp), max(temp)
	;xyouts, 0.60, 0.86, "os_x_zaxis", size=1.2, color=oplotit, /normal, charthick=3.0

	; ------------------------------
        ; are old stars aligned with z?
        ;oplotit= 150
        ;temp= process_histogram(a_t0spin1_x_osg1, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
        ;print, min(temp), max(temp)
        ;xyouts, 0.60, 0.90, "t0spin1_x_osg1", size=1.2, color=oplotit, /normal, charthick=3.0
        ;oplotit= 50
        ;temp= process_histogram(a_t0spin2_x_osg2, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
        ;print, min(temp), max(temp)
        ;xyouts, 0.60, 0.86, "t0spin2_x_osg2", size=1.2, color=oplotit, /normal, charthick=3.0

	; are old stars aligned with galaxy initial spins ?
        ;oplotit= 100
        ;temp= process_histogram(a_t0spin1_x_nsg1, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
        ;print, min(temp), max(temp)
        ;xyouts, 0.60, 0.82, "t0spin1_x_nsg1", size=1.2, color=oplotit, /normal, charthick=3.0
        ;oplotit= 0
        ;temp= process_histogram(a_t0spin2_x_nsg2, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
        ;print, min(temp), max(temp)
        ;xyouts, 0.60, 0.78, "t0spin2_x_nsg2", size=1.2, color=oplotit, /normal, charthick=3.0


        ; ------------------------------
        ; are the inner new stars aligned with initial galaxy spins ?
        ;oplotit= 150
        ;temp= process_histogram(a_t0spin1_x_nslt5g1, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
        ;print, min(temp), max(temp)
        ;xyouts, 0.60, 0.90, "t0spin1_x_nslt5g1", size=1.2, color=oplotit, /normal, charthick=3.0
        ;oplotit= 50
        ;temp= process_histogram(a_t0spin2_x_nslt5g2, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
        ;print, min(temp), max(temp)
        ;xyouts, 0.60, 0.86, "t0spin2_x_nslt5g2", size=1.2, color=oplotit, /normal, charthick=3.0

        ; are inner aligned new stars aligned with galaxy spins?
        ;oplotit= 100
        ;temp= process_histogram(a_t0spin1_x_nsal, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
        ;print, min(temp), max(temp)
        ;xyouts, 0.60, 0.82, "t0spin1_x_nsal", size=1.2, color=oplotit, /normal, charthick=3.0
        ;oplotit= 0
        ;temp= process_histogram(a_t0spin2_x_nsal, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
        ;print, min(temp), max(temp)
        ;xyouts, 0.60, 0.78, "t0spin2_x_nsal", size=1.2, color=oplotit, /normal, charthick=3.0

	; ------------------------------
        ; is the gas disk aligned with anything?
        oplotit= 150
        temp= process_histogram(a_t0spin1_x_gd, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
        print, min(temp), max(temp)
        xyouts, 0.60, 0.90, "t0spin1_x_gd", size=1.2, color=oplotit, /normal, charthick=3.0
        oplotit= 50
        temp= process_histogram(a_t0spin2_x_gd, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
        print, min(temp), max(temp)
        xyouts, 0.60, 0.86, "t0spin2_x_gd", size=1.2, color=oplotit, /normal, charthick=3.0
        oplotit= 100
        temp= process_histogram(a_nslt5_x_gd, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
        print, min(temp), max(temp)
        xyouts, 0.60, 0.82, "nslt5_x_gd", size=1.2, color=oplotit, /normal, charthick=3.0




endif



;  collisionless
; -----------------

;do_h2= 1
do_h2= 0
if do_h2 eq 1 then begin
	fruns= ['cvc3vc3h','cvc3vc3b','cvc3vc3c','cvc3vc3d','cvc3vc3e', $
        	'cvc3vc3f','cvc3vc3g','cvc3vc3i','cvc3vc3j','cvc3vc3k', $
        	'cvc3vc3l','cvc3vc3m','cvc3vc3n','cvc3vc3o','cvc3vc3p']

        jfruns= '/raid4/tcox/collisionless/'+fruns(*)


        fload_all_alignments, jfruns, $
                a_nslt5_x_os, $
                a_nslt5_x_gd, $
                a_os_x_zaxis, $
                a_t0spin1_x_osg1, $
                a_t0spin2_x_osg2, $
                a_t0spin1_x_nsg1, $
                a_t0spin2_x_nsg2, $
                a_t0spin1_x_gd, $
                a_t0spin2_x_gd, $
		a_t0spin1_x_nslt5g1, $
		a_t0spin2_x_nslt5g2, $
		a_t0spin1_x_nsal, $
		a_t0spin2_x_nsal



        ; print whatever histograms we want
        ;-----------------------------------

        ; are old stars aligned with z?
        oplotit= 50
        temp= process_histogram(a_os_x_zaxis, xmax=xmax, xmin=xmin, levels=9, oplotit=oplotit, mannorm=1)
        print, min(temp), max(temp)
        xyouts, 0.60, 0.90, "os_x_zaxis (c)", size=1.2, color=oplotit, /normal, charthick=3.0




endif




; -----------------------------------------------------------------------------



; -------------
;  Done
; -------------

device, /close



end







