Function fload_newstars_j, typeofam, center=center, vcom=vcom, total=total

;  typeofam according to the following
;  note: this procedure return the total
;        angular momentum, while the 
;        procedure fload_dm_angmo returns
;        the ang mo per particle
; -------------------------------------
;  0 - total
;  1 - x component
;  2 - y component
;  3 - z component


    COMMON GalaxyHeader
    COMMON NewStarData
    COMMON Center

    if keyword_set(center) then begin
	c= center
    endif else begin
	c= com
    endelse


    if not keyword_set(vcom) then begin
	vcom= [0.0, 0.0, 0.0]
    endif


        ; are there stars?
        ; -----------------
        if npart(4) GT 0 then begin
          allmasses = mstars
	  allxs = xstars
	  allys = ystars
	  allzs = zstars
          allvxs = vxstars
          allvys = vystars
          allvzs = vzstars
        endif else begin
	  print, "NO NEWSTARS, sorry!"
	  return, [0]
	endelse




    CASE typeofam OF

	0: BEGIN
	   ; total angular momentum
	   ; ----------------------
	   jx= allmasses*((allys-c(1))*(allvzs-vcom(2))- (allzs-c(2))*(allvys-vcom(1)))
	   jy= -allmasses*((allxs-c(0))*(allvzs-vcom(2))- (allzs-c(2))*(allvxs-vcom(0)))
	   jz= allmasses*((allxs-c(0))*(allvys-vcom(1))- (allys-c(1))*(allvxs-vcom(0)))

	   n_j= n_elements(jx)

	   jx= total(jx)
	   jy= total(jy)
	   jz= total(jz)

	   jtot= sqrt(jx*jx + jy*jy + jz*jz)

	   print, n_j," particles have J_total= ", jtot, "    (",jx,",",jy,",",jz,")"


	   ; what is it without centering??
	   ; ---------------------------------
           ;jx= allmasses*(allys*(allvzs-vcom(2))- allzs*(allvys-vcom(1)))
           ;jy= -allmasses*(allxs*(allvzs-vcom(2))- allzs*(allvxs-vcom(0)))
           ;jz= allmasses*(allxs*(allvys-vcom(1))- allys*(allvxs-vcom(0)))
           ;jx= total(jx)
           ;jy= total(jy)
           ;jz= total(jz)
           ;print, "J_total= ",sqrt(jx*jx + jy*jy + jz*jz)

	   return, jtot
	END

	1: BEGIN
           jx= allmasses*((allys-c(1))*(allvzs-vcom(2))- (allzs-c(2))*(allvys-vcom(1)))
           return, total(jx)
	END

	2: BEGIN
           jy= - allmasses*((allxs-c(0))*(allvzs-vcom(2))- (allzs-c(2))*(allvxs-vcom(0)))
           return, total(jy)
	END

	3: BEGIN
           jz= allmasses*((allxs-c(0))*(allvys-vcom(1))- (allys-c(1))*(allvxs-vcom(0)))
           return, total(jz)
	END

        44: BEGIN
           ; specific total angular momentum
           ; -------------------------------
           jx= allmasses*((allys-c(1))*(allvzs-vcom(2))- (allzs-c(2))*(allvys-vcom(1)))
           jy= -allmasses*((allxs-c(0))*(allvzs-vcom(2))- (allzs-c(2))*(allvxs-vcom(0)))
           jz= allmasses*((allxs-c(0))*(allvys-vcom(1))- (allys-c(1))*(allvxs-vcom(0)))

	   n_j= n_elements(jx)

           jx= total(jx)
           jy= total(jy)
           jz= total(jz)
	   tm=total(allmasses)

           jtot= sqrt(jx*jx + jy*jy + jz*jz)

           print, n_j," particles have J_specf= ", jtot/tm, "    (",jx/tm,",",jy/tm,",",jz/tm,")"

           return, jtot/tm
        END

        10: BEGIN
           ; total angular momentum
           ; ----------------------
           jx= allmasses*((allys-c(1))*(allvzs-vcom(2)) - (allzs-c(2))*(allvys-vcom(1)))
           jy= -allmasses*((allxs-c(0))*(allvzs-vcom(2)) - (allzs-c(2))*(allvxs-vcom(0)))
           jz= allmasses*((allxs-c(0))*(allvys-vcom(1)) - (allys-c(1))*(allvxs-vcom(0)))

           jtot= sqrt(jx*jx + jy*jy + jz*jz)
           return, jtot
        END

        11: BEGIN
           jx= allmasses*((allys-c(1))*(allvzs-vcom(2)) - (allzs-c(2))*(allvys-vcom(1)))
           return, jx
        END

        12: BEGIN
           jy= - allmasses*((allxs-c(0))*(allvzs-vcom(2)) - (allzs-c(2))*(allvxs-vcom(0)))
           return, jy
        END

        13: BEGIN
           jz= allmasses*((allxs-c(0))*(allvys-vcom(1)) - (allys-c(1))*(allvxs-vcom(0)))
           return, jz
        END


	;   specific angular momentum
	; ------------------------------
        21: BEGIN
           jx= ((allys-c(1))*(allvzs-vcom(2)) - (allzs-c(2))*(allvys-vcom(1)))
           return, jx
        END

        22: BEGIN
           jy= - ((allxs-c(0))*(allvzs-vcom(2)) - (allzs-c(2))*(allvxs-vcom(0)))
           return, jy
        END

        23: BEGIN
           jz= ((allxs-c(0))*(allvys-vcom(1)) - (allys-c(1))*(allvxs-vcom(0)))
           return, jz
        END

    ELSE: print, "I don't know what to do now?"

    ENDCASE


end


