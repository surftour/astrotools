Function fload_allstars_j, typeofam, center=center, vcom=vcom, total=total

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
    COMMON DiskData
    COMMON BulgeData
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


    if npart(2)+npart(3)+npart(4) lt 1 then begin
	print, "No stars: halt!!!"
	return, [0]
    endif

        ; is there a disk?
        ; -----------------
        if npart(2) GT 0 then begin
          allmasses = [mdisk]
	  allxs = [xdisk]
	  allys = [ydisk]
	  allzs = [zdisk]
          allvxs = [vxdisk]
          allvys = [vydisk]
          allvzs = [vzdisk]
        endif

        ; is there a bulge?
        ; -----------------
        if npart(3) GT 0 then begin
	  if npart(2) gt 0 then begin
            allmasses = [allmasses, mbulge]
	    allxs = [allxs, xbulge]
	    allys = [allys, ybulge]
	    allzs = [allzs, zbulge]
            allvxs = [allvxs, vxbulge]
            allvys = [allvys, vybulge]
            allvzs = [allvzs, vzbulge]
	  endif else begin
            allmasses = [mbulge]
	    allxs = [xbulge]
	    allys = [ybulge]
	    allzs = [zbulge]
            allvxs = [vxbulge]
            allvys = [vybulge]
            allvzs = [vzbulge]
	  endelse
        endif

	; is there new stars?
	; --------------------
        if npart(4) GT 0 then begin
          if ((npart(3) GT 0) or (npart(2) GT 0)) then begin
            allmasses = [allmasses, mstars]
            allxs = [allxs, xstars]
            allys = [allys, ystars]
            allzs = [allzs, zstars]
            allvxs = [allvxs, vxstars]
            allvys = [allvys, vystars]
            allvzs = [allvzs, vzstars]
          endif else begin
            allmasses = [mstars]
            allxs = [xstars]
            allys = [ystars]
            allzs = [zstars]
            allvxs = [vxstars]
            allvys = [vystars]
            allvzs = [vzstars]
          endelse
	endif





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

        33: BEGIN
           ; total angular momentum
           ; ----------------------
           jx= allmasses*((allys-c(1))*(allvzs-vcom(2))- (allzs-c(2))*(allvys-vcom(1)))
           jy= -allmasses*((allxs-c(0))*(allvzs-vcom(2))- (allzs-c(2))*(allvxs-vcom(0)))
           jz= allmasses*((allxs-c(0))*(allvys-vcom(1))- (allys-c(1))*(allvxs-vcom(0)))

           ; n_j= n_elements(jx)

           jtot= sqrt(jx*jx + jy*jy + jz*jz)

           ; print, n_j," particles have J_total= ", jtot, "    (",jx,",",jy,",",jz,")"

           return, jtot
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
        20: BEGIN
           jx= ((allys-c(1))*(allvzs-vcom(2)) - (allzs-c(2))*(allvys-vcom(1)))
           jy= -((allxs-c(0))*(allvzs-vcom(2)) - (allzs-c(2))*(allvxs-vcom(0)))
           jz= ((allxs-c(0))*(allvys-vcom(1)) - (allys-c(1))*(allvxs-vcom(0)))

           jtot= sqrt(jx*jx + jy*jy + jz*jz)
           return, jtot
        END

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


