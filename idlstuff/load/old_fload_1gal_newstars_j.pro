Function fload_1gal_newstars_j, typeofam, startid, numpart, center=center, vcom=vcom, idtofollow=idtofollow

;  typeofam according to the following
; -------------------------------------
;  0 - total, i.e. one number, |J|
;  1 - x component, one number, |Jx|
;  2 - y component, one number, |Jy|
;  3 - z component, one number, |Jz|
; 10 - total for each particle, i.e. returns |J| for each particle
; 11 - |Jx| for each particle
; 12 - |Jy| for each particle
; 13 - |Jz| for each particle
; 44 - total specific, i.e. one number, |J|/Mtot, so 0 divided by total mass
; 87 - specific Jz of each particle
; 88 - specific of each particle, return |J|/m for each particle
; 99 - total |J| for individual particle, used for idtofollow


    COMMON GalaxyHeader
    COMMON OtherData
    COMMON Center
    COMMON NewStarData
    COMMON ParentID


    if npart(4) le 0 then return, [0]

    if keyword_set(center) then begin
        c= center
    endif else begin
        c= com
    endelse


    if not keyword_set(vcom) then begin
        vcom= [0.0, 0.0, 0.0]
    endif


    ; default
    x = xstars
    y = ystars
    z = zstars
    vx= vxstars
    vy= vystars
    vz= vzstars
    m= mstars



    ; now instead of grabbing new star id number's
    ;  we'll grab their parent id's
    if flag_stargens eq 1 and flag_parentid eq 1  then begin
        stars_ids= [parentid]
    endif else begin
        nsstartid= npart(0)+npart(1)+npart(2)+npart(3)
        stars_ids= id(nsstartid:N-1)
    endelse


    ; hey we'll just grab galaxy i's particles, then carry on
    if keyword_set(numpart) then begin

	idx = where((stars_ids ge long(startid)) and (stars_ids lt long(startid)+long(numpart)))

	x = xstars(idx)
	y = ystars(idx)
	z = zstars(idx)
	vx= vxstars(idx)
        vy= vystars(idx)
        vz= vzstars(idx)
        m= mstars(idx)
    endif



    ; used to follow individual particles, used with typeofam= 99
    if keyword_set(idtofollow) then begin
        stars_ids= id(0:npart(0)-1)
        idx= where(stars_ids EQ idtofollow)
        if idx(0) EQ -1 or n_elements(idx) NE 1 then begin
                print, "Hey, what's up with idx?  printing"
                print, idx
                return, [0,0,0]
        endif

	x= xstars(idx)
	y= ystars(idx)
	z= zstars(idx)
        vx= vxstars(idx)
        vy= vystars(idx)
        vz= vzstars(idx)
        m= mstars(idx)
    endif





    CASE typeofam OF

	0: BEGIN
	   ; total angular momentum - return jtot (i.e. one number)
	   ; --------------------------------------------------------
	   jx= m*((y-c(1))*(vz-vcom(2)) - (z-c(2))*(vy-vcom(1)))
	   jy= -m*((x-c(0))*(vz-vcom(2)) - (z-c(2))*(vx-vcom(0)))
	   jz= m*((x-c(0))*(vy-vcom(1)) - (y-c(1))*(vx-vcom(0)))

	   jx= total(jx)
	   jy= total(jy)
	   jz= total(jz)

	   jtot= sqrt(jx*jx + jy*jy + jz*jz)
	   return, total(jtot)
	END

	1: BEGIN
           jx= m*((y-c(1))*(vz-vcom(2)) - (z-c(2))*(vy-vcom(1)))
           return, total(jx)
	END

	2: BEGIN
           jy= - m*((x-c(0))*(vz-vcom(2)) - (z-c(2))*(vx-vcom(0)))
           return, total(jy)
	END

	3: BEGIN
           jz= m*((x-c(0))*(vy-vcom(1)) - (y-c(1))*(vx-vcom(0)))
           return, total(jz)
	END

	44: BEGIN
           ; specific total angular momentum - return j (i.e. one number)
           ; ---------------------------------------------------------------
           jx= m*((y-c(1))*(vz-vcom(2)) - (z-c(2))*(vy-vcom(1)))
           jy= -m*((x-c(0))*(vz-vcom(2)) - (z-c(2))*(vx-vcom(0)))
           jz= m*((x-c(0))*(vy-vcom(1)) - (y-c(1))*(vx-vcom(0)))

           jx= total(jx)
           jy= total(jy)
           jz= total(jz)
	   tm= total(m)

           jtot= sqrt(jx*jx + jy*jy + jz*jz)
           return, total(jtot)/tm
	END

	45: BEGIN
	   jz= m*((x-c(0))*(vy-vcom(1)) - (y-c(1))*(vx-vcom(0)))
	   tm= total(m)
	   return, total(jz)/tm
	END

        10: BEGIN
           ; total angular momentum -  returns jtot(n)
           ; --------------------------------------------
           jx= m*((y-c(1))*(vz-vcom(2)) - (z-c(2))*(vy-vcom(1)))
           jy= -m*((x-c(0))*(vz-vcom(2)) - (z-c(2))*(vx-vcom(0)))
           jz= m*((x-c(0))*(vy-vcom(1)) - (y-c(1))*(vx-vcom(0)))

           jtot= sqrt(jx*jx + jy*jy + jz*jz)
           return, jtot
        END

        11: BEGIN
           jx= m*((y-c(1))*(vz-vcom(2)) - (z-c(2))*(vy-vcom(1)))
           return, jx
        END

        12: BEGIN
           jy= - m*((x-c(0))*(vz-vcom(2)) - (z-c(2))*(vx-vcom(0)))
           return, jy
        END

        13: BEGIN
           jz= m*((x-c(0))*(vy-vcom(1)) - (y-c(1))*(vx-vcom(0)))
           return, jz
        END

        87: BEGIN
           ; specific z angular momentum - return j(n)
           ; -----------------------------------------------
           jz= (x-c(0))*(vy-vcom(1)) - (y-c(1))*(vx-vcom(0))
           return, jz
        END

        88: BEGIN
           ; specific total angular momentum - return j(n)
           ; -----------------------------------------------
           jx= (y-c(1))*(vz-vcom(2)) - (z-c(2))*(vy-vcom(1))
           jy= -(x-c(0))*(vz-vcom(2)) + (z-c(2))*(vx-vcom(0))
           jz= (x-c(0))*(vy-vcom(1)) - (y-c(1))*(vx-vcom(0))

           jtot= sqrt(jx*jx + jy*jy + jz*jz)
           return, jtot
        END

	99: BEGIN
	   jx= m*((y-c(1))*vz - (z-c(2))*vy)
           jy= -m*((x-c(0))*vz- (z-c(2))*vx)
           jz= m*((x-c(0))*vy- (y-c(1))*vx)

	   jtot= sqrt(jx*jx + jy*jy + jz*jz)
           return, total(jtot)
	END

    ELSE: print, "I don't know what to do now?"

    ENDCASE


end


