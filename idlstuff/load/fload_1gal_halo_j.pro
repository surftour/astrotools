Function fload_1gal_halo_j, typeofam, halo_startid, halo_numpart, center=center, vcom=vcom, idtofollow=idtofollow

;  typeofam according to the following
; -------------------------------------
;  0 - total
;  1 - x component
;  2 - y component
;  3 - z component


    COMMON GalaxyHeader
    COMMON HaloData
    COMMON OtherData
    COMMON Center


    if keyword_set(center) then begin
        c= center
    endif else begin
        c= com
    endelse


    if not keyword_set(vcom) then begin
        vcom= [0.0, 0.0, 0.0]
    endif

    ; default
    x = xhalo
    y = yhalo
    z = zhalo
    vx= vxhalo
    vy= vyhalo
    vz= vzhalo
    m= mhalo


    ; hey we'll just grad galaxy i's particles, then carry on
    if keyword_set(halo_numpart) then begin
        hid = id(npart(0):npart(0)+npart(1)-1)
        idx = where((hid ge long(halo_startid)) and (hid lt long(halo_startid)+long(halo_numpart)))

        x = xhalo(idx)
        y = yhalo(idx)
        z = zhalo(idx)
        vx= vxhalo(idx)
        vy= vyhalo(idx)
        vz= vzhalo(idx)
        m= mhalo(idx)
    endif


    ; used to follow individual particles, used with typeofam= 99
    if keyword_set(idtofollow) then begin
	sids= npart(0)
	eids= npart(0)+npart(1)-1
        halo_ids= id(sids:eids)
        idx= where(halo_ids EQ idtofollow)
        if idx(0) EQ -1 or n_elements(idx) NE 1 then begin
                print, "Hey, what's up with idx?  printing"
                print, idx
                return, [0,0,0]
        endif

	x= xhalo(idx)
	y= yhalo(idx)
	z= zhalo(idx)
        vx= vxhalo(idx)
        vy= vyhalo(idx)
        vz= vzhalo(idx)
        m= mhalo(idx)
    endif



    CASE typeofam OF

	0: BEGIN
	   ; total angular momentum
	   ; ----------------------
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
           ; specific total angular momentum
           ; -------------------------------
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
           ; (individual) total angular momentum
           ; -------------------------------------
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
           jz= ((x-c(0))*(vy-vcom(1)) - (y-c(1))*(vx-vcom(0)))
	   ;tm= total(m)
           return, jz
        END

        88: BEGIN
           ; specific total angular momentum
           ; -------------------------------
           jx= ((y-c(1))*(vz-vcom(2)) - (z-c(2))*(vy-vcom(1)))
           jy= -((x-c(0))*(vz-vcom(2)) - (z-c(2))*(vx-vcom(0)))
           jz= ((x-c(0))*(vy-vcom(1)) - (y-c(1))*(vx-vcom(0)))

	   tm= 1.0
           ;tm= total(m)

           jtot= sqrt(jx*jx + jy*jy + jz*jz)
           return, jtot/tm
        END

	; NOTE:  this one is for follow id but will probably be eliminated
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


