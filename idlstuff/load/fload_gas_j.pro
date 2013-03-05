Function fload_gas_j, typeofam, center=center, vcom=vcom, idtofollow=idtofollow

;  typeofam according to the following
; -------------------------------------
;  0 - total
;  1 - x component
;  2 - y component
;  3 - z component


    COMMON GalaxyHeader
    COMMON GasData
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


    ; used to follow individual particles, used with typeofam= 99
    if keyword_set(idtofollow) then begin
        gas_ids= id(0:npart(0)-1)
        idx= where(gas_ids EQ idtofollow)
        if idx(0) EQ -1 or n_elements(idx) NE 1 then begin
                print, "Hey, what's up with idx?  printing"
                print, idx
                return, [0,0,0]
        endif

        x= xgas(idx)
        y= ygas(idx)
        z= zgas(idx)
        vx= vxgas(idx)-vcom(0)
        vy= vygas(idx)-vcom(1)
        vz= vzgas(idx)-vcom(2)
        m= mgas(idx)
    endif


    CASE typeofam OF

        0: BEGIN
           ; total angular momentum
           ; ----------------------
           jx= mgas*((ygas-c(1))*(vzgas-vcom(2)) - (zgas-c(2))*(vygas-vcom(1)))
           jy= -mgas*((xgas-c(0))*(vzgas-vcom(2)) - (zgas-c(2))*(vxgas-vcom(0)))
           jz= mgas*((xgas-c(0))*(vygas-vcom(1)) - (ygas-c(1))*(vxgas-vcom(0)))

           jx= total(jx)
           jy= total(jy)
           jz= total(jz)

           jtot= sqrt(jx*jx + jy*jy + jz*jz)
           return, total(jtot)
        END

        1: BEGIN
           jx= mgas*((ygas-c(1))*(vzgas-vcom(2)) - (zgas-c(2))*(vygas-vcom(1)))
           return, total(jx)
        END

        2: BEGIN 
           jy= - mgas*((xgas-c(0))*(vzgas-vcom(2)) - (zgas-c(2))*(vxgas-vcom(0)))
           return, total(jy)
        END

        3: BEGIN
           jz= mgas*((xgas-c(0))*(vygas-vcom(1)) - (ygas-c(1))*(vxgas-vcom(0)))
           return, total(jz)
        END

        44: BEGIN
           ; specific total angular momentum
           ; -------------------------------
           jx= mgas*((ygas-c(1))*(vzgas-vcom(2)) - (zgas-c(2))*(vygas-vcom(1)))
           jy= -mgas*((xgas-c(0))*(vzgas-vcom(2)) - (zgas-c(2))*(vxgas-vcom(0)))
           jz= mgas*((xgas-c(0))*(vygas-vcom(1)) - (ygas-c(1))*(vxgas-vcom(0)))

           jx= total(jx)
           jy= total(jy)
           jz= total(jz)
           tm= total(mgas)

           jtot= sqrt(jx*jx + jy*jy + jz*jz)
           return, total(jtot)/tm
        END

        10: BEGIN
           ; total angular momentum
           ; ----------------------
           jx= mgas*((ygas-c(1))*(vzgas-vcom(2)) - (zgas-c(2))*(vygas-vcom(1)))
           jy= -mgas*((xgas-c(0))*(vzgas-vcom(2)) - (zgas-c(2))*(vxgas-vcom(0)))
           jz= mgas*((xgas-c(0))*(vygas-vcom(1)) - (ygas-c(1))*(vxgas-vcom(0)))

           jtot= sqrt(jx*jx + jy*jy + jz*jz)
           return, jtot
        END

        11: BEGIN
           jx= mgas*((ygas-c(1))*(vzgas-vcom(2)) - (zgas-c(2))*(vygas-vcom(1)))
           return, jx
        END

        12: BEGIN
           jy= - mgas*((xgas-c(0))*(vzgas-vcom(2)) - (zgas-c(2))*(vxgas-vcom(0)))
           return, jy
        END

        13: BEGIN
           jz= mgas*((xgas-c(0))*(vygas-vcom(1)) - (ygas-c(1))*(vxgas-vcom(0)))
           return, jz
        END

        ;   specific angular momentum
        ; ------------------------------
        21: BEGIN
           jx= ((ygas-c(1))*(vzgas-vcom(2)) - (zgas-c(2))*(vygas-vcom(1)))
           return, jx
        END

        22: BEGIN
           jy= - ((xgas-c(0))*(vzgas-vcom(2)) - (zgas-c(2))*(vxgas-vcom(0)))
           return, jy
        END

        23: BEGIN
           jz= ((xgas-c(0))*(vygas-vcom(1)) - (ygas-c(1))*(vxgas-vcom(0)))
           return, jz
        END



        88: BEGIN
           ; specific total angular momentum
           ; -------------------------------
           jx= mgas*((ygas-c(1))*(vzgas-vcom(2)) - (zgas-c(2))*(vygas-vcom(1)))
           jy= -mgas*((xgas-c(0))*(vzgas-vcom(2)) - (zgas-c(2))*(vxgas-vcom(0)))
           jz= mgas*((xgas-c(0))*(vygas-vcom(1)) - (ygas-c(1))*(vxgas-vcom(0)))

           tm= total(mgas)

           jtot= sqrt(jx*jx + jy*jy + jz*jz)
           return, jtot/tm
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


