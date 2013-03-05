;------------------------------------------------------------------------
;
;    Determine Local Group Gas Map
;
;
;
;
;------------------------------------------------------------------------









;------------------------------------------------------------------------
;
;     Density Map from Earth
;     ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro lg_earthmap, frun, snapnum


frun= "/raid4/tcox/localgroup/v3"
snapnum= 28

if not keyword_set(frun) then begin
   print, "  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    "
   print, "  "
   print, "  "
   print, "     Need to use hidl for this procedure"
   print, "  "
   print, "  "
   print, "  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    "
   print, "  "
   print, "lg_earthmap, frun, snapnum "
   print, "  "
   print, "  * filenames are set from within:"
   print, "              columnden.eps"
   print, "              temp.eps"
   print, "              thermalSZ.eps"
   print, "              kineticSZ.eps"
   print, "  "
   return
endif




; ------------------------------------------------------------------------

if not keyword_set(frun) then frun= "/raid4/tcox/localgroup/v3"
if not keyword_set(snapnum) then snapnum= 28


if (fload_snapshot_bh(frun, snapnum)) then begin
   print, "PROBLEM: opening ",frun,snapnum
   return
endif



; /raid4/tcox/localgroup/v2
g1_sid= 1L
g1_npart= 500000L
g2_sid= 500001L
g2_npart= 800000L

; what time is it?
Ttime= fload_time(1)

bhid= fload_blackhole_id(1)
print, "Blackhole ID's : ", bhid

if fload_npart(5) gt 1 then begin
	bh_center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid[0])
	bh_center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid[1])
	;bh_center1= fload_blackhole_xyz('xyz',idtofollow=bhid[0])
	;bh_center2= fload_blackhole_xyz('xyz',idtofollow=bhid[1])
	print, "center1= ", bh_center1
	print, "center2= ", bh_center2
endif else begin
	;center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
	bh_center1= fload_blackhole_xyz('xyz',idtofollow=bhid)
	bh_center2= bh_center1
	print, "center1= ", bh_center1
endelse




; Earth - follow one star particle
; ----------------------------------------
			    ; these angles are pre-rotation
			    ; ------------------------------
;star_to_follow= 488933L    ; both are right on top of each other
;star_to_follow= 490923L    ; -20, 10
;star_to_follow= 495138L    ; -30, 0
;star_to_follow= 484475L    ; -45, 5
;star_to_follow= 490252L    ; 60 degreee sep
;star_to_follow= 490609L    ; 150 degree sep
;star_to_follow= 487194L    ; -150 degree sep
;star_to_follow= 487816L    ; -60, 0
;star_to_follow= 487751L    ; -90, 0
;star_to_follow= 4152L      ; -145, 5
;star_to_follow= 487755L      ; -138, -2
;star_to_follow= 488771L      ; -130, -2
star_to_follow= 484414L      ; -110, -2
;star_to_follow= 4268L      ; -   - sounds good, but there are 2 particles with this ID
;star_to_follow= 2870L      ; -   - sounds good, but there are 2 particles with this ID

center1= fload_allstars_xyz('dummy',center=[0,0,0],idtofollow=star_to_follow)
center2= center1
center1velocity= fload_allstars_v('dummy',idtofollow=star_to_follow)
center2velocity= center1velocity

print, "star_to_follow= ", star_to_follow
print, "star position= ", center1
print, "star velocity= ", center1velocity
print, " "
print, " dist. between bh's and star"
print, " bh1-star= ", sqrt(total((bh_center1-center1)*(bh_center1-center1)))
print, " bh2-star= ", sqrt(total((bh_center2-center1)*(bh_center2-center1)))



;
;  Find MW/M31 angles
; ----------------------------
if fload_npart(5) eq 1 then begin
	print, "ONLY 1 Black Hole - calling it MW"
        MW_r= bh_center2 - center1
        print, "MW r_dir= ", MW_r
        R= sqrt(MW_r[0]*MW_r[0] + MW_r[1]*MW_r[1])
        MW_theta= (180.0 / !PI) * acos(MW_r[2]/R)
        MW_phi= (180.0 / !PI) * atan(MW_r[1],MW_r[0])
        if MW_theta le 90.0 then MW_theta= -1.0*(90.0-MW_theta) else MW_theta= (MW_theta-90.0)
        print, "MW theta= ", MW_theta
        print, "     phi= ", MW_phi
endif else begin
   if bhid[0] eq 1300002L then M31_r= bh_center1 - center2
   if bhid[0] eq 500001L then MW_r= bh_center1 - center2
   if bhid[1] eq 500001L then MW_r= bh_center2 - center2
   if bhid[1] eq 1300002L then M31_r= bh_center2 - center2

	print, "M31 r_dir= ", M31_r
	R= sqrt(M31_r[0]*M31_r[0] + M31_r[1]*M31_r[1])
	M31_theta= (180.0 / !PI) * acos(M31_r[2]/R)
	M31_phi= (180.0 / !PI) * atan(M31_r[1],M31_r[0])
	if M31_theta le 90.0 then M31_theta= 90.0 - M31_theta else M31_theta= -1.0*(M31_theta-90.0)
	print, "M31 theta= ", M31_theta
	print, "      phi= ", M31_phi

	print, "MW r_dir= ", MW_r
	R= sqrt(MW_r[0]*MW_r[0] + MW_r[1]*MW_r[1])
	MW_theta= (180.0 / !PI) * acos(MW_r[2]/R)
	MW_phi= (180.0 / !PI) * atan(MW_r[1],MW_r[0])
	if MW_theta le 90.0 then MW_theta= -1.0*(90.0-MW_theta) else MW_theta= (MW_theta-90.0)
	print, "MW theta= ", MW_theta
	print, "     phi= ", MW_phi
endelse


; load data
; ------------



	Ngas= fload_npart(0)
	;Ntheta= 8L
	;Nphi= 8L
	;Ntheta= 16L
	;Nphi= 16L
	;Ntheta= 32L
	;Nphi= 32L
	;Ntheta= 256L
	;Nphi= 256L
	;Ntheta= 1024L
	;Nphi= 1024L

	; new HealPix stuff
	Nside= 16L
	;Nside= 64L
	Ntheta= Nside
	Nphi= 12*Nside

	; snap=0
	;earth_x= -561.92
	;earth_y= -6.41    +   5.0    ; put earth off the BH a bit
	;earth_y= -6.41
	;earth_z= 0.0    + 10.0
	; snap=4
	;earth_x= -560.28
	;earth_y= -25.62    +   5.0    ; put earth off the BH a bit
	;earth_y= -25.62    +   15.0    ; put earth off the BH a bit
	;earth_y= -25.62
	;earth_z= 0.0  ; + 10.0

	;earth_x= center2[0]
	;earth_y= center2[1]
	;earth_z= center2[2]
	earth_x= 0.0
	earth_y= 0.0
	earth_z= 0.0

        ; fields to pass
        Coord = fltarr(13,Ngas)

        ; gas
        ;Coord(0,0:Ngas-1) = fload_gas_xyz('x', center=[0,0,0])
        ;Coord(1,0:Ngas-1) = fload_gas_xyz('y', center=[0,0,0])
        ;Coord(2,0:Ngas-1) = fload_gas_xyz('z', center=[0,0,0])
        x = fload_gas_xyz('x', center=center2)
        y = fload_gas_xyz('y', center=center2)
        z = fload_gas_xyz('z', center=center2)
        vx = fload_gas_v('x', comvel=center2velocity)
        vy = fload_gas_v('y', comvel=center2velocity)
        vz = fload_gas_v('z', comvel=center2velocity)

	x_new= x
        y_new= y
        z_new= z
	vx_new= vx
        vy_new= vy
        vz_new= vz
	; ------------
	; rotate so that MW is at 0,0
	;
	rotate_MW_to00 = 1
	;rotate_MW_to00 = 0
	if rotate_MW_to00 eq 1 then begin
	    MW_theta= (!PI/180.0)*MW_theta
	    MW_phi=(!PI/180.0)*MW_phi
	    x_new= x*(cos(MW_theta)*cos(MW_phi)) + y*(cos(MW_theta)*sin(MW_phi)) - z*sin(MW_theta)
            y_new= -x*sin(MW_phi) + y*cos(MW_phi)
            z_new= x*(sin(MW_theta)*cos(MW_phi)) + y*(sin(MW_theta)*sin(MW_phi)) + z*cos(MW_theta)
	    vx_new= vx*(cos(MW_theta)*cos(MW_phi)) + vy*(cos(MW_theta)*sin(MW_phi)) - vz*sin(MW_theta)
            vy_new= -vx*sin(MW_phi) + vy*cos(MW_phi) 
            vz_new= vx*(sin(MW_theta)*cos(MW_phi)) + vy*(sin(MW_theta)*sin(MW_phi)) + vz*cos(MW_theta)
	    ;
	    if fload_npart(5) gt 1 then begin
		newM31_r= M31_r
		newM31_r[0]= M31_r[0]*(cos(MW_theta)*cos(MW_phi)) + M31_r[1]*(cos(MW_theta)*sin(MW_phi)) - M31_r[2]*sin(MW_theta)
		newM31_r[1]= -M31_r[0]*sin(MW_phi) + M31_r[1]*cos(MW_phi)
		newM31_r[2]= M31_r[0]*(sin(MW_theta)*cos(MW_phi)) + M31_r[1]*(sin(MW_theta)*sin(MW_phi)) + M31_r[2]*cos(MW_theta)
		R= sqrt(newM31_r[0]*newM31_r[0] + newM31_r[1]*newM31_r[1])
		M31_theta= (180.0 / !PI) * acos(newM31_r[2]/R)
		M31_phi= (180.0 / !PI) * atan(newM31_r[1],newM31_r[0])
		if M31_theta le 90.0 then M31_theta= 90.0 - M31_theta else M31_theta= -1.0*(M31_theta-90.0)
		print, "new M31 theta,phi= ",M31_theta, M31_phi
	    endif
	    ;
	    newMW_r= MW_r
	    newMW_r[0]= MW_r[0]*(cos(MW_theta)*cos(MW_phi)) + MW_r[1]*(cos(MW_theta)*sin(MW_phi)) - MW_r[2]*sin(MW_theta)
	    newMW_r[1]= -MW_r[0]*sin(MW_phi) + MW_r[1]*cos(MW_phi)
	    newMW_r[2]= MW_r[0]*(sin(MW_theta)*cos(MW_phi)) + MW_r[1]*(sin(MW_theta)*sin(MW_phi)) + MW_r[2]*cos(MW_theta)
	    R= sqrt(newMW_r[0]*newMW_r[0] + newMW_r[1]*newMW_r[1])
	    MW_theta= (180.0 / !PI) * acos(newMW_r[2]/R)
	    MW_phi= (180.0 / !PI) * atan(newMW_r[1],newMW_r[0])
	    if MW_theta le 90.0 then MW_theta= 90.0 - MW_theta else MW_theta= -1.0*(MW_theta-90.0)
	    print, "new MW theta,phi= ", MW_theta, MW_phi
	endif



	; rotate non-MW particles so that M31 is at -20
	; ---------------------------------------------
	rotate_non_MW= 1
	;rotate_non_MW= 0
	if rotate_non_MW eq 1 then begin
	    r_MW= sqrt(x_new*x_new + y_new*y_new + z_new*z_new)
	    idx= where(r_MW gt 200.0)

	    rottheta= -1.0*(!PI/180.0)* 20.0
            rotphi=0.0
            x_new(idx)= x(idx)*(cos(rottheta)*cos(rotphi)) + y(idx)*(cos(rottheta)*sin(rotphi)) - z(idx)*sin(rottheta)
            y_new(idx)= -x(idx)*sin(rotphi) + y(idx)*cos(rotphi)
            z_new(idx)= x(idx)*(sin(rottheta)*cos(rotphi)) + y(idx)*(sin(rottheta)*sin(rotphi)) + z(idx)*cos(rottheta)
            vx_new(idx)= vx(idx)*(cos(rottheta)*cos(rotphi)) + vy(idx)*(cos(rottheta)*sin(rotphi)) - vz(idx)*sin(rottheta)
            vy_new(idx)= -vx(idx)*sin(rotphi) + vy(idx)*cos(rotphi)
            vz_new(idx)= vx(idx)*(sin(rottheta)*cos(rotphi)) + vy(idx)*(sin(rottheta)*sin(rotphi)) + vz(idx)*cos(rottheta)
            ;
            if fload_npart(5) gt 1 then begin
                newM31_r= M31_r
                newM31_r[0]= M31_r[0]*(cos(rottheta)*cos(rotphi)) + M31_r[1]*(cos(rottheta)*sin(rotphi)) - M31_r[2]*sin(rottheta)
                newM31_r[1]= -M31_r[0]*sin(rotphi) + M31_r[1]*cos(rotphi)
                newM31_r[2]= M31_r[0]*(sin(rottheta)*cos(rotphi)) + M31_r[1]*(sin(rottheta)*sin(rotphi)) + M31_r[2]*cos(rottheta)
                R= sqrt(newM31_r[0]*newM31_r[0] + newM31_r[1]*newM31_r[1])
                M31_theta= (180.0 / !PI) * acos(newM31_r[2]/R)
                M31_phi= (180.0 / !PI) * atan(newM31_r[1],newM31_r[0])
		print, "pre-new M31 theta,phi= ", M31_theta, M31_phi
                if M31_theta le 90.0 then M31_theta= 90.0 - M31_theta else M31_theta= -1.0*(M31_theta-90.0)
                print, "new M31 theta,phi= ",M31_theta, M31_phi
		M31_theta= -1.0*M31_theta
                print, "new M31 theta,phi= ",M31_theta, M31_phi
            endif
	endif


	; boost the entire group with respect to the CMB
	; ----------------------------------------------
	group_v_boost= 1
	;group_v_boost= 0
	if group_v_boost eq 1 then begin
	   boost_v_tot= 300.0
	   boost_theta= (90.0 - 48.25) * !PI / 180.0
	   boost_phi= (263.85 - 180.0) * !PI / 180.0
	   boost_v_x= boost_v_tot * sin(boost_theta) * cos(boost_phi)
	   boost_v_y= boost_v_tot * sin(boost_theta) * sin(boost_phi)
	   boost_v_z= boost_v_tot * cos(boost_theta)
	   print, "Add velocity boost: ", boost_v_tot
	   print, "   boost_x= ", boost_v_x
	   print, "   boost_y= ", boost_v_y
	   print, "   boost_z= ", boost_v_z
print, 'old 0: ', vx_new[0], vy_new[0], vz_new[0]
print, 'old 1: ', vx_new[1], vy_new[1], vz_new[1]
print, 'old 2: ', vx_new[2], vy_new[2], vz_new[2]
	   vx_new= vx_new + boost_v_x
	   vy_new= vy_new + boost_v_y
	   vz_new= vz_new + boost_v_z
	endif



        Coord(0,0:Ngas-1) = x_new
        Coord(1,0:Ngas-1) = y_new
        Coord(2,0:Ngas-1) = z_new
        Coord(3,0:Ngas-1) = vx_new
        Coord(4,0:Ngas-1) = vy_new
        Coord(5,0:Ngas-1) = vz_new
print, '0: ', Coord[3,0], Coord[4,0], Coord[5,0]
print, '1: ', Coord[3,1], Coord[4,1], Coord[5,1]
print, '2: ', Coord[3,2], Coord[4,2], Coord[5,2]
        Coord(6,0:Ngas-1) = fload_gas_mass(1)
        Coord(7,0:Ngas-1) = fload_gas_u(1)
        Coord(8,0:Ngas-1) = fload_gas_rho(1)
        Coord(9,0:Ngas-1) = fload_gas_hsml(1)
        Coord(10,0:Ngas-1) = fload_gas_numh(1)
        Coord(11,0:Ngas-1) = fload_gas_nume(1)
        Coord(12,0:Ngas-1) = fload_gas_metallicity(1)


        los_NH= fltarr(Ntheta*Nphi)
        los_Z= fltarr(Ntheta*Nphi)
        los_Temp= fltarr(Ntheta*Nphi)
        los_ThermalSZ= fltarr(Ntheta*Nphi)
        los_KineticSZ= fltarr(Ntheta*Nphi)
        los_xrayL= fltarr(Ntheta*Nphi)

        print, "PASSING: "
        print, "Ngas= ", Ngas
        print, "Ntheta= ", Ntheta
        print, "Nphi= ", Nphi
        print, "earth_x= ", earth_x 
        print, "earth_y= ", earth_y
        print, "earth_z= ", earth_z
        help, Coord


	; Call our All-in-one C-program
	; ----------------------------------
        ;S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/EarthLOS_GasProp/getnh.so', $
        S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/testEarthLOS/getnh.so', $
                'getnh', $
                Ngas, $
                Ntheta, $
                Nphi, $
                earth_x, $
                earth_y, $
                earth_z, $
                Coord, $
                los_NH, $
                los_Z, $
		los_Temp, $
		los_ThermalSZ, $
		los_KineticSZ, $
		los_xrayL)

	    print, "max/min NH= ",max(los_NH), min(los_NH)
	    print, "max/min Z= ",max(los_Z), min(los_Z)
	    print, "max/min Temp= ",max(los_Temp), min(los_Temp)
	    print, "max/min ThermalSZ= ",max(los_ThermalSZ), min(los_ThermalSZ)
	    print, "max/min KineticSZ= ",max(los_KineticSZ), min(los_KineticSZ)
	    print, "max/min xrayL= ",max(los_xrayL), min(los_xrayL)


	rho= fload_gas_rho(1)
	u= fload_gas_u(1)
	idx=where(rho gt 0.000854924)
	if idx(0) ne -1 then u(idx)= 100.0
	Coord(7,0:Ngas-1)= u


        temp1= fltarr(Ntheta*Nphi)
        temp2= fltarr(Ntheta*Nphi)
        los_Temp_nodisk= fltarr(Ntheta*Nphi)
        los_ThermalSZ_nodisk= fltarr(Ntheta*Nphi)
        los_KineticSZ_nodisk= fltarr(Ntheta*Nphi)
        los_xrayL_nodisk= fltarr(Ntheta*Nphi)


	; Now, call program again, only without disk hot gas
        ;S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/EarthLOS_GasProp/getnh.so', $
        S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/testEarthLOS/getnh.so', $
                'getnh', $
                Ngas, $
                Ntheta, $
                Nphi, $
                earth_x, $
                earth_y, $
                earth_z, $
                Coord, $
                temp1, $
                temp2, $
                los_Temp_nodisk, $
                los_ThermalSZ_nodisk, $
                los_KineticSZ_nodisk, $
                los_xrayL_nodisk)

            print, "max/min Temp= ",max(los_Temp_nodisk), min(los_Temp_nodisk)
            print, "max/min ThermalSZ= ",max(los_ThermalSZ_nodisk), min(los_ThermalSZ_nodisk)
            print, "max/min KineticSZ= ",max(los_KineticSZ_nodisk), min(los_KineticSZ_nodisk)
            print, "max/min xrayL= ",max(los_xrayL_nodisk), min(los_xrayL_nodisk)



	;
	;  N_e
	; -------
	; trap for really low NH values
	idx= where(los_NH lt 1.0e+10)
	if idx(0) ne -1 then los_NH(idx)= 1.0e+10

	QuantityToProject= los_NH

	;process_data_to_epsfile, QuantityToProject, MW_phi, MW_theta, M31_phi, M31_theta, $
	;			Ntheta, Nphi, $
	;			cbtitlename= '!6Log N!DH!N (cm!E-2!N)', filename='columndensity.eps'


	mollview, QuantityToProject, ps='plot_nh.ps', /graticule, /online, title='N!De!N (cm!E2!N)'


	;
	;  Z
	; ------
	; trap for zero metallicity (make it small instead)
	idx= where(los_Z le 0.0)
	if idx(0) ne -1 then begin
	   other_idx= where(los_Z gt 0.0)
	   min_Z= min(los_Z(other_idx))
	   los_Z(idx)= min_Z/10.0
	   ;los_Z(idx)= 1.0e-5
	endif

	QuantityToProject= los_Z

	;process_data_to_epsfile, QuantityToProject, MW_phi, MW_theta, M31_phi, M31_theta, $
	;			Ntheta, Nphi, $
	;			cbtitlename= 'Log Z (Z!D!9n!6!N)', filename='metallicitymap.eps'

	mollview, QuantityToProject, ps='plot_z.ps', /graticule, /online, title='Metallicity (Z!D!9n!6!N)'


	;
	;  Temperature
	; --------------
	; Do deltaT/T
	;meanT= mean(los_Temp)
	;QuantityToProject= (los_Temp - meanT) - 1.0
        idx= where(los_Temp le 0.0)
        if idx(0) ne -1 then los_Temp(idx)= 1.0e+2

	QuantityToProject= los_Temp

	;process_data_to_epsfile, QuantityToProject, MW_phi, MW_theta, M31_phi, M31_theta, $
	;			Ntheta, Nphi, $
	;			cbtitlename= 'Log Temperature (K)', filename='temperature.eps'

	mollview, QuantityToProject, ps='plot_t.ps', /graticule, /online, title='Temperature (K)'



        ;
        ;  Temperature (no disk)
        ; -----------------------
        ; Do deltaT/T
        ;meanT= mean(los_Temp)
        ;QuantityToProject= (los_Temp - meanT) - 1.0
        idx= where(los_Temp le 0.0)
        if idx(0) ne -1 then los_Temp(idx)= 1.0e+2

        QuantityToProject= los_Temp_nodisk

        ;process_data_to_epsfile, QuantityToProject, MW_phi, MW_theta, M31_phi, M31_theta, $
        ;                        Ntheta, Nphi, $
        ;                        cbtitlename= 'Log Temperature (K)', filename='temperature_nd.eps'


	mollview, QuantityToProject, ps='plot_t_nd.ps', /graticule, /online, title='Temperature - no disk (K)'



	;
	;  Thermal SZ
	; --------------
	c= 2.9979e+10
	m_e= 9.10953e-28
	k= 1.3806e-16
	thompson = 6.65245d-25
	thermalSZ_factor= k*thompson/(m_e * c * c)
	print, "thermalSZ_factor= ",thermalSZ_factor

	QuantityToProject= 2.0 * thermalSZ_factor * los_ThermalSZ
	print, "Thermal SZ max/min = ", max(QuantityToProject), min(QuantityToProject)

	;process_data_to_epsfile, QuantityToProject, MW_phi, MW_theta, M31_phi, M31_theta, $
	;			Ntheta, Nphi, $
	;			Cl=Cl, $
	;			cbtitlename= 'Log !7D!6T / T', filename='thermalSZ.eps'

	;cl_t= Cl		; save the thermal Cl's

	mollview, QuantityToProject, ps='plot_tsz.ps', /graticule, /online, title='Thermal SZ'
	write_fits_map, 'lg_tsz.fits', QuantityToProject, coordsys='G', ordering='ring'




        ;
        ;  Thermal SZ  (without the disk)
        ; --------------------------------
        c= 2.9979e+10
        m_e= 9.10953e-28
        k= 1.3806e-16
        thompson = 6.65245d-25 
        thermalSZ_factor= k*thompson/(m_e * c * c)
        ;print, "thermalSZ_factor= ",thermalSZ_factor

        QuantityToProject= 2.0 * thermalSZ_factor * los_ThermalSZ_nodisk
        print, "Thermal SZ (no disk) max/min = ", max(QuantityToProject), min(QuantityToProject)

        ;process_data_to_epsfile, QuantityToProject, MW_phi, MW_theta, M31_phi, M31_theta, $
        ;                        Ntheta, Nphi, $
        ;                        Cl=Cl, $
        ;                        cbtitlename= 'Log !7D!6T / T', filename='thermalSZ_nd.eps'

        ;cl_t_nodisk= Cl                ; save the thermal Cl's

	mollview, QuantityToProject, ps='plot_tsz_nd.ps', /graticule, /online, title='Thermal SZ (no disk)'
	write_fits_map, 'lg_tsz_nd.fits', QuantityToProject, coordsys='G', ordering='ring'



	;
	;  Kinetic SZ
	; --------------
	kineticSZ_factor= thompson/c
	print, "kineticSZ_factor= ",kineticSZ_factor

	QuantityToProject= kineticSZ_factor * los_KineticSZ 
	print, "Kinetic SZ max/min = ", max(QuantityToProject), min(QuantityToProject)

	; units
	;QuantityToProject= QuantityToProject * 1.0e8

	;process_data_to_epsfile, QuantityToProject, MW_phi, MW_theta, M31_phi, M31_theta, $
	;			Ntheta, Nphi, $
	;			Cl=Cl, $
	;			cbtitlename= '10!E-8!N !7D!6T / T', filename='kineticSZ.eps'

	;cl_k= Cl * 1.0d-16		; Cl's for kinetic

	mollview, QuantityToProject, ps='plot_ksz.ps', /graticule, /online, title='Kinetic SZ'
	write_fits_map, 'lg_ksz.fits', QuantityToProject, coordsys='G', ordering='ring'




        ;
        ;  X-ray Luminosity (no disk)
        ; ----------------------------
        xrayL_factor= 1.0 * 1.989d+33
        ;print, "xrayL_factor= ", xrayL_factor

        QuantityToProject= xrayL_factor * los_xrayL_nodisk
	QuantityToProject= alog10(QuantityToProject)
        print, "xray Luminosity (no disk) max/min = ", max(QuantityToProject), min(QuantityToProject)

        mollview, QuantityToProject, ps='plot_xrayl_nd.ps', /graticule, /online, title='Xray Background (no disk)'
        ;write_fits_map, 'lg_xrayl.fits', QuantityToProject, coordsys='G', ordering='ring'


        ;
        ;  X-ray Luminosity 
        ; -------------------
        xrayL_factor= 1.0 * 1.989d+33
        ;print, "xrayL_factor= ", xrayL_factor

        QuantityToProject= xrayL_factor * los_xrayL
	QuantityToProject= alog10(QuantityToProject)
        print, "xray Luminosity max/min = ", max(QuantityToProject), min(QuantityToProject)

        mollview, QuantityToProject, ps='plot_xrayl.ps', /graticule, /online, title='Xray Background'
        ;write_fits_map, 'lg_xrayl.fits', QuantityToProject, coordsys='G', ordering='ring'




	; 
	;  Print C_l's 
	; --------------
	;process_cls, cl_t, cl_k, cl_t_nodisk



end







;
;
;=====================================================================








;
;
;=====================================================================

pro process_data_to_epsfile, QuantityToProject, MW_phi, MW_theta, M31_phi, M31_theta, $
				Ntheta, Nphi, $
				Cl=Cl, $
				cbtitlename= cbtitlename, filename=filename



initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize= 18
;setup_plot_stuff, 'ps', filename=filename, colortable=4



	; the los quantities _NH and los_Z arrays are returned such
	; that los_NH[theta,phi] but it's in one long
	; string, hence we need to deconvolve it so
	; that matches our expected data array
	data= fltarr(Ntheta,Nphi)
	idx= 0L
	for i=1,Ntheta do begin
	  for j=1,Nphi do begin
		lon= j-1
		lat= i-1
		;lon= i-1
		;lat= j-1
		data[lon,lat]= QuantityToProject[idx]
		idx=idx+1
	  endfor
	endfor



; rotate such that bh1 is at theta,phi of 0,0
;rotated_data= data
;rotated_data[0:10,*]= data[21:31,*]
;rotated_data[11:31,*]= data[0:20,*]
;data= rotated_data     ; need to rotate MW/M31 angles too







; process the moments of data
; -------------------------------
;calc_moments= 0
calc_moments= 1
if calc_moments eq 1 then begin

	area= 0.0
	dtheta= !PI/Ntheta
	dphi= 2.*!PI/Nphi

	l_max= 100


	
	;a00= 0.0
	;a10= 0.0 & a11= 0.0 & a1n1= 0.0
	;a20= 0.0 & a22= 0.0 & a21= 0.0 & a2n1= 0.0 & a2n2= 0.0
	;a30= 0.0
	;a40= 0.0
	;a50= 0.0
	;a60= 0.0
	;a70= 0.0
	;a80= 0.0
	a0= complexarr(1)
	a1= complexarr(3)
	a2= complexarr(5)
	a3= complexarr(7)
	a4= complexarr(9)
	a5= complexarr(11)
	a6= complexarr(13)
	a7= complexarr(15)
	a8= complexarr(17)

	moment_data= data
	;moment_data= data-mean(data)
	;moment_data= data-mean(data)-0.636556   ; substract mean of data and sin(theta)
	;moment_data= (data-mean(data))/mean(data)

	for j=1,Ntheta-1 do begin
	   for i=1,Nphi-1 do begin
	;for i=0,255 do begin
	   ;for j=0,255 do begin
		;
		; assume long, lat are as above
		;   N= 256
		;  
		;theta=  !PI*(j/255.0)   ; in radians
		;phi=  2.0*!PI*(i/255.0) - !PI     ; in radians
		theta=  !PI*(float(j)/Ntheta)   ; in radians
		phi=  2.0*!PI*(float(i)/Nphi) - !PI     ; in radians

		area= area + sin(theta)*dtheta*dphi

		a0[0]= a0[0] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,0,0,/double)

                a1[0]= a1[0] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,1,0,/double)
                a1[1]= a1[1] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,1,1,/double)
                a1[2]= a1[2] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,1,-1,/double)

                a2[0]= a2[0] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,2,0,/double)
                a2[1]= a2[1] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,2,2,/double)
                a2[2]= a2[2] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,2,1,/double)
                a2[3]= a2[3] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,2,-1,/double)
                a2[4]= a2[4] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,2,-2,/double)


		a3[0]= a3[0] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,3,3,/double)
		a3[1]= a3[1] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,3,2,/double)
		a3[2]= a3[2] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,3,1,/double)
		a3[3]= a3[3] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,3,0,/double)
		a3[4]= a3[4] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,3,-1,/double)
		a3[5]= a3[5] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,3,-2,/double)
		a3[6]= a3[6] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,3,-3,/double)

		a4[0]= a4[0] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,4,4,/double)
		a4[1]= a4[1] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,4,3,/double)
		a4[2]= a4[2] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,4,2,/double)
		a4[3]= a4[3] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,4,1,/double)
		a4[4]= a4[4] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,4,0,/double)
		a4[5]= a4[5] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,4,-1,/double)
		a4[6]= a4[6] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,4,-2,/double)
		a4[7]= a4[7] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,4,-3,/double)
		a4[8]= a4[8] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,4,-4,/double)

		a5[0]= a5[0] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,5,5,/double)
		a5[1]= a5[1] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,5,4,/double)
		a5[2]= a5[2] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,5,3,/double)
		a5[3]= a5[3] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,5,2,/double)
		a5[4]= a5[4] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,5,1,/double)
		a5[5]= a5[5] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,5,0,/double)
		a5[6]= a5[6] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,5,-1,/double)
		a5[7]= a5[7] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,5,-2,/double)
		a5[8]= a5[8] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,5,-3,/double)
		a5[9]= a5[9] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,5,-4,/double)
		a5[10]= a5[10] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,5,-5,/double)

		a6[0]= a6[0] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,6,/double)
		a6[1]= a6[1] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,5,/double)
		a6[2]= a6[2] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,4,/double)
		a6[3]= a6[3] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,3,/double)
		a6[4]= a6[4] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,2,/double)
		a6[5]= a6[5] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,1,/double)
		a6[6]= a6[6] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,0,/double)
		a6[7]= a6[7] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,-1,/double)
		a6[8]= a6[8] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,-2,/double)
		a6[9]= a6[9] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,-3,/double)
		a6[10]= a6[10] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,-4,/double)
		a6[11]= a6[11] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,-5,/double)
		a6[12]= a6[12] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,6,-6,/double)

		a7[0]= a7[0] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,7,/double)
		a7[1]= a7[1] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,6,/double)
		a7[2]= a7[2] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,5,/double)
		a7[3]= a7[3] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,4,/double)
		a7[4]= a7[4] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,3,/double)
		a7[5]= a7[5] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,2,/double)
		a7[6]= a7[6] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,1,/double)
		a7[7]= a7[7] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,0,/double)
		a7[8]= a7[8] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,-1,/double)
		a7[9]= a7[9] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,-2,/double)
		a7[10]= a7[10] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,-3,/double)
		a7[11]= a7[11] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,-4,/double)
		a7[12]= a7[12] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,-5,/double)
		a7[13]= a7[13] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,-6,/double)
		a7[14]= a7[14] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,7,-7,/double)

		a8[0]= a8[0] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,8,/double)
		a8[1]= a8[1] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,7,/double)
		a8[2]= a8[2] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,6,/double)
		a8[3]= a8[3] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,5,/double)
		a8[4]= a8[4] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,4,/double)
		a8[5]= a8[5] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,3,/double)
		a8[6]= a8[6] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,2,/double)
		a8[7]= a8[7] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,1,/double)
		a8[8]= a8[8] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,0,/double)
		a8[9]= a8[9] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,-1,/double)
		a8[10]= a8[10] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,-2,/double)
		a8[11]= a8[11] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,-3,/double)
		a8[12]= a8[12] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,-4,/double)
		a8[13]= a8[13] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,-5,/double)
		a8[14]= a8[14] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,-6,/double)
		a8[15]= a8[15] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,-7,/double)
		a8[16]= a8[16] + moment_data[i,j]*dtheta*dphi*sin(theta)*spher_harm(theta,phi,8,-8,/double)

	   endfor
	endfor

        c0= (1.0/(4.0*!PI))*(1./(2.*0. + 1.)) * total(a0*conj(a0))
        c1= (1.0/(4.0*!PI))*(1./(2.*1. + 1.)) * total(a1*conj(a1))
        c2= (1.0/(4.0*!PI))*(1./(2.*2. + 1.)) * total(a2*conj(a2))
        c3= (1.0/(4.0*!PI))*(1./(2.*3. + 1.)) * total(a3*conj(a3))
        c4= (1.0/(4.0*!PI))*(1./(2.*4. + 1.)) * total(a4*conj(a4))
        c5= (1.0/(4.0*!PI))*(1./(2.*5. + 1.)) * total(a5*conj(a5))
        c6= (1.0/(4.0*!PI))*(1./(2.*6. + 1.)) * total(a6*conj(a6))
        c7= (1.0/(4.0*!PI))*(1./(2.*7. + 1.)) * total(a7*conj(a7))
        c8= (1.0/(4.0*!PI))*(1./(2.*8. + 1.)) * total(a8*conj(a8))

print, "------------------------------------"
print, "area= ", area, " comapre to: ",4.0*!PI
print, "  "
print, "monopole, c0 =", c0
print, "dipole, c1 =", c1
print, "quadrupole, c2 =", c2
print, "c3=", c3
print, "c4=", c4
print, "c5=", c5
print, "c6=", c6
print, "c7=", c7
print, "c8=", c8
print, "------------------------------------"


Cl= [c0, c1, c2, c3, c4, c5, c6, c7, c8]

endif



; rescale data to 0-255
; ------------------------

;	data[longitude,latitude]
help, data


print, "data max/min= ", max(data), min(data)

linearplot= 0
if min(data) lt 0.0 then begin
	print, "data is negative: adjusting"
	;data= data - min(data) + 0.001
	linearplot= 1
endif

ma= max(data)
mi= min(data)
;mi= 0.0
;ma= 1.0e+22
;mi= 1.0e+20

print, "clipping at max/min= ", ma, mi

idx= where(data gt ma)
if idx(0) ne -1 then data(idx)= ma
idx= where(data lt mi)
if idx(0) ne -1 then data(idx)= mi

ncolors= 255

; log scaling
Pic= (alog10(data)-alog10(mi))/alog10(ma/mi) * (ncolors-3) + 2
; linear scaling
if linearplot eq 1 then Pic= (data-mi)/(ma-mi) * (ncolors-3) + 2

idx= where(Pic ge 256)
if idx(0) ne -1 then Pic(idx)=255
idx=where(Pic le 0)
if idx(0) ne -1 then Pic(idx)=1


;invertcolors= 1
;if invertcolors eq 1 then begin
;	Pic=256-Pic                          ; invert color table
;	idx= where(Pic EQ 254)               ; set background to white
;endif else begin
;	idx= where(Pic EQ 2)
;endelse

;if idx(0) ne -1 then Pic(idx)= 1       ; this is setting it to white
;if idx(0) ne -1 then Pic(idx)= 0       ; this sets background to black




help, Pic



; ------------------------------------------------------------------------





; setup the map projection
; ------------------------
x0= 0.0
x1= 0.89
y0= 0.0
y1= 1.0

!p.position= [x0, y0, x1, y1]

;map_set, MW_theta, MW_phi, $
map_set, 0, 0, $
;map_set, 0, 180, $       ; puts the center at 0, 180 (instead of the default 0,0)
;	projections
;	------------
;       /aitoff, $       ; kind of eggy
;       /mercator, $     ; square, cuts off some
;       /cylindrical, $  ; square and stretches at top
;       /miller, $       ; cuts off +/- 10 degrees
	/hammer, $       ; fits it all onto circle (similar to aitoff?)
;	/orthographic, $ ; circle, does only half
;	/mollweide, $    ; not that much different that aitoff
;
;	other
;	-----
;       limit=[latmin,lonmin,latmax,lonmax]
;       /isotropic, $
;       /grid, $
;       /label, $
        color= 0


;
; get image ready to overplot on the map
;
image=map_image(Pic,xx,yy,/bilinear, compress=1)
;image=map_image(Pic,xx,yy,/bilinear)
;image=map_image(Pic,xx,yy, $
;               compress=1, $
;               latmin=latmin, lonmin=lonmin, $
;               latmax=latmax, lonmax=lonmax)
;image= Pic
idx=where(image eq 0)
if idx(0) ne -1 then image(idx)= 1
help, image
help, xx, yy

;
; display warped image on the map, and do
; it at the proper position (xx,yy)
;
;tv, image, xx, yy
tv, image, x0, y0, xsize= (x1-x0), ysize= (y1-y0), /normal

; continents
; ------------
;map_continents, color=0
;map_continents, color=150, /fill

; grid
; -----
;map_grid
;map_grid, latdel=10, londel=10, /label, /horizon
map_grid, latdel=15, londel=30, /label, color= 1



; ----------------------------------------------------------------------



; overplot MW and M31 centers?
; -----------------------------
;plots, MW_phi, MW_theta, psym=7, color= 150, symsize=1.7
;plots, M31_phi, M31_theta, psym=7, color= 220, symsize=1.7

if fload_npart(5) gt 1 then begin
   xyouts, MW_phi, MW_theta, 'MW', color= 1, charthick=4.7, size=1.2
   xyouts, M31_phi, M31_theta, 'M31', color= 1, charthick=4.7, size=1.2
endif else begin
   xyouts, M31_phi, M31_theta, 'BH', color= 1, charthick=4.7, size=1.2
endelse



; ----------------------------------------------------------------------

; print random information ?

;xyouts, 0.72, 0.83, '40% gas', /normal, size=1.5, charthick=3.0, color= 0
;xyouts, 0.68, 0.83, 'collisionless', /normal, size=1.5, charthick=3.0, color= 0

xyouts, 0.03, 0.92, fload_timelbl(1,2), /normal, size= 1.5, color= 0, charthick=4.5



; ----------------------------------------------------------------------


colorbar= 1
if keyword_set(colorbar) then begin
	invertcolors= 0    ; set above
	;ma=MaxDensXY      ; these are set above
	;mi=MaxDensXY/DynRange

	bar= REPLICATE(1B, 10)#BINDGEN(256)
	if invertcolors eq 1 then bar= 255-bar
	idx= where(bar eq 0 or bar eq 1)
	if idx(0) ne -1 then begin
	   ;print, "idx has ", n_elements(idx), idx(0)
	   ;colornumber= idx(0)-1
	   ;if colornumber lt 0 then colornumber= 11   ; if it's bar[*,0] then set to bar[*,1]
	   bar(idx)= 2   ; bar(colornumber)
	endif

	barwidth= 0.03

	;tv, bar, x1+0.01, y0, xsize=barwidth, ysize=(y1-y0), /normal
	tv, bar, x1, y0+0.02, xsize=barwidth, ysize=(y1-y0-0.04), /normal

	!p.position= [x1,y0+0.02,x1+barwidth,y1-0.02]

	if linearplot then plotmax= ma else plotmax= alog10(ma)
	if linearplot then plotmin= mi else plotmin= alog10(mi)
	plot, [0],[0], psym=3, xrange=[0,1], yrange=[plotmin,plotmax], $
	   /normal, ticklen=0.3, charthick= 2.5, $
	   /noerase, color= 0, xstyle=1, ystyle=1, $
	   xthick=4.0, ythick=4.0, $
	   xticks=1, xtickformat='(a1)', ytickformat='(a1)', $
	   /nodata

	axis, yaxis=1, yrange=[plotmin,plotmax], ystyle=1, $
	   ;xcharsize=1.50, ycharsize=1.50, charthick=3.0, $
	   /normal, ticklen= 0.3, charthick= 2.5, $
	   ytitle= cbtitlename
endif



; ----------------------------------------------------------------------

device, /close



end

















; =========================================================================


pro test_data, junk

; testing purposes
; ------------------
do_test= 0
;do_test= 1
if do_test eq 1 then begin
	data= dist(256)
	;;data= bytscl(sin(dist(400)/10))
	for lon= 0,255 do begin
	  for lat= 0,255 do begin
		;ilat= max(lon,lat)
		;ilat= 250.0*sin(!PI*(lat/255.0))
		;if ilat gt 0.0 then ilat= sqrt(ilat)
		;if ilat gt 0.0 then ilat= ilat*ilat

		;ilat= cos(!PI*(lat/255.0))  + 0.2
		ilat= cos(!PI*(lat/255.0))
		ilat= ilat*ilat

		data[lon,lat]= ilat
	  endfor
	endfor
endif

; rescale to between 250 and 0
;data= 250.0*data/max(data)
;help, data




	QuantityToProject= data

	process_data_to_epsfile, QuantityToProject, MW_phi, MW_theta, M31_phi, M31_theta, $
				cbtitlename= 'Test', filename='testdata.eps'



end









;==========================================================











;
;
;------------------------------------------------------------------------

pro testmap, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "test, junk"
   print, "  "
   return
endif

filename='test.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;--------------------------------------
;--------------------------------------



; for us with /azim projection
;map_set, 32, -100, /azim, limit=[10,-130,55,-70], /grid, /cont, /label

; for city
map_set, /mercator, /grid, /continent, limit=[10,-130,60,-70]

; lat/lon of santa cruz
lat= [34.0]
lon= [-119.4]
city=['Santa Cruz, CA']
plots,lon,lat,psym=4,symsize=1.4,color=150
xyouts,lon,lat,city,color=50, charthick=2, charsize=1.5,align=0.5



; --------------------

; print random information ?

;xyouts, 0.72, 0.83, '40% gas', /normal, size=1.5, charthick=3.0, color= 0
;xyouts, 0.68, 0.83, 'collisionless', /normal, size=1.5, charthick=3.0, color= 0


; --------------------

device, /close



end







;==========================================================








;
;
;-------------------------------------------------------
pro determine_earth, junk


if not keyword_set(frun) then frun= "/raid4/tcox/localgroup/v3"
if not keyword_set(snapnum) then snapnum= 28


if (fload_snapshot_bh(frun, snapnum)) then begin
   print, "PROBLEM: opening ",frun,snapnum
   return
endif



; /raid4/tcox/localgroup/v2
g1_sid= 1L
g1_npart= 500000L
g2_sid= 500001L
g2_npart= 800000L

; what time is it?
Ttime= fload_time(1)

bhid= fload_blackhole_id(1) 
print, "Blackhole ID's : ", bhid

if fload_npart(5) gt 1 then begin
        bh_center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid[0])
        bh_center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid[1])
        ;bh_center1= fload_blackhole_xyz('xyz',idtofollow=bhid[0])
        ;bh_center2= fload_blackhole_xyz('xyz',idtofollow=bhid[1])
        print, "center1= ", bh_center1
        print, "center2= ", bh_center2
endif else begin
        ;center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
        bh_center1= fload_blackhole_xyz('xyz',idtofollow=bhid)
        bh_center2= bh_center1
        print, "center1= ", bh_center1
endelse




; Earth - which particle
; ----------------------------------------
r=fload_allstars_xyz('r',center=bh_center1)
x=fload_allstars_xyz('x',center=bh_center1)
y=fload_allstars_xyz('y',center=bh_center1)
z=fload_allstars_xyz('z',center=bh_center1)

id=fload_allstars_ids(1)

idx=where((r gt 5.95) and (r lt 6.0))
help, idx

stop
print, x(idx)
print, y(idx)
print, id(idx)


star_to_follow= 480422L
center1= fload_allstars_xyz('dummy',center=[0,0,0],idtofollow=star_to_follow)
center2= center1 
center1velocity= fload_allstars_v('dummy',idtofollow=star_to_follow)
center2velocity= center1velocity

print, "star_to_follow= ", star_to_follow
print, "star position= ", center1
print, "star velocity= ", center1velocity
print, " "
print, " dist. between bh's and star"
print, " bh1-star= ", sqrt((bh_center1-center1)*(bh_center1-center1))
print, " bh2-star= ", sqrt((bh_center2-center1)*(bh_center2-center1))



;
;  Find MW/M31 angles
; ----------------------------
M31_r= bh_center1 - center1
print, "M31 r_dir= ", M31_r 
R= sqrt(M31_r[0]*M31_r[0] + M31_r[1]*M31_r[1])
M31_theta= (180.0 / !PI) * acos(M31_r[2]/R)
M31_phi= (180.0 / !PI) * atan(M31_r[1],M31_r[0])
if M31_theta le 90.0 then M31_theta= 90.0 - M31_theta else M31_theta= -1.0*(M31_theta-90.0)
print, "M31 theta= ", M31_theta
print, "      phi= ", M31_phi


if fload_npart(5) gt 1 then begin
        MW_r= bh_center2 - center1
        print, "MW r_dir= ", MW_r 
        R= sqrt(MW_r[0]*MW_r[0] + MW_r[1]*MW_r[1])
        MW_theta= (180.0 / !PI) * acos(MW_r[2]/R)
        MW_phi= (180.0 / !PI) * atan(MW_r[1],MW_r[0])
        if MW_theta le 90.0 then MW_theta= -1.0*(90.0-MW_theta) else MW_theta= (MW_theta-90.0)
        print, "MW theta= ", MW_theta
        print, "     phi= ", MW_phi
endif
                                                                                                       


end






;==========================================================









;
;  Print Angular power spectrum
; ---------------------------------
pro process_cls, cl_t, cl_k, cl_t_nodisk


; Print thie mess up
; -------------------
filename='apowerspec.eps' 

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4

x0= 0.18 & x1= 0.99
y0= 0.15 & y1= 0.99

yaxistitle='!12l(l+1)!6 C!D!12l!6!N / 2!7p!6 (!7l!3K!E2!N)'
yaxistitle='(!12l(l+1)!6 C!D!12l!6!N / 2!7p!6)!E1/2!N (!7l!3K)'
;ymax= 1.0e-9
;ymin= 1.0e-18
ymax= 1.0e+4
ymin= 1.0e-3

xaxistitle='!12l!6'
xmax= 12.0
xmin= 0.0

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        ;xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, /ylog, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase


;lls= [0.1,1,2,3,4,5,6,7,8]
lls= [2,3,4,5,6,7,8]

; take out monopole and dipole
cl_t= cl_t[2:8]
cl_k= cl_k[2:8]
cl_t_nodisk= cl_t_nodisk[2:8]

cl_t= lls * (lls + 1.0) / (2.0 * !PI) * cl_t
cl_k= lls * (lls + 1.0) / (2.0 * !PI) * cl_k
cl_t_nodisk= lls * (lls + 1.0) / (2.0 * !PI) * cl_t_nodisk


; change to micro-K squared
cl_t = cl_t * 2.8 * 1.0e+12
cl_k = cl_k * 2.8 * 1.0e+12
cl_t_nodisk = cl_t_nodisk * 2.8 * 1.0e+12


cl_t= sqrt(cl_t)
cl_k= sqrt(cl_k)
cl_t_nodisk= sqrt(cl_t_nodisk)

oplot, lls, cl_t, psym=-2, color=50, thick= 3.0, symsize= 1.5
oplot, lls, cl_t_nodisk, psym=-2, color=100, thick= 3.0, symsize= 1.5
oplot, lls, cl_k, psym=-5, color=150, thick= 3.0, symsize= 1.5


xyouts, 0.38, 0.87, 'thermal', /normal, charthick=3.0, size=1.5, color=50
xyouts, 0.38, 0.83, 'thermal (remove disk)', /normal, charthick=3.0, size=1.5, color=100
xyouts, 0.38, 0.79, 'kinetic', /normal, charthick=3.0, size=1.5, color=150

overplot_WMAP_cls, 1, /justmicroK


; -------------
;  Done
; -------------

device, /close



end





;==========================================================




;
;  Print Angular power spectrum
; ---------------------------------
pro plot_powerspectrum, junk


filename='apowerspec.eps' 

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4

x0= 0.18 & x1= 0.99
y0= 0.15 & y1= 0.99

;yaxistitle='!12l(l+1)!6 C!D!12l!6!N / 2!7p!6 (!7l!3K!E2!N)'
;ymax= 1.0e-9
;ymin= 1.0e-18

yaxistitle='(!12l(l+1)!6 C!D!12l!6!N / 2!7p!6)!E1/2!N (!7l!3K)'
ymax= 1.0e+4
ymin= 1.0e-3

xaxistitle='!6Multipole moment (!12l!6)'
xmax= 80.0
xmin= 1.0

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, /ylog, $
        ;xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata


; --------------------

; kinetic
; --------
fits2cl, clar, "/home/tcox/CMBStuff/lg/lg_ksz_Cl.fits", multipoles=mps
power= 1.0/(2.*!PI) * clar[2:1000]*mps[2:1000]*(mps[2:1000] + 1)
lls= mps[2:1000]
power= sqrt(power)*1.0e+6   ; converts to microK
oplot, lls, power, psym=-3, color= 150


; thermal
; --------
fits2cl, clar, "/home/tcox/CMBStuff/lg/lg_tsz_Cl.fits", multipoles=mps
power= 1.0/(2.*!PI) * clar[2:1000]*mps[2:1000]*(mps[2:1000] + 1)
lls= mps[2:1000]
power= sqrt(power)*1.0e+6   ; converts to microK
oplot, lls, power, psym=-3, color= 50


; thermal (no disk)
; ------------------
fits2cl, clar, "/home/tcox/CMBStuff/lg/lg_tsz_nd_Cl.fits", multipoles=mps
power= 1.0/(2.*!PI) * clar[2:1000]*mps[2:1000]*(mps[2:1000] + 1)
lls= mps[2:1000]
power= sqrt(power)*1.0e+6   ; converts to microK
oplot, lls, power, psym=-3, color= 100


; --------------------

xyouts, 0.38, 0.87, 'thermal', /normal, charthick=3.0, size=1.5, color=50
xyouts, 0.38, 0.83, 'thermal (remove disk)', /normal, charthick=3.0, size=1.5, color=100
xyouts, 0.38, 0.79, 'kinetic', /normal, charthick=3.0, size=1.5, color=150

overplot_WMAP_cls, 1, /justmicroK


; Tegmark's cleaned map
fits2cl, clar, "/home/tcox/CMBStuff/test/qaz_anafast_Cl.fits", multipoles=mps
power= 1.0/(2.*!PI) * clar[2:1000]*mps[2:1000]*(mps[2:1000] + 1)
lls= mps[2:1000]
power= sqrt(power)*1.0e+3   ; converts to microK
oplot, lls, power, psym=-3, color=200


; -------------
;  Done
; -------------

device, /close



end






;==========================================================






;
;   Overplot the WMAP angular power spectrum
; -------------------------------------------
pro overplot_WMAP_cls, junk, justmicroK=justmicroK


wmapfile= '/home/tcox/CMBStuff/WMAP_TTangularpowerspec.txt'


;
;# Wilkinson Microwave Anisotropy Probe (WMAP) year 1 data release.
;# WMAP One-Year Binned Combined TT Power Spectrum, version 1.1 (Oct 2003)
;# Reference = WMAP Explanatory Supplement: http://lambda.gsfc.nasa.gov/
;# Column 1 = mean multipole moment l for the bin
;# Column 2 = smallest l contributing to the bin
;# Column 3 = largest l contributing to the bin
;# Column 4 = mean value of TT power spectrum (=  l(l+1)/2pi * C_l) in the bin,
;#             units = uK^2
;# Column 5 = 'Error' for binned value, as computed from diagonal terms of the
;#             Fisher matrix, units = uK^2.
;#             Included only as an estimate of their magnitude.  The
;#             multipole moments are slightly coupled, so a correct
;#             treatment of errors requires use of the entire Fisher matrix.
;# Column 6 = portion of column5 error attributed to measurement errors,
;#             units = uK^2.
;# Column 7 = portion of column5 error attributed to cosmic variance,
;#             assuming the best-fit running index model described in
;#             "First Year WMAP Results: Determination of Cosmological
;#              Parameters", Spergel et al. (2003). Units = uK^2.
;#
;# Revision History:
;#  v1  : initial release
;#  v1p1: inserted 2 new columns containing the lower and upper l bin boundaries
;#        after column 1 and re-numbered the columns.  Contents of v1
;#        columns remain unchanged.  31 Oct 2003.
;#
;    2     2     2    123.3820    762.6369       3.3813     759.2556
;    3     3     3    611.7750    608.1737       4.0262     604.1475
;    4     4     5   1006.6580    331.9451       3.1432     328.8178
;    6     6     7    763.1490    255.8872       3.4429     252.4549
;    9     8    11    828.1140    143.4565       2.9464     140.5284
;
; etc...
;

spawn, "wc "+wmapfile,result
lines=long(result)
datalines=lines(0)-26
meanl= fltarr(datalines)
minl= fltarr(datalines)
maxl= fltarr(datalines)
TT_Cl= fltarr(datalines)
TT_Cl_err= fltarr(datalines)

openr, 1, wmapfile
junk=''

; read the header
for i=0,25 do begin
	readf, 1, junk
endfor

; read the data
for i=0,lines(0)-27 do begin
        readf, 1, junk
        ;print, junk
        tempjunk= strsplit(junk,/extract,count=count)
        meanl(i)= float(tempjunk(0))
        minl(i)= float(tempjunk(1))  
        maxl(i)= float(tempjunk(2))  
        TT_Cl(i)= float(tempjunk(3))
        TT_Cl_err(i)= float(tempjunk(4))
endfor

close, 1


; make it just micro-K
if keyword_set(justmicroK) then begin
	TT_Cl= sqrt(TT_Cl)
endif


; overplot the spectrum
symsize= 1.0
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0

oplot, meanl, TT_Cl, psym=-8, color=0, thick= 3.0

xyouts, 0.38, 0.92, 'WMAP TT Power Spectrum', /normal, charthick=3.0, size=1.5, color=0



end









;==========================================================================



pro HealPIX_Stuff, junk



; must use hidl
read_fits_map, "CMBStuff/tegmark_cleaned_map.fits", thismap, hdr, exthdr
mollview, thismap, ps='plot.ps', /graticule, /online, max=0.2, min=-0.2


; downgrade the image to a lower resolution
;
ud_grade, thismap, outmap, nside_out=32, order_in='ring'



write_fits_map, 'CMBStuff/tegmark_cleaned_map_lowres.fits', outmap, coordsys='G', ordering='ring'




; --------------------
;  other random info
; --------------------
;
;  anafast - performs the fourier analysis, i.e., it calculates the C_l's by
;            integrating over the entire sphere (unless a mask is provided)
;            up to a specificed l.
;
;
;  to run anafast (in /home/tcox/CMBStuff/test):
;    ../anafast.scr 1000 ../tegmark_cleaned_map.fits test_cl.dat test_alm.fits
;
;
;
;  synfast - produces a map based upon a power specturm
;
;  ud_grade - can up/down grade a HealPix map
;



; idl procedure to read in the alm fits file generated by
; Matias' program anafast.scr
fits2alm, idx, almarray, "CMBStuff/test/test_alm.fits"

fits2cl, clar, "CMBStuff/test/qaz_anafast_Cl.fits", multipoles=mps
plot, mps[2:1000], clar[2:1000]*mps[2:1000]*(mps[2:1000] + 1), psym=-3, /xlog


end




pro gen_healpix_angles, nside

if not keyword_set(nside) then nside= 16

; generate the angles in the "Ring" orientation
;
; nside:   this is the number of rings
; pixels:  equal to 12.*nside*nside
;    pixel_number can be an array, or a
;    single pixel number
; theta:   angle, in radians, of the pixel (can be an array)
; phi:     angle, in radians, of the pixel (can be an array)
;
npixels= 12.*nside*nside
pixeln= findgen(npixels)

pix2ang_ring, nside, pixeln, theta, phi



; ----------------------------------------
; write to file

openw, 1, 'healpixangles.txt', ERROR=err

printf, 1, "#   healpixangles.txt"     
printf, 1, "#       Ring orientation, with nside= "+string(nside)
printf, 1, "#               "
printf, 1, "# pixle  theta       phi"
printf, 1, "# (Gyr)  (rad)     (rad)"
for i=0L,npixels-1 do begin
        printf, 1, FORMAT= '(F6.1," ",F8.5,"  ",F8.5)', $
                pixeln[i], theta[i], phi[i]
endfor
close, 1    


            
end












