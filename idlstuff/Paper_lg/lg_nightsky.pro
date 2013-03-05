;------------------------------------------------------------------------
;
;
;
;     Night Sky (or more accurately, a stellar mass Density Map)
;        from a candidate Earth
;     -----------------------------------------------------------
;     -----------------------------------------------------------
;
;
;
;------------------------------------------------------------------------

;
;  v1
;--------------------------
pro lg_nightsky_v1, junk

frun= "/raid4/tcox/localgroup/v7"
snapnum= 29
filename= 'nightsky_nh.eps'
imgname= 'nightsky_nh.png'

if not keyword_set(junk) then begin
   print, "  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    "
   print, "  "
   print, "  "
   print, "     Need to use hidl for this procedure"
   print, "  "
   print, "  "
   print, "  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    "
   print, "  "
   print, "  lg_nightsky, junk"
   print, "  "
   print, "  static variables:"
   print, "       frun= ", frun
   print, "    snapnum= ", snapnum
   print, "   filename= ", filename
   print, "  "
   print, "  "
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


;
; standard /raid4/tcox/localgroup/v* (2-19)
;
g1_sid= 1L
g1_npart= 500001L
mw_bhid= 500001L
g2_sid= 500002L
g2_npart= 800001L
m31_bhid= 1300002L


;frun= "/raid4/tcox/localgroup/bhires"
;mw_bhid= 817955
;m31_bhid= 2145299


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
;star_to_follow= 484414L      ; -110, -2
star_to_follow= 484400L      ; ?
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

	Nstars= fload_npart(2)+fload_npart(3)+fload_npart(4)

	; new HealPix stuff (read it in)
	healpix= '/home2/tcox/CMBStuff/HealPix_RingPixels/healpixangles.txt'
	openr, 1, healpix
	junk=''
	readf, 1, junk
	readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & Nside= long(tempjunk(5))
	close,1
	;Nside= 16L
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
        ;Coord = fltarr(13,Ngas)
        Coord = fltarr(13,Nstars)

        ; gas
        ;Coord(0,0:Ngas-1) = fload_gas_xyz('x', center=[0,0,0])
        ;Coord(1,0:Ngas-1) = fload_gas_xyz('y', center=[0,0,0])
        ;Coord(2,0:Ngas-1) = fload_gas_xyz('z', center=[0,0,0])
        ;x = fload_gas_xyz('x', center=center2)
        ;y = fload_gas_xyz('y', center=center2)
        ;z = fload_gas_xyz('z', center=center2)
        ;vx = fload_gas_v('x', comvel=center2velocity)
        ;vy = fload_gas_v('y', comvel=center2velocity)
        ;vz = fload_gas_v('z', comvel=center2velocity)

        ; stars
        x = fload_allstars_xyz('x', center=center2)
        y = fload_allstars_xyz('y', center=center2)
        z = fload_allstars_xyz('z', center=center2)
        vx = fload_allstars_v('x', comvel=center2velocity)
        vy = fload_allstars_v('y', comvel=center2velocity)
        vz = fload_allstars_v('z', comvel=center2velocity)

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
	;group_v_boost= 1
	group_v_boost= 0
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
	   vx_new= vx_new + boost_v_x
	   vy_new= vy_new + boost_v_y
	   vz_new= vz_new + boost_v_z
	endif



        Coord(0,0:Nstars-1) = x_new
        Coord(1,0:Nstars-1) = y_new
        Coord(2,0:Nstars-1) = z_new
        Coord(3,0:Nstars-1) = vx_new
        Coord(4,0:Nstars-1) = vy_new
        Coord(5,0:Nstars-1) = vz_new
print, '0: ', Coord[3,0], Coord[4,0], Coord[5,0]
print, '1: ', Coord[3,1], Coord[4,1], Coord[5,1]
print, '2: ', Coord[3,2], Coord[4,2], Coord[5,2]
        Coord(6,0:Nstars-1) = fload_allstars_mass(1)
        Coord(7,0:Nstars-1) = x_new*0.0 + 1.0
        Coord(8,0:Nstars-1) = x_new*0.0 + 1.0
        Coord(9,0:Nstars-1) = x_new*0.0 + 20.0   ; hsml
        Coord(10,0:Nstars-1) = x_new*0.0 + 1.0
        Coord(11,0:Nstars-1) = x_new*0.0 + 1.0
        Coord(12,0:Nstars-1) = x_new*0.0 + 1.0


        los_NH= fltarr(Ntheta*Nphi)

        print, "PASSING: "
        print, "Nstars= ", Nstars
        print, "Ntheta= ", Ntheta
        print, "Nphi= ", Nphi
        print, "earth_x= ", earth_x 
        print, "earth_y= ", earth_y
        print, "earth_z= ", earth_z
        help, Coord


	; Call our All-in-one C-program
	; ----------------------------------
        ;S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/EarthLOS_GasProp/getnh.so', $
        ;S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/testEarthLOS/getnh.so', $
        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/EarthLOS_stars/getnh.so', $
                'getnh', $
                Nstars, $
                Ntheta, $
                Nphi, $
                earth_x, $
                earth_y, $
                earth_z, $
                Coord, $
                los_NH)

	    print, "max/min NH= ",max(los_NH), min(los_NH)



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


	lg_mollview, QuantityToProject, ps=filename, /graticule, /online, title='N!De!N (cm!E2!N)', $
					png=imgname, /crop


end







;
;
;=====================================================================


pro doit, junk
  ;lg_nightsky_v2, 1, snapnum= 29, filename='nightsky_today_hr.eps', imgname='nightsky_today_hr.png'
;  lg_nightsky_v2, 1, snapnum= 46, filename='nightsky_first_hr.eps', imgname='nightsky_first_hr.png'
;  lg_nightsky_v2, 1, snapnum= 49, filename='nightsky_post1_hr.eps', imgname='nightsky_post1_hr.png'
;  lg_nightsky_v2, 1, snapnum= 57, filename='nightsky_second_hr.eps', imgname='nightsky_second_hr.png'
;  lg_nightsky_v2, 1, snapnum= 59, filename='nightsky_post2_hr.eps', imgname='nightsky_post2_hr.png'
;  lg_nightsky_v2, 1, snapnum= 62, filename='nightsky_final_hr.eps', imgname='nightsky_final_hr.png'
;  lg_nightsky_v2, 1, snapnum= 80, filename='nightsky_milkomeda_hr.eps', imgname='nightsky_milkomeda_hr.png'

;  lg_nightsky_v2, 1, snapnum= 29, filename='nightsky_today.eps', imgname='nightsky_today.png'
;  lg_nightsky_v2, 1, snapnum= 46, filename='nightsky_first.eps', imgname='nightsky_first.png'
;  lg_nightsky_v2, 1, snapnum= 49, filename='nightsky_post1.eps', imgname='nightsky_post1.png'
;  lg_nightsky_v2, 1, snapnum= 57, filename='nightsky_second.eps', imgname='nightsky_second.png'
;  lg_nightsky_v2, 1, snapnum= 59, filename='nightsky_post2.eps', imgname='nightsky_post2.png'
;  lg_nightsky_v2, 1, snapnum= 62, filename='nightsky_final.eps', imgname='nightsky_final.png'
;  lg_nightsky_v2, 1, snapnum= 80, filename='nightsky_milkomeda.eps', imgname='nightsky_milkomeda.png'


;lg_nightsky_v2, 1, snapnum= 29, filename='nightsky_today_32.eps', imgname='nightsky_today_32.png', $
;				o_filename='nightsky_today_32_o.eps', o_imgname='nightsky_today_32_o.png', $
;				h_filename='nightsky_today_32_h.eps', h_imgname='nightsky_today_32_h.png'


lg_nightsky_v2, 1, snapnum= 29, basename='nightsky_today_256'
lg_nightsky_v2, 1, snapnum= 46, basename='nightsky_first_256'
lg_nightsky_v2, 1, snapnum= 49, basename='nightsky_post1_256'
lg_nightsky_v2, 1, snapnum= 57, basename='nightsky_second_256'
lg_nightsky_v2, 1, snapnum= 59, basename='nightsky_post2_256'
lg_nightsky_v2, 1, snapnum= 62, basename='nightsky_final_256'
lg_nightsky_v2, 1, snapnum= 80, basename='nightsky_milkomeda_256'


end



;
;  v2
;--------------------------
;pro lg_nightsky_v2, junk, snapnum=snapnum, $
;			filename=filename, imgname=imgname, $
;			o_filename=o_filename, o_imgname=o_imgname, $
;			h_filename=h_filename, h_imgname=h_imgname
pro lg_nightsky_v2, junk, snapnum=snapnum, basename=basename

filename=basename+'.eps'
imgname=basename+'.png'
o_filename=basename+'_o.eps'
o_imgname=basename+'_o.png'
h_filename=basename+'_h.eps'
h_imgname=basename+'_h.png'


frun= "/raid4/tcox/localgroup/v7"
if not keyword_set(snapnum) then snapnum= 29
if not keyword_set(filename) then filename= 'nightsky_nh4.eps'
if not keyword_set(imgname) then imgname= 'nightsky_nh4.png'

if not keyword_set(junk) then begin
   print, "  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    "
   print, "  "
   print, "  "
   print, "     Need to use hidl for this procedure"
   print, "  "
   print, "  "
   print, "  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    "
   print, "  "
   print, "  lg_nightsky_v2, junk"
   print, "  "
   print, "  static variables:"
   print, "       frun= ", frun
   print, "    snapnum= ", snapnum
   print, "   filename= ", filename
   print, "  "
   print, "  "
   print, "  "
   return
endif




; ------------------------------------------------------------------------



if (fload_snapshot_bh(frun, snapnum)) then begin
   print, "PROBLEM: opening ",frun,snapnum
   return
endif


;
; standard /raid4/tcox/localgroup/v* (2-19)
;
g1_sid= 1L
g1_npart= 500001L
mw_bhid= 500001L
g2_sid= 500002L
g2_npart= 800001L
m31_bhid= 1300002L


;frun= "/raid4/tcox/localgroup/bhires"
;mw_bhid= 817955
;m31_bhid= 2145299


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
;star_to_follow= 484414L      ; -110, -2
star_to_follow= 484400L      ; ?
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
	Nstars= fload_npart(2)+fload_npart(3)+fload_npart(4)

	; new HealPix stuff (read it in)
	read_healpix_angles, Nside, pixeln, theta, phi
	Ntheta= Nside
	Nphi= 12*Nside

	earth_x= 0.0
	earth_y= 0.0
	earth_z= 0.0

        ; stars
        x = fload_allstars_xyz('x', center=center2)
        y = fload_allstars_xyz('y', center=center2)
        z = fload_allstars_xyz('z', center=center2)
        vx = fload_allstars_v('x', comvel=center2velocity)
        vy = fload_allstars_v('y', comvel=center2velocity)
        vz = fload_allstars_v('z', comvel=center2velocity)

	x_new= x
        y_new= y
        z_new= z
	vx_new= vx
        vy_new= vy
        vz_new= vz

	mass= fload_allstars_mass(1)
	print, "mass min/max= ", min(mass), max(mass)


	; ------------
	; rotate so that MW is at 0,0
	;
	;rotate_MW_to00 = 1
	rotate_MW_to00 = 0
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
	;rotate_non_MW= 1
	rotate_non_MW= 0
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



	;
        ;  Project the Stellar Particles
	; -------------------------------
	los_stars= pixeln
	los_stars(*)= 1.0e-10

	for i= 0L, Nstars-1 do begin
	   this_R= sqrt(x_new[i]*x_new[i] + y_new[i]*y_new[i] + z_new[i]*z_new[i])
	   if this_R gt 1.0e-3 then begin
		this_theta= acos(z_new[i]/this_R)
		this_phi= atan(y_new[i], x_new[i]) + !PI

		; find nearest theta/phi value
		theta_idx= find_closest(this_theta, theta)
		phi_idx= find_closest(this_phi, phi)

		; now find the ONE
		if n_elements(theta_idx) ge n_elements(phi_idx) then begin
			diff= abs(phi(theta_idx) - this_phi)
			final_idx= theta_idx(where(diff eq min(diff)))
		endif else begin
			diff= abs(theta(phi_idx) - this_theta)
			final_idx= phi_idx(where(diff eq min(diff)))
		endelse

		if (i mod 10000) eq 0 then print, "i= ", i

		; put mass somewhere
		los_stars(final_idx) = los_stars(final_idx) + mass(i)
	   endif
	endfor
	print, "i= ", i

	print, "max/min Stellar Density= ",max(los_stars), min(los_stars)

	print, "total mass= ", total(mass)
	print, "area of each element ()= ", 4.*!PI/n_elements(los_stars)
	print, "total(los_stars)= ", total(los_stars)
	print, "total area (based upon los_stars)= ", total(los_stars * 4.*!PI/n_elements(los_stars))
	print, total(los_stars) * 4. * !PI


	;
	;  N_e
	; -------
	; trap for really low stellar density values
	;idx= where(los_stars lt 1.0e+10)
	;if idx(0) ne -1 then los_stars(idx)= 1.0e+10

	QuantityToProject= los_stars

	;process_data_to_epsfile, QuantityToProject, MW_phi, MW_theta, M31_phi, M31_theta, $
	;			Ntheta, Nphi, $
	;			cbtitlename= '!6Log N!DH!N (cm!E-2!N)', filename='columndensity.eps'


	lg_mollview, QuantityToProject, ps=filename, /graticule, /online, $
					title='N!Dstars!N (10!E10!N M!D!9n!6!N)', $
					png=imgname, /crop, /log, colt=0, pxsize=2400


	lg_orthview, QuantityToProject, ps=h_filename, /graticule, /online, $
					title='N!Dstars!N (10!E10!N M!D!9n!6!N)', $
					png=h_imgname, /crop, /log, colt=0, pxsize=2400, $
					/half_sky


	lg_orthview, QuantityToProject, ps=o_filename, /graticule, /online, $
					title='N!Dstars!N (10!E10!N M!D!9n!6!N)', $
					png=o_imgname, /crop, /log, colt=0, pxsize=2400


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



; /raid4/tcox/localgroup/v3
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








; ----------------------------------------
; read file
pro read_healpix_angles, nside, pixeln, theta, phi

healpix= '/home2/tcox/CMBStuff/HealPix_RingPixels/healpixangles.txt'
print, "reading file= ", healpix
openr, 1, healpix
junk=''
readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & Nside= long(tempjunk(5))
print, "Nside= ", Nside
readf, 1, junk
readf, 1, junk
readf, 1, junk

npixels= 12.*nside*nside
print, "npixels= ", npixels
pixeln= findgen(npixels)
theta= findgen(npixels)
phi= findgen(npixels)

for i=0L,npixels-1 do begin
        ;readf, 1, FORMAT= '(F6.1," ",F8.5,"  ",F8.5)', $
        ;        pixeln[i], theta[i], phi[i]
	readf, 1, junk
	tempjunk=strsplit(junk,/extract,count=count)
	pixeln[i]= i
	theta[i]= float(tempjunk(1))
	phi[i]= float(tempjunk(2))
endfor
close, 1    


            
end









function find_closest, this_theta, theta

	temptheta= theta

	diff= abs(theta - this_theta)

	idx= where(diff eq min(diff), nidx)

	return, idx

end





