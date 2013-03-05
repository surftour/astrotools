pro process_rotation, x_orig, y_orig, z_orig, rotate_theta, rotate_phi, x_new, y_new, z_new

	if not keyword_set(rotate_theta) then rotate_theta= 0.0
	if not keyword_set(rotate_phi) then rotate_phi= 0.0

        print, "*** process_rotation ***"
        print, "rotate: by (theta,phi)=", rotate_theta,rotate_phi

        ; transform to radians
	theta= 0.0
	phi= 0.0
        if strpos('0123456789',strmid(strcompress(string(rotate_theta),/remove_all),0,1)) gt 0 then theta= !PI*float(rotate_theta)/180.0
        if strpos('0123456789',strmid(strcompress(string(rotate_phi),/remove_all),0,1)) gt 0 then phi= !PI*float(rotate_phi)/180.0

	x= x_orig
	y= y_orig
	z= z_orig

        ; rotate the coordinates
        ; ------------------------



	if strmid(string(rotate_phi),0,2) eq "xy" or strmid(string(rotate_theta),0,2) eq "xy" then begin
		x_new= x
		y_new= y
		z_new= z
		return
	endif

	if strmid(string(rotate_phi),0,2) eq "xz" or strmid(string(rotate_theta),0,2) eq "xz" then begin
		x_new= x
		y_new= z
		z_new= y
		return
	endif

	if strmid(string(rotate_phi),0,2) eq "yz" or strmid(string(rotate_theta),0,2) eq "yz" then begin
		x_new= y
		y_new= z
		z_new= x
		return
	endif



        ; 1.
        ; this transformation is 1st rotating around the y-axis such
        ; that +z -> +x (theta), and then 2nd rotating around the z-axis such
        ; that +y -> +x (phi).
        ;
	;
        ;x_new=  x*cos(theta)*cos(phi)  - y*sin(phi) + z*cos(phi)*sin(theta)
        ;y_new=  x*cos(theta)*sin(phi)  + y*cos(phi) + z*sin(phi)*sin(theta)
        ;z_new= -x*sin(theta)                        + z*cos(theta)
    


        ; 2. alternatively, this is the inverse (or not exactly) rotation of the above, 
        ;    where we rotate around z-axis (by phi) and then y-axis (by theta)
	;
	;
        ;x_new=  x*cos(theta)*cos(phi) + y*cos(theta)*sin(phi) + z*sin(theta)
        ;y_new= -x*sin(theta)          + y*cos(theta)
        ;z_new=  x*sin(theta)*cos(phi) - y*sin(theta)*sin(phi) + z*cos(theta)



    
        ; 2.5 first around z axis (phi) - clockwise
        ; then around y axis (theta)  - counter clockwise
	;
	;
        x_new=  x*cos(theta)*cos(phi)  + y*(cos(theta)*sin(phi)) - z*sin(theta)
        y_new= -x*sin(phi)             + y*cos(phi)
        z_new=  x*sin(theta)*cos(phi)  + y*(sin(theta)*sin(phi)) + z*cos(theta)



    
        ; 3. first around x axis (theta)
        ; then around z axis (phi)
	;
	;
        ;x_new=  x*cos(phi)   + y*(sin(phi)*cos(theta)) + z*(sin(phi)*sin(theta))
        ;y_new= -x*sin(phi)   + y*(cos(phi)*cos(theta)) + z*(cos(phi)*sin(theta))
        ;z_new=               - y*sin(theta)            + z*cos(theta)




        ; 4. first around y axis (theta)
        ; then around x axis (phi)
	;
	;
        ;x_new=  x*cos(theta)                           + z*sin(theta)
        ;y_new= -x*(sin(phi)*sin(theta))  + y*cos(phi)  + z*(cos(theta)*sin(phi))
        ;z_new= -x*(cos(phi)*sin(theta))  - y*sin(phi)  + z*(cos(theta)*cos(phi))




        ; 5. first around z axis (theta) in counter-clockwise sense
        ; then around x axis (phi) in clockwise sense
	;
	; this effectively rotates to the standard spherical coordinates theta, phi
	;
	; i'll caution, however, that i've had confusions in the past about the
	; order of my matrix operators, so double check this!!!
	;
	;
	;  R_x(phi) * R_z2(theta)
	;
        ;x_new=  x*cos(theta)          + y*sin(theta)
        ;y_new= -x*cos(phi)*sin(theta) + y*cos(phi)*cos(theta) + z*sin(phi)
        ;z_new=  x*sin(phi)*sin(theta) - y*sin(phi)*cos(theta) + z*cos(phi)
	;
	;  R_x(-phi) * R_z2(theta)
	;
        ;x_new=  x*cos(theta)          - y*sin(theta)
        ;y_new=  x*cos(phi)*sin(theta) + y*cos(phi)*cos(theta) - z*sin(phi)
        ;z_new=  x*sin(phi)*sin(theta) + y*sin(phi)*cos(theta) + z*cos(phi)
	;
	;
	;
	;
	; inverse of above?
	;
        ;x_new=  x*cos(theta) + y*cos(phi)*sin(theta) + z*sin(phi)*sin(theta)
        ;y_new= -x*sin(theta) + y*cos(phi)*cos(theta) + z*sin(phi)*cos(theta)
        ;z_new=               - y*sin(phi)            + z*cos(phi)



end

