
; overplot the image
; --------------------
pro process_kin_image_and_overplot, x, y, z, mass, xlen, $
			x0, y0, x1, y1, $
			hsml=hsml, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			center=center, showbhs=showbhs, $
			show_re=show_re, R_e=R_e, $
			show_xlen=show_xlen, $
			fitstoo=fitstoo, pixels=pixels, $
			msg1=msg1, msg2=msg2, filename=filename, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			ReturnedMap=ReturnedMap, $
			msg=msg, show_maj_slit=show_maj_slit, show_min_slit=show_min_slit, $
			xslit=xslit, yslit=yslit, $
			fit_a=fit_a, fit_b=fit_b, fit_ellip=fit_ellip, a4diva=a4diva, $
			pa=pa, HalfMassSB=HalfMassSB, colortable=colortable, $
			FindMassFactor=FindMassFactor

	m= mass
	if keyword_set(center) then orig_center= center else orig_center=[0,0,0]
	if not keyword_set(colortable) then colortable= 0


        print, "====================================="
        print, "    process_kin_image_and_overplot   "
        print, "====================================="



	x_orig= x
	y_orig= y
	z_orig= z
	center= orig_center
	rotate_phi_orig= rotate_phi
	rotate_theta_orig= rotate_theta
        contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
			hsml=hsml, $
                        filename=filename, fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, $
                        pixels=pixels, zthickness=zthickness, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        crude=crude, center=center, /compute_halfmasscnt, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, domaxrng=domaxrng, $
                        NxNImage=NxNImage, HalfMassSB=HalfMassSB, FindMassFactor=FindMassFactor
        
        Pic= NxNImage
	ReturnedMap= NxNImage

	pixels= (size(NxNImage))[1]
	print, "pixels= ", pixels

    ;print, "center= ", center
    ;print, "center (original)= ", orig_center
    center= [0,0,0]

    xmin= -xlen+center[0]
    xmax=  xlen+center[0]
    ymin= -xlen+center[1]
    ymax=  xlen+center[1]


    fload_newcolortable, colortable


    tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


    !p.position=[x0,y0,x1,y1]
    !p.ticklen=0.03

                ; creates axes and plot style
                ; ---------------------------
                plot,[0],[0], psym=3,  $
                      xrange=[xmin,xmax], $
                      yrange=[ymin,ymax], $
                      /NOERASE, $
                      color= 0, $
                      xcharsize=0.01, $
		      ycharsize=0.01, $
                      xthick=4.0, $
		      ythick=4.0, $
                      charthick=3.0, $
		      ;xticks=10, $
		      ;yticks=10, $
		      xtickformat='(a1)', ytickformat='(a1)', $
                      xstyle=1, ystyle=1, $
                      ;xtitle=xtit, $
                      ;ytitle=ytit, $ ;/normal, $
                      /nodata


     ;  key
     ; ------

	produce_key= 0
	;produce_key= 1
	if produce_key eq 1 then begin
	     ; white box
	     bx0=xmin+1.0
	     bx1=xmin+6.0
	     by0=ymax-1.0
	     by1=ymax-3.6
	     polyfill, [bx0,bx0,bx1,bx1], [by0,by1,by1,by0], color= 1
	     oplot, [bx0,bx0,bx1,bx1,bx0], [by0,by1,by1,by0,by0] , color= 0, thick=3.0

	     ;xyouts, x0+0.04, y1-0.08, fload_timelbl(1,2,/noteq), /normal, size= 1.7, charthick=3.0, color= 0
	     xyouts, x0+0.02, y1-0.05, msg1, /normal, size= 1.3, charthick=4.0, color= 0
	     xyouts, x0+0.02, y1-0.10, msg2, /normal, size= 1.3, charthick=4.0, color= 0
	endif



                ; -----------------------------
                ;  X marks the black hole spot
                ; -----------------------------
                if keyword_set(showbhs) and fload_npart(5) ge 1 then begin

			fload_newcolortable, 4

                        bhpsym= 7
                        bhsymsize= 2.0
                        bhcolor= 220

                        n_bh= fload_npart(5)
                        bhid= fload_blackhole_id(1)
                        bhid1= bhid[0]
                        bh_center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
                        
                        bh_new= bh_center_1 - orig_center
                        if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin
                           bh_c= bh_new
                           theta= !PI*rotate_theta/180.0
                           phi= !PI*rotate_phi/180.0
                           ; 4. rotate y (theta) and then x (phi)
                           ;bh_new[0]= bh_c[0]*cos(theta) + bh_c[2]*sin(theta)
                           ;bh_new[1]= -bh_c[0]*(sin(phi)*sin(theta)) + bh_c[1]*cos(phi) + bh_c[2]*(cos(theta)*sin(phi))
                           ;bh_new[2]= -bh_c[0]*(cos(phi)*sin(theta)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(theta)*cos(phi))
                           ; 3. rotate x (theta) and then z (phi)
                           ;bh_new[0]= bh_c[0]*cos(phi) + bh_c[1]*(sin(phi)*cos(theta)) + bh_c[2]*(sin(phi)*sin(theta))
                           ;bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*(cos(phi)*cos(theta)) + bh_c[2]*(cos(phi)*sin(theta))
                           ;bh_new[2]= -bh_c[1]*sin(theta) + bh_c[2]*cos(theta)
			   ; 1. rotate y (theta) and then z (phi)
			   ;bh_new[0]= bh_c[0]*(cos(theta)*cos(phi)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(phi)*sin(theta))
			   ;bh_new[1]= bh_c[0]*(cos(theta)*sin(phi)) + bh_c[1]*cos(phi) + bh_c[2]*(sin(phi)*sin(theta))
			   ;bh_new[2]= -bh_c[0]*sin(theta) + bh_c[2]*cos(theta)
	        	   ; 2.5 first around z axis (phi) - clockwise
	        	   ; then around y axis (theta)  - counter clockwise
	        	   bh_new[0]= bh_c[0]*(cos(theta)*cos(phi)) + bh_c[1]*(cos(theta)*sin(phi)) - bh_c[2]*sin(theta)
	        	   bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*cos(phi)
	        	   bh_new[2]= bh_c[0]*(sin(theta)*cos(phi)) + bh_c[1]*(sin(theta)*sin(phi)) + bh_c[2]*cos(theta)
                        endif
                        
                        ; for the moment, just plot xy black hole position
                        oplot, [bh_new[0]], [bh_new[1]], psym=bhpsym, thick=6.0, $
                                        color=bhcolor, symsize=bhsymsize
                        
                        if n_bh gt 1 then begin
                                bhid2= bhid[1]
                                bh_center_2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
                                bh_new= bh_center_2 - orig_center
                                if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin
                                   bh_c= bh_new
                                   theta= !PI*rotate_theta/180.0
                                   phi= !PI*rotate_phi/180.0
                                   ; 4. rotate y (theta) and then x (phi)
                                   ;bh_new[0]= bh_c[0]*cos(theta) + bh_c[2]*sin(theta)
                                   ;bh_new[1]= -bh_c[0]*(sin(phi)*sin(theta)) + bh_c[1]*cos(phi) + bh_c[2]*(cos(theta)*sin(phi))
                                   ;bh_new[2]= -bh_c[0]*(cos(phi)*sin(theta)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(theta)*cos(phi))
                                   ; 3. rotate x (theta) and then z (phi)
                                   ;bh_new[0]= bh_c[0]*cos(phi) + bh_c[1]*(sin(phi)*cos(theta)) + bh_c[2]*(sin(phi)*sin(theta))
                                   ;bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*(cos(phi)*cos(theta)) + bh_c[2]*(cos(phi)*sin(theta))
                                   ;bh_new[2]= -bh_c[1]*sin(theta) + bh_c[2]*cos(theta)
				   ; 1. rotate y (theta) and then z (phi)
				   ;bh_new[0]= bh_c[0]*(cos(theta)*cos(phi)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(phi)*sin(theta))
				   ;bh_new[1]= bh_c[0]*(cos(theta)*sin(phi)) + bh_c[1]*cos(phi) + bh_c[2]*(sin(phi)*sin(theta))
				   ;bh_new[2]= -bh_c[0]*sin(theta) + bh_c[2]*cos(theta)
				   ; 2.5 first around z axis (phi) - clockwise
				   ; then around y axis (theta)  - counter clockwise
				   bh_new[0]= bh_c[0]*(cos(theta)*cos(phi)) + bh_c[1]*(cos(theta)*sin(phi)) - bh_c[2]*sin(theta)
				   bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*cos(phi)
				   bh_new[2]= bh_c[0]*(sin(theta)*cos(phi)) + bh_c[1]*(sin(theta)*sin(phi)) + bh_c[2]*cos(theta)
                                endif
                                
                                print, "original center= ",orig_center
                                print, "Blackhole #1: ",bhid1, bh_center_1
                                print, "Blackhole #2: ",bhid2, bh_center_2
                                
                                oplot, [bh_new[0]], [bh_new[1]], psym=bhpsym, thick=6.0, $                                  
                                        color=bhcolor, symsize=bhsymsize
                                
                        endif else begin
                                print, "original center= ",orig_center
                                print, "Blackhole: ", bhid1, bh_center_1
                        endelse
                
                endif

		fload_newcolortable, 4

		if keyword_set(pa) then pa_orig= pa else pa_orig= -1

		; overplot contours, and ellipse
		;oplot_contours, Pic, x0, y0, x1, y1, hmsb=HalfMassSB, xlen=xlen, pa=pa
		oplot_contours, Pic, x0, y0, x1, y1, hmsb=HalfMassSB, $
					;/fitellipse, xlen=xlen, pa=pa, $      ; don't think i need to fitellipse
					xlen=xlen, pa=pa, $
					/fitdiskyboxy, $
					fit_a=fit_a, fit_b=fit_b, $
					fit_ellip=fit_ellip, a4diva=a4diva, pixels= pixels


		if pa_orig ne -1 then pa= -1.0*pa_orig   ; minus sign so manual setting
							 ; goes clockwise


                ; -------------------------------
                ;  Draw Box (slit, for instance)
                ; -------------------------------
                ;show_maj_slit= 1
                if keyword_set(show_maj_slit) and not keyword_set(pa) then begin
                   slitcolor= 1
                   slitthick= 4.0
                   oplot, [-xslit, -xslit], [-yslit,yslit], psym=-3, color= slitcolor, thick= slitthick
                   oplot, [-xslit, xslit], [yslit,yslit], psym=-3, color= slitcolor, thick= slitthick
                   oplot, [xslit, xslit], [-yslit,yslit], psym=-3, color= slitcolor, thick= slitthick
                   oplot, [-xslit, xslit], [-yslit,-yslit], psym=-3, color= slitcolor, thick= slitthick
                endif

		print, "xyslit= ",xslit, yslit


		; rotate just x and y so that it's
		; along the major axis position angle (pa)
		; ----------------------------------------
		if keyword_set(pa) and keyword_set(show_maj_slit) then begin
		   print, "slit pa= ",pa
		   ; version 0
		   ; actually the pa passed is clockwise rotation, so we
		   ; need to flips this around
		   ;if pa gt 0.0 then pa= 180.0 - pa else pa= -pa
		   ; version 1 (i think i'm fixing things here??)
		   ;pa= 90.0 - pa
		   ; version 2 (newer diskboxy pa is counter-clockwise rotation)
		   pa= 90.0 + pa
		   if pa gt 180.0 then pa= pa - 180.0    ; doesn't change slit position, but does change + side
		   par = pa * !PI/180.0
		   print, "new pa= ", pa

		   xtr= xslit*cos(par) - yslit*sin(par)       ; this choice of signs rotates
		   ytr= xslit*sin(par) + yslit*cos(par)       ; counter-clockwise

		   xbr= xslit*cos(par) - (-yslit)*sin(par)
		   ybr= xslit*sin(par) +(-yslit)*cos(par)

		   xtl= (-xslit)*cos(par) - yslit*sin(par)
		   ytl= (-xslit)*sin(par) + yslit*cos(par)

		   xbl= (-xslit)*cos(par) - (-yslit)*sin(par)
		   ybl= (-xslit)*sin(par) + (-yslit)*cos(par)

		   ; need axes to be reset

    		   !p.position=[x0,y0,x1,y1]
    		   !p.ticklen=0.03

                   ; creates axes and plot style
                   ; ---------------------------
                   plot,[0],[0], psym=3,  $
                      xrange=[xmin,xmax], $
                      yrange=[ymin,ymax], $
                      /NOERASE, $
                      color= 0, $
                      xcharsize=0.01, $
		      ycharsize=0.01, $
                      xthick=4.0, $
		      ythick=4.0, $
                      charthick=3.0, $
		      ;xticks=10, $
		      ;yticks=10, $
		      xtickformat='(a1)', ytickformat='(a1)', $
                      xstyle=1, ystyle=1, $
                      ;xtitle=xtit, $
                      ;ytitle=ytit, $ ;/normal, $
                      /nodata


                   ;slitcolor= 1    ; white
                   slitcolor= 240     ; yellow
                   slitthick= 4.0
                   oplot, [xtr,xbr,xbl,xtl,xtr], [ytr,ybr,ybl,ytl,ytr], psym=-3, color= slitcolor, thick= slitthick

		   ;
		   ; show the minor axis slit
                   if keyword_set(show_min_slit) then oplot, -[ytr,ybr,ybl,ytl,ytr], [xtr,xbr,xbl,xtl,xtr], psym=-3, color= slitcolor, thick= slitthick

		endif


	if keyword_set(show_re) then begin
	       ; 2.5 first around z axis (phi) - clockwise
	       ; then around y axis (theta)  - counter clockwise
		;if keyword_set(rotate_theta) then theta= !PI*rotate_theta/180.0 else theta= 0.0
		;if keyword_set(rotate_phi) then phi= !PI*rotate_phi/180.0 else phi= 0.0
		;xin= xin-orig_center[0]
		;yin= yin-orig_center[1]
		;zin= zin-orig_center[2]
		;x_new= xin*(cos(theta)*cos(phi)) + yin*(cos(theta)*sin(phi)) - zin*sin(theta)
		;y_new= -xin*sin(phi) + yin*cos(phi)

		; rotate
		; ---------
		rotate_phi= rotate_phi_orig
		rotate_theta= rotate_theta_orig
		center= orig_center
		if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin

		        print, "rotate: moving to center", center
		        xin=x_orig-center[0]
		        yin=y_orig-center[1]
		        zin=z_orig-center[2]

		        print, "rotate: set center = [0,0,0]"
		        center= [0,0,0]

		        process_rotation, xin, yin, zin, rotate_theta, rotate_phi, x_new, y_new, z_new
		        ;process_rotation, vxin, vyin, vzin, rotate_theta, rotate_phi, vx_new, vy_new, vz_new
		endif else begin
			x_new= x_orig - center[0]
			y_new= y_orig - center[1]
			z_new= z_orig - center[2]
		endelse




		; actually figure out R_e
		;R_e= 2.4384
		;R_e= process_2d_halfmass(x_new,y_new,m)
		R_e= process_2d_halfmass(x_new,y_new,m,FindMassFactor=FindMassFactor)
		;Sig= process_2d_sigma(x_new,y_new,vz_new,R_e*0.5)

		phi = dindgen(101)*2D*!dpi/100
		; data coord
		re_x = R_e * cos(phi)
		re_y = R_e * sin(phi)

		;re_x= re_x + center[0]
		;re_y= re_y + center[1]

		oplot, re_x, re_y, color= 50, thick=4.0, linestyle=0, psym=-3
		relbl= strcompress(string(R_e),/remove_all)
		relbl= '!6R!De!N='+strmid(relbl,0,4)+' kpc/h'
		xyouts, x1-0.13, y0+0.05, relbl, /normal, size=1.4, color= 0, charthick= 3.0

	endif else begin
		R_e= 0.0
		;Sig= 0.0
	endelse



	if keyword_set(show_xlen) then begin
		xlenlbl= strcompress(string(2.0*xlen),/remove_all)
		if (2.0*xlen) ge 1.0 then digs= 1
		if (2.0*xlen) ge 10.0 then digs= 2
		if (2.0*xlen) ge 100.0 then digs= 3
		xlenlbl = strmid(xlenlbl,0,digs)        ; T=0.x (4+digits after decimal)
		xlenlbl= '!94!6'+xlenlbl+' kpc/h !96!6'
		xyouts, x1-0.10, y1-0.07, xlenlbl, /normal, size= 1.2, charthick=3.0, color= 0
	endif


end





; -----------------------------------------------------------------------------







