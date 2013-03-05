

; does the work for each panel
; -------------------------------
pro img_do_one_panel, frun, snapnum, xlen, $
	x0, y0, xs, ys, $
	colorbar=colorbar, draw_circle=draw_circle, $
	origcolortable=origcolortable, $
	center=center, $
	plottype=plottype, $
	fitstoo=fitstoo, arepo=arepo, $
	xthickness=xthickness, ythickness=ythickness, $
	pixels=pixels, zthickness=zthickness, $
	crude=crude, msg=msg, plotpts=plotpts, $
	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	nolabels=nolabels, pubstyle=pubstyle, particlesonly=particlesonly, $
	showbhs=showbhs, set_maxden=set_maxden, set_dynrng=set_dynrng, $
	track_to_draw=track_to_draw, colortable=colortable


	if not keyword_set(origcolortable) then origcolortable= 4

	if keyword_set(center) then begin
        	orig_center= center
	endif else begin
        	orig_center=[0,0,0]
        	center=[0,0,0]
	endelse




	;  get images
	; -------------
	;center=[0,0,0]
	img_produce_one_image, frun, snapnum, xlen, xz=xz, yz=yz, $
		xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
		pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, crude=crude, $
		center=center, $
		ReturnImage=ReturnImage, ReturnTimeLbl= ReturnTimeLbl, $
		arepo=arepo, plottype=plottype

	center_used= center

	Pic= ReturnImage
	TimeLbl= ReturnTimeLbl



	tv, Pic, x0, y0, xsize=xs, ysize=ys, /normal


	x1= x0 + xs
	y1= y0 + ys
	oplot_one_axes, x0, y0, x1, y1, TimeLbl, xlen    ;, /add_j_vector


	;snaplbl= '('+strcompress(string(snapnum),/remove_all)+')'
        ;xyouts, x0+0.03, y0+ys-0.05, snaplbl, /normal, size= 1.0, charthick=3.0, color= 0


	; generate white axes if desired
	;axis, xaxis=0, xstyle=1,ystyle=1,color=1,charsize=0.001, xthick=2.0,ythick=2.0
	;axis, xaxis=1, xstyle=1,ystyle=1,color=1,charsize=0.001, xthick=2.0,ythick=2.0
	;axis, yaxis=0, xstyle=1,ystyle=1,color=1,charsize=0.001, xthick=2.0,ythick=2.0
	;axis, yaxis=1, xstyle=1,ystyle=1,color=1,charsize=0.001, xthick=2.0,ythick=2.0




; ====================
;
;   Extras  (do search for number repeated)
;
;	1. plot points
;	2. draw centers  (several versions, some old)
;	3. show BHs
;	4. overplot IDs
;	5. draw various shapes (circle, rectangle, box)
;	6. add a J vector
;	7. overplot a contour of some sort
;	8. add a colorbar to the plot
;	9. add velocity vector field
;	10. do disky/boxy calc.
;
;
; ====================



		;============================================================================
		;  1111111111111111111


                ; plot the actual points
		; ------------------------
                if keyword_set(plotpts) or keyword_set(particlesonly) then begin
                  if keyword_set(yz) then begin
                        x= y
                        y= z
                  endif
   
                  if keyword_set(xz) then y= z
                  ;oplot, x, y, psym=5, thick=2.0, color= 50
		  ;oplot, x, y, psym=2, thick=2.0, color= 150
		  ;oplot, x, y, psym=7, thick=2.0, color= 100
		  ;oplot, x, y, psym=1, thick=2.0, color= 0
		  oplot, x, y, psym=3, color= 0
                endif



		;============================================================================
		;  22222222222222222222
		;
		;  Draw centers of some sort
		;  
		;




		; OLD OLD OLD OLD OLD
                ; --------------------------
                ;   Draw track
                ; --------------------------
                if keyword_set(track_to_draw) then begin

                        frun= fload_frun(1)
                        read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=frun+"/centers.txt"

                        ; hard-coding xy projection
                        xtrack1= transpose(cp_cen1(0,*))
                        ytrack1= transpose(cp_cen1(1,*))
                        xtrack2= transpose(cp_cen2(0,*))
                        ytrack2= transpose(cp_cen2(1,*))

                        idx= where(cp_time lt fload_time(1))

                        if idx(0) ne -1 then xtrack1= xtrack1(idx)
                        if idx(0) ne -1 then ytrack1= ytrack1(idx)
                        oplot, xtrack1, ytrack1, color= 0, thick=3.0, linestyle= 2, psym=-3

                        if idx(0) ne -1 then xtrack2= xtrack2(idx)
                        if idx(0) ne -1 then ytrack2= ytrack2(idx)
                        oplot, xtrack2, ytrack2, color= 150, thick=3.0, linestyle= 2, psym=-3

                        ;xtrack= track_to_draw[*,Axis1]
                        ;ytrack= track_to_draw[*,Axis2]
                        ;idx= where((xtrack ne 0) or (ytrack ne 0))
                        ;if idx(0) ne -1 then begin
                        ;  xtrack= xtrack(idx)
                        ;  ytrack= ytrack(idx)
                        ;  oplot, xtrack, ytrack, color= 0, thick=3.0, linestyle= 2, psym=-3
                        ;endif
                endif



		; OLD OLD OLD OLD OLD
                ; -----------------------------
                ;  Draw center path, similar to
                ;  tracktodraw
                ; -----------------------------
                ;drawcenters= 1
                drawcenters= 0
                if keyword_set(drawcenters) then begin

                        if total(abs(center_used)) gt 0.0 then center= center_used else center= [0,0,0]
                        print, "center_used= ", center_used
                        print, "center= ", center

                        ; warning, need to load process_directory for this
                        read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=fload_frun(1)+"/centers.txt"

                        cp_cen1_x= cp_cen1[0,*] - center[0]
                        cp_cen1_y= cp_cen1[1,*] - center[1]
                        cp_cen1_z= cp_cen1[2,*] - center[2]
                        idx= where(cp_time le fload_time(1))
                        cp_cen1_x= cp_cen1_x(idx)
                        cp_cen1_y= cp_cen1_y(idx)
                        oplot, cp_cen1_x, cp_cen1_y, psym=-3, linestyle=1, color= 50, thick=2.0

                        cp_cen2_x= cp_cen2[0,*] - center[0]
                        cp_cen2_y= cp_cen2[1,*] - center[1]
                        cp_cen2_z= cp_cen2[2,*] - center[2]
                        idx= where(cp_time le fload_time(1))
                        cp_cen2_x= cp_cen2_x(idx)
                        cp_cen2_y= cp_cen2_y(idx)
                        ;oplot, cp_cen2_x, cp_cen2_y, psym=-3, linestyle=1, color= 150, thick=4.0


                        ; double up baby
                        ;read_centerpositions, cp_time, cp_cen1, cp_cen2, filename="Data/nodrag/centers.txt"
                        ;cp_cen2_x= cp_cen2[0,*]
                        ;cp_cen2_y= cp_cen2[1,*]
                        ;idx= where(cp_time le fload_time(1))
                        ;cp_cen2_x= cp_cen2_x(idx)
                        ;cp_cen2_y= cp_cen2_y(idx)
                        ;oplot, cp_cen2_x, cp_cen2_y, psym=-3, linestyle=0, color= 100, thick=4.0


                endif




                ; -----------------------------
                ;
                ; Print the center separation
                ;
                ; -----------------------------
                if keyword_set(show_rsep) and fload_npart(5) eq 2 then begin

                        bhid= fload_blackhole_id(1)

                        ; bh 1
                        bhidlbl1= strcompress(string(bhid[0]),/remove_all)
                        read_center, ttime1, center1, filename=fload_frun(1)+"/centers_bh_"+bhidlbl1+".txt"
                        idx= where(ttime1 ge fload_time(1)-0.01)
                        bhcen1= center1(*,idx(0))

                        ; bh 2
                        bhidlbl2= strcompress(string(bhid[1]),/remove_all)
                        read_center, ttime2, center2, filename=fload_frun(1)+"/centers_bh_"+bhidlbl2+".txt"
                        idx= where(ttime2 ge fload_time(1)-0.01)
                        bhcen2= center2(*,idx(0))


                        rdiff= bhcen1-bhcen2
                        rsep= sqrt(rdiff(0)*rdiff(0) + rdiff(1)*rdiff(1) + rdiff(2)*rdiff(2))
                        rprojsep= sqrt(rdiff(0)*rdiff(0) + rdiff(1)*rdiff(1))

                        rsep=rsep/0.7
                        rprojsep=rprojsep/0.7
                        print, "separations: ", rsep, rprojsep

                        rseplbl= strcompress(string(rsep),/remove_all)
                        if (rsep) ge 1.0 then digs= 1
                        if (rsep) ge 10.0 then digs= 2
                        rseplbl = strmid(rseplbl,0,digs+2)        ; T=0.x (4+digits after decimal)
                        rseplbl= '!7D!6R= '+rseplbl+' kpc'
                        xyouts, 0.07, 0.90, rseplbl, /normal, size= 1.4, charthick=3.0, color= 0

                        rsepprojlbl= strcompress(string(rprojsep),/remove_all)
                        if (rsep) ge 1.0 then digs= 1
                        if (rsep) ge 10.0 then digs= 2
                        rsepprojlbl = strmid(rsepprojlbl,0,digs+2)        ; T=0.x (4+digits after decimal)
                        rsepprojlbl= '!7D!6R!Dproj!N= '+rsepprojlbl+' kpc'
                        xyouts, 0.07, 0.85, rsepprojlbl, /normal, size= 1.4, charthick=3.0, color= 0



                endif



		;============================================================================
		;  33333333333333333
		;
		;
		;  show BH's
		;
		;



		; -----------------------------
		;  X marks the black hole spot
		; -----------------------------
		;showbhs= 0
		;showbhs= 1
		if keyword_set(showbhs) and fload_npart(5) ge 1 then begin

			fload_newcolortable, 4

			; symbol
                        usersym,cos(findgen(49)/49*2*!pi),sin(findgen(49)/49*2*!pi), thick=4.0, /fill
                        bhpsym= 8
                        ;bhpsym= 7
                        bhsymsize= 2.0
                        ;bhcolor= 220
                        bhcolor= 150

                        n_bh= fload_npart(5)
                        bhid= fload_blackhole_id(1)


                        print, "original center= ",orig_center

                        ; loop through BH's and mark their location
                        ; ----------------------------------------
                        for i=0, n_bh-1 do begin

                                bhid1= bhid[i]
                                ;bh_center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
                                bh_center_1= fload_blackhole_xyz('xyz',center=orig_center,idtofollow=bhid1)
                                bh_new= bh_center_1

                                if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin

                                   ; new rotation procedure
                                   process_rotation, bh_c[0], bh_c[1], bh_c[2], rotate_theta, rotate_phi, bh_new[0], bh_new[1], bh_new[2]

                                   ;bh_c= bh_new
                                   ;theta= !PI*rotate_theta/180.0
                                   ;phi= !PI*rotate_phi/180.0
                                   ; rotate y (theta) and then x (phi)
                                   ;bh_new[0]= bh_c[0]*cos(theta) + bh_c[2]*sin(theta)
                                   ;bh_new[1]= -bh_c[0]*(sin(phi)*sin(theta)) + bh_c[1]*cos(phi) + bh_c[2]*(cos(theta)*sin(phi))
                                   ;bh_new[2]= -bh_c[0]*(cos(phi)*sin(theta)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(theta)*cos(phi))
                                   ; rotate x (theta) and then z (phi)
                                   ;bh_new[0]= bh_c[0]*cos(phi) + bh_c[1]*(sin(phi)*cos(theta)) + bh_c[2]*(sin(phi)*sin(theta))
                                   ;bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*(cos(phi)*cos(theta)) + bh_c[2]*(cos(phi)*sin(theta))
                                   ;bh_new[2]= -bh_c[1]*sin(theta) + bh_c[2]*cos(theta)
                                endif

                                if keyword_set(xz) then bh_new[1]= bh_new[2]
                                if keyword_set(yz) then begin
                                        bh_new[0]= bh_new[1]
                                        bh_new[1]= bh_new[2]
                                endif

                                ; for the moment, just plot xy black hole position
                                oplot, [bh_new[0]], [bh_new[1]], psym=bhpsym, thick=6.0, $
                                        color=bhcolor, symsize=bhsymsize

                                print, "Blackhole #"+strcompress(string(i),/remove_all)+": ", bhid1, bh_center_1

                        endfor

			fload_newcolortable, origcolortable

		endif



		;============================================================================
		;   4444444444444444444
		;
		;
		;  Various ID overplotting
		;
		;

                ;
                ;  overplot from a list 
                ;
                ; -----------------------
                overplot_idlist= 0
                ;overplot_idlist= 1
                if overplot_idlist eq 1 then begin
                        ;contour_add_idlist, 'dog', center=orig_center, $
                        ;        rotate_theta=rotate_theta, rotate_phi=rotate_phi
                        frun= fload_frun(1)
                        center= orig_center
                        ;contour_add_idlist, frun+'/tidal_tail_idlist.txt', center=center, $
                        ;        rotate_theta=rotate_theta, rotate_phi=rotate_phi
                        contour_add_idlist, frun+'/idlist_rem_rgt4.txt', center=center, $
                        ;contour_add_idlist, frun+'/tidal_tail_idlist.txt', center=center, $
                        ;contour_add_idlist, frun+'/idlist_020_egt0.txt', center=center, $
                        ;contour_add_idlist, frun+'/idlist_020_e1.txt', center=center, $
                        ;contour_add_idlist, frun+'/idlist_020_e2.txt', center=center, $
                        ;contour_add_idlist, frun+'/idlist_020_e3.txt', center=center, $
                                rotate_theta=rotate_theta, rotate_phi=rotate_phi

                        fload_newcolortable, colortable
                endif



		;
                ; overplot the new stars, if we want
		;
                ; ------------------------------------
                ;do_new_star_overplotting= 1
                do_new_star_overplotting= 0
                if do_new_star_overplotting eq 1 then begin
                        print, "getting to overplot"
                        ;overplot_all_newstars, frun, x0, xs, y0, ys, xmin, xmax, ymin, ymax, xlen, $
                        ;               snap1, snap2, snap3, snap4, snap5, snap6, xz=xz, $
                        oplot_newstars, frun, snapnum, x0, y0, x1, y1, xmin, xmax, ymin, ymax, xlen
                endif



                ;
                ;   Add Young Stars
                ;
                ; -----------------------
                overplot_youngstars= 0
                ;overplot_youngstars= 1
                if overplot_youngstars eq 1 then begin
                        fload_newcolortable, 1
                        oplot_contour_add_youngstars, 'junk', center=orig_center, $
                                rotate_theta=rotate_theta, rotate_phi=rotate_phi
                endif





		;============================================================================
		;    55555555555555
		;
		;
		;    Add an angular momentum vector
		;
		;
        	;
        	; ----------------------------
        	if keyword_set(add_j_vector) then begin

                	j_cold= fload_gas_cold_j(1)
                	;j_cold= j_cold-orig_center

                	; rescale it
                	j_cold= 20.0*j_cold/sqrt(total(j_cold*j_cold))

                	; rotate?
                	if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin
                        	; transform to radians
                        	theta= !PI*rotate_theta/180.0
                        	phi= !PI*rotate_phi/180.0

                        	x= j_cold[0]
                        	y= j_cold[1]
                        	z= j_cold[2]
                        	; first around y axis (theta)
                        	; then around x axis (phi)
                        	x_new= x*cos(theta) + z*sin(theta)
                        	y_new= -x*(sin(phi)*sin(theta)) + y*cos(phi) + z*(cos(theta)*sin(phi))
                        	z_new= -x*(cos(phi)*sin(theta)) - y*sin(phi) + z*(cos(theta)*cos(phi))
                        	j_cold[0]= x_new
                        	j_cold[1]= y_new
                        	j_cold[2]= z_new
                	endif

                	;arrow, 0.0, 0.0, j_cold[0], j_cold[1], /data, color= 0, thick=4.0
        	        arrow, 0.0, 0.0, j_cold[0], j_cold[2], /data, color= 0, thick=4.0
	        endif




		;============================================================================
		;    6666666666666
		;
		;
		;    Draw various shapes
		;
		;



                ; ---------------------
                ;  Draw Circle
                ; ---------------------
                ;draw_circle=1
                ;draw_circle=0
                if keyword_set(draw_circle) then begin
                        R_e= 3.0
   
                        phi = dindgen(101)*2D*!dpi/100
                        ; data coord
                        re_x = R_e * cos(phi)
                        re_y = R_e * sin(phi)

                        loadct, 4
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0  ; define colors 0 and 1 as black and white

                        oplot, re_x, re_y, color= 220, thick=4.0, linestyle=0, psym=-3


			loadct, 1          ; blue/white
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0
                        
                endif



                ; ---------------------
                ;  Draw Rectangle
                ; ---------------------
                ;draw_square=1
                ;draw_square=0
                if keyword_set(draw_square) then begin

                        bs= 6.0 * 0.7
   
                        loadct, 4
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0  ; define colors 0 and 1 as black and white

                        oplot, [-bs,bs], [bs,bs], color= 210, thick=4.0, linestyle=1, psym=-3
                        oplot, [-bs,bs], [-bs,-bs], color= 210, thick=4.0, linestyle=1, psym=-3
                        oplot, [-bs,-bs], [-bs,bs], color= 210, thick=4.0, linestyle=1, psym=-3
                        oplot, [bs,bs], [-bs,bs], color= 210, thick=4.0, linestyle=1, psym=-3


			loadct, 0          ; black/white
			;loadct, 1          ; blue/white
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0
                        
                endif



                ; -------------------------------
                ;  Draw Box (slit, for instance)
                ; -------------------------------
                ;drawbox= 1
                if keyword_set(drawbox) then begin
                   slitcolor= 1
                   slitthick= 4.0
                   oplot, [-xslit, -xslit], [-yslit,yslit], psym=-3, color= slitcolor, thick= slitthick
                   oplot, [-xslit, xslit], [yslit,yslit], psym=-3, color= slitcolor, thick= slitthick
                   oplot, [xslit, xslit], [-yslit,yslit], psym=-3, color= slitcolor, thick= slitthick
                   oplot, [-xslit, xslit], [-yslit,-yslit], psym=-3, color= slitcolor, thick= slitthick
                endif




		;============================================================================
		;  777777777777
		;
		;
		;  Overplot a contour of some kind
		;
		;
		;
		;



	        ; ---------------------------
	        ;  Contour Plot
	        ;
	        ; Do a contour plot of 
	        ; this mess.
	        ; ---------------------------
	        ;docontour= 1
	        docontour= 0
	        if docontour EQ 1 then begin

	                oplot_contours, Pic, x0, y0, x1, y1, $
        	                hmsb=HalfMassSB, xlen=xlen , /fitellipse, pa=pa

  	      endif



	        ; ---------------------------
	        ;  Contour Plot (of gas)
	        ;
	        ; ---------------------------
	        docontour= 1
	        ;docontour= 0
	        if docontour EQ 1 then begin
	
	                center= orig_center

        	        oplot_contours_comp, 'gas', x0, y0, x1, y1, $
	                        center=center, $
                	        xlen=xlen, pixels=pixels, $
				rotate_theta=rotate_theta, rotate_phi=rotate_phi

        	endif




        	; --------------------------------
        	;  Contour Plot (of dark matter)
        	;
        	; --------------------------------
        	;docontour= 1
        	docontour= 0
        	if docontour EQ 1 then begin
        	        center= orig_center

        	        ;oplot_contours_comp, 'dm', x0, y0, x1, y1, $
        	        oplot_contours_comp, 'special', x0, y0, x1, y1, $
        	                center=center, $
        	                xlen=xlen, pixels=pixels
        	endif




        	; --------------------------------
        	;  Add X-ray contours
		;   (the following script doesn't exist yet, but
		;    but would be easy to adapt from the above
		;    and old_crap/contour_xrays)
        	;
        	; --------------------------------
        	;docontour= 1
        	docontour= 0
        	if docontour EQ 1 then begin
        	        center= orig_center

        	        ;oplot_contours_comp, 'dm', x0, y0, x1, y1, $
        	        oplot_contours_xrays, 'special', x0, y0, x1, y1, $
        	                center=center, $
        	                xlen=xlen, pixels=pixels
        	endif






		;============================================================================
		;    888888888888888
		;
		;
		;   Add a color bar somewhere
		;
		;


                ; -----------------------
                ;  ColorBar #1
                ; 
                ; The original colorbar.
                ; Oriented on the right
                ; side of the plot, keyword
                ; also affects the plot
                ; positioning above.
                ; ----------------
                if keyword_set(colorbar) then begin

                invertcolors= 1    ; warning: this should match the value in contour_makepic
                MaxDensXY= 5.0e-2
                DynRange=1.0e3
                ma=MaxDensXY
                mi=MaxDensXY/DynRange

                  bar= REPLICATE(1B, 10)#BINDGEN(256)
                  if invertcolors eq 1 then bar= 255-bar
                  idx= where(bar eq 0 or bar eq 1)
                  if idx(0) ne -1 then begin
                        ;print, "idx has ", n_elements(idx), idx(0)
                        ;colornumber= idx(0)-1
                        ;if colornumber lt 0 then colornumber= 11     ; if it's bar[*,0] then set to bar[*,1]
                        bar(idx)= 2   ; bar(colornumber)
                  endif
                  barwidth= 0.04

                  tv, bar, x1+0.01, y0, xsize=barwidth, ysize=(y1-y0), /normal

                  !p.position= [x1+0.01,y0,x1+0.01+barwidth,y1]

                  plot, [0],[0], psym=3, xrange=[0,1], yrange=[alog10(mi),alog10(ma)], $
                        /normal, $
                        /noerase, color= 1, xstyle=1, ystyle=1, $
                        xthick=4.0, ythick=4.0, $
                        xticks=1, xtickformat='(a1)', ytickformat='(a1)', $
                        /nodata

                  axis, yaxis=1, yrange=[alog10(mi)+4.0,alog10(ma)+4.0], ystyle=1, $
                        ;xcharsize=1.50, ycharsize=1.50, charthick=3.0, $
                        /normal, $
                        ytitle= 'Log !4R!D!3gas!N(M!D!9n!3!N pc!E-2!N)'   ; +4.0 above comes from converting to these units
                endif




                ; ----------------------
                ;  ColorBar #2
                ;
                ; This one is on the bottom,
                ; horizontally, and inside
                ; the plotted area (for
                ; nolabels keyword)
                ; ----------------
		colorbar2= 0
                if keyword_set(colorbar2) then begin

                  ; bar parameters
                  ;bar_x0= 0.62
                  ;bar_y0= 0.12
                  ;bar_x1= 0.95
                  ;bar_y1= 0.17
                  bar_x0= x0+0.20
                  bar_y0= y0+0.032
                  bar_x1= x0+0.30
                  bar_y1= y0+0.047

                  ; setup bar 
                  invertcolors= 1
                  bar= REPLICATE(1B, 10)#BINDGEN(256)
                  if invertcolors eq 1 then bar= 255-bar
                  idx= where(bar eq 1) 
                  if idx(0) ne -1 then begin
                        print, "idx has ", n_elements(idx), idx(0)
                        colornumber= idx(0)-1
                        if colornumber lt 0 then colornumber= 11     ; if it's bar[*,0] then set to bar[*,1]
                        bar(idx)= bar(colornumber)
                        ;bar(idx)= 2
                  endif

                  ; do this if it's horizontal bar
                  bar= transpose(bar)

                  tv, bar, bar_x0, bar_y0, xsize=(bar_x1-bar_x0), ysize=(bar_y1-bar_y0), /normal

                  !p.position= [bar_x0,bar_y0,bar_x1,bar_y1]

                  ma= 1.0e+1
                  mi= 1.0e-4
		  numxticks= 5

                  ; puts axes
                  plot, [0],[0], psym=3, yrange=[0,1], xrange=[alog10(mi)+4.0,alog10(ma)+4.0], $
                        /normal, /noerase, color= 0, xstyle=1, ystyle=1, $
                        xthick=3.0, ythick=3.0, yticks= 1, $
                        xticks=numxticks, xtickformat='(a1)', ytickformat='(a1)', $
                        /nodata, xticklen=0.22

                  ; puts title and ticks on right-hand side
                  ;axis, yaxis=1, yrange=[alog10(mi)+4.0,alog10(ma)+4.0], ystyle=1, $                  
                  ;     xcharsize=1.50, ycharsize=1.50, charthick=3.0, $
                  ;      /normal, $
                  ;      ytitle= 'Log !4R!D!3gas!N(M!D!9n!3!N pc!E-2!N)'   ; +4.0 above comes from converting to these units

                  ; puts title on top
                  axis, xaxis=1, xrange=[alog10(mi)+4.0,alog10(ma)+4.0], xstyle=1, color= 0, $
                        xcharsize=0.95, charthick=3.0, xthick= 3.0, /normal, $
                        xticks=numxticks, xtickformat='(a1)', xticklen=0.08, $
                        ;xtitle= 'Log !4R!D!3gas!N(M!D!9n!3!N pc!E-2!N)'
                        xtitle= 'Log !4R!D!3 !N(M!D!9n!3!N pc!E-2!N)'
                        ; +4.0 above comes from converting to these units
                  ; puts ticks on bottom
                  axis, xaxis=0, xrange=[alog10(mi)+4.0,alog10(ma)+4.0], xstyle=1, $
                        xticks=numxticks,xcharsize=0.95, ycharsize=0.95, charthick=3.0, xthick= 3.0, $
                        /normal, color= 0, xticklen=0.08
                  ; y axes
                  axis, yaxis=1, yrange=[0,1], ystyle=1, ythick= 3.0, color= 0, /normal, $
                        yticks=1, ytickformat='(a1)'
                  ; puts ticks on bottom
                  axis, yaxis=0, yrange=[0,1], ystyle=1, ythick= 3.0, color= 0, /normal,  $
                        yticks=1, ytickformat='(a1)'

                endif




                ; ----------------------
                ;  ColorBar #3
                ;
                ; This one is on the bottom,
                ; horizontally, and outside
                ; the plotted area (this one
                ; expects the ysize = 14 cm)
                ; ----------------
                colorbar3= 0
                if keyword_set(colorbar3) then begin

                  ; bar parameters
                  bar_x0= 0.05
                  bar_y0= 0.085
                  bar_x1= 0.95
                  bar_y1= 0.135

                  ; setup bar 
                  invertcolors= 1
                  bar= REPLICATE(1B, 10)#BINDGEN(256)
                  if invertcolors eq 1 then bar= 255-bar
                  idx= where(bar eq 1)
                  if idx(0) ne -1 then begin
                        print, "idx has ", n_elements(idx), idx(0)
                        colornumber= idx(0)-1
                        if colornumber lt 0 then colornumber= 11     ; if it's bar[*,0] then set to bar[*,1]
                        bar(idx)= bar(colornumber)
                        ;bar(idx)= 2
                  endif

                  ; do this if it's horizontal bar
                  bar= transpose(bar)

                  tv, bar, bar_x0, bar_y0, xsize=(bar_x1-bar_x0), ysize=(bar_y1-bar_y0), /normal

                  !p.position= [bar_x0,bar_y0,bar_x1,bar_y1]

                  ma= 1.0e-1
                  mi= 1.0e-4

                  ; puts axes
                  plot, [0],[0], psym=3, yrange=[0,1], xrange=[alog10(mi),alog10(ma)], $
                        /normal, $
                        /noerase, color= 0, xstyle=1, ystyle=1, $
                        xthick=3.0, ythick=3.0, yticks= 1, $
                        xticks=3, xtickformat='(a1)', ytickformat='(a1)', $
                        /nodata, xticklen=0.02

                  ; puts title and ticks on right-hand side
                  ;axis, yaxis=1, yrange=[alog10(mi)+4.0,alog10(ma)+4.0], ystyle=1, $            
                  ;     xcharsize=1.50, ycharsize=1.50, charthick=3.0, $
                  ;      /normal, $
                  ;      ytitle= 'Log !4R!D!3gas!N(M!D!9n!3!N pc!E-2!N)'   ; +4.0 above comes from converting to these units

                  ; puts title on top
                  axis, xaxis=0, xrange=[alog10(mi)+4.0,alog10(ma)+4.0], xstyle=1, color= 0, $
                        xcharsize=1.20, ycharsize=1.20, charthick=3.0, xthick= 3.0, /normal, xticks=3, $
                        xtitle= 'Log !4R!D!3gas!N(M!D!9n!3!N pc!E-2!N)'   ; +4.0 above comes from converting to these units
                  ; puts ticks on bottom
                  axis, xaxis=1, xrange=[alog10(mi)+4.0,alog10(ma)+4.0], xstyle=1, $
                        /normal, xticks=3, xthick= 3.0, color= 0, xticklen=0.02, xtickformat='(a1)'
                  ; y axis - right
                  axis, yaxis=1, yrange=[0,1], ystyle=1, ythick= 3.0, color= 0, /normal, $
                        yticks=1, ytickformat='(a1)'
                  ; y axis - left 
                  axis, yaxis=0, yrange=[0,1], ystyle=1, ythick= 3.0, color= 0, /normal,  $
                        yticks=1, ytickformat='(a1)'

                endif





		;============================================================================
		;   999999999
		;
		;
		;   Put Velocity Vectors On
		;
		; ----------------------------
		add_vel_vectors= 0
		;add_vel_vectors= 1
		if add_vel_vectors eq 1 then begin

			oplot_velocity_vectors, "gas"

		endif




		

		;============================================================================
		; 1010101010101010101010
		;
		;
		;  Fit isophote & calc deviations
        	;
        	; --------------------------------
        	;calc_isophote_boxydisky= 1
        	calc_isophote_boxydisky= 0
        	if calc_isophote_boxydisky eq 1 then begin

           	 Thresh= 256-HalfMassSB
           	 SurfaceBrightness= Pic
           	 Scale= 2.0*xlen

           	 Xfit=fltarr(360)
           	 YFit=fltarr(360)

           	 Xcon=fltarr(360)
           	 Ycon=fltarr(360)

           	 AFit=0.0
           	 BFit=0.0
           	 X0Fit=0.0
           	 Y0Fit=0.0
           	 PhiFit=0.0
           	 a4Fit=0.0


           	 spawn, 'echo $HOME', result
           	 homedir= strcompress(result,/remove_all)
           	 libfile= homedir+'/Tools/C-Routines_for_IDL/ComputeIsoPhotFitIDL/isophotfit.so'
           	 S = CALL_EXTERNAL(libfile, $
                      'fit_isophot', $
                      PIXELS, $
                      SurfaceBrightness,$
                      Thresh,$
                      Scale, $
                      Xfit, Yfit,$
                      PhiFit, AFit, BFit, X0Fit, Y0Fit, a4Fit, $
                      Xcon, Ycon)

           	; disky/boxy
           	print, "a4/a= ", a4Fit/AFit

           	; normal coordinates
           	;e_x = x0 + (x1-x0)*1./480.0 * Xfit
           	;e_y = y0 + (y1-y0)*1./480.0 * Yfit
           	;plots, e_x, e_y, color= 150, thick= 4.0, /normal

           	loadct, 4
           	tvlct,r,g,b,/get
           	v1=[0,255]
        	   v2=[0,255]
         	  v3=[0,255]
           	tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white

           	plot,[0],[0],psym=3,xrange=[-xlen,xlen],yrange=[-xlen, xlen],xstyle=1+4,ystyle=1+4, $
                /noerase, /nodata

	           oplot, [Xfit, xfit(0)], [Yfit, yfit(0)], color=150, thick=4.0




	        endif





; -------------
;  Done
; -------------



end




