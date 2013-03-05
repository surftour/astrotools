
;===================================================================================================



; overplot the velocity image
; -----------------------------
pro process_kin_velimage_and_overplot, x, y, z, vx, vy, vz, mass, xlen, $
			x0, y0, x1, y1, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			center=center, showbhs=showbhs, $
			msg1=msg1, msg2=msg2, fitstoo=fitstoo, filename=filename, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			pa=pa, show_maj_slit=show_maj_slit, pixels=pixels, $
			msg=msg, drawbox=drawbox, xslit=xslit, yslit=yslit, $
			velocitymap=velocitymap, dispersionmap=dispersionmap, $
			ContourMap=ContourMap, HalfMassSB=HalfMassSB, $
			noxaxistitle=noxaxistitle, colorbar_onright=colorbar_onright

	m= mass
	if keyword_set(center) then orig_center= center else orig_center=[0,0,0]

	print, "====================================="
	print, "  process_kin_velimage_and_overplot  "
	print, "====================================="

	vx_orig= vx
	vy_orig= vy
	vz_orig= vz

	if keyword_set(dispersionmap) then begin
		if not keyword_set(set_maxden) then set_maxden= 200.0
		if not keyword_set(set_dynrng) then set_dynrng= 100.0

		vx= vx_orig
		vy= vy_orig
		vz= vz_orig
		center= orig_center

       		contour_makevelpic, x, y, z, vx, vy, vz, m, xlen, xz=xz, yz=yz, $
                        filename=filename, fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, $
                        pixels=pixels, zthickness=zthickness, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        center=center, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
			dispersionmap=dispersionmap, $
                        NxNImage=NxNImage

		colorbarcolor= 33
		fload_newcolortable, colorbarcolor
        
	endif


        if keyword_set(velocitymap) then begin
                if not keyword_set(set_maxden) then set_maxden= 150.0
                if not keyword_set(set_dynrng) then set_dynrng= 2.0*set_maxden

		vx= vx_orig
		vy= vy_orig
		vz= vz_orig
		center= orig_center

                contour_makevelpic, x, y, z, vx, vy, vz, m, xlen, xz=xz, yz=yz, $
                        filename=filename, fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, $
                        pixels=pixels, zthickness=zthickness, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        center=center, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        velocitymap=velocitymap, $
                        NxNImage=NxNImage

		colorbarcolor= 25
		fload_newcolortable, colorbarcolor

        endif


        
       	Pic= NxNImage

	pixels= (size(NxNImage))[1]
	print, "pixels= ", pixels


    if n_elements(Pic) le 0 then begin
	print, "  "
	print, "  WARNING: no map set"
	print, "  "
	return
    endif


    ;print, "center= ", center
    ;print, "center (original)= ", orig_center
    center= [0,0,0]

    xmin= -xlen+center[0]
    xmax=  xlen+center[0]
    ymin= -xlen+center[1]
    ymax=  xlen+center[1]


    tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


    !p.position=[x0,y0,x1,y1]
    !p.ticklen=0.03

    ; creates axes and plot style
    ; ---------------------------
    if x0 lt 0.2 then begin    ; in this case, print x tick labels
	plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, /nodata, color= 0,  $
		      xstyle=1, ystyle=1, xtitle=xtit, ytitle=ytit, $
                      xcharsize=0.01, ycharsize=0.01, xthick=4.0, ythick=4.0, charthick=3.0, $
		      ;xticks=10, $
		      ;yticks=10, $
		      ytickformat='(a1)'
    endif else begin
	plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, /nodata, color= 0,  $
		      xstyle=1, ystyle=1, xtitle=xtit, ytitle=ytit, $
                      xcharsize=0.01, ycharsize=0.01, xthick=4.0, ythick=4.0, charthick=3.0, $
		      ;xticks=10, $
		      ;yticks=10, $
		      xtickformat='(a1)', ytickformat='(a1)'
    endelse


     ;xyouts, x0+0.04, y1-0.08, fload_timelbl(1,2,/noteq), /normal, size= 1.7, charthick=3.0, color= 0
     xyouts, x0+0.02, y1-0.05, msg1, /normal, size= 1.3, charthick=4.0, color= 0
     xyouts, x0+0.02, y1-0.10, msg2, /normal, size= 1.3, charthick=4.0, color= 0




                ; -----------------------------
                ;  X marks the black hole spot
                ; -----------------------------
                if keyword_set(showbhs) and fload_npart(5) ge 1 then begin

                        loadct, 4          ; blue/green/red/yellow (std)
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0

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



		; ---------------------------
		;  Stellar Density Contour
		; ---------------------------
		oplot_contours, ContourMap, x0, y0, x1, y1, hmsb=HalfMassSB, xlen=xlen, pixels= pixels

		; now turn this off!
		;ContourMap= 0   - bad, we loose our map information
		;if keyword_set(ContourMap) then begin
                ;    levels = 12
                ;    ;step = (Max(XrayPic) - Min(XrayPic)) / levels
                ;    step = (256)/levels
                ;    ;userLevels = IndGen(levels) * step + Min(XrayPic)
                ;    userLevels = IndGen(levels) * step
                ;    ;userLevels = [0,userLevels[2:5]]
                ;    ;userLevels = [userLevels[3],userLevels[5:11]]
                ;    userLevels = userLevels[2:11]
                ;    print, "userLevels= ",userLevels
;
		;    ; parse map
		;    FixedMapIdx= where(ContourMap eq 1)
		;print, "Max ContourMap= ", max(ContourMap)
		;    FixedMap= ContourMap * (256.0/max(ContourMap))
		;    FixedMap= 256-FixedMap
		;    if FixedMapIdx(0) ne -1 then FixedMap(FixedMapIdx)= 1

		;    ; smooth the map
		;    ;NewMap= min_curve_surf(ContourMap, /regular)
		;    NewMap= smooth(FixedMap, 10, /edge_truncate)

                ;    ; load contours into variables (this doesn't print)
                ;    ;contour, ContourMap, path_xy=xy_contours, path_info=info, levels=userLevels, $
                ;    contour, NewMap, path_xy=xy_contours, path_info=info, levels=userLevels, $
                ;                        min_value=2, /PATH_DATA_COORDS

                    ;loadct, 4
                    ;tvlct,r1,g1,b1,/get
                    ;v12=[0,255]
                    ;v22=[0,255]
                    ;v32=[0,255]
                    ;tvlct,v12,v22,v32,0        ; define colors 0 and 1 as black and white

                ;    print, "N_contours= ", n_elements(info)
                
                ;    for i=0, (n_elements(info)-1) do begin
                ;        cnt_index= [indgen(info(i).N),0]
                ;        xycnt= xy_contours(*,info(i).offset+cnt_index)
                ;        n_level= info(i).level
;print, "i= ",i,"  n_level= ",n_level
                ;        xycnt[0,*]= xycnt[0,*]*(x1-x0)/480.0 + x0
                ;        xycnt[1,*]= xycnt[1,*]*(y1-y0)/480.0 + y0
                ;        idx= where(xycnt[0,*] eq x0)
                ;        if idx(0) eq -1 then begin
                ;           ;if n_level eq 0 then plots, xycnt, /normal, color= 150, thick=4.0
                ;           ;if n_level eq 1 then plots, xycnt, /normal, color= 150, thick=3.0
                ;           ;if n_level eq 2 then plots, xycnt, /normal, color= 150, thick=2.0
                ;           ;if n_level eq 3 then plots, xycnt, /normal, color= 150, thick=1.0, linestyle= 1
                ;           ;if n_level eq 4 then plots, xycnt, /normal, color= 130
		;	   ;if (n_level ge 0 and n_level le 4) then plots, xycnt, /normal, color= 0, thick=2.0
		;	   if n_level gt 0 then plots, xycnt, /normal, color= 0, thick=2.0
                ;           ;plots, xy_contours(*,info(i).offset+cnt_index), /normal, color= 150
                ;        endif
                ;    endfor

		;endif

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

                ;print, "xyslit= ",xslit, yslit



                ; rotate just x and y so that it's
                ; along the major axis position angle (pa)
                ; ----------------------------------------
                if keyword_set(pa) and keyword_set(show_maj_slit) then begin
                   print, "velocity map slit pa= ",pa
                   ; actually the pa passed is clockwise rotation, so we
                   ; need to flips this around
                   ;if pa gt 0.0 then pa= 180.0 - pa else pa= -pa
                   par = pa * !PI/180.0

                   xtr= xslit*cos(par) - yslit*sin(par)       ; this choice of signs rotates
                   ytr= xslit*sin(par) + yslit*cos(par)       ; counter-clockwise

                   xbr= xslit*cos(par) - (-yslit)*sin(par)
                   ybr= xslit*sin(par) +(-yslit)*cos(par)

                   xtl= (-xslit)*cos(par) - yslit*sin(par)
                   ytl= (-xslit)*sin(par) + yslit*cos(par)

                   xbl= (-xslit)*cos(par) - (-yslit)*sin(par)
                   ybl= (-xslit)*sin(par) + (-yslit)*cos(par)

                   ; need axes to be reset

    		   fload_newcolortable, 4
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



                ; -----------------------
                ;  ColorBar #1
                ; 
                ; The original colorbar.
                ; Oriented on the right
                ; side of the plot, keyword
                ; also affects the plot
                ; positioning above.
                ; ----------------
		;colorbar= 1
                ;if keyword_set(colorbar) then begin


                if keyword_set(colorbar_onright) then begin

    		fload_newcolortable, colorbarcolor


                invertcolors= 1    ; warning: this should match the value in contour_makepic
		if keyword_set(dispersionmap) then begin
		   ;set_maxden= 200.0
		   ;set_dynrng= 100.0
		endif
		if keyword_set(velocitymap) then begin
		   ;set_maxden= 120.0
		   ;set_dynrng= 2.0*set_maxden
		endif
                MaxDensXY= set_maxden
                DynRange= set_dynrng
                ma=MaxDensXY
                mi=MaxDensXY-DynRange

                  bar= REPLICATE(1B, 10)#BINDGEN(256)
                  if invertcolors eq 1 then bar= 255-bar
                  idx= where(bar eq 0 or bar eq 1)
                  if idx(0) ne -1 then begin
                        ;print, "idx has ", n_elements(idx), idx(0)
                        ;colornumber= idx(0)-1
                        ;if colornumber lt 0 then colornumber= 11     ; if it's bar[*,0] then set to bar[*,1]
                        bar(idx)= 2   ; bar(colornumber)
                  endif
                  barwidth= 0.025

                  tv, bar, x1+0.01, y0+0.04, xsize=barwidth, ysize=(y1-y0-0.08), /normal

                  !p.position= [x1+0.01,y0+0.04,x1+0.01+barwidth,y1-0.04]
		  !p.ticklen= 0.20

                  plot, [0],[0], psym=3, xrange=[0,1], yrange=[mi,ma], $
                        /normal, $
                        /noerase, color= 0, xstyle=1, ystyle=1, $
                        xthick=4.0, ythick=4.0, $
                        xticks=1, xtickformat='(a1)', ytickformat='(a1)', $
                        /nodata

		  if keyword_set(dispersionmap) then ytit= '!7r!6 (km s!E-1!N)'
		  if keyword_set(velocitymap) then ytit= '!6Velocity (km s!E-1!N)' 
                  axis, yaxis=1, yrange=[mi,ma], ystyle=1, /normal, $
                        ;ycharsize=0.80, charthick=3.0, ytitle= ytit   ;, $ 
                        ycharsize=1.1, charthick=3.0, ytitle= ytit   ;, $ 
                        ;ytickformat='(a1)', ymargin= [1,6]
                  ;axis, yaxis=1, yrange=[mi,ma], ystyle=1, /normal, $
                        ;ycharsize=0.80, charthick=2.0
                endif




end










;#########################################################################
