pro oplot_velocity_vectors, comptype, xlen, xz=xz, yz=yz, $
			center=center, msg=msg, plotpts=plotpts, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta



		; velocity field variables
		;min_velocity= 50.0
		min_velocity= 0.0
		;draw_a_key= 1
		draw_a_key= 0
        	GRIDPTS= 20L
		zthickness= 3.0
        	normlen= 2.0
        	VNorm= 200.0  ; normalize to 200  (then divide by 10, so 1 unit is 10 km/s on graph
        	;VNorm= 100.0

		if keyword_set(center) then center=center else center= [0,0,0]


                xmin= -xlen
                xmax=  xlen
                ymin= -xlen
                ymax=  xlen

		; default is gas
		;if comptype eq "gas" then begin
                	x= fload_gas_xyz('x',center=[0,0,0])
                	y= fload_gas_xyz('y',center=[0,0,0])
                	z= fload_gas_xyz('z',center=[0,0,0])
                	m= fload_gas_mass(1)
                	hsml= fload_gas_hsml(1)
                	vx= fload_gas_v('x')
                	vy= fload_gas_v('y')
                	vz= fload_gas_v('z')
		;endif

		;vx= fload_gas_v('x')
		;vy= fload_gas_v('y')
		;vz= fload_gas_v('z')
		;vx= fload_allstars_v('x')
		;vy= fload_allstars_v('y')
		;vz= fload_allstars_v('z')

		if keyword_set(xz) then begin
			vy=vz
			z=-y
		endif

		;
		; unsure which one this is
		;
		;contour_makevelpic, x, y, z, vx, vy, vz, m, xlen, $
		contour_velmap, x, y, z, vx, vy, vz, m, xlen, $
                        ;hsml=hsml, $
                        pixels=GRIDPTS, zthickness=zthickness, $
                        ;pixels=pixels, zthickness=zthickness, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        center=center, $
                        VelxMap= VelxMap, $
                        VelyMap= VelyMap


		; this needed resetting in one instance
		;!p.position=[x0,y0,x1,y1]
                ;!p.ticklen=0.03
                ;plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                ;      xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
                ;      xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'




	        ; draw velocity field
	        ; ------------------------
	        ;pix_int= fix(XYPIXELS/(gridpts+1))
	        pix_int= 1.0
	        vx= fltarr(2)
	        vy= fltarr(2)

		do_manual_vectors= 1
		if do_manual_vectors eq 1 then begin
	          for i=0,gridpts-1 do begin
	            for j=0,gridpts-1 do begin

			; arrow head size
			;thishsize= -0.25
			thishsize= -0.425

	                ; physical location
	                vx[0]= (xmax-xmin)*(i+0.5)/gridpts + xmin
	                vy[0]= (ymax-ymin)*(j+0.5)/gridpts + ymin

	                ; add on velocity at that point
	                ;vx[1]= vx[0] + VelxMap[(i+1)*pix_int,(j+1)*pix_int]/VNorm
	                ;vy[1]= vy[0] + VelyMap[(i+1)*pix_int,(j+1)*pix_int]/VNorm
	                vx[1]= vx[0] + VelxMap[i,j]/VNorm
	                vy[1]= vy[0] + VelyMap[i,j]/VNorm

			if vx[1] gt xmax then begin
				vx[1]= xmax   & thishsize= 0
			endif
			if vx[1] lt xmin then begin
				vx[1]= xmin   & thishsize= 0
			endif
			if vy[1] gt ymax then begin
				vy[1]= ymax   & thishsize= 0
			endif
			if vy[1] lt ymin then begin
				vy[1]= ymin   & thishsize= 0
			endif

	                vellen= sqrt((vx[1]-vx[0])*(vx[1]-vx[0]) + (vy[1]-vy[0])*(vy[1]-vy[0]))
	                ;if vellen gt 0.0 then begin
	                if vellen gt (min_velocity/VNorm) then begin
	                    ;arrow, vx[0], vy[0], vx[1], vy[1], /data, color= 1
	                    ;arrow, vx[0], vy[0], vx[1], vy[1], /data, color= 0, hsize=thishsize
	                    arrow, vx[0], vy[0], vx[1], vy[1], /data, color= 0, hsize=thishsize, thick=2.0
	                endif

	            endfor
	          endfor
		endif else begin
		  dpt= (xmax-xmin)/GRIDPTS/2.0
		  xs= (xmax-xmin)*findgen(GRIDPTS)/(GRIDPTS) + xmin + dpt
		  VELOVECT, VelxMap, VelyMap, xs, xs, /overplot, color= 0 ;[, X, Y] [, COLOR=index] [, MISSING=value [, /DOTS]] [, LENGTH=value] [, /OVERPLOT]
		endelse


	        ;  Draw a key
        	; -----------------
        	if draw_a_key eq 1 then begin

        	   ; big images.
        	   ;bx0=xmin+0.9
        	   ;bx1=xmin+7.1
        	   ;by0=ymax-1.0
        	   ;by1=ymax-3.1
        	   bx0=xmin+0.9
        	   bx1=xmin+17.1
        	   by0=ymin+10.0
        	   by1=ymin+1.1
        	   polyfill, [bx0,bx0,bx1,bx1], [by0,by1,by1,by0], color= 1
        	   oplot, [bx0,bx0,bx1,bx1,bx0], [by0,by1,by1,by0,by0] , color= 0, thick=3.0

        	   ;xyouts, bx0+1.2, by1-4.8, fload_timelbl(0.7,2,/round_off), /data, size= 1.3, charthick=4.0, color= 0
		   windlbl= '!7g!6=2.0'
		   ;windlbl= '!7g!6=0.5'
		   ;windlbl= '!7g!6=0.05'
        	   ;xyouts, bx0+0.2, by0-0.7, windlbl, /data, size= 1.3, charthick=4.0, color= 0
        	   xyouts, bx0+0.5, by0-1.7, windlbl, /data, size= 1.3, charthick=4.0, color= 0
		   ;windlbl= '!8v!6!Dw!N=837 km s!E-1!N'
		   windlbl= '!8v!6!Dw!N=105 km s!E-1!N'
        	   ;xyouts, bx0+0.2, by0-1.7, windlbl, /data, size= 1.3, charthick=4.0, color= 0
        	   xyouts, bx0+0.5, by0-7.7, windlbl, /data, size= 1.3, charthick=4.0, color= 0

        	   ; add a little key
        	   ;keylen= 200.0
        	   ;keylen= keylen/VNorm

        	   ;arrow, bx0+1.0, by1-2.0, bx0+1.0+keylen, by1-2.0, /data, thick=3.0, hthick=2.0, color= 0
        	   ;xyouts, bx0+1.0, by1-7.0, '200 km s!E-1!N', /data, size= 1.3, charthick=4.0, color= 0
        	   ;xyouts, bx0+1.0, by1-12.0, fload_timelbl(h,2), /data, size= 1.3, charthick=4.0, color= 0
        	endif




end


