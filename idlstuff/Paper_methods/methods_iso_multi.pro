pro produce_one_image, frun, snapnum, xlen, sendto, $
			center=center, $
			colorbar=colorbar, $
			crude=crude, $
			filename=filename, $
			fitstoo=fitstoo, $
			loadedsnap=loadedsnap, $
			msg=msg, $
			nolabels=nolabels, $
			particlesonly=particlesonly, $
			pixels=pixels, $
			pubstyle=pubstyle, $
			plotpts=plotpts, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			showbhs=showbhs, $
			startid=startid, numpart=numpart, $
			thumbnail=thumbnail, $
			track_to_draw= track_to_draw, $
			xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
			xz=xz, yz=yz, $
			ReturnImage=ReturnImage, $
			ReturnTimeLbl= ReturnTimeLbl

if not keyword_set(xlen) then xlen=100.0
if not keyword_set(sendto) then sendto='ps'
if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(filename) then filename='contour.eps'
if not keyword_set(frun) then begin
   print, "  "
   print, " PROBLEM: produce_one_image"
   print, "  "
   print, "  "
   return
endif



; --------------------------------
;  Load variables for smoothing
; --------------------------------
    if not keyword_set(loadedsnap) then begin
      if (fload_snapshot(frun, snapnum)) then begin
      ;if (fload_snapshot_bh(frun, snapnum)) then begin
        print, "PROBLEM: opening file"
        return
      endif
    endif


;---------------------------------
;  Gas Particles
;---------------------------------

    npart= fload_npart(0)
    N= long(npart)

    ; default
    ;x= (1./0.7)*fload_gas_xyz('x',center=[0,0,0])
    ;y= (1./0.7)*fload_gas_xyz('y',center=[0,0,0])
    ;z= (1./0.7)*fload_gas_xyz('z',center=[0,0,0])
    ;m= (1./0.7)*fload_gas_mass(1)
    x= fload_gas_xyz('x',center=[0,0,0])
    y= fload_gas_xyz('y',center=[0,0,0])
    z= fload_gas_xyz('z',center=[0,0,0])
    m= fload_gas_mass(1)
    hsml= fload_gas_hsml(1)
    ;hsml= 0.2 + m*0.0


    ;; 1 galaxy
    ;if keyword_set(numpart) and keyword_set(startid) then begin
    ;	x= fload_1gal_gas_xyz('x',startid,numpart,center=[0,0,0])
    ;	y= fload_1gal_gas_xyz('y',startid,numpart,center=[0,0,0])
    ;	z= fload_1gal_gas_xyz('z',startid,numpart,center=[0,0,0])
    ;	m= fload_1gal_gas_mass(startid,numpart) - fload_1gal_gas_mfs(startid,numpart)
    ;endif



;--------------------------------
;  Stellar Particles
;--------------------------------

    ;ndisk= fload_npart(2)                               ; disk particles
    ;nbulge= fload_npart(3)                              ; bulge particles
    ;nstars= fload_npart(4)
    ;npart= long(ndisk) + long(nbulge) + long(nstars)
    ;print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars
    
    ;x= fload_allstars_xyz('x',center=[0,0,0])
    ;y= fload_allstars_xyz('y',center=[0,0,0])
    ;z= fload_allstars_xyz('z',center=[0,0,0])
    ;m= fload_allstars_mass(1)
    ;;hsml= fload_allstars_hsml(1)
    ;hsml= 0.2 + m*0.0


    ;x= x - center[0]
    ;y= y - center[1]
    ;z= z - center[2]


;--------------------------------
;   New Stellar Particles
;--------------------------------

    ;x= fload_newstars_xyz('x',center=[0,0,0])
    ;y= fload_newstars_xyz('y',center=[0,0,0])
    ;z= fload_newstars_xyz('z',center=[0,0,0])
    ;m= fload_newstars_mass(1)




; determine_center
; -----------------
;fancy_center= 1
fancy_center= 0
if fancy_center eq 1 then begin

	;center= [0,0,0]
	center_bh= fload_center_alreadycomp(1)

	; two black holes
	if fload_npart(5) gt 1 then begin
	        ;bhid= fload_blackhole_id(1)
	        ;bhid1= bhid[0]
	        ;bhid2= bhid[1]
	        ;center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
	        ;print, "Blackhole : ", bhid1, center1
	        ;center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
	        ;print, "Blackhole : ", bhid2, center2
	        ;center_bh= center1
	        ;;center_bh= center2
	endif

	; one black hole
	if fload_npart(5) eq 1 then begin
	        ;bhid= fload_blackhole_id(1)
	        ;center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
	        ;print, "Blackhole : ", bhid, center_bh
	endif

	center= center_bh
	print, "using this center= ", center

endif



; ----------------------------------------------
;  Call our generic contour plotting procedure
; ----------------------------------------------

;set_maxden= 10.0    ; 10^4 msolar/pc2
set_maxden= 1.0e-1
set_dynrng= 1.0e+3
;set_maxden= 0.5e-1    ; 10^4 msolar/pc2
;set_dynrng= 1.0e+4

;pixels= 256L


contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
	hsml=hsml, $
	filename=filename, fitstoo=fitstoo, $
	xthickness=xthickness, ythickness=ythickness, $
	pixels=pixels, zthickness=zthickness, $
	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	crude=crude, center=center, $
	set_maxden= set_maxden, $
        set_dynrng= set_dynrng, $
	NxNImage=NxNImage

	ReturnImage= NxNImage




	ReturnTimeLbl= fload_timelbl(1,3,/noteq)
; -------------
;  Done
; -------------



end















pro produce_multiplot, junk

;
;  This is specifically to make a six panel image
;
;
;
;


if not keyword_set(junk) then begin
   print, "  "
   print, "produce_multiplot, junk"
   print, "    (set parameters"
   print, "         from within)"
   print, "  "
   print, "  "
   return
endif


; --------------------------------
;  Center this on something other
;  than 0,0,0
; --------------------------------
center=[0.0,0.0,0.0]
msg=' '
;xlen= 500.0
xlen= 50.0
;xlen= 25.0
sendto= 'ps'



filename= 'isopanelfig.eps'

frun1= "execute/Sbc11i4-u55" & msg1= '!8n0low!6' &  snap1= 10
frun2= "execute/Sbc11i4-u50" & msg2= '!8n0med!6' &  snap2= 10
frun3= "execute/Sbc11i4-u51" & msg3= '!8n0high!6' &  snap3= 10

frun4= "execute/Sbc11i4-u56" & msg4= '!8n1low!6' &  snap4= 10
frun5= "execute/Sbc11i4-u53" & msg5= '!8n1med!6' &  snap5= 10
frun6= "execute/Sbc11i4-u52" & msg6= '!8n1high!6' &  snap6= 10

frun7= "execute/Sbc11i4-u57" & msg7= '!8n2low!6' & snap7= 10
frun8= "/data/tcox/Sbc11i4-u4" & msg8= '!8n2med!6' &  snap8= 10
frun9= "execute/Sbc11i4-u43" & msg9= '!8n2high!6' &  snap9= 10








;--------------------------------------
;  Now plot this mess (to postscript)
;--------------------------------------


        initialize_plotinfo, 1

	;setup_plot_stuff, sendto, filename=filename, colortable= 4
        ;setup_plot_stuff, sendto, filename=filename, colortable= 4, newxsize=18, newysize=24
        ;setup_plot_stuff, sendto, filename=filename, colortable= 3, newxsize=18, newysize=24
        ;setup_plot_stuff, sendto, filename=filename, colortable= 1, newxsize=18, newysize=24
        setup_plot_stuff, sendto, filename=filename, colortable= 1, newxsize=18, newysize=18
        ;setup_plot_stuff, sendto, filename=filename, colortable= 0, newxsize=18, newysize=24


	x0= 0.008
	xs= 0.328   ; assumes 3 panels
	x1= x0+xs
	x2= x0+xs+xs
	x3= x0+xs+xs+xs

        y0= 0.008
	ys= 0.328   ; assumes 3 panels
        y1= y0+ys
	y2= y0+ys+ys
	y3= y0+ys+ys+ys



	;  place the image(s) down
	; ----------------------

	do_panel, frun1, snap1, xlen, x0, y2, xs, ys, msg=msg1
	do_panel, frun2, snap2, xlen, x1, y2, xs, ys, msg=msg2
	do_panel, frun3, snap3, xlen, x2, y2, xs, ys, msg=msg3

	do_panel, frun4, snap4, xlen, x0, y1, xs, ys, msg=msg4
	do_panel, frun5, snap5, xlen, x1, y1, xs, ys, msg=msg5
	do_panel, frun6, snap6, xlen, x2, y1, xs, ys, msg=msg6

	do_panel, frun7, snap7, xlen, x0, y0, xs, ys, msg=msg7
	do_panel, frun8, snap8, xlen, x1, y0, xs, ys, msg=msg8
	do_panel, frun9, snap9, xlen, x2, y0, xs, ys, msg=msg9



	; print msg?
	; -----------
        ;if keyword_set(msg) then begin
        ;        ;xyouts, 0.07, 0.91, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
	;	xyouts, 0.07, 0.83, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
        ;endif



        ; done, close this up
        ; --------------------
	device, /close






; -------------
;  Done
; -------------



end







; does the work for each panel
; -------------------------------
pro do_panel, frun, snapnum, xlen, $
	x0, y0, xs, ys, $
        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	colorbar=colorbar, draw_circle=draw_circle, $
	msg=msg



	;  get images
	; -------------
	center=[0,0,0]
	produce_one_image, frun, snapnum, xlen, xz=xz, yz=yz, $
		xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
		pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
		crude=crude, center=center, ReturnImage=ReturnImage, ReturnTimeLbl= ReturnTimeLbl

	Pic= ReturnImage
	TimeLbl= ReturnTimeLbl

	if keyword_set(msg) then TimeLbl= ''


	tv, Pic, x0, y0, xsize=xs, ysize=ys, /normal


	x1= x0 + xs
	y1= y0 + ys
	overplot_axes, x0, y0, x1, y1, TimeLbl, xlen    ;, /add_j_vector


	;snaplbl= '('+strcompress(string(snapnum),/remove_all)+')'
        ;xyouts, x0+0.03, y0+ys-0.05, snaplbl, /normal, size= 1.0, charthick=3.0, color= 0


	if keyword_set(msg) then begin
		xyouts, x0+0.02, y0+ys-0.04, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
		;print, "printing msg= ", msg
	endif

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


		if keyword_set(track_to_draw) then begin
			xtrack= track_to_draw[*,Axis1]
			ytrack= track_to_draw[*,Axis2]
			idx= where((xtrack ne 0) or (ytrack ne 0))
			if idx(0) ne -1 then begin
			  xtrack= xtrack(idx)
			  ytrack= ytrack(idx)
                          oplot, xtrack, ytrack, color= 0, thick=3.0, linestyle= 2, psym=-5
			endif
		endif


	    ; generate white axes if desired
	    ;axis, xaxis=0, xstyle=1,ystyle=1,color=1,charsize=0.001, xthick=2.0,ythick=2.0
            ;axis, xaxis=1, xstyle=1,ystyle=1,color=1,charsize=0.001, xthick=2.0,ythick=2.0
            ;axis, yaxis=0, xstyle=1,ystyle=1,color=1,charsize=0.001, xthick=2.0,ythick=2.0
            ;axis, yaxis=1, xstyle=1,ystyle=1,color=1,charsize=0.001, xthick=2.0,ythick=2.0


		; -----------------------------
		;  X marks the black hole spot
		; -----------------------------
		showbhs= 0
		if keyword_set(showbhs) and fload_npart(5) ge 1 then begin

                        loadct, 4          ; blue/green/red/yellow (std)
                	tvlct,r,g,b,/get
                	v1=[0,255]
                	v2=[0,255]
                	v3=[0,255]
                	tvlct,v1,v2,v3,0

			n_bh= fload_npart(5)
			bhid= 0.0; fload_blackhole_id(1)


			if n_bh gt 1 then begin

				bhid1= bhid[0]
				bh_center_1= 0.0; fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
				orig_bh_center_1= bh_center_1
				bh_center_1= bh_center_1 - bh_center_1
				oplot, [bh_center_1[0]], [bh_center_1[1]], psym=7, thick=6.0, color=150, symsize=2.0
				print, "Blackhole #1: ",bhid1, bh_center_1

				bhid2= bhid[1]
				bh_center_2= 0.0; fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
				bh_center_2= bh_center_2 - orig_bh_center_1
				oplot, [bh_center_2[0]], [bh_center_2[1]], psym=7, thick=6.0, color=150, symsize=2.0
				print, "Blackhole #2: ",bhid2, bh_center_2
				
			endif else begin
				bhid= fload_blackhole_id(1)
				bh_center_1= 0.0; fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
				bh_center_1= bh_center_1 - bh_center_1
				oplot, [bh_center_1[0]], [bh_center_1[1]], psym=7, thick=6.0, color=150, symsize=2.0
				print, "Blackhole: ", bhid, bh_center_1
			endelse

			loadct, 1          ; blue/white
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0


		endif




        	; ----------------------------
        	;
        	; Add angular momentum vector.
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



                ; -----------------------------
                ;  Draw center path, similar to
		;  tracktodraw
                ; -----------------------------
		;drawcenters= 1
		if keyword_set(drawcenters) then begin

			; warning, need to load process_directory for this
			read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=fload_frun(1)+"/centers.txt"
			cp_cen2_x= cp_cen2[0,*]
			cp_cen2_y= cp_cen2[1,*]
			cp_cen2_z= cp_cen2[2,*]
			idx= where(cp_time le fload_time(1))
			cp_cen2_x= cp_cen2_x(idx)
			cp_cen2_y= cp_cen2_y(idx)

			oplot, cp_cen2_x, cp_cen2_y, psym=-3, linestyle=0, color= 150, thick=4.0


			; double up baby
			;read_centerpositions, cp_time, cp_cen1, cp_cen2, filename="Data/nodrag/centers.txt"
			;cp_cen2_x= cp_cen2[0,*]
                        ;cp_cen2_y= cp_cen2[1,*]
			;idx= where(cp_time le fload_time(1))
			;cp_cen2_x= cp_cen2_x(idx)
                        ;cp_cen2_y= cp_cen2_y(idx)
                        ;oplot, cp_cen2_x, cp_cen2_y, psym=-3, linestyle=0, color= 100, thick=4.0


		endif




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






        ; do a contour plot of this mess
        ; ---------------------------------
	docontour= 0
        if docontour EQ 1 then begin
;               for i=0, (n_elements(info)-1) do begin
;                       cnt_index= [indgen(info(i).N),0]
;                       plots, xy_contours(*,info(i).offset+cnt_index), /normal, color= 0
;               endfor

		levels = 6
		step = (Max(Pic) - Min(Pic)) / levels
		userLevels = IndGen(levels) * step + Min(Pic)

		; load contours into variables (this doesn't print)
		contour, Pic, path_xy=xy_contours, path_info=info

                for i=0, (n_elements(info)-1) do begin
                        cnt_index= [indgen(info(i).N),0]
                        plots, xy_contours(*,info(i).offset+cnt_index), /normal, color= 0
                endfor
        endif



	; overplot the new stars, if we want
	; ------------------------------------
	;do_new_star_overplotting= 1
	do_new_star_overplotting= 0
	if do_new_star_overplotting eq 1 then begin
		print, "getting to overplot"
		overplot_all_newstars, frun, x0, xs, y0, ys, xmin, xmax, ymin, ymax, xlen, $
				snap1, snap2, snap3, snap4, snap5, snap6, xz=xz, $
				panel6_seen_edgeon=panel6_seen_edgeon
	endif





                ; ----------------------
                ;  ColorBar #2
                ;
                ; This one is on the bottom,
                ; horizontally, and inside
                ; the plotted area (for
                ; nolabels keyword)
                ; ----------------
                if keyword_set(colorbar) then begin

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




	; print msg?
	; -----------
        ;if keyword_set(msg) then begin
        	;xyouts, x0+0.03, y0+ys-0.05, snaplbl, /normal, size= 1.0, charthick=3.0, color= 0
		;xyouts, x0+(0.07*xs), y0-(0.12*xs), msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
	;	xyouts, x0+0.03, y0-0.08, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
        ;endif




; -------------
;  Done
; -------------



end




pro overplot_all_newstars, frun, x0, xs, y0, ys, xmin, xmax, ymin, ymax, xlen, $
			snap1, snap2, snap3, snap4, snap5, snap6, xz=xz, $
			panel6_seen_edgeon=panel6_seen_edgeon

	x1= x0+xs
	x2= x0+xs+xs
	x3= x0+xs+xs+xs

        y1= y0+ys
	y2= y0+ys+ys


        ;  1
        ; ---
	overplot_newstars, frun, snap1, x0, y1, x1, y2, xlen, xmin, xmax, ymin, ymax, $
                                xz=xz, /do_newstar_contour

        ;  2
        ; ---
	overplot_newstars, frun, snap2, x1, y1, x2, y2, xlen, xmin, xmax, ymin, ymax, $
                                xz=xz, /do_newstar_contour

        ;  3
        ; ---
	overplot_newstars, frun, snap3, x2, y1, x3, y2, xlen, xmin, xmax, ymin, ymax, $
                                xz=xz, /do_newstar_contour

        ;  4
        ; ---
	overplot_newstars, frun, snap4, x0, y0, x1, y1, xlen, xmin, xmax, ymin, ymax, $
                                xz=xz, /do_newstar_contour

        ;  5
        ; ---
	overplot_newstars, frun, snap5, x1, y0, x2, y1, xlen, xmin, xmax, ymin, ymax, $
				xz=xz, /do_newstar_contour

        ;  6
        ; ---
	if panel6_seen_edgeon eq 1 then xz= 1
	overplot_newstars, frun, snap6, x2, y0, x3, y1, xlen, xmin, xmax, ymin, ymax, $
				xz=xz, /do_newstar_contour


end






;-----------------------------------------
; Overplot New Stars
;-----------------------------------------
pro overplot_newstars, frun, snapnum, x0, y0, x1, y1, xlen, xmin, xmax, ymin, ymax, $
			do_newstar_contour=do_newstar_contour, $
			rotate_phi=rotate_phi, $
			rotate_theta=rotate_theta, $
			center=center, $
                        fitstoo=fitstoo, $
			xz=xz, yz=yz, $
                        xthickness=xthickness, $
			ythickness=ythickness, $
                        zthickness=zthickness, $
                        pixels=pixels

	center= [0,0,0]

	; --------------------------------
	;  Load variables for smoothing
	; --------------------------------
	if not keyword_set(loadedsnap) then begin
	   ;if (fload_snapshot(frun, snapnum)) then begin
	   if (fload_snapshot_bh(frun, snapnum)) then begin
		print, "PROBLEM: opening file"
		return
	   endif
	endif


        npart= fload_npart(4)                         ;newstars  particles
        print, "Ntot= ",npart
        N= long(npart)

        x= fload_newstars_xyz('x',center=[0,0,0])
        y= fload_newstars_xyz('y',center=[0,0,0])
        z= fload_newstars_xyz('z',center=[0,0,0])
        m= fload_newstars_mass(1)


        loadct, 4               ; blue/green/red/yellow (std)
        tvlct,rr,gg,bb,/get
        vv1=[0,255]
        vv2=[0,255]
        vv3=[0,255]
        tvlct,vv1,vv2,vv3,0


	; plot axes
	; ----------
        !p.position=[x0,y0,x1,y1] & !p.ticklen=0.03
        plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.00, ycharsize=0.00, xstyle=1, ystyle=1, /normal, /nodata, $
                      xthick=0.0, ythick=0.0, xtickformat='(a1)', ytickformat='(a1)'

        ; plot the actual points
	; ------------------------
        if keyword_set(yz) then begin
              x= y
              y= z
        endif
   
        if keyword_set(xz) then y= z


	; actually overplot the new stars
	; ----------------------------------
	if not keyword_set(do_newstar_contour) then begin
	  ;!p.position=[x0,y0,x1,y1]
	  ;plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	  ;	/NOERASE, color= 0, xcharsize=0.01, ycharsize=0.01, $
	  ;	xthick=0.01, ythick=0.01, xstyle=1, ystyle=1, /normal, $
	  ;	xtickformat='(a1)', ytickformat='(a1)', /nodata

          ;oplot, x, y, psym=5, thick=2.0, color= 50
	  ;oplot, x, y, psym=2, thick=2.0, color= 150
	  ;oplot, x, y, psym=7, thick=2.0, color= 100
	  ;oplot, x, y, psym=1, thick=2.0, color= 0
	  oplot, x, y, psym=3, color= 150
	endif else begin

	        contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
                        fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, $
                        pixels=pixels, zthickness=zthickness, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        crude=crude, center=center, $
                        NxNImage=NxNImage

	        NewStarPic= NxNImage

                levels = 6
                ;step = (Max(XrayPic) - Min(XrayPic)) / levels
                step = (256)/levels
                ;userLevels = IndGen(levels) * step + Min(XrayPic)
                userLevels = IndGen(levels) * step
                ;userLevels = [0,userLevels[2:5]]
                userLevels = userLevels[2:5]
                ;userLevels = [userLevels[1],userLevels[5]]
                print, "userLevels= ",userLevels
                
                ; invert Pic, so that high is central bright regions
                ;NewStarPic= 255-NewStarPic
                
                ; load contours into variables (this doesn't print)
                ;contour, NewStarPic, path_xy=xy_contours, path_info=info
                ;contour, NewStarPic, path_xy=xy_contours, path_info=info, levels=userLevels
                contour, NewStarPic, path_xy=xy_contours, path_info=info, levels=userLevels, $
                                        min_value=2, /PATH_DATA_COORDS 
                                        
                print, "N_contours= ", n_elements(info)

                for i=0, (n_elements(info)-1) do begin
                        cnt_index= [indgen(info(i).N),0]
                        xycnt= xy_contours(*,info(i).offset+cnt_index)
                        n_level= info(i).level
                        xycnt[0,*]= xycnt[0,*]*(x1-x0)/480.0 + x0
                        xycnt[1,*]= xycnt[1,*]*(y1-y0)/480.0 + y0
                        ;if n_level eq 0 then plots, xycnt, /normal, color= 150, thick=4.0
                        ;if n_level eq 1 then plots, xycnt, /normal, color= 150, thick=2.0
                        ;if n_level eq 2 then plots, xycnt, /normal, color= 150, linestyle= 1
                        ;if n_level eq 3 then plots, xycnt, /normal, color= 130, linestyle= 1
                        ;if n_level eq 4 then plots, xycnt, /normal, color= 130
                        ;plots, xy_contours(*,info(i).offset+cnt_index), /normal, color= 150
			;if n_level gt 0 then plots, xycnt, /normal, color= 150
			;plots, xycnt, /normal, color= 150
			idx= where(xycnt[0,*] eq x0)     ; get rid of box at x0,x1,y0,y1
			;if idx(0) eq -1 then plots, xycnt, /normal, color= 150
			if idx(0) eq -1 then plots, xycnt, /normal, color= 220
                endfor
	endelse


; -------------
;  Done
; -------------


end






pro overplot_axes, xx0, yy0, xx1, yy1, thetimelbl, xlen, $
			extramsg=extramsg, add_j_vector=add_j_vector, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta

        !p.position=[xx0,yy0,xx1,yy1]
	!p.ticklen=0.03

	xmin= -xlen
	xmax=  xlen
	ymin= -xlen
	ymax=  xlen

        plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
                      xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'

        xyouts, xx0+0.02, yy1-0.03, '!6'+thetimelbl, /normal, size= 1.2, charthick=3.0, color= 0

	if keyword_set(extramsg) then begin
		xyouts, xx0+0.02, yy1-0.10, extramsg, /normal, size= 1.2, charthick=3.0, color= 0
	endif


end



