pro contour_makeplot, x, y, z, m, hsml, xlen, xz=xz, yz=yz, $
			filename=filename, thumbnail=thumbnail, fitstoo=fitstoo, $
			xthickness=xthickness, ythickness=ythickness, colorbar=colorbar, $
			pixels=pixels, zthickness=zthickness, $
			crude=crude, center=center, msg=msg, plotpts=plotpts, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			pubstyle=pubstyle, particlesonly=particlesonly, $
			showbhs=showbhs, set_maxden=set_maxden, set_dynrng=set_dynrng, $
			show_rsep=show_rsep, colortable=colortable, $
			track_to_draw=track_to_draw, nolabels=nolabels, $
			drawbox=drawbox, xslit=xslit, yslit=yslit


;
;  This is our general contour plotting procedure, note that this will
; make matters less confusing since we can simply plot everything from here now.
;
;
;
;
;
;
;
;
;


if not keyword_set(xlen) then xlen=100.0
if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(colortable) then colortable= 4
if not keyword_set(x) then begin
   print, "  "
   print, "PROBLEM: contour_makeplot"
   print, "  "
   print, "  "
   print, "  "
   return
endif


if keyword_set(center) then begin
	orig_center= center
endif else begin
	orig_center=[0,0,0]
	center=[0,0,0]
endelse




;--------------------------------------
;  Make image
;--------------------------------------

if not keyword_set(particlesonly) then begin

;  old code to make the NxN image
; ---------------------------------
use_old_smoothing= 1
if use_old_smoothing eq 1 then begin
	contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
			hsml=hsml, $
                        filename=filename, fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, $
                        pixels=pixels, zthickness=zthickness, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        crude=crude, center=center, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			NxNImage=NxNImage, HalfMassSB=HalfMassSB

	Pic= NxNImage
	PIXELS= n_elements(Pic[0,*])
endif

;  brant's surface brightness code
; ----------------------------------
;use_brants_sbcode= 0
;if use_brants_sbcode eq 1 then begin
;	N= n_elements(x)
;	npixels= 480
;	Pic= brants_surfacebrightness_map(N,x,y,z,m,hsml,xlen, $
;				npixels=npixels, $
;                                rotate_phi=rotate_phi, $
;				rotate_theta=rotate_theta, $
;				set_maxden=set_maxden, $
;				set_dynrng=set_dynrng, $
;				fitstoo=fitstoo, $
;				filename=filename, $
;                                center=center)
;endif


; particlesonly
endif




;-------------------------------------------------------------------------



    ; default is xy axis
    ; ------------------
    xtit="!8x!3 (kpc/h)"
    ytit="!8y!3 (kpc/h)"
    ;xtit="x (kpc)"
    ;ytit="y (kpc)"

    xmin= -xlen+center[0]
    xmax=  xlen+center[0]
    ymin= -xlen+center[1]
    ymax=  xlen+center[1]

    ; default is xy axis
    ; ------------------
    Axis1=0L    ; select x-axis
    Axis2=1L    ; select y-axis

    if keyword_set(xz) then begin
        xtit="x (kpc)"
        ytit="z (kpc)"
	ymin= -xlen+center[2]
	ymax=  xlen+center[2]
	Axis2=2L
    endif

    if keyword_set(yz) then begin
        xtit="y (kpc)"
        ytit="z (kpc)"
	xmin= -xlen+center[1]
	xmax=  xlen+center[1]
	ymin= -xlen+center[2]
	ymax=  xlen+center[2]
	Axis1=1L
	Axis2=2L
    endif




; -------------------------------------
;  Now What do we do with the image?
; -------------------------------------



        ; get position information
        ; --------------------------
;       plot, [-10000],[0], col=col, row=row, position=cpos, /noplot


	; -----------------------------
        ;   Send it to postscript
        ; -----------------------------

        initialize_plotinfo, 1

	setup_plot_stuff, 'ps', filename=filename, colortable= colortable
	;setup_plot_stuff, 'ps', filename=filename, colortable= 4   ; colors
	;setup_plot_stuff, 'ps', filename=filename, colortable= 3   ; red
	;setup_plot_stuff, 'ps', filename=filename, colortable= 1   ; blue
	;setup_plot_stuff, 'ps', filename=filename, colortable= 0   ; b/w
        ;setup_plot_stuff, 'ps', filename=filename, colortable= 4, imgxsize=14
	;setup_plot_stuff, 'ps', filename=filename, colortable= 0, imgxsize=14


	if keyword_set(nolabels) then begin
            x0= 0.01
            x1= 0.99
            y0= 0.01
            y1= 0.99
	endif else begin
	    x0= 0.20
            if keyword_set(colorbar) then x1= 0.85 else x1= 0.97
            y0= 0.15
            y1= 0.97
	endelse


        !p.position=[x0,y0,x1,y1]
        !p.ticklen=0.03

	;if keyword_set(colorbar) then pictxsize= 14 else pictxsize= 13


	;  place the image down
	; ----------------------
        if not keyword_set(particlesonly) then begin
		tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal
	endif
        ;print, "xxxxxxxxxxxxxxxxxxxxxxxxxx"
        ;print, "xxxxxxxxxxxxxxxxxxxxxxxxxx"
        ;print, "    WARNING !!!!!!  "
        ;print, "xxxxxxxxxxxxxxxxxxxxxxxxxx"
        ;print, "xxxxxxxxxxxxxxxxxxxxxxxxxx"
        ;print, "xxxxxxxxxxxxxxxxxxxxxxxxxx"
        ;print, "   pixel image has been"
        ;print, "     diabled."
        ;print, "xxxxxxxxxxxxxxxxxxxxxxxxxx"
        ;print, "xxxxxxxxxxxxxxxxxxxxxxxxxx"
        ;print, "xxxxxxxxxxxxxxxxxxxxxxxxxx"


        ; creates axes and plot style
        ; ---------------------------
	if keyword_set(nolabels) then begin
		plot,[0],[0], psym=3, xstyle=1, ystyle=1, color= 0, /noerase, $
                      xrange=[xmin,xmax], yrange=[ymin,ymax], $
                      xcharsize=0.01, ycharsize=0.01, $
                      xthick=4.0, ythick=4.0, /normal, /nodata, $
                      xtickformat='(a1)', ytickformat='(a1)'
	endif else begin
        	plot,[0],[0], psym=3, xstyle=1, ystyle=1, color= 0, /noerase, $
        	      xrange=[xmin,xmax], yrange=[ymin,ymax], $
		      xcharsize=1.50, ycharsize=1.50, $
		      xthick=4.0, ythick=4.0, /normal, /nodata, $
		      charthick=3.0, $
        	      xtitle=xtit, ytitle=ytit
	endelse


	; generate white axes if desired
	;axis, xaxis=0, xstyle=1,ystyle=1,color=1,charsize=0.001, xthick=2.0,ythick=2.0
        ;axis, xaxis=1, xstyle=1,ystyle=1,color=1,charsize=0.001, xthick=2.0,ythick=2.0
        ;axis, yaxis=0, xstyle=1,ystyle=1,color=1,charsize=0.001, xthick=2.0,ythick=2.0
        ;axis, yaxis=1, xstyle=1,ystyle=1,color=1,charsize=0.001, xthick=2.0,ythick=2.0



	; -----------------------------------
	;  Puts std information in upper left
	; -----------------------------------
	if not keyword_set(nolabels) then begin
	   if keyword_set(pubstyle) then begin
		xyouts, 0.23, 0.88, fload_timelbl(1,2), /normal, size= 1.7, charthick=3.0, color= 0
	   endif else begin
		; plot other information
		; ------------------------
               	xyouts, 0.22, 0.88, fload_fid(1), /normal, size= 1.7, charthick=3.0, color= 0
               	xyouts, 0.22, 0.83, fload_timelbl(1,3), /normal, size= 1.7, charthick=3.0, color= 0

		if keyword_set(msg) then begin
		   xyouts, 0.65, 0.88, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
		endif
	   endelse
	endif else begin
	     if keyword_set(pubstyle) then begin
		;xyouts, 0.07, 0.90, fload_timelbl(1,2,/noteq), /normal, size= 1.7, charthick=3.0, color= 0
		;xyouts, 0.07, 0.10, fload_timelbl(1,2,/noteq), /normal, size= 1.7, charthick=3.0, color= 0
		;xyouts, 0.07, 0.10, fload_timelbl(0.7,2,/noteq), /normal, size= 1.7, charthick=3.0, color= 0
		xyouts, 0.07, 0.10, fload_timelbl(0.7,2,/noteq)+' Gyr', /normal, size= 1.7, charthick=3.0, color= 0

		;;xyouts, 0.07, 0.83, 'side view', /normal, size= 1.7, charthick=3.0, color= 0
		;xyouts, 0.07, 0.83, fload_fid('firsttwo'), /normal, size= 1.7, charthick=3.0, color= 0
                xyouts, 0.07, 0.92, fload_fid(1), /normal, size= 1.7, charthick=3.0, color= 0

		;xyouts, 0.07, 0.83, '!6no black hole', /normal, size= 1.7, charthick=3.0, color= 0

		;sidelen= 2.0 * xlen
		sidelen= 2.0 * xlen / 0.7
	        xlenlbl= strcompress(string(sidelen),/remove_all)
	        if (sidelen) ge 1.0 then digs= 1
	        if (sidelen) ge 10.0 then digs= 2
	        if (sidelen) ge 100.0 then digs= 3
	        xlenlbl = strmid(xlenlbl,0,digs)        ; T=0.x (4+digits after decimal)
	        ;xlenlbl = xlenlbl+' kpc/h '
	        xlenlbl = xlenlbl+' kpc '
	        xlenlbl= '!94!6'+xlenlbl+'!96!6'
	        xyouts, 0.07, 0.05, xlenlbl, /normal, size= 1.4, charthick=3.0, color= 0



	     endif else begin
	        ; plot other information
                ; ------------------------
                xyouts, 0.07, 0.88, fload_fid(1), /normal, size= 1.7, charthick=3.0, color= 0
                xyouts, 0.07, 0.82, fload_timelbl(1,4), /normal, size= 1.7, charthick=3.0, color= 0
	     endelse

             if keyword_set(msg) then begin
                ;xyouts, 0.07, 0.91, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
		;xyouts, 0.07, 0.83, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
             endif
	endelse




		; -------------------------------
		;
                ;   Plot the actual points!!!
		;
                ; -------------------------------
                if keyword_set(plotpts) or keyword_set(particlesonly) then begin
                  if keyword_set(yz) then begin
                        x= y
                        y= z
                  endif

                  if keyword_set(xz) then y= z

                  usecen= fload_center_alreadycomp(1)


                  xtoplot= x-usecen[0]
                  ytoplot= y-usecen[1]

                  rad= sqrt(xtoplot*xtoplot + ytoplot*ytoplot)
                  ;idx= where(rad gt 20.0)
                  ;x= x(idx)
                  ;y= y(idx)

                  ;oplot, x, y, psym=5, thick=2.0, color= 50
                  ;oplot, x, y, psym=2, thick=2.0, color= 150
                  ;oplot, x, y, psym=7, thick=2.0, color= 100
                  ;oplot, x, y, psym=1, thick=2.0, color= 0
                  ;oplot, x, y, psym=3, color= 0
                  oplot, x, y, psym=3, color= 50

                endif



; ====================
;
;   Extras
;
; ====================

		; --------------------------
		;   Draw track
		; --------------------------
                if keyword_set(track_to_draw) then begin

			frun= fload_frun(1)


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



		; ---------------------
		;  Draw Circle
		; ---------------------
		;draw_re_circle=1
		draw_re_circle=0
		if draw_re_circle eq 1 then begin
			R_e= 2.4384

			phi = dindgen(101)*2D*!dpi/100
			; data coord
			re_x = R_e * cos(phi)
			re_y = R_e * sin(phi)

			re_x= re_x + center[0]
			re_y= re_y + center[1]

			loadct, 4
			tvlct,r,g,b,/get
			v1=[0,255]
			v2=[0,255]
			v3=[0,255]
			tvlct,v1,v2,v3,0  ; define colors 0 and 1 as black and white

			oplot, re_x, re_y, color= 220, thick=4.0, linestyle=0, psym=-3
			xyouts, 0.22, 0.68, 'R!De!N=2.4', /normal, size=1.3, color= 220, charthick= 4.0
			
		endif


                ; -----------------------------
                ;  Draw center path, similar to
		;  tracktodraw
                ; -----------------------------
		;drawcenters= 1
		drawcenters= 0
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



		; -----------------------------
		;  Black Hole X's
		;
		; Mark the position of the black
		; hole with a red X.
		; -----------------------------
		if keyword_set(showbhs) and fload_npart(5) ge 1 then begin

			fload_newcolortable, 4

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

			fload_newcolortable, colortable

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





                ; -----------------------
                ;
                ;   Add Young Stars
                ; --------------------
                ; 
                ;
                ; -----------------------
                overplot_youngstars= 0
                ;overplot_youngstars= 1
                if overplot_youngstars eq 1 then begin 
			fload_newcolortable, 1
                        contour_add_youngstars, 'junk', center=orig_center, $
                                rotate_theta=rotate_theta, rotate_phi=rotate_phi
                endif






		; -----------------------
		;
		;   ID's
		; --------------
		; 
		;
		; -----------------------
		overplot_idlist= 0
		;overplot_idlist= 1
		if overplot_idlist eq 1 then begin
			;contour_add_idlist, 'dog', center=orig_center, $
                        ;        rotate_theta=rotate_theta, rotate_phi=rotate_phi
			frun= fload_frun(1)
			center= orig_center
			contour_add_idlist, frun+'/idlist_rem_rgt4.txt', center=center, $
			;contour_add_idlist, frun+'/tidal_tail_idlist.txt', center=center, $
			;contour_add_idlist, frun+'/idlist_020_egt0.txt', center=center, $
			;contour_add_idlist, frun+'/idlist_020_e1.txt', center=center, $
			;contour_add_idlist, frun+'/idlist_020_e2.txt', center=center, $
			;contour_add_idlist, frun+'/idlist_020_e3.txt', center=center, $
                                rotate_theta=rotate_theta, rotate_phi=rotate_phi

			fload_newcolortable, 0
		endif




		; --------------------------------------------------------------------



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
                  bar_x0= 0.62
                  bar_y0= 0.12
                  bar_x1= 0.95
                  bar_y1= 0.17

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
                  axis, xaxis=1, xrange=[alog10(mi)+4.0,alog10(ma)+4.0], xstyle=1, color= 0, $
                        xcharsize=1.50, ycharsize=1.50, charthick=3.0, xthick= 3.0, /normal, $
                        xticks=1, xtickformat='(a1)', $
                        xtitle= 'Log !4R!D!3gas!N(M!D!9n!3!N pc!E-2!N)'   ; +4.0 above comes from converting to these units
                  ; puts ticks on bottom
                  axis, xaxis=0, xrange=[alog10(mi)+4.0,alog10(ma)+4.0], xstyle=1, $
                        xticks=3,xcharsize=1.20, ycharsize=1.20, charthick=3.0, xthick= 3.0, $
                        /normal, color= 0, xticklen=0.02
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
        ;docontour= 1
        docontour= 0
        if docontour EQ 1 then begin

		center= orig_center

                oplot_contours_comp, 'gas', x0, y0, x1, y1, $
			center=center, $
                        xlen=xlen, pixels=pixels

		; default
		;x= fload_gas_xyz('x',center=[0,0,0])
		;y= fload_gas_xyz('y',center=[0,0,0])
		;z= fload_gas_xyz('z',center=[0,0,0])
		;m= fload_gas_mass(1)
		;hsml= fload_gas_hsml(1)

		;center= orig_center

		;contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
                ;        hsml=hsml, $
                ;        pixels=pixels, $
                ;        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                ;        crude=crude, center=center, $
                ;        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                ;        NxNImage=NxNImage, HalfMassSB=HalfMassSB

		;GasPic= NxNImage


                ;oplot_contours, GasPic, x0, y0, x1, y1, $
                ;        xlen=xlen, pixels=pixels


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



	; ----------------------------
	;
	; Add angular momentum vector.
	;
        ; ----------------------------
        add_j_vector= 0
        if add_j_vector eq 1 then begin

                j_cold= fload_gas_cold_j(1)
                ;j_cold= j_cold-orig_center
		; rescale it
                j_cold= 20.0*j_cold/sqrt(total(j_cold*j_cold))

        	; transform to radians
        	;theta= !PI*rotate_theta/180.0
        	;phi= !PI*rotate_phi/180.0

        	;x= j_cold[0]
        	;y= j_cold[1]
        	;z= j_cold[2]
        	; first around y axis (theta)
        	; then around x axis (phi)
        	;x_new= x*cos(theta) + z*sin(theta)
        	;y_new= -x*(sin(phi)*sin(theta)) + y*cos(phi) + z*(cos(theta)*sin(phi))
        	;z_new= -x*(cos(phi)*sin(theta)) - y*sin(phi) + z*(cos(theta)*cos(phi))
        	;j_cold[0]= x_new
        	;j_cold[1]= y_new
        	;j_cold[2]= z_new

                ;arrow, 0.0, 0.0, j_cold[0], j_cold[1], /data, color= 0, thick=4.0
                arrow, 0.0, 0.0, j_cold[0], j_cold[2], /data, color= 0, thick=4.0
        endif




	; ----------------------------
	;
	;   Put Velocity Vectors On
	;
	; ----------------------------
	add_vel_vectors= 0
	;add_vel_vectors= 1
	if add_vel_vectors eq 1 then begin

		; velocity field variables
		min_velocity= 0.0
		;draw_a_key= 1
		draw_a_key= 0
        	GRIDPTS= 20L
		zthickness= 3.0
        	normlen= 2.0
        	VNorm= 200.0  ; normalize to 200  (then divide by 10, so 1 unit is 10 km/s on graph


		;vx= fload_gas_v('x')
		;vy= fload_gas_v('y')
		;vz= fload_gas_v('z')
		vx= fload_allstars_v('x')
		vy= fload_allstars_v('y')
		vz= fload_allstars_v('z')

		if keyword_set(xz) then begin
			vy=vz
			z=-y
		endif

		contour_velmap, x, y, z, vx, vy, vz, m, xlen, $
                        ;hsml=hsml, $
                        pixels=GRIDPTS, zthickness=zthickness, $
                        ;pixels=pixels, zthickness=zthickness, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        center=center, $
                        VelxMap= VelxMap, $
                        VelyMap= VelyMap



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

	endif




	; --------------------
	; done, close this up
	; --------------------
	device, /close




end


