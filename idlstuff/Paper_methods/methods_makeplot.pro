pro methods_makeplot, x, y, z, m, xlen, sendto, xz=xz, yz=yz, $
			filename=filename, thumbnail=thumbnail, fitstoo=fitstoo, $
			xthickness=xthickness, ythickness=ythickness, colorbar=colorbar, $
			pixels=pixels, zthickness=zthickness, $
			crude=crude, center=center, msg=msg, plotpts=plotpts, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			pubstyle=pubstyle, particlesonly=particlesonly, $
			track_to_draw=track_to_draw, nolabels=nolabels


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
if not keyword_set(sendto) then sendto='x'
if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(x) then begin
   print, "  "
   print, "PROBLEM: contour_makeplot, x, y, z, m, xlen, sendto, /xz, /yz,"
   print, "              filename=filename, /thumbnail, /fitstoo, /colorbar, "
   print, "              pixels=pixels, zthickness=zthickness,"
   print, "              xthickness=xthickness, ythickness=ythickness"
   print, "  "
   print, "  "
   print, "  "
   return
endif




don0toolow= 1




;--------------------------------------
;  Make image
;--------------------------------------

if not keyword_set(particlesonly) then begin
	contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
                        fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, colorbar=colorbar, $
                        pixels=pixels, zthickness=zthickness, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        crude=crude, center=center, $
			NxNImage=NxNImage

	Pic= NxNImage
endif



    ; default is xy axis
    ; ------------------
    ;xtit="!8x!3 (kpc/h)"
    ;ytit="!8y!3 (kpc/h)"
    xtit="x (kpc)"
    ytit="y (kpc)"

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
        if sendto EQ 'ps' and not keyword_set(nolabels) then begin


            initialize_plotinfo, 1

            setup_plot_stuff, sendto, filename=filename, colortable= 4


	    x0= 0.18
            if keyword_set(colorbar) then x1= 0.85 else x1= 0.95

            y0= 0.15
            y1= 0.95

            !p.position=[x0,y0,x1,y1]
            !p.ticklen=0.03

	    ;if keyword_set(colorbar) then pictxsize= 14 else pictxsize= 13


	    ;  place the image down
	    ; ----------------------
            if not keyword_set(particlesonly) then begin
		tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal
	    endif

                ; creates axes and plot style
                ; ---------------------------
                plot,[0],[0], psym=3,  $
                      xrange=[xmin,xmax], $
                      yrange=[ymin,ymax], $
                      /NOERASE, $
                      color= 0, $
		      xcharsize=1.50, ycharsize=1.50, $
		      xthick=4.0, ythick=4.0, $
		      charthick=3.0, $
                      xstyle=1, ystyle=1, $
                      xtitle=xtit, $
                      ytitle=ytit, /normal, $
                      /nodata



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

                ; create color bar
                ; ----------------
                if keyword_set(colorbar) then begin
                  bar= REPLICATE(1B, 10)#BINDGEN(256)
		  if invertcolors eq 1 then bar= 255-bar
		  idx= where(bar eq 1)
		  if idx(0) ne -1 then begin
			print, "idx has ", n_elements(idx), idx(0)
			colornumber= idx(0)-1
			if colornumber lt 0 then colornumber= 11     ; if it's bar[*,0] then set to bar[*,1]
			bar(idx)= bar(colornumber)
		  endif
                  barwidth= 0.04

                  tv, bar, x1+0.01, y0, xsize=barwidth, ysize=(y1-y0), /normal

                  !p.position= [x1+0.01,y0,x1+0.01+barwidth,y1]

		  ma= 5.0e-2
		  mi= 5.0e-5
                  plot, [0],[0], psym=3, xrange=[0,1], yrange=[alog10(mi),alog10(ma)], $
                        /normal, $
                        /noerase, color= 1, xstyle=1, ystyle=1, $
			xthick=4.0, ythick=4.0, $
                        xticks=1, xtickformat='(a1)', ytickformat='(a1)', $
                        /nodata

                  axis, yaxis=1, yrange=[alog10(mi)+4.0,alog10(ma)+4.0], ystyle=1, $
			xcharsize=1.50, ycharsize=1.50, charthick=3.0, $
                        /normal, $
                        ytitle= 'Log !4R!D!3gas!N(M!D!9n!3!N/pc!E2!N)'   ; +4.0 above comes from converting to these units
                endif


		if keyword_set(pubstyle) then begin
			xyouts, 0.23, 0.88, fload_timelbl(1,2), /normal, size= 1.7, charthick=3.0, color= 0
		endif else begin
			; plot other information
			; ------------------------
                	xyouts, 0.22, 0.88, fload_fid(1), /normal, size= 1.7, charthick=3.0, color= 0
                	xyouts, 0.22, 0.83, fload_timelbl(1,2), /normal, size= 1.7, charthick=3.0, color= 0

			if keyword_set(msg) then begin
			   xyouts, 0.78, 0.88, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
			endif
		endelse

	        device, /close

        endif 





        ; ---------------------------------------------------------------------
        ;   Send it to postscript - and fill the plotted area with the figure
        ; ---------------------------------------------------------------------
        if sendto EQ 'ps' and keyword_set(nolabels) then begin


	    initialize_plotinfo, 1

	    if not keyword_set(don0toolow) then begin
		setup_plot_stuff, sendto, filename=filename, colortable= 4
		;setup_plot_stuff, sendto, filename=filename, colortable= 0
	    endif else begin
		setup_plot_stuff, sendto, filename=filename, colortable= 4, newysize=14.0
		;setup_plot_stuff, sendto, filename=filename, colortable= 0, newysize=14.0
	    endelse


            x0= 0.01
            x1= 0.99
            if keyword_set(don0toolow) then y0= 0.15 else y0= 0.01
            y1= 0.99
            !p.position=[x0,y0,x1,y1]


            ;  place the image down
            ; ----------------------
            tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

            ; creates axes and plot style
            ; ---------------------------
            plot,[0],[0], psym=3,  $
                      xrange=[xmin,xmax], $
                      yrange=[ymin,ymax], $
                      /NOERASE, $
                      color= 0, $
                      xcharsize=0.01, ycharsize=0.01, $
                      xthick=4.0, ythick=4.0, $
                      xstyle=1, ystyle=1, $
		      xtickformat='(a1)', ytickformat='(a1)', $
                      /normal, $
                      /nodata


	     if keyword_set(pubstyle) then begin
		;xyouts, 0.07, 0.90, fload_timelbl(1,2,/noteq), /normal, size= 1.7, charthick=3.0, color= 0
		xyouts, 0.07, 0.90, fload_timelbl(1,2), /normal, size= 1.7, charthick=3.0, color= 0
		;xyouts, 0.07, 0.83, 'side view', /normal, size= 1.7, charthick=3.0, color= 0
		;xyouts, 0.07, 0.83, fload_fid('firsttwo'), /normal, size= 1.7, charthick=3.0, color= 0
	     endif else begin
	        ; plot other information
                ; ------------------------
                xyouts, 0.07, 0.88, fload_fid(1), /normal, size= 1.7, charthick=3.0, color= 0
                xyouts, 0.77, 0.83, fload_timelbl(1,2), /normal, size= 1.7, charthick=3.0, color= 0
	     endelse

             if keyword_set(msg) then begin
                ;xyouts, 0.07, 0.91, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
		xyouts, 0.07, 0.83, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
             endif


                if keyword_set(track_to_draw) then begin
                        xtrack= track_to_draw[*,Axis1]
                        ytrack= track_to_draw[*,Axis2]
                        idx= where((xtrack ne 0) or (ytrack ne 0))
                        if idx(0) ne -1 then begin
                          xtrack= xtrack(idx)
                          ytrack= ytrack(idx)
                          oplot, xtrack, ytrack, color= 0, thick=3.0, linestyle= 2, psym=-3
                        endif
                endif

                ; create color bar
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
		  ;	xcharsize=1.50, ycharsize=1.50, charthick=3.0, $
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


                ; create color bar
                ; ----------------
		colorbar3= 1
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
		  ;	xcharsize=1.50, ycharsize=1.50, charthick=3.0, $
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




	     ; done, close this up
	     ; --------------------
	     device, /close


	endif










	; ------------------------------------
	;   Send it to the screen
	; ------------------------------------
	if sendto EQ 'x' then begin


	        set_plot, 'x'
        	thetitle = "density contour  ->  "+fload_fid(1)
        	!p.background= getcolor('white')
        	window,11,xsize=600,ysize=600, title=thetitle

		;tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

                ; creates axes and plot style
                ; ---------------------------
                plot,[0],[0], psym=3,  $
                      xrange=[xmin,xmax], $
                      yrange=[ymin,ymax], $
		      /noerase, $
                      color= 0, $
                      xstyle=1, ystyle=1, $
                      charsize=1.4, $
                      xthick=1.7, ythick=1.7, $
                      xtitle=xtit, $
                      ytitle=ytit, /normal, $
                      /nodata


		; plot the actual points
		; -------------------------
		if keyword_set(plotpts) then begin
	          if keyword_set(yz) then begin
	                x= y
	                y= z
	          endif

	          if keyword_set(xz) then y= z
		  oplot, x, y, psym=1, color= getcolor('red')
		endif

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




        ; time on graph
        ; --------------
	if sendto EQ 'x' then begin
		xyouts, 0.22, 0.88, fload_fid(1), /normal, size= 1.7, color= 1
		xyouts, 0.22, 0.83, fload_timelbl(1,2), /normal, size= 1.7, color= 1
	endif




	; -------------------------------------------
	;    image files
	; -------------------------------------------
	if (sendto EQ 'jpg') or (sendto EQ 'gif') or (sendto EQ 'png') then begin
           if not keyword_set(filename) then begin
                if sendto EQ 'jpg' then filename= 'contour.jpg' else filename= 'contour.gif'
                ans= ''
                read, ans, PROMPT='eps filename ['+filename+']:'
                if strlen(ans) GT 0 then filename= ans
           endif

           loadct, 4
           tvlct,r,g,b,/get
           v1=[0,255]
           v2=[0,255]
           v3=[0,255]
           tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)

           tvlct,r1,g1,b1,/get     ; reget new color table with the black and white

	endif



	; write a jpg file (this currently isn't working, something with the color table)
	; ----------------
	if sendto EQ 'jpg' then begin

		write_jpeg, filename, Pic, true=1

	endif


        ; write a gif file if that's what we desire
        ; -----------------------------------------
        if sendto EQ 'gif' then begin

		if keyword_set(thumbnail) then begin
			n= n_elements(Pic[0,*])
			Pic[n-2:n-1,*]= 0
			Pic[0:1,*]= 0
			Pic[*,n-2:n-1]= 0
			Pic[*,0:1]= 0

			Pic= congrid(Pic, 200, 200)
		endif

		bPic= byte(Pic)
                write_gif, filename, bPic, r1, g1, b1
        endif


	if sendto eq 'png' then begin
		write_png, filename, Pic, r1, g1, b1
	endif

	; ----------------------------------
	; use the z buffer
	; ----------------------------------
	; this is a nice way to combine just the pixel image
	; and words or these types of things
	;
	if sendto eq 'z' then begin
           if not keyword_set(filename) then begin
                filename= 'contour.png'
		;filename= 'contour.eps'
                ans= ''
                read, ans, PROMPT='eps filename ['+filename+']:'
                if strlen(ans) GT 0 then filename= ans
           endif

           set_plot, 'z'
           device, set_resolution= [500,500], set_character_size= [12,18]

           loadct, 4
           tvlct,r,g,b,/get
           v1=[0,255]
           v2=[0,255]
           v3=[0,255]
           tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)

           tvlct,r1,g1,b1,/get     ; reget new color table with the black and white


		tv, Pic, 0, 0, xsize=1.0, ysize=1.0, /normal
		xyouts, 0.05, 0.90, fload_timelblnot(1), /normal, color= 1, charthick=3.0, size=1.5

		; draw 10 kpc line
		barlength= 10.0/(xmax-xmin)    ; in normalized coordinates
		; can we do plot in normalized coordinates?

		plottedimage= tvrd()

		; send to ps file - strange white on the edge
		;set_plot, 'ps'
        	;device, filename= filename,/color,bits_per_pixel=8
		;device, xsize=10, ysize= 10
		;tv, plottedimage
		;device, /close


		; send it to png file
		;plottedimage= reverse(plottedimage, 1)
		write_png, filename, plottedimage, r1, g1, b1
	endif



; -------------
;  Done
; -------------



end


