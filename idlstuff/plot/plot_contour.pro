;-------------------------------------------------------
;-------------------------------------------------------
;
;  general plotting procedure
;
;
;-------------------------------------------------------
;-------------------------------------------------------
pro plot_contour, runinfo, plottype, panels=panels, $
			filename=filename, $
			arepo=arepo, $
			evensampling=evensampling, $
                        center=center, $
                        colorbar=colorbar, $
                        colortable=colortable, $
                        crude=crude, $
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
			threeprojection=threeprojection, $
                        thumbnail=thumbnail, $
                        track_to_draw= track_to_draw, $
                        use_calc_center=use_calc_center, $
                        use_bary_center=use_bary_center, $
                        xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
                        xz=xz, yz=yz, $
			xxlen=xxlen


if not keyword_set(runinfo) then begin
   print, "  "
   print, "plot_contour, runinfo, plottype, ..... (see below for options) panels=panels, "
   print, "    /evensampling, /arepo, "
   print, "    xxlen=xxlen "
   print, "  "
   print, " keyword [default]       - description "
   print, " crude [off]             - make a crude pixel map rather than c-program one"
   print, " center [0,0,0]          - zoom in on center"
   print, " colorbar [off]          - prints color bar and scale on right side"
   print, " contour [off]           - Can do, but need to turn it on within"
   print, " evensampling [off]      - will automatically grab regularly spaced images and do a 3x4 plot"
   print, " filename [contour.eps]  - same file as"
   print, " fitstoo [off]           - write a fits file"
   print, " loadedsnap [off]        - snap is already load it, routine does not"
   print, " msg [gas]               - print a message on figure"
   print, " nolabels [off]          - suppress labels, so that the figure fills the plotted area"
   print, " numpart [all]           - with startid, can select to print certain particles"
   print, " panels ['3x4']          - make a 3 x 4 panel plot, requires single quotes, can't do >10 dims"
   print, " particlesonly [off]     - plot the particles rather than the density map"
   print, " plotpts [off]           - plot the actual data points, can be prohibitive if there are many"
   print, " pixels [480]            - "
   print, " pubstyle [off]          - omit run name and msg printing"
   print, " rotate_phi [0]          - rotate the image by phi, about z-axis"
   print, " rotate_theta [0]        - rotate the image by theta, about y-axis"
   print, " showbhs [off]           - put an x where the black hole is"
   print, " snapnum [100]           - snapshot number, reads snapshot_<snapnum>"
   print, " startid [all]           - with numpart, can select to print certain particles"
   print, " track_to_draw [off]     - will add the path as found in centers.txt"
   print, " use_calc_center [off]   - uses the center automatically calculated when a snap is opened"
   print, " use_bary_center [off]   - uses the baryonic center "
   print, " xxlen [100.0]            - length of each side is 2x(xxlen)"
   print, " x|y|zthickness [off]    - takes slice of this thickness centered at center"
   print, " xz [off]                - project along xz, default is xy"
   print, " yz [off]                - project along yz, default is xy"
   print, "  "
   print, "  "
   print, "  "
   return
endif


; do some checking, if runinfo is only a string, then
; create a structure with it and print a single plot
;

n_images= (size(runinfo))[1]    ; not sure if this accomplishes the above

s= size(runinfo)
if s[s[0]+1] ne 8 then begin
	print, " "
	print, " NOTICE: runinfo is not a structure, creating one"
	print, " "
endif



;
;
;

if not keyword_set(filename) then filename='panelfig.eps'

if not keyword_set(panels) then panels='none'

if not keyword_set(colortable) then colortable= 4


center=[0.0,0.0,0.0]


; default is stars
if not keyword_set(plottype) then plottype= "allstars"


; universal xlen
if keyword_set(xxlen) then runinfo[*].xlen= xxlen

;
; 3x4 panel which shows regular frequency of snapshots
;
if keyword_set(evensampling) then begin

	; determine the number of
	; snapshots in frun directory
	spawn, "/bin/ls "+frun+"/snapshot* | wc ",result
	nsnaps=long(result[0])

	snaparr= (indgen(12) + 1) * long(nsnaps/12)
	idx=where(snaparr ge nsnaps)
	if idx(0) ne -1 then snaparr(idx)= nsnaps-1
	snaparr[11]= nsnaps-1

	; now put this in the structure
	if n_images ne 1 then begin
		print, " "
		print, " I'm confused the the dimensionality of the structure, I'm stopping. "
		print, " "
		return
	endif

	panels='3x4'

endif



; 
; 3x1 panel which shows all three projections of a snapshot
;
if keyword_set(threeprojection) then begin

        ; now put this in the structure
        if n_images ne 1 and n_images ne 3 then begin
                print, " "
                print, " n_images ne (1 or 3) - "
                print, " I'm confused about the dimensionality of the structure, I'm stopping. "
                print, " "
                return
        endif

	;
	;  make the appropriate 3 element galaxy information structure
	;     (might crap out if rotate_theta is already in there)
	;
	if where((tag_names(runinfo)) eq "ROTATE_THETA") lt 0 then begin
		tempinfo= replicate(create_struct(runinfo[0], 'rotate_theta', ''), 3)
	endif else begin
		tempinfo= replicate(runinfo[0], 3)
	endelse

	tempinfo[*].frun= runinfo[0].frun
	tempinfo[*].snapnum= runinfo[0].snapnum
	tempinfo[*].xlen= runinfo[0].xlen
	tempinfo[0].rotate_theta= "xy"
	tempinfo[1].rotate_theta= "xz"
	tempinfo[2].rotate_theta= "yz"

	runinfo= tempinfo

        panels='3x1'
	n_images= (size(runinfo))[1]    ; not sure if this accomplishes the above
endif




;===========================================================================



;--------------------------------------
;  Set plot dimensions
;--------------------------------------


;default is single image

xdim= 1
ydim= 1

xtit="!8x!3 (kpc/h)"
ytit="!8y!3 (kpc/h)"
;xtit="x (kpc)"
;ytit="y (kpc)"

x0= 0.18
xs= 0.80
y0= 0.15
ys= 0.83

if keyword_set(colorbar) then xs= 0.68

if keyword_set(nolabels) then begin
	x0= 0.01
	xs= 0.98
	y0= 0.01
	ys= 0.98
endif



if strpos(panels,'x') gt 0 then begin

	; assume it's always xdir x ydir and both are single digits
	xdim= long(strmid(panels,0,1))
	ydim= long(strmid(panels,2,1))
	ndim= xdim * ydim

	if n_images ne ndim then begin
		print, " "
		print, " Problem, stopping: "
		print, "      panels= ", panels
		print, "      ndim= ", ndim
		print, "      n_images= ", n_images
		print, " "
		return
	endif

	;newxsize= 6.0 * xdim
	;newysize= 6.0 * ydim
	newxsize= 10.0 * (xdim^(1.5*abs(xdim-ydim)/(xdim+ydim)))
	newysize= 10.0 * (ydim^(1.5*abs(xdim-ydim)/(xdim+ydim)))
	;newxsize= 14.0 * (xdim^(3./10.))
	;newysize= 14.0 * (ydim^(3./10.))

        x0= 0.01
        xs= (1.0 - 2.0*x0) / xdim

        y0= 0.01
        ys= (1.0 - 2.0*y0) / ydim

endif


;--------------------------------------
;  Now plot this mess (to postscript)
;--------------------------------------


        initialize_plotinfo, 1

        setup_plot_stuff, 'ps', filename=filename, colortable= colortable, newxsize=newxsize, newysize=newysize


	paneli= 0

	for yy= 0, ydim-1 do begin
	  for xx= 0, xdim-1 do begin

        	print, " "
        	print, "paneli, xx, yy= ", paneli, xx, yy
        	print, " "

		xo= x0 + ((paneli mod xdim ) * xs)
		yo= y0 + ((ydim-1-yy ) * ys)

        	print, "xo, xs= ", xo, xs
        	print, "yo, ys= ", yo, ys

		frun= runinfo[paneli].frun
		snapnum= runinfo[paneli].snapnum
		xlen= runinfo[paneli].xlen

        	print, "frun= ", frun
        	print, "snapnum= ", snapnum
		print, " "

		if where((tag_names(runinfo)) eq "ROTATE_PHI") ge 0 then rotate_phi= runinfo[paneli].rotate_phi
		if where((tag_names(runinfo)) eq "ROTATE_THETA") ge 0 then rotate_theta= runinfo[paneli].rotate_theta
			
		if keyword_set(xxlen) then xlen=xxlen

		img_do_one_panel, frun, snapnum, xlen, $
				xo, yo, xs, ys, $
        			plottype=plottype, $
        			draw_circle=draw_circle, $
        			origcolortable=origcolortable, $
        			arepo=arepo, $
				fitstoo=fitstoo, $
                        	xthickness=xthickness, ythickness=ythickness, colorbar=colorbar, $
                        	pixels=pixels, zthickness=zthickness, $
                        	crude=crude, center=center, msg=msg, plotpts=plotpts, $
                        	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        	nolabels=nolabels, pubstyle=pubstyle, particlesonly=particlesonly, $
                        	showbhs=showbhs, set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        	track_to_draw=track_to_draw, colortable=colortable


		paneli= paneli + 1

	  endfor
	endfor




	; print msg?
	; -----------
	if not keyword_set(msg) then msg= fload_fid(1)
	xyouts, 0.03, y0+0.055, msg, /normal, size= 1.2, charthick=3.0, color= 0   ; 0=black, 1=white
	;xyouts, 0.07, 0.83, msg, /normal, size= 1.2, charthick=3.0, color= 0   ; 0=black, 1=white


	; put time units on panel 1
	; has to match up with the xyouts in plot/oplot_one_axes.pro
	; --------------------------
       	xyouts, x0+0.13, y0+ys-0.06, '!6Gyr', /normal, size= 1.2, charthick=3.0, color= 0
       	;xyouts, x0+0.065, y0+ys-0.065, '!6Gyr/h', /normal, size= 1.2, charthick=3.0, color= 0


	; put side length on panel 1
	; (by default in lower-left of panel 0,0)
	; -----------------------------------------
	xlenarr=[xlen,xlen]
        xlenlbl= strcompress(string(2.0*xlenarr[0]),/remove_all)
        if (2.0*xlenarr[0]) ge 1.0 then digs= 1
        if (2.0*xlenarr[0]) ge 10.0 then digs= 2
        if (2.0*xlenarr[0]) ge 100.0 then digs= 3
        xlenlbl = strmid(xlenlbl,0,digs)        ; T=0.x (4+digits after decimal)
        xlenlbl= '!94!6'+xlenlbl+'!96!6'
        ;xyouts, x0+xs-0.06, y0+0.05, xlenlbl, /normal, size= 1.2, charthick=3.0, color= 0
        xyouts, x0+0.04, y0+0.05, xlenlbl, /normal, size= 1.0, charthick=2.0, color= 0




        ; done, close this up
        ; --------------------
	print, " "
	print, " plot saved to: ", filename
	print, " "
	device, /close






; -------------
;  Done
; -------------



end






