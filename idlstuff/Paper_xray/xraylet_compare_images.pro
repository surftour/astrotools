pro produce_one_image, frun, snapnum, xlen, sendto, $
			center=center, $
			colorbar=colorbar, $
			crude=crude, $
			filename=filename, $
			fitstoo=fitstoo, $
			loadedsnap=loadedsnap, $
			loadstars=loadstars, $
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
      ;if (fload_snapshot(frun, snapnum)) then begin
      if (fload_snapshot_bh(frun, snapnum)) then begin
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
    x= fload_gas_xyz('x',center=[0,0,0])
    y= fload_gas_xyz('y',center=[0,0,0])
    z= fload_gas_xyz('z',center=[0,0,0])
    m= fload_gas_mass(1)


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

if keyword_set(loadstars) then begin
    ndisk= fload_npart(2)                               ; disk particles
    nbulge= fload_npart(3)                              ; bulge particles
    nstars= fload_npart(4)
    npart= long(ndisk) + long(nbulge) + long(nstars)
    N= long(npart)
    print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars
    
    x= fload_allstars_xyz('x',center=[0,0,0])
    y= fload_allstars_xyz('y',center=[0,0,0])
    z= fload_allstars_xyz('z',center=[0,0,0])
    m= fload_allstars_mass(1)
endif




;--------------------------------
;   New Stellar Particles
;--------------------------------

    ;x= fload_newstars_xyz('x',center=[0,0,0])
    ;y= fload_newstars_xyz('y',center=[0,0,0])
    ;z= fload_newstars_xyz('z',center=[0,0,0])
    ;m= fload_newstars_mass(1)



; ----------------------------------------------
;  Call our generic contour plotting procedure
; ----------------------------------------------

set_maxden= 1.0    ; 10^4 msolar/pc2
;set_dynrng= 1.0e+6
set_dynrng= 1.0e+5
;set_maxden= 0.5e-1    ; 10^4 msolar/pc2
;set_dynrng= 1.0e+4

if keyword_set(loadstars) then set_dynrng= 1.0e+4


contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
	filename=filename, fitstoo=fitstoo, $
	xthickness=xthickness, ythickness=ythickness, $
	pixels=pixels, zthickness=zthickness, $
	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	crude=crude, center=center, $
	set_maxden= set_maxden, $
	set_dynrng= set_dynrng, $
	NxNImage=NxNImage

	ReturnImage= NxNImage




	;ReturnTimeLbl= fload_timelbl(1,2,/noteq)
	ReturnTimeLbl= fload_timelbl(0.7,2,/noteq)
; -------------
;  Done
; -------------



end









pro produce_xray_image, frun, snapnum, xlen, sendto, $
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
      ;if (fload_snapshot(frun, snapnum)) then begin
      if (fload_snapshot_bh(frun, snapnum)) then begin
        print, "PROBLEM: opening file"
        return
      endif
    endif


;---------------------------------
;  Gas Particles
;---------------------------------

    ; default
    x= fload_gas_xyz('x',center=[0,0,0])
    y= fload_gas_xyz('y',center=[0,0,0])
    z= fload_gas_xyz('z',center=[0,0,0])
    ;m= fload_gas_mass(1)
    xray= fload_gas_xray_luminosity(1,/diffuse_hotgas)   ; returns in ergs sec^-1
print, "L_x= ",total(xray)
    idxx= where(xray gt 0.0)
    xray= xray(idxx)/1.0e+30
    x= x(idxx)
    y= y(idxx)
    z= z(idxx)

    ;set_maxden= 1.0e+10     ; 1e40
    set_maxden= 1.0e+9     ; 1e39
    ;set_maxden= 1.0e+12     ; 1e42
    set_dynrng= 1.0e+6      ; 1e34
    ;set_dynrng= 1.0e+6      ; 1e37
    ;set_dynrng= 1.0e+8      ; 1e34
    ;set_dynrng= 1.0e+3      ; 1e34

    ; this is our special one that totals
    contour_makepic_total, x, y, z, xray, xlen, xz=xz, yz=yz, $
                        filename=filename, fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, colorbar=colorbar, $
                        pixels=pixels, zthickness=zthickness, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        center=center, $  ; /crude,  $
                        NxNImage=NxNImage, $
                        set_maxden= set_maxden, $
                        set_dynrng= set_dynrng      ;, /var_already_inlog


	ReturnImage= NxNImage


; -------------
;  Done
; -------------



end










;========================================
;
;   Primary callable program
;
;========================================
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
xlen= 35.0
filename='xraypanel.eps'
sendto= 'ps'

;xz= 1

frun1= "/raid4/tcox/vc3bvc3b"
frun2= "/raid4/tcox/vc3bvc3b_no"

snapnum= 20


;--------------------------------------
;  Galaxy Number One
;--------------------------------------

; gas image
; -----------
produce_one_image, frun1, snapnum, xlen, xz=xz, yz=yz, $
	xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
	pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	crude=crude, center=center, ReturnImage=ReturnImage, ReturnTimeLbl= ReturnTimeLbl

GasPic1= ReturnImage
TimeLbl1= ReturnTimeLbl

; xray image
; -----------
produce_xray_image, frun1, snapnum, xlen, xz=xz, yz=yz, /loadedsnap, $
	xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
	pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	crude=crude, center=center, ReturnImage=ReturnImage

XrayPic1= ReturnImage

; stellar image
; -------------
produce_one_image, frun1, snapnum, xlen, xz=xz, yz=yz, /loadedsnap, /loadstars, $
	xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
	pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	crude=crude, center=center, ReturnImage=ReturnImage, ReturnTimeLbl= ReturnTimeLbl

StarsPic1= ReturnImage




;--------------------------------------
;  Galaxy Number Two
;--------------------------------------

center=[4.0,0.0,4.0]

; gas image
; -----------
produce_one_image, frun2, snapnum, xlen, xz=xz, yz=yz, $
	xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
	pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	crude=crude, center=center, ReturnImage=ReturnImage, ReturnTimeLbl= ReturnTimeLbl

GasPic2= ReturnImage

; xray image
; -----------
produce_xray_image, frun2, snapnum, xlen, xz=xz, yz=yz, /loadedsnap, $
	xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
	pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	crude=crude, center=center, ReturnImage=ReturnImage

XrayPic2= ReturnImage

; stellar image
; --------------
produce_one_image, frun2, snapnum, xlen, xz=xz, yz=yz, /loadedsnap, /loadstars, $
	xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
	pixels=pixels, rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	crude=crude, center=center, ReturnImage=ReturnImage, ReturnTimeLbl= ReturnTimeLbl

StarsPic2= ReturnImage





;--------------------------------------
;  Now plot this mess
;--------------------------------------

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

	; -----------------------------
        ;   Send it to postscript
        ; -----------------------------


        initialize_plotinfo, 1

	;setup_plot_stuff, sendto, filename=filename, colortable= 4
        ;setup_plot_stuff, sendto, filename=filename, colortable= 4, newxsize=18, newysize=24
        ;setup_plot_stuff, sendto, filename=filename, colortable= 1, newxsize=30, newysize=20
        setup_plot_stuff, sendto, filename=filename, colortable= 0, newxsize=21, newysize=14


	x0= 0.005
	xs= 0.33   ; assumes 3 panels
	x1= x0+xs
	x2= x0+xs+xs
	x3= x0+xs+xs+xs

        y0= 0.005
	ys= 0.495   ; assumes 2 panels
        y1= y0+ys
	y2= y0+ys+ys



	;  place the image(s) down
	; ----------------------

	; stars
	tv, StarsPic1, x0, y0, xsize=xs, ysize=ys, /normal
	tv, StarsPic2, x0, y1, xsize=xs, ysize=ys, /normal


	loadct, 1
	tvlct,r,g,b,/get
	v1=[0,255]
	v2=[0,255]
	v3=[0,255]
	tvlct,v1,v2,v3,0

	; gas
	tv, GasPic1, x1, y0, xsize=xs, ysize=ys, /normal
	tv, GasPic2, x1, y1, xsize=xs, ysize=ys, /normal

	loadct, 3
	tvlct,r,g,b,/get
	v1=[0,255]
	v2=[0,255]
	v3=[0,255]
	tvlct,v1,v2,v3,0

	; xrays
	tv, XrayPic1, x2, y0, xsize=xs, ysize=ys, /normal
	tv, XrayPic2, x2, y1, xsize=xs, ysize=ys, /normal



        ; creates axes and plot style
        ; ---------------------------

	;  first row
	; -----------
	overplot_axes, x0, y1, x1, y2, TimeLbl1, xlen
	overplot_axes, x1, y1, x2, y2, ' ', xlen
	overplot_axes, x2, y1, x3, y2, ' ', xlen

	;  second row
	; -------------
	overplot_axes, x0, y0, x1, y1, ' ', xlen
	overplot_axes, x1, y0, x2, y1, ' ', xlen
	overplot_axes, x2, y0, x3, y1, ' ', xlen




	; overplot xray contours
	print, "getting to overplot"
	overplot_xraycontour, XrayPic1, x2, y0, x3, y1, xlen
	overplot_xraycontour, XrayPic2, x2, y1, x3, y2, xlen



	; print msg?
	; -----------
        if keyword_set(msg) then begin
                ;xyouts, 0.07, 0.91, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
		xyouts, 0.07, 0.83, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
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

                ; create color bar
                ; ----------------
		invertcolors= 1    ; warning: this should match the value in contour_makepic
		MaxDensXY= 5.0e-2
		DynRange=1.0e3
		ma=MaxDensXY
		mi=MaxDensXY/DynRange
                if keyword_set(colorbar) then begin
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

			n_bh= fload_npart(5)
			;bhid= fload_blackhole_id(1)
			;bhid1= bhid[0]
			;bh_center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)

			; for the moment, just plot xy black hole position
                        ;oplot, [bh_center_1[0]], [bh_center_1[1]], psym=7, thick=6.0, color=150, symsize=2.0

			;if n_bh gt 1 then begin
			;	bhid2= bhid[1]
			;	bh_center_2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
			;	print, "Blackhole #1: ",bhid1, bh_center_1
			;	print, "Blackhole #2: ",bhid2, bh_center_2

			;	oplot, [bh_center_2[0]], [bh_center_2[1]], psym=7, thick=6.0, color=150, symsize=2.0
				
			;endif else begin
			;	print, "Blackhole: ", bhid1, bh_center_1
			;endelse

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
                  bar_y0= y3+0.032
                  bar_x1= x0+0.30
                  bar_y1= y3+0.047

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
                  plot, [0],[0], psym=3, yrange=[0,1], xrange=[alog10(mi)+4.0,alog10(ma)+4.0], $
                        /normal, /noerase, color= 0, xstyle=1, ystyle=1, $
                        xthick=3.0, ythick=3.0, yticks= 1, $
                        xticks=3, xtickformat='(a1)', ytickformat='(a1)', $
                        /nodata, xticklen=0.22
                
                  ; puts title and ticks on right-hand side
                  ;axis, yaxis=1, yrange=[alog10(mi)+4.0,alog10(ma)+4.0], ystyle=1, $                  
                  ;     xcharsize=1.50, ycharsize=1.50, charthick=3.0, $
                  ;      /normal, $
                  ;      ytitle= 'Log !4R!D!3gas!N(M!D!9n!3!N pc!E-2!N)'   ; +4.0 above comes from converting to these units

                  ; puts title on top
                  axis, xaxis=1, xrange=[alog10(mi)+4.0,alog10(ma)+4.0], xstyle=1, color= 0, $
                        xcharsize=0.95, charthick=3.0, xthick= 3.0, /normal, $
                        xticks=3, xtickformat='(a1)', xticklen=0.08, $
                        ;xtitle= 'Log !4R!D!3gas!N(M!D!9n!3!N pc!E-2!N)'
                        xtitle= 'Log !4R!D!3 !N(M!D!9n!3!N pc!E-2!N)'
			; +4.0 above comes from converting to these units
                  ; puts ticks on bottom
                  axis, xaxis=0, xrange=[alog10(mi)+4.0,alog10(ma)+4.0], xstyle=1, $
                        xticks=3,xcharsize=0.95, ycharsize=0.95, charthick=3.0, xthick= 3.0, $
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



        ; done, close this up
        ; --------------------
	if sendto eq 'ps' then device, /close






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
	   if (fload_snapshot(frun, snapnum)) then begin
	   ;if (fload_snapshot_bh(frun, snapnum)) then begin
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
                        xthickness=xthickness, ythickness=ythickness, colorbar=colorbar, $
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







;-----------------------------------------
; Overplot New Stars
;-----------------------------------------
pro overplot_xraycontour, imagemap, x0, y0, x1, y1, xlen


        ;loadct, 4               ; blue/green/red/yellow (std)
        loadct, 0               ; blue/green/red/yellow (std)
        tvlct,rr,gg,bb,/get
        vv1=[0,255]
        vv2=[0,255]
        vv3=[0,255]
        tvlct,vv1,vv2,vv3,0


	; plot axes
	; ----------
	xmax= xlen
	xmin=-xlen
	ymax= xlen
	ymin=-xlen
        !p.position=[x0,y0,x1,y1] & !p.ticklen=0.03
        plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.00, ycharsize=0.00, xstyle=1, ystyle=1, /normal, /nodata, $
		      xticks= 1, yticks= 1, $
                      xthick=0.0, ythick=0.0, xtickformat='(a1)', ytickformat='(a1)'


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
                contour, imagemap, path_xy=xy_contours, path_info=info, levels=userLevels, $
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
			;if idx(0) eq -1 then plots, xycnt, /normal, color= 220
                        if idx(0) eq -1 then begin

			   ; do the following if colortable=4
                           ;if n_level eq 0 then plots, xycnt, /normal, color= 150, thick=5.0
                           ;if n_level eq 1 then plots, xycnt, /normal, color= 150, thick=4.0
                           ;if n_level eq 2 then plots, xycnt, /normal, color= 150, thick=3.0
                           ;if n_level eq 3 then plots, xycnt, /normal, color= 150, thick=2.0

			   ; do the following if colortable=0
                           if n_level eq 0 then plots, xycnt, /normal, color= 0, thick=5.0
                           if n_level eq 1 then plots, xycnt, /normal, color= 0, thick=4.0
                           if n_level eq 2 then plots, xycnt, /normal, color= 0, thick=3.0
                           if n_level eq 3 then plots, xycnt, /normal, color= 0, thick=2.0

                           ;if n_level eq 3 then plots, xycnt, /normal, color= 150, thick=1.0, linestyle= 1
                           ;if n_level eq 4 then plots, xycnt, /normal, color= 130
                           ;plots, xy_contours(*,info(i).offset+cnt_index), /normal, color= 150
                        endif
                endfor


; -------------
;  Done
; -------------


end











pro overplot_axes, xx0, yy0, xx1, yy1, thetimelbl, xlen, $
			extramsg=extramsg

        !p.position=[xx0,yy0,xx1,yy1]
	!p.ticklen=0.03

	xlen= xlen/0.7
	xmin= -xlen
	xmax=  xlen
	ymin= -xlen
	ymax=  xlen

        plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
		      xticks= 10, yticks=10, $
                      xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'

        xyouts, xx0+0.04, yy1-0.06, thetimelbl, /normal, size= 1.2, charthick=3.0, color= 0

	if keyword_set(extramsg) then begin
		xyouts, xx0+0.04, yy1-0.10, extramsg, /normal, size= 1.2, charthick=3.0, color= 0
	endif


end




