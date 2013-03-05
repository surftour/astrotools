;-----------------------------------------------------------
;
;
;
;
;
;-----------------------------------------------------------
pro disk_thickness, junk, filename=filename

if not keyword_set(junk) then begin
   print, "  "
   print, "  "
   print, "disk_thickness, junk"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename = 'disk_thickness.eps'

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename=filename, colortable=4
; should I try and send a different size for this?



;-------------
;  Plot axes 
;-------------

xaxistitle= "Radius (kpc)"
xmax = 15.0
xmin = 0.0

yaxistitle= "Thickness, z!D1/2!N (kpc)"
ymax = 2.0
ymin = 0.0



!p.position= [0.20, 0.15, 0.97, 0.98]

plot, [10.0], [10.0], psym=-3, linestyle=0, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata





;-----------------------
;  Put data on graph
;-----------------------2


;do_1_to_8_0gas= 1
do_1_to_8_0gas= 0
if do_1_to_8_0gas eq 1 then begin

        ;  initial
        ; --------
        load_and_plot_one_thickness, 'data/minor/min_30_fg0.0', 0, xmax, xmin, pcolor= 0, plotthick= 8.0, plotlinest= 1, rotateby= 30.0
        xyouts, 0.45, 0.90, 'initial disk (dotted)', color= 0, charthick=2.0, size=1.50, /normal

        ;  final
        ; --------
        fload_newcolortable, 1
        load_and_plot_one_thickness, 'data/minor/min_30_fg0.0', 200, xmax, xmin, pcolor= 50, plotthick= 4.0, rotateby= 30.0
        xyouts, 0.45, 0.85, 'final (fgas= 0)', color= 50, charthick=3.0, size=1.10, /normal
        load_and_plot_one_thickness, 'data/minor/min_30', 200, xmax, xmin, pcolor= 100, plotthick= 4.0, rotateby= 30.0
        xyouts, 0.45, 0.81, 'final (fgas= 20%)', color= 100, charthick=3.0, size=1.10, /normal
        load_and_plot_one_thickness, 'data/minor/min_30_fg0.4', 200, xmax, xmin, pcolor=150, plotthick=4.0, rotateby= 30.0
        xyouts, 0.45, 0.77, 'final (fgas= 40%)', color= 150, charthick=3.0, size=1.10, /normal
endif



do_1_to_8_40gas= 1
;do_1_to_8_40gas= 0
if do_1_to_8_40gas eq 1 then begin

	;  initial
	; --------
        load_and_plot_one_thickness, 'data/minor/Sbfg0.4Imfg0.4_000', 0, xmax, xmin, pcolor= 0, plotthick= 8.0, plotlinest= 1
        xyouts, 0.25, 0.90, 'initial disk (dotted)', color= 0, charthick=2.0, size=1.50, /normal
        
	;  final
	; --------
	load_and_plot_one_thickness, 'data/minor/Sbfg0.4Imfg0.4_000', 120, xmax, xmin, pcolor= 0, plotthick= 8.0
        xyouts, 0.25, 0.85, 'final disk (solid)', color= 0, charthick=2.0, size=1.50, /normal

	load_and_plot_one_thickness, 'data/minor/Sbfg0.4Imfg0.4_000', 120, xmax, xmin, pointt=1, /do_old_stars
        xyouts, 0.25, 0.81, 'old stars', color= 150, charthick=3.0, size=1.10, /normal
	;load_and_plot_one_thickness, 'data/minor/Sbfg0.4Imfg0.4_000', 120, xmax, xmin, pointt=6, /do_new_stars
        ;xyouts, 0.28, 0.80, 'new', color= 20, charthick=2.0, size=1.10, /normal
	load_and_plot_one_thickness, 'data/minor/Sbfg0.4Imfg0.4_000', 120, xmax, xmin, pointt=3, do_postburst= 1.9
        xyouts, 0.25, 0.77, 'stars formed post merger', color= 50, charthick=3.0, size=1.10, /normal
endif


;do_1_to_8_80gas= 1
do_1_to_8_80gas= 0
if do_1_to_8_80gas eq 1 then begin

	;  initial
	; --------
        load_and_plot_one_thickness, 'data/minor/Sbfg0.8Imfg0.8_000', 0, xmax, xmin, pcolor= 0, plotthick= 8.0, plotlinest= 1
        xyouts, 0.25, 0.90, 'initial disk (dotted)', color= 0, charthick=2.0, size=1.50, /normal
        
	;  final
	; --------
	load_and_plot_one_thickness, 'data/minor/Sbfg0.8Imfg0.8_000', 120, xmax, xmin, pcolor= 0, plotthick= 8.0
        xyouts, 0.25, 0.85, 'final disk (solid)', color= 0, charthick=2.0, size=1.50, /normal

	load_and_plot_one_thickness, 'data/minor/Sbfg0.8Imfg0.8_000', 120, xmax, xmin, pointt=1, /do_old_stars
        xyouts, 0.25, 0.81, 'old stars', color= 150, charthick=3.0, size=1.10, /normal
	;load_and_plot_one_thickness, 'data/minor/Sbfg0.8Imfg0.8_000', 120, xmax, xmin, pointt=6, /do_new_stars
        ;xyouts, 0.28, 0.80, 'new', color= 20, charthick=2.0, size=1.10, /normal
	load_and_plot_one_thickness, 'data/minor/Sbfg0.8Imfg0.8_000', 120, xmax, xmin, pointt=3, do_postburst= 1.9
        xyouts, 0.25, 0.77, 'stars formed post merger', color= 50, charthick=3.0, size=1.10, /normal
endif



;do_iso= 1
do_iso= 0
if do_iso eq 1 then begin

        ;  stad
        ; --------
	fload_newcolortable, 1

	; no gas disk
        ;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.0',  0, xmax, xmin, plotthick= 1.0, pcolor=  20
        ;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.0', 10, xmax, xmin, plotthick= 2.0, pcolor=  50
        ;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.0', 20, xmax, xmin, plotthick= 3.0, pcolor=  80
        ;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.0', 30, xmax, xmin, plotthick= 4.0, pcolor= 110
        ;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.0', 40, xmax, xmin, plotthick= 5.0, pcolor= 140
        ;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.0', 60, xmax, xmin, plotthick= 6.0, pcolor= 170
        ;;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.0', 80, xmax, xmin, plotthick= 7.0, pcolor= 200
        ;xyouts, 0.28, 0.85, 'isolated disk (without gas)', color=50, charthick=2.0, size=1.50, /normal

	; 40% gas
        ;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.4',  0, xmax, xmin, plotthick= 1.0, pcolor=  20
        ;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.4', 10, xmax, xmin, plotthick= 2.0, pcolor=  50
        ;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.4', 20, xmax, xmin, plotthick= 3.0, pcolor=  80
        ;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.4', 30, xmax, xmin, plotthick= 4.0, pcolor= 110
        ;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.4', 40, xmax, xmin, plotthick= 5.0, pcolor= 140
        ;load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.4', 60, xmax, xmin, plotthick= 6.0, pcolor= 170
        ;xyouts, 0.28, 0.85, 'isolated disk (initially 40% gas)', color=50, charthick=2.0, size=1.50, /normal

	; 40% gas
        load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.4',  0, xmax, xmin, plotthick= 1.0, pcolor=  20
        load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.4', 60, xmax, xmin, plotthick= 6.0, pcolor= 170
        load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.4', 60, xmax, xmin, plotthick= 6.0, pcolor= 170, plotlinest= 1, /do_old_stars
        load_and_plot_one_thickness, 'data/minor/iso_Sb_fg0.4', 60, xmax, xmin, plotthick= 6.0, pcolor= 170, plotlinest= 2, /do_new_stars
        xyouts, 0.28, 0.85, 'isolated disk (initially 40% gas)', color=50, charthick=2.0, size=1.50, /normal


endif







;load_and_plot_seq, "/raid4/tcox/ds/vc3vc3e_2/", 1, r_kenn=r_kenn



;--------------------------------------
device, /close


end









;==========================================================================




pro load_and_plot_one_thickness, frun, snapnum, xmax, xmin, pointt=pointt, $
				plotlinest=plotlinest, $
				plotthick=plotthick, $
				pstyle=pstyle, $
				pcolor=pcolor, $
				do_old_stars=do_old_stars, $
				do_new_stars=do_new_stars, $
				do_postburst=do_postburst, $
				rotateby=rotateby

	if not keyword_set(rotateby) then rotateby= 0.0

	bins = 40

	ok=fload_snapshot_bh(frun,snapnum,/nopot_in_snap)
	bhid= fload_blackhole_id(1)
        ;bhid= bhid[0]
        ;bhid= bhid[1]
        ;bhid= 200001L
        ;bhid= 280002L   ; used for z3/b4e
        ;bhid= 400002L   ; used for ds/vc3vc3e_2
        bhid= long(min(bhid))   ; uses minimum bhid - then you don't need to know the BH ID
        center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
        print, "Blackhole ID: ", bhid
        print, "Blackhole center: ", center_bh


	x= fload_allstars_xyz('x') / 0.7
	y= fload_allstars_xyz('y') / 0.7
	z= fload_allstars_xyz('z') / 0.7
	c= fload_allstars_mass(1) / 0.7

	if keyword_set(do_old_stars) then begin
		x= fload_oldstars_xyz('x') / 0.7
		y= fload_oldstars_xyz('y') / 0.7
		z= fload_oldstars_xyz('z') / 0.7
		c= fload_oldstars_mass(1) / 0.7
	endif

	if keyword_set(do_new_stars) then begin
		x= fload_newstars_xyz('x') / 0.7
		y= fload_newstars_xyz('y') / 0.7
		z= fload_newstars_xyz('z') / 0.7
		c= fload_newstars_mass(1) / 0.7
	endif

	if keyword_set(do_postburst) then begin
		x= fload_newstars_xyz('x') / 0.7
		y= fload_newstars_xyz('y') / 0.7
		z= fload_newstars_xyz('z') / 0.7
		m= fload_newstars_mass(1) / 0.7
		age= fload_newstars_age(1)

		idx= where(age gt do_postburst)
		x= x(idx)
		y= y(idx)
		z= z(idx)
		c= m(idx)
	endif

	; rotate so that disk is in x-y plane
	;check_rotation, 1
	print, "rotateby= ", rotateby
	rotateby= rotateby * (!PI / 180.)
	print, "rotateby= ", rotateby
        x_new= cos(rotateby)*x - sin(rotateby)*z
        y_new= y
        z_new= sin(rotateby)*x + cos(rotateby)*z

	a= sqrt(x_new*x_new + y_new*y_new)
	b= z_new

	process_thickness, a, b, bins, xmax, xmin, $
                        radius, $
                        thckness, $
                        y_weighting=c

	; cheat
	thckness= thckness * 0.8

	if keyword_set(pointt) then begin
		psymsize= 1.0
		plotthick= 2.0
		plotlinest= 0    ; didn't work that will with other linestyles
		select_thispoint, pointt, thispsym, thiscolor
	endif else begin
		psymsize= 1.0
		if not keyword_set(pstyle) then thispsym= 3 else thispsym= pstyle
		if not keyword_set(plotthick) then plotthick= 2.0
		if not keyword_set(plotlinest) then plotlinest= 0
		if not keyword_set(pcolor) then thiscolor= 0 else thiscolor= pcolor
	endelse

	oplot, radius, thckness, psym=-thispsym, linestyle= plotlinest, $
				color= thiscolor, thick= plotthick, $
				symsize= psymsize



end




;-----------------------------------------------------------
;
;
;
;
;
;-----------------------------------------------------------
pro disk_temp, junk, filename=filename

if not keyword_set(junk) then begin
   print, "  "
   print, "  "
   print, "disk_temp, junk"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename = 'disk_temp.eps'

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename=filename, colortable=4
; should I try and send a different size for this?



;-------------
;  Plot axes 
;-------------

xaxistitle= "Radius (kpc)"
xmax = 15.0
xmin = 0.0

yaxistitle= "!7r!3!Dz!N (km/sec)"
ymax = 220.0
ymin = 0.0



!p.position= [0.20, 0.15, 0.97, 0.98]

plot, [10.0], [10.0], psym=-3, linestyle=0, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata





;-----------------------
;  Put data on graph
;-----------------------2

;do_1_to_8_0gas= 1
do_1_to_8_0gas= 0
if do_1_to_8_0gas eq 1 then begin

        ;  initial
        ; --------
        load_and_plot_one_temp, 'data/minor/min_30_fg0.0', 0, xmax, xmin, pcolor= 0, plotthick= 8.0, plotlinest= 1, rotateby= 30.0
        xyouts, 0.45, 0.90, 'initial disk (dotted)', color= 0, charthick=2.0, size=1.50, /normal

        ;  final
        ; --------
	fload_newcolortable, 1
        load_and_plot_one_temp, 'data/minor/min_30_fg0.0', 200, xmax, xmin, pcolor= 50, plotthick= 4.0, rotateby= 30.0
        xyouts, 0.45, 0.85, 'final (fgas= 0)', color= 50, charthick=3.0, size=1.10, /normal
        load_and_plot_one_temp, 'data/minor/min_30', 200, xmax, xmin, pcolor= 100, plotthick= 4.0, rotateby= 30.0
        xyouts, 0.45, 0.81, 'final (fgas= 20%)', color= 100, charthick=3.0, size=1.10, /normal
        load_and_plot_one_temp, 'data/minor/min_30_fg0.4', 200, xmax, xmin, pcolor=150, plotthick=4.0, rotateby= 30.0
        xyouts, 0.45, 0.77, 'final (fgas= 40%)', color= 150, charthick=3.0, size=1.10, /normal
endif


do_1_to_8_40gas= 1
;do_1_to_8_40gas= 0
if do_1_to_8_40gas eq 1 then begin

	;  initial
	; --------
        load_and_plot_one_temp, 'data/minor/Sbfg0.4Imfg0.4_000', 0, xmax, xmin, pcolor= 0, plotthick= 8.0, plotlinest= 1
        xyouts, 0.45, 0.90, 'initial disk (dotted)', color= 0, charthick=2.0, size=1.50, /normal
        
	;  final
	; --------
	load_and_plot_one_temp, 'data/minor/Sbfg0.4Imfg0.4_000', 120, xmax, xmin, pcolor= 0, plotthick= 8.0
        xyouts, 0.45, 0.85, 'final disk (solid)', color= 0, charthick=2.0, size=1.50, /normal

	load_and_plot_one_temp, 'data/minor/Sbfg0.4Imfg0.4_000', 120, xmax, xmin, pointt=1, /do_old_stars
        xyouts, 0.45, 0.81, 'old stars', color= 150, charthick=3.0, size=1.10, /normal
	;load_and_plot_one_temp, 'data/minor/Sbfg0.4Imfg0.4_000', 120, xmax, xmin, pointt=6, /do_new_stars
        ;xyouts, 0.45, 0.80, 'new', color= 20, charthick=2.0, size=1.10, /normal
	load_and_plot_one_temp, 'data/minor/Sbfg0.4Imfg0.4_000', 120, xmax, xmin, pointt=3, do_postburst= 1.9
        xyouts, 0.45, 0.77, 'stars formed post merger', color= 50, charthick=3.0, size=1.10, /normal
endif


;do_1_to_8_80gas= 1
do_1_to_8_80gas= 0
if do_1_to_8_80gas eq 1 then begin

	;  initial
	; --------
        load_and_plot_one_temp, 'data/minor/Sbfg0.8Imfg0.8_000', 0, xmax, xmin, pcolor= 0, plotthick= 8.0, plotlinest= 1
        xyouts, 0.45, 0.90, 'initial disk (dotted)', color= 0, charthick=2.0, size=1.50, /normal
        
	;  final
	; --------
	load_and_plot_one_temp, 'data/minor/Sbfg0.8Imfg0.8_000', 120, xmax, xmin, pcolor= 0, plotthick= 8.0
        xyouts, 0.45, 0.85, 'final disk (solid)', color= 0, charthick=2.0, size=1.50, /normal

	load_and_plot_one_temp, 'data/minor/Sbfg0.8Imfg0.8_000', 120, xmax, xmin, pointt=1, /do_old_stars
        xyouts, 0.45, 0.81, 'old stars', color= 150, charthick=3.0, size=1.10, /normal
	;load_and_plot_one_temp, 'data/minor/Sbfg0.8Imfg0.8_000', 120, xmax, xmin, pointt=6, /do_new_stars
        ;xyouts, 0.45, 0.80, 'new', color= 20, charthick=2.0, size=1.10, /normal
	load_and_plot_one_temp, 'data/minor/Sbfg0.8Imfg0.8_000', 120, xmax, xmin, pointt=3, do_postburst= 1.9
        xyouts, 0.45, 0.77, 'stars formed post merger', color= 50, charthick=3.0, size=1.10, /normal

endif



;do_iso= 1
do_iso= 0
if do_iso eq 1 then begin

        ;  stad
        ; --------
	fload_newcolortable, 1

	; no gas disk
        ;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.0',  0, xmax, xmin, plotthick= 1.0, pcolor=  20
        ;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.0', 10, xmax, xmin, plotthick= 2.0, pcolor=  50
        ;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.0', 20, xmax, xmin, plotthick= 3.0, pcolor=  80
        ;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.0', 30, xmax, xmin, plotthick= 4.0, pcolor= 110
        ;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.0', 40, xmax, xmin, plotthick= 5.0, pcolor= 140
        ;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.0', 60, xmax, xmin, plotthick= 6.0, pcolor= 170
        ;;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.0', 80, xmax, xmin, plotthick= 7.0, pcolor= 200
        ;xyouts, 0.28, 0.85, 'isolated disk (without gas)', color=50, charthick=2.0, size=1.50, /normal

	; 40% gas
        ;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.4',  0, xmax, xmin, plotthick= 1.0, pcolor=  20
        ;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.4', 10, xmax, xmin, plotthick= 2.0, pcolor=  50
        ;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.4', 20, xmax, xmin, plotthick= 3.0, pcolor=  80
        ;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.4', 30, xmax, xmin, plotthick= 4.0, pcolor= 110
        ;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.4', 40, xmax, xmin, plotthick= 5.0, pcolor= 140
        ;load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.4', 60, xmax, xmin, plotthick= 6.0, pcolor= 170
        ;xyouts, 0.28, 0.85, 'isolated disk (initially 40% gas)', color=50, charthick=2.0, size=1.50, /normal

	; 40% gas
        load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.4',  0, xmax, xmin, plotthick= 1.0, pcolor=  20
        load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.4', 60, xmax, xmin, plotthick= 6.0, pcolor= 170
        load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.4', 60, xmax, xmin, plotthick= 6.0, pcolor= 170, plotlinest= 1, /do_old_stars
        load_and_plot_one_temp, 'data/minor/iso_Sb_fg0.4', 60, xmax, xmin, plotthick= 6.0, pcolor= 170, plotlinest= 2, /do_new_stars
        xyouts, 0.28, 0.85, 'isolated disk (initially 40% gas)', color=50, charthick=2.0, size=1.50, /normal


endif







;load_and_plot_seq, "/raid4/tcox/ds/vc3vc3e_2/", 1, r_kenn=r_kenn



;--------------------------------------
device, /close


end









;==========================================================================




pro load_and_plot_one_temp, frun, snapnum, xmax, xmin, pointt=pointt, $
				plotlinest=plotlinest, $
				plotthick=plotthick, $
				pstyle=pstyle, $
				pcolor=pcolor, $
				do_old_stars=do_old_stars, $
				do_new_stars=do_new_stars, $
				do_postburst=do_postburst, $
				rotateby=rotateby

	if not keyword_set(rotateby) then rotateby= 0.0
	bins = 40

	ok=fload_snapshot_bh(frun,snapnum,/nopot_in_snap)
	x= fload_allstars_xyz('x') / 0.7
	y= fload_allstars_xyz('y') / 0.7
	z= fload_allstars_xyz('z') / 0.7
	vx= fload_allstars_v('x')
	vy= fload_allstars_v('y')
	vz= fload_allstars_v('z')
	c= fload_allstars_mass(1) / 0.7

	if keyword_set(do_old_stars) then begin
		x= fload_oldstars_xyz('x') / 0.7
		y= fload_oldstars_xyz('y') / 0.7
		z= fload_oldstars_xyz('z') / 0.7
		vx= fload_oldstars_v('x')
		vy= fload_oldstars_v('y')
		vz= fload_oldstars_v('z')
		c= fload_oldstars_mass(1) / 0.7
	endif

	if keyword_set(do_new_stars) then begin
		x= fload_newstars_xyz('x') / 0.7
		y= fload_newstars_xyz('y') / 0.7
		z= fload_newstars_xyz('z') / 0.7
		vx= fload_newstars_v('x')
		vy= fload_newstars_v('y')
		vz= fload_newstars_v('z')
		c= fload_newstars_mass(1) / 0.7
	endif

	if keyword_set(do_postburst) then begin
		x= fload_newstars_xyz('x') / 0.7
		y= fload_newstars_xyz('y') / 0.7
		z= fload_newstars_xyz('z') / 0.7
		vx= fload_newstars_v('x')
		vy= fload_newstars_v('y')
		vz= fload_newstars_v('z')
		m= fload_newstars_mass(1) / 0.7
		age= fload_newstars_age(1)

		idx= where(age gt do_postburst)
		x= x(idx)
		y= y(idx)
		z= z(idx)
		vx= vx(idx)
		vy= vy(idx)
		vz= vz(idx)
		c= m(idx)
	endif


	; rotate so that disk is in x-y plane
	rotateby= rotateby * (!PI / 180.)
        x_new= cos(rotateby)*x - sin(rotateby)*z
        y_new= y
        z_new= sin(rotateby)*x + cos(rotateby)*z

        vx_new= cos(rotateby)*vx - sin(rotateby)*vz
        vy_new= vy
        vz_new= sin(rotateby)*vx + cos(rotateby)*vz

	a= sqrt(x_new*x_new + y_new*y_new)
	b= vz_new


	process_disktemp, a, b, bins, xmax, xmin, $
                        radius, $
                        temp, $
                        y_weighting=c


	if keyword_set(pointt) then begin
		psymsize= 1.0
		plotthick= 2.0
		plotlinest= 0    ; didn't work that will with other linestyles
		select_thispoint, pointt, thispsym, thiscolor
	endif else begin
		psymsize= 1.0
		if not keyword_set(pstyle) then thispsym= 3 else thispsym= pstyle
		if not keyword_set(plotthick) then plotthick= 2.0
		if not keyword_set(plotlinest) then plotlinest= 0
		if not keyword_set(pcolor) then thiscolor= 0 else thiscolor= pcolor
	endelse

	oplot, radius, temp, psym=-thispsym, linestyle= plotlinest, $
				color= thiscolor, thick= plotthick, $
				symsize= psymsize



end





pro check_rotation, junk

           com=fload_center(1)
           r= fload_allstars_xyz('r', center=com)
           comvel= fload_all_comvel(1, center=com, justcenter=10.0)
           jx= fload_allstars_j(21, center=com, vcom=comvel)   ; specific j_x
           jy= fload_allstars_j(22, center=com, vcom=comvel)   ; specific j_y
           jz= fload_allstars_j(23, center=com, vcom=comvel)   ; specific j_z

           ; angular momentum within some radius
           radius= 5.0
           idx= where(r lt radius)
           if idx(0) eq -1 then stop
           jtot_x= total(jx(idx))
           jtot_y= total(jy(idx))
           jtot_z= total(jz(idx))
           jtot= sqrt(jtot_x*jtot_x + jtot_y*jtot_y + jtot_z*jtot_z)
           nx= jtot_x/jtot
           ny= jtot_y/jtot
           nz= jtot_z/jtot

           n_hat= [nx, ny, nz]
           n_hat_tot= sqrt(nx*nx + ny*ny + nz*nz)

           theta= acos(nz)*180./!PI
           phi= atan(ny,nx)*180./!PI

           print, "angular momentum vector= ", nx, ny, nz
           print, "              n_hat_tot= ", n_hat_tot
           print, "                  theta= ", theta
           print, "                    phi= ", phi


end
