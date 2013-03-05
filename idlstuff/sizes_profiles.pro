
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Profiles
;   -------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------


pro profs_bycomp, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "profs_bycomp, junk"
   print, "  "
   return
endif

filename='profiles_bycomp.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 40

yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"

xaxistitle="!6[R (kpc)]!E1/4!N"
xmax = 2.4
xmin = 0.4
ymax = 6.0
ymin = -0.2




; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, /nodata, $
        ;/ylog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle




; -----------------------------------------------

; original profile

;process_and_plot_sd, "data1/remergers/iso_b3e", 10, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=0, thislinest= 2
;xyouts, 0.7, 0.88, "original (b3e)", /normal, charthick=3.0, size=1.2, color= 0


; -----------------------------------------------


;frun= "data1/remergers/b3eb3e2" & msg="b3e"
;frun= "data1/remergers/b3esph" & msg="sph"
frun= "data1/remergers/b3ev4" & msg="v4"
;frun= "data1/remergers/b3ev3" & msg="v3"
;frun= "data1/remergers/b3ev2" & msg="v2"

ok= fload_snapshot_bh(frun,0,/nopot_in_snap,/skip_center)
init_bhids= fload_blackhole_id(1,n_bh=n_bh)
print, "Blackhole ID's: ", init_bhids

process_and_plot_sd, frun, -99, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=0
xyouts, 0.7, 0.82, fload_frun(frun,/runidonly), /normal, charthick=3.0, size=1.2, color= 0


;
; galaxy 1 (b3e)
; ---------------
startid= long(1)
numpart= init_bhids[0]
x=fload_1gal_allstars_xyz('x', startid, numpart)
y=fload_1gal_allstars_xyz('y', startid, numpart)
z=fload_1gal_allstars_xyz('z', startid, numpart)
m=fload_1gal_allstars_mass(startid, numpart)

process_and_plot_sd, "loaded", 0, "special", xmin, xmax, bins, linecolor=150, $
			special_m=m, special_x=x, special_y=y, special_z=z, /x_is_devac

xyouts, 0.7, 0.78, '!6b3e', /normal, charthick=3.0, size= 1.2, color= 150



;
; galaxy 2
; -----------
startid= init_bhids[0]+1
numpart= init_bhids[1]-init_bhids[0]
x=fload_1gal_allstars_xyz('x', startid, numpart)
y=fload_1gal_allstars_xyz('y', startid, numpart)
z=fload_1gal_allstars_xyz('z', startid, numpart)
m=fload_1gal_allstars_mass(startid, numpart)

process_and_plot_sd, "loaded", 0, "special", xmin, xmax, bins, linecolor=100, $
                        special_m=m, special_x=x, special_y=y, special_z=z, /x_is_devac

xyouts, 0.7, 0.74, msg, /normal, charthick=3.0, size= 1.2, color= 100
xyouts, 0.7, 0.74, msg+"-total", /normal, charthick=3.0, size= 1.2, color= 100


; new stars
x=fload_1gal_newstars_xyz('x', startid, numpart)
y=fload_1gal_newstars_xyz('y', startid, numpart)
z=fload_1gal_newstars_xyz('z', startid, numpart)
m=fload_1gal_newstars_mass(startid, numpart)

process_and_plot_sd, "loaded", 0, "special", xmin, xmax, bins, linecolor=50, thislinest= 2, $
                        special_m=m, special_x=x, special_y=y, special_z=z, /x_is_devac

xyouts, 0.7, 0.70, msg+"-new stars", /normal, charthick=3.0, size= 1.2, color= 50


; old stars
x=fload_1gal_oldstars_xyz('x', startid, numpart)
y=fload_1gal_oldstars_xyz('y', startid, numpart)
z=fload_1gal_oldstars_xyz('z', startid, numpart)
m=fload_1gal_oldstars_mass(startid, numpart)

process_and_plot_sd, "loaded", 0, "special", xmin, xmax, bins, linecolor=200, thislinest= 1, $
                        special_m=m, special_x=x, special_y=y, special_z=z, /x_is_devac

xyouts, 0.7, 0.66, msg+"-old stars", /normal, charthick=3.0, size= 1.2, color= 200



; print extras
; -------------
if not keyword_set(smoothlen) then smoothlen= 0.1/0.7 else smoothlen= smoothlen/0.7
x=[smoothlen,smoothlen]
x=x^(0.25)
y=[ymin,ymax]
oplot, x, y, linestyle=1, color= 0




;  done
; -------
device, /close



end





;================================================================================
;================================================================================





; -------------------------------
;
;  Overplot devac profile for one
;  galaxy (and one snap).
;
; -------------------------------
pro do_one_galaxy_devac, frun, snapnum, linecolor=linecolor, msg=msg, $
					fitdevac=fitdevac, $
					fitsersic=fitsersic, $
					includeerr=includeerr, $
					frommanyproj=frommanyproj

; make sure these are same as 
; calling routine
;
bins = 40

xmax = 2.4
;xmax = 4.2
xmin = 0.4
;xmin = 0.5
;ymax = 6.0
;ymin = -1.0

	;
	; total (all star) profile
	; --------------------------

	;process_and_plot_sd, x, y, z, m, xmin, xmax, bins, linecolor=linecolor, /calcfit
	process_and_plot_sd, frun, snapnum, "allstars", xmin, xmax, bins, linecolor=linecolor, $
					fitdevac=fitdevac, /x_is_devac, $
					fitsersic=fitsersic, $
					includeerr=includeerr, $
					frommanyproj=frommanyproj, $
					h= 0.7

	if keyword_set(msg) then begin
	  ;xyouts, 0.75, 0.90, msg, size=1.5, color=linecolor, /normal, charthick=4.0
	endif

end





;================================================================================





pro profs_devac, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: profs_devac"
   print, "  "
   return
endif

filename= 'profiles_devac.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
;setup_plot_stuff, 'ps', filename=filename, colortable=2
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 40

xmax = 2.4
xmin = 0.4
ymax = 6.0
ymin = -0.2


xaxistitle= "!6[R (kpc)]!E1/4!N"
yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"


if not keyword_set(snapnum) then snapnum= 25

; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
	;/ylog, $
	color= 0, $
	xrange=[xmin,xmax], $
	yrange=[ymin,ymax], $
	xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata


; -----------------------------------------------



;process_and_plot_sd, "data1/remergers/iso_b3e", 1, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=0, thislinest= 2
process_and_plot_sd, "data1/remergers/iso_b3e", 0, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=0, thislinest= 2, /includeerr
xyouts, 0.7, 0.82, "b3e", /normal, charthick=3.0, size=1.2, color= 0

;process_and_plot_sd, "data1/remergers/iso_sph", 1, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=50
process_and_plot_sd, "data1/remergers/iso_sph", 0, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=50, /includeerr
xyouts, 0.7, 0.78, "sph", /normal, charthick=3.0, size=1.2, color= 50

;process_and_plot_sd, "data1/remergers/iso_v4x3", 1, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=150
process_and_plot_sd, "data1/remergers/iso_v4x3", 0, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=150, /includeerr
xyouts, 0.7, 0.74, "v4x3", /normal, charthick=3.0, size=1.2, color= 150
;process_and_plot_sd, "data1/remergers/iso_v3x2", 1, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=100
process_and_plot_sd, "data1/remergers/iso_v3x2", 0, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=100, /includeerr
xyouts, 0.7, 0.70, "v3x2", /normal, charthick=3.0, size=1.2, color= 100
;process_and_plot_sd, "data1/remergers/iso_v2v2", 1, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=200
process_and_plot_sd, "data1/remergers/iso_v2v2", 0, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=200, /includeerr
xyouts, 0.7, 0.66, "v2v2", /normal, charthick=3.0, size=1.2, color= 200


; -----------------------------------------------


;process_and_plot_sd, "data1/remergers/iso_b3e", 10, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=0, thislinest= 2
;process_and_plot_sd, "data1/remergers/iso_b3e", 10, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=0, thislinest= 2, /includeerr
;xyouts, 0.7, 0.88, "b3e", /normal, charthick=3.0, size=1.2, color= 0

;process_and_plot_sd, "data1/remergers/b3eb3e2", -99, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=0
;xyouts, 0.7, 0.82, "b3eb3e2", /normal, charthick=3.0, size=1.2, color= 0
;process_and_plot_sd, "data1/remergers/b3esph", -99, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=50
;xyouts, 0.7, 0.78, "b3esph", /normal, charthick=3.0, size=1.2, color= 50
;process_and_plot_sd, "data1/remergers/b3ev4", -99, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=150
;;process_and_plot_sd, "data1/remergers/b3ev4", -99, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=150, /includeerr
;xyouts, 0.7, 0.74, "b3ev4", /normal, charthick=3.0, size=1.2, color= 150 
;process_and_plot_sd, "data1/remergers/b3ev3", -99, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=100
;xyouts, 0.7, 0.70, "b3ev3", /normal, charthick=3.0, size=1.2, color= 100
;process_and_plot_sd, "data1/remergers/b3ev2", -99, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=200
;xyouts, 0.7, 0.66, "b3ev2", /normal, charthick=3.0, size=1.2, color= 200
   

; -----------------------------------------------


;fload_newcolortable, 1
;xyouts, 0.64, 0.82, 'T= ', /normal, charthick=3.0, size=1.2, color= 10
;xyouts, 0.76, 0.82, 'Gyr', /normal, charthick=3.0, size=1.2, color= 10
;process_and_plot_sd, "data1/remergers/iso_b3e", 10, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=2
;xyouts, 0.3, 0.88, "b3e (ErrTolIntAccuracy= 0.005)", /normal, charthick=3.0, size=1.2, color= 0
;process_and_plot_sd, "data1/remergers/iso_b3e_2", 1, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=2
;xyouts, 0.3, 0.88, "b3e_2 (ErrTolIntAccuracy= 0.0005)", /normal, charthick=3.0, size=1.2, color= 0
;xyouts, 0.7, 0.82, fload_timelbl(0.7,1,/noteq), /normal, charthick=3.0, size=1.2, color= 2

;process_and_plot_sd, "data1/remergers/iso_b3e", 70, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=40
;process_and_plot_sd, "data1/remergers/iso_b3e_2", 10, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=40
;xyouts, 0.7, 0.78, fload_timelbl(0.7,1,/noteq), /normal, charthick=3.0, size=1.2, color= 40

;process_and_plot_sd, "data1/remergers/iso_b3e", 140, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=80
;process_and_plot_sd, "data1/remergers/iso_b3e_2", 20, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=80
;xyouts, 0.7, 0.74, fload_timelbl(0.7,1,/noteq), /normal, charthick=3.0, size=1.2, color= 80

;process_and_plot_sd, "data1/remergers/iso_b3e", 210, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=120
;process_and_plot_sd, "data1/remergers/iso_b3e_2", 30, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=120
;xyouts, 0.7, 0.70, fload_timelbl(0.7,1,/noteq), /normal, charthick=3.0, size=1.2, color= 120

;process_and_plot_sd, "data1/remergers/iso_b3e", 280, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=180
;process_and_plot_sd, "data1/remergers/iso_b3e_2", 40, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=180
;xyouts, 0.7, 0.66, fload_timelbl(0.7,1,/noteq), /normal, charthick=3.0, size=1.2, color= 180

;process_and_plot_sd, "data1/remergers/iso_b3e", 350, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=220
;process_and_plot_sd, "data1/remergers/iso_b3e_2", 50, "allstars", xmin, xmax, bins, /x_is_devac, linecolor=220
;xyouts, 0.7, 0.62, fload_timelbl(0.7,1,/noteq), /normal, charthick=3.0, size=1.2, color= 220


; -----------------------------------------------


; print extras
; -------------
if not keyword_set(smoothlen) then smoothlen= 0.1/0.7 else smoothlen= smoothlen/0.7
x=[smoothlen,smoothlen]
x=x^(0.25)
y=[ymin,ymax]
oplot, x, y, linestyle=1, color= 0

;smoothlen= 0.4
;x=[smoothlen,smoothlen]
;x=x^(0.25)
;y=[ymin,ymax]
;oplot, x, y, linestyle=1, color= 0




; done
; -----
device, /close


end







;=================================================================================
;=================================================================================
;
;
;
;
;=================================================================================
;=================================================================================
;
;
;
;
;=================================================================================
;=================================================================================



; --------------------------------
;  6 panel devac profile figure
;
; --------------------------------
pro plot_devacouleur_6, junk, smoothlen=smoothlen, filename=filename, snapnum=snapnum

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: plot_devacouleur_6, junk"
   print, "  "
   return
endif

if not keyword_set(smoothlen) then smoothlen= 0.1
if not keyword_set(filename) then filename= 'devac.eps'


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1, newxsize=30, newysize=20
;setup_plot_stuff, 'ps', filename=filename, colortable=2, newxsize=30, newysize=20
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=30, newysize=20



; -----------------------------
; set up constants
; -----------------------------

bins = 50

xmax = 2.4
;xmax = 4.2
xmin = 0.5
ymax = 5.3
ymin = -1.0


xaxistitle= "!6[R (kpc)]!E1/4!N"
yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"


frun1= "/raid4/tcox/gfs/vc3vc3u_k"   & msg1='1.0'
frun2= "/raid4/tcox/gfs/vc3vc3v_k"   & msg2='0.8'
frun3= "/raid4/tcox/gfs/vc3vc3w_k"   & msg3='0.6'
frun4= "/raid4/tcox/gfs/vc3vc3x2_k"   & msg4='0.4'
frun5= "/raid4/tcox/gfs/vc3vc3y_k"   & msg5='0.2'
frun6= "/raid4/tcox/gfs/vc3vc3z_k"   & msg6='0.05'
cfrun= "/raid4/tcox/cvc3vc3k"

if not keyword_set(snapnum) then snapnum= 30

x0= 0.08
xs= 0.30   ; assumes 3 panels
x1= x0+xs
x2= x0+xs+xs
x3= x0+xs+xs+xs
    
y0= 0.08
ys= 0.45   ; assumes 2 panels
y1= y0+ys
y2= y0+ys+ys



; smoothing length
; -----------------
smoothlen= 0.2
x=[smoothlen,smoothlen]
xsmooth=x^(0.25)
ysmooth=[ymin,ymax]


; ------------------
;    1
; ------------------
!p.position= [x0, y1, x1, y2]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	xtickformat='(a1)', ytitle=yaxistitle, /noerase

;do_one_galaxy_devac, frun1, snapnum, linecolor=150, msg=fload_getid(frun1)
do_one_galaxy_devac, frun1, snapnum, linecolor=150, /oneproj
do_one_galaxy_devac, cfrun, snapnum, linecolor= 0, /oneproj

oplot, xsmooth, ysmooth, linestyle=1, color= 0

oplot, [1.2,1.4], [4.6,4.6], psym=-3, color= 0, thick=6.0, linestyle= 1
xyouts, 1.50, 4.5, 'collisionless', /data, size= 1.5, color= 0, charthick=5.0

oplot, [1.2,1.4], [4.1,4.1], psym=-3, color= 150, thick=8.0, linestyle= 0
xyouts, 1.50, 4.0, 'f= ', /data, size= 1.5, color= linecolor, charthick=5.0
xyouts, 1.80, 4.0, msg1, /data, size= 1.5, color= linecolor, charthick=5.0


; ------------------
;    2
; ------------------
!p.position= [x1, y1, x2, y2]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase

do_one_galaxy_devac, frun2, snapnum, linecolor=150, /oneproj
;do_one_galaxy_devac, frun2, snapnum, linecolor= 150, /fitsersic, /includeerr
do_one_galaxy_devac, cfrun, snapnum, linecolor= 0, /oneproj

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 1.80, 4.0, msg2, /data, size= 1.5, color= linecolor, charthick=5.0


; ------------------
;    3
; ------------------
!p.position= [x2, y1, x3, y2]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase

do_one_galaxy_devac, frun3, snapnum, linecolor=150, /oneproj
do_one_galaxy_devac, cfrun, snapnum, linecolor= 0, /oneproj

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 1.80, 4.0, msg3, /data, size= 1.5, color= linecolor, charthick=5.0


; ------------------
;    4
; ------------------
!p.position= [x0, y0, x1, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /noerase

do_one_galaxy_devac, frun4, snapnum, linecolor=150, /oneproj
do_one_galaxy_devac, cfrun, snapnum, linecolor= 0, /oneproj

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 1.80, 4.0, msg4, /data, size= 1.5, color= linecolor, charthick=5.0


; ------------------
;    5
; ------------------
!p.position= [x1, y0, x2, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytickformat='(a1)', /noerase

do_one_galaxy_devac, frun5, snapnum, linecolor=150, /oneproj
do_one_galaxy_devac, cfrun, snapnum, linecolor= 0, /oneproj

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 1.80, 4.0, msg5, /data, size= 1.5, color= linecolor, charthick=5.0


; ------------------
;    6
; ------------------
!p.position= [x2, y0, x3, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytickformat='(a1)', /noerase

do_one_galaxy_devac, frun6, snapnum, linecolor=150, /oneproj
do_one_galaxy_devac, cfrun, snapnum, linecolor= 0, /oneproj

oplot, xsmooth, ysmooth, linestyle=1, color= 0

xyouts, 1.80, 4.0, msg6, /data, size= 1.5, color= linecolor, charthick=5.0





; done
; -----
device, /close


end





;================================================================================
;================================================================================
;================================================================================





; --------------------------------
;  Profiles
;
; --------------------------------
pro profs, junk, smoothlen=smoothlen, filename=filename

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: profs"
   print, "  "
   return
endif

if not keyword_set(smoothlen) then smoothlen= 0.1
if not keyword_set(filename) then filename= 'profiles.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 40

;xmax = 50.0
xmax = 33.0    ; now it matches 2.4^4.0 (i.e., the devac plot range)
;xmax = 20.0
xmin = 0.0
ymax = 6.0
;ymax = 4.5
;ymax = 3.5
;ymin = 1.5
ymin = -0.2
;ymin = -1.0


xaxistitle= "!6R (kpc)"
yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"
;yaxistitle= "!4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"



; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
	;/ylog, $
	color= 0, $
	xrange=[xmin,xmax], $
	yrange=[ymin,ymax], $
	xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata


; -----------------------------------------------


;process_and_plot_sd, "data1/remergers/iso_b3e", 10, "allstars", xmin, xmax, bins, linecolor=0, thislinest= 2
;xyouts, 0.2, 0.48, "b3e", /normal, charthick=3.0, size=1.2, color= 0

;process_and_plot_sd, "data1/remergers/iso_sph", 10, "allstars", xmin, xmax, bins, linecolor=50
;xyouts, 0.2, 0.44, "sph", /normal, charthick=3.0, size=1.2, color= 50

;process_and_plot_sd, "data1/remergers/iso_v2v2", 10, "allstars", xmin, xmax, bins, linecolor=100
;xyouts, 0.2, 0.40, "v2v2", /normal, charthick=3.0, size=1.2, color= 100
;process_and_plot_sd, "data1/remergers/iso_v3x2", 10, "allstars", xmin, xmax, bins, linecolor=150
;xyouts, 0.2, 0.36, "v3x2", /normal, charthick=3.0, size=1.2, color= 150
;process_and_plot_sd, "data1/remergers/iso_v4x3", 10, "allstars", xmin, xmax, bins, linecolor=200
;xyouts, 0.2, 0.32, "v4x3", /normal, charthick=3.0, size=1.2, color= 200


; -----------------------------------------------


process_and_plot_sd, "data1/remergers/iso_b3e", 10, "allstars", xmin, xmax, bins, linecolor=0, thislinest= 2
xyouts, 0.7, 0.88, "b3e", /normal, charthick=3.0, size=1.2, color= 0

process_and_plot_sd, "data1/remergers/b3eb3e2", -99, "allstars", xmin, xmax, bins, linecolor=0
xyouts, 0.7, 0.82, "b3eb3e2", /normal, charthick=3.0, size=1.2, color= 0
process_and_plot_sd, "data1/remergers/b3esph", -99, "allstars", xmin, xmax, bins, linecolor=50
xyouts, 0.7, 0.78, "b3esph", /normal, charthick=3.0, size=1.2, color= 50
process_and_plot_sd, "data1/remergers/b3ev4", -99, "allstars", xmin, xmax, bins, linecolor=150
xyouts, 0.7, 0.74, "b3ev4", /normal, charthick=3.0, size=1.2, color= 150
process_and_plot_sd, "data1/remergers/b3ev3", -99, "allstars", xmin, xmax, bins, linecolor=100
xyouts, 0.7, 0.70, "b3ev3", /normal, charthick=3.0, size=1.2, color= 100
process_and_plot_sd, "data1/remergers/b3ev2", -99, "allstars", xmin, xmax, bins, linecolor=200
xyouts, 0.7, 0.66, "b3ev2", /normal, charthick=3.0, size=1.2, color= 200


; -----------------------------------------------


; print extras
; -------------
if not keyword_set(smoothlen) then smoothlen= 0.1/0.7 else smoothlen= smoothlen/0.7
x=[smoothlen,smoothlen]
;x=x^(0.25)
y=[ymin,ymax]
oplot, x, y, linestyle=1, color= 0

;smoothlen= 0.4
;x=[smoothlen,smoothlen]
;x=x^(0.25)
;y=[ymin,ymax]
;oplot, x, y, linestyle=1, color= 0




; done
; -----
device, /close


end


