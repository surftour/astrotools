
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Detect Shells in Phase Space
;     -------------------------------------------
;     We'll do this by plotting v_r versus r,
;     and essentially pick it out of phase space.
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro doshells, frun, snapnum, filename=filename

if not keyword_set(frun) then begin
   print, "  "
   print, "doshells, frun, snapnum, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='shells.eps'
if not keyword_set(snapnum) then snapnum=25


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; radius, linear scale
xaxistitle= "r (h!E-1!Nkpc)"
xmax = 120.0
xmin =  0.0
;xmin =  0.1

; velocity
yaxistitle= "!6V!Dr!N (km s!E-1!N)"
ymax = 500.0
ymin = -500.0



; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.98
!p.position= [0.18, 0.15, 0.98, 0.98]




; ------------------
; compute histogram
; ------------------

;NxNImage= process_one_shellphasespace_map(frun, snapnum, xmin, xmax, ymin, ymax)

;tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal



!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
	;/ylog, $
	;/xlog, $
	color= 0, $
	xrange=[xmin,xmax], $
	yrange=[ymin,ymax], $
	xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata, /noerase


; -----------------------------------------------

; raw points if we want em
;

;process_one_shellphasespace, frun, snapnum, linecolor=50, /ditchmostbound
;process_one_shellphasespace, frun, snapnum, linecolor=150, /ditchsmallr
process_one_shellphasespace, frun, snapnum, linecolor=50, /ditchsmallr
;process_one_shellphasespace, frun, snapnum, linecolor=100

; -----------------------------------------------


;xyouts, 0.75, 0.90, fload_timelbl(1,2), size=1.2, color=0, /normal, charthick=4.0
;xyouts, 0.55, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0
;xyouts, 0.75, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0



; print extras
; -------------

; zero
x=[xmin,xmax]
y=x*0.0
oplot, x, y, linestyle= 2, color= 0



; done
; -----
device, /close


end






; --------------------------------
;  process one snaps phase space
; --------------------------------
pro process_one_shellphasespace, frun, snapnum, xmin, xmax, linecolor=linecolor, msg=msg, $
					ditchmostbound=ditchmostbound, $
					ditchsmallr=ditchsmallr


ok=fload_snapshot_bh(frun,snapnum)
radius= fload_allstars_xyz('r')
print, "radius   max/min= ", max(radius), min(radius)
comvel= fload_all_comvel(1)
vel_radial= fload_allstars_v('r',comvel=comvel)
print, "radial velocity   max/min= ", max(vel_radial), min(vel_radial)
;vel_theta= fload_allstars_v('theta',comvel=comvel)
;vel_phi= fload_allstars_v('phi',comvel=comvel)
;vel_tan= fload_allstars_v('tan',comvel=comvel)

if keyword_set(ditchmostbound) then begin

        ;e= fload_gas_energy()
        ke= fload_allstars_energy(1,/kinetic)
        pe= fload_allstars_energy(1,/potential)
        e= ke + pe

        ; least bound
	;
	Energy_level= 0.5
	;
	Ns= n_elements(radius)
	print, "N_allstars= ", Ns
	e_idx= sort(e)
	sorted_e= e(e_idx)
	radius= radius(e_idx)
	vel_radial= vel_radial(e_idx)

	print, "energy max/min= ", max(e), min(e)
	print, "Energy_level= ", Energy_level
	print, "Energy at level, e[",Ns*0.8,"] = ", e[Ns*0.8]
	radius= radius[Ns*0.8:Ns-1]
	vel_radial= vel_radial[Ns*0.8:Ns-1]
	print, "N_radius= ", n_elements(radius)

endif


if keyword_set(ditchsmallr) then begin

	min_r= 10.0
	print, "min_r= ", min_r
	idx= where(radius gt min_r)
	radius= radius(idx)
	vel_radial= vel_radial(idx)
	
endif

oplot, radius, vel_radial, psym=7, color=linecolor, thick=1.0, symsize= 0.2



end







; --------------------------------
;  process one snaps phase space
;
;  ***  MAP version  ***
;
; --------------------------------
function process_one_shellphasespace_map, frun, snapnum, xmin, xmax, ymin, ymax


ok=fload_snapshot_bh(frun,snapnum)
radius= fload_allstars_xyz('r')
print, "radius   max/min= ", max(radius), min(radius)
comvel= fload_all_comvel(1)
vel_radial= fload_allstars_v('r',comvel=comvel)
print, "radial velocity   max/min= ", max(vel_radial), min(vel_radial)
;vel_theta= fload_allstars_v('theta',comvel=comvel)
;vel_phi= fload_allstars_v('phi',comvel=comvel)
;vel_tan= fload_allstars_v('tan',comvel=comvel)



; ------------------
; compute histogram
; ------------------

bins= 500

contour_makegeneralpic, radius, vel_radial, xmax, xmin, ymax, ymin, $
                                pixels= bins, $
                                NxNImage=NxNImage



return, NxNImage





end




;=================================================================





pro shell_hist, frun, snapnum, filename=filename

if not keyword_set(frun) then begin
   print, "  "
   print, "shell_hist, frun, snapnum, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='shellhist.eps'
if not keyword_set(snapnum) then snapnum=25


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; radius, linear scale
xaxistitle= "r (h!E-1!Nkpc)"
xmax = 120.0
xmin =  0.0
;xmin =  0.1

; number (histogram)
;yaxistitle= "!6V!Dr!N (km s!E-1!N)"
yaxistitle= ' '
ymax = 1.2
;ymin = 0.0
ymin = 0.00001



; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
	/ylog, $
	;/xlog, $
	color= 0, $
	xrange=[xmin,xmax], $
	yrange=[ymin,ymax], $
	xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata, /noerase


; -----------------------------------------------


; ------------------
; compute histogram
; ------------------


ok=fload_snapshot_bh(frun,snapnum)
radius= fload_allstars_xyz('r')
print, "radius   max/min= ", max(radius), min(radius)
temp= process_histogram(radius, xmax=xmax, xmin=xmin, levels=100, oplotit=50)
print, min(temp), max(temp)

; -----------------------------------------------


xyouts, 0.75, 0.90, fload_timelbl(1,2), size=1.2, color=0, /normal, charthick=4.0
;xyouts, 0.55, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0
xyouts, 0.75, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0



; print extras
; -------------


; done
; -----
device, /close


end






