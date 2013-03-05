pro determine_orbit, frun

if not keyword_set(frun) then begin
	print, "  "
	print, "determine_orbit, frun"
	print, "  "
	return
endif

	xlen= 1000.0
	sendto= 'ps'


	; read centers
	; -------------
	read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=frun+"/centers.txt"

	cp_time= transpose(cp_time)
	cp_cen1= transpose(cp_cen1)
	cp_cen2= transpose(cp_cen2)
	comvel_1= 0.0*cp_cen1
	comvel_2= 0.0*cp_cen2


        ; get startup info
        ; ------------------
	ok=fload_snapshot_bh(frun,0)
	;ok=fload_snapshot_bh(frun,0,/nopot_in_snap)

        ;startid= 200001
        ;numpart= 200000
        ;startid= 1L
        ;numpart= 140000L
        startid= 1L
        numpart= 500000L
	thiscenter= cp_cen1[0,*]
	comvel_1[0,*]= fload_1gal_comvel(startid,numpart,center=thiscenter,rfact=1.0)
	print, fload_1gal_center(startid,numpart)
	m1= fload_1gal_mass(startid,numpart)
	print, "Galaxy 1 --------------"
	print, "mass= ",m1
	print, "center= ",thiscenter
	print, "comvel= ", comvel_1[0,*]

        ;startid= 1
        ;numpart= 200000
        ;startid= 140001L
        ;numpart= 220000L
        startid= 500001L
        numpart= 800000L
	thiscenter= cp_cen2[0,*]
	comvel_2[0,*]= fload_1gal_comvel(startid,numpart,center=thiscenter,rfact=1.0)
	print, fload_1gal_center(startid,numpart)
	m2= fload_1gal_mass(startid,numpart)
	print, "Galaxy 2 --------------"
	print, "mass= ",m2
	print, "center= ",thiscenter
	print, "comvel= ", comvel_2[0,*]


	; transform the non-h units
	cp_time= cp_time/0.7
	cp_cen1= cp_cen1/0.7
	cp_cen2= cp_cen2/0.7
	


;	plot_mergerseparation_vs_time, cp_time, cp_cen1, cp_cen1, cp_cen1, cp_cen2, cp_cen2, cp_cen2, $
;			comvel_1, comvel_2, m1, m2, $
;			xlen=xlen, $
;			frun, sendto, filename=frun+"/csep.eps"
	plot_mergerseparation_vs_time, cp_time, cp_cen2, cp_cen2, cp_cen2, cp_cen1, cp_cen1, cp_cen1, $
			comvel_2, comvel_1, m2, m1, $
			xlen=xlen, $
			frun, sendto, filename=frun+"/csep.eps"

end















;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Merger Seperation vs. time
;     -------------------------------------------
;  plot the distance between the two interacting galaxies as a function
;  of time, do this both using all the mass to find the center and
;  just the baryonic mass
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro plot_mergerseparation_vs_time, time, center_1, center_b_1, center_d_1, center_2, center_b_2, center_d_2, $
			comvel_1, comvel_2, M1, M2, $
			xlen=xlen, $
                        frun, sendto, filename=filename



if not keyword_set(sendto) then sendto='ps'
if not keyword_set(xlen) then xlen=100.0
if n_elements(time) lt 2 then begin
   print, "  "
   print, "PROBLEM: plot_mergerseparation_vs_time"
   print, "  "
   return
endif


initialize_plotinfo, 1

setup_plot_stuff, sendto, filename=filename, colortable=4



; -------------------------------
;  Set Parameters
; -------------------------------

G= 43007.1     ; use gadget units

xmin = -xlen
xmax = xlen
ymin = -xlen
ymax = xlen

zmin = -3
zmax = 3

timemin = 0.0
timemax = max(time)

Mtot= M1 + M2
mu= M1*M2/Mtot

x = center_1[*,0] - center_2[*,0]
y = center_1[*,1] - center_2[*,1]
z = center_1[*,2] - center_2[*,2]

rdiff = sqrt(x*x + y*y + z*z)

vx= comvel_1[*,0] - comvel_2[*,0]
vy= comvel_1[*,1] - comvel_2[*,1]
vz= comvel_1[*,2] - comvel_2[*,2]

velrel= sqrt(vx*vx + vy*vy + vz*vz)




; Baryons and Dark Matter Seperately
; ------------------------------------
x_b = center_b_1[*,0] - center_b_2[*,0]
y_b = center_b_1[*,1] - center_b_2[*,1]
z_b = center_b_1[*,2] - center_b_2[*,2]

rdiff_b = sqrt(x_b*x_b + y_b*y_b + z_b*z_b)



x_d = center_d_1[*,0] - center_d_2[*,0]
y_d = center_d_1[*,1] - center_d_2[*,1]
z_d = center_d_1[*,2] - center_d_2[*,2]

rdiff_d = sqrt(x_d*x_d + y_d*y_d + z_d*z_d)









; initial separation and angle
; ----------------------------
r_s= rdiff[0]
phi_s= acos(x[0]/r_s)    ; in radians



; initial velcities
; ------------------
;phidot_s= (vx[0]/r_s)*sin(phi_s)  + (vy[0]/r_s)*cos(phi_s)
phidot_s= (vx[0]/r_s)*(y[0]/r_s) - (vy[0]/r_s)*(x[0]/r_s)
;rdot_s= vx[0]*cos(phi_s) + vy[0]*sin(phi_s)
rdot_s= vx[0]*x[0]/r_s + vy[0]*y[0]/r_s

v2= rdot_s*rdot_s + r_s*r_s*phidot_s*phidot_s



;constants of the orbit
p= phidot_s*phidot_s*r_s*r_s*r_s*r_s/(G*Mtot)

l= mu*r_s*r_s*phidot_s    ; total angular momentum
if l lt 0 then l=-l

;ecc= (p/r_s - 1.0)/cos(phi_s)   ; eccentricity
ecc= (p/r_s - 1.0)*r_s/x[0]

;ecc_another= (phidot_s*rdot_s*r_s*r_s)/(G*Mtot*sin(phi_s))

KEradial= 0.5*mu*rdot_s*rdot_s
KE0= 0.5*mu*velrel[0]*velrel[0]
PE0= -G*M1*M2/r_s
Energy= KE0 + PE0

Energy_another= (ecc*ecc - 1.0) * G*M1*M2 / (2.0*p)

a= -G*M1*M2/(2.0*Energy)       ; semi-major axis

rmin= p/(1+ecc)

;stop

; ----------------------------------
;   Try this plot thingy
; ----------------------------------

; plot center separation vs. time

!p.position= [0.18, 0.15, 0.95, 0.95]


plot, [1.0], [1.0], psym=-3, $
        xrange=[timemin,timemax], $
	yrange=[0,ymax], $
        xstyle=1, ystyle=1, $
        color= 0, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle="!6Time (Gyr)", $
        ;ytitle="Center Separation (kpc)", $
        ytitle="!6MW/M31 Separation (kpc)", $
        /nodata


if (sendto EQ 'ps') then begin
	oplot, time, rdiff_b, psym=-3, color= 50, thick=3.0, linestyle= 1
	oplot, time, rdiff_d, psym=-3, color= 100, thick=3.0, linestyle= 1
	oplot, time, rdiff, psym=-3, color= 150, thick=5.0
endif else begin
	oplot, time, rdiff, psym=-3, color= getcolor('red')
endelse


;xyouts, 0.65, 0.88, fload_fid(1), /normal, size= 1.7, color= 0, charthick= 3.0


plotcurrentsep= 0
plotcurrentsep= 1
if plotcurrentsep eq 1 then begin
	currsep= 700.0
	x= [min(time), max(time)]
	y= time*0.0 + currsep
	oplot, x, y, psym=-3, color= 0, thick= 2.0, linestyle= 1
endif




plotextras= 0
;plotextras= 1
if plotextras eq 1 then begin

	; print initial separation
	rslbl='R!Dstart!N='+strcompress(string(r_s),/remove_all)
	rslbl=strmid(rslbl,0,16)
	;rslbl= rslbl +' kpc/h'
	xyouts, 0.70, 0.82, rslbl, /normal, size= 1.2, color=0, charthick=2.5

	rminlbl='R!Dmin!N='+strcompress(string(rmin),/remove_all)
	rminlbl=strmid(rminlbl,0,13)
	xyouts, 0.70, 0.78, rminlbl, /normal, size= 1.2, color=0, charthick=2.5

	plbl='p='+strcompress(string(p),/remove_all)
	plbl=strmid(plbl,0,6)
	xyouts, 0.70, 0.74, plbl, /normal, size= 1.2, color=0, charthick=2.5

	ecclbl='e='+strcompress(string(ecc),/remove_all)
	ecclbl=strmid(ecclbl,0,6)
	xyouts, 0.70, 0.70, ecclbl, /normal, size= 1.2, color=0, charthick=2.5

	; angular momentum
	llbl='l='+strcompress(string(l),/remove_all)
	llbl=strmid(llbl,0,10)
	xyouts, 0.70, 0.66, llbl, /normal, size= 1.2, color=0, charthick=2.5

	; energy
	;elbl='E='+strcompress(string(Energy),/remove_all)
	;elbl=strmid(elbl,0,16)
	;xyouts, 0.70, 0.66, elbl, /normal, size= 1.2, color=0

	elbl='KE/PE='+strcompress(string(-KE0/PE0),/remove_all)
	elbl=strmid(elbl,0,10)
	xyouts, 0.70, 0.62, elbl, /normal, size= 1.2, color=0, charthick=2.5

	elbl='KE!Drad!N/KE='+strcompress(string(KEradial/KE0),/remove_all)
	elbl=strmid(elbl,0,17)
	xyouts, 0.70, 0.58, elbl, /normal, size= 1.2, color=0, charthick=2.5

	elbl='logKE='+strcompress(string(alog10(KE0)),/remove_all)
	elbl=strmid(elbl,0,16)
	xyouts, 0.70, 0.54, elbl, /normal, size= 1.2, color=0, charthick=2.5

	; semi-major axis
	;albl='a='+strcompress(string(a),/remove_all)
	;albl=strmid(albl,0,6)
	;xyouts, 0.70, 0.62, albl, /normal, size= 1.2, color=0


endif




; -------------------------------

if (sendto EQ 'ps') then device, /close




end









;----------------------------------
;
;
;
;----------------------------------
pro show_centers_time, junk


if not keyword_set(junk) then begin
	print, "  "
	print, "  show_centers_time, junk"
	print, "  "
endif


filename='centers_eh.eps'

timemin= 0.0
;timemin= 1.0
timemax= 1.2
;timemax= 2.0
;timemax= 4.0
;timemax= 4.25
;timemax= 4.6
;timemax= 13.0
;timemax= 25.0

;ymax= 1000.0
;ymax= 200.0
ymax= 120.0
;ymax= 1.4
ymin= 0.1
;ymin= 0


; ----------------------------------
;   Try this plot thingy
; ----------------------------------

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename=filename, colortable= 4



!p.position= [0.19, 0.15, 0.97, 0.98]

plot, [1.0], [1.0], psym=-3, $
        xrange=[timemin,timemax], $
	yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
	ytickformat='exp_label', /ylog, $
        color= 0, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle="!6Time (Gyr)", $
        ;ytitle="Center Separation (kpc)", $
        ytitle="Center Separation (!8h!6!E-1!Nkpc)", $
        ;ytitle="MW/M31 Separation (kpc)", $
        ;ytitle="!6MW/M31 Separation (Mpc)", $
        /nodata



; go through files
; ------------------

x0= 0.65
y0= 0.90

;process_onesim_centers, "/raid4/tcox/vc3vc3", 4, x0, y0
;process_onesim_centers, "/raid4/tcox/vc3rem_vc3", 6, x0, y0-0.03
;process_onesim_centers, "/raid4/tcox/vc3rem_vc3rem", 3, x0, y0-0.06
;process_onesim_centers, "/raid4/tcox/vc3rem_vc2", 2, x0, y0-0.09
;process_onesim_centers, "/raid4/tcox/vc3rem_vc1", 7, x0, y0-0.12
;process_onesim_centers, "/raid4/tcox/vc3rem_vc1a", 1, x0, y0-0.15
;process_onesim_centers, "/raid4/tcox/vc3rem_vc1b", 5, x0, y0-0.18
;process_onesim_centers, "/raid4/tcox/vc3rem_vc1c", 8, x0, y0-0.21

; -------------------------------

;process_onesim_centers, "/raid4/tcox/z5/jyDfg0.4_1", 2, 0.3, 0.35, /ylog, msg='tilted'
process_onesim_centers, "/raid4/tcox/z5/jyDfg0.4_3", 3, 0.4, 0.88, /ylog, msg='retrograde-retrograde'

;process_onesim_centers, "/raid4/tcox/vc3vc3e_2", 2, 0.3, 0.35, /ylog, msg='tilted'
;process_onesim_centers, "/raid4/tcox/vc3vc3h_2", 3, 0.3, 0.30, /ylog, msg='prograde-prograde'
;xyouts, 0.3, 0.25, '!7D!6t!Dsnap!N= 7 Myr', charthick=2.0, size=1.2, color= 0, /normal

x=[timemin,timemax]
y=[7.0,7.0]
oplot, x, y, psym=-3, linestyle=1, thick=3.0, color= 0

; -------------------------------

;process_onesim_centers, "/raid4/tcox/sbw/sb7", 1, x0, y0-0.03
;process_onesim_centers, "/raid4/tcox/sbw/sb8", 2, x0, y0-0.06
;process_onesim_centers, "/raid4/tcox/sbw/sb9", 3, x0, y0-0.09
;process_onesim_centers, "/raid4/tcox/sbw/sb10", 4, x0, y0-0.12
;process_onesim_centers, "/raid4/tcox/sbw/sb14", 5, x0, y0-0.15
;process_onesim_centers, "/raid4/tcox/sbw/sb17", 6, x0, y0-0.18

; -------------------------------

;process_onesim_centers, "/raid4/tcox/localgroup/v4", 2, x0, y0
;process_onesim_centers, "/raid4/tcox/localgroup/v5", 1, x0, y0-0.03
;process_onesim_centers, "/raid4/tcox/localgroup/v5_noigm", 2, x0, y0-0.08
;process_onesim_centers, "/raid4/tcox/localgroup/v6", 3, x0, y0-0.06

;x=[timemin,timemax]
;y=[714.0,714.0]
;y=y/1000.0
;oplot, x, y, psym=-3, linestyle=1, thick=3.0, color= 0

;xyouts, 0.57, 0.60, 'separation today', charthick=2.0, size=1.2, color= 0, /normal

; -------------------------------

device, /close


end








; do the dirty work
; --------------------
pro process_onesim_centers, frun, pointselection, x0, y0, ylog=ylog, msg=msg



;  open (or filled) square
; --------------------------
if pointselection eq 1 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 8
        symcolor= 150
endif

;  x's
; -----
if pointselection eq 2 then begin
        symsel= 7
        symcolor= 100
endif

;  open circle
; --------------
if pointselection eq 3 then begin
        symsize= 1.0
        ;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
        symsel= 8
        symcolor= 50
endif

;  open triangle
; ----------------
if pointselection eq 4 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
        usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
        symsel= 8
        symcolor= 0
endif

;  asterisk
; ----------
if pointselection eq 5 then begin
        symsel= 2
        symcolor= 200
endif

;  filled square
; --------------------------
if pointselection eq 6 then begin
        symsize= 1.0
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 8
        symcolor= 20
endif

;  diamond
; ----------
if pointselection eq 7 then begin
        symsel= 4
        symcolor= 120
endif

;  triangle (open)
; -----------------
if pointselection eq 8 then begin
        symsel= 5
        symcolor= 180
endif

;  plus sign
; ----------
if pointselection eq 9 then begin
        symsel= 1
        symcolor= 170
endif


	; read centers
	; -------------
	read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=frun+"/centers.txt"

	cp_time= transpose(cp_time)
	cp_cen1= transpose(cp_cen1)
	cp_cen2= transpose(cp_cen2)

x = cp_cen1[*,0] - cp_cen2[*,0]
y = cp_cen1[*,1] - cp_cen2[*,1]
z = cp_cen1[*,2] - cp_cen2[*,2]

rdiff = sqrt(x*x + y*y + z*z)

;rdiff= rdiff/0.7
;rdiff= rdiff/1000.0
cp_time= cp_time/0.7


if keyword_set(ylog) then begin
idx= where(rdiff lt 1.0e-6)
if idx(0) ne -1 then rdiff(idx)= 1.0e-6
endif

;oplot, time, rdiff_b, psym=-3, color= 50, thick=3.0, linestyle= 1
;oplot, time, rdiff_d, psym=-3, color= 100, thick=3.0, linestyle= 1
oplot, cp_time, rdiff, psym=-symsel, color= symcolor, thick=3.0

if keyword_set(msg) then begin
	xyouts, x0, y0, msg, /normal, size= 1.1, color= symcolor, charthick= 3.0
endif else begin
	xyouts, x0, y0, fload_getid(frun), /normal, size= 1.1, color= symcolor, charthick= 3.0
endelse




end





; ------------------------------------------------------------------


pro velocity_info, frun

	if not keyword_set(frun) then begin
		print, " "
		print, " velocity_info, frun  (need centers.txt and cenvel.txt) "
		print, " "
		return
	endif

	sendto= 'ps'

	; read centers
	; -------------
	read_centerpositions, center_time, center_1, center_2, filename=frun+"/centers.txt"

	center_time= transpose(center_time)
	center_1= transpose(center_1)
	center_2= transpose(center_2)

	center_time= center_time/0.7
	center_1= center_1/0.7
	center_2= center_2/0.7

	x = center_1[*,0] - center_2[*,0]
	y = center_1[*,1] - center_2[*,1]
	z = center_1[*,2] - center_2[*,2]

	rdiff = sqrt(x*x + y*y + z*z)



	; read com velocity
	; ------------------
	;read_centerpositions, comvel_time, comvel_1, comvel_2, filename=frun+"/comvel.txt"
	read_centerpositions, comvel_time, comvel_1, comvel_2, filename=frun+"/cenvel.txt"

	comvel_time= transpose(comvel_time)
	comvel_1= transpose(comvel_1)
	comvel_2= transpose(comvel_2)


	vx= comvel_1[*,0] - comvel_2[*,0]
	vy= comvel_1[*,1] - comvel_2[*,1]
	vz= comvel_1[*,2] - comvel_2[*,2]

	velrel= sqrt(vx*vx + vy*vy + vz*vz)

	vrad_x= vx*x/rdiff 
	vrad_y= vy*y/rdiff 
	vrad_z= vz*z/rdiff 

	vrad= sqrt(vrad_x*vrad_x + vrad_y*vrad_y + vrad_z*vrad_z)
	vtan = sqrt(velrel*velrel - vrad*vrad)

	;for i=0,n_elements(center_time)-1 do begin
	;	print, center_time[i], rdiff[i], velrel[i], vrad[i], vtan[i]
	;endfor


; ----------------------------------
;   Try this plot thingy
; ----------------------------------

;timemin= 0.0
;timemax= max(center_time)
timemin= 1.0
timemax= 2.0

;ymax= 500.0
ymax= 250.0
;ymax= 120.0

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename='velocity_info.eps', colortable= 4



!p.position= [0.18, 0.15, 0.98, 0.99]

plot, [1.0], [1.0], psym=-3, $
        xrange=[timemin,timemax], $
	yrange=[0,ymax], $
        xstyle=1, ystyle=1, $
        color= 0, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        ;xtitle="!6Time (h!E-1!N Gyr)", $
        xtitle="!6Time (Gyr)", $
        ytitle="!6Velocity (km sec!E-1!N)", $
        /nodata


oplot, center_time, velrel, psym=-3, color= 0, thick=5.0
oplot, center_time, vrad, psym=-3, color= 150, thick=5.0, linestyle= 1
oplot, center_time, vtan, psym=-3, color= 50, thick=5.0, linestyle= 2


xyouts, 0.62, 0.93, 'total relative', /normal, size= 1.2, color= 0, charthick= 3.0
xyouts, 0.65, 0.89, 'velocity', /normal, size= 1.2, color= 0, charthick= 3.0
xyouts, 0.62, 0.84, 'radial', /normal, size= 1.2, color= 150, charthick= 3.0
xyouts, 0.62, 0.79, 'tangential', /normal, size= 1.2, color= 50, charthick= 3.0


x= [5.4,5.4]
y=[0,300.0]
oplot, x, y, psym=-3, color= 0, thick=3.0, linestyle= 1

; -------------------------------

device, /close



end
