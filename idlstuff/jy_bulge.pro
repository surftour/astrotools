;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Determine Disk and Bulge Components
;   -----------------------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro do_one, junk

do_one_set, "/raid4/tcox/minor/min_0", 185
;do_one_set, "/raid4/tcox/minor/min_30", 200
;do_one_set, "/raid4/tcox/minor/min_30_fg0.0", 200
;do_one_set, "/raid4/tcox/minor/min_30_fg0.4", 200
;do_one_set, "/raid4/tcox/minor/min_150", 220
;do_one_set, "/raid4/tcox/minor/iso_Sb", 120

end



pro do_one_set, frun, snapnum

	jze, frun, snapnum
	ecirchist, frun, snapnum, /loadedsnap
	parts_4, frun, snapnum, /loadedsnap
	parts_6_cnt, frun, snapnum, /loadedsnap

end



;=====================================================================



pro jze, frun, snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "jze, frun, snapnum"
   print, "  "
   return
endif

filename=frun+'/jze.eps'
;frun= "/raid4/tcox/minor/min_30"
;snapnum=200
;frun= "/raid4/tcox/minor/iso_Sb"
;snapnum=120


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; specific bining energy
xaxistitle= "!6Specific Binding Energy (10!E5!N km!E2!N s!E-2!N)"
xmax =  2.3
;xmax =  0.3
xmin = -8.0
;xmin = -1.4

; angular momentum
yaxistitle= 'J!Dz!N (10!E3!N kpc km s!E-1!N)'
ymax =  50.0
;ymax =  10.0
;ymin = -10.0
ymin = -50.0



; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.95
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
	color= 0, xcharsize=1.50, ycharsize=1.50, $
	xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle


; -----------------------------------------------



ok= fload_snapshot_bh(frun,snapnum)




; -------------------------
;  All Stars
; -------------------------

print, " "
print, " ALL STARS "
print, " "

; specific binding energy
; -------------------------
etot= fload_allstars_energy(1, /specific)
etot=etot/1.0e+5
print, "energy  max/min=", max(etot), min(etot)

; angular momentum (z)
; -------------------------
jznew= return_allstars_jznew(1, n_hat=n_hat)
print, "n_hat= ", n_hat
n_hat_allstars=n_hat
jznew= jznew/1000.0
print, "jznew  max/min=", max(jznew), min(jznew)


oplot, etot, jznew, psym=3, color= 50


xyouts, 0.20, 0.32, 'all stars', size=1.2, color=50, /normal, charthick=2.0


; determine maximum j_z
bins= 200
e_min= min(e_tot) ;-2.0
e_max= max(e_tot) ;0.8
determine_max_jz, bins, e_min, e_max, etot, jznew, e_e=e_e, jz_e=jz_e
oplot, e_e, jz_e, psym=-3, color= 200, thick= 4.0

oplot, [-1.30,-1.20], [-6,-6], psym=-3, color=200, thick=2.0
xyouts, 0.20, 0.28, 'J!Dz!N(E)', size=1.2, color=200, /normal, charthick=2.0





; -------------------------
;  New Stars
; -------------------------

etot= fload_newstars_energy(1, /specific)
etot=etot/1.0e+5
print, "energy  max/min=", max(etot), min(etot)

jznew= return_newstars_jznew(1, n_hat=n_hat_allstars)
jznew= jznew/1000.0
print, "jznew  max/min=", max(jznew), min(jznew)


oplot, etot, jznew, psym=3, color= 150

xyouts, 0.20, 0.24, 'new stars', size=1.2, color=150, /normal, charthick=2.0





; -------------------------
;  Gas
; -------------------------

etot= fload_gas_energy(1, /specific)
etot=etot/1.0e+5
print, "gas energy  max/min=", max(etot), min(etot)

jznew= return_gas_jznew(1, n_hat=n_hat_allstars)
jznew= jznew/1000.0
print, "gas jznew  max/min=", max(jznew), min(jznew)


oplot, etot, jznew, psym=3, color= 100

xyouts, 0.20, 0.20, 'gas', size=1.2, color=100, /normal, charthick=2.0







; print extras
; -------------
xyouts, 0.20, 0.40, frun, size=1.2, color=0, /normal, charthick=4.0
xyouts, 0.20, 0.36, fload_timelbl(0.7,2), size=1.2, color=0, /normal, charthick=4.0




; done
; -------------
device, /close


end







;====================================================================================







pro ecirchist, frun, snapnum, loadedsnap=loadedsnap

if not keyword_set(frun) then begin
   print, "  "
   print, "ecirchist, frun, snapnum"
   print, "  "
   return
endif

filename=frun+'/ecirchist.eps'
;frun= "/raid4/tcox/minor/min_30"
;snapnum=200
;frun= "/raid4/tcox/minor/iso_Sb"
;snapnum=120


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; radius
xaxistitle= "!6J!Dz!N / J!Dcirc!N"
xmax =  1.2
xmin = -1.0

; number (histogram)
yaxistitle= ' '
ymax = 1.2
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.06
y0= 0.15
x1= 0.98
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
	color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'


; -----------------------------------------------



if not keyword_set(loadedsnap) then ok= fload_snapshot_bh(frun,snapnum)



; -------------------------
;  All Stars
; -------------------------


; specific binding energy
etot= fload_allstars_energy(1, /specific)
etot=etot/1.0e+5
print, "energy  max/min=", max(etot), min(etot)

; angular momentum (z)
jznew= return_allstars_jznew(1, n_hat=n_hat)
print, "n_hat= ", n_hat
n_hat_allstars=n_hat
jznew= jznew/1000.0
print, "jznew  max/min=", max(jznew), min(jznew)

; determine maximum j_z
bins= 200
;e_min= -2.0
;e_max= 0.8
e_min= min(e_tot)
e_max= max(e_tot)
determine_max_jz, bins, e_min, e_max, etot, jznew, e_e=e_e, jz_e=jz_e



;---------------------------------

; determine the e_j = j_z / j_circ(E)

;map the particle's energy to an above bin
e_bin= bins * (etot - e_min) / (e_max - e_min)


e_idx= long(e_bin)
jz_circ_e= jz_e(e_idx)
;idx=where(jz_circ_e eq 0)
;if idx(0) ne -1 then jz_circ_e(idx)= max(jz_circ_e)

e_j= jznew / jz_circ_e


; plot histogram
; ------------------
temp= process_histogram(e_j, xmax=xmax, xmin=xmin, levels=100, oplotit=50, normalization=normalization)
orig_normalization=normalization
print, min(temp), max(temp)


xyouts, 0.10, 0.82, 'all stars', size=1.2, color=50, /normal, charthick=2.0


   
; -------------------------
;  New Stars
; -------------------------

etot= fload_newstars_energy(1, /specific)
etot=etot/1.0e+5
print, "ns energy  max/min=", max(etot), min(etot)

jznew= return_newstars_jznew(1, n_hat=n_hat_allstars)
jznew= jznew/1000.0
print, "ns jznew  max/min=", max(jznew), min(jznew)

e_bin= bins * (etot - e_min) / (e_max - e_min)
e_idx= long(e_bin)
jz_circ_e= jz_e(e_idx)
e_j= jznew / jz_circ_e

; plot histogram
; ------------------
temp= process_histogram(e_j, xmax=xmax, xmin=xmin, levels=100, oplotit=150, mannorm=orig_normalization)
print, min(temp), max(temp)


xyouts, 0.10, 0.78, 'new stars', size=1.2, color=150, /normal, charthick=2.0



; -------------------------
;  Gas
; -------------------------

etot= fload_gas_energy(1, /specific)
etot=etot/1.0e+5
print, "gas energy  max/min=", max(etot), min(etot)

jznew= return_gas_jznew(1, n_hat=n_hat_allstars)
jznew= jznew/1000.0
print, "gas jznew  max/min=", max(jznew), min(jznew)

e_bin= bins * (etot - e_min) / (e_max - e_min)
e_idx= long(e_bin)
jz_circ_e= jz_e(e_idx)
e_j= jznew / jz_circ_e


; plot histogram
; ------------------
temp= process_histogram(e_j, xmax=xmax, xmin=xmin, levels=100, oplotit=100, mannorm=orig_normalization)
print, min(temp), max(temp)

xyouts, 0.10, 0.74, 'gas', size=1.2, color=100, /normal, charthick=2.0





;---------------------------------
; another method to determin
; the circular j is via the mass
; interior to radius.
;
; make binned mass data
; m_r= fltarr(bins)
;
; circular v for each particle
;mtot= fload_all_mass(1)
;vcirc= 0.
;
; --------------------------------





; print extras
; -------------
xyouts, 0.10, 0.90, frun, size=1.2, color=0, /normal, charthick=4.0
xyouts, 0.10, 0.86, fload_timelbl(0.7,2), size=1.2, color=0, /normal, charthick=4.0




; done
; -----
device, /close


end





;====================================================================================
;
;
;         ---------------------------------------
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |    all stars     |    thin disk     |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         ---------------------------------------
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |    other         |    gas           |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         ---------------------------------------
;
;
;


pro parts_4, frun, snapnum, loadedsnap=loadedsnap

if not keyword_set(frun) then begin
   print, "  "
   print, "parts_4, frun, snapnum"
   print, "  "
   return
endif

filename=frun+'/components.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=20, newysize=20



; -----------------------------
; set up constants
; -----------------------------

xlen=  50.0

xmax = +xlen
xmin = -xlen
ymax = +xlen
ymin = -xlen


x0= 0.02
x1= 0.50
x2= 0.98
y0= 0.02
y1= 0.50
y2= 0.98




if not keyword_set(loadedsnap) then ok= fload_snapshot_bh(frun,snapnum)





;------------------------------------------------------------------
;
;  Now split things up by component
;
;


; -------------------------
;  All Stars
; -------------------------

; x and z are already fixed to all star quantities

; specific binding energy
etot= fload_allstars_energy(1, /specific)
etot=etot/1.0e+5
print, "energy  max/min=", max(etot), min(etot)

; angular momentum (z)
jznew= return_allstars_jznew(1, n_hat=n_hat)
print, "n_hat= ", n_hat
n_hat_allstars=n_hat
jznew= jznew/1000.0
print, "jznew  max/min=", max(jznew), min(jznew)

; determine maximum j_z
bins= 200
;e_min= -2.0
;e_max= 0.8
e_min= min(e_tot)
e_max= max(e_tot)
determine_max_jz, bins, e_min, e_max, etot, jznew, e_e=e_e, jz_e=jz_e

; determine the e_j = j_z / j_circ(E)
;map the particle's energy to an above bin
e_bin= bins * (etot - e_min) / (e_max - e_min)
e_idx= long(e_bin)
jz_circ_e= jz_e(e_idx)
e_j= jznew / jz_circ_e


x= fload_allstars_xyz('x')
y= fload_allstars_xyz('y')
z= fload_allstars_xyz('z')

; now rotate so we're edgeon
rotate_theta= acos(n_hat[2])*180./!PI
;rotate_theta= atan(sqrt(n_hat[0]*n_hat[0] + n_hat[1]*n_hat[1])/n_hat[2])*180./!PI
rotate_phi= atan(n_hat[1]/n_hat[0])*180./!PI
process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

x_all= x_new
y_all= y_new
z_all= z_new

idx= where(e_j ge 0.85)
x_thin= x_new(idx)
z_thin= z_new(idx)

idx= where(e_j lt 0.85)
x_other= x_new(idx)
z_other= z_new(idx)



; ----------
;  1111111  - top left
; ----------

!p.position= [x0, y1, x1, y2]
!p.ticklen=0.03

plot, [0], [0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, /normal, $
        xtickformat='(a1)', ytickformat='(a1)', xcharsize=0.01, ycharsize=0.01

oplot, x_all, z_all, psym=3, color= 0

xyouts, 0.35, 0.94, 'all stars', size=1.2, color=0, /normal, charthick=2.0

; print extras
; -------------
xyouts, 0.05, 0.56, frun, size=1.2, color=0, /normal, charthick=4.0
xyouts, 0.05, 0.52, fload_timelbl(0.7,2), size=1.2, color=0, /normal, charthick=4.0


; ----------
;  2222222
; ----------

!p.position= [x1, y1, x2, y2]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

oplot, x_thin, z_thin, psym=3, color= 50

xyouts, 0.85, 0.94, 'thin disk', size=1.2, color=50, /normal, charthick=2.0


; ----------
;  3333333
; ----------

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

oplot, x_other, z_other, psym=3, color= 150

xyouts, 0.35, 0.46, 'other', size=1.2, color=150, /normal, charthick=2.0


; ----------
;  4444444
; ----------

!p.position= [x1, y0, x2, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

x= fload_gas_xyz('x')
y= fload_gas_xyz('y')
z= fload_gas_xyz('z')

process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

x= x_new
y= y_new
z= z_new

oplot, x, z, psym=3, color= 100

xyouts, 0.85, 0.46, 'gas', size=1.2, color=100, /normal, charthick=2.0


; done
; -----
device, /close



end






;====================================================================================
;   PLOT FORMAT
;
;
;         ---------------------------------------
;         |                  |                  |
;         |                  |                  |
;         |                  |     primary      |
;         |    all stars     |    old stars     |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         ---------------------------------------
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |    primary       |    satellite     |
;         |   new stars      |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         ---------------------------------------
;
;


pro parts_4_cnt, frun, snapnum, loadedsnap=loadedsnap

if not keyword_set(frun) then begin
   print, "  "
   print, "parts_6_cnt, frun, snapnum"
   print, "  "
   return
endif

filename=frun+'/components4_v2.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=20, newysize=20



; -----------------------------
; set up constants
; -----------------------------

xlen=  50.0 

xmax = +xlen
xmin = -xlen
ymax = +xlen
ymin = -xlen


x0= 0.02
x1= 0.50
x2= 0.98
y0= 0.02
y1= 0.50 
y2= 0.98 





;set_maxden= 10.0
set_maxden= 1.0    ; 10^4 msolar/pc2
;set_dynrng= 1.0e+5
;set_dynrng= 1.0e+6
;set_maxden= 0.5e-1
set_dynrng= 1.0e+6




if not keyword_set(loadedsnap) then ok= fload_snapshot_bh(frun,snapnum)



;------------------------------------------------------------------
;
;  Now split things up by component
;
;


; -------------------------
;  All Stars
; -------------------------

; x and z are already fixed to all star quantities

; specific binding energy
etot= fload_allstars_energy(1, /specific)
etot=etot/1.0e+5
print, "energy  max/min=", max(etot), min(etot)

; angular momentum (z)
jznew= return_allstars_jznew(1, n_hat=n_hat)
print, "n_hat= ", n_hat
n_hat_allstars=n_hat
jznew= jznew/1000.0
print, "jznew  max/min=", max(jznew), min(jznew)

; determine maximum j_z
bins= 200
;e_min= -2.0
;e_max= 0.8
e_min= min(etot)
e_max= max(etot)
determine_max_jz, bins, e_min, e_max, etot, jznew, e_e=e_e, jz_e=jz_e

; determine the e_j = j_z / j_circ(E)
;map the particle's energy to an above bin
e_bin= bins * (etot - e_min) / (e_max - e_min)
e_idx= long(e_bin)
jz_circ_e= jz_e(e_idx)
e_j= jznew / jz_circ_e


x= fload_allstars_xyz('x')
y= fload_allstars_xyz('y')
z= fload_allstars_xyz('z')
m= fload_allstars_mass(1)
;hsml= fload_allstars_hsml(1)
hsml= m*0.0 + 0.2

; now rotate so we're edgeon
rotate_theta= acos(n_hat[2])*180./!PI
;rotate_theta= atan(sqrt(n_hat[0]*n_hat[0] + n_hat[1]*n_hat[1])/n_hat[2])*180./!PI
rotate_phi= atan(n_hat[1]/n_hat[0])*180./!PI
process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new


x_all= x_new
y_all= -y_new
z_all= z_new
m_all= m

;---------------------------------------------------------

; high resolution runs
startid_primary= 1
npart_primary= 1600001
startid_sat= 1600002
npart_sat= 170001


;---------------------------------------------------------

;x= fload_1gal_oldstars_xyz('x', startid_primary, npart_primary)
;y= fload_1gal_oldstars_xyz('y', startid_primary, npart_primary)
;z= fload_1gal_oldstars_xyz('z', startid_primary, npart_primary)
;m= fload_1gal_oldstars_mass(startid_primary, npart_primary)
xo= fload_1gal_oldstars_xyz('x', startid_primary, npart_primary)
yo= fload_1gal_oldstars_xyz('y', startid_primary, npart_primary)
zo= fload_1gal_oldstars_xyz('z', startid_primary, npart_primary)
mo= fload_1gal_oldstars_mass(startid_primary, npart_primary)
;
;process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new
;
;x_pri_old= x_new
;y_pri_old= -y_new
;z_pri_old= z_new
;m_pri_old= m
;hsml_pri_old= m*0.0 + 0.2

;---------------------------------------------------------

;x= fload_1gal_newstars_xyz('x', startid_primary, npart_primary)
;y= fload_1gal_newstars_xyz('y', startid_primary, npart_primary)
;z= fload_1gal_newstars_xyz('z', startid_primary, npart_primary)
;m= fload_1gal_newstars_mass(startid_primary, npart_primary)
xn= fload_1gal_newstars_xyz('x', startid_primary, npart_primary)
yn= fload_1gal_newstars_xyz('y', startid_primary, npart_primary)
zn= fload_1gal_newstars_xyz('z', startid_primary, npart_primary)
mn= fload_1gal_newstars_mass(startid_primary, npart_primary)
an= fload_1gal_newstars_age(startid_primary, npart_primary)
stop
;

idx= where(an lt 2.9)
x= [xo, xn(idx)]
y= [yo, yn(idx)]
z= [zo, zn(idx)]
m= [mo, mn(idx)]

process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

x_pri_old= x_new
y_pri_old= -y_new
z_pri_old= z_new
m_pri_old= m
hsml_pri_old= m*0.0 + 0.2

idx= where(an ge 2.9)
x= xn(idx)
y= yn(idx)
z= zn(idx)
m= mn(idx)

process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

x_pri_new= x_new
y_pri_new= -y_new
z_pri_new= z_new
m_pri_new= m
hsml_pri_new= m*0.0 + 0.2

;---------------------------------------------------------

x= fload_1gal_allstars_xyz('x', startid_sat, npart_sat)
y= fload_1gal_allstars_xyz('y', startid_sat, npart_sat)
z= fload_1gal_allstars_xyz('z', startid_sat, npart_sat)
m= fload_1gal_allstars_mass(startid_sat, npart_sat)

process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

x_sat= x_new
y_sat= -y_new
z_sat= z_new
m_sat= m
hsml_sat= m*0.0 + 0.2

;---------------------------------------------------------





; ----------
;  1111111  - top left
; ----------

contour_makepic, x_all, z_all, y_all, m_all, xlen, xz=xz, yz=yz, hsml=hsml, $
        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

tv, NxNImage, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal

!p.position= [x0, y1, x1, y2]
!p.ticklen=0.03

plot, [0], [0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, /normal, $
        xtickformat='(a1)', ytickformat='(a1)', xcharsize=0.01, ycharsize=0.01

;oplot, x_all, z_all, psym=3, color= 0

xyouts, 0.35, 0.92, 'all stars', size=1.4, color=0, /normal, charthick=3.0

; print extras
; -------------
;xyouts, 0.05, 0.68, frun, size=1.2, color=0, /normal, charthick=2.0
;xyouts, 0.05, 0.705, fload_timelbl(0.7,2), size=1.2, color=0, /normal, charthick=2.0


; ----------
;  2222222
; ----------

contour_makepic, x_pri_old, z_pri_old, y_pri_old, m_pri_old, xlen, xz=xz, yz=yz, hsml=hsml_pri_old, $
        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

tv, NxNImage, x1, y1, xsize=(x2-x1), ysize=(y2-y1), /normal

!p.position= [x1, y1, x2, y2]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

;oplot, x_thin, z_thin, psym=3, color= 50

xyouts, 0.66, 0.92, 'primary: old stars', size=1.4, color=0, /normal, charthick=3.0


; ----------
;  3333333
; ----------

contour_makepic, x_pri_new, z_pri_new, y_pri_new, m_pri_new, xlen, xz=xz, yz=yz, hsml=hsml_pri_new, $
        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

;oplot, x_other, z_other, psym=3, color= 120

xyouts, 0.21, 0.44, 'primary: post-merger stars', size=1.4, color=0, /normal, charthick=3.0


; ----------
;  4444444
; ----------

contour_makepic, x_sat, z_sat, y_sat, m_sat, xlen, xz=xz, yz=yz, hsml=hsml_sat, $
        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

tv, NxNImage, x1, y0, xsize=(x2-x1), ysize=(y1-y0), /normal


!p.position= [x1, y0, x2, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

;oplot, x_nonrot, z_nonrot, psym=3, color= 150

xyouts, 0.77, 0.44, 'satellite', size=1.4, color=0, /normal, charthick=3.0



; done
; -----
device, /close



end









;====================================================================================
;   PLOT FORMAT
;
;
;      -----------------------------------------------------------------------------
;      |                  |                  |                  |                  |
;      |                  |                  |                  |                  |
;      |                  |                  |   primary:       |                  |
;      |    initial disk  |    all stars     |   old stars      |    satellite     |
;      |                  |                  |                  |                  |
;      |                  |                  |                  |                  |
;      |                  |                  |                  |                  |
;      |                  |                  |                  |                  |
;      |                  |                  |                  |                  |
;      -----------------------------------------------------------------------------
;
;


pro parts_1x4_cnt, frun, snapnum, loadedsnap=loadedsnap

if not keyword_set(frun) then begin
   print, "  "
   print, "parts_1x4_cnt, frun, snapnum"
   print, "  "
   return
endif

filename=frun+'/components1x4.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=32, newysize=8



; -----------------------------
; set up constants
; -----------------------------

xlen=  50.0 

xmax = +xlen
xmin = -xlen
ymax = +xlen
ymin = -xlen


x0= 0.006667
x1= 0.253333
x2= 0.50
x3= 0.746667
x4= 0.993333

y0= 0.02
y1= 0.98 





;set_maxden= 10.0
set_maxden= 1.0    ; 10^4 msolar/pc2
;set_dynrng= 1.0e+5
;set_dynrng= 1.0e+6
;set_maxden= 0.5e-1
set_dynrng= 1.0e+6







;------------------------------------------------------------------
;
;  Now split things up by component
;
;


ok= fload_snapshot_bh("data/minor/iso_Sb",10)

x_init= fload_allstars_xyz('x')
y_init= fload_allstars_xyz('y')
z_init= fload_allstars_xyz('z')
m_init= fload_allstars_mass(1)
;hsml= fload_allstars_hsml(1)
hsml_init= m_init*0.0 + 0.2






; -------------------------
;  All Stars
; -------------------------

if not keyword_set(loadedsnap) then ok= fload_snapshot_bh(frun,snapnum)

; x and z are already fixed to all star quantities

; specific binding energy
etot= fload_allstars_energy(1, /specific)
etot=etot/1.0e+5
print, "energy  max/min=", max(etot), min(etot)

; angular momentum (z)
jznew= return_allstars_jznew(1, n_hat=n_hat)
print, "n_hat= ", n_hat
n_hat_allstars=n_hat
jznew= jznew/1000.0
print, "jznew  max/min=", max(jznew), min(jznew)

; determine maximum j_z
bins= 200
;e_min= -2.0
;e_max= 0.8
e_min= min(etot)
e_max= max(etot)
determine_max_jz, bins, e_min, e_max, etot, jznew, e_e=e_e, jz_e=jz_e

; determine the e_j = j_z / j_circ(E)
;map the particle's energy to an above bin
e_bin= bins * (etot - e_min) / (e_max - e_min)
e_idx= long(e_bin)
jz_circ_e= jz_e(e_idx)
e_j= jznew / jz_circ_e


x= fload_allstars_xyz('x')
y= fload_allstars_xyz('y')
z= fload_allstars_xyz('z')
m= fload_allstars_mass(1)
;hsml= fload_allstars_hsml(1)
hsml= m*0.0 + 0.2

; now rotate so we're edgeon
rotate_theta= acos(n_hat[2])*180./!PI
;rotate_phi= atan(n_hat[0],n_hat[1])*180./!PI
rotate_phi= atan(n_hat[1],n_hat[0])*180./!PI

;
;
;  needs 2.5 rotation
;
;
process_rotation, n_hat[0], n_hat[1], n_hat[2], rotate_theta, rotate_phi, xx, yy, zz
print, "n_hat= ", n_hat
print, "xx, yy, zz= ", xx, yy, zz
;stop
process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new


x_all= x_new
y_all= -y_new
z_all= z_new
m_all= m

;---------------------------------------------------------

; high resolution runs
startid_primary= 1
npart_primary= 1600001
startid_sat= 1600002
npart_sat= 170001


;---------------------------------------------------------

;x= fload_1gal_oldstars_xyz('x', startid_primary, npart_primary)
;y= fload_1gal_oldstars_xyz('y', startid_primary, npart_primary)
;z= fload_1gal_oldstars_xyz('z', startid_primary, npart_primary)
;m= fload_1gal_oldstars_mass(startid_primary, npart_primary)
xo= fload_1gal_oldstars_xyz('x', startid_primary, npart_primary)
yo= fload_1gal_oldstars_xyz('y', startid_primary, npart_primary)
zo= fload_1gal_oldstars_xyz('z', startid_primary, npart_primary)
mo= fload_1gal_oldstars_mass(startid_primary, npart_primary)
;
;process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new
;
;x_pri_old= x_new
;y_pri_old= -y_new
;z_pri_old= z_new
;m_pri_old= m
;hsml_pri_old= m*0.0 + 0.2

;---------------------------------------------------------

;x= fload_1gal_newstars_xyz('x', startid_primary, npart_primary)
;y= fload_1gal_newstars_xyz('y', startid_primary, npart_primary)
;z= fload_1gal_newstars_xyz('z', startid_primary, npart_primary)
;m= fload_1gal_newstars_mass(startid_primary, npart_primary)
xn= fload_1gal_newstars_xyz('x', startid_primary, npart_primary)
yn= fload_1gal_newstars_xyz('y', startid_primary, npart_primary)
zn= fload_1gal_newstars_xyz('z', startid_primary, npart_primary)
mn= fload_1gal_newstars_mass(startid_primary, npart_primary)
an= fload_1gal_newstars_age(startid_primary, npart_primary)
;stop
;

idx= where(an lt 2.9)
x= [xo, xn(idx)]
y= [yo, yn(idx)]
z= [zo, zn(idx)]
m= [mo, mn(idx)]

process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

x_pri_old= x_new
y_pri_old= -y_new
z_pri_old= z_new
m_pri_old= m
hsml_pri_old= m*0.0 + 0.2

idx= where(an ge 2.9)
x= xn(idx)
y= yn(idx)
z= zn(idx)
m= mn(idx)

process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

x_pri_new= x_new
y_pri_new= -y_new
z_pri_new= z_new
m_pri_new= m
hsml_pri_new= m*0.0 + 0.2

;---------------------------------------------------------

x= fload_1gal_allstars_xyz('x', startid_sat, npart_sat)
y= fload_1gal_allstars_xyz('y', startid_sat, npart_sat)
z= fload_1gal_allstars_xyz('z', startid_sat, npart_sat)
m= fload_1gal_allstars_mass(startid_sat, npart_sat)

process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

x_sat= x_new
y_sat= -y_new
z_sat= z_new
m_sat= m
hsml_sat= m*0.0 + 0.2

;---------------------------------------------------------




; ----------
;  1111111
; ----------

contour_makepic, x_init, z_init, y_init, m_init, xlen, xz=xz, yz=yz, hsml=hsml_init, $
        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

!p.position= [x0, y0, x1, y1]
!p.ticklen=0.03

plot, [0], [0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, /normal, $
        xtickformat='(a1)', ytickformat='(a1)', xcharsize=0.01, ycharsize=0.01

xyouts, 0.02, 0.87, '!6primary evolved in isolation', size=1.4, color=0, /normal, charthick=3.0





; ----------
;  2222222 
; ----------

contour_makepic, x_all, z_all, y_all, m_all, xlen, xz=xz, yz=yz, hsml=hsml, $
        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

tv, NxNImage, x1, y0, xsize=(x2-x1), ysize=(y1-y0), /normal

!p.position= [x1, y0, x2, y1]
!p.ticklen=0.03

plot, [0], [0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, /normal, $
        xtickformat='(a1)', ytickformat='(a1)', xcharsize=0.01, ycharsize=0.01

xyouts, 0.29, 0.87, '1:8 merger remnant', size=1.4, color=0, /normal, charthick=3.0
xyouts, 0.40, 0.10, '(all stars)', size=1.4, color=0, /normal, charthick=3.0






; ----------
;  3333333
; ----------

contour_makepic, x_pri_new, z_pri_new, y_pri_new, m_pri_new, xlen, xz=xz, yz=yz, hsml=hsml_pri_new, $
        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

tv, NxNImage, x2, y0, xsize=(x3-x2), ysize=(y1-y0), /normal


!p.position= [x2, y0, x3, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

;xyouts, 0.29, 0.87, '1:8 merger remnant', size=1.4, color=0, /normal, charthick=3.0
xyouts, 0.56, 0.10, '(post-merger stars)', size=1.4, color=0, /normal, charthick=3.0






; ----------
;  4444444
; ----------

contour_makepic, x_sat, z_sat, y_sat, m_sat, xlen, xz=xz, yz=yz, hsml=hsml_sat, $
        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

tv, NxNImage, x3, y0, xsize=(x4-x3), ysize=(y1-y0), /normal

!p.position= [x3, y0, x4, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

xyouts, 0.90, 0.10, '(satellite)', size=1.4, color=0, /normal, charthick=3.0




; done
; -----
device, /close



end
















;====================================================================================
;
;
;         ---------------------------------------
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |    all stars     |    thin disk     |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         ---------------------------------------
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |    mod disk      |    non-rot       |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         ---------------------------------------
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |    gas disk      |    gas halo      |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         ---------------------------------------
;
;
;


pro parts_6, frun, snapnum, loadedsnap=loadedsnap

if not keyword_set(frun) then begin
   print, "  "
   print, "parts_6, frun, snapnum"
   print, "  "
   return
endif

filename=frun+'/components6.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=20, newysize=30



; -----------------------------
; set up constants
; -----------------------------

xlen=  50.0

xmax = +xlen
xmin = -xlen
ymax = +xlen
ymin = -xlen


x0= 0.02
x1= 0.50
x2= 0.98

y0= 0.01
y1= 0.336667
y2= 0.663333
y3= 0.99




if not keyword_set(loadedsnap) then ok= fload_snapshot_bh(frun,snapnum)





;------------------------------------------------------------------
;
;  Now split things up by component
;
;


; -------------------------
;  All Stars
; -------------------------

; x and z are already fixed to all star quantities

; specific binding energy
etot= fload_allstars_energy(1, /specific)
etot=etot/1.0e+5
print, "energy  max/min=", max(etot), min(etot)

; angular momentum (z)
jznew= return_allstars_jznew(1, n_hat=n_hat)
print, "n_hat= ", n_hat
n_hat_allstars=n_hat
jznew= jznew/1000.0
print, "jznew  max/min=", max(jznew), min(jznew)

; determine maximum j_z
bins= 200
;e_min= -2.0
;e_max= 0.8
e_min= min(e_tot)
e_max= max(e_tot)
determine_max_jz, bins, e_min, e_max, etot, jznew, e_e=e_e, jz_e=jz_e

; determine the e_j = j_z / j_circ(E)
;map the particle's energy to an above bin
e_bin= bins * (etot - e_min) / (e_max - e_min)
e_idx= long(e_bin)
jz_circ_e= jz_e(e_idx)
e_j= jznew / jz_circ_e


x= fload_allstars_xyz('x')
y= fload_allstars_xyz('y')
z= fload_allstars_xyz('z')

; now rotate so we're edgeon
rotate_theta= acos(n_hat[2])*180./!PI
;rotate_theta= atan(sqrt(n_hat[0]*n_hat[0] + n_hat[1]*n_hat[1])/n_hat[2])*180./!PI
rotate_phi= atan(n_hat[1]/n_hat[0])*180./!PI
process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

x_all= x_new
y_all= y_new
z_all= z_new

idx= where(e_j ge 0.85)
x_thin= x_new(idx)
z_thin= z_new(idx)

idx= where((e_j lt 0.85) and (e_j gt 0.25))
x_other= x_new(idx)
z_other= z_new(idx)

idx= where(e_j le 0.25)
x_nonrot= x_new(idx)
z_nonrot= z_new(idx)



; ----------
;  1111111  - top left
; ----------

!p.position= [x0, y2, x1, y3]
!p.ticklen=0.03

plot, [0], [0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, /normal, $
        xtickformat='(a1)', ytickformat='(a1)', xcharsize=0.01, ycharsize=0.01

oplot, x_all, z_all, psym=3, color= 0

xyouts, 0.35, 0.95, 'all stars', size=1.2, color=0, /normal, charthick=2.0

; print extras
; -------------
xyouts, 0.05, 0.72, frun, size=1.2, color=0, /normal, charthick=4.0
xyouts, 0.05, 0.68, fload_timelbl(0.7,2), size=1.2, color=0, /normal, charthick=4.0


; ----------
;  2222222
; ----------

!p.position= [x1, y2, x2, y3]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

oplot, x_thin, z_thin, psym=3, color= 50

xyouts, 0.85, 0.95, 'thin disk', size=1.2, color=50, /normal, charthick=2.0


; ----------
;  3333333
; ----------

!p.position= [x0, y1, x1, y2]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

oplot, x_other, z_other, psym=3, color= 120

xyouts, 0.35, 0.60, 'other', size=1.2, color=120, /normal, charthick=2.0


; ----------
;  4444444
; ----------

!p.position= [x1, y1, x2, y2]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

oplot, x_nonrot, z_nonrot, psym=3, color= 150

xyouts, 0.85, 0.60, 'non-rotating', size=1.2, color=150, /normal, charthick=2.0


; ----------
;  5555555
; ----------

etot= fload_gas_energy(1, /specific)
etot=etot/1.0e+5
print, "gas energy  max/min=", max(etot), min(etot)

jznew= return_gas_jznew(1, n_hat=n_hat_allstars)
jznew= jznew/1000.0
print, "gas jznew  max/min=", max(jznew), min(jznew)

e_bin= bins * (etot - e_min) / (e_max - e_min)
e_idx= long(e_bin)
jz_circ_e= jz_e(e_idx)
e_j= jznew / jz_circ_e


x= fload_gas_xyz('x')
y= fload_gas_xyz('y')
z= fload_gas_xyz('z')

process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

idx= where(e_j ge 0.85)
x_disk= x_new(idx)
z_disk= z_new(idx)

idx= where(e_j lt 0.85)
x_halo= x_new(idx)
z_halo= z_new(idx)



!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

oplot, x_disk, z_disk, psym=3, color= 100

xyouts, 0.85, 0.46, 'gas: thin disk', size=1.2, color=100, /normal, charthick=2.0


; ----------
;  5555555
; ----------

!p.position= [x1, y0, x2, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

oplot, x_halo, z_halo, psym=3, color= 100

xyouts, 0.85, 0.46, 'gas: thin disk', size=1.2, color=100, /normal, charthick=2.0


; done
; -----
device, /close



end





;====================================================================================
;   PLOT FORMAT
;
;  contour version of the above
;
;         ---------------------------------------
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |    all stars     |    thin disk     |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         ---------------------------------------
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |    mod disk      |    non-rot       |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         ---------------------------------------
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |    gas disk      |    gas halo      |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         |                  |                  |
;         ---------------------------------------
;
;
;


pro parts_6_cnt, frun, snapnum, loadedsnap=loadedsnap, $
				msg1=msg1, msg2=msg2

if not keyword_set(frun) then begin
   print, "  "
   print, "parts_6_cnt, frun [, snapnum]"
   print, "  "
   return
endif

if not keyword_set(snapnum) then begin
   print, "  "
   print, "seting snapnum"

   relaxationtime= 0.4

   spawn, "/bin/ls "+frun+"/snap* | wc ",result
   totnumsnap=long(result[0])-1
   snapnum= totnumsnap

   bhmtime= fload_bh_mergertime(frun)

   for i=0,totnumsnap do begin
            ok= fload_snapshot_bh(frun,i,/header_only)

            ; grab snapnum when system is "relaxed"
            if fload_time(1) gt (bhmtime+relaxationtime) then begin
                snapnum= i
                break
            endif
   endfor

   print, "snapnum= ", snapnum
   print, "  "
endif

filename=frun+'/components6.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=20, newysize=30



; -----------------------------
; set up constants
; -----------------------------

xlen=  50.0

xmax = +xlen
xmin = -xlen
ymax = +xlen
ymin = -xlen


x0= 0.02
x1= 0.50
x2= 0.98

y0= 0.01
y1= 0.336667
y2= 0.663333
y3= 0.99


;set_maxden= 10.0
set_maxden= 1.0    ; 10^4 msolar/pc2
;set_dynrng= 1.0e+5
;set_dynrng= 1.0e+6
;set_maxden= 0.5e-1
set_dynrng= 1.0e+6




if not keyword_set(loadedsnap) then ok= fload_snapshot_bh(frun,snapnum)



;------------------------------------------------------------------
;
;  Now split things up by component
;
;


; -------------------------
;  All Stars
; -------------------------

; x and z are already fixed to all star quantities

; specific binding energy
etot= fload_allstars_energy(1, /specific)
etot=etot/1.0e+5
print, "energy  max/min=", max(etot), min(etot)

; angular momentum (z)
jznew= return_allstars_jznew(1, n_hat=n_hat)
print, "n_hat= ", n_hat
n_hat_allstars=n_hat
jznew= jznew/1000.0
print, "jznew  max/min=", max(jznew), min(jznew)

; determine maximum j_z
bins= 200
;e_min= -2.0
;e_max= 0.8
e_min= min(e_tot)
e_max= max(e_tot)
determine_max_jz, bins, e_min, e_max, etot, jznew, e_e=e_e, jz_e=jz_e

; determine the e_j = j_z / j_circ(E)
;map the particle's energy to an above bin
e_bin= bins * (etot - e_min) / (e_max - e_min)
e_idx= long(e_bin)
jz_circ_e= jz_e(e_idx)
e_j= jznew / jz_circ_e


x= fload_allstars_xyz('x')
y= fload_allstars_xyz('y')
z= fload_allstars_xyz('z')
m= fload_allstars_mass(1)
;hsml= fload_allstars_hsml(1)
hsml= m*0.0 + 0.2

; now rotate so we're edgeon
rotate_theta= acos(n_hat[2])*180./!PI
;rotate_theta= atan(sqrt(n_hat[0]*n_hat[0] + n_hat[1]*n_hat[1])/n_hat[2])*180./!PI
rotate_phi= atan(n_hat[1]/n_hat[0])*180./!PI
process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

x_all= x_new
y_all= -y_new
z_all= z_new

idx= where(e_j ge 0.85)
x_thin= x_new(idx)
y_thin= y_new(idx)
z_thin= z_new(idx)
m_thin= m(idx)
hsml_thin= hsml(idx)

idx= where((e_j lt 0.85) and (e_j gt 0.25))
x_other= x_new(idx)
y_other= y_new(idx)
z_other= z_new(idx)
m_other= m(idx)
hsml_other= hsml(idx)

idx= where(e_j le 0.25)
x_nonrot= x_new(idx)
y_nonrot= y_new(idx)
z_nonrot= z_new(idx)
m_nonrot= m(idx)
hsml_nonrot= hsml(idx)



; ----------
;  1111111  - top left
; ----------

contour_makepic, x_all, z_all, y_all, m, xlen, xz=xz, yz=yz, hsml=hsml, $
        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

tv, NxNImage, x0, y2, xsize=(x1-x0), ysize=(y3-y2), /normal

!p.position= [x0, y2, x1, y3]
!p.ticklen=0.03

plot, [0], [0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, /normal, $
        xtickformat='(a1)', ytickformat='(a1)', xcharsize=0.01, ycharsize=0.01

;oplot, x_all, z_all, psym=3, color= 0

xyouts, 0.35, 0.95, 'all stars', size=1.4, color=0, /normal, charthick=3.0

; print extras
; -------------
xyouts, 0.05, 0.68, frun, size=1.2, color=0, /normal, charthick=2.0
xyouts, 0.05, 0.705, fload_timelbl(0.7,2), size=1.2, color=0, /normal, charthick=2.0


; ----------
;  2222222
; ----------

contour_makepic, x_thin, z_thin, y_thin, m_thin, xlen, xz=xz, yz=yz, hsml=hsml_thin, $
        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

tv, NxNImage, x1, y2, xsize=(x2-x1), ysize=(y3-y2), /normal

!p.position= [x1, y2, x2, y3]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

;oplot, x_thin, z_thin, psym=3, color= 50

xyouts, 0.82, 0.95, 'thin disk', size=1.4, color=0, /normal, charthick=3.0

if keyword_set(msg1) then xyouts, 0.55, 0.715, msg1, size=1.8, color=0, /normal, charthick=3.0
if keyword_set(msg2) then xyouts, 0.55, 0.685, msg2, size=1.8, color=0, /normal, charthick=3.0

; ----------
;  3333333
; ----------

contour_makepic, x_other, z_other, y_other, m_other, xlen, xz=xz, yz=yz, hsml=hsml_other, $
        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

tv, NxNImage, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal


!p.position= [x0, y1, x1, y2]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

;oplot, x_other, z_other, psym=3, color= 120

xyouts, 0.21, 0.62, 'moderate-rotating', size=1.4, color=0, /normal, charthick=3.0


; ----------
;  4444444
; ----------

contour_makepic, x_nonrot, z_nonrot, y_nonrot, m_nonrot, xlen, xz=xz, yz=yz, hsml=hsml_nonrot, $
        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

tv, NxNImage, x1, y1, xsize=(x2-x1), ysize=(y2-y1), /normal


!p.position= [x1, y1, x2, y2]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

;oplot, x_nonrot, z_nonrot, psym=3, color= 150

xyouts, 0.77, 0.62, 'non-rotating', size=1.4, color=0, /normal, charthick=3.0

b_d_ratio= total(m_nonrot)/(total(m_thin)+total(m_other))

b_d_lbl= strcompress(string(b_d_ratio),/remove_all)
b_d_lbl= 'B/D= ' + strmid(b_d_lbl,0,4)
xyouts, 0.77, 0.59, b_d_lbl, size=1.4, color=0, /normal, charthick=3.0

b_t_ratio= total(m_nonrot)/(total(m_thin)+total(m_other)+total(m_nonrot))

b_t_lbl= strcompress(string(b_t_ratio),/remove_all)
b_t_lbl= 'B/T= ' + strmid(b_t_lbl,0,4)
xyouts, 0.77, 0.56, b_t_lbl, size=1.4, color=0, /normal, charthick=3.0


; ----------
;  5555555
; ----------

n_gas= fload_npart(0)

if n_gas gt 0 then begin
	etot= fload_gas_energy(1, /specific)
	etot=etot/1.0e+5
	print, "gas energy  max/min=", max(etot), min(etot)

	jznew= return_gas_jznew(1, n_hat=n_hat_allstars)
	jznew= jznew/1000.0
	print, "gas jznew  max/min=", max(jznew), min(jznew)

	e_bin= bins * (etot - e_min) / (e_max - e_min)
	e_idx= long(e_bin)
	jz_circ_e= jz_e(e_idx)
	e_j= jznew / jz_circ_e


	x= fload_gas_xyz('x')
	y= fload_gas_xyz('y')
	z= fload_gas_xyz('z')
	m= fload_gas_mass(1)
	hsml= fload_gas_hsml(1)

	process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

	idx= where(e_j ge 0.85)
	x_disk= x_new(idx)
	y_disk= -y_new(idx)
	z_disk= z_new(idx)
	m_disk= m(idx)
	hsml_disk= hsml(idx)

	idx= where(e_j lt 0.85)
	x_halo= x_new(idx)
	y_halo= -y_new(idx)
	z_halo= z_new(idx)
	m_halo= m(idx)
	hsml_halo= hsml(idx)
endif


fload_newcolortable, 1

; ----------
;  5555555
; ----------

if n_gas gt 0 then begin
	contour_makepic, x_disk, z_disk, y_disk, m_disk, xlen, xz=xz, yz=yz, hsml=hsml_disk, $
        	pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

	tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal
endif else begin
	xyouts, 0.18, 0.19, 'NO GAS', size=1.4, color=0, /normal, charthick=3.0
endelse

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

;oplot, x_disk, z_disk, psym=3, color= 100

xyouts, 0.28, 0.29, 'gas: thin disk', size=1.4, color=0, /normal, charthick=3.0


; ----------
;  666666
; ----------

if n_gas gt 0 then begin
	contour_makepic, x_halo, z_halo, y_halo, m_halo, xlen, xz=xz, yz=yz, hsml=hsml_halo, $
	        pixels=pixels, set_maxden= set_maxden, set_dynrng= set_dynrng, NxNImage=NxNImage

	tv, NxNImage, x1, y0, xsize=(x2-x1), ysize=(y1-y0), /normal
endif else begin
	xyouts, 0.18, 0.19, 'NO GAS', size=1.4, color=0, /normal, charthick=3.0
endelse

!p.position= [x1, y0, x2, y1]

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

;oplot, x_halo, z_halo, psym=3, color= 100

xyouts, 0.81, 0.29, 'gas: other', size=1.4, color=0, /normal, charthick=3.0


; done
; -----
device, /close



end







;====================================================================================




;
;  This procedure will calculate the bulge/disk mass,
;  and determine its composition.
;
;----------------------------------------
pro compute_bd, frun, snapnum, loadedsnap=loadedsnap


if not keyword_set(frun) then begin
	print, " "
	print, " compute_bd, frun, snapnum, loadedsnap=loadedsnap"
	print, " "
	return
endif

relaxationtime= 0.4

; determine last snapshot (if not provided)
if not keyword_set(snapnum) then begin

	spawn, "/bin/ls "+frun+"/snap* | wc ",result
	totnumsnap=long(result[0])-1
	snapnum=long(result[0])-4

	bhmtime= fload_bh_mergertime(frun)

	for i=0,totnumsnap do begin
	    ok= fload_snapshot_bh(frun,i,/header_only)
	    if fload_time(1) gt (bhmtime+relaxationtime) then begin
		snapnum= i
		break
	    endif
	endfor
endif



;---------------------------------------------------------

; assume it's a major merger and grab
ok= fload_snapshot_bh(frun,0,/header_only)

startid_primary= long(1)
npart_primary= long(total(fload_npart(99))/2.)
startid_sat= npart_primary+1
npart_sat= npart_primary

; high resolution runs
;startid_primary= 1
;npart_primary= 1600001
;startid_sat= 1600002
;npart_sat= 170001


; low resolution runs
;startid_primary= 1
;npart_primary= 533333
;startid_sat= 533334
;npart_sat= 66666    ; Im
;npart_sat= 286067   ; Sc
;npart_sat= 130208   ; Sd
;npart_sat= 533333   ; Sb


; std vc3 runs
;startid_primary= 1
;npart_primary= 200001
;startid_sat= 200002
;npart_sat= 200001

; bs runs
;startid_primary= 1
;npart_primary= 550001
;startid_sat= 550002
;npart_sat= 550001


;---------------------------------------------------------


; now grab snapshot to process
if not keyword_set(loadedsnap) then ok= fload_snapshot_bh(frun,snapnum)



;----------------------------------------------------------
;
;  Now split things up by component
;
;


; -------------------------
;  All Stars
; -------------------------

; x and z are already fixed to all star quantities

; specific binding energy
etot= fload_allstars_energy(1, /specific)
etot=etot/1.0e+5
print, "energy  max/min=", max(etot), min(etot)

; angular momentum (z)
jznew= return_allstars_jznew(1, n_hat=n_hat)
print, "n_hat= ", n_hat
n_hat_allstars=n_hat
jznew= jznew/1000.0
print, "jznew  max/min=", max(jznew), min(jznew)

; determine maximum j_z
bins= 200
;e_min= -2.0
;e_max= 0.8
e_min= min(etot)
e_max= max(etot)
determine_max_jz, bins, e_min, e_max, etot, jznew, e_e=e_e, jz_e=jz_e

; determine the e_j = j_z / j_circ(E)
;map the particle's energy to an above bin
e_bin= bins * (etot - e_min) / (e_max - e_min)
e_idx= long(e_bin)
jz_circ_e= jz_e(e_idx)
e_j= jznew / jz_circ_e


ids= fload_allstars_ids(1)
m= fload_allstars_mass(1)
print, "Total Stellar Mass=          ", total(m)

idx= where(e_j ge 0.85)
m_thin= m(idx)
ids_thin= ids(idx)
print, "Thin Disk Mass=              ", total(m_thin)

idx= where((e_j lt 0.85) and (e_j gt 0.25))
m_other= m(idx)
ids_other= ids(idx)
print, "Other Disk Mass=             ", total(m_other)

idx= where(e_j le 0.25)
m_nonrot= m(idx)
ids_nonrot= ids(idx)
print, "Non-rotating component Mass= ", total(m_nonrot)



b_d_ratio= total(m_nonrot)/(total(m_thin)+total(m_other))



;---------------------------------------------------------

;---------------------------------------------------------

print, "------"
g1idx= where((ids_thin ge startid_primary) and (ids_thin lt (startid_primary+npart_primary))) 
g2idx= where((ids_thin ge startid_sat) and (ids_thin lt (startid_sat+npart_sat))) 
print, "Thin Primary/Sat= ", total(m_thin(g1idx)), total(m_thin(g2idx))

g1idx= where((ids_other ge startid_primary) and (ids_other lt (startid_primary+npart_primary))) 
g2idx= where((ids_other ge startid_sat) and (ids_other lt (startid_sat+npart_sat))) 
print, "Other Primary/Sat= ", total(m_other(g1idx)), total(m_other(g2idx))

g1idx= where((ids_nonrot ge startid_primary) and (ids_nonrot lt (startid_primary+npart_primary))) 
g2idx= where((ids_nonrot ge startid_sat) and (ids_nonrot lt (startid_sat+npart_sat))) 
print, "Non-rot Primary/Sat= ", total(m_nonrot(g1idx)), total(m_nonrot(g2idx))


print, "------"



; ------------
;
;  Primary
;


if npart_primary gt 0 then begin 

	; old stars
        etot= fload_1gal_oldstars_energy(startid_primary, npart_primary, /specific) & etot=etot/1.0e+5
        jznew= return_oldstars_jznew(startid_primary, npart=npart_primary, n_hat=n_hat_allstars) & jznew= jznew/1000.0
	g1osm= fload_1gal_oldstars_mass(startid_primary,npart_primary)

        e_bin= bins * (etot - e_min) / (e_max - e_min)
        e_idx= long(e_bin)
        jz_circ_e= jz_e(e_idx)
        e_j= jznew / jz_circ_e

        idx= where(e_j ge 0.85)
        g1_os_m_thin= g1osm(idx)

	idx= where((e_j lt 0.85) and (e_j gt 0.25))
	g1_os_m_other= g1osm(idx)

        idx= where(e_j le 0.25)
        g1_os_m_nonrot= g1osm(idx)


        ; new stars
        etot= fload_1gal_newstars_energy(startid_primary, npart_primary, /specific) & etot=etot/1.0e+5
        jznew= return_newstars_jznew(startid_primary, npart=npart_primary, n_hat=n_hat_allstars) & jznew= jznew/1000.0
	g1nsm= fload_1gal_newstars_mass(startid_primary,npart_primary)

        e_bin= bins * (etot - e_min) / (e_max - e_min)
        e_idx= long(e_bin)
        jz_circ_e= jz_e(e_idx)
        e_j= jznew / jz_circ_e

        idx= where(e_j ge 0.85)
        if idx(0) ne -1 then g1_ns_m_thin= g1nsm(idx) else g1_ns_m_thin= 0.0

	idx= where((e_j lt 0.85) and (e_j gt 0.25))
	if idx(0) ne -1 then g1_ns_m_other= g1nsm(idx) else g1_ns_m_other= 0.0

        idx= where(e_j le 0.25)
        g1_ns_m_nonrot= g1nsm(idx)

endif




; --------------
;
;  Satellite
;


if npart_sat gt 0 then begin

        ; old stars
        etot= fload_1gal_oldstars_energy(startid_sat, npart_sat, /specific) & etot=etot/1.0e+5
        jznew= return_oldstars_jznew(startid_sat, npart=npart_sat, n_hat=n_hat_allstars) & jznew= jznew/1000.0
        g2osm= fload_1gal_oldstars_mass(startid_sat,npart_sat)

        e_bin= bins * (etot - e_min) / (e_max - e_min)
        e_idx= long(e_bin)
        jz_circ_e= jz_e(e_idx)
        e_j= jznew / jz_circ_e

        idx= where(e_j ge 0.85)
        g2_os_m_thin= g2osm(idx)

        idx= where((e_j lt 0.85) and (e_j gt 0.25))
        g2_os_m_other= g2osm(idx)

        idx= where(e_j le 0.25)
        g2_os_m_nonrot= g2osm(idx)


        ; new stars
        etot= fload_1gal_newstars_energy(startid_sat, npart_sat, /specific) & etot=etot/1.0e+5
        jznew= return_newstars_jznew(startid_sat, npart=npart_sat, n_hat=n_hat_allstars) & jznew= jznew/1000.0
        g2nsm= fload_1gal_newstars_mass(startid_sat,npart_sat)

        e_bin= bins * (etot - e_min) / (e_max - e_min)
        e_idx= long(e_bin)
        jz_circ_e= jz_e(e_idx)
        e_j= jznew / jz_circ_e

        idx= where(e_j ge 0.85)
        if idx(0) ne -1 then g2_ns_m_thin= g2nsm(idx) else g2_ns_m_thin= 0.0

        idx= where((e_j lt 0.85) and (e_j gt 0.25))
	if idx(0) ne -1 then g2_ns_m_other= g2nsm(idx) else g2_ns_m_other= 0.0

        idx= where(e_j le 0.25)
        g2_ns_m_nonrot= g2nsm(idx)

endif





;---------------------------------------------------------------------------


print, "running totals"
print, "total mass=  ", total(g1osm)+ total(g1nsm)+ total(g2osm)+ total(g2nsm)
print, "total thin=  ", total(g1_os_m_thin)+ total(g1_ns_m_thin)+ total(g2_os_m_thin)+ total(g2_ns_m_thin)
print, "total other= ", total(g1_os_m_other)+ total(g1_ns_m_other)+ total(g2_os_m_other)+ total(g2_ns_m_other)
print, "total nonr=  ", total(g1_os_m_nonrot)+ total(g1_ns_m_nonrot)+ total(g2_os_m_nonrot)+ total(g2_ns_m_nonrot)




; -----------------
;  write to a file
; -----------------
spawn, 'mkdir '+frun+'/bulgedisk', result
print, 'mkdir '+frun+'/bulgedisk'
print, result
exts='0000'
exts=exts+strcompress(string(snapnum),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)
openw, 1, frun+'/bulgedisk/'+exts+'.txt', ERROR=err

printf, 1, "#   bulgedisk/xxx.txt information from snap= ", snapnum
printf, 1, "# "
printf, 1, "# "
printf, 1, "# "
printf, 1, "# Total "
printf, 1, total(m),"       # total stellar mass"
printf, 1, total(m_thin),"       # total thin disk mass"
printf, 1, total(m_other),"       # total mod-rotating, thick, or other mass"
printf, 1, total(m_nonrot),"       # total non-rotating stellar mass"
printf, 1, "# "
printf, 1, "# Primary : Old Stellar Component "
printf, 1, total(g1osm),"       # total stellar mass"
printf, 1, total(g1_os_m_thin),"       # total thin disk mass"
printf, 1, total(g1_os_m_other),"       # total mod-rotating, thick, or other mass"
printf, 1, total(g1_os_m_nonrot),"       # total non-rotating stellar mass"
printf, 1, "# Primary : New Stellar Component "
printf, 1, total(g1nsm),"       # total stellar mass"
printf, 1, total(g1_ns_m_thin),"       # total thin disk mass"
printf, 1, total(g1_ns_m_other),"       # total mod-rotating, thick, or other mass"
printf, 1, total(g1_ns_m_nonrot),"       # total non-rotating stellar mass"
printf, 1, "# "
printf, 1, "# Satellite : Old Stellar Component "
printf, 1, total(g2osm),"       # total stellar mass"
printf, 1, total(g2_os_m_thin),"       # total thin disk mass"
printf, 1, total(g2_os_m_other),"       # total mod-rotating, thick, or other mass"
printf, 1, total(g2_os_m_nonrot),"       # total non-rotating stellar mass"
printf, 1, "# Satellite : New Stellar Component "
printf, 1, total(g2nsm),"       # total stellar mass"
printf, 1, total(g2_ns_m_thin),"       # total thin disk mass"
printf, 1, total(g2_ns_m_other),"       # total mod-rotating, thick, or other mass"
printf, 1, total(g2_ns_m_nonrot),"       # total non-rotating stellar mass"

close, 1


; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            "
print, "---------------------------------"
print, "---------------------------------"






end



;====================================================================================
;
;
;   Read the above txt file
;
;
;
;
;pro read_bulge_file, filename, m, m_thin, m_other, m_nonrot, $
pro read_bulge_file, frun, snapnum, m, m_thin, m_other, m_nonrot, $
			g1osm, g1_os_m_thin, g1_os_m_other, g1_os_m_nonrot, $
			g1nsm, g1_ns_m_thin, g1_ns_m_other, g1_ns_m_nonrot, $
			g2osm, g2_os_m_thin, g2_os_m_other, g2_os_m_nonrot, $
			g2nsm, g2_ns_m_thin, g2_ns_m_other, g2_ns_m_nonrot


exts='0000'
exts=exts+strcompress(string(snapnum),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)
spawn, "/bin/ls "+frun+"/bulgedisk/"+exts+".txt", result
if strlen(result) eq 0 then return
filename=result


openr, 1, filename
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & m= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & m_thin= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & m_other= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & m_nonrot= float(tempjunk(0))
readf, 1, junk
readf, 1, junk
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1osm= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1_os_m_thin= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1_os_m_other= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1_os_m_nonrot= float(tempjunk(0))
readf, 1, junk
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1nsm= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1_ns_m_thin= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1_ns_m_other= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g1_ns_m_nonrot= float(tempjunk(0))
readf, 1, junk
readf, 1, junk
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2osm= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2_os_m_thin= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2_os_m_other= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2_os_m_nonrot= float(tempjunk(0))
readf, 1, junk
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2nsm= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2_ns_m_thin= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2_ns_m_other= float(tempjunk(0))
readf, 1, junk & tempjunk= strsplit(junk,/extract,count=count) & g2_ns_m_nonrot= float(tempjunk(0))

close, 1




end








;====================================================================================
;
;
;   Returns the z-component of the angular momentum for every particle.
;   Where z is defined by the total angular momentum within 5 kpc.
;
;
;

function return_allstars_jznew, junk, radius=radius, $
				n_hat=n_hat, $
				npart=npart


if not keyword_set(radius) then radius=5.0


; -------------------------
; angular momentum (z)
; -------------------------

; positions
; ----------
if keyword_set(npart) then begin
	r= fload_1gal_allstars_xyz('r', startid, npart)
	jx= fload_1gal_allstars_j(21,startid=junk,numpart=npart)   ; specific j_x
	jy= fload_1gal_allstars_j(22,startid=junk,numpart=npart)   ; specific j_y
	jz= fload_1gal_allstars_j(23,startid=junk,numpart=npart)   ; specific j_z
endif else begin
	;x= fload_allstars_xyz('x')
	;y= fload_allstars_xyz('y')
	;z= fload_allstars_xyz('z')
	r= fload_allstars_xyz('r')
	comvel= fload_all_comvel(1)
	print, "comvel= ", comvel
   
	; total j's
	;jx= fload_allstars_j(11)   ; j_x
	;jy= fload_allstars_j(12)   ; j_y
	;jz= fload_allstars_j(13)   ; j_z
	jx= fload_allstars_j(21)   ; specific j_x
	jy= fload_allstars_j(22)   ; specific j_y
	jz= fload_allstars_j(23)   ; specific j_z
endelse

; angular momentum within some radius
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

print, "angular momentum vector= ", nx, ny, nz
print, "              n_hat_tot= ", n_hat_tot
print, "                  theta= ", acos(nz)*180./!PI

; determine jz (i.e., w.r.t. n-hat)
jznew= jx*nx + jy*ny + jz*nz


print, "all stars j_z max/min= ",  max(jznew), min(jznew)

return, jznew




end






;====================================================================================




function return_newstars_jznew, junk, radius=radius, $
				n_hat=n_hat, $
				npart= npart


if not keyword_set(radius) then radius=5.0


; -------------------------
if not keyword_set(n_hat) then begin

	; positions
	; ----------
	r= fload_allstars_xyz('r')
	comvel= fload_all_comvel(1)
	print, "comvel= ", comvel
   
	; total j's
	;jx= fload_allstars_j(11)   ; j_x
	;jy= fload_allstars_j(12)   ; j_y
	;jz= fload_allstars_j(13)   ; j_z
	jx= fload_allstars_j(21)   ; specific j_x
	jy= fload_allstars_j(22)   ; specific j_y
	jz= fload_allstars_j(23)   ; specific j_z

	; angular momentum within some radius
	idx= where(r lt radius)
	if idx(0) eq -1 then stop
	jtot_x= total(jx(idx))
	jtot_y= total(jy(idx))
	jtot_z= total(jz(idx))
	jtot= sqrt(jtot_x*jtot_x + jtot_y*jtot_y + jtot_z*jtot_z)
	nx= jtot_x/jtot
	ny= jtot_y/jtot
	nz= jtot_z/jtot

	print, "angular momentum vector= ", nx, ny, nz
	print, "                  theta= ", acos(nz)*180./!PI

endif else begin
	nx= n_hat[0]
	ny= n_hat[1]
	nz= n_hat[2]
endelse


; now we have the unit vector for the total angular momentum

if keyword_set(npart) then begin
	jx= fload_1gal_newstars_j(21,startid=junk,numpart=npart)   ; specific j_x
	jy= fload_1gal_newstars_j(22,startid=junk,numpart=npart)   ; specific j_y
	jz= fload_1gal_newstars_j(23,startid=junk,numpart=npart)   ; specific j_z
endif else begin
	jx= fload_newstars_j(21)   ; specific j_x
	jy= fload_newstars_j(22)   ; specific j_y
	jz= fload_newstars_j(23)   ; specific j_z
endelse


; determine jz (i.e., w.r.t. n-hat)
jznew= jx*nx + jy*ny + jz*nz


print, "new stars j_z max/min= ", max(jznew), min(jznew)

return, jznew



end






;====================================================================================





function return_oldstars_jznew, junk, radius=radius, $
				n_hat=n_hat, $
				npart= npart


if not keyword_set(radius) then radius=5.0


; -------------------------
if not keyword_set(n_hat) then begin

	; positions
	; ----------
	r= fload_allstars_xyz('r')
	comvel= fload_all_comvel(1)
	print, "comvel= ", comvel
   
	; total j's
	;jx= fload_allstars_j(11)   ; j_x
	;jy= fload_allstars_j(12)   ; j_y
	;jz= fload_allstars_j(13)   ; j_z
	jx= fload_allstars_j(21)   ; specific j_x
	jy= fload_allstars_j(22)   ; specific j_y
	jz= fload_allstars_j(23)   ; specific j_z

	; angular momentum within some radius
	idx= where(r lt radius)
	if idx(0) eq -1 then stop
	jtot_x= total(jx(idx))
	jtot_y= total(jy(idx))
	jtot_z= total(jz(idx))
	jtot= sqrt(jtot_x*jtot_x + jtot_y*jtot_y + jtot_z*jtot_z)
	nx= jtot_x/jtot
	ny= jtot_y/jtot
	nz= jtot_z/jtot

	print, "angular momentum vector= ", nx, ny, nz
	print, "                  theta= ", acos(nz)*180./!PI

endif else begin
	nx= n_hat[0]
	ny= n_hat[1]
	nz= n_hat[2]
endelse


; now we have the unit vector for the total angular momentum

if keyword_set(npart) then begin
	jx= fload_1gal_oldstars_j(21,startid=junk,numpart=npart)   ; specific j_x
	jy= fload_1gal_oldstars_j(22,startid=junk,numpart=npart)   ; specific j_y
	jz= fload_1gal_oldstars_j(23,startid=junk,numpart=npart)   ; specific j_z
endif else begin
	jx= fload_oldstars_j(21)   ; specific j_x
	jy= fload_oldstars_j(22)   ; specific j_y
	jz= fload_oldstars_j(23)   ; specific j_z
endelse


; determine jz (i.e., w.r.t. n-hat)
jznew= jx*nx + jy*ny + jz*nz


print, "old stars j_z max/min= ", max(jznew), min(jznew)

return, jznew



end





;====================================================================================






function return_gas_jznew, junk, radius=radius, $
				n_hat=n_hat


if not keyword_set(radius) then radius=5.0


; -------------------------
if not keyword_set(n_hat) then begin

	; positions
	; ----------
	r= fload_allstars_xyz('r')
	comvel= fload_all_comvel(1)
	print, "comvel= ", comvel
   
	; total j's
	;jx= fload_allstars_j(11)   ; j_x
	;jy= fload_allstars_j(12)   ; j_y
	;jz= fload_allstars_j(13)   ; j_z
	jx= fload_allstars_j(21)   ; specific j_x
	jy= fload_allstars_j(22)   ; specific j_y
	jz= fload_allstars_j(23)   ; specific j_z

	; angular momentum within some radius
	idx= where(r lt radius)
	if idx(0) eq -1 then stop
	jtot_x= total(jx(idx))
	jtot_y= total(jy(idx))
	jtot_z= total(jz(idx))
	jtot= sqrt(jtot_x*jtot_x + jtot_y*jtot_y + jtot_z*jtot_z)
	nx= jtot_x/jtot
	ny= jtot_y/jtot
	nz= jtot_z/jtot

	print, "angular momentum vector= ", nx, ny, nz
	print, "                  theta= ", acos(nz)*180./!PI

endif else begin
	nx= n_hat[0]
	ny= n_hat[1]
	nz= n_hat[2]
endelse


; now we have the unit vector for the total angular momentum

jx= fload_gas_j(21)   ; specific j_x
jy= fload_gas_j(22)   ; specific j_y
jz= fload_gas_j(23)   ; specific j_z


; determine jz (i.e., w.r.t. n-hat)
jznew= jx*nx + jy*ny + jz*nz


print, "gas j_z max/min= ", max(jznew), min(jznew)

return, jznew



end




;====================================================================================





pro determine_max_jz, bins, e_min, e_max, etot, jznew, $
				e_e=e_e, jz_e=jz_e


; make some binned data
;bins= 200
;e_min= -2.0
;e_max= 0.8
e_i= (e_max - e_min) / bins

; jz at specific binding energy
jz_e= fltarr(bins)
e_e= fltarr(bins)
for i=0, bins-1 do begin
        e_low= i * e_i + e_min
        e_high= e_low + e_i
        e_e[i]= 0.5 * (e_low + e_high)

        e_in= where((etot ge e_low) and (etot lt e_high))
        if e_in(0) ne -1 then begin
                jznew_in= jznew(e_in)
                jznew_in_sorted= jznew_in(sort(jznew_in))
                n_idx= n_elements(jznew_in)

                jz_e[i]= mean(jznew_in_sorted[long(n_idx*.9):n_idx-1])
        endif
endfor






end





;====================================================================================






;====================================================================================



