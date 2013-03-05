;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   Procedures related to the Tensor Virial Theorem
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------




; Calculates the Velocity, and Velocity Dispersion Tensors
;   (or, really, both the ordered and dispersion components
;    of the kinetic energy tensor)
; ---------------------------------------------------------------

;pro calc_velfield, junk
pro calc_velfield, frun


if not keyword_set(frun) then begin
	print, "  "
	print, " calc_velfield, frun"
	print, "  "
	return
endif

;frun= "/raid4/tcox/vc3vc3e"
;snapnum= 25
;snapnum= 31
;ok=fload_snapshot_bh(frun,snapnum)

spawn, "/bin/ls "+frun+"/snap* | wc ",result                                    
lastsnapshot=long(result[0])-1                                                  
ok=fload_snapshot_bh(frun,lastsnapshot)         

x= fload_allstars_xyz('x')
y= fload_allstars_xyz('y')
z= fload_allstars_xyz('z')
m= fload_allstars_mass(1)
ids= fload_allstars_ids(1)
comvel= fload_all_comvel(1)
vx= fload_allstars_v('x',comvel=comvel)
vy= fload_allstars_v('y',comvel=comvel)
vz= fload_allstars_v('z',comvel=comvel)

N= long(n_elements(x))

;Coord=fltarr(3,N)
;Coord=fltarr(6,N)
Coord=fltarr(7,N)
Masses= fltarr(N)
;Id= lonarr(N)
Velx= fltarr(N)
Vely= fltarr(N)
Velz= fltarr(N)

Coord(0,*)= x
Coord(1,*)= y
Coord(2,*)= z
Masses(*)= m
;Id(*)= ids
;Velx(*)= vx
;Vely(*)= vy
;Velz(*)= vz
Coord(3,*)= vx
Coord(4,*)= vy
Coord(5,*)= vz
Coord(6,*)= ids


DesNgb=96L
Hmax= 100.0



if(N gt 0) then begin

        VelEllipsoid= fltarr(6,N)      ;  full velocity field

        S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/VelocityEllipsoid/velocityellipsoid.so', $
                'velocityellipsoid', $
		N, $
                Coord, $
                Masses, $
                ;Velx, $
                ;Vely, $
                ;Velz, $
                ;Id, $
		DesNgb, $
		Hmax, $
                VelEllipsoid)
endif else begin
        print,'No stars, no problem.'
        return
endelse

help, VelEllipsoid

AvgVx= VelEllipsoid[0,*]
AvgDispx= VelEllipsoid[1,*]
AvgVy= VelEllipsoid[2,*]
AvgDispy= VelEllipsoid[3,*]
AvgVz= VelEllipsoid[4,*]
AvgDispz= VelEllipsoid[5,*]

;print, "--------------------------"
;print, "PI_xx   ", total(AvgDispx)
;print, "PI_yy   ", total(AvgDispy)
;print, "PI_zz   ", total(AvgDispz)
;print, "--------------------------"

openw, 1, frun+'/tvt.txt', ERROR=err
printf, 1, "#   "
printf, 1, "# Tensor Virial Theorem   "
printf, 1, "#   "
printf, 1, "#   "
printf, 1, "#   "

sigx2= vx-AvgVx
printf, 1, "PI_xx    ", total(m*sigx2*sigx2)

sigy2= vy-AvgVy
printf, 1, "PI_yy    ", total(m*sigy2*sigy2)

sigz2= vz-AvgVz
printf, 1, "PI_zz    ", total(m*sigz2*sigz2)

printf, 1, "PI_xy    ", total(m*sigx2*sigy2)
printf, 1, "PI_xz    ", total(m*sigx2*sigz2)
printf, 1, "PI_yz    ", total(m*sigy2*sigz2)


printf, 1, "T_xx    ", total(m*AvgVx*AvgVx)
printf, 1, "T_yy    ", total(m*AvgVy*AvgVy)
printf, 1, "T_zz    ", total(m*AvgVz*AvgVz)
printf, 1, "T_xy    ", total(m*AvgVx*AvgVy)
printf, 1, "T_xz    ", total(m*AvgVx*AvgVz)
printf, 1, "T_yz    ", total(m*AvgVy*AvgVz)

close, 1

;print, "--------------------------"

end





; =============================================================================





pro read_tvt_file, frun, pi_xx, pi_yy, pi_zz, $
			pi_xy, pi_xz, pi_yz, $
			t_xx, t_yy, t_zz, $
			t_xy, t_xz, t_yz


tvtfile= '/raid4/tcox/'+frun+'/tvt.txt'

openr, 1, tvtfile
junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

pi_xx= 0.0 & pi_yy= 0.0 & pi_zz= 0.0
pi_xy= 0.0 & pi_xz= 0.0 & pi_yz= 0.0
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & pi_xx= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & pi_yy= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & pi_zz= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & pi_xy= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & pi_xz= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & pi_yz= float(tempjunk(1))

t_xx= 0.0 & t_yy= 0.0 & t_zz= 0.0
t_xy= 0.0 & t_xz= 0.0 & t_yz= 0.0
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_xx= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_yy= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_zz= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_xy= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_xz= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_yz= float(tempjunk(1))

close, 1


end






; =============================================================================




; Calculates the Moment of Inertia Tensor
; ---------------------------------------------------------------

;pro calc_MoI, junk
pro calc_MoI, frun, snapnum


if not keyword_set(frun) then begin
	print, "  "
	print, " calc_MoI, frun, snapnum"
	print, "  "
	return
endif

;frun= "/raid4/tcox/vc3vc3e"
if not keyword_set(snapnum) then snapnum= 31

ok=fload_snapshot_bh(frun,snapnum)

x= fload_allstars_xyz('x')
y= fload_allstars_xyz('y')
z= fload_allstars_xyz('z')
m= fload_allstars_mass(1)
etot= fload_allstars_energy(1)
print, "energy  max/min=", max(etot), min(etot)


;e_sorted= etot(sort(etot))
;nstars= n_elements(m)
;halfenergy= etot[long(nstars/2.0)]
;halfenergy= etot[long(nstars/5.0)]
halfenergy= median(etot)
;halfenergy= median(etot) * 2.0
;halfenergy= median(etot) * 3.0
;halfenergy= max(etot)

print, "halfenergy= ", halfenergy
mostboundidx= where(etot lt halfenergy)
help, etot
help, mostboundidx

x= x(mostboundidx) 
y= y(mostboundidx) 
z= z(mostboundidx) 
m= m(mostboundidx) 


;print, "--------------------------"
;print, "--------------------------"

openw, 1, frun+'/MoI.txt', ERROR=err
printf, 1, "#   "
printf, 1, "# (most bound particles) Moment of Inertia"
printf, 1, "#  snap: ", snapnum
printf, 1, "#   "
printf, 1, "#   "

I_xx= total(m*x*x)
I_yy= total(m*y*y)
I_zz= total(m*z*z)
I_xy= total(m*x*y)
I_xz= total(m*x*z)
I_yz= total(m*y*z)
printf, 1, "I_xx    ", I_xx
printf, 1, "I_yy    ", I_yy
printf, 1, "I_zz    ", I_zz
printf, 1, "I_xy    ", I_xy
printf, 1, "I_xz    ", I_xz
printf, 1, "I_yz    ", I_yz

I_array= [[I_xx, I_xy, I_xz], $
	  [I_xy, I_yy, I_yz], $
	  [I_xz, I_yz, I_zz]]

Evs= EIGENQL(I_array, /double)

print, "Eigenvalues= ", Evs

; this is the strange Gonzalez-Garcia/van Albada (astro-ph/0506014)
; method, which doesn't seem to match anybody else's methods
; 
;a_over_b= sqrt((Evs(0)+Evs(1)-Evs(2))/(Evs(0)+Evs(2)-Evs(1)))
;print, Evs(0)+Evs(1)-Evs(2)
;print, Evs(0)+Evs(2)-Evs(1)
;print, Evs(1)+Evs(2)-Evs(0) 
;a_over_c= sqrt((Evs(0)+Evs(1)-Evs(2))/(Evs(1)+Evs(2)-Evs(0)))
;printf, 1, "a/b     ", a_over_b
;printf, 1, "a/c     ", a_over_c



; std. barnes '92, Hernquist '92/'93, and Springel 00 methods
bb= sqrt(Evs(1)/Evs(0))
cc= sqrt(Evs(2)/Evs(0))
TT= (1.0 - bb*bb)/(1.0 - cc*cc)
printf, 1, "b     ", bb
printf, 1, "c     ", cc
printf, 1, "T     ", TT


close, 1

;print, "--------------------------"

end






; --------------------------------------------------------------------------






pro read_moi_file, frun, i_xx, i_yy, i_zz, $
			i_xy, i_xz, i_yz, $
			bb, cc, TT


tvtfile= '/raid4/tcox/'+frun+'/MoI.txt'

openr, 1, tvtfile
junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

i_xx= 0.0 & i_yy= 0.0 & i_zz= 0.0
i_xy= 0.0 & i_xz= 0.0 & i_yz= 0.0
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & i_xx= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & i_yy= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & i_zz= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & i_xy= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & i_xz= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & i_yz= float(tempjunk(1))

bb= 0.0
cc= 0.0
TT= 0.0
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & bb= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & cc= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & TT= float(tempjunk(1))

close, 1


end







; ============================================================================
; ============================================================================
; ============================================================================





pro shapes, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "shapes, junk"
   print, "  "
   return
endif

filename='shape.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;--------------------------------------
;--------------------------------------


;xaxistitle= '1-c'
xaxistitle= '!6b'
xmax= 1.0
xmin= 0.3
;xmin= 0.0

;yaxistitle= '1-b'
yaxistitle= '!6c'
ymax = 1.0
ymin = 0.3
;ymin = 0.0

;----------------------------------------
; Generate plot
;----------------------------------------

x0= 0.15
y0= 0.15
x_size= 0.80
y_size= 0.80

!p.position=[x0, y0, x0+x_size,y0+y_size]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, $
        xrange=[xmin,xmax],$
        yrange=[ymin,ymax],$
        color=0, $
        xcharsize=1.5, ycharsize=1.5, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;xticks= 14, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /noerase, $
	/nodata


; --------------------

;do_regs= 1
do_regs= 0
if do_regs eq 1 then begin
   process_one_shape, 'vc3vc3b', pttype=1
   process_one_shape, 'vc3vc3c', pttype=1
   process_one_shape, 'vc3vc3d', pttype=1
   process_one_shape, 'vc3vc3e', pttype=1
   process_one_shape, 'vc3vc3f', pttype=1
   process_one_shape, 'vc3vc3g', pttype=1
   process_one_shape, 'vc3vc3h', pttype=1
   process_one_shape, 'vc3vc3i', pttype=1
   process_one_shape, 'vc3vc3j', pttype=1
   process_one_shape, 'vc3vc3k', pttype=1
   process_one_shape, 'vc3vc3l', pttype=1
   process_one_shape, 'vc3vc3m', pttype=1
   process_one_shape, 'vc3vc3n', pttype=1
   process_one_shape, 'vc3vc3o', pttype=1
   process_one_shape, 'vc3vc3p', pttype=1
   xyouts, 0.3, 0.80, '40% gas', /normal, charthick=3.0, size=1.5, color=150
endif

; --------------------

do_cs= 1
;do_cs= 0
if do_cs eq 1 then begin
   process_one_shape, 'collisionless/cvc3vc3b', pttype=2
   process_one_shape, 'collisionless/cvc3vc3c', pttype=2
   process_one_shape, 'collisionless/cvc3vc3d', pttype=2
   process_one_shape, 'collisionless/cvc3vc3e', pttype=2
   process_one_shape, 'collisionless/cvc3vc3f', pttype=2
   process_one_shape, 'collisionless/cvc3vc3g', pttype=2
   process_one_shape, 'collisionless/cvc3vc3h', pttype=2
   process_one_shape, 'collisionless/cvc3vc3i', pttype=2
   process_one_shape, 'collisionless/cvc3vc3j', pttype=2
   process_one_shape, 'collisionless/cvc3vc3k', pttype=2
   process_one_shape, 'collisionless/cvc3vc3l', pttype=2
   process_one_shape, 'collisionless/cvc3vc3m', pttype=2
   process_one_shape, 'collisionless/cvc3vc3n', pttype=2
   process_one_shape, 'collisionless/cvc3vc3o', pttype=2
   process_one_shape, 'collisionless/cvc3vc3p', pttype=2
   xyouts, 0.3, 0.80, 'dissipationless', /normal, charthick=3.0, size=1.5, color=50
endif

; --------------------


; straight line a=c
x=findgen(50)*(1.0/40.)
y=x
oplot, x, y, psym=-3, color=0, thick=3.0


; --------------------

posidx= where(x lt 1.0)
bb=x(posidx)
TT= 1.0/3.0
cc= sqrt(1.0 - (1.0 - bb*bb)/TT)
oplot, bb, cc, psym=-3, color=0, thick=1.0, linestyle= 1

TT= 2.0/3.0
cc= sqrt(1.0 - (1.0 - bb*bb)/TT)
oplot, bb, cc, psym=-3, color=0, thick=1.0, linestyle= 2

; --------------------


;xyouts, 0.3, 0.80, '40% gas', /normal, charthick=3.0, size=1.5, color=150
;xyouts, 0.3, 0.74, 'dissipationless', /normal, charthick=3.0, size=1.5, color=50


device, /close


end





;------------------------------------------------------------------------

pro process_one_shape, frun, pttype=pttype


read_moi_file, frun, i_xx, i_yy, i_zz, $
			i_xy, i_xz, i_yz, $
			bb, cc, TT

justletter= '!6'+strmid(frun,strlen(frun)-1,1)

;if pttype eq 1 then oplot, [1.0-cc], [1.0-bb], psym= 5, color= 150, symsize=1.2
if pttype eq 1 then begin
	;oplot, [bb], [cc], psym= 5, color= 150, symsize=1.2
	xyouts, bb, cc, justletter, color= 150, charthick= 3, size=1.5
endif

;if pttype eq 2 then oplot, [1.0-cc], [1.0-bb], psym= 2, color= 50, symsize=1.2
if pttype eq 2 then begin
	;oplot, [bb], [cc], psym= 2, color= 50, symsize=1.2
	xyouts, bb, cc, justletter, color= 50, charthick= 3, size=1.5
endif


end










;===================================================================================




pro triax_hist, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "triax_hist, junk"
   print, "  "
   return
endif

filename='triax.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;--------------------------------------
;--------------------------------------


;xaxistitle= '1-c'
xaxistitle= '!6Triaxiality, T'
xmax= 1.0
xmin= 0.0

;yaxistitle= '1-b'
;yaxistitle= 'c'
ymax = 1.2
ymin = 0.0

;----------------------------------------
; Generate plot
;----------------------------------------

x0= 0.02
y0= 0.15
x_size= 0.96
y_size= 0.83

!p.position=[x0, y0, x0+x_size,y0+y_size]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, $
        xrange=[xmin,xmax],$
        yrange=[ymin,ymax],$
        color=0, $
        xcharsize=1.5, ycharsize=1.5, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xticks= 10, $
        yticks= 1, $
	ytickformat='(a1)', $
	xtickname=['0',' ','0.2',' ','0.4',' ','0.6',' ','0.8',' ','1'], $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /noerase, $
	/nodata


; --------------------

read_moi_file, 'vc3vc3b', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [TT]
read_moi_file, 'vc3vc3c', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3d', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3e', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3f', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3g', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3h', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3i', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3j', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3k', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3l', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3m', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3n', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3o', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]
read_moi_file, 'vc3vc3p', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
alltts= [alltts, TT]

; --------------------

read_moi_file, 'collisionless/cvc3vc3b', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [TT]
read_moi_file, 'collisionless/cvc3vc3c', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3d', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3e', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3f', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3g', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3h', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3i', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3j', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3k', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3l', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3m', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3n', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3o', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]
read_moi_file, 'collisionless/cvc3vc3p', i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
allctts= [allctts, TT]

; --------------------


levels= 10.0
step= (xmax-xmin)/(levels)
bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
;bins=[bins,xmax]
bins=bins+(step*0.5)


; std histogram
hist_alltts= histogram(alltts, binsize=step, max=xmax, min=xmin)


; collisionless histogram
hist_allctts= histogram(allctts, binsize=step, max=xmax, min=xmin)

; make sure it extends to 0
bins= [0,bins]
hist_alltts= [hist_alltts(0),hist_alltts]
hist_allctts= [hist_allctts(0),hist_allctts]


normalization= float(max([hist_alltts,hist_allctts]))
print, "histogram normalization= ", normalization
hist_alltts= hist_alltts/normalization
hist_allctts= hist_allctts/normalization

oplot, bins, hist_alltts, psym=10, color=150, thick=4.0
;oplot, bins, hist_allctts, psym=10, color=50, thick=4.0
oplot, bins, hist_allctts, psym=10, color=50, thick=8.0

; fill in histogram
; ------------------
nbins= bins+0.05          ; make x coord
nbins[0]= 0.0
nbins=[nbins,nbins]
nbins= nbins(sort(nbins)) 

ntts= fltarr(2.*levels + 2) 
nctts= fltarr(2.*levels + 2) 
for i=1,levels do begin
   ntts[2.*i-1]= hist_alltts[i]
   ntts[2.*i]= hist_alltts[i]
   nctts[2.*i-1]= hist_allctts[i]
   nctts[2.*i]= hist_allctts[i]
endfor

; collisionless
;polyfill, nbins, nctts, /data, color= 50, /line_fill, linestyle=0, $
;				thick=3.0

; 40% gas
polyfill, nbins, ntts, /data, color= 150, /line_fill, linestyle=0, $
				thick=3.0, orientation=45.0



xyouts, 0.62, 0.85, '40% gas', /normal, charthick=3.0, size=1.5, color=150
xyouts, 0.62, 0.79, 'dissipationless', /normal, charthick=3.0, size=1.5, color=50


device, /close


end









;===================================================================================







pro delta_hist, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "delta_hist, junk"
   print, "  "
   return
endif

filename='delta.eps'
;filename='tk_hist.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;--------------------------------------
;--------------------------------------


xaxistitle= '!6Anisotropy, !7d!6'
;xaxistitle= 'T/K'
xmax= 1.0
xmin= 0.0

ymax = 1.2
ymin = 0.0

;----------------------------------------
; Generate plot
;----------------------------------------

x0= 0.02
y0= 0.15
x_size= 0.96
y_size= 0.83

!p.position=[x0, y0, x0+x_size,y0+y_size]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, /nodata, $
        color= 0, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, charthick=3.0, $
        xticks= 10, yticks= 1, ytickformat='(a1)', xtitle=xaxistitle, ytitle=yaxistitle, $
	xtickname=['0',' ','0.2',' ','0.4',' ','0.6',' ','0.8',' ','1']


; --------------------

read_tvt_file, 'vc3vc3b', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)
alltks=  (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)
read_tvt_file, 'vc3vc3c', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3d', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3e', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3f', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3g', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3h', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3i', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3j', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3k', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3l', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3m', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3n', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3o', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'vc3vc3p', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allds=  [allds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
alltks=  [alltks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]


; --------------------

read_tvt_file, 'cvc3vc3b', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)
allctks=  (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)
read_tvt_file, 'cvc3vc3c', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3d', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3e', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3f', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3g', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3h', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3i', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3j', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3k', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3l', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3m', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3n', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3o', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]
read_tvt_file, 'cvc3vc3p', pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
allcds=  [allcds, 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)]
allctks=  [allctks, (t_xx + t_yy + t_zz) / ( t_xx + t_yy + t_zz + pi_xx + pi_yy + pi_zz)]


; --------------------



levels= 10.0


; two arrays to histogram

; deltas
alltts= allds
allctts= allcds

; t/k
alltts= alltks
allctts= allctks


xyouts, 0.67, 0.85, '40% gas', /normal, charthick=3.0, size=1.5, color=150
xyouts, 0.67, 0.79, 'dissipationless', /normal, charthick=3.0, size=1.5, color=50



; --------------------


step= (xmax-xmin)/(levels)
bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
;bins=[bins,xmax]
bins=bins+(step*0.5)


; std histogram
hist_alltts= histogram(alltts, binsize=step, max=xmax, min=xmin)


; collisionless histogram
hist_allctts= histogram(allctts, binsize=step, max=xmax, min=xmin)

; make sure it extends to 0
bins= [0,bins]
hist_alltts= [hist_alltts(0),hist_alltts]
hist_allctts= [hist_allctts(0),hist_allctts]


normalization= float(max([hist_alltts,hist_allctts]))
print, "histogram normalization= ", normalization
hist_alltts= hist_alltts/normalization
hist_allctts= hist_allctts/normalization

oplot, bins, hist_alltts, psym=10, color=150, thick=4.0
;oplot, bins, hist_allctts, psym=10, color=50, thick=4.0
oplot, bins, hist_allctts, psym=10, color=50, thick=8.0

; fill in histogram
; ------------------
nbins= bins+(step*0.5)          ; make x coord
nbins[0]= 0.0
nbins=[nbins,nbins]
nbins= nbins(sort(nbins)) 

ntts= fltarr(2.*levels + 2) 
nctts= fltarr(2.*levels + 2) 
for i=1,levels do begin
   ntts[2.*i-1]= hist_alltts[i]
   ntts[2.*i]= hist_alltts[i]
   nctts[2.*i-1]= hist_allctts[i]
   nctts[2.*i]= hist_allctts[i]
endfor


; collisionless
;polyfill, nbins, nctts, /data, color= 50, /line_fill, linestyle=0, $
;				thick=3.0

; 40% gas
polyfill, nbins, ntts, /data, color= 150, /line_fill, linestyle=0, $
				thick=3.0, orientation=45.0


device, /close


end









;=================================================================================








pro doit, junk


print, " -- Collisional -- "

doit_it, 'vc3_iso'
doit_it, 'vc3vc3b'
doit_it, 'vc3vc3c'
doit_it, 'vc3vc3d'
doit_it, 'vc3vc3e'
doit_it, 'vc3vc3f'
doit_it, 'vc3vc3g'
doit_it, 'vc3vc3h'
doit_it, 'vc3vc3i'
doit_it, 'vc3vc3j'
doit_it, 'vc3vc3k'
doit_it, 'vc3vc3l'
doit_it, 'vc3vc3m'
doit_it, 'vc3vc3n'
doit_it, 'vc3vc3o'
doit_it, 'vc3vc3p'

print, " -- Collisionless -- "

doit_it, 'cvc3vc3b'
doit_it, 'cvc3vc3c'
doit_it, 'cvc3vc3d'
doit_it, 'cvc3vc3e'
doit_it, 'cvc3vc3f'
doit_it, 'cvc3vc3g'
doit_it, 'cvc3vc3h'
doit_it, 'cvc3vc3i'
doit_it, 'cvc3vc3j'
doit_it, 'cvc3vc3k'
doit_it, 'cvc3vc3l'
doit_it, 'cvc3vc3m'
doit_it, 'cvc3vc3n'
doit_it, 'cvc3vc3o'
doit_it, 'cvc3vc3p'


end




pro doit_it, frun

  read_tvt_file, frun, pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
  ;print, 2.0*(t_xx + t_yy + t_zz) / (pi_xx + pi_yy + pi_zz)
  ;print, pi_xx/pi_zz, pi_yy/pi_zz
  print, 2.0*pi_zz/(pi_xx+pi_yy), pi_xx/pi_yy

end






;===============================================================



pro virial_eq, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "virial_eq, junk"
   print, "  "
   return
endif

filename='virialeq.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;--------------------------------------
;--------------------------------------


xaxistitle= '!6I!Dzz!N / I!Dxy!N'
xmax= 1.0
xmin= 0.0

yaxistitle= '!6K!Dzz!N / K!Dxy!N'
ymax = 1.0
ymin = 0.0

;----------------------------------------
; Generate plot
;----------------------------------------

x0= 0.15
y0= 0.15
x_size= 0.80
y_size= 0.80

!p.position=[x0, y0, x0+x_size,y0+y_size]


plot, [0], [0], psym=3,xstyle=1,ystyle=1, $
        xrange=[xmin,xmax],$
        yrange=[ymin,ymax],$
        color=0, $
        xcharsize=1.5, ycharsize=1.5, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;xticks= 14, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /noerase, $
	/nodata



; --------------------

do_regs= 1
;do_regs= 0
if do_regs eq 1 then begin
   process_one_virialeq, 'vc3vc3b', pttype=1
   process_one_virialeq, 'vc3vc3c', pttype=1
   process_one_virialeq, 'vc3vc3d', pttype=1
   process_one_virialeq, 'vc3vc3e', pttype=1
   process_one_virialeq, 'vc3vc3f', pttype=1
   process_one_virialeq, 'vc3vc3g', pttype=1
   process_one_virialeq, 'vc3vc3h', pttype=1
   process_one_virialeq, 'vc3vc3i', pttype=1
   process_one_virialeq, 'vc3vc3j', pttype=1
   process_one_virialeq, 'vc3vc3k', pttype=1
   process_one_virialeq, 'vc3vc3l', pttype=1
   process_one_virialeq, 'vc3vc3m', pttype=1
   process_one_virialeq, 'vc3vc3n', pttype=1
   process_one_virialeq, 'vc3vc3o', pttype=1
   process_one_virialeq, 'vc3vc3p', pttype=1
   xyouts, 0.23, 0.80, '40% gas', /normal, charthick=3.0, size=1.5, color=150
endif

; --------------------

do_cs= 1
;do_cs= 0
if do_cs eq 1 then begin
   process_one_virialeq, 'collisionless/cvc3vc3b', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3c', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3d', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3e', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3f', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3g', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3h', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3i', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3j', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3k', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3l', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3m', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3n', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3o', pttype=2
   process_one_virialeq, 'collisionless/cvc3vc3p', pttype=2
   xyouts, 0.23, 0.85, 'dissipationless', /normal, charthick=3.0, size=1.5, color=50
endif

; --------------------





device, /close


end






;------------------------------------------------------------------------

pro process_one_virialeq, frun, pttype=pttype


read_tvt_file, frun, pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz

read_moi_file, frun, i_xx, i_yy, i_zz, $
			i_xy, i_xz, i_yz, $
			bb, cc, TT


k_zz= t_zz + 0.5*pi_zz
k_xx= t_xx + 0.5*pi_xx
k_yy= t_yy + 0.5*pi_yy

k_xy= sqrt(k_xx*k_xx + k_yy*k_yy)

Kzzxy= k_zz / k_xy
bc= bb/cc



bc= i_zz / sqrt(i_xx*i_xx + i_yy*i_yy)


print, frun, Kzzxy, bc


justletter= '!6'+strmid(frun,strlen(frun)-1,1)

;if pttype eq 1 then oplot, [1.0-cc], [1.0-bb], psym= 5, color= 150, symsize=1.2
if pttype eq 1 then begin
	;oplot, [bb], [cc], psym= 5, color= 150, symsize=1.2
	xyouts, bc, Kzzxy, justletter, color= 150, charthick= 3, size=1.5
endif

;if pttype eq 2 then oplot, [1.0-cc], [1.0-bb], psym= 2, color= 50, symsize=1.2
if pttype eq 2 then begin
	;oplot, [bb], [cc], psym= 2, color= 50, symsize=1.2
	xyouts, bc, Kzzxy, justletter, color= 50, charthick= 3, size=1.5
endif


end








; =============================================================================









;------------------------------------------------------------------------


; Calculates the Velocity, and Velocity Dispersion Tensors
;   (or, really, both the ordered and dispersion components
;    of the kinetic energy tensor)
;   as a function of radius
; ---------------------------------------------------------------

pro tvt_r, frun, snapnum= snapnum


if not keyword_set(frun) then begin
	print, "  "
	print, " tvt_r, frun, snapnum= snapnum"
	print, "  "
	return
endif

;frun= "/raid4/tcox/vc3vc3e"
;snapnum= 25
if not keyword_set(snapnum) then snapnum= 30
;ok=fload_snapshot_bh(frun,snapnum)

; this grabs the last snapshot
;spawn, "/bin/ls "+frun+"/snap* | wc ",result                                    
;lastsnapshot=long(result[0])-1                                                  
;ok=fload_snapshot_bh(frun,lastsnapshot)         



; all stars
; -----------
r_as= fload_allstars_xyz('r')

x= fload_allstars_xyz('x')
y= fload_allstars_xyz('y')
z= fload_allstars_xyz('z')
m= fload_allstars_mass(1)
ids= fload_allstars_ids(1)
comvel= fload_all_comvel(1)
vx= fload_allstars_v('x',comvel=comvel)
vy= fload_allstars_v('y',comvel=comvel)
vz= fload_allstars_v('z',comvel=comvel)

calc_tvt, x, y, z, m, vx, vy, vz, ids, $
		AvgVx=AvgVx_as, AvgDispx=AvgDispx_as, $
		AvgVx=AvgVy_as, AvgDispy=AvgDispy_as, $
		AvgVx=AvgVz_as, AvgDispz=AvgDispz_as


process_and_oplot_PIdivT, r_as, AvgVx_as, AvgDispx_as, $
				AvgVy_as, AvgDispy_as, $
				AvgVz_as, AvgDispz_as



; old stars
; -----------
r_as= fload_oldstars_xyz('r')

x= fload_oldstars_xyz('x')
y= fload_oldstars_xyz('y')
z= fload_oldstars_xyz('z')
m= fload_oldstars_mass(1)
ids= fload_oldstars_ids(1)
vx= fload_oldstars_v('x',comvel=comvel)   ; use total comvel
vy= fload_oldstars_v('y',comvel=comvel)
vz= fload_oldstars_v('z',comvel=comvel)

calc_tvt, x, y, z, m, vx, vy, vz, ids, $
                AvgVx=AvgVx_os, AvgDispx=AvgDispx_os, $
                AvgVx=AvgVy_os, AvgDispy=AvgDispy_os, $
                AvgVx=AvgVz_os, AvgDispz=AvgDispz_os


process_and_oplot_PIdivT, r_os, AvgVx_os, AvgDispx_os, $
				AvgVy_os, AvgDispy_os, $
				AvgVz_os, AvgDispz_os




; new stars
; -----------
r_as= fload_newstars_xyz('r')

x= fload_newstars_xyz('x')
y= fload_newstars_xyz('y')
z= fload_newstars_xyz('z')
m= fload_newstars_mass(1)
ids= fload_newstars_ids(1)
vx= fload_newstars_v('x',comvel=comvel)    ; use total comvel
vy= fload_newstars_v('y',comvel=comvel)
vz= fload_newstars_v('z',comvel=comvel)

calc_tvt, x, y, z, m, vx, vy, vz, ids, $
                AvgVx=AvgVx_ns, AvgDispx=AvgDispx_ns, $
                AvgVx=AvgVy_ns, AvgDispy=AvgDispy_ns, $
                AvgVx=AvgVz_ns, AvgDispz=AvgDispz_ns


process_and_oplot_PIdivT, r_ns, AvgVx_ns, AvgDispx_ns, $
				AvgVy_ns, AvgDispy_ns, $
				AvgVz_ns, AvgDispz_ns

end




; =============================================================================




pro calc_tvt, x, y, z, m, vx, vy, vz, ids, $
		AvgVx=AvgVx, AvgDispx=AvgDispx, $
		AvgVx=AvgVy, AvgDispy=AvgDispy, $
		AvgVx=AvgVz, AvgDispz=AvgDispz


N= long(n_elements(x))

Coord=fltarr(7,N)
Masses= fltarr(N)
Velx= fltarr(N)
Vely= fltarr(N)
Velz= fltarr(N)

Coord(0,*)= x
Coord(1,*)= y
Coord(2,*)= z
Masses(*)= m
Coord(3,*)= vx
Coord(4,*)= vy
Coord(5,*)= vz
Coord(6,*)= ids


DesNgb=96L
Hmax= 100.0



if(N gt 0) then begin

        VelEllipsoid= fltarr(6,N)      ;  full velocity field

        S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/VelocityEllipsoid/velocityellipsoid.so', $
                'velocityellipsoid', $
		N, $
                Coord, $
                Masses, $
                ;Velx, $
                ;Vely, $
                ;Velz, $
                ;Id, $
		DesNgb, $
		Hmax, $
                VelEllipsoid)
endif else begin
        print,'No stars, no problem.'
        return
endelse


AvgVx= VelEllipsoid[0,*]
AvgDispx= VelEllipsoid[1,*]
AvgVy= VelEllipsoid[2,*]
AvgDispy= VelEllipsoid[3,*]
AvgVz= VelEllipsoid[4,*]
AvgDispz= VelEllipsoid[5,*]


end





; =============================================================================





pro process_and_oplot_PIdivT, r_ns, AvgVx_ns, AvgDispx_ns, $
				AvgVy_ns, AvgDispy_ns, $
				AvgVz_ns, AvgDispz_ns








end








