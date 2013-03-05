;compile kernel_2d_functions.pro
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;	xray_map.pro
;
;	Reads in snapshots w/ black holes
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Set some x device options
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
LCTNUM = 0

;loadct,lctnum
;set_plot,'x'
;device,decompose=0,true_color=24,retain=2
f = 2.0*!Pi*findgen(16)/17.0
usersym,cos(f),sin(f),/fill

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make a grid for the grid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
grav_smoothing = 0.1
DynRange=1.0e4
AXIS0 = 0
AXIS1 = 1
XPIXELS = 400
YPIXELS = 400
;XMIN    = -15.0*.470
;YMIN    = -15.0*.470
;XMAX    = 15.0*.470
;YMAX    = 15.0*.470
XMIN    = -10.0;*.470
YMIN    = -10.0;*.470
XMAX    = 10.0;*.470
YMAX    = 10.0;*.470
lx    = XMAX-XMIN
ly    = YMAX-YMIN
dx    = (XMAX-XMIN)/float(XPIXELS)
dy    = (YMAX-YMIN)/float(YPIXELS)
grid = dblarr(XPIXELS,YPIXELS)
grid_soft   = dblarr(XPIXELS,YPIXELS)
grid_hard   = dblarr(XPIXELS,YPIXELS)
grid_medium = dblarr(XPIXELS,YPIXELS)
gsum = dblarr(XPIXELS,YPIXELS)
gsum_soft  = dblarr(XPIXELS,YPIXELS)
gsum_hard  = dblarr(XPIXELS,YPIXELS)
gsum_medium= dblarr(XPIXELS,YPIXELS)
grid(*,*) = 0.0D
grid_soft(*,*) = 0.0D
grid_hard(*,*) = 0.0D
grid_medium(*,*) = 0.0D
gsum(*,*) = 0.0D
gsum_soft(*,*) = 0.0D
gsum_hard(*,*) = 0.0D
gsum_medium(*,*) = 0.0D
image  = dblarr(XPIXELS,YPIXELS)
image(*,*) = 0.0D
p_x   = fltarr(XPIXELS)
p_y   = fltarr(XPIXELS)
p_x(*) = dx*findgen(XPIXELS)
p_y(*) = dy*findgen(YPIXELS)

n_table = 1000
kernel_table =fltarr(n_table)


for i=0,n_table-1 do kernel_table(i) =  kernel_2d(float(i)/float(n_table-1))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Read in snapshot
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fdir = "A3/"
snapbase = "snapshot"


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Select snapshot to read in
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

cnum   = ' '
read,cnum,prompt='Enter snapshot number: '
NUM   = LONG(cnum)
exts='000'
exts=exts+strcompress(string(num),/remove_all) ; sets exts to snap number
exts=strmid(exts,strlen(exts)-3,3)

fsnap = fdir+snapbase +"_"+exts	;snapshot filename


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; set some units if you need them
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
UnitLength_in_cm =        3.085678d21        ;  1.0 kpc 
UnitMass_in_g    =        1.989d43 ;  1.0e10 solar masses 
UnitVelocity_in_cm_per_s = 1d5 ;  1 km/sec  
UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s 
UnitDensity_in_cgs= UnitMass_in_g/ UnitLength_in_cm^3
UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ UnitTime_in_s^2
UnitEnergy_in_cgs= UnitMass_in_g * UnitLength_in_cm^2 / UnitTime_in_s^2
GRAVITY   = 6.672d-8
BOLTZMANN = 1.3806d-16
PROTONMASS = 1.6726e-24
G=GRAVITY/ UnitLength_in_cm^3 * UnitMass_in_g * UnitTime_in_s^2
Xh=0.76D  ; mass fraction of hydrogen
gamma= 5.0/3


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Set variables to read header info
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

time = 0.0D
redshift = 0.0D
BoxSize = 0.0D
HubbleParam = 0.0D
Omega0 = 0.0D
OmegaLambda = 0.0D
flag_sfr= 0L
flag_feedback= 0L
npartTotal =lonarr(6)
flag_cooling= 0L
flag_multiphase= 0L
flag_stellarage= 0L
flag_sfrhistogram= 0L
num_files= 0L
seed = 0L
npart = lonarr(6)
massarr = dblarr(6)
npartTotal = lonarr(6)
slack_array = lonarr(21)


;;;;;;;;;;;;;;;;;;;;;;;;;;;
; open snapshot
;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,1,fsnap,/f77_unformatted


;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read header
;;;;;;;;;;;;;;;;;;;;;;;;;;;
readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartTotal,flag_cooling,num_files,BoxSize,Omega0,OmegaLambda,HubbleParam,flag_multiphase,flag_stellarage,flag_sfrhistogram,slack_array
time  = float(time)
HubbleParam = float(HubbleParam)

print,'npart',npart
print,'massaray',massarr
print,"a = ", time
print,"z = ",redshift

;;;;;;;;;;;;;;;;;;;;;;;;;
; set numbers of particles
;;;;;;;;;;;;;;;;;;;;;;;;;

NGas  =npart(0)
NHalo =npart(1)
NDisk =npart(2)
NBulge=npart(3)
NStars=npart(4)
NBH   =npart(5)
N=npart(0)+npart(1)+npart(2)+npart(3)+npart(4)+npart(5)


;;;;;;;;;;;;;;;;;;;;;;;;;
; make data arrays
;;;;;;;;;;;;;;;;;;;;;;;;;

pos		=fltarr(3,N)				;positions
vel		=fltarr(3,N) 				;velocities
id		=lonarr(N)				;particle ids
if (npart(0)+npart(4)+npart(5) gt 0) then gsbhmass	=fltarr(npart(0)+npart(4)+npart(5))	;gas, stellar, and black hole masses
if npart(0) gt 0 then begin
u		=fltarr(npart(0))			;gas internal energy
rho		=fltarr(npart(0))			;gas density
nel		=fltarr(npart(0))			;gas electron density
nh		=fltarr(npart(0))			;gas hydrogen density
hsmltmp		=fltarr(npart(0))			;gas smoothing length
sfr		=fltarr(npart(0))			;gas star formation rate
endif
age 		=fltarr(N)
if npart(4) gt 0 then begin				;if stellar particles
agetmp		=fltarr(npart(4))			;stellar ages
endif

zm		=fltarr(N)
if npart(0)+npart(4) gt 0 then zmtmp		=fltarr(npart(0)+npart(4))		;gas and stellar metallicities
hsml		=fltarr(N)				;particle smoothing lengths
pot		=fltarr(N)				;particle potentials
mass 		= fltarr(N)				;masses
types 		= lonarr(N)				;types




;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Read in data 
;
;;;;;;;;;;;;;;;;;;;;;;;;;;

readu,1,pos
readu,1,vel
readu,1,id 
if npart(0)+npart(4)+npart(5) gt 0 then readu,1,gsbhmass
if npart(0) gt 0 then begin
readu,1,u
readu,1,rho
readu,1,nel
readu,1,nh
readu,1,hsmltmp
readu,1,sfr
endif
if npart(4) gt 0 then readu,1,agetmp
if npart(0) + npart(4) gt 0 then readu,1,zmtmp
close,1


;;;;;;;;;;;;;;;;;;;;;;;;;;
; set particle masses
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;
;gas mass
;;;;;;;;;;;;;;;;;;;;;;;;;;

if npart(0) gt 0 then mass(0:npart(0)-1) = gsbhmass(0:npart(0)-1)


;;;;;;;;;;;;;;;;;;;;;;;;;;
;hrdm mass
;;;;;;;;;;;;;;;;;;;;;;;;;;
if npart(1) gt 0 then mass(npart(0):npart(0)+npart(1)-1) = float(massarr(1))

;;;;;;;;;;;;;;;;;;;;;;;;;;
;disk mass
;;;;;;;;;;;;;;;;;;;;;;;;;;
if npart(2) gt 0 then mass(npart(0)+npart(1):npart(0)+npart(1)+npart(2)-1) = float(massarr(2))

;;;;;;;;;;;;;;;;;;;;;;;;;;
;bulge mass
;;;;;;;;;;;;;;;;;;;;;;;;;;
if npart(3) gt 0 then mass(npart(0)+npart(1)+npart(2):npart(0)+npart(1)+npart(2)+npart(3)-1) = float(massarr(3))


;;;;;;;;;;;;;;;;;;;;;;;;;;
;stellar mass
;;;;;;;;;;;;;;;;;;;;;;;;;;
if npart(4) gt 0 then mass(npart(0)+npart(1)+npart(2)+npart(3):npart(0)+npart(1)+npart(2)+npart(3)+npart(4)-1) = gsbhmass(npart(0):npart(0)+npart(4)-1)

;;;;;;;;;;;;;;;;;;;;;;;;;;
;black hole masses
;;;;;;;;;;;;;;;;;;;;;;;;;;
if npart(5) gt 0 then mass(npart(0)+npart(1)+npart(2)+npart(3)+npart(4):npart(0)+npart(1)+npart(2)+npart(3)+npart(4)+npart(5)-1) = gsbhmass(npart(0)+npart(4):npart(0)+npart(4)+npart(5) -1)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set mass units to solar masses
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
mass(*) = mass(*)*10.0^10.0
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set particle types
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if npart(0) gt 0 then types(0:npart(0)-1)=0
if npart(1) gt 0 then types(npart(0):npart(0)+npart(1)-1)=1
if npart(2) gt 0 then types(npart(0)+npart(1):npart(0)+npart(1)+npart(2)-1)=2
if npart(3) gt 0 then types(npart(0)+npart(1)+npart(2):npart(0)+npart(1)+npart(2)+npart(3)-1)=3
if npart(4) gt 0 then types(npart(0)+npart(1)+npart(2)+npart(3):npart(0)+npart(1)+npart(2)+npart(3)+npart(4)-1)=4
if npart(5) gt 0 then types(npart(0)+npart(1)+npart(2)+npart(3)+npart(4):npart(0)+npart(1)+npart(2)+npart(3)+npart(4)+npart(5)-1)=5

Nhalo  = 0L
Nbulge = 0L
Ndisk  = 0L
Ngas   = 0L
Nstars = 0L
Nbh    = 0L

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find particle type indicies
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

gas_index  = where(types eq 0,Ngas)
hrdm_index = where(types eq 1,Nhalo)
disk_index = where(types eq 2,Ndisk)
bulge_index = where(types eq 3,Nbulge)
star_index = where(types eq 4,Nstars)
dbs_index = where(types eq 2 or types eq 3 or types eq 4,Ndbs)
sg_index   = where(types eq 4 or types eq 0,Nsg)
bh_index   = where(types eq 5,Nbh)


print,'Ngas   :',Ngas
print,'Nhalo  :',Nhalo
print,'Ndisk  :',Ndisk
print,'Nbulge :',Nbulge
print,'Nstars :',Nstars
print,'Nbh    :',Nbh


;;;;;;;;;
;
; Adjust metallicities and ages
;
;;;;;;;;;
zm(*)=0.0
age(*)=0.0
hsml(*)=grav_smoothing
if nstars gt 0 then age(star_index) = (Time-agetmp)/HubbleParam
if nstars gt 0 then zm(star_index) = zmtmp(npart(0):npart(0)+npart(4)-1)
if ngas   gt 0 then zm(gas_index) = zmtmp(0:npart(0)-1)
if ngas   gt 0 then hsml(gas_index) = hsmltmp(*)



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make a kernel, w_h, that we can 
; scale and distribute within pixels
; for the pseudo-Monte Carlo 
; integration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
iseed = 3028
N_SAMPLE = 1000
n_s = long(4.0*N_SAMPLE/(!PI))
x_s = fltarr(2,n_s)
x_s(*,*) = 2.0*randomu(iseed,2,n_s)-1.0
i_h = where( sqrt( x_s(0,*)^2 + x_s(1,*)^2)  lt 1.0, n_h)
x_h = fltarr(2,n_h)
w_h = fltarr(n_h)
x_h(*,*) = x_s(*,i_h)
for i=0,n_h-1 do w_h(i) =  kernel_table(float(n_table-1)*sqrt(  x_h(0,i)^2 + x_h(1,i)^2  ))





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;       Find cold_fractions
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

tstar = 4.5
tSN   = 3.0e8
fEVP  = 3000.0
fSN   = 0.1
comoving = 0

print,tstar,tSN,fEVP,fSN,hubbleparam,time,comoving

if(ngas gt 0) then begin

        gas_temperature = fltarr(Ngas)  ;gas tmp
        cold_volume_fraction = fltarr(Ngas)     ;cold_volume_fraction
        cold_mass_fraction = fltarr(Ngas)       ;cold_mass_fraction

        S = CALL_EXTERNAL('/home/brant/code/idl/gas_temperature/gas_temperature', $
                'gas_temperature', $
                Ngas, $
                tstar, $
                TSN, $
                fEVP, $
                fSN, $
                hubbleparam, $
                time, $
                comoving, $
                rho, $
                nel, $
                u, $
                gas_temperature, $
                cold_volume_fraction,$
                cold_mass_fraction,/F_VALUE)


endif 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;       Find luminosities
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

nlum = ngas
if(Nlum gt 0) then begin

        soft_xray_lum = fltarr(Nlum)
        hard_xray_lum = fltarr(Nlum)
        medium_xray_lum = fltarr(Nlum)
        xray_lum = fltarr(Nlum)
        soft_xray_lum(*) = 0.0
        hard_xray_lum(*) = 0.0
        medium_xray_lum(*) = 0.0
        xray_lum(*) = 0.0
        ;soft_xray_lum(*) = 1.0
        ;hard_xray_lum(*) = 1.0
        ;medium_xray_lum(*) = 1.0
        ;xray_lum(*) = 1.0
        rho_nh  = fltarr(Nlum)
        hot_nel = fltarr(Nlum)
        hot_rnh = fltarr(Nlum)
        hot_mass   = fltarr(Nlum)
        z_gas_solar = fltarr(Nlum)

        rho_nh = rho(0:nlum-1)*UnitDensity_in_cgs*0.76/(1.67e-24)
        rho_nh = rho_nh(*)*hubbleparam*hubbleparam*(1.0-cold_mass_fraction(0:nlum-1))/(1.0-cold_volume_fraction(0:nlum-1))
        z_gas_solar(*) = zm(gas_index(0:nlum-1))/0.02
        isupersolar = where(z_gas_solar(*) lt 1.0e-6,nss)
        z_gas_solar(isupersolar) = 1.0e-6
        hot_rnh(*) = rho_nh(0:nlum-1)
        hot_nel(*) = rho_nh(*)*nel(0:nlum-1)
        hot_mass(*) = (1.0-cold_mass_fraction(0:nlum-1))*mass(gas_index(0:nlum-1))

        S = CALL_EXTERNAL('/home/brant/code/idl/raymond_smith_3_color/raymond_smith_3_color', $
                'raymond_smith', $
                Nlum, $
                redshift, $
                z_gas_solar, $
                hot_nel, $
                hot_rnh, $
                hot_mass, $
                gas_temperature(0:nlum-1), $
                hsml(0:nlum-1), $
                soft_xray_lum, $
                medium_xray_lum, $
                hard_xray_lum, $
                /F_VALUE)
	xray_lum(*) = soft_xray_lum(*)+medium_xray_lum(*)+hard_xray_lum(*)

        ;LUMINOSITIES are in h^-1 units

endif else begin
        print,'No gas, no gaseous x-ray luminosities!'
endelse

openw,2,'xray_lums.'+exts+'.txt'
printf,2,Nlum
for ni =0,Nlum-1 do begin
	printf,2,soft_xray_lum(ni),medium_xray_lum(ni),hard_xray_lum(ni)
endfor
close,2

;;;;;;;;;;;;;
;show kernel
;;;;;;;;;;;;;
;nkt = 400L
;ktest = fltarr(nkt,nkt)
;ktest(*,*) = 0.0
;for i=0,nkt-1 do begin
;	for j=0,nkt-1 do begin
;		r =sqrt( (float(i) - float(nkt)/2.0)^2 + (float(j) - float(nkt)/2.0)^2 )/(float(nkt)/2.0)
;		ktest(i,j) = kernel_2d(r)
;	endfor
;endfor
;window,2,xsize=nkt,ysize=nkt
;cols = 255.0
;ktmax = max(ktest)
;ktmin = ktmax/dynrange
;min_index = where(ktest lt ktmin,nmin)
;if nmin gt 0 then ktest(min_index)=ktmin
;image = ((alog10(ktest)-alog10(ktmin))/alog10(ktmax/ktmin))*(cols-254.9)+254.9
;tvscl,alog10(image)


;;;;;;;;;;;;;;;;;;;;;
; SELECT PARTICLES 
; TO SMOOTH!!!!!!!!!!
;;;;;;;;;;;;;;;;;;;;;

n_smooth = ngas
i_smooth = gas_index


;;;;;;;;;
;reposition on BH or cm
;;;;;;;;;

cm = fltarr(3)
if nbh gt 0 then begin
	for i=0,2 do cm(i)  = total(mass(bh_index)*pos(i,bh_index))/total(mass(bh_index))
endif else begin
	for i=0,2 do cm(i)  = total(mass(i_smooth)*pos(i,i_smooth))/total(mass(i_smooth))
endelse
print,'cm : ',cm
for i=0,2 do pos(i,*)  = pos(i,*)-cm(i)

;;;;;;;;;;;;;;;;;;;;;
; bin particles
;;;;;;;;;;;;;;;;;;;;;

l_x = fltarr(n_smooth)
u_x = fltarr(n_smooth)
l_y = fltarr(n_smooth)
u_y = fltarr(n_smooth)
i_x = lonarr(n_smooth)
i_y = lonarr(n_smooth)
t_x = lonarr(n_smooth)
t_y = lonarr(n_smooth)
edge_flag = lonarr(n_smooth)
edge_flag(*) = 0
i=0l
while i lt n_smooth do  begin
	l_x(i) = pos(AXIS0,i_smooth(i)) - hsml(i_smooth(i))  - xmin 
	l_y(i) = pos(AXIS1,i_smooth(i)) - hsml(i_smooth(i))  - ymin 
	u_x(i) = pos(AXIS0,i_smooth(i)) + hsml(i_smooth(i))  - xmin 
	u_y(i) = pos(AXIS1,i_smooth(i)) + hsml(i_smooth(i))  - ymin 
	i_x(i) = long(float(XPIXELS)*l_x(i)/lx)
	i_y(i) = long(float(YPIXELS)*l_y(i)/ly)
	t_x(i) = long(float(XPIXELS)*u_x(i)/lx)
	t_y(i) = long(float(YPIXELS)*u_y(i)/ly)
	i = i+1
endwhile


;;;;;;;;;;;;;;;;;;;
;
; Limit smoothing to
; particle kernel -
; image overlap
;
;;;;;;;;;;;;;;;;;;;

i_plot = where(  (t_x ge 0 and i_x le XPIXELS-1) and (t_y ge 0 and i_y  le YPIXELS-1), n_plot)
print,'Particles contributing to image: ',n_plot
p = where(t_x(i_plot(*)) ge XPIXELS, nb)
if nb gt 0 then begin
	t_x(i_plot(p(*))) = XPIXELS-1
	edge_flag(i_plot(p(*))) = 1
endif
p = where(i_x(i_plot(*)) lt 0, nb)
if nb gt 0 then begin
	i_x(i_plot(p(*))) = 0
	edge_flag(i_plot(p(*))) = 1
endif
p = where(t_y(i_plot(*)) ge YPIXELS, nb)
if nb gt 0 then begin
	t_y(i_plot(p(*))) = YPIXELS-1
	edge_flag(i_plot(p(*))) = 1
endif

p = where(i_y(i_plot(*)) lt 0, nb)
if nb gt 0 then begin
	i_y(i_plot(p(*))) = 0
	edge_flag(i_plot(p(*))) = 1
endif


;window,3,xsize=400,ysize=400
;plot,pos(0,i_smooth(i_plot)),pos(1,i_smooth(i_plot)),psym=3,xrange=[xmin,xmax],xstyle=1,yrange=[ymin,ymax],ystyle=1
bh_pos = fltarr(3,nbh)

for j=0,2 do bh_pos(j,*) = pos(j,bh_index)
;oplot,bh_pos(0,*),bh_pos(1,*),psym=1 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SMOOTHING ROUTINE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


i=0L
x_s_k = 0.0
y_s_k = 0.0
x_a_k = fltarr(n_h)
y_a_k = fltarr(n_h)
r_s_k = 0.0
h_s_k = 0.0
w_s_k = 0.0
while i lt n_plot do begin
	lum_soft   = soft_xray_lum(i_plot(i))
	lum_medium = medium_xray_lum(i_plot(i))
	lum_hard   = hard_xray_lum(i_plot(i))
	lum =mass(i_smooth(i_plot(i)))
	h_s_k = hsml(i_smooth(i_plot(i)))
	if( (i eq 20000) or (i eq 50000) or (i eq 1) or (i eq 10000) or (i eq 12000) or (i eq 5000) or (i eq 0) or (i eq 1000) or (i eq 2000) or (i eq 500) or (i eq 250) or (i eq 100) or (i eq 10)) then print,i,t_x(i_plot(i))-i_x(i_plot(i)),t_y(i_plot(i))-i_y(i_plot(i)),hsml(i_smooth(i_plot(i)))
	if(t_x(i_plot(i))-i_x(i_plot(i)) le 4) then begin
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; SMOOTHING ROUTINE - psuedo-Monte Carlo
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		x_a_k(*) = h_s_k*x_h(0,*)+pos(AXIS0,i_smooth(i_plot(i)))-xmin
		y_a_k(*) = h_s_k*x_h(1,*)+pos(AXIS1,i_smooth(i_plot(i)))-ymin
		for j=i_x(i_plot(i)),t_x(i_plot(i)) do begin
			for k=i_y(i_plot(i)),t_y(i_plot(i)) do begin
				i_w = where( (x_a_k ge p_x(j)) and (x_a_k lt p_x(j)+dx) and (y_a_k ge p_y(k)) and (y_a_k lt p_y(k)+dy), n_i_w)
				if n_i_w gt 0 then begin
					gsum(j,k) = gsum(j,k)+!PI*lum*(total(w_h(i_w))/float(n_h))/(dx*dy)
					gsum_soft(j,k) = gsum_soft(j,k)+!PI*lum_soft*(total(w_h(i_w))/float(n_h))/(dx*dy)
					gsum_hard(j,k) = gsum_hard(j,k)+!PI*lum_hard*(total(w_h(i_w))/float(n_h))/(dx*dy)
					gsum_medium(j,k) = gsum_medium(j,k)+!PI*lum_medium*(total(w_h(i_w))/float(n_h))/(dx*dy)
					if edge_flag(i_plot(i)) eq 0 then begin
						;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
						; particle doesn't over lap with edge
						;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
						w_s_k = w_s_k+!PI*(total(w_h(i_w))/float(n_h))
					endif else begin
						;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
						; particle overlaps with edge, no correction
						;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
						w_s_k = 1.0
					endelse
				endif
			endfor
		endfor
	endif else begin
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; SMOOTHING ROUTINE - summation normalization
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		for j=i_x(i_plot(i)),t_x(i_plot(i)) do begin
			for k=i_y(i_plot(i)),t_y(i_plot(i)) do begin
				x_s_k = pos(AXIS0,i_smooth(i_plot(i)))-(p_x(j)+dx/2.0) - xmin
				y_s_k = pos(AXIS1,i_smooth(i_plot(i)))-(p_y(k)+dy/2.0) - ymin
				r_s_k = sqrt(x_s_k^2 +  y_s_k^2)/h_s_k
				ui = long(float(n_table-1)*r_s_k)
				if ui le (n_table-1)  then begin
					gsum(j,k) =  gsum(j,k)+lum*(kernel_table(ui)/(h_s_k*h_s_k))
					gsum_soft(j,k) =  gsum_soft(j,k)+lum_soft*(kernel_table(ui)/(h_s_k*h_s_k))
					gsum_hard(j,k) =  gsum_hard(j,k)+lum_hard*(kernel_table(ui)/(h_s_k*h_s_k))
					gsum_medium(j,k) =  gsum_medium(j,k)+lum_medium*(kernel_table(ui)/(h_s_k*h_s_k))
					if edge_flag(i_plot(i)) eq 0 then begin
						;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
						; particle doesn't over lap with edge
						;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
						w_s_k = w_s_k + kernel_table(ui)*dx*dy/(h_s_k*h_s_k)
					endif else begin
						;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
						; particle overlaps with edge, no correction
						;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
						w_s_k = 1.0
					endelse
				endif
			endfor
		endfor
	endelse	

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; SMOOTHING ROUTINE - normalize particle SD
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if w_s_k eq 0.0 then w_s_k = 1.0
	gsum(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i))) = gsum(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i)))/w_s_k
	gsum_soft(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i))) = gsum_soft(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i)))/w_s_k
	gsum_hard(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i))) = gsum_hard(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i)))/w_s_k
	gsum_medium(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i))) = gsum_medium(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i)))/w_s_k

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; SMOOTHING ROUTINE - add particle SD to grid
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	grid(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i))) = grid(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i))) + gsum(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i)))
	grid_soft(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i))) = grid_soft(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i))) + gsum_soft(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i)))
	grid_hard(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i))) = grid_hard(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i))) + gsum_hard(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i)))
	grid_medium(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i))) = grid_medium(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i))) + gsum_medium(i_x(i_plot(i)):t_x(i_plot(i)),i_y(i_plot(i)):t_y(i_plot(i)))
	gsum(*,*) = 0.0
	gsum_medium(*,*) = 0.0
	gsum_hard(*,*) = 0.0
	gsum_soft(*,*) = 0.0
	w_s_k = 0.0

	i= i+1
endwhile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SMOOTHING ROUTINE - DONE!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; normalize image
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ps_width = 6.0
i_width  = 5.0
xi = ((ps_width-i_width)/2.0)/ps_width
xaa = ((ps_width-i_width)/2.0)/ps_width + i_width/ps_width
yi = ((ps_width-i_width)/2.0)/ps_width
yaa = ((ps_width-i_width)/2.0)/ps_width + i_width/ps_width
ya = i_width/ps_width
xa = i_width/ps_width

cols = 255.0
image = grid
imax = max(grid)
imin = imax/DynRange
min_index = where(image lt imin,nmin)
if nmin gt 0 then image(min_index)=imin
image = ((alog10(image)-alog10(imin))/alog10(imax/imin))*(cols-2) +2

set_plot,'ps'
device,filename='gas.eps',/color,bits_per_pixel=8
;device,xsize=5.0,ysize=5.0,/inches
device,xsize=ps_width,ysize=ps_width,/inches
loadct,lctnum
tvlct,r,g,b,/get
v1=[255,0]
v2=[255,0]
v3=[255,0]
tvlct,v1,v2,v3,0
tv, image, /normal, xi, yi, xsize = xa, ysize=ya
plot,[0],[0],xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=1,ystyle=1,/noerase,color=1,position=[xi,yi,xaa,yaa],thick=3,xtitle='h!U-1!N kpc',ytitle='h!U-1!N kpc',font=1
oplot,bh_pos(0,*),bh_pos(1,*),psym=1,color=255,symsize=1.5,thick=1.5; ,position=[xi,yi,xaa,yaa]
device,/close
set_plot,'x'

cols = 255.0
image_soft = grid_soft
imax = max(grid_soft)
imin = imax/DynRange
min_index = where(image_soft lt imin,nmin)
if nmin gt 0 then image_soft(min_index)=imin
image_soft = ((alog10(image_soft)-alog10(imin))/alog10(imax/imin))*(cols-2) +2

set_plot,'ps'
device,filename='soft_xray_sb.eps',/color,bits_per_pixel=8
;device,xsize=5.0,ysize=5.0,/inches
device,xsize=ps_width,ysize=ps_width,/inches
loadct,lctnum
tvlct,r,g,b,/get
v1=[255,0]
v2=[255,0]
v3=[255,0]
tvlct,v1,v2,v3,0
tv, image_soft, /normal, xi, yi, xsize = xa, ysize=ya
plot,[0],[0],xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=1,ystyle=1,/noerase,color=1,position=[xi,yi,xaa,yaa],thick=3,xtitle='h!U-1!N kpc',ytitle='h!U-1!N kpc',font=1
oplot,bh_pos(0,*),bh_pos(1,*),psym=1,color=255,symsize=1.5,thick=1.5; ,position=[xi,yi,xaa,yaa]
device,/close
set_plot,'x'

cols = 255.0
image_medium = grid_medium
imax = max(grid_medium)
imin = imax/DynRange
min_index = where(image_medium lt imin,nmin)
if nmin gt 0 then image_medium(min_index)=imin
image_medium = ((alog10(image_medium)-alog10(imin))/alog10(imax/imin))*(cols-2) +2

set_plot,'ps'
device,filename='medium_xray_sb.eps',/color,bits_per_pixel=8
;device,xsize=5.0,ysize=5.0,/inches
device,xsize=ps_width,ysize=ps_width,/inches
loadct,lctnum
tvlct,r,g,b,/get
v1=[255,0]
v2=[255,0]
v3=[255,0]
tvlct,v1,v2,v3,0
tv, image_medium, /normal, xi, yi, xsize = xa, ysize=ya
plot,[0],[0],xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=1,ystyle=1,/noerase,color=1,position=[xi,yi,xaa,yaa],thick=3,xtitle='h!U-1!N kpc',ytitle='h!U-1!N kpc',font=1
oplot,bh_pos(0,*),bh_pos(1,*),psym=1,color=255,symsize=1.5,thick=1.5; ,position=[xi,yi,xaa,yaa]
device,/close
set_plot,'x'

cols = 255.0
image_hard = grid_hard
imax = max(grid_hard)
imin = imax/DynRange
min_index = where(image_hard lt imin,nmin)
if nmin gt 0 then image_hard(min_index)=imin
image_hard = ((alog10(image_hard)-alog10(imin))/alog10(imax/imin))*(cols-2) +2

set_plot,'ps'
device,filename='hard_xray_sb.eps',/color,bits_per_pixel=8
;device,xsize=5.0,ysize=5.0,/inches
device,xsize=ps_width,ysize=ps_width,/inches
loadct,lctnum
tvlct,r,g,b,/get
v1=[255,0]
v2=[255,0]
v3=[255,0]
tvlct,v1,v2,v3,0
tv, image_hard, /normal, xi, yi, xsize = xa, ysize=ya
plot,[0],[0],xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=1,ystyle=1,/noerase,color=1,position=[xi,yi,xaa,yaa],thick=3,xtitle='h!U-1!N kpc',ytitle='h!U-1!N kpc',font=1
oplot,bh_pos(0,*),bh_pos(1,*),psym=1,color=255,symsize=1.5,thick=1.5; ,position=[xi,yi,xaa,yaa]
device,/close
set_plot,'x'

image = fltarr(3,xpixels,ypixels)
image(0,*,*) = image_soft
image(1,*,*) = image_medium 
image(2,*,*) = image_hard
;tv,image,true=1
;filename='image.'+exts+'.eps'
;set_plot,'ps'
;device,filename=filename,/encapsulated,/color,bits=8
;device,/inches,xsize = xsiz,ysize=ysiz
;loadct,0
;tv, image, /inches, xsize=xres, ysize=yres,true=1
;device,/close
;set_plot,'x' 

set_plot,'ps'
device,filename='xray_map.eps',/color,bits_per_pixel=8
device,xsize=ps_width,ysize=ps_width,/inches
loadct,0
tv, image, /normal, xi, yi, xsize=xa, ysize=ya,true=1
plot,[0],[0],xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=1,ystyle=1,/noerase,color=1,position=[xi,yi,xaa,yaa],thick=3,xtitle='h!U-1!N kpc',ytitle='h!U-1!N kpc',font=1
oplot,bh_pos(0,*),bh_pos(1,*),psym=1,color=255,symsize=1.5,thick=1.5; ,position=[xi,yi,xaa,yaa]
device,/close

end
