;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;	colors.bh.pro
;
;	Reads in snapshots w/ black holes
;	Assigns luminosities colors to 
;	stellar particles
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Set some x device options
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
set_plot,'x'
device,decompose=0,true_color=24,retain=2


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Read in snapshot
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fdir = "/home/tcox/Data/aA1/"
snapbase = "snap"


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

;Luminosity index legend
; 0 = Sloan u
; 1 = Sloan g
; 2 = Sloan r
; 3 = Sloan i
; 4 = Sloan z
; 5 = Cousins J
; 6 = Cousins H
; 7 = Cousins K

s_UBVRIJHK = fltarr(8)
s_UBVRIJHK(0) = 5.66; //U
s_UBVRIJHK(1) = 5.47; //B
s_UBVRIJHK(2) = 4.82; //V
s_UBVRIJHK(3) = 4.28; //R
s_UBVRIJHK(4) = 3.94; //I
s_UBVRIJHK(5) = 3.64; //J ?
s_UBVRIJHK(6) = 3.44; //H ?
s_UBVRIJHK(7) = 3.33; //K

s_ugrizJHK = fltarr(8)
s_ugrizJHK(0) = 6.2789; //u
s_ugrizJHK(1) = 4.9489; //g
s_ugrizJHK(2) = 4.44964; //r
s_ugrizJHK(3) = 4.34644; //i
s_ugrizJHK(4) = 4.3592; //z
s_ugrizJHK(5) = 3.64; //J ?
s_ugrizJHK(6) = 3.44; //H ?
s_ugrizJHK(7) = 3.33; //K


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


print,'npart',npart
print,'massaray',massarr
print,"a = ", time
print,"z = ",redshift
print,"flag_sfr = ",flag_sfr
print,"flag_feedback= ",flag_feedback
print,"npartTotal = ",npartTotal
print,"flag_cooling= ",flag_cooling
print,"numfiles = ",num_files
print,'BoxSize= ',BoxSize
print,'Omega_m = ',Omega0
print,'Omega_L = ',OmegaLambda
print,'h = ',HubbleParam
print,"flag_multiphase= ",flag_multiphase
print,"flag_stellarage= ",flag_stellarage

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
gsbhmass	=fltarr(npart(0)+npart(4)+npart(5))	;gas, stellar, and black hole masses
u		=fltarr(npart(0))			;gas internal energy
rho		=fltarr(npart(0))			;gas density
nel		=fltarr(npart(0))			;gas electron density
nh		=fltarr(npart(0))			;gas hydrogen density
hsml		=fltarr(npart(0))			;gas smoothing length
sfr		=fltarr(npart(0))			;gas star formation rate

if npart(4) gt 0 then begin				;if stellar particles
age		=fltarr(npart(4))			;stellar ages
endif

zm		=fltarr(npart(0)+npart(4))		;gas and stellar metallicities
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
readu,1,gsbhmass
readu,1,u
readu,1,rho
readu,1,nel
readu,1,nh
readu,1,hsml
readu,1,sfr
if npart(4) gt 0 then readu,1,age
readu,1,zm
close,1


;;;;;;;;;;;;;;;;;;;;;;;;;;
; set particle masses
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;
;gas mass
;;;;;;;;;;;;;;;;;;;;;;;;;;

mass(0:npart(0)-1) = gsbhmass(0:npart(0)-1)


;;;;;;;;;;;;;;;;;;;;;;;;;;
;hrdm mass
;;;;;;;;;;;;;;;;;;;;;;;;;;
mass(npart(0):npart(0)+npart(1)-1) = float(massarr(1))

;;;;;;;;;;;;;;;;;;;;;;;;;;
;lrdm mass
;;;;;;;;;;;;;;;;;;;;;;;;;;
if npart(2) gt 0 then mass(npart(0)+npart(1):npart(0)+npart(1)+npart(2)-1) = float(massarr(2))

;;;;;;;;;;;;;;;;;;;;;;;;;;
;bldm mass
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
types(0:npart(0)-1)=0
types(npart(0):npart(0)+npart(1)-1)=1
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
lrdm_index = where(types eq 2,Ndisk)
bldm_index = where(types eq 3,Nbulge)
stars_index = where(types eq 4,Nstars)
bh_index   = where(types eq 5,Nbh)


print,'Ngas   :',Ngas
print,'Nhalo  :',Nhalo
print,'Ndisk  :',Ndisk
print,'Nbulge :',Nbulge
print,'Nstars :',Nstars
print,'Nbh    :',Nbh


new_mass= mass(stars_index)
new_z= zm(npart(0):npart(0)+npart(4)-1)

if(nstars gt 0) then begin

	luminosities	= fltarr(8,Nstars)	;luminosities
	
	S = CALL_EXTERNAL('colors', $
		'colors', $
		Nstars, $
		redshift, $
		new_mass, $
		age, $
		new_z, $
		luminosities)
endif else begin
	print,'No stars, no luminosities!'
endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Plot u-band luminosities for fun
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

window,3,xsize=400,ysize=400
plot,age,luminosities(0,*),psym=1,/ylog,ytitle='L!Lu,sun!N',xtitle='Age [Gyr]',charsize=1.3

end
