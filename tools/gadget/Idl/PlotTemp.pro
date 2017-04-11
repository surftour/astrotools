
Base=   "~/bfd/mcdonald/run128/"

SnapBase="snap"

Num= 10



UnitLength_in_cm         =        3.085678d21 ;  1.0 kpc 
UnitMass_in_g            =        1.989d43 ;  1.0e10 solar masses 
UnitVelocity_in_cm_per_s =        1d5 ;  1 km/sec 
HubbleParam              =        0.7

GRAVITY    = 6.672d-8
BOLTZMANN  = 1.3806d-16
PROTONMASS = 1.6726e-24

UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s 
UnitDensity_in_cgs= UnitMass_in_g/ UnitLength_in_cm^3
UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ UnitTime_in_s^2
UnitEnergy_in_cgs= UnitMass_in_g * UnitLength_in_cm^2 / UnitTime_in_s^2
G=GRAVITY/ UnitLength_in_cm^3 * UnitMass_in_g * UnitTime_in_s^2

Xh=0.76D                        ; mass fraction of hydrogen

Gamma= 5.0/3



exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)
f=base + snapbase+"_"+exts
f=strcompress(f,/remove_all)


npart=lonarr(6)	
massarr=dblarr(6)
time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedback=0L
npartall=lonarr(6)	
bytesleft= 96
filler=intarr(bytesleft/2)


openr,1,f ,/f77_unformatted
readu,1, npart, massarr, time, redshift, flag_sfr, flag_feedback, $
  npartall, flag_cooling, num_files, BoxSize, Omega0, OmegaLambda, HubbleParam, filler
flag_cooling =0L
num_files = 0L
BoxSize = 0.0D
Omega0 = 0.0D
OmegaLambda = 0.0D
HubbleParam = 0.0D


print, npart, massarr
print, "Time= ", time
print, "Redshift= ", redshift
N=total(npart)
pos=fltarr(3,N)
vel=fltarr(3,N)
id=lonarr(N)
readu, 1, pos
readu, 1, vel
readu, 1, id


ind=where((npart gt 0) and (massarr eq 0))
if ind(0) ne -1 then begin
    Nm= total(npart(ind))
    mass=fltarr(Nm)	
    readu,1,mass
endif

NGas=  npart(0)
NDM =  npart(1)
NStars=npart(4)

if Ngas gt 0 then begin
    u=fltarr(Ngas)
    readu,1,u
    rho=fltarr(Ngas)
    readu,1,rho
    Nelec=fltarr(Ngas)
    readu,1 ,Nelec
    NH0=fltarr(Ngas)
    readu,1 ,NH0
    hsml=fltarr(Ngas)
    readu,1,hsml
endif
close,1

;; convert velocities to peculiar velocity

Vel = Vel * sqrt(Time)

;; extract coordinates and velocities of gas and dark matter particles

if NGas gt 0 then begin
    xgas=pos(0,0:Ngas-1)
    ygas=pos(1,0:Ngas-1)
    zgas=pos(2,0:Ngas-1)
    vxgas=vel(0,0:Ngas-1)
    vygas=vel(1,0:Ngas-1)
    vzgas=vel(2,0:Ngas-1)
endif

if Ndm gt 0 then begin
    xdm= pos(0,0+Ngas:Ndm+Ngas-1)
    ydm= pos(1,0+Ngas:Ndm+Ngas-1)
    zdm= pos(2,0+Ngas:Ndm+Ngas-1)
    vxdm= vel(0,0+Ngas:Ndm+Ngas-1)
    vydm= vel(1,0+Ngas:Ndm+Ngas-1)
    vzdm= vel(2,0+Ngas:Ndm+Ngas-1)
endif



;; Compute temperature in Kelvin

MeanWeight= 4.0/(3*Xh+1+4*Xh*Nelec) * PROTONMASS
Temp = MeanWeight/BOLTZMANN * (gamma-1) * U * UnitEnergy_in_cgs/ UnitMass_in_g

;; Compute mean density

rho_mean_gas= 0.04 * 3*0.1^2 /(8*!DPI*G) 



ind=where(randomu(seed,Ngas) lt 0.1)

plot,  rho(ind)/rho_mean_gas, temp(ind),  psym=3  , /xlog, /ylog


end

















