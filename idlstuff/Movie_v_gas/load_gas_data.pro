pro load_gas_data,  Num, Path, Snapbase, N, CurrentTime, Pos, Vel, Hsml, Mass, Temp, Id


UnitLength_in_cm =        3.085678d21 ;  1.0 kpc 
UnitMass_in_g    =        1.989d43 ;  1.0e10 solar masses 
UnitVelocity_in_cm_per_s = 1d5  ;  1 km/sec 
UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s 
UnitDensity_in_cgs= UnitMass_in_g/ UnitLength_in_cm^3
UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ UnitTime_in_s^2
UnitEnergy_in_cgs= UnitMass_in_g * UnitLength_in_cm^2 / UnitTime_in_s^2
GRAVITY   = 6.672d-8
BOLTZMANN = 1.3806d-16
PROTONMASS = 1.6726e-24
G=GRAVITY/ UnitLength_in_cm^3 * UnitMass_in_g * UnitTime_in_s^2
Xh=0.76D                        ; mass fraction of hydrogen
HubbleParam= 0.65
gamma= 5.0/3



exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)
f= path + snapbase+ "_"+exts
f=strcompress(f,/remove_all)

npart=lonarr(6)	
massarr=dblarr(6)
time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedback=0L
npartall=lonarr(6)	
bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4
la=intarr(bytesleft/2)


openr,1,f,/f77_unformatted
readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,la

CurrentTime = Time

N= npart(0)

Pos= fltarr(3,N)
Vel= fltarr(3,N)
Id = lonarr(N)
Mass=fltarr(N)	

Hsml=fltarr(N)
U=fltarr(N)
Nelec=fltarr(N)


readu,1, Pos
readu,1, Vel
readu,1, Id
readu,1, Mass


readu,1, U

readu,1  ; ,rho
readu,1, Nelec
readu,1  ;,NH0

readu,1, Hsml

readu,1  ;; , Sfr

close,1


MeanWeight= 4.0/(3*Xh+1+4*Xh*Nelec) * PROTONMASS


Temp = MeanWeight/BOLTZMANN * (gamma-1) * U * UnitEnergy_in_cgs/ UnitMass_in_g


end
