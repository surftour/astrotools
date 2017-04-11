
N  = 10000L

fout = "seed_glass.dat"



npart=lonarr(6)	
massarr=dblarr(6)
time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedback=0L
npartall=lonarr(6)	
flag_cooling= 0L
num_files= 1L  
bytesleft=120
la=intarr(bytesleft/2)


BoxSize = 100.0D * 1000.0D ; 100 Mpc/h box

npart(1) = N
npartall(1) = N


UnitLength_in_cm =        3.085678d21        ;  1.0 kpc
UnitMass_in_g    =        1.989d43 ;  1.0e10 solar masses
UnitVelocity_in_cm_per_s = 1d5 ;  1 km/sec
UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s
GRAVITY   = 6.672d-8
G=GRAVITY/ UnitLength_in_cm^3 * UnitMass_in_g * UnitTime_in_s^2
HUBBLE    =   3.2407789d-18	

H = HUBBLE * UnitTime_in_s

Omega = 1.0D

mm = Omega * 3 * H * H /(8*!DPI*G) * BoxSize^3 / N

massarr(1)=  mm

pos= fltarr(3, N)

for n=0L, N-1 do begin
    pos(0, n)  = randomu(seed) * BoxSize
    pos(1, n)  = randomu(seed) * BoxSize
    pos(2, n)  = randomu(seed) * BoxSize
endfor

vel= fltarr(3, N)
id = lonarr(N)


openw,1,fout,/f77_unformatted
writeu,1, npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,flag_cooling,num_files,BoxSize,la
writeu,1, pos
writeu,1, vel
writeu,1, id
close,1

end


