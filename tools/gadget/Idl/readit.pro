
f="../../snap_000"

npart=lonarr(6)	
massarr=dblarr(6)
time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedback=0L
npartall=lonarr(6)	
bytesleft= 136
la=intarr(bytesleft/2)


openr,1,f,/f77_unformatted
readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,la
print,npart,massarr
print,time,redshift

NGas=  npart(0)

N=  NGas

pos=fltarr(3,N)
vel=fltarr(3,N)
id=lonarr(N)
readu,1, pos
readu,1, vel
readu,1, id
ind=where((npart gt 0) and (massarr eq 0))
if ind(0) ne -1 then begin
  Nm= total(npart(ind))
  mass=fltarr(Nm)	        ; masses for variable mass particles (usually gas+stars)
  readu,1,mass
endif

u=fltarr(Ngas)
readu,1,u                       ; internal energy per unit mass
rho=fltarr(Ngas)
readu,1,rho                     ; comoving gas density

if flag_sfr gt 0 then begin
  Nelec=fltarr(Ngas)
  readu,1 ,Nelec                ; gas electron abundance relative to hydrogen
  NH0=fltarr(Ngas)
  readu,1 ,NH0                  ; neutral hydrogen abundance relative to hydrogen
endif

hsml=fltarr(Ngas)
readu,1,hsml                    ; SPH smoothing length

close,1

t =  -2147475680L
t =  2323L

ind=where(id eq t)


i = ind(0)


print, hsml(i)

x = Pos(0,i)
y = Pos(1,i)
z = Pos(2,i)

r = sqrt((pos(0,*)-x)^2 + (pos(1,*)-y)^2 + (pos(2,*)-z)^2)

ind =where(r lt 20.0)

plot, r(ind), psym=4

rr= r(ind)


rr=rr(sort(rr))

cc = indgen(n_elements(rr))

plot, rr, cc

end





 









