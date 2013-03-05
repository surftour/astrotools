
Num = 1023

Base= "/gpfs/mpa/vrs/Billennium/halo_C02/C02_400/"
SnapBase="snap_C02_400"

GrNr = 1

;--------------------------------------------------

if num lt 1000 then begin
   exts='000'
   exts=exts+strcompress(string(Num),/remove_all)
   exts=strmid(exts,strlen(exts)-3,3)
endif else begin
   exts='0000'
   exts=exts+strcompress(string(Num),/remove_all)
   exts=strmid(exts,strlen(exts)-4,4)
endelse


skip = 0L
skip_sub = 0L
fnr = 0L

repeat begin

    f = Base + "/groups_" + exts +"/subhalo_tab_"+exts +"."+strcompress(string(fnr),/remove_all)

    Ngroups = 0L
    TotNgroups = 0L
    Nids = 0L
    TotNids = 0LL
    NTask = 0L
    Nsubgroups = 0L
    TotNsubgroups = 0L

    openr,1,f
    readu,1, Ngroups, TotNgroups, Nids, TotNids, NTask, Nsubgroups, TotNsubgroups

    if fnr eq 0 then begin
        GroupLen = lonarr(TotNgroups)
        GroupOffset = lonarr(TotNgroups)
        GroupMass = fltarr(TotNgroups)
        GroupPos = fltarr(3, TotNgroups)
        Group_M_Mean200 = fltarr(TotNgroups)
        Group_R_Mean200 = fltarr(TotNgroups)
        Group_M_Crit200 = fltarr(TotNgroups)
        Group_R_Crit200 = fltarr(TotNgroups)
        Group_M_TopHat200 = fltarr(TotNgroups)
        Group_R_TopHat200 = fltarr(TotNgroups)
        GroupContaminationCount = lonarr(TotNgroups)
        GroupContaminationMass = fltarr(TotNgroups)
        GroupNsubs = lonarr(TotNgroups)
        GroupFirstSub = lonarr(TotNgroups)

        SubhaloLen = lonarr(TotNsubgroups)
        SubhaloOffset = lonarr(TotNsubgroups)
        SubhaloParent = lonarr(TotNsubgroups)
        SubhaloMass = fltarr(TotNsubgroups)
        SubhaloPos = fltarr(3, TotNsubgroups)
        SubhaloVel = fltarr(3, TotNsubgroups)
        SubhaloCM = fltarr(3, TotNsubgroups)
        SubhaloSpin = fltarr(3, TotNsubgroups)
        SubhaloVelDisp = fltarr(TotNsubgroups)
        SubhaloVmax = fltarr(TotNsubgroups)
        SubhaloVmaxRad = fltarr(TotNsubgroups)
        SubhaloHalfmassRad = fltarr(TotNsubgroups)
        SubhaloIDMostbound = lonarr(TotNsubgroups)
        SubhaloGrNr = lonarr(TotNsubgroups)
    endif

    if Ngroups gt 0 then begin
        
        locLen = lonarr(Ngroups)
        locOffset = lonarr(Ngroups)
        locMass = fltarr(Ngroups)
        locPos = fltarr(3, Ngroups)
        loc_M_Mean200 = fltarr(Ngroups)
        loc_R_Mean200 = fltarr(Ngroups)
        loc_M_Crit200 = fltarr(Ngroups)
        loc_R_Crit200 = fltarr(Ngroups)
        loc_M_TopHat200 = fltarr(Ngroups)
        loc_R_TopHat200 = fltarr(Ngroups)
        locContaminationCount = lonarr(Ngroups)
        locContaminationMass = fltarr(Ngroups)
        locNsubs = lonarr(Ngroups)
        locFirstSub = lonarr(Ngroups)

        readu,1, loclen 
        readu,1, locOffset
        readu,1, locMass
        readu,1, locPos
        readu,1, loc_M_Mean200
        readu,1, loc_R_Mean200
        readu,1, loc_M_Crit200
        readu,1, loc_R_Crit200
        readu,1, loc_M_TopHat200
        readu,1, loc_R_TopHat200
        readu,1, locContaminationCount
        readu,1, locContaminationMass
        readu,1, locNsubs
        readu,1, locFirstSub

        GroupLen(skip:skip+Ngroups-1) = locLen(*)
        GroupOffset(skip:skip+Ngroups-1) = locOffset(*)
        GroupMass(skip:skip+Ngroups-1) = locMass(*)
        GroupPos(*, skip:skip+Ngroups-1) = locPos(*,*)
        Group_M_Mean200(skip:skip+Ngroups-1) = loc_M_Mean200(*)
        Group_R_Mean200(skip:skip+Ngroups-1) = loc_R_Mean200(*)
        Group_M_Crit200(skip:skip+Ngroups-1) = loc_M_Crit200(*)
        Group_R_Crit200(skip:skip+Ngroups-1) = loc_R_Crit200(*)
        Group_M_TopHat200(skip:skip+Ngroups-1) = loc_M_TopHat200(*)
        Group_R_TopHat200(skip:skip+Ngroups-1) = loc_R_TopHat200(*)
        GroupContaminationCount(skip:skip+Ngroups-1) = locContaminationCount(*) 
        GroupContaminationMass(skip:skip+Ngroups-1) = locContaminationMass(*)
        GroupNsubs(skip:skip+Ngroups-1) = locNsubs(*)
        GroupFirstsub(skip:skip+Ngroups-1) = locFirstSub(*)

        skip+= Ngroups
    endif


    if Nsubgroups gt 0 then begin
        
        locLen = lonarr(Nsubgroups)
        locOffset = lonarr(Nsubgroups)
        locParent = lonarr(Nsubgroups)
        locMass = fltarr(Nsubgroups)
        locPos = fltarr(3, Nsubgroups)
        locVel = fltarr(3, Nsubgroups)
        locCM = fltarr(3, Nsubgroups)
        locSpin= fltarr(3, Nsubgroups)
        locVelDisp = fltarr(Nsubgroups)
        locVmax = fltarr(Nsubgroups)
        locVmaxRad = fltarr(Nsubgroups)
        locHalfMassRad = fltarr(Nsubgroups)
        locIDMostBound = lonarr(Nsubgroups)
        locGrNr = lonarr(Nsubgroups)

        readu,1, loclen 
        readu,1, locOffset
        readu,1, locParent
        readu,1, locMass
        readu,1, locPos
        readu,1, locVel
        readu,1, locCM
        readu,1, locSpin
        readu,1, locVelDisp
        readu,1, locVmax
        readu,1, locVmaxRad
        readu,1, locHalfMassRad
        readu,1, locIDMostBound
        readu,1, locGrNr

        SubhaloLen(skip_sub:skip_sub+Nsubgroups-1) = locLen(*)
        SubhaloOffset(skip_sub:skip_sub+Nsubgroups-1) = locOffset(*)
        SubhaloParent(skip_sub:skip_sub+Nsubgroups-1) = locParent(*)
        SubhaloMass(skip_sub:skip_sub+Nsubgroups-1) = locMass(*)
        SubhaloPos(*, skip_sub:skip_sub+Nsubgroups-1) = locPos(*,*)
        SubhaloVel(*, skip_sub:skip_sub+Nsubgroups-1) = locVel(*,*)
        SubhaloCM(*, skip_sub:skip_sub+Nsubgroups-1) = locCM(*,*)
        SubhaloSpin(*, skip_sub:skip_sub+Nsubgroups-1) = locSpin(*,*)
        SubhaloVeldisp(skip_sub:skip_sub+Nsubgroups-1) = locVeldisp(*)
        SubhaloVmax(skip_sub:skip_sub+Nsubgroups-1) = locVmax(*)
        SubhaloVmaxRad(skip_sub:skip_sub+Nsubgroups-1) = locVmaxRad(*)
        SubhaloHalfmassRad(skip_sub:skip_sub+Nsubgroups-1) = locHalfmassRad(*)
        SubhaloIDMostBound(skip_sub:skip_sub+Nsubgroups-1) = locIDMostBound(*)
        SubhaloGrNr(skip_sub:skip_sub+Nsubgroups-1) = locGrNr(*)

        skip_sub+= Nsubgroups
    endif

    close, 1

    fnr++

endrep until fnr eq NTask

print
print, "TotNgroups   =", TotNgroups
print, "TotNsubgroups=", TotNsubgroups
print
print, "Largest group of length ", GroupLen(0)," has", GroupNsubs(0)," substructures"
print



skip = 0L
fnr = 0L

repeat begin

    f = Base + "/groups_" + exts +"/subhalo_ids_"+exts +"."+strcompress(string(fnr),/remove_all)

    Ngroups = 0L
    TotNgroups = 0L
    Nids = 0L
    TotNids = 0LL
    NTask = 0L
    Offset = 0L


    openr,1,f
    readu,1, Ngroups, TotNgroups, Nids, TotNids, NTask, Offset

    if fnr eq 0 then begin
        IDs = lonarr(TotNids)
    endif

    if Nids gt 0 then begin
        
        locIDs = lonarr(Nids)
        readu,1, locIDs

        IDs(skip:skip+Nids-1) = locIDs(*)
        skip+= Nids

    endif

    close, 1

    fnr++

endrep until fnr eq NTask

print
print, "Total IDs in grps=", TotNids
print











f= Base + "/snapdir_"+exts+"/"+ snapbase+"_"+exts +".0"
f= Strcompress(f, /remove_all)

falt = Base + "/"+ snapbase+"_"+exts 
falt= Strcompress(falt, /remove_all)


npart=lonarr(6)		
massarr=dblarr(6)
time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedback=0L
npartall=lonarr(6)	
flag_cooling= 0L
Nsubfiles = 0L
BoxSize = 0.0D

openr,1,f,/f77_unformatted, error = err
if err ne 0 then begin
  openr,1,falt,/f77_unformatted, error = err
endif

readu,1, npart, massarr, time, redshift, flag_sfr, flag_feedback, $
         npartall, flag_cooling, Nsubfiles, BoxSize 
close,1


print, npartall
print,time,redshift
print, boxsize


count=0L
for subfile=0, Nsubfiles-1 do begin

    f= Base + "/snapdir_"+exts+"/"+ snapbase+"_"+exts + "."+ string(subfile)
    f= Strcompress(f, /remove_all)
  
    if Nsubfiles le 1 then begin
       f=base + "/"+ snapbase+"_"+exts 
       f=strcompress(f,/remove_all)
    endif

    openr,1,f,/f77_unformatted

    readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall
    print,npart,massarr
    print
    N=  npart(0) + npart(1) + npart(2) + npart(3) + npart(4) + npart(5)
    pos1 = dblarr(3,N)
    vel1 = dblarr(3,N)
    readu,1, pos1
    readu,1, vel1

    id1 = lonarr(N)
    readu,1, ID1

    mass1 = dblarr(N)
    readu,1, mass1

    close,1

    if subfile eq 0 then begin
       Pos= dblarr(3, npartall(0) + npartall(1) + npartall(2) + npartall(3) + npartall(4) + npartall(5))
       Vel= dblarr(3, npartall(0) + npartall(1) + npartall(2) + npartall(3) + npartall(4) + npartall(5))
       ID= lonarr(npartall(0) + npartall(1) + npartall(2) + npartall(3) + npartall(4) + npartall(5))
       mass= dblarr(npartall(0) + npartall(1) + npartall(2) + npartall(3) + npartall(4) + npartall(5))

    endif

    Pos(*, count:count+N-1)= Pos1(*,*)
    Vel(*, count:count+N-1)= Vel1(*,*)
    mass(count:count+N-1)= mass1(*)
    id(count:count+N-1)= id1(*)
    count=count+N
endfor


ind = sort(ID)

idlist = IDs(GroupOffset(GrNr):GroupOffset(GrNr)+GroupLen(GrNr)-1)

ind_halo = sort(idlist)

poslist = fltarr(3, GroupLen(GrNr))


  j = 0L

  for i = 0L, GroupLen(GrNr)-1 do begin
    while ID(ind(j)) lt idlist(ind_halo(i)) do begin
       j++
    endwhile

    poslist(*,ind_halo(i)) = Pos(*, ind(j))
  endfor

window, xsize=600,ysize=600

plot, poslist(0,*) - GroupPos(0,GrNr), poslist(1,*) - GroupPos(1,GrNr), psym=3


r=[255,   0,   0, 255,   0, 255]
g=[  0, 255,   0, 255, 255, 0  ]
b=[  0,   0, 255,   0, 255, 255]

co = r + g*256L + b*256L^2

for snr = 0, GroupNsubs(GrNr)-1 do begin

   start = SubhaloOffset(GroupFirstSub(GrNr)+snr) - GroupOffset(GrNr)
   len = SubhaloLen(GroupFirstSub(GrNr)+snr)

   if snr eq 0 then symnr=3 else symnr = 4

   print, snr, start, len

   oplot, poslist(0,start:start+len-1) - GroupPos(0,GrNr), $
          poslist(1,start:start+len-1) - GroupPos(1,GrNr), psym=symnr, color=co(snr mod n_elements(co))

endfor





end
