
Num = 127

Base= "/gpfs/mpa/vrs/Billennium/halo_C02/C02_200/"
SnapBase="snap_C02_200"

SubNr = 0  ; select the subhalo number

LONGIDs = 0   ; set to 1 for 64bit IDs


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
        if LONGIDs then begin
           SubhaloIDMostbound = lon64arr(TotNsubgroups)
        endif else begin
           SubhaloIDMostbound = lonarr(TotNsubgroups)
        endelse
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
        if LONGIDs then begin
           locIDMostBound = lon64arr(Nsubgroups)
        endif else begin
           locIDMostBound = lonarr(Nsubgroups)
        endelse
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





f= Base + "/snapdir_"+exts+"/"+ snapbase+"_subidorder_"+exts +".0"
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





len = SubhaloLen(SubNr)
offset = SubhaloOffset(SubNr)


print, "desired substructure has length=", len

partmass = SubhaloMass(0) / SubhaloLen(0)


Pos = dblarr(3, len)
Vel = dblarr(3, len)
if LONGIDs then begin
   ID  = lon64arr(len)
endif else begin
   ID  = lonarr(len)
endelse

Hsml = fltarr(len)
Dens = fltarr(len)
VelDisp = fltarr(len)


found=0L

count=0LL
for subfile=0, Nsubfiles-1 do begin

    f= Base + "/snapdir_"+exts+"/"+ snapbase+"_subidorder_"+exts + "."+ string(subfile)
    f= Strcompress(f, /remove_all)
  
    if Nsubfiles le 1 then begin
       f=base + "/"+ snapbase+"_"+exts 
       f=strcompress(f,/remove_all)
    endif

    print, f

    dummy = 0L
    dummy2 = 0L

    openr,1,f

    readu,1,dummy
    readu,1,npart, massarr
    skip_lun, 1, dummy - 6*(4L+8L)
    readu,1,dummy2
    print,dummy, dummy2
    


    if offset lt (count + npart(1)) then begin
       
       nskip = offset - count
       nrem  =  count + npart(1) - offset
       if len gt nrem then begin
          n_to_read = nrem
       endif else begin
          n_to_read = len
       endelse


       if n_to_read gt 0 then begin

          tmp_pos = dblarr(3, n_to_read)
          tmp_vel = dblarr(3, n_to_read)
          if LONGIDs then begin
             tmp_id  = lon64arr(n_to_read)
          endif else begin
             tmp_id  = lonarr(n_to_read)
          endelse
          tmp_hsml = fltarr(n_to_read)
          tmp_dens = fltarr(n_to_read)
          tmp_veldisp = fltarr(n_to_read)


          readu,1,dummy
          skip_lun, 1, 3*8L * nskip 
          readu,1, tmp_pos
          skip_lun, 1, dummy - 3*8L * (nskip+n_to_read) 
          readu,1,dummy2
          print,dummy, dummy2

          readu,1,dummy
          skip_lun, 1, 3*8L * nskip 
          readu,1, tmp_vel
          skip_lun, 1, dummy - 3*8L * (nskip+n_to_read) 
          readu,1,dummy2
          print,dummy, dummy2

          if LONGIDs then size = 8L else size=4L

          readu,1,dummy
          skip_lun, 1, size * nskip 
          readu,1, tmp_id
          skip_lun, 1, dummy - size * (nskip+n_to_read) 
          readu,1,dummy2
          print,dummy, dummy2


          ; if there are masses stored, skip them
          ind = where((npart gt 0) and (massarr eq 0))
          if ind(0) ne -1 then begin
             readu,1,dummy
             skip_lun, 1, dummy
             readu,1,dummy2
             print, dummy, dummy2
          endif

          ; now hsml, density and veldisp

          size = 4L

          readu,1,dummy
          skip_lun, 1, size * nskip 
          readu,1, tmp_hsml
          skip_lun, 1, dummy - size * (nskip+n_to_read) 
          readu,1,dummy2
          print,dummy, dummy2

          readu,1,dummy
          skip_lun, 1, size * nskip 
          readu,1, tmp_dens
          skip_lun, 1, dummy - size * (nskip+n_to_read) 
          readu,1,dummy2
          print,dummy, dummy2

          readu,1,dummy
          skip_lun, 1, size * nskip 
          readu,1, tmp_veldisp
          skip_lun, 1, dummy - size * (nskip+n_to_read) 
          readu,1,dummy2
          print,dummy, dummy2


          Pos(*, found:found+n_to_read-1) = tmp_pos(*,*)
          Vel(*, found:found+n_to_read-1) = tmp_vel(*,*)
          ID(found:found+n_to_read-1) = tmp_id(*)

          Hsml(found:found+n_to_read-1) = tmp_hsml(*)
          Dens(found:found+n_to_read-1) = tmp_dens(*)
          VelDisp(found:found+n_to_read-1) = tmp_veldisp(*)

          found += n_to_read

       endif

       offset += n_to_read
       len -= n_to_read

    endif

    close,1

    count += npart(1)

endfor



plot, Pos(0,*) - SubhaloPos(0,SubNr), Pos(1,*) - SubhaloPos(1,SubNr), psym = 3

;r = sqrt( (Pos(0,*) - SubhaloPos(0,SubNr))^2 + (Pos(1,*) - SubhaloPos(1,SubNr))^2 + (Pos(2,*) - SubhaloPos(2,SubNr))^2)

;window, xsize=600,ysize=600

;plot, poslist(0,*) - GroupPos(0,GrNr), poslist(1,*) - GroupPos(1,GrNr), psym=3


;r=[255,   0,   0, 255,   0, 255]
;g=[  0, 255,   0, 255, 255, 0  ]
;b=[  0,   0, 255,   0, 255, 255]


end
