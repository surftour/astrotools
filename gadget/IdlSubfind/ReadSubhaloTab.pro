

FLAG_Group_VelDisp = 0  ; Set this to one if the SO-properties computed by SUBFIND for
                        ; the FOF halos contain velocity dispersions


Num = 1023

Base= "/ptmp/vrs/Billennium/C02_400/"


;--------------------------------------------------

if num ge 1000 then begin
   exts='0000'
   exts=exts+strcompress(string(Num),/remove_all)
   exts=strmid(exts,strlen(exts)-4,4)
endif else begin
   exts='000'
   exts=exts+strcompress(string(Num),/remove_all)
   exts=strmid(exts,strlen(exts)-3,3)
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
        Group_VelDisp_Mean200 = fltarr(TotNgroups)
        Group_VelDisp_Crit200 = fltarr(TotNgroups)
        Group_VelDisp_TopHat200 = fltarr(TotNgroups)
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
        loc_VelDisp_Mean200 = fltarr(Ngroups)
        loc_VelDisp_Crit200 = fltarr(Ngroups)
        loc_VelDisp_TopHat200 = fltarr(Ngroups)
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
        if FLAG_Group_VelDisp ne 0 then begin
          readu,1, loc_VelDisp_Mean200
          readu,1, loc_VelDisp_Crit200
          readu,1, loc_VelDisp_TopHat200
        endif
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
        Group_VelDisp_Mean200(skip:skip+Ngroups-1) = loc_VelDisp_Mean200(*)
        Group_VelDisp_Crit200(skip:skip+Ngroups-1) = loc_VelDisp_Crit200(*)
        Group_VelDisp_TopHat200(skip:skip+Ngroups-1) = loc_VelDisp_TopHat200(*)
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

end
