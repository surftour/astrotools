
Num = 0  

Base= "/afs/mpa/home/snhess/tmp/P-Gadget3/SBcluster"

sfrflag = 0
bhflag = 0

;--------------------------------------------------

exts='000'
exts=exts+strcompress(string(Num),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)


skip = 0L
fnr = 0L

repeat begin

    f = Base + "/groups_" + exts +"/group_tab_"+exts +"."+strcompress(string(fnr),/remove_all)

    Ngroups = 0L
    TotNgroups = 0L
    Nids = 0L
    TotNids = 0LL
    NTask = 0L


    openr,1,f
    readu,1, Ngroups, TotNgroups, Nids, TotNids, NTask

    if fnr eq 0 then begin
        GroupLen = lonarr(TotNgroups)
        GroupOffset = lonarr(TotNgroups)
        GroupMass = fltarr(TotNgroups)
        GroupCM = fltarr(3, TotNgroups)
        GroupVel= fltarr(3, TotNgroups)
        GroupLenType = lonarr(6, TotNgroups)
        GroupLenMass = fltarr(6, TotNgroups)
        if sfrflag then begin
            GroupSfr = fltarr(TotNgroups)
        endif
        if bhflag then begin
            GroupBHMass = fltarr(TotNgroups)
            GroupBHMdot = fltarr(TotNgroups)
        endif
    endif

    if Ngroups gt 0 then begin
        
        locLen = lonarr(Ngroups)
        locOffset = lonarr(Ngroups)
        locMass = fltarr(Ngroups)
        locCM = fltarr(3, Ngroups)
        locVel= fltarr(3, Ngroups)
        locLenType = lonarr(6, Ngroups)
        locLenMass = fltarr(6, Ngroups)
        if sfrflag then begin
            locSfr = fltarr(Ngroups)
        endif
        if bhflag then begin
            locBHMass = fltarr(Ngroups)
            locBHMdot = fltarr(Ngroups)
        endif

        readu,1, loclen 
        readu,1, locOffset
        readu,1, locMass
        readu,1, locCM
        readu,1, locVel
        readu,1, locLenType
        readu,1, locLenMass
        if sfrflag then begin
            readu,1, locSfr
        endif
        if bhflag then begin
            readu,1, locBHMass
            readu,1, locBHMdot
        endif

        GroupLen(skip:skip+Ngroups-1) = locLen(*)

        GroupOffset(skip:skip+Ngroups-1) = locOffset(*)
        GroupMass(skip:skip+Ngroups-1) = locMass(*)
        GroupCM(*, skip:skip+Ngroups-1) = locCM(*,*)
        GroupVel(*, skip:skip+Ngroups-1) = locVel(*,*)
        GroupLenType(*, skip:skip+Ngroups-1) = locLenType(*,*)
        GroupLenMass(*, skip:skip+Ngroups-1) = locLenMass(*,*)
        if sfrflag then begin
            GroupSfr(skip:skip+Ngroups-1) = locSfr(*)
        endif
        if bhflag then begin
            GroupBHMass(skip:skip+Ngroups-1) = locBHMass(*)
            GroupBHMdot(skip:skip+Ngroups-1) = locBHMdot(*)
        endif
        
        skip+= Ngroups

    endif

    close, 1

    fnr++

endrep until fnr eq NTask

print
print, "TotNgroups   =", TotNgroups
print, "Largest      =", GroupLen(0)
print
print, "Total in grps=", TotNids
print

end
