
Num = 0 

Base= "/afs/mpa/home/snhess/P-Gadget3/SBcluster"

GrNr = 0

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
            lovSfr = fltarr(Ngroups)
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
print, "GroupLen("+string(GrNr)+") =", GroupLen(GrNr)
print


Id_List = lonarr(GroupLen(GrNr))


goff = GroupOffset(GrNr)
glen =  GroupLen(GrNr)
skip = 0L
fnr = 0L

repeat begin

    f = Base + "/groups_" + exts +"/group_ids_"+exts +"."+strcompress(string(fnr),/remove_all)

    Ngroups = 0L
    TotNgroups = 0L
    Nids = 0L
    TotNids = 0LL
    NTask = 0L
    Offset = 0L


    openr,1,f
    readu,1, Ngroups, TotNgroups, Nids, TotNids, NTask, Offset


    if (goff lt (Offset + Nids)) and (goff ge Offset) then begin
        point_lun, 1, goff - offset
        
;;    print,"forward in file", fnr," by ", goff-offset

        if (Offset + Nids) lt (goff + glen) then begin
            len = Offset + Nids - goff
        endif else begin
            len = glen
        endelse

        if len gt 0 then begin
            locIDs = lonarr(len)
            readu,1,locIDs 

;;    print, "read ", len, " particles in file", fnr

            Id_List(skip:skip+len-1) = locIDs(*)
            skip += len
            glen -= len
            goff += len 
        endif
    endif
    close,1

    fnr++

endrep until fnr eq NTask


end
