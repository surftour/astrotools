
BaseDir = "/afs/rzg/bc-b/vrs/tmp/"   ; The output-directory of the simulation
SnapBase= "snap_S1"                   ; The basename of the snapshot files

Num = 0                            ; The number of the snapshot we look at

;;; The object-file of the compile C-library for accessing the group catalogue

ObjectFile = "/afs/rzg/bc-b/vrs/tmp/P-Gadget2/IdlGroupLib/idlgrouplib.so"

Num=long(Num)
exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)
Outputdir  = Basedir + "/"



;;;;;;;;;; First, we get the number of groups

Ngroups = CALL_EXTERNAL(ObjectFile, $
                       'get_total_number_of_groups', /UL_VALUE, $
                        OutputDir, $
                        Num)

print, "Number of groups in the catalogue: ", Ngroups

;;;;;;;;;; Now we load the group catalogue


GroupLen = lonarr(Ngroups)
GroupFileNr = lonarr(Ngroups)
GroupNr = lonarr(Ngroups)
GroupTypeLen = lonarr(6, Ngroups)
GroupTypeMass = dblarr(6, Ngroups)
GroupCenter = fltarr(3, Ngroups)
GroupSfr = fltarr(Ngroups)

Ngroups = CALL_EXTERNAL(ObjectFile, $
                       'get_group_catalogue', /UL_VALUE, $
                        OutputDir, $
                        Num, $
                        GroupLen, $
                        GroupFileNr, $
                        GroupNr, $
                        GroupTypeLen, $
                        GroupTypeMass, $
                        GroupCenter, $
			GroupSfr)

;;; Now let's plot a simple cumulative mass function


Count= lindgen(Ngroups) + 1

plot, GroupLen, Count, /xlog, /ylog


Gas_Count= GroupTypeLen(0,*)
DM_Count= GroupTypeLen(1,*)
Star_Count= GroupTypeLen(4,*)
BH_Count= GroupTypeLen(5,*)



end
