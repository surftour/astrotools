;=========================================
;
;  ID's for KDC's
;
;=========================================

;
; Write a file which lists the ids
;  of star particles (new or old)
;  that satisfy some criteria
;-------------------------------
pro g_id, frun, snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "g_id, frun, snapnum"
   return
endif

;ok= fload_snapshot(frun,snapnum)
ok= fload_snapshot_bh(frun,snapnum)


center= fload_center_alreadycomp(1)
comvel= fload_comvel(1,center=center)

; all particles
;r= fload_all_xyz('r')
;m= fload_all_mass(1)

; load new star info
;nsids= fload_newstars_id(1)
;r= fload_newstars_xyz('r')
;jz= fload_newstars_j(23)

; load old star info
;osids= fload_oldstars_id(1)
;r= fload_oldstars_xyz('r')
;jz= fload_oldstars_j(23)


allstar_id= fload_allstars_ids(1)
allstar_r= fload_allstars_xyz('r')
allstar_v_tot= fload_allstars_v('tot',comvel=comvel)
allstar_v_phi= fload_allstars_v('phi',comvel=comvel)

allstar_totjx= fload_allstars_j(1,vcom=comvel)
allstar_totjy= fload_allstars_j(2,vcom=comvel)
allstar_totjz= fload_allstars_j(3,vcom=comvel)
print, "total j= ", allstar_totjx, allstar_totjy, allstar_totjz


allstar_jx= fload_allstars_j(21,vcom=comvel)
allstar_jy= fload_allstars_j(22,vcom=comvel)
allstar_jz= fload_allstars_j(23,vcom=comvel)
print, "total j= ", total(allstar_jx), total(allstar_jy), total(allstar_jz)
allstar_jtot= fload_allstars_j(20,vcom=comvel)


idx= where(allstar_r lt 5.0)
print, n_elements(allstar_r), "  total stellar particles"
print, n_elements(idx), "  particles are within 5.0 kpc/h"

if idx(0) ne -1 then begin

    as_id= allstar_id(idx)

    ; velocity
    ; ----------
    as_vratio= abs(allstar_v_phi(idx) / allstar_v_tot(idx))

    ; this shouldn't be large
    idxx= where(as_vratio gt 1.0)
    if idxx(0) ne -1 then print, "BIG PROBLEM: number of as_vratio values greater than 1= ", n_elements(idxx)

    ; grab circular orbits??
    ;idxx= where(as_vratio gt 0.95)
    ;idlist= as_id(idxx)

    ; angular momentum
    ; ------------------
    as_jx= allstar_jx(idx)
    as_jy= allstar_jy(idx)
    as_jz= allstar_jz(idx)
    as_jtot= allstar_jtot(idx)

    tjx= total(as_jx)
    tjy= total(as_jy)
    tjz= total(as_jz)
    print, "total j (<5)= ", tjx, tjy, tjz

    ; get direction of total j
    tjnorm= sqrt(tjx*tjx + tjy*tjy + tjz*tjz)
    n_x= tjx/tjnorm
    n_y= tjy/tjnorm
    n_z= tjz/tjnorm

    j_dot_n= as_jx*n_x + as_jy*n_y + as_jz*n_z
    thetadiff= acos(j_dot_n/as_jtot)*180./!PI

    ; grab aligned angular momenta??
    ;idxx= where(thetadiff lt 10.0)
    idxx= where(thetadiff lt 5.0)
    ;idxx= where(thetadiff lt 3.0)
    idlist= as_id(idxx)

endif


print, "Particles that meet criteria= ", n_elements(idlist)

; write id list
; ----------------------------------------
id_filename= frun+"/id_list_test.txt"
print, ' writing star ids to '+id_filename
write_id_list, idlist, id_filename, cmt=cmt



end






;==================================
;  Write the ID list
;==================================
pro write_id_list, idlist, idfilename, cmt=cmt

if not keyword_set(cmt) then cmt=' '

; ----------------------------
;  Write id list to a file
;
get_lun, unit
openw, unit, idfilename

printf, unit, '#  file: '+idfilename
printf, unit, "#  list of new star ids "
printf, unit, "#  "
printf, unit, '#  '+cmt
printf, unit, "#  "
for i=0L,n_elements(idlist)-1 do printf, unit, idlist[i]

close, unit

end





;==================================
;  Load the ID list
;==================================
function load_id_list, idfilename


spawn, 'wc '+idfilename, result
lines= long(result)
idlist= lonarr(lines-5)

get_lun, unit
openr, unit, idfilename

textjunk= ''
readf, unit, textjunk
readf, unit, textjunk
readf, unit, textjunk
readf, unit, textjunk
readf, unit, textjunk
readf, unit, idlist
close, unit


return, idlist


end





;==============================================
;
;
;==============================================


