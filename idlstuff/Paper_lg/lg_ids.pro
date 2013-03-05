pro do_the_following, junk

   ;
   ; generate the ID list, all
   ; the details need to be 
   ; set manually
   ;
   generate_idlist, "/raid4/tcox/localgroup/v4", 140



end






;=========================================
;
;  Do the grunt work of figuring out
; what the id's are of certain particles
;
;=========================================

;
; Write a file that lists the
;  ids of all gas particles which
;  become unbound from the system
;-------------------------------
;pro generate_idlist, junk
pro generate_idlist, frun, snapnum

;if not keyword_set(junk) then begin
if not keyword_set(frun) then begin
   print, "  "
   ;print, "generate_idlist, junk"
   print, "generate_idlist, frun, snapnum"
   return
endif



; -----------------
; -----------------
;frun= "/raid4/tcox/localgroup/v4"
;snapnum= 34
; -----------------
;mw_sid= 1L
;mw_npart= 500001L
mw_bhid= 500001L
; -----------------
;m31_sid= 500002L
;m31_npart= 800001L
m31_bhid= 1300002L
; -----------------
; -----------------
;frun= "/raid4/tcox/localgroup/bhires"
;mw_bhid= 817955
;m31_bhid= 2145299


idlistname= frun+"/earth_idlist.txt"

snaplbl= strcompress(string(snapnum),/remove_all)

ok= fload_snapshot_bh(frun,snapnum)

;bhid= fload_blackhole_id(1)
;bhid1= bhid[0]
;bhid2= bhid[1]
;print, "Blackhole ID's: ", bhid1, bhid2

mw_bh_xyz= fload_blackhole_xyz('xyz',idtofollow=mw_bhid,center=[0,0,0])
print, "MW blackhole ID= ", mw_bhid
print, " position = ", mw_bh_xyz

star_ids= fload_allstars_ids(1)
;r= fload_allstars_xyz('r', center=mw_bh_xyz)
rxy= fload_allstars_xyz('rxy', center=mw_bh_xyz)

idlist= [0]



; take stars in small radial range
; --------------------------------------
earthr= 8.0
dr= 0.5
;dr= 0.05

earthrlbl= strcompress(string(earthr),/remove_all)
drlbl= strcompress(string(dr),/remove_all)

;idx= where((r ge 8.0) and (r lt 8.5))
idx= where((rxy ge (earthr-dr)) and (rxy lt (earthr+dr)))
if idx(0) ne -1 then idlist= star_ids(idx)
print, 'there are ', n_elements(idx), ' star particles at the earth radius'
print, 'writing ', n_elements(idlist), ' ids'


;  Write id list to a file
; ----------------------------
cmt1= 'snap= '+snaplbl
cmt2= 'assumed earth radius= '+earthrlbl+'  +/- '+drlbl
id_filename= idlistname
print, ' writing comments'
print, '     '+cmt1
print, '     '+cmt2
print, ' to '+id_filename
write_id_list, idlist, id_filename, cmt1=cmt1, cmt2=cmt2



end




; ---------------------------------------------------------------------




;==================================
;  Write the ID list
;==================================
pro write_id_list, idlist, idfilename, cmt1=cmt1, cmt2=cmt2

if not keyword_set(cmt) then cmt=' '

; ----------------------------
;  Write id list to a file
;
openw, 1, idfilename

printf, 1, '#  star ids, file: '+idfilename
printf, 1, "#  "
printf, 1, '#  '+cmt1
printf, 1, '#  '+cmt2
printf, 1, "#  "
for i=0,n_elements(idlist)-1 do printf, 1, idlist[i]

close, 1

end






; ---------------------------------------------------------------------





;==================================
;  Load the ID list
;==================================
function load_id_list, idfilename


spawn, 'wc '+idfilename, result
lines= long(result)
idlist= lonarr(lines-5)

openr, 1, idfilename

textjunk= ''
readf, 1, textjunk
readf, 1, textjunk
readf, 1, textjunk
readf, 1, textjunk
readf, 1, textjunk
readf, 1, idlist
close, 1


return, idlist


end






; ------------------------------------------------------------------------










