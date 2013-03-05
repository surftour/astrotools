pro do_the_following, junk

   ;
   ; generate the ID list, all
   ; the details need to be 
   ; set manually
   ;

   generate_idlist, "data/z3/b5e", 75

end






;=========================================
;
;  Do the grunt work of figuring out
; what the id's are of certain particles
;
;=========================================





;
; Write a file that lists the
;  ids of all stars in the tidal tail
;  (where tidal tail has a "flexible"
;   definition)
;-------------------------------
pro generate_idlist, frun, snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "generate_idlist, frun, snapnum"
   return
endif


; -----------------
; -----------------


ok= fload_snapshot_bh(frun,snapnum,/nopot_in_snap)


bhid= fload_blackhole_id(1)
print, bhid
bhid_max= long(max(bhid))   ;   assume there is only 2 galaxies
bhid_min= long(min(bhid))   ; 
bhid= bhid_min              ; uses minimum bhid - then you don't need to know the BH ID
bhidlbl= strcompress(string(bhid), /remove_all)
center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
centerlbl= string(center_bh[0])+','+string(center_bh[1])+','+string(center_bh[2])


print, "Blackhole ID: ", bhid
print, "Blackhole center: ", center_bh

star_ids= fload_allstars_ids(1)
x= fload_allstars_xyz('x', center=center_bh)
y= fload_allstars_xyz('y', center=center_bh)
z= fload_allstars_xyz('z', center=center_bh)

;r= sqrt(x*x + y*y + z*z)
r= sqrt(x*x + y*y)


;
; grab large r particles
;

idx= where((r gt 4.0))
if idx(0) eq -1 then stop
star_ids= star_ids(idx)
x= x(idx)
y= y(idx)
z= z(idx)
r= r(idx)


;---------------------------------------

idlist= long(star_ids)


print, 'there are ', n_elements(idx), ' star particles at the tail'
print, 'writing ', n_elements(idlist), ' ids'




;  Write id list to a file
; ----------------------------
idlistname= frun+"/idlist_rem_rgt4.txt"

snaplbl= strcompress(string(snapnum),/remove_all)
cmt1= 'snap= '+snaplbl+'  bhid= '+bhidlbl+'    center= '+centerlbl
cmt2= 'all particle with r greater than 4: '+strcompress(string(n_elements(idlist)),/remove_all)
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
for i=0L,n_elements(idlist)-1 do printf, 1, idlist[i]

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










