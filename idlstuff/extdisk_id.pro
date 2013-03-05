pro do_the_following, junk

   ;
   ; generate the ID list, all
   ; the details need to be 
   ; set manually
   ;

   generate_idlist, "data2/Sbc_10x_wBH", 0,  0.0,  10.0, filename="idlist_r_10.txt"
   generate_idlist, "data2/Sbc_10x_wBH", 0, 10.0,  20.0, filename="idlist_r_20.txt"
   generate_idlist, "data2/Sbc_10x_wBH", 0, 20.0,  30.0, filename="idlist_r_30.txt"
   generate_idlist, "data2/Sbc_10x_wBH", 0, 30.0,  50.0, filename="idlist_r_50.txt"
   generate_idlist, "data2/Sbc_10x_wBH", 0, 50.0, 100.0, filename="idlist_r_100.txt"

end






;=========================================
;
;  Do the grunt work of figuring out
;  what are particle id's at certain
;  radii, and write a file
;
;=========================================



pro generate_idlist, frun, snapnum, rmin, rmax, filename=filename

if not keyword_set(filename) then filename= 'idlist.txt'

if not keyword_set(frun) then begin
   print, "  "
   print, "generate_idlist, frun, snapnum, rmin, rmax, filename=filename"
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

ids= fload_gas_ids(1)
x= fload_gas_xyz('x', center=center_bh)
y= fload_gas_xyz('y', center=center_bh)
z= fload_gas_xyz('z', center=center_bh)

r= sqrt(x*x + y*y + z*z)


;---------------------------------------


;
; grab particles at specific radii
;

idx= where((r gt rmin) and (r lt rmax))
if idx(0) eq -1 then return
star_ids= star_ids(idx)
x= x(idx)
y= y(idx)
z= z(idx)
r= r(idx)


;---------------------------------------

idlist= long(ids)

print, 'there are '+string(n_elements(idx))+' particles between '+string(rmin)+' and '+string(rmax)
print, 'writing ', n_elements(idlist), ' ids'



;  Write id list to a file
; ----------------------------

snaplbl= strcompress(string(snapnum),/remove_all)
cmt1= 'snap= '+snaplbl+'  bhid= '+bhidlbl+'    center= '+centerlbl
cmt2= 'particle between: rmin= '+strcompress(string(rmin),/remove_all)+', rmax= '+strcompress(string(rmax),/remove_all)

print, ' writing comments'
print, '     '+cmt1
print, '     '+cmt2
print, ' to file: '+frun+filename
write_id_list, idlist, frun+filename, cmt1=cmt1, cmt2=cmt2



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










