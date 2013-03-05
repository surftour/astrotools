pro do_the_following, junk

   ;
   ; generate the ID list, all
   ; the details need to be 
   ; set manually
   ;
   generate_idlist, "data/altsf/vc3vc3e_N2", 0



end






;=========================================
;
;  Do the grunt work of figuring out
; what the id's are of certain particles
;
;=========================================





;
; Write a file that lists the
;  ids of all stars at 8+/- some radius
;-------------------------------
pro generate_idlist, frun, snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "generate_idlist, frun, snapnum"
   return
endif


; -----------------
; -----------------


idlistname= frun+"/earth_idlist.txt"

snaplbl= strcompress(string(snapnum),/remove_all)

ok= fload_snapshot_bh(frun,snapnum)

bhid= fload_blackhole_id(1)
print, bhid
;bhid= bhid[0]
;bhid= bhid[1]
;bhid= 200001L
;bhid= 280002L   ; used for z3/b4e
;bhid= 400002L   ; used for ds/vc3vc3e_2
bhid= long(max(bhid))   ; uses minimum bhid - then you don't need to know the BH ID
center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)

print, "Blackhole ID: ", bhid
print, "Blackhole center: ", center_bh

star_ids= fload_allstars_ids(1)
rxy= fload_allstars_xyz('rxy', center=center_bh)

idlist= [0]



; take stars in small radial range
; --------------------------------------
earthr= 8.0
;dr= 0.5
dr= 0.1
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




;
; Write a file that lists the
;  ids of all stars at 8+/- some radius
;-------------------------------
pro generate_idlist_2, frun, snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "generate_idlist_2, frun, snapnum"
   return
endif


; -----------------
; -----------------


idlistname= frun+"/earth_idlist.txt"

snaplbl= strcompress(string(snapnum),/remove_all)

ok= fload_snapshot_bh(frun,snapnum)

bhid= fload_blackhole_id(1)
print, bhid
;bhid= bhid[0]
;bhid= bhid[1]
;bhid= 200001L
;bhid= 280002L   ; used for z3/b4e
;bhid= 400002L   ; used for ds/vc3vc3e_2
bhid= long(max(bhid))   ; uses minimum bhid - then you don't need to know the BH ID
bhidlbl= strcompress(string(bhid), /remove_all)
center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
centerlbl= string(center_bh[0])+','+string(center_bh[1])+','+string(center_bh[2])

print, "Blackhole ID: ", bhid
print, "Blackhole center: ", center_bh

star_ids= fload_allstars_ids(1)
x= fload_allstars_xyz('x', center=center_bh)
y= fload_allstars_xyz('y', center=center_bh)
z= fload_allstars_xyz('z', center=center_bh)
rxy= fload_allstars_xyz('rxy', center=center_bh)




; take stars in small radial range
; --------------------------------------
earthr= 8.0
dr= 0.3
rpts= 360

earthrlbl= strcompress(string(earthr),/remove_all)
rlbl= strcompress(string(rpts),/remove_all)

;idx= where((r ge 8.0) and (r lt 8.5))
idx= where((rxy ge (earthr-dr)) and (rxy lt (earthr+dr)))
if idx(0) eq -1 then stop

ex= x(idx)
ey= y(idx)
ez= z(idx)
eids= star_ids(idx)

idlist= lonarr(rpts)

for i=0, rpts-1 do begin
	thisradian = i * (2*!PI/rpts)
	thisx= earthr * cos(thisradian)
	thisy= earthr * sin(thisradian)
	thisr= sqrt((ex-thisx)*(ex-thisx) + (ey-thisy)*(ey-thisy))
	eidx= where(thisr le min(thisr))
	idlist[i]= eids(eidx)
endfor

print, 'there are ', n_elements(idx), ' star particles at the earth radius'
print, 'writing ', n_elements(idlist), ' ids'




;  Write id list to a file
; ----------------------------
cmt1= 'snap= '+snaplbl+'  bhid= '+bhidlbl+'    center= '+centerlbl
cmt2= 'assumed earth radius= '+earthrlbl+'  numberofparticles= '+rlbl
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










