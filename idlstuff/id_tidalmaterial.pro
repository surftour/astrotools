pro do_the_following, junk

   ;
   ; generate the ID list, all
   ; the details need to be 
   ; set manually
   ;


   ;generate_idlist, "data1/tests/Sb10xSb10xhs_e_256", 200


   ;generate_idlist, "data/altsf/vc3vc3e_N2", 50

   generate_idlist, "data1/tides/4x_pro", 61

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

;bhid= bhid[0]
;bhid= bhid[1]
;bhid= 200001L
;bhid= 280002L   ; used for z3/b4e
;bhid= 400002L   ; used for ds/vc3vc3e_2
;bhid= long(max(bhid))   ; uses minimum bhid - then you don't need to know the BH ID
bhid= long(min(bhid))   ; uses minimum bhid - then you don't need to know the BH ID
center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)

bhidlbl= strcompress(string(bhid), /remove_all)
center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
centerlbl= string(center_bh[0])+','+string(center_bh[1])+','+string(center_bh[2])


print, "Blackhole ID: ", bhid
print, "Blackhole center: ", center_bh

star_ids= fload_allstars_ids(1)
x= fload_allstars_xyz('x', center=center_bh)
y= fload_allstars_xyz('y', center=center_bh)
z= fload_allstars_xyz('z', center=center_bh)

r= sqrt(x*x + y*y + z*z)


;
; grab just one galaxy's ID's
;

idx= where((star_ids gt 1) and (star_ids le bhid_min))
;idx= where((star_ids gt bhid_min) and (star_ids le bhid_max))
if idx(0) eq -1 then stop
star_ids= star_ids(idx)
x= x(idx)
y= y(idx)
z= z(idx)
r= r(idx)



;-----------------------------------------
;
; take stars outside of some radius
;
;fiducial_disk_r= 8.0 & fiducial_disk_r_lbl= strcompress(string(fiducial_disk_r),/remove_all)
fiducial_disk_r= 15.0 & fiducial_disk_r_lbl= strcompress(string(fiducial_disk_r),/remove_all)

idx= where(r ge fiducial_disk_r)
if idx(0) eq -1 then stop

ttr= r(idx)
ttx= x(idx)
tty= y(idx)
ttz= z(idx)
ttids= star_ids(idx)





;-----------------------------------------
;
; do cut along an axis (if desired)
;

do_other_cut= 1
;do_other_cut= 0
if do_other_cut gt 0 then begin

	;
	; outside arm
	;
	;idx_oc= where((ttx lt -10.0) or (tty lt -5.0))
	;idx_oc= where((ttx lt -8.0) or (tty lt -8.0))
	;idx_oc= where(ttx lt -5.0) & print, "cut at x<-5"
	;idx_oc= where(tty lt -8.0) & print, "cut at y<-8"
	;idx_oc= where(ttx gt 0.0)
	idx_oc= where(ttx lt 0.0)

	;
	; inside arm
	;
	;idx_oc= where((ttx gt -5.0) and (tty gt 10.0))

	if idx_oc(0) ne -1 then begin
		ttr= ttr(idx_oc)
		ttx= ttx(idx_oc)
		tty= tty(idx_oc)
		ttz= ttz(idx_oc)
		ttids= ttids(idx_oc)
	endif
endif






;---------------------------------------

idlist= long(ttids)


print, 'there are ', n_elements(idx), ' star particles at the tail'
print, 'writing ', n_elements(idlist), ' ids'




;  Write id list to a file
; ----------------------------
idlistname= frun+"/tidal_tail_idlist.txt"

snaplbl= strcompress(string(snapnum),/remove_all)
cmt1= 'snap= '+snaplbl+'  bhid= '+bhidlbl+'    center= '+centerlbl
cmt2= 'fiducial_disk_r= '+fiducial_disk_r_lbl+'  numberofparticles= '+strcompress(string(n_elements(idlist)),/remove_all)
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










