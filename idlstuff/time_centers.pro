pro time_centers, frun, $
		init_bhids=init_bhids, $
		init_gal1= init_gal1, $
		init_gal2= init_gal2



if not keyword_set(frun) then begin
	print, " "
	print, " time_centers, frun, init_bhids=init_bhids, init_gal1=init_gal1, init_gal2= init_gal2"
	print, " "
	print, " "
	return
endif


;
; use BH information to determine
; the center
;
ok= fload_snapshot_bh(frun,0,/nopot_in_snap,/skip_center)
nsnaps= fload_frun_nsnaps(frun)
init_bhids= fload_blackhole_id(1,n_bh=n_bh)
print, "Blackhole ID's: ", init_bhids


; ----------------------------------------
; ----------------------------------------
;
; write headers to the files
;


for bhi=0, n_bh-1 do begin

        bhlbl= strcompress(string(init_bhids[bhi]),/remove_all)

        ; center file
        ; -------------------------------
        ;write_line_to_file, frun+'/center_bh_'+bhlbl+'.txt', comment="#  center for bhid= "+bhlbl
        ;write_line_to_file, frun+'/center_bh_'+bhlbl+'.txt', comment="#  "
        ;write_line_to_file, frun+'/center_bh_'+bhlbl+'.txt', comment="#  "
        ;write_line_to_file, frun+'/center_bh_'+bhlbl+'.txt', comment="# snap   time    --------------- dispersions (in km/sec) ------------------ "
        ;write_line_to_file, frun+'/center_bh_'+bhlbl+'.txt', comment="#  num  (Gyr)          xy         xz         yz        avg        err       "

endfor




; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------


snapi= 0
;snapi= 1202

repeat begin

	print, "--------------------------------------"


	;  open snapshot
	; ---------------
	print, "snapi= ", snapi
	if snapi ne 0 then begin
	  if snapi ge 1000 then begin
		 ok=fload_snapshot_bh(frun,snapi,/nopot_in_snap,/do_four,/skip_center)
	  endif else begin
		 ok=fload_snapshot_bh(frun,snapi,/nopot_in_snap,/skip_center)
	  endelse
	endif



        ;  now grab info
        ; ---------------
        if ok eq -1 then begin
                ; just skip this snapshot, but goto the next one
        endif else begin

	    ; what time is it?
	    time= fload_time(1)


	    ;
	    ; ---------------------------------------------------
	    for bhi= 0, n_bh-1 do begin

		if bhi eq 0 then sid= long(1) else sid= init_bhids[bhi-1]+1
		if bhi eq 0 then npart= init_bhids[bhi] else npart= init_bhids[bhi]-init_bhids[bhi-1]

		if (bhi eq 0) and keyword_set(init_gal1) and n_elements(init_gal1) eq 2 then begin
			sid= long(init_gal1[0])
			npart= long(init_gal1[1])
		endif

		if (bhi eq 1) and keyword_set(init_gal2) and n_elements(init_gal2) eq 2 then begin
			sid= long(init_gal2[0])
			npart= long(init_gal2[1])
		endif

		print, "========================================"
		print, "bhi= ", bhi
		print, "sid= ", sid
		print, "npart= ", npart

		bhlbl= strcompress(string(init_bhids[bhi]),/remove_all)

		; straight-up 
		com= fload_1gal_com(sid,npart)
		write_line_to_file, frun+"/com_"+bhlbl+".txt", [time,com]

		; density method
		center_den= fload_1gal_center(sid,npart)
		;center_den= fload_1gal_baryonic_center(sid,npart)
		;center_den= fload_1gal_dm_center(sid,npart)
		write_line_to_file, frun+"/centers_den_"+bhlbl+".txt", [time,center_den]

		; center-of-mass velocity, and central density velocity
		;   both of these need the particle numbers to be set
		;   correctly.
		comvel_den= fload_1gal_comvel(sid,npart,center=center_den,rfact=1.0)
		write_line_to_file, frun+"/comvel_den_"+bhlbl+".txt", [time,comvel_den]
		cenvel_den= fload_comvel(1,center=center_den,rfact=0.001)
		write_line_to_file, frun+"/cenvel_den_"+bhlbl+".txt", [time,cenvel_den]


		; now, use BH positions
		bhidx= where(init_bhids[bhi] eq fload_blackhole_id(1))
		;if fload_npart(5) ge 1 then begin
		if bhidx(0) ne -1 then begin
			;print, "in bh center routine"
			;print, "fload_npart(5)= ", fload_npart(5)
			;print, "idtofollow= ", init_bhids[bhi]
			;print, "fload_blackhole_ids(1)= ", fload_blackhole_id(1)
			center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=init_bhids[bhi])
			write_line_to_file, frun+"/centers_bh_"+bhlbl+".txt", [time, center_bh]

			; center-of-mass velocity, and central density velocity
			;   both of these need the particle numbers to be set
			;   correctly.
			comvel_bh= fload_1gal_comvel(sid,npart,center=center_bh,rfact=1.0)
			write_line_to_file, frun+"/comvel_bh_"+bhlbl+".txt", [time,comvel_bh]
			cenvel_bh= fload_comvel(1,center=center_bh,rfact=0.001)
			write_line_to_file, frun+"/cenvel_bh_"+bhlbl+".txt", [time,cenvel_bh]
		endif else begin
			write_line_to_file, frun+"/centers_bh_"+bhlbl+".txt", [time, -1, -1, -1] 
			write_line_to_file, frun+"/comvel_bh_"+bhlbl+".txt", [time, -1, -1, -1]
			write_line_to_file, frun+"/cenvel_bh_"+bhlbl+".txt", [time, -1, -1, -1]
		endelse
	    endfor
	endelse

	snapi= snapi + 1

        if snapi gt (nsnaps-1) then goto, done

endrep until (snapi gt 5999)


done:

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



; ----------------------------------------


;time_relative_diffs, frun


; ----------------------------------------



end




;===============================================================================




;===============================================================================



pro find_rdiffs, frun, time=time

        ; read centers
        ; -------------
        read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=frun+"/centers.txt"

        cp_time= transpose(cp_time)
        center_1= transpose(cp_cen1)
        center_2= transpose(cp_cen2)

	x = center_1[*,0] - center_2[*,0]
	y = center_1[*,1] - center_2[*,1]
	z = center_1[*,2] - center_2[*,2]

	rdiff = sqrt(x*x + y*y + z*z)

	idx=where(rdiff lt 15.0)

	if idx(0) ne -1 then begin
	   for i=0, (n_elements(idx)-1) do begin
		print, "snap= ", i, cp_time(idx(i)), rdiff(idx(i))
		if rdiff(idx(i)) le 0.0 then begin
		  print, "bh's have merged"
		  break
		endif
	   endfor
	endif


end






pro find_rdiff, fbase, bhlbl1, bhlbl2, rdiff=rdiff, time=time

	;
	; gal 1
	;
	read_center, time1, center1, filename=fbase+bhlbl1+".txt"

	;
	; gal 2
	;
	read_center, time2, center2, filename=fbase+bhlbl2+".txt"


	time= transpose(time1)
	center1= transpose(center1)
	center2= transpose(center2)

	x = center1[*,0] - center2[*,0]
	y = center1[*,1] - center2[*,1]
	z = center1[*,2] - center2[*,2]
	rdiff= sqrt(x*x + y*y + z*z)

	; check for disappearing BH's
	idx= where((center1(*,0) eq -1.0000) and (center1(*,1) eq -1.0000) and (center1(*,2) eq -1.0000))
	if idx(0) ne -1 then rdiff(idx)= 0.00
	idx= where((center2(*,0) eq -1.0000) and (center2(*,1) eq -1.0000) and (center2(*,2) eq -1.0000))
	if idx(0) ne -1 then rdiff(idx)= 0.00

end





;===============================================================================
;
;
;

pro time_relative_diffs, frun

if not keyword_set(frun) then begin
        print, " "
        print, " time_relative_diffs, frun"
        print, " "
        print, "  * requires center files"
        print, " "
        return
endif


;
; get bh ids
;
ok= fload_snapshot_bh(frun,0,/nopot_in_snap)
init_bhids= fload_blackhole_id(1)
init_bhids= long(init_bhids(sort(init_bhids)))
n_bh= n_elements(init_bhids)
print, "Blackhole ID's: ", init_bhids


;
; gal 1
;
bhlbl1= strcompress(string(init_bhids[0]),/remove_all)

;
; gal 2
;
if n_bh gt 1 then begin
   bhlbl2= strcompress(string(init_bhids[1]),/remove_all)
endif else begin
	print, " "
	print, "              *****  PROBLEM  ****** "
	print, " "
	print, "   there is only one BH in this run - stopping "
	print, " "
	print, " "
	return
endelse



find_rdiff, frun+"/centers_den_", bhlbl1, bhlbl2, rdiff=rdiff_den, time=time
find_rdiff, frun+"/centers_bh_", bhlbl1, bhlbl2, rdiff=rdiff_bh

find_rdiff, frun+"/cenvel_den_", bhlbl1, bhlbl2, rdiff=veldiff_cenden
find_rdiff, frun+"/cenvel_bh_", bhlbl1, bhlbl2, rdiff=veldiff_cenbh

find_rdiff, frun+"/comvel_den_", bhlbl1, bhlbl2, rdiff=veldiff_comden
find_rdiff, frun+"/comvel_bh_", bhlbl1, bhlbl2, rdiff=veldiff_combh



;
; write to file(s)
;
filename= frun+"/centers_den_reldiff.txt"
openw, 1, filename
printf, 1, frun
printf, 1, n_elements(time)
data= fltarr(4,n_elements(time))
data[0,*]=time
data[1,*]=rdiff_den
data[2,*]=veldiff_cenden
data[3,*]=veldiff_comden
printf, 1, data
close, 1

filename= frun+"/centers_bh_reldiff.txt"
openw, 1, filename
printf, 1, frun
printf, 1, n_elements(time)
data= fltarr(4,n_elements(time))
data[0,*]=time
data[1,*]=rdiff_bh
data[2,*]=veldiff_cenbh
data[3,*]=veldiff_combh
printf, 1, data
close, 1



end






; ----------------------------------------











