pro toomre_centers, frun


if not keyword_set(frun) then begin
	print, " "
	print, " toomre_centers, frun"
	print, " "
	print, " "
	return
endif



; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------


snapi= 0
done= 0

repeat begin

	print, "--------------------------------------"


	;  open snapshot
	; ---------------
	print, "snapi= ", snapi
	if snapi ge 1000 then begin
		 ok=fload_snapshot_bh(frun,snapi,/nopot_in_snap,/do_four,/skip_center)
	endif else begin
		 ok=fload_snapshot_bh(frun,snapi,/nopot_in_snap,/skip_center)
	endelse


	;  now grab info
	; ---------------
	if ok eq -1 then begin
		done= 1
	endif else begin

		; what time is it?
		time= fload_time(1)


		x0= fload_halo_xyz('x',center=[0,0,0])
		y0= fload_halo_xyz('y',center=[0,0,0])
		z0= fload_halo_xyz('z',center=[0,0,0])
		vx0= fload_halo_v('x',center=[0,0,0])
		vy0= fload_halo_v('y',center=[0,0,0])
		vz0= fload_halo_v('z',center=[0,0,0])
		m0= fload_halo_mass(1)
		ids0= fload_halo_id(1)
		n_h= fload_npart(1)


		;
		; ---------------------------------------------------
		for hi= 0, n_h-1 do begin

			hlbl= strcompress(string(ids0[hi]),/remove_all)

			;print, "========================================"
			print, "hi= ", hi


			;  position
			; ----------
			hlbl= strcompress(string(ids0[hi]),/remove_all)
			filename= frun+"/center_"+hlbl+".txt"
			if not file_test(filename) then begin
				openw, 1, filename
				printf, 1, filename
				printf, 1, -1
			endif else begin
				openu, 1, filename, /append
			endelse
			printf, 1, FORMAT= '(F6.3,"  ",3(F10.3,"  "))', time, x0[hi], y0[hi], z0[hi]
			close, 1


			;  velocity
			; ----------
                        filename= frun+"/cenvel_"+hlbl+".txt"
                        if not file_test(filename) then begin
                                openw, 1, filename
                                printf, 1, filename
                                printf, 1, -1
                        endif else begin
                                openu, 1, filename, /append
                        endelse
                        printf, 1, FORMAT= '(F6.3,"  ",3(F10.3,"  "))', time, vx0[hi], vy0[hi], vz0[hi]
                        close, 1
		endfor
	endelse

	snapi= snapi + 1

endrep until (done eq 1)
        



print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end






;===============================================================================







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


end




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
find_rdiff, frun+"/center_", '1', '233', rdiff=rdiff, time=time
find_rdiff, frun+"/cenvel_", '1', '233', rdiff=veldiff, time=time


;
; write to file(s)
;
filename= frun+"/centers_reldiff.txt"
openw, 1, filename
printf, 1, frun
printf, 1, n_elements(time)
data= fltarr(4,n_elements(time))
data[0,*]=time
data[1,*]=rdiff
data[2,*]=veldiff
data[3,*]=veldiff
printf, 1, data
close, 1


end



