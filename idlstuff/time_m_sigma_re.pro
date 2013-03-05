; ---------------------------------------------------------------------------
;
;   Find the sizes and dispersions as a function of time.
;
;
; ---------------------------------------------------------------------------


;=========================================================
;
;  Main program here.
;
;=========================================================

pro time_m_sigma_re, frun, first_n_last_only=first_n_last_only


if not keyword_set(frun) then begin
	print, " "
	print, " m_sigma_time, frun"
	print, " "
	print, " "
	return
endif



nsnaps= fload_frun_nsnaps(frun)
ok= fload_snapshot_bh(frun,0,/nopot_in_snap,/skip_center)
init_bhids= fload_blackhole_id(1,n_bh=n_bh)
n_bh_init= n_bh
print, "Blackhole ID's: ", init_bhids


if keyword_set(first_n_last_only) then begin
	lastsnapnum= nsnaps
	nsnaps= 2
endif



; ----------------------------------------
; ----------------------------------------
;
; write headers to the files
;


for bhi=0, n_bh-1 do begin

	bhlbl= strcompress(string(init_bhids[bhi]),/remove_all)

	; sigma file
	; -------------------------------
	write_line_to_file, frun+'/sigma_'+bhlbl+'.txt', comment="#  sigma for bhid= "+bhlbl
	write_line_to_file, frun+'/sigma_'+bhlbl+'.txt', comment="#  "
	write_line_to_file, frun+'/sigma_'+bhlbl+'.txt', comment="#  "
	write_line_to_file, frun+'/sigma_'+bhlbl+'.txt', comment="# snap   time    --------------- dispersions (in km/sec) ------------------ "
	write_line_to_file, frun+'/sigma_'+bhlbl+'.txt', comment="#  num  (Gyr)          xy         xz         yz        avg        err       "

	; R_e file
	; -------------------------------
	write_line_to_file, frun+'/R_e_'+bhlbl+'.txt', comment="#  half-mass radius for bhid= "+bhlbl
	write_line_to_file, frun+'/R_e_'+bhlbl+'.txt', comment="#  "
	write_line_to_file, frun+'/R_e_'+bhlbl+'.txt', comment="#  "
	write_line_to_file, frun+'/R_e_'+bhlbl+'.txt', comment="# snap   time    ---- R_e (in kpc/h) ----"
	write_line_to_file, frun+'/R_e_'+bhlbl+'.txt', comment="#  num  (Gyr)          avg        err    "

endfor




; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------
snapi= 0

for i=0,nsnaps-1 do begin

	print, "--------------------------------------"

	if keyword_set(first_n_last_only) and i ne 0 then snapi= lastsnapnum-1

	; open snapshot
	repeat begin
	    ;ok=fload_snapshot_bh(frun,snapi)
	    ok=fload_snapshot_bh(frun,snapi,/nopot_in_snap)
	    ;ok=fload_snapshot(frun,snapi)
	    snapi= snapi+1
	    if snapi gt nsnaps then goto, done
	endrep until (ok eq 0)

	; what time is it?
	time= fload_time(1)

	; grab this snaps BH info
	bhids= fload_blackhole_id(1,n_bh=n_bh)

	c_sep= 0
	if n_bh gt 1 then begin
		c1= fload_center_alreadycomp(time,denorbh='den',bhlbl=strcompress(string(bhids[0]),/remove_all))
		c2= fload_center_alreadycomp(time,denorbh='den',bhlbl=strcompress(string(bhids[1]),/remove_all))
		c_sep= sqrt(total((c2-c1)*(c2-c1)))
	endif

	;
	; loop through this snaps BHs
	;
	for bhi= 0, n_bh-1 do begin

		;
		; get center of this galaxy - from the file (if possible)
		;
		;center= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=init_bhid[bhi])
		;center= fload_center_alreadycomp(time,denorbh='bh',bhlbl=strcompress(string(bhids[bhi]),/remove_all))
		center= fload_center_alreadycomp(time,denorbh='den',bhlbl=strcompress(string(bhids[bhi]),/remove_all))     ; there's no -1's here

		
		; all luminous material
		; -----------------------
		sx= fload_allstars_xyz('x', center=[0,0,0])
		sy= fload_allstars_xyz('y', center=[0,0,0])
		sz= fload_allstars_xyz('z', center=[0,0,0])
             
		svx= fload_allstars_v('x')
		svy= fload_allstars_v('y') 
		svz= fload_allstars_v('z')

		smass= fload_allstars_mass(1)


		; calculate sigma and half-mass radius
		; -------------------------------------
		if c_sep le 2.0 then begin
			Reff= process_2d_halfmass(sx, sy, smass, center=center)
			sigma_avg= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 1.0*Reff, $
					avg_Reff=avg_Reff, center= center, sigmaerr=sigmaerr, refferr=refferr)
		endif else begin
			if bhi eq 0 then startid= long(1) else startid= bhids[bhi-1]+1
			if bhi eq 0 then numpart= bhids[bhi] else numpart= bhids[bhi]-bhids[bhi-1]+1
			thisgalmass= total(fload_1gal_allstars_mass(startid,numpart))
			Reff= process_2d_halfmass(sx, sy, smass, center=center, fixed_mass=thisgalmass*0.5)
			sigma_avg= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 1.0*Reff, $
					avg_Reff=avg_Reff, center= center, sigmaerr=sigmaerr, refferr=refferr, fixed_mass=thisgalmass*0.5)
		endelse
		sigma_xy= process_2d_sigma(sx, sy, svz, 1.0*Reff, center= center)
		sigma_xz= process_2d_sigma(sx, sz, svy, 1.0*Reff, center= [center[0],center[2]])
		sigma_yz= process_2d_sigma(sy, sz, svx, 1.0*Reff, center= [center[1],center[2]])
		
		sigma_err= sigmaerr
		re_avg= avg_Reff
		re_err= refferr
		print, "effective radius= ",Reff," (xy) ", avg_Reff, " (avg)"


		; write to files
		; ----------------
		write_line_to_file, frun+'/sigma_'+strcompress(string(bhids[bhi]),/remove_all)+'.txt', [time, sigma_xy, sigma_xz, sigma_yz, sigma_avg, sigma_err], col1txt="  "+fload_finfo_exts(1)
		write_line_to_file, frun+'/R_e_'+strcompress(string(bhids[bhi]),/remove_all)+'.txt', [time, re_avg, re_err], col1txt="  "+fload_finfo_exts(1)
	endfor

endfor
        
; ----------------------------------------
; ----------------------------------------

done:

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end






;=========================================================
;
;  Alternate version for sigma investigation
;
;=========================================================

pro time_many_sigmas, frun



if not keyword_set(frun) then begin
	print, " "
	print, " time_many_sigmas, frun"
	print, " "
	print, " "
	return
endif



; this assumes they are ordered
;  0 through


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])



time= fltarr(nsnaps)

bmass= fltarr(nsnaps)


; all stars
; ----------
st_sigma_1= fltarr(nsnaps)
st_sigma_2= fltarr(nsnaps)
st_sigma_3= fltarr(nsnaps)
st_sigma_4= fltarr(nsnaps)
st_sigma_5= fltarr(nsnaps)
st_sigma_6= fltarr(nsnaps)



; gas
; -----
g_sigma_1= fltarr(nsnaps)
g_sigma_2= fltarr(nsnaps)
g_sigma_3= fltarr(nsnaps)
g_sigma_4= fltarr(nsnaps)
g_sigma_5= fltarr(nsnaps)
g_sigma_6= fltarr(nsnaps)



; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------
snapi= 0

for i=0,nsnaps-1 do begin

	print, "--------------------------------------"

	; open snapshot
	;ok=fload_snapshot(frun,snapnum)
	repeat begin
	    ok=fload_snapshot_bh(frun,snapi)
	    snapi= snapi+1
	endrep until (ok eq 0)

	; what time is it?
	time[i]= fload_time(1)


	; get black hole mass
	bmass[i]= total(fload_blackhole_mass(1))


	; get blackhole id's
	if i eq 0 then begin
	     bhid= fload_blackhole_id(1)
	     bhid1= bhid[0]
	     bhid2= bhid[1]
	     print, "Blackhole ID's: ", bhid1, bhid2
	endif


	; get centers, for the moment we do this by hand, but could
	; generalize this in the future to search centers.txt file
	twobhs= 0
	if fload_npart(5) gt 1 then twobhs= 1

	if twobhs eq 1 then begin
             ;center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=170001)
             ;center_2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=340002)
             center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
             center_2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
	     rdiff= (center_1[0]-center_2[0])*(center_1[0]-center_2[0])  +  $
	            (center_1[1]-center_2[1])*(center_1[1]-center_2[1])  +  $
	            (center_1[2]-center_2[2])*(center_1[2]-center_2[2])
	     rdiff= sqrt(rdiff)
	     if rdiff lt 5.0 then twobhs= 0
	     print, "center 1: ", center_1
	     print, "center 2: ", center_2, "  (rdiff= ",rdiff,")"
        endif else begin
             bhid= fload_blackhole_id(1)
             center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
             center_2= center_1
	     print, "center: ", center_1
        endelse


	; all luminous material
	; -----------------------
	sx= fload_allstars_xyz('x', center=[0,0,0])
	sy= fload_allstars_xyz('y', center=[0,0,0])
	sz= fload_allstars_xyz('z', center=[0,0,0])
             
	svx= fload_allstars_v('x')
	svy= fload_allstars_v('y') 
	svz= fload_allstars_v('z')

	smass= fload_allstars_mass(1)

        ; gaseous material
        ; -----------------------
        gx= fload_gas_xyz('x', center=[0,0,0])
        gy= fload_gas_xyz('y', center=[0,0,0])
        gz= fload_gas_xyz('z', center=[0,0,0])
        gvx= fload_gas_v('x')
        gvy= fload_gas_v('y')
        gvz= fload_gas_v('z')
        gmass= fload_gas_mass(1)

        ; cold gas
        ; -----------
        ;cgx= fload_gas_xyz('x', center=[0,0,0])
        ;cgy= fload_gas_xyz('y', center=[0,0,0])
        ;cgz= fload_gas_xyz('z', center=[0,0,0])
        ;cgvx= fload_gas_v('x')
        ;cgvy= fload_gas_v('y')
        ;cgvz= fload_gas_v('z')
        ;cgmass= fload_gas_mass(1)


	;if twobhs eq 0 then begin
	   ; one center
	   ; -----------
           st_sigma_1[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 1.0, /usefixedRe, $
					avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)
           st_sigma_2[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 2.0, /usefixedRe, $
					avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)
           st_sigma_3[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 3.0, /usefixedRe, $
					avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)
           st_sigma_4[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 4.0, /usefixedRe, $
					avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)
           st_sigma_5[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 5.0, /usefixedRe, $
					avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)
           st_sigma_6[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 6.0, /usefixedRe, $
					avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)


	   ; gas
	   g_sigma_1[i]= process_2d_sigma_avg(gx, gy, gz, gvx, gvy, gvz, gmass, 1.0, /usefixedRe, $
                                        avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)
	   g_sigma_2[i]= process_2d_sigma_avg(gx, gy, gz, gvx, gvy, gvz, gmass, 2.0, /usefixedRe, $
                                        avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)
	   g_sigma_3[i]= process_2d_sigma_avg(gx, gy, gz, gvx, gvy, gvz, gmass, 3.0, /usefixedRe, $
                                        avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)
	   g_sigma_4[i]= process_2d_sigma_avg(gx, gy, gz, gvx, gvy, gvz, gmass, 4.0, /usefixedRe, $
                                        avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)
	   g_sigma_5[i]= process_2d_sigma_avg(gx, gy, gz, gvx, gvy, gvz, gmass, 5.0, /usefixedRe, $
                                        avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)
	   g_sigma_6[i]= process_2d_sigma_avg(gx, gy, gz, gvx, gvy, gvz, gmass, 6.0, /usefixedRe, $
                                        avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)

	;endif else begin
	   ; gal 1
	   ;Reff= process_2d_halfmass(sx, sy, smass, center=center_1, /do_quarter_mass)
	   ;as_sigma_xy[i]= process_2d_sigma(sx, sy, svz, 1.0*Reff, center= center_1)
	   ;as_sigma_xz[i]= process_2d_sigma(sx, sz, svy, 1.0*Reff, center= [center_1[0],center_1[2]])
	   ;as_sigma_yz[i]= process_2d_sigma(sy, sz, svx, 1.0*Reff, center= [center_1[1],center_1[2]])
	   ;as_sigma_avg[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 1.0*Reff, avg_Reff=avg_Reff, $
	;			center= center_1, sigmaerr=sigmaerr, refferr=refferr, /do_quarter_mass)

	   ;as_sigma_err[i]= sigmaerr
           ;as_re_avg[i]= avg_Reff
           ;as_re_err[i]= refferr
	   ;print, "#1: effective radius= ",Reff," (xy) ", avg_Reff, " (avg)"

	   ; gas
	   ;g_sigma_avg[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 1.0*avg_Reff, avg_Reff=avg_Reff, $
           ;                     center= center_1, sigmaerr=sigmaerr, refferr=refferr, /usefixedRe, /do_quarter_mass)

           ;g_sigma_err[i]= sigmaerr
           ;g_re_avg[i]= avg_Reff
           ;g_re_err[i]= refferr


	   ; gal 2
	   ;Reff= process_2d_halfmass(sx, sy, smass, center=center_2, /do_quarter_mass)
	   ;as_2_sigma_xy[i]= process_2d_sigma(sx, sy, svz, 1.0*Reff, center= center_2)
	   ;as_2_sigma_xz[i]= process_2d_sigma(sx, sz, svy, 1.0*Reff, center= [center_2[0],center_2[2]])
	   ;as_2_sigma_yz[i]= process_2d_sigma(sy, sz, svx, 1.0*Reff, center= [center_2[1],center_2[2]])
	   ;as_2_sigma_avg[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 1.0*Reff, avg_Reff=avg_Reff, $
	;			center= center_2, sigmaerr=sigmaerr, refferr=refferr, /do_quarter_mass)
	   ;as_2_sigma_err[i]= sigmaerr

	   ;as_2_re_avg[i]= avg_Reff
	   ;as_2_re_err[i]= refferr

	   ;print, "#2: effective radius= ",Reff," (xy) ", avg_Reff, " (avg)"
	;endelse

endfor
        
; ----------------------------------------


; print gas sigma information to file
; -------------------------------------
openw, 1, frun+'/sigma_gas.txt', ERROR=err

printf, 1, "#   sigma_gas.txt"
printf, 1, "#   all velocity dispersion measurements in units of km/sec"
printf, 1, "#        "
printf, 1, "# time      1 kpc      2 kpc      3 kpc      4 kpc      5 kpc      6 kpc"
printf, 1, "#(Gyr)        avg        avg        avg        avg        avg        avg"
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"  ",6(F9.4,"  "))', $
                time[i], g_sigma_1[i], g_sigma_2[i], g_sigma_3[i], g_sigma_4[i], $
                        g_sigma_5[i], g_sigma_6[i]
endfor
close, 1



; print star sigma information to file
; --------------------------------------
openw, 1, frun+'/sigma_stars.txt', ERROR=err

printf, 1, "#   sigma_stars.txt"
printf, 1, "#   all velocity dispersion measurements in units of km/sec"
printf, 1, "#   "
printf, 1, "# time      1 kpc      2 kpc      3 kpc      4 kpc      5 kpc      6 kpc"
printf, 1, "#(Gyr)        avg        avg        avg        avg        avg        avg"
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"  ",6(F9.4,"  "))', $
                	time[i], st_sigma_1[i], st_sigma_2[i], st_sigma_3[i], st_sigma_4[i], $
                        st_sigma_5[i], st_sigma_6[i]
endfor
close, 1





; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end





;=========================================================
;
;  3D Re
;
;=========================================================

pro time_3D_re, frun


if not keyword_set(frun) then begin
	print, " "
	print, " time_3D_re, frun"
	print, " "
	print, " "
	return
endif



; this assumes they are ordered
;  0 through x


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])



;  time sequence of quantities
; ------------------------------
time= fltarr(nsnaps)
as_re= fltarr(nsnaps)
as_2_re= fltarr(nsnaps)
g_re= fltarr(nsnaps)
g_2_re= fltarr(nsnaps)



; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------
snapi= 0

for i=0,nsnaps-1 do begin

	print, "--------------------------------------"

	if keyword_set(first_n_last_only) and i ne 0 then snapi= lastsnapnum-1

	; open snapshot
	repeat begin
	    ;ok=fload_snapshot_bh(frun,snapi)
	    ok=fload_snapshot_bh(frun,snapi,/nopot_in_snap)
	    ;ok=fload_snapshot(frun,snapi)
	    snapi= snapi+1
	endrep until (ok eq 0)

	; what time is it?
	time[i]= fload_time(1)


	; get blackhole id's
	bhid= fload_blackhole_id(1)
	n_bh= n_elements(bhid)
	sortbymass= bsort(bhm,/reverse)
	bhm= bhm(sortbymass)
	bhid= bhid(sortbymass)
	for ii=0, n_bh-1 do begin
		print, "Blackhole ID  "+strcompress(string(bhid[ii]),/remove_all)+"    Mass= ",bhm[ii]
		;if ii ge 1 then twobhs= 1
	endfor


	; get centers, for the moment we do this by hand, but could
	; generalize this in the future to search centers.txt file
	;twobhs= 0
	;if fload_npart(5) gt 1 then twobhs= 1

	;if twobhs eq 1 then begin
	if n_bh ge 2 then begin
             ;center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=170001)
             ;center_2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=340002)
             ;center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
             ;center_2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
             center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid[0])
             center_2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid[1])
	     rdiff= (center_1[0]-center_2[0])*(center_1[0]-center_2[0])  +  $
	            (center_1[1]-center_2[1])*(center_1[1]-center_2[1])  +  $
	            (center_1[2]-center_2[2])*(center_1[2]-center_2[2])
	     rdiff= sqrt(rdiff)
	     if rdiff lt 5.0 then twobhs= 0
	     print, "bh 1: ", bhid[0]
	     print, "bh 2: ", bhid[1]
	     print, "center 1: ", center_1
	     print, "center 2: ", center_2, "  (rdiff= ",rdiff,")"
        endif else begin
             ;bhid= fload_blackhole_id(1)
             ;center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
             center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid[0])
             center_2= center_1
	     print, "center: ", center_1
        endelse



xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



	;if twobhs eq 0 then begin
	if n_bh eq 1 then begin
	   ; one center
           Reff= process_2d_halfmass(sx, sy, smass, center=center_1)    ; Hey, this procedure is in process_directory
           as_sigma_xy[i]= process_2d_sigma(sx, sy, svz, 1.0*Reff, center= center_1)
           as_sigma_xz[i]= process_2d_sigma(sx, sz, svy, 1.0*Reff, center= [center_1[0],center_1[2]])
           as_sigma_yz[i]= process_2d_sigma(sy, sz, svx, 1.0*Reff, center= [center_1[1],center_1[2]])
           as_sigma_avg[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 1.0*Reff, $
					avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)

	   as_sigma_err[i]= sigmaerr
	   as_re_avg[i]= avg_Reff
	   as_re_err[i]= refferr
	   print, "effective radius= ",Reff," (xy) ", avg_Reff, " (avg)"

	   ; gas
	   g_sigma_avg[i]= process_2d_sigma_avg(gx, gy, gz, gvx, gvy, gvz, gmass, 1.0*avg_Reff, /usefixedRe, $
                                        avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)
           g_sigma_err[i]= sigmaerr
           g_re_avg[i]= avg_Reff
           g_re_err[i]= refferr

	endif else begin
	   ; gal 1
	   Reff= process_2d_halfmass(sx, sy, smass, center=center_1, /do_quarter_mass)
	   as_sigma_xy[i]= process_2d_sigma(sx, sy, svz, 1.0*Reff, center= center_1)
	   as_sigma_xz[i]= process_2d_sigma(sx, sz, svy, 1.0*Reff, center= [center_1[0],center_1[2]])
	   as_sigma_yz[i]= process_2d_sigma(sy, sz, svx, 1.0*Reff, center= [center_1[1],center_1[2]])
	   as_sigma_avg[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 1.0*Reff, avg_Reff=avg_Reff, $
				center= center_1, sigmaerr=sigmaerr, refferr=refferr, /do_quarter_mass)

	   as_sigma_err[i]= sigmaerr
           as_re_avg[i]= avg_Reff
           as_re_err[i]= refferr
	   print, "#1: effective radius= ",Reff," (xy) ", avg_Reff, " (avg)"

	   ; gas
	   g_sigma_avg[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 1.0*avg_Reff, avg_Reff=avg_Reff, $
                                center= center_1, sigmaerr=sigmaerr, refferr=refferr, /usefixedRe, /do_quarter_mass)

           g_sigma_err[i]= sigmaerr
           g_re_avg[i]= avg_Reff
           g_re_err[i]= refferr


	   ; gal 2
	   Reff= process_2d_halfmass(sx, sy, smass, center=center_2, /do_quarter_mass)
	   as_2_sigma_xy[i]= process_2d_sigma(sx, sy, svz, 1.0*Reff, center= center_2)
	   as_2_sigma_xz[i]= process_2d_sigma(sx, sz, svy, 1.0*Reff, center= [center_2[0],center_2[2]])
	   as_2_sigma_yz[i]= process_2d_sigma(sy, sz, svx, 1.0*Reff, center= [center_2[1],center_2[2]])
	   as_2_sigma_avg[i]= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 1.0*Reff, avg_Reff=avg_Reff, $
				center= center_2, sigmaerr=sigmaerr, refferr=refferr, /do_quarter_mass)
	   as_2_sigma_err[i]= sigmaerr

	   as_2_re_avg[i]= avg_Reff
	   as_2_re_err[i]= refferr

	   print, "#2: effective radius= ",Reff," (xy) ", avg_Reff, " (avg)"
	endelse

endfor
        
; ----------------------------------------


; print 3D R_e information to file
; -------------------------------
openw, 1, frun+'/R_e_3D.txt', ERROR=err

printf, 1, "#   R_e_3D.txt"
printf, 1, "#   effective radius in kpc h-1"
printf, 1, "#   "
printf, 1, "#  time    -- stars --   -- stars --   --  gas  --   --  gas  -- " 
printf, 1, "#(Gyr/h)   --    1  --   --    2  --   --   1   --   --   2   --"
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F8.3,"   ",4(F9.4,"   "))', $
				time[i], as_re_avg[i], as_re_err[i], $
					as_2_re_avg[i], as_2_re_err[i], $
					g_re_avg[i], g_re_err[i]
endfor
close, 1





; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end





;--------------------------------------------------------------------------
;==========================================================================
;--------------------------------------------------------------------------




;
;
;  single snapshot information
;
;  needs to be relaxed remnant (one center)
;

pro onesnap_m_sigma_re, frun, snapnum


if not keyword_set(frun) then begin
	print, " "
	print, " onesnap_m_sigma_re, frun, snapnum"
	print, " "
	print, " "
	return
endif



print, "--------------------------------------"

;ok=fload_snapshot_bh(frun,snapnum)
ok=fload_snapshot(frun,snapnum)

; what time is it?
time= fload_time(1)


center_1= fload_center_alreadycomp(1)




	; all luminous material
	; -----------------------
	sx= fload_allstars_xyz('x', center=[0,0,0])
	sy= fload_allstars_xyz('y', center=[0,0,0])
	sz= fload_allstars_xyz('z', center=[0,0,0])
             
	svx= fload_allstars_v('x')
	svy= fload_allstars_v('y') 
	svz= fload_allstars_v('z')

	smass= fload_allstars_mass(1)

        ; gaseous material
        ; -----------------------
	if fload_npart(0) gt 0 then begin
        	gx= fload_gas_xyz('x', center=[0,0,0])
        	gy= fload_gas_xyz('y', center=[0,0,0])
        	gz= fload_gas_xyz('z', center=[0,0,0])
        	gvx= fload_gas_v('x')
        	gvy= fload_gas_v('y')
        	gvz= fload_gas_v('z')
        	gmass= fload_gas_mass(1)
	endif else begin
        	gx= [0.0]
        	gy= [0.0]
        	gz= [0.0]
        	gvx= [0.0]
        	gvy= [0.0]
        	gvz= [0.0]
        	gmass= [0.0]
	endelse



	   ; one center
           Reff= process_2d_halfmass(sx, sy, smass, center=center_1)    ; Hey, this procedure is in process_directory
           as_sigma_xy= process_2d_sigma(sx, sy, svz, 1.0*Reff, center= center_1)
           as_sigma_xz= process_2d_sigma(sx, sz, svy, 1.0*Reff, center= [center_1[0],center_1[2]])
           as_sigma_yz= process_2d_sigma(sy, sz, svx, 1.0*Reff, center= [center_1[1],center_1[2]])
           as_sigma_avg= process_2d_sigma_avg(sx, sy, sz, svx, svy, svz, smass, 1.0*Reff, $
					avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)

	   as_sigma_err= sigmaerr
	   as_re_avg= avg_Reff
	   as_re_err= refferr
	   print, "effective radius= ",Reff," (xy) ", avg_Reff, " (avg)"

	   ; gas
	   g_sigma_avg= process_2d_sigma_avg(gx, gy, gz, gvx, gvy, gvz, gmass, 1.0*avg_Reff, /usefixedRe, $
                                        avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)
           g_sigma_err= sigmaerr
           g_re_avg= avg_Reff
           g_re_err= refferr


	; all baryons
	bx= [sx, gx]
	by= [sy, gy]
	bz= [sz, gz]
	bvx= [svx, gvx]
	bvy= [svy, gvy]
	bvz= [svz, gvz]
	bmass= [smass, gmass]

	Reff= process_2d_halfmass(bx, by, bmass, center=center_1)    ; Hey, this procedure is in process_directory
	ab_sigma_xy= process_2d_sigma(bx, by, bvz, 1.0*Reff, center= center_1)
	ab_sigma_xz= process_2d_sigma(bx, bz, bvy, 1.0*Reff, center= [center_1[0],center_1[2]])
	ab_sigma_yz= process_2d_sigma(by, bz, bvx, 1.0*Reff, center= [center_1[1],center_1[2]])
	ab_sigma_avg= process_2d_sigma_avg(bx, by, bz, bvx, bvy, bvz, bmass, 1.0*Reff, $
					avg_Reff=avg_Reff, center= center_1, sigmaerr=sigmaerr, refferr=refferr)

	ab_sigma_err= sigmaerr
	ab_re_avg= avg_Reff
	ab_re_err= refferr
	print, "effective radius= ",Reff," (xy) ", avg_Reff, " (avg)"



; print sigma information to file
; -------------------------------
openw, 1, frun+'/sigma_one.txt', ERROR=err

printf, 1, "#   sigma.txt"
printf, 1, "#   all velocity dispersion measurements in units of km/sec"
printf, 1, "#        "
printf, 1, "# time    -------------------- center 1 ---------------------  ------- gas --------"
printf, 1, "#(Gyr)         xy         xz         yz        avg        err        avg        err"
printf, 1, FORMAT= '(F6.3,"  ",7(F9.4,"  "))', $
                time, as_sigma_xy, as_sigma_xz, as_sigma_yz, as_sigma_avg, as_sigma_err, $
			g_sigma_avg, g_sigma_err
close, 1



; print total baryonic information to file
; -----------------------------------------
openw, 1, frun+'/sigma_allbs_one.txt', ERROR=err

printf, 1, "#   sigma.txt for all the baryons"
printf, 1, "#   all velocity dispersion measurements in units of km/sec"
printf, 1, "#        "
printf, 1, "# time    -------------------- center 1 ---------------------  -------- R_e -------"
printf, 1, "#(Gyr)         xy         xz         yz        avg        err        avg        err"
printf, 1, FORMAT= '(F6.3,"  ",7(F9.4,"  "))', $
                time, ab_sigma_xy, ab_sigma_xz, ab_sigma_yz, ab_sigma_avg, ab_sigma_err, $
			ab_re_avg, ab_re_err
close, 1



; print R_e information to file
; -------------------------------
openw, 1, frun+'/R_e_one.txt', ERROR=err

printf, 1, "#   R_e.txt"
printf, 1, "#   effective radius in kpc h-1"
printf, 1, "#   "
printf, 1, "# time   ---- center 1 -----   ------- gas --------"
printf, 1, "#(Gyr)        avg        err        avg        err"
printf, 1, FORMAT= '(F6.3,"  ",4(F9.4,"  "))', $
				time, as_re_avg, as_re_err, $
					g_re_avg, g_re_err
close, 1





; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end











;--------------------------------------------------------------------------
;==========================================================================
;--------------------------------------------------------------------------




