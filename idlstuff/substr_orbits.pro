pro centers, frun

if not keyword_set(frun) then begin
	print, " "
	print, " centers, frun (all else set within)"
	print, " "
endif


; this assumes they are ordered
;  0 through


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])

if not keyword_set(startsnap) then startsnap= 0
if not keyword_set(endsnap) then endsnap= nsnaps


time= fltarr(nsnaps)
center=       fltarr(nsnaps,3)
comvel=       fltarr(nsnaps,3)
cenvel=       fltarr(nsnaps,3)


; std vc3c mergers
gal_sid= 1L
gal_npart= 2150001L

satnum='000'
sat_sid= 2150002L
sat_npart= 200001L




; center method
;use_bh_positions = 1
use_bh_positions = 0
use_density_peak = 1
;use_density_peak = 0


; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=startsnap,endsnap-1 do begin

	print, "--------------------------------------"


	; open snapshot
	;ok=fload_snapshot(frun,i)
	ok=fload_snapshot_bh(frun,i,/nopot_in_snap)
	;ok=fload_snapshot_bh(frun,i)


	; what time is it?
	time[i]= fload_time(1)


	if use_density_peak eq 1 then begin
	   center[i,*]= fload_1gal_center(sat_sid,sat_npart)
	   ;center[i,*]= fload_1gal_baryonic_center(sat_sid,sat_npart)
	   ;center[i,*]= fload_1gal_dm_center(sat_sid,sat_npart)
	endif


	; instead use BH positions
	; ------------------------
	if use_bh_positions eq 1 then begin

           ; get blackhole id's
           if i eq 0 then begin
        	bhid= fload_blackhole_id(1)
        	bhid1= bhid[0]
        	bhid2= bhid[1]
        	print, "Blackhole ID's: ", bhid1, bhid2
           endif

	   if fload_npart(5) gt 1 then begin
		center[i,*]= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
	   endif else begin
		bhid= fload_blackhole_id(1)
		center[i,*]= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
		;center[i,*]= fload_blackhole_xyz('xyz',idtofollow=bhid)
	   endelse
	endif



	; center-of-mass velocity, and central density velocity
	;   both of these need the particle numbers to be set
	;   correctly.

	comvel[i,*]= fload_1gal_comvel(sat_sid,sat_npart,center=center[i,*],rfact=1.0)
	;comvel[i,*]= fload_comvel(1,center=center[i,*],rfact=0.1)

	;cenvel[i,*]= fload_1gal_comvel(sat_sid,sat_npart,center=center[i,*],rfact=0.05)
	cenvel[i,*]= fload_comvel(1,center=center[i,*],rfact=0.001)

endfor
        


;===============================================================================



; ----------------------------------------

; write centers to file
write_centerpositions_1, time, center, filename=frun+"/sat_"+satnum+"_center.txt"

; write center-of-mass velocities to file
write_centerpositions_1, time, comvel, filename=frun+"/sat_"+satnum+"_comvel.txt"

; write center velocity to file (may be different than comvel)
write_centerpositions_1, time, cenvel, filename=frun+"/sat_"+satnum+"_cenvel.txt"

; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end




;===============================================================================



pro write_centerpositions_1, time, center, filename=filename



openw, 1, filename

; galaxy 1
printf, 1, 'satellite'

lines= n_elements(time)
printf, 1, lines

data= fltarr(4,lines)
data[0,*]= time
c1= transpose(center)
data[1:3,*]= c1
printf, 1, data

close, 1


end








;===============================================================================



pro read_centerpositions_1, time, center, filename=filename
        
        
        
; read this with the following
openr, 1, filename
junk= ''
lines= 0

; galaxy 1
readf, 1, junk
readf, 1, lines
data= fltarr(4,lines)
readf, 1, data  
time= data[0,*] 
gal1_xyz= data[1:3,*]
center= gal1_xyz
                
close, 1   
           

end        











;===============================================================================



pro find_rdiffs, frun

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



;
