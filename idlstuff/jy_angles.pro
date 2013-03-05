;---------------------------------------------------------------
;
;
;
;---------------------------------------------------------------
pro jang, junk


frun= '/raid4/tcox/z5/jyCfg0.4'
snapnum= 35

desired_separation= 7.0   ; in kpc
tolerance= 0.7   ; in kpc


; ---------------------
; note, the above is in physical
; units, and below we'll do things
; in gadget units, so divide by h
desired_separation= desired_separation / 0.7
tolerance= tolerance / 0.7


; open the file up
;   and load variables
; ---------------------

ok= fload_snapshot_bh(frun, snapnum, /nopot_in_snap)


; get "center" - here it's the BH positions
if fload_npart(5) gt 1 then begin
	bhid= fload_blackhole_id(1)
	bhid1= bhid[0]
	bhid2= bhid[1]
	print, "Blackhole ID's: ", bhid1, bhid2

	center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
	center_2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
	;center_1= fload_blackhole_xyz('xyz',idtofollow=bhid1)
	;center_2= fload_blackhole_xyz('xyz',idtofollow=bhid2)
endif else begin
	print, " "
	print, " "
	print, "The Black Holes have merged - this is not the program you're looking for."
	print, " "
	print, " "
	return
endelse




; --------------------------------------------------------

;  read the list of theta's and phi's

anglefile= '/home/tcox/unitsphereangles.txt'
;anglefile= '/home/tcox/unitsphereangles.txt.test'

spawn, "wc "+anglefile,result
lines=long(result)
Nangles=lines(0)-5
if Nangles GT 0 then angle_data= fltarr(2,Nangles)

openr, 1, anglefile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, angle_data
close, 1


theta= angle_data[0,*]
phi= angle_data[1,*]


; --------------------------------------------------------
;
;   Loop through the angles, project and calculate separation
;

listof= fltarr(Nangles) & listof(*)= -1
listof_theta= theta
listof_phi= phi


for i= 0, Nangles-1 do begin

	rotate_phi= phi[i]
	rotate_theta= theta[i]

	if rotate_phi lt 1.0e-6 then rotate_phi= 1.0e-6
	if rotate_theta lt 1.0e-6 then rotate_theta= 1.0e-6

	; ------------------------------------------------------------------

	x= [center_1[0]]
	y= [center_1[1]]
	z= [center_1[2]]
        process_rotation, x, y, z, rotate_theta, rotate_phi, x1_new, y1_new, z1_new

	x= [center_2[0]]
	y= [center_2[1]]
	z= [center_2[2]]
        process_rotation, x, y, z, rotate_theta, rotate_phi, x2_new, y2_new, z2_new

	; ------------------------------------------------------------------

	xdiff= x1_new - x2_new
	ydiff= y1_new - y2_new
	projectedseparation= sqrt(xdiff*xdiff  + ydiff*ydiff)

	; do our check
	if (abs(desired_separation - projectedseparation) le tolerance) then listof[i]= +1
	print, desired_separation, projectedseparation
	;print, desired_separation, projectedseparation, tolerance

	; ------------------------------------------------------------------

endfor


; ------------------------------------------------------------------


snaplbl= '000'+strcompress(string(snapnum),/remove_all)
snaplbl= strmid(snaplbl,strlen(snaplbl)-3,3)
filename=frun+'/angles_with_sep_'+snaplbl+'.txt'
   


openw, 1, filename, ERROR=err
        
printf, 1, "#   "
printf, 1, "#  List of angles with projected separation of 7 kpc."
printf, 1, "#   "
printf, 1, "#  theta      phi     "
printf, 1, "#  (deg)     (deg)    "

; are there any close ones?
idx= where(listof gt 0)
Nlistof= n_elements(idx)
if idx(0) ne -1 then begin
	print, Nlistof," out of ", Nangles," are within range"
	listof_theta= listof_theta(idx)
	listof_phi= listof_phi(idx)
	for i=0,Nlistof-1 do begin
		printf, 1, FORMAT= '(2(F8.3,"  "))', $
                	listof_theta[i], listof_phi[i]
	endfor
endif else begin
	printf, 1,  " "
	printf, 1,  "no angles are within range"
	print, "no angles are within range"
endelse


close, 1




; ------------------------------------------------------------------



end









; =================================================================
; =================================================================





