




pro sinfo, frun, snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "sinfo, frun, snapnum"
   print, "  "
   print, "  "
   return
endif





print, "--------------------------------------"


;ok=fload_snapshot_bh(frun,snapnum)
ok=fload_snapshot_bh(frun,snapnum, /nopot_in_snap)

; what time is it?
;time[si]= fload_time(1)


print, " ---------------------"
print, "  DETERMINE: center   "
print, " ---------------------"


; -------------------------------
;
;
; determine_center
;
;   (make sure its the same as
;      sph info - it is)
; -------------------------------
center= [0,0,0]
center_bh= fload_center_alreadycomp(1)

; two black holes
if fload_npart(5) gt 1 then begin
        bhid= fload_blackhole_id(1)
        bhid1= bhid[0]
        bhid2= bhid[1]
        print, "Blackhole ID's: ", bhid1, bhid2
        center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
        center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
        center_bh= center1
        ;center_bh= center2
endif

; one black hole
if fload_npart(5) eq 1 then begin
        bhid= fload_blackhole_id(1)
        center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
endif




; grab luminosities
; ----------------------
ndisk= fload_npart(2)				; disk particles
nbulge= fload_npart(3)				; bulge particles
nstars= fload_npart(4)
npart= long(ndisk) + long(nbulge) + long(nstars)
print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars
TTime= float(fload_time(1))
N= npart
m= 1.0e+10*fload_allstars_mass(1)
age=fload_allstars_age(1)
age=float(TTime-age)
zmets=fload_allstars_z(1)

x= fload_allstars_xyz('x', center=center_bh)
y= fload_allstars_xyz('y', center=center_bh)
z= fload_allstars_xyz('z', center=center_bh)





; get the luminosities
;  - in units of solar luminosities
print, "load luminosities"
load_all_stellar_luminosities, N, TTime, m, age, zmets, $
	Lum_Bol, Lum_U, Lum_B, Lum_V, Lum_R, Lum_I, Lum_J, Lum_H, Lum_K, $
	Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z   ;, $
	;/notsolar

	; trap for any NaN's
	;idx=where(finite(Lum_B) eq 0)
	idx=where(finite(Lum_K) eq 0)
	if idx(0) ne -1 then begin
		Lum_Bol(idx)= 100.0 & Lum_U(idx)= 100.0 & Lum_B(idx)= 100.0 & Lum_V(idx)= 100.0 & Lum_R(idx)= 100.0
		Lum_I(idx)= 100.0 & Lum_J(idx)= 100.0 & Lum_H(idx)= 100.0 & Lum_K(idx)= 100.0
		Lum_sdss_u(idx)= 100.0 & Lum_sdss_g(idx)= 100.0 & Lum_sdss_r(idx)= 100.0
		Lum_sdss_i(idx)= 100.0 & Lum_sdss_z(idx)= 100.0
	print, "Fixed ",n_elements(idx), " luminosities when trapping for NaN values."
	endif

	print, "Bolometric luminosity= ", total(Lum_Bol)
	print, "Total K band luminosity= ", total(Lum_K)
	print, "Total B band luminosity= ", total(Lum_B)




print, "--------------------------------------"
print, "  order by radius"


radius= x*x + y*y + z*z
sortr= sort(radius)
x= x(sortr)
y= y(sortr)
z= z(sortr)
Lum_Bol= Lum_Bol(sortr)


print, "--------------------------------------"




; ------------------------------
;  Done - now write text file
; ------------------------------
print, "Writing to file: starinfo.txt"
fadd= ''
openw, 1, 'starinfo'+fadd+'.txt', ERROR=err

printf, 1, "#   starinfo.txt (+fadd, maybe)"
printf, 1, "#   frun=  "+frun
printf, 1, "#   time (Gyr/h)=  "+fload_timelbl(1,4,/noteq)
printf, 1, "#   star          x           y           z        luminosity"
printf, 1, "#  number      (kpc/h)     (kpc/h)     (kpc/h)        (solar) "
for i=0L,npart-1 do begin
        printf, 1, FORMAT= '("   ",F7.0,"   ",3(F9.4,"   "),F12.1)', $
                i, x[i], y[i], z[i], Lum_Bol[i]
endfor
close, 1


print, " "


end








; ================================================================================








