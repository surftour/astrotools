


pro testit, junk

do_one_angle, "data/ds/vc3vc3i", 0.0, 0.0, snapnum=19, filename="test.eps", xlen=7.0, $
					R_e= R_e, $
                                        Vrot_maj=Vrot_maj, avgVrot_maj=avgVrot_maj, $
                                        Vrot_min=Vrot_min, avgVrot_min=avgVrot_min, $
                                        pa=pa, Sig=sig, $
                                        fit_a=fit_a, fit_b=fit_b, $
                                        fit_ellip=fit_ellip, a4diva=a4diva


print, "R_e= ", R_e
print, "Vrot_maj= ", Vrot_maj
print, "avgVrot_maj= ", avgVrot_maj
print, "Vrot_min= ", Vrot_min
print, "avgVrot_min= ", avgVrot_min
print, "pa= ", pa
print, "Sig= ", Sig
print, "fit_a= ", fit_a
print, "fit_b= ", fit_b
print, "a4diva= ", fit_ellip
print, "a4diva= ", a4diva


end



; -----------------------------------------------------------------------------



pro test_1, junk

pa= 0.0

;
; this test didn't yield that much pertinent, excet that measuring the inner versus outer
; kinematics can be a tricky business.  in general, inner measurements give more velocity,
; but there is a large scatter.  in the end, i liked the tack of taking two separate
; images with different fields of view and measure the kinematics from each separate.
;
;
;


do_one_angle, "data/ds/vc3vc3l", 0.0, 0.0, snapnum=20, filename="temp_48.eps", xlen= 8.0, pa=pa, reff=4.0, $
;do_one_angle, "data/ds/vc3vc3l", 0.0, 0.0, snapnum=20, filename="temp_480.eps", xlen= 8.0, pa=pa, reff=4.0, $
;do_one_angle, "data/ds/vc3vc3l", 0.0, 0.0, snapnum=20, filename="temp_kin.eps", xlen= 8.0, pa=pa, $
;do_one_angle, "data/ds/vc3vc3l", 0.0, 0.0, snapnum=20, filename="temp_kin2.eps", xlen= 15.0, pa=pa, $
;do_one_angle, "data/ds/vc3vc3l", 0.0, 0.0, snapnum=20, filename="temp_kin_inner.eps", xlen= 15.0, pa=pa, /innerkinematics, $
;do_one_angle, "data/ds/vc3vc3l", 0.0, 0.0, snapnum=20, filename="temp_kin_outer.eps", xlen= 15.0, pa=pa, /outerkinematics, $
;do_one_angle, "data/ds/vc3vc3l", 0.0, 0.0, snapnum=20, filename="temp_kin_masswt.eps", xlen= 3.5, pa=pa, /masswtkinematics, $
;do_one_angle, "data/ds/vc3vc3l", 0.0, 0.0, snapnum=20, filename="temp_kin_smallr.eps", xlen= 3.5, pa=pa, $
;do_one_angle, "data/ds/vc3vc3l", 0.0, 0.0, snapnum=20, filename="temp_kin_smallr2.eps", xlen= 4.5, pa=pa, $
;do_one_angle, "data/ds/vc3vc3l", 0.0, 0.0, snapnum=20, filename="temp_kin_smallr_inner.eps", xlen= 4.5, pa=pa, /innerkinematics, $
					R_e= R_e, $
                                        Vrot_maj=Vrot_maj, avgVrot_maj=avgVrot_maj, $
                                        Vrot_min=Vrot_min, avgVrot_min=avgVrot_min, $
                                        Sig=sig, $
                                        fit_a=fit_a, fit_b=fit_b, $
                                        fit_ellip=fit_ellip, a4diva=a4diva


print, "R_e= ", R_e
print, "Vrot_maj= ", Vrot_maj
print, "avgVrot_maj= ", avgVrot_maj
print, "Vrot_min= ", Vrot_min
print, "avgVrot_min= ", avgVrot_min
print, "pa= ", pa
print, "Sig= ", Sig
print, "fit_a= ", fit_a
print, "fit_b= ", fit_b
print, "fit_ellip= ", fit_ellip
print, "a4diva= ", a4diva


end




; -----------------------------------------------------------------------------



pro test_all_3, junk

pa= 0.0

kinnum= 48   ; for the moment, this is the pixel number

istring= strcompress(string(kinnum), /remove_all)

;frun= "data/ds/vc3vc3l"
frun= "data/ds/vc3vc3f"
fitsdir= "./"

do_one_angle, frun, 0.0, 0.0, snapnum=20, filename="test_"+istring+".eps", xlen= 8.0, pa=pa, reff=4.0, $
					R_e= R_e, $
                                        Vrot_maj=Vrot_maj, avgVrot_maj=avgVrot_maj, $
                                        Vrot_min=Vrot_min, avgVrot_min=avgVrot_min, $
                                        Sig=sig, $
                                        fit_a=fit_a, fit_b=fit_b, $
                                        fit_ellip=fit_ellip, a4diva=a4diva


print, "R_e= ", R_e
print, "Vrot_maj= ", Vrot_maj
print, "avgVrot_maj= ", avgVrot_maj
print, "Vrot_min= ", Vrot_min
print, "avgVrot_min= ", avgVrot_min
print, "pa= ", pa
print, "Sig= ", Sig
print, "fit_a= ", fit_a
print, "fit_b= ", fit_b
print, "fit_ellip= ", fit_ellip
print, "a4diva= ", a4diva



do_lambda_profile_4, frun, kinnum, snapnum=20, filename="test_lambda_"+istring+".eps", fitsdir=fitsdir, /loadedsnap, $
			sdfitsf="test_"+istring+".fits", vfitsf="test_"+istring+"_vel.fits", dfitsf="test_"+istring+"_dis.fits"


;do_isophote_profile, frun, kinnum, snapnum=20, filename="test_isoprof_"+istring+".eps", fitsdir=fitsdir, /loadedsnap, $
;			sdfitsf="test_"+istring+".fits", vfitsf="test_"+istring+"_vel.fits", dfitsf="test_"+istring+"_dis.fits"




end



; -----------------------------------------------------------------------------




pro testvel, junk


frun= "data/isolated/d3"

do_one_angle, frun, 0.0, 0.0, snapnum= 5, filename="test.eps", xlen= 10.0, reff=3.0
do_lambda_profile_4, frun, '0', snapnum= 5, filename="test_lambda.eps", /loadedsnap, $
				sdfitsf="test_sd.fits", $
				vfitsf= "test_vel.fits", $
				dfitsf= "test_dis.fits"

do_one_angle, frun, 90.0, 0.0, snapnum= 5, filename="test_90.eps", xlen= 10.0, reff=3.0
do_lambda_profile_4, frun, '0', snapnum= 5, filename="test_90_lambda.eps", /loadedsnap, $
				sdfitsf="test_90_sd.fits", $
				vfitsf= "test_90_vel.fits", $
				dfitsf= "test_90_dis.fits"

do_one_angle, frun, 90.0, 90.0, snapnum= 5, filename="test_9090.eps", xlen= 10.0, reff=3.0
do_lambda_profile_4, frun, '0', snapnum= 5, filename="test_9090_lambda.eps", /loadedsnap, $
				sdfitsf="test_9090_sd.fits", $
				vfitsf= "test_9090_vel.fits", $
				dfitsf= "test_9090_dis.fits"


end






; -----------------------------------------------------------------------------




pro lorentest, junk


frun= "data/ds/vc3vc3i"
snapnum= 19
;
; long axis view
;do_one_angle, frun, -17.7, 46.4, snapnum= snapnum, filename=frun+"/snap_019_slitmap_loren.eps", xlen= 10.0, reff=3.0
;
; intermediate axis view
;do_one_angle, frun, -167.6, 47.9, snapnum= snapnum, filename=frun+"/snap_019_slitmap_loren.eps", xlen= 6.0, reff=3.0
;do_one_angle, frun, -167.6, 47.9, snapnum= 21, filename=frun+"/snap_021_slitmap_loren.eps", xlen= 6.0, reff=3.0
;do_one_angle, frun, -167.6, 47.9, snapnum= 23, filename=frun+"/snap_023_slitmap_loren.eps", xlen= 6.0, reff=3.0
;do_one_angle, frun, -167.6, 47.9, snapnum= 25, filename=frun+"/snap_025_slitmap_loren.eps", xlen= 6.0, reff=3.0
;do_one_angle, frun, -167.6, 47.9, snapnum= 27, filename=frun+"/snap_027_slitmap_loren.eps", xlen= 6.0, reff=3.0
;do_one_angle, frun, -167.6, 47.9, snapnum= 29, filename=frun+"/snap_029_slitmap_loren.eps", xlen= 6.0, reff=3.0
;do_one_angle, frun, -167.6, 47.9, snapnum= 31, filename=frun+"/snap_031_slitmap_loren.eps", xlen= 6.0, reff=3.0
;
; short axis view
;do_one_angle, frun, 87.6, 74.3, snapnum= snapnum, filename=frun+"/snap_019_slitmap_loren.eps", xlen= 10.0, reff=3.0
;
;
;
;do_lambda_profile_4, frun, '0', snapnum= snapnum, filename=frun+"/snap_019_lambda.eps", /loadedsnap, $
;                                sdfitsf=frun+"/snap_019_slitmap_loren_sd.fits", $
;                                vfitsf= frun+"/snap_019_slitmap_loren_vel.fits", $
;                                dfitsf= frun+"/snap_019_slitmap_loren_dis.fits"


frun= "data/collisionless/cvc3vc3i"
snapnum= 23
;
; long axis view
;do_one_angle, frun, 37.7, 85.4, snapnum= snapnum, filename=frun+"/snap_023_slitmap_loren.eps", xlen= 10.0, reff=3.0
;
; intermediate axis view
;do_one_angle, frun, -55.7, 55.9, snapnum= snapnum, filename=frun+"/snap_023_slitmap_loren.eps", xlen= 10.0, reff=3.0
do_one_angle, frun, -55.7, 55.9, snapnum= 25, filename=frun+"/snap_025_slitmap_loren.eps", xlen= 10.0, reff=3.0
do_one_angle, frun, -55.7, 55.9, snapnum= 27, filename=frun+"/snap_027_slitmap_loren.eps", xlen= 10.0, reff=3.0
do_one_angle, frun, -55.7, 55.9, snapnum= 29, filename=frun+"/snap_029_slitmap_loren.eps", xlen= 10.0, reff=3.0
do_one_angle, frun, -55.7, 55.9, snapnum= 31, filename=frun+"/snap_031_slitmap_loren.eps", xlen= 10.0, reff=3.0
;
; short axis view
;do_one_angle, frun, -45.7, 146.1, snapnum= snapnum, filename=frun+"/snap_023_slitmap_loren.eps", xlen= 10.0, reff=3.0
;
;
;
;do_lambda_profile_4, frun, '0', snapnum= snapnum, filename=frun+"/snap_023_lambda.eps", /loadedsnap, $
;                                sdfitsf=frun+"/snap_023_slitmap_loren_sd.fits", $
;                                vfitsf= frun+"/snap_023_slitmap_loren_vel.fits", $
;                                dfitsf= frun+"/snap_023_slitmap_loren_dis.fits"




end








; -----------------------------------------------------------------------------
; -----------------------------------------------------------------------------
; -----------------------------------------------------------------------------
; -----------------------------------------------------------------------------
; -----------------------------------------------------------------------------
; -----------------------------------------------------------------------------



pro process_kinematics, frun, snapnum, xlen=xlen, $
			small_fov=small_fov, $
			med_fov=med_fov, $
			large_fov=large_fov, $
			outerkinematics=outerkinematics, $
			innerkinematics=innerkinematics, $
			masswtkinematics=masswtkinematics, $
			skipkinemetry=skipkinemetry



if not keyword_set(frun) then begin
        print, "  "
        print, " process_kinematics, frun, snapnum"
        print, "  "
        return
endif

if not keyword_set(snapnum) then begin

        ;print, "  "
        ;print, " process_kinematics, frun, snapnum"
        ;print, "  "
        ;return

	bhmtime= fload_bh_mergertime(frun)
        if bhmtime le 0 then begin
                ;print, " "
                ;print, " WARNING: can't determine bhmtime, setting to be 2.0"
                ;print, " "
                ;bhmtime= 2.0

                print, " "
                print, " WARNING: can't determine bhmtime, using centers.txt"
                print, " "
		bhmtime= fload_mergertime(frun)
        endif

        spawn, "/bin/ls "+frun+"/snapshot* | wc ",result
        totnumsnap=long(result[0])-1

	relaxationtime= 0.7  ; Gyr/h

	bhmtime_snaptime= -1
	bhmtime_snapnum= -1

        for i=0,totnumsnap do begin
            ;ok= fload_snapshot_bh(frun,i,/header_only)
            ok= fload_snapshot_bh(frun,i,/header_only,/arepo)

	    bhnum= fload_npart(5)
	    if i eq 0 then bhnum_lastsnap= bhnum

	    if bhnum ne bhnum_lastsnap then begin
		bhmtime_snaptime= fload_time(1)
		bhmtime_snapnum= i
	    endif

            ; grab snapnum when system is "relaxed"
            if (fload_time(1) gt (bhmtime+relaxationtime)) or (i eq totnumsnap) then begin
		snapnum= i
		break
	    endif

	    bhnum_lastsnap= bhnum

	endfor

	print, " "
	print, "         Timing Information "
	print, " ----------------------------------- "
	print, " fload_bh_mergertime= ", bhmtime
	print, " bhmtime_snaptime= ", bhmtime_snaptime
	print, " bhmtime_snapnum= ", bhmtime_snapnum
	print, " analyze at "+string(bhmtime_snaptime)+" + "+string(relaxationtime)+" Gyr/h later"
	print, " snapnum= ", snapnum
	print, " "


	if bhnum eq 2 then begin
		print, " "
		print, " PROBLEM: BH's never merged"
		print, " "
		;return
	endif

endif


ok=fload_snapshot_bh(frun, snapnum, /nopot_in_snap)
;ok=fload_snapshot_bh(frun, snapnum, /nopot_in_snap, /arepo)


; setup directory information
fitsdir= frun+'/snap_'+fload_finfo_exts(1)+'_info/'

; first, save copy of previous kinematics file
;spawn, "mv -i "+frun+"/kinematics.txt.old "+frun+"/kinematics.txt.1"
;spawn, "mv -f "+frun+'/snap_'+fload_finfo_exts(1)+'_info/ '+frun+'/snap_'+fload_finfo_exts(1)+'_info_2/'
;spawn, "mv -f "+frun+"/kinematics.txt "+frun+"/kinematics.txt.old"
;spawn, "mv -i "+fitsdir+" "+frun+'/snap_'+fload_finfo_exts(1)+'_info_2/'


; find R_e & set xlen for all images
reff= fload_allstars_3dhalfmassr(1)
xlen_orig= long(reff * 2.25 + 0.5)      ; inner regions
;xlen_orig= long(reff * 7.5 + 0.5)


if keyword_set(small_fov) then begin
	fitsdir= frun+'/snap_'+fload_finfo_exts(1)+'_small_fov/'
	;xlen_orig= long(reff + 0.5)
	xlen_orig= long(reff * 1.5)


	spawn, "mv -f "+frun+'/kinematics_small_fov.txt '+fitsdir
	spawn, "mv -f "+fitsdir+" "+frun+'/snap_'+fload_finfo_exts(1)+'_small_fov_trial1/'
;stop
endif

if keyword_set(large_fov) then begin
	fitsdir= frun+'/snap_'+fload_finfo_exts(1)+'_large_fov/'
	xlen_orig= long(reff * 7.5 + 0.5)
endif


if keyword_set(med_fov) then begin
	fitsdir= frun+'/snap_'+fload_finfo_exts(1)+'_med_fov/'
	xlen_orig= long(reff * 3.0 + 0.5)
endif


if keyword_set(outerkinematics) then begin
	xlen_orig= long(reff * 7.0)
endif


; now, make directory
print, "making directory ", fitsdir
spawn, "mkdir "+fitsdir


; --------------------------------------------------------

;  read the list of theta's and phi's

;anglefile= '/n/home/tcox/unitsphereangles.txt'
;anglefile= '/n/home/tcox/Documents/Data_UnitSphere/unitsphereangles.txt'
anglefile= '/n/home/tcox/Documents/Data_UnitSphere/unitsphereangles.txt.10'
;anglefile= '/n/home/tcox/Documents/Data_UnitSphere/unitsphereangles.txt.190'

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



; --------------------------------------------------------
;
;   Loop through the angles, project and calculate
;  the rotation, dispersion, size, ellipticity, and disky/boxy - ness
;

Re_i= fltarr(Nangles)
vrot_i_maj= fltarr(Nangles)
vrotavg_i_maj= fltarr(Nangles)
vrot_i_min= fltarr(Nangles)
vrotavg_i_min= fltarr(Nangles)
sig_i= fltarr(Nangles)
PA_i= fltarr(Nangles)
fit_a_i= fltarr(Nangles)
fit_b_i= fltarr(Nangles)
ellip_i= fltarr(Nangles)
a4diva_i= fltarr(Nangles)

;lambda_Re= fltarr(Nangles)
;lambda_5Re= fltarr(Nangles)
;
;ell_Re= fltarr(Nangles)
;posang_Re= fltarr(Nangles)
;a4a_Re= fltarr(Nangles)
;ell_5Re= fltarr(Nangles)
;posang_5Re= fltarr(Nangles)
;a4a_5Re= fltarr(Nangles)


for i= 0, Nangles-1 do begin

        print, " ------------------------------------------------------------------ "
	print, "    processing: "+string(i)
	print, "    theta= "+string(theta[i])
	print, "    phi= "+string(phi[i])
        print, " ------------------------------------------------------------------ "

        rotate_phi= phi[i]
        rotate_theta= theta[i]


	istring= strcompress(string(i), /remove_all)
        ;thisfilename= fitsdir+'/kinmap_'+istring+'.eps'    ; kind of funky becuase
        thisfilename= fitsdir+'/'+istring+'_slitmap.eps'    ; kind of funky becuase
                                                        ; contour_makepic assumes
                                                        ; the ending is eps


	if not keyword_set(xlen) then xlen= float(xlen_orig)

	pa= 0.0     ; need this so it resets the position angle for each slit

	do_one_angle, frun, rotate_theta, rotate_phi, snapnum=snapnum, filename=thisfilename, /loadedsnap, $
					R_e= R_e, reff=reff, $
					xlen= xlen, $
					outerkinematics=outerkinematics, $
					innerkinematics=innerkinematics, $
					masswtkinematics=masswtkinematics, $
					small_fov=small_fov, med_fov=med_fov, large_fov=large_fov, $
					Vrot_maj=Vrot_maj, avgVrot_maj=avgVrot_maj, $
					Vrot_min=Vrot_min, avgVrot_min=avgVrot_min, $
					pa=pa, Sig=sig, $
					fit_a=fit_a, fit_b=fit_b, $
					fit_ellip=fit_ellip, a4diva=a4diva, $
					fitsdir=fitsdir, istring=istring


	Re_i[i]= R_e
	vrot_i_maj[i]= Vrot_maj
	vrotavg_i_maj[i]= avgVrot_maj
	vrot_i_min[i]= Vrot_min
	vrotavg_i_min[i]= avgVrot_min
	PA_i[i]= pa
	Sig_i[i]= Sig
	fit_a_i[i]= fit_a
	fit_b_i[i]= fit_b
	ellip_i[i]= fit_ellip
	a4diva_i[i]= a4diva


	;------------------------------------------------


	kinnum= i
	;thisfilename= fitsdir+'/kinlambda_'+istring+'.eps'
	thisfilename= fitsdir+'/'+istring+'_lambda.eps'

	do_lambda_profile_4, frun, kinnum, snapnum=snapnum, filename=thisfilename, fitsdir=fitsdir, /loadedsnap
					;lam_Re=lam_Re, lam_5Re=lam_5Re


	;lambda_Re[i]= lam_Re
	;lambda_5Re[i]= lam_5Re


	;------------------------------------------------


	kinnum= i
	;thisfilename= fitsdir+'/isoprof_'+istring+'.eps'
	thisfilename= fitsdir+'/'+istring+'_isophotes.eps'

	do_isophote_profile, frun, kinnum, snapnum=snapnum, filename=thisfilename, fitsdir=fitsdir, /loadedsnap
					;ellip_Re=ellip_Re, pa_Re=pa_Re, a4diva_Re=a4diva_Re, $
					;ellip_5Re=ellip_5Re, pa_5Re=pa_5Re, a4diva_5Re=a4diva_5Re


	;ell_Re[i]= ellip_Re
	;posang_Re[i]= pa_Re
	;a4a_Re[i]= a4diva_Re
	;ell_5Re[i]= ellip_5Re
	;posang_5Re[i]= pa_5Re
	;a4a_5Re[i]= a4diva_5Re


	;------------------------------------------------

	;
	; now do a much smaller pixel version
	;

        thisfilename= fitsdir+'/'+istring+'_slitmap_48.eps'
        pa= 0.0

        do_one_angle, frun, rotate_theta, rotate_phi, snapnum=snapnum, filename=thisfilename, /loadedsnap, $
                                        R_e= R_e, reff=reff, $
                                        xlen= xlen, $
                                        outerkinematics=outerkinematics, $
                                        innerkinematics=innerkinematics, $
                                        masswtkinematics=masswtkinematics, $
                                        small_fov=small_fov, med_fov=med_fov, large_fov=large_fov, $
					pixels= 48L


	;------------------------------------------------


	if not keyword_set(skipkinemetry) then kinematics_kinemetry, frun, nbins= 30, istring=istring, fitsdir=fitsdir


	;------------------------------------------------


	epsfilename= fitsdir+'/'+istring+'_rprofs.eps'

	plot_compiled_radial_data, fitsdir, istring, filename=epsfilename


endfor


; ------------------------------------------------------------------


filename= frun+'/kinematics.txt'
if keyword_set(outerkinematics) then filename= frun+'/kinematics_outer.txt'
if keyword_set(innerkinematics) then filename= frun+'/kinematics_inner.txt'
if keyword_set(masswtkinematics) then filename= frun+'/kinematics_masswt.txt'
if keyword_set(small_fov) then filename= frun+'/kinematics_small_fov.txt'
if keyword_set(med_fov) then filename= frun+'/kinematics_med_fov.txt'
if keyword_set(large_fov) then filename= frun+'/kinematics_large_fov.txt'

openw, 1, filename, ERROR=err
        
printf, 1, "#   kinematics.txt, "+frun+" snap= "+string(snapnum)
printf, 1, "#    (all quantites are calculated for the half-mass isophote)"
printf, 1, "#                                                                                          -- Major Axis --    -- Minor Axis --  "
printf, 1, "#  theta      phi     R_e         a         b                                   Sig_0     V_rot    <V_rot>    V_rot    <V_rot>  "
printf, 1, "#  (deg)    (deg)  (kpc/h)   (kpc/h)   (kpc/h)     ellip       pa      a4/a    (km/s)    (km/s)     (km/s)   (km/s)     (km/s)  "

for i=0,Nangles-1 do begin 
        printf, 1, FORMAT= '(2(F8.3," "),4(F8.4,"  "),1(F8.2,"  "),1(F8.4,"  "),5(F8.2,"  "))', $
                theta[i], phi[i], $
		Re_i[i], $
		fit_a_i[i], fit_b_i[i], ellip_i[i], PA_i[i], a4diva_i[i], $
		sig_i[i], $
		vrot_i_maj[i], vrotavg_i_maj[i], vrot_i_min[i], vrotavg_i_min[i]
endfor 
close, 1



; ------------------------------------------------------------------
;
;
;filename= frun+'/kinematics_lambda.txt'
;if keyword_set(outerkinematics) then filename= frun+'/kinematics_lambda_outer.txt'
;if keyword_set(innerkinematics) then filename= frun+'/kinematics_lambda_inner.txt'
;if keyword_set(masswtkinematics) then filename= frun+'/kinematics_lambda_masswt.txt'
;if keyword_set(small_fov) then filename= frun+'/kinematics_lambda_small_fov.txt'
;if keyword_set(large_fov) then filename= frun+'/kinematics_lambda_large_fov.txt'
;
;openw, 1, filename, ERROR=err
;
;printf, 1, "#   kinematics_lambda.txt, "+frun+" snap= "+string(snapnum)
;printf, 1, "#   "
;printf, 1, "#                                  "
;printf, 1, "#  theta      phi   Lambda     Ellip       PA   100*a4/a    Lambda     Ellip       PA   100*a4/a"
;printf, 1, "#  (deg)    (deg)   (      ---      1R_e      ---      )    (      ---      5R_e      -----    )"
;
;for i=0,Nangles-1 do begin 
;        printf, 1, FORMAT= '(2(F8.3," "),8(F8.4,"  "))', $
;                theta[i], phi[i], $
;                lambda_Re[i], ell_Re[i], posang_Re[i], a4a_Re[i], $
;		lambda_5Re[i], ell_5Re[i], posang_5Re[i], a4a_5Re[i]
;endfor
;close, 1
;
;
;
; ------------------------------------------------------------------



end




;##############################################################################################
;##############################################################################################
;##############################################################################################





;=============================================================
;
; grid: 3x1
;
; 1) slit info
; 2) project gas or star density
; 3) velocity field
;
;=============================================================


pro do_one_angle, frun, rt, rp, pa=pa, snapnum=snapnum, loadedsnap=loadedsnap, $
					filename=filename, fitsdir=fitsdir, istring=istring, $
					xlen=xlen, $
					R_e= R_e, reff=reff, $
					outerkinematics=outerkinematics, $
					innerkinematics=innerkinematics, $
					masswtkinematics=masswtkinematics, $
					pixels=pixels, $
					small_fov=small_fov, $
					med_fov=med_fov, $
					large_fov=large_fov, $
					Vrot_maj=Vrot_maj, avgVrot_maj=avgVrot_maj, $
					Vrot_min=Vrot_min, avgVrot_min=avgVrot_min, $
					Sig=sig, $
					fit_a=fit_a, fit_b=fit_b, $
					fit_ellip=fit_ellip, a4diva=a4diva


if not keyword_set(frun) then begin
        print, " "
        print, " do_one_angle, frun, rotate_theta, rotate_phi,"
        print, "                          pa=pa, snapnum=snapnum,"
        print, "                          filename=filename"
        print, " "
        return
endif


if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(snapnum) then begin
        print, " "
        print, " do_one_angle, frun, rotate_theta, rotate_phi,"
        print, "                          pa=pa, snapnum=snapnum,"
        print, "                          filename=filename"
        print, " "
        return
endif



; --------------------------
;  Setup stuff
; --------------------------

if not keyword_set(filename) then filename='slit_ho.eps'

if not keyword_set(xlen) then xlen= 10.0

; view 0
if keyword_set(rt) then rotate_theta= rt else rotate_theta= 0.0
if keyword_set(rp) then rotate_phi= rp else rotate_phi= 0.0

if keyword_set(pa) then pa_fixed= pa else pa_fixed= -1




; slit dimensions
; -----------------
slit_bins = 50
; default
;slit_width = xlen*0.06
;slit_len = xlen*0.9
; trial 
slit_width = xlen*0.06
slit_len = xlen*0.9

; plot range
slitmax = xlen*1.1

vmax= 250.0
vscale= 200.0

; image size
;xlen= 1.2*slit_len
;xlen= 8.0


if not keyword_set(pixels) then pixels= 480L     ; standard size of images
;if not keyword_set(pixels) then pixels= 48L      ; useful for the kinemetry stuff
pixels= long(pixels)


fitstoo= 1
;fitstoo= 0


; start it up
; -------------


x0= 0.08
xs= (0.90-x0)/3.0   ; assumes 3 panels

x1= x0+xs
x2= x0+xs+xs
x3= x0+xs+xs+xs
x4= x0+xs+xs+xs+xs
    
y0= 0.02
y1= 0.98



; ---------------
;  Do it Up
; ---------------

if not keyword_set(loadedsnap) then ok=fload_snapshot_bh(frun,snapnum)


if keyword_set(masswtkinematics) then begin
	domaxrng= 1
endif

if keyword_set(outerkinematics) then begin
	;
	; for the moment, we'll keep original slit
	;
	;FindMassFactor= 0.75
endif



; get print stuff ready
; ----------------------
initialize_plotinfo, 1
        
;setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=21, newysize=10
;setup_plot_stuff, 'ps', filename=filename, colortable= 1, newxsize=21, newysize=10
setup_plot_stuff, 'ps', filename=filename, colortable= 0, newxsize=32, newysize=10
print, "writing: ", filename


; what center do we use?
; ------------------------

; option 1  (use what's computed at end of snapshot loader)
center= fload_center_alreadycomp(1)
print, "center (alreadycomp)= ", center

; option 2  (use BH position)
bhid= fload_blackhole_id(1)
;stop
;if bhid[0] eq 1 then begin 
if bhid[0] gt 0 then begin 
	bhid= long(max(bhid))   ; uses minimum bhid - then you don't need to know the BH ID
	center= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
	print, "Blackhole ID: ", bhid
	print, "Blackhole center (used as center): ", center
endif


orig_center= center



; Do Stars
;===========
;do_stars= 0
do_stars= 1
if do_stars eq 1 then begin
	x= fload_allstars_xyz('x', center=[0,0,0])
	y= fload_allstars_xyz('y', center=[0,0,0])
	z= fload_allstars_xyz('z', center=[0,0,0])
	vx= fload_allstars_v('x')
	vy= fload_allstars_v('y')
	vz= fload_allstars_v('z')
	mass= fload_allstars_mass(1)

	r2= (x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2])
	r= sqrt(r2)
	;hsml= 0.1 + 1.0*(r/reff)
	hsml= 0.1 + 0.5*(r/reff)^2

	;idx= where(r gt sqrt(2)*xlen)
	idx= where(hsml gt 1.5*xlen)
	hsml(idx)= 1.5*xlen

	;if keyword_set(small_fov) then hsml(*)= 0.2
endif


do_oldstars= 0
;do_oldstars= 1
if do_oldstars eq 1 then begin
        x= fload_oldstars_xyz('x', center=[0,0,0])
        y= fload_oldstars_xyz('y', center=[0,0,0])
        z= fload_oldstars_xyz('z', center=[0,0,0])
        vx= fload_oldstars_v('x')
        vy= fload_oldstars_v('y')
        vz= fload_oldstars_v('z')
        mass= fload_oldstars_mass(1)

        r2= (x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2])
        r= sqrt(r2)
        ;hsml= 0.1 + 1.0*(r/reff)
        hsml= 0.1 + 0.5*(r/reff)^2

        ;idx= where(r gt sqrt(2)*xlen)
        idx= where(hsml gt 1.5*xlen)
        hsml(idx)= 1.5*xlen
endif



do_newstars= 0
;do_newstars= 1
if do_newstars eq 1 then begin
        x= fload_newstars_xyz('x', center=[0,0,0])
        y= fload_newstars_xyz('y', center=[0,0,0])
        z= fload_newstars_xyz('z', center=[0,0,0])
        vx= fload_newstars_v('x')
        vy= fload_newstars_v('y')
        vz= fload_newstars_v('z')
        mass= fload_newstars_mass(1)

        r2= (x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2])
        r= sqrt(r2)
        ;hsml= 0.1 + 1.0*(r/reff)
        hsml= 0.1 + 0.5*(r/reff)^2

        ;idx= where(r gt sqrt(2)*xlen)
        idx= where(hsml gt 1.5*xlen)
        hsml(idx)= 1.5*xlen
endif




; Do Gas
;===========
do_gas= 0
;do_gas= 1
if do_gas eq 1 then begin
        x= fload_gas_xyz('x', center=[0,0,0])
        y= fload_gas_xyz('y', center=[0,0,0])
        z= fload_gas_xyz('z', center=[0,0,0])
        vx= fload_gas_v('x')
        vy= fload_gas_v('y')
        vz= fload_gas_v('z')
        mass= fload_gas_mass(1)
	hsml= fload_gas_hsml(1)
endif

        

;  middle - projected surface density + slit
; ---------------------------------------------
x_orig= x
y_orig= y
z_orig= z
center= orig_center
process_kin_image_and_overplot, x_orig, y_orig, z_orig, mass, xlen, x1, y0, x2, y1, $
	hsml=hsml, $
	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	center=center, /show_re, R_e=R_e, $ ;/showbhs, $
	set_maxden=10.0, set_dynrng=1.0e5, pixels=pixels, $
	fitstoo=fitstoo, filename=filename, $
	;msg1='stars', msg2=' ', /drawbox, xslit=slit_len, yslit=0.5*slit_width, $
	;msg1=' ', msg2=' ', /show_maj_slit, /show_min_slit, $
	msg1=' ', msg2=' ', /show_maj_slit, /show_xlen, $
	xslit=slit_len, yslit=0.5*slit_width, $
	;msg1=' ', msg2=' ', xslit=slit_len, yslit=0.5*slit_width, $
	fit_a=fit_a, fit_b=fit_b, fit_ellip=fit_ellip, a4diva=a4diva, $
	ReturnedMap=ReturnedMap, pa=pa, HalfMassSB=HalfMassSB, $
	FindMassFactor=FindMassFactor



if pa_fixed ne -1 then pa= pa_fixed


;  left - velocity field within slit
; ------------------------------------------
x_orig= x
y_orig= y
z_orig= z
vx_orig= vx
vy_orig= vy
vz_orig= vz
center= orig_center
process_kin_slit_and_overplot, x_orig, y_orig, z_orig, vx_orig, vy_orig, vz_orig, mass, $
        x0, y0, x1, y1, $
        center=center, /ditchsigma, $
        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        slit_len=slit_len, slit_width=slit_width, /show_slitlen, $
        slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
	Rin=0.5*R_e, Sig=Sig, $
	Vrot_maj=Vrot_maj, avgVrot_maj=avgVrot_maj, $
	Vrot_min=Vrot_min, avgVrot_min=avgVrot_min, $
        pa=pa, fitsdir=fitsdir, istring=istring, $
	outerkinematics=outerkinematics, $
	innerkinematics=innerkinematics, $
	masswtkinematics=masswtkinematics



; labels, if we want 'em
thetalbl= strcompress(string(rotate_theta),/remove_all)
digs= 3
if rotate_theta gt 100.0 then digs= 5
if rotate_theta gt 10.0 then digs= 4
if rotate_theta lt 0.0 then digs= 4
if rotate_theta lt 10.0 then digs= 5
if rotate_theta lt 100.0 then digs= 6
thetalbl= '!7h!6='+strmid(thetalbl,0,digs)
xyouts, x0+0.02, y1-0.08, thetalbl, /normal, size=1.4, color= 0, charthick= 3.0

philbl= strcompress(string(rotate_phi),/remove_all)
digs= 3
if rotate_phi gt 100.0 then digs= 5
if rotate_phi gt 10.0 then digs= 4
if rotate_phi lt 0.0 then digs= 4
if rotate_phi lt 10.0 then digs= 5
if rotate_phi lt 100.0 then digs= 6
philbl= '!7u!6='+strmid(philbl,0,digs)
xyouts, x0+0.02, y1-0.14, philbl, /normal, size=1.4, color= 0, charthick= 3.0






;  right - velocity field
; ----------------------------
x_orig= x
y_orig= y
z_orig= z
vx_orig= vx
vy_orig= vy
vz_orig= vz
center= orig_center

set_maxden= 150.0
set_dynrng= 2.0*set_maxden

process_kin_velimage_and_overplot, x_orig, y_orig, z_orig, vx_orig, vy_orig, vz_orig, mass, $
	xlen, x2, y0, x3, y1, /velocitymap, /dispersionmap, $   ; both are set - will do both fits, but only print velmap
	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	center=center, $ ;/showbhs, $
	set_maxden=set_maxden, set_dynrng=set_dynrng, pixels=pixels, $
	;pa=pa, /show_maj_slit, /show_min_slit, $
	pa=pa, /show_maj_slit, $
	xslit=slit_len, yslit=0.5*slit_width, $
	fitstoo=fitstoo, filename=filename, $
	msg1=' ', msg2=' ', $ ; /drawbox, xslit=slit_len, yslit=0.5*slit_width, $
	ContourMap=ReturnedMap, HalfMassSB=HalfMassSB, /colorbar_onright








; any extras?
; -------------

;xyouts, 0.02, 0.45, fload_timelbl(1,2,/noteq), /normal, size= 1.33, charthick=3.0, color= 0


device, /close

; -------------
;  Done
; -------------


end



; -----------------------------------------------------------------------------



