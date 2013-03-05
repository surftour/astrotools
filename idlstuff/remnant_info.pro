pro remnant_info, frun


;frun= "data/Sbc201a-u4"
frun= "data/"+frun
;frun= "data/Sbc201b"
;frun= "data/Sbc201b_wBH"

snapnum= 30
;snapnum= 25

;ok=fload_snapshot(frun,snapnum)
ok=fload_snapshot_bh(frun,snapnum)


ngas= fload_npart(0)


; ----------------------------------------
; This script calls all our procedures
; which calculate and display remnant
; information
; ----------------------------------------

        
initialize_plotinfo, 1
 


plot_devacouleur_profile, frun, sendto, smoothlen=smoothlen, filename=frun+"/mrem.eps"

plot_density_profile_dm, frun, sendto, smoothlen=smoothlen, filename=frun+"/dmprofile.eps"

plot_density_profile_b, frun, sendto, smoothlen=smoothlen, filename=frun+"/bprofile.eps"

plot_density_profile_ratio, frun, sendto, smoothlen=smoothlen, filename=frun+"/proratio.eps"

if ngas gt 0 then begin
   plot_density_profile_gas, frun, sendto, smoothlen=smoothlen, filename=frun+"/gprofile.eps"
endif



plot_Re_sigma_random, frun, sendto


plot_surface_density, frun, sendto, filename=frun+"/surden.eps"


plot_circular_velocity, frun, sendto, filename=frun+"/vc.eps"


plot_sigma_center, frun, sendto, filename=frun+"/csigma.eps"


; lay slits across merger remnant 
plot_sigma_slit, /xy, frun, sendto, $
	slitfilename= frun+"/slitxy-baryon.eps", $
	sigfilename= frun+"/slitxy-sigma.eps", $
	velfilename= frun+"/slitxy-velocity.eps", $
	vsigfilename= frun+"/slitxy-vsig.eps"

plot_sigma_slit, /xz, frun, sendto, $
        slitfilename= frun+"/slitxz-baryon.eps", $
        sigfilename= frun+"/slitxz-sigma.eps", $
        velfilename= frun+"/slitxz-velocity.eps", $
	vsigfilename= frun+"/slitxz-vsig.eps"

plot_sigma_slit, /yz, frun, sendto, $
        slitfilename= frun+"/slityz-baryon.eps", $
        sigfilename= frun+"/slityz-sigma.eps", $
        velfilename= frun+"/slityz-velocity.eps", $
	vsigfilename= frun+"/slityz-vsig.eps"


if ngas gt 0 then begin
   plot_radialvelocity_vs_r, frun, sendto, filename=frun+"/vel-r.eps"
endif


if ngas gt 0 then begin
   plot_gascoolingtime_histogram, frun, sendto, filename=frun+"/gcoolt.eps"
endif




print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end














