; ---------------------------------------------------------------------------
;
;  Determine sigma (and Reff) for N random lines of sight (used below)
;
; ---------------------------------------------------------------------------
function process_2d_sigma_avg, x, y, z, vx, vy, vz, mass, reff, $
			avg_Reff=avg_Reff, center=center, do_quarter_mass=do_quarter_mass, $
			fixed_mass=fixed_mass, $
			usefixedRe=usefixedRe, sigmaerr=sigmaerr, refferr=refferr


n_projections= 12 
;n_projections= 100
;n_projections= 144
;n_projections= 1000
seed= 104L

allthetas= fltarr(n_projections+1)
allphis= fltarr(n_projections+1)
allsigs= fltarr(n_projections+1)
allReffs= fltarr(n_projections+1)


if keyword_set(center) then c=center else c=[0,0,0]


for i=0,n_projections do begin

        rdphi= randomu(seed)
        rdtheta= randomu(seed)

        ;allthetas[i]= rdtheta*!PI
        allthetas[i]= acos(rdtheta)   ; the previous line isn't random
        allphis[i]= rdphi*2*!PI

        ; in radians
        theta= allthetas[i]
        phi= allphis[i]


	; rotate
        rot_x= (x-c[0])*(cos(theta)*cos(phi)) + (y-c[1])*(cos(theta)*sin(phi)) + (z-c[2])*sin(theta)
        rot_y= -(x-c[0])*sin(theta) + (y-c[1])*cos(theta)
        ;rot_z= x*(sin(theta)*cos(phi)) - y*(sin(theta)*sin(phi)) + z*cos(theta)
        ;rot_vx= vx*(cos(theta)*cos(phi)) + vy*(cos(theta)*sin(phi)) + vz*sin(theta)
        ;rot_vy= -vx*sin(theta) + vy*cos(theta)
        rot_vz= vx*(sin(theta)*cos(phi)) - vy*(sin(theta)*sin(phi)) + vz*cos(theta)

	; center should already be taken into consideration
	reff_avg= -1.0
	if keyword_set(do_quarter_mass) then reff_avg= process_2d_halfmass(rot_x, rot_y, mass, center=[0,0,0], /do_quarter_mass)
	if keyword_set(fixed_mass) then reff_avg= process_2d_halfmass(rot_x, rot_y, mass, center=[0,0,0], fixed_mass=fixed_mass)
	if reff_avg lt 0 then reff_avg= process_2d_halfmass(rot_x, rot_y, mass, center=[0,0,0])

	re_touse= reff_avg
	if keyword_set(usefixedRe) then re_touse= reff
        sigma= process_2d_sigma(rot_x, rot_y, rot_vz, re_touse)

	; save
        allsigs[i]= sigma
        allReffs[i]= reff_avg

        
        if (i mod 100) eq 0 then print, "i= ",i

endfor  
        

sigma_moment= moment(allsigs)
sigma_random= sigma_moment[0]
sigmaerr= sqrt(sigma_moment[1])

print, "average velocity dispersion= ", sigma_random, "  (+/- ",sqrt(sigma_moment[1]),")"

reff_moment= moment(allReffs)
avg_Reff= reff_moment[0]
refferr= sqrt(reff_moment[1])

print, "average effective radius=    ", avg_Reff, "  (+/- ",sqrt(reff_moment[1]),")"

return, sigma_random

end









; ---------------------------------------------------------------------------

