;; returns  (k^3/(2pi^2))*P(k,z=0)  -- in CDM, no post-recombination 
;;   changes to the linear power spectrum, so make the redshift-dependent 
;;   version by just taking this * (D1(z)/D1(z=0))^2

function power_spectrum, k
	COMMON COSMOLOGICAL_PARAMETERS
	forward_function master_transfer_function
	
	T = master_transfer_function(k)
	;; k will be in 1/Mpc, H0 = h * 100 km/s/Mpc, so c = 3x10^5 km/s
	H0 = 100. * HUBBLE
	c = 2.998d5
	delta_H = 1.0	;; amplitude of perturbations on the horizon scale today -- 
					;;   just leave it at 1, since the normalization we care about
					;;   will be the ultimate sigma norm, normalized with sigma_8
	n_squig = SPECTRAL_INDEX - 1.
	delta_H = 1.94d-5 * (OMEGA_MATTER^(-0.785 - 0.05*alog(OMEGA_MATTER))) * $
				EXP(-0.95*n_squig - 0.169*n_squig*n_squig)				
					
	scale_fac = (c*k/H0)^(3.+SPECTRAL_INDEX)
	return, delta_H*delta_H*scale_fac*T*T	
end
