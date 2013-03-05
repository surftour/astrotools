;; called by sigma_of_Mhalo_z0.pro
function sigma_of_Mhalo_z0_deriv, ln_k
	COMMON COSMOLOGICAL_PARAMETERS
	COMMON sigma_of_Mhalo_z0_TO_PASS

	k = EXP(ln_k)
	x = R_of_M_to_pass * k	;; R must be in Mpc
	window_part = (abs(real_space_window(x)))^2
	return, power_spectrum(k) * window_part
end

