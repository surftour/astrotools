;; transfer functions from EH99 used to compute power spectra
;;
;; for now, we ignore massive neutrinos, since only upper limits to 
;;   Omega_nu exist -- in this case, Dcb = Dcbnu = D1, and 
;;   Tcb(q,z) = Tmaster(q) (post-recomination transfer function)
;;
function master_transfer_function, k
	;; k must be in inverse Mpc !!!
	COMMON COSMOLOGICAL_PARAMETERS
	Omega_nu = 0.
	Omega_M  = OMEGA_MATTER
	Omega_b  = OMEGA_BARYON
	theta_cmb= THETA_CMB
	h = HUBBLE
	
	
	f_nu = Omega_nu/Omega_M
	f_nub = Omega_b/Omega_M
	N_nu = 3.
	f_c = (Omega_M-Omega_b-Omega_nu)/Omega_M
	f_cb = (Omega_M-Omega_nu)/Omega_M
	zeq = 2.50d4 * Omega_M *h*h* (theta_cmb^(-4))
	b1 = 0.313 * ((Omega_M*h*h)^(-0.419)) * (1. + 0.607*((Omega_M*h*h)^(0.674)))
	b2 = 0.238 * ((Omega_M*h*h)^(0.223))	
	z_d = 1291. * (((Omega_M*h*h)^(0.251))/(1. + 0.659*((Omega_M*h*h)^(0.828)))) * $
					(1. + b1*((Omega_b*h*h)^b2))
	y_d = (1. + zeq)/(1. + z_d)
	p_c = (1./4.)*(5. - SQRT(1. + 24.*f_c))
	p_cb= (1./4.)*(5. - SQRT(1. + 24.*f_cb))

	alpha_nu = (f_c/f_cb) * ((5. - 2.*(p_c + p_cb))/(5. - 4.*p_cb)) * $
		((1. - 0.553*f_nub + 0.126*(f_nub^3))/(1. - 0.193*SQRT(f_nu*N_nu) + 0.169*f_nu*(N_nu^(0.2)))) * $
		((1. + y_d)^(p_cb - p_c)) * $
		(1. + ((p_c - p_cb)/2.)*(1. + 1./((3. - 4.*p_c)*(7. - 4.*p_cb)))/(1. + y_d))
	
	s = 44.5 * alog(9.83/(Omega_M*h*h)) / SQRT(1. + 10.*((Omega_b*h*h)^(3./4.))) ;; dimensional, Mpc
	Gamma_eff = Omega_M*h*h*(SQRT(alpha_nu) + (1.-SQRT(alpha_nu))/(1.+(0.43*k*s)^4))
	q_eff = k * theta_cmb*theta_cmb/Gamma_eff	;; dimensional! k must be in Mpc^-1 !
	beta_c = 1./(1. - 0.949*f_nub)
	L = alog(exp(1.) + 1.84*beta_c*SQRT(alpha_nu)*q_eff)
	C = 14.4 + 325./(1. + 60.5 * (q_eff^(1.11)))

	T_sup = L / (L + C*q_eff*q_eff)
	B_k = 1.	;; again, ignoring the massive neutrino corrections for now
	return, T_sup*B_k
end
