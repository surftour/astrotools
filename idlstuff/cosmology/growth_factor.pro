function growth_factor, z, BRUTE_FORCE=BRUTE_FORCE
	COMMON COSMOLOGICAL_PARAMETERS
	;; growth factor D(z) normalized to D(z=0)=1
	forward_function growth_factor_deriv
	
	if (keyword_set(BRUTE_FORCE)) then begin
	D_0 = 1.0
	a = 1./(1.+z)
	prefac = SQRT(Omega_LAMBDA*a*a*a + Omega_MATTER)/(a^(3./2.))
	intfac = QROMO('growth_factor_deriv',0.,a)
	a0 = 1.
	prefac_0 = SQRT(Omega_LAMBDA*a0*a0*a0 + Omega_MATTER)/(a0^(3./2.))
	intfac_0 = QROMO('growth_factor_deriv',0.,a0)
	return, (prefac*intfac)/(prefac_0*intfac_0)	
	endif

;; alternatively, could adopt the 
;; use the Carroll et al. 1992 approximation, which is basically the same to ~1%
;;
;	if (keyword_set(BRUTE_FORCE)) then begin
	Omz = omega_matter_z(z)
	Olz = omega_lambda_z(z)
	gz  = (5./2.)*omz/(omz^(4./7.) - olz + (1.+omz/2.)*(1.+olz/70.))
	Omz = OMEGA_MATTER
	Olz = OMEGA_LAMBDA
	g0  = (5./2.)*omz/(omz^(4./7.) - olz + (1.+omz/2.)*(1.+olz/70.))
	Dz  = (gz/g0)/(1.+z)
	return, Dz
;	endif
end
