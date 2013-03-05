pro test_sig
	COMMON COSMOLOGICAL_PARAMETERS

	M0 = 8.+0.25*findgen(33)
	s0 = 0.*M0
	for i=0,n_elements(M0)-1 do s0[i]=sigma_of_Mhalo_z0(M0[i])      
	plot,M0,s0,/ylog,xstyle=1,xrange=[MIN(M0),MAX(M0)]                                  
		m00 = alog10((2.78d11 * OMEGA_MATTER)*4.*!PI/3. * (8.)^3)  - alog10(0.7)

		oplot,[2.,20.],[SIGMA_8,SIGMA_8],linestyle=1            
		oplot,[m00,m00],[0.00001,100.],linestyle=1
end

pro compare_dift_croom
	h = 0.7
	z = [0.526, 0.804, 1.026, 1.225, 1.413, 1.579, 1.745, 1.921, 2.131, 2.475]
	s = [ 5.73,  3.94,  4.76,  5.52,  5.28,  4.87,  6.25,  6.39,  8.00,  8.81]
	ds= [ 0.79,  1.00,  0.97,  0.98,  0.98,  0.95,  0.83,  0.98,  0.99,  0.98]
	b = [ 1.13,  1.49,  1.71,  2.31,  2.32,  2.24,  2.17,  2.91,  3.53,  4.24]
	db= [ 0.18,  0.21,  0.24,  0.23,  0.27,  0.30,  0.35,  0.35,  0.38,  0.53]
	M = [ 0.82,  2.09,  2.31,  5.76,  3.69,  2.05,  1.15,  3.05,  4.46,  4.78] * 1.0d12
	Mp= [ 1.55,  2.18,  2.23,  2.90,  2.24,  1.61,  1.24,  1.85,  2.20,  2.68] * 1.0d12
	Mm= [ 0.47,  1.30,  1.37,  2.21,  1.62,  1.07,  0.72,  1.34,  1.68,  1.99] * 1.0d12

	plot,z,s,/nodata,xrange=[0.,5.],yrange=[0.,15.],xstyle=1,ystyle=1
		
	;; their measured s0
	plotsym,8,/fill
	oploterror,z,s,ds,PSYM=8
	
	;; direct from bias 
	r0 = 5. * (b*growth_factor(z))^(2./1.8)
	plotsym,0
	oploterror,z,r0,ds,PSYM=8

	
	;; compare bias
	f  = load_cosmological_parameters(1)
	om = omega_matter_z(z)
	xi_rat = (s/4.)^(1.8)/(growth_factor(z)^2)
	b_est  = SQRT(xi_rat - 4.*(om^1.2)/45.) - (om^0.6)/3.
	
	xi_rat_real = b_est^2
	r_est = (xi_rat_real*(growth_factor(z)^2))^(1./1.8) * 4.
	plotsym,3,1.2,/fill
	oplot,z+0.03,r_est,PSYM=8
	
	r_est = (((b + (om^(0.6)/3.))^2 + (4./45.)*(om^(1.2)))*(growth_factor(z)^2))^(1./1.8) * 4.
	plotsym,3,1.2
	oplot,z-0.03,r_est,PSYM=8

	b0 = 0.*b
	z00 = [0., 0.5, 1., 2., 3., 4., 5.]
	om  = omega_matter_z(z00)
	for i=0,n_elements(z00)-1 do b0[i] = bias_of_mhalo(12.+2.*alog10(1.+z00[i]),z00[i])
	r_est = (((b0 + (om^(0.6)/3.))^2 + (4./45.)*(om^(1.2)))*(growth_factor(z00)^2))^(1./1.8) * 4.
	plotsym,3,2.
	oplot,z00-0.03,r_est,PSYM=-7,SYMSIZE=2.

	for i=0,n_elements(z00)-1 do b0[i] = bias_of_mhalo(12.+alog10(3.),z00[i])
	r_est = (((b0 + (om^(0.6)/3.))^2 + (4./45.)*(om^(1.2)))*(growth_factor(z00)^2))^(1./1.8) * 4.
	plotsym,0
	oplot,z00-0.03,r_est,PSYM=-8,SYMSIZE=2.

	for i=0,n_elements(z00)-1 do b0[i] = bias_of_mhalo(14.+alog10(0.7),z00[i])
	r_est = (((b0 + (om^(0.6)/3.))^2 + (4./45.)*(om^(1.2)))*(growth_factor(z00)^2))^(1./1.8) * 4.
	plotsym,8
	oplot,z00-0.03,r_est,PSYM=-8,SYMSIZE=2.


	;plot,z,b
	;	oploterror,z,b,db,PSYM=8
	;	oplot,z,b_est,PSYM=7

end

