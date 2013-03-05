function func_c, x

   COMMON NFWParams

	;rho_0= 0.00105243/2.77550e-08
	if rho_0 lt 1.0 then rho0= rho_0/2.77550e-08 else rho0= rho_0


	;x= is the guess at the concentration
	return, x[0]*x[0]*x[0]/(alog(1.0+x[0]) - (x[0]/(1.0+x[0]))) - (3.0*rho0/105.0)

end
