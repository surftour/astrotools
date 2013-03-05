function func_nfw, x, p

  COMMON NFWParams, rho_0, r_s

	rho_0= p[0]
	r_s= p[1]

	;if p[0] lt 0.0 then p[0]= -p[0]*100
	;if p[1] lt 0.0 then p[1]= -p[1]

	return, p[0]*((x/p[1])^(-1))*((1+(x/p[1]))^(-2))

end
