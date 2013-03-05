function func_exponential, x, p

	; 
	; p[0] = Mass
	; p[1] = R_d
	;
	Sigma_0= p[0] / (2. * !PI * p[1] * p[1])
	return, Sigma_0 * exp(-x/p[1])

end


