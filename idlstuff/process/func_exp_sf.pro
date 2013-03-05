function func_exp_sf, x, p

	; 
	; p[0] = Norm
	; p[1] = tau
	; p[2] = + constant
	;
	return, p[0] * exp(-x/p[1]) + p[2]

end


