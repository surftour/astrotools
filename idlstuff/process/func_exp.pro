function func_exp, x, p

	; 
	; p[0] = Norm
	; p[1] = tau
	;
	return, p[0] * exp(-x/p[1])

end


