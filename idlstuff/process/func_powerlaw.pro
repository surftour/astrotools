function func_powerlaw, x, p

	; 
	; p[0] = Normalization
	; p[1] = power
	;
	return, p[0] * (x^p[1])

end


