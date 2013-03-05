function func_hernquist_mass, x, p

	; 
	; p[0] = Total Mass
	; p[1] = scale radius, a
	;
	return, p[0] * x * x / ((x + p[1])*(x + p[1]))

end


