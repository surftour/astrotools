function func_hernquist, x, p

	; 
	; p[0] = Total Mass
	; p[1] = scale radius, a
	;
	return, (p[0] / (2.0*!PI)) * (p[1] / x) / ((x + p[1])*(x + p[1])*(x + p[1]))

end


