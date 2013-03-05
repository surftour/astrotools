function func_linear, x, p

	; 
	; p[0] = a, the intercept
	; p[1] = m, the slope
	;
	; a + m*x
	;
	return, p[0] + (x * p[1])

end


