function func_sersic, x, p

	; note, that we are in real space
	; p[0] = I (R_e)
	; p[1] = R_e
	; p[2] = n
	;
	; b = dimensionless constant fixed
	;     such that p[1] is R_e, i.e.,
	;     the projected radius enclosing
	;     half of the total luminosity
	;
	b = 2 * p[2] - 0.33 + 4.0/(405.*p[2]) + 46.0/(25515*p[2]*p[2])
	power= -1.0 * alog10(exp(1)) * b * ((x/p[1])^(1.0/p[2]) - 1.0)
	return, p[0]*10^(power)

end


