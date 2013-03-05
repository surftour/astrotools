function func_devac, x, p

	; note, that we are in real space
	; p[0] = I (R_e)
	; p[1] = R_e
	power= -3.331 * ((x/p[1])^(0.25) - 1.0)
	return, p[0]*10^(power)

end


