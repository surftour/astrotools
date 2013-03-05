function func_gaussian_pluszeropt, x, p

	; p[0] = Normalization, i.e., the SFR_max
	; p[1] = Width of the dist.; sigma
	; p[2] = mean of the dist., i.e., time of SFR_max
	; p[3] = zero point, i.e. allow it to be non-zero
	return, p[0]*exp(-(x-p[2])^2/(2*(p[1]^2))) + p[3]

end


