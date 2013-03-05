function func_gaussian, x, p

	; p[0] = Normalization, i.e., the SFR_max
	; p[1] = Width of the dist.; sigma
	; p[2] = mean of the dist., i.e., time of SFR_max
	return, p[0]*exp(-(x-p[2])^2/(2*(p[1]^2)))

end


