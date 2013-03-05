function minor_func_samsblp, x, p

	;x= is the mass ratio

	y= p[0]*alog10(x) + p[1]

	return, y

end
