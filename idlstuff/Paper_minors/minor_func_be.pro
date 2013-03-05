function minor_func_be, x, p

	;x= is the mass ratio
	return, p[0]*( x / ( 1 + x))

end
