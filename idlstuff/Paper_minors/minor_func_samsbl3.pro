function minor_func_samsbl3, x, p

	;x= is the mass ratio
	xcutoff= 0.09
	;xcutoff= 0.07

	y= x*0.0

	idx= where(x gt xcutoff)
	if idx(0) ne -1 then begin
		y(idx)= p[1]*((x(idx) - xcutoff)^(p[0]))
	endif

	return, y

end
