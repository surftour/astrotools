function minor_func_samsbl2, x, p

	;x= is the mass ratio

	y= x*0.0

	idx= where(x gt p[2])
	if idx(0) ne -1 then begin
		y(idx)= p[1]*((x(idx) - p[2])^(p[0]))
	endif

	return, y

end
