;----------------------------------------------------------
;   Actually solve for the Cooling rate
;----------------------------------------------------------
function fload_lambdan_array, temp

COMMON TabulatedCoolingData


	idx= where(temp eq 0)
	if idx(0) EQ -1 then LogTemp= alog10(temp) else return, 0

	n= n_elements(temp)
	i= intarr(n)
	LambdaN= fltarr(n)


        ; get index of cooling table
        ; -----------------------------
        i= fix((LogTemp-coolingtemp[0])/0.05)

	; trim
	; -----
	idx_low= where(LogTemp LT coolingtemp[0])
	idx_high= where(LogTemp GE coolingtemp[90])

	if idx_low(0) ne -1 then i(idx_low)= 0
	if idx_high(0) ne -1 then i(idx_high)= 89

	; now actually map interpolate
	; -----------------------------
	f= (LogTemp - coolingtemp[i])/0.05
	if idx_low(0) ne -1 then f(idx_low)= 0

	LambdaN= coolingLambdaN[i]*(1-f) + coolingLambdaN[i+1]*f

	LambdaN= 10^(LambdaN)

	if idx_low(0) ne -1 then LambdaN(idx_low)= 0

	return, LambdaN

end


