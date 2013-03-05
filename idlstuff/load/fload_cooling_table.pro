;----------------------------------------------------------
;   Load the Cooling table from 
;----------------------------------------------------------
function fload_cooling_table, Metallicity


COMMON TabulatedCoolingData, coolingtemp, coolingLambdaN

	if not keyword_set(Metallicity) then Metallicity= 'zero'
	;Metallicity= '-00'


	; get cooling data
	; ---------------------
	coolingfile= '/home/tcox/Tools/CoolingFiles/m'+Metallicity+'.cie.noheader'
	coolingdata= read_ascii(coolingfile)

	coolingtemp= coolingdata.field01[0,*]
	coolingLambdaN= coolingdata.field01[4,*]

	print, "Loaded Data from", coolingfile

	return, Metallicity
end



