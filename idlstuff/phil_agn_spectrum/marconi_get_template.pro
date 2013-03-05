;; call to load template marconi et al. spectrum
function marconi_get_template, Lbol
	COMMON TEMPLATE_SPECTRA_BLOCK
	flot_index = (alog10(Lbol)-7.0)/(0.1)
	Lbol_index = fix(flot_index)
	if (Lbol_index LT 0) then begin
		print, 'Bolometric luminosity is too low - using minimum (10^7 Lsolar) spectrum'
		Lbol_index = 0
		return, template_nuLnu_list[0,*]
	end
	if (Lbol_index GT template_N_Lbol-2) then begin
		print, 'Bolometric luminosity is too high - using maximum (10^17 Lsolar) spectrum'
		Lbol_index = template_N_Lbol-1
		return, template_nuLnu_list[template_N_Lbol-1,*]
	end
	dy = template_nuLnu_list[Lbol_index+1,*] - template_nuLnu_list[Lbol_index,*]
	dx = flot_index - Lbol_index
	return, (template_nuLnu_list[Lbol_index,*] + dx*dy)
end

