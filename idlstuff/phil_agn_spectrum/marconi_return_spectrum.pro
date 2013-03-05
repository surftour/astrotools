;; function to return the marconi spectrum in log(nuLnu), in solar 
;;;   luminosities, as a function of input frequency (in Hz) and 
;;;	  bolometric luminosity (Lbol, in solar)
;;;

COMMON TEMPLATE_SPECTRA_BLOCK,  template_freq_list, template_nuLnu_list, $
			template_Lbol_list, template_N_Lbol, template_N_freq

function marconi_return_spectrum, freq_list, Lbol
	COMMON TEMPLATE_SPECTRA_BLOCK
	forward_function marconi_get_template
	marconi_load_template
	L = marconi_get_template(Lbol)
	return, INTERPOL(L,template_freq_list,alog10(freq_list))
end

