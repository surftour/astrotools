;; load template marconi et al. agn spectra
pro marconi_load_template
	COMMON TEMPLATE_SPECTRA_BLOCK
  	OPENR,1,return_idl_routines_homedir(0)+'/agn_spectrum/marconi_template_spectra.dat'
  	template_N_Lbol = 0
  	template_N_freq = 0
  	READU,1,template_N_freq,template_N_Lbol
  	template_freq_list = fltarr(template_N_freq)
  	READU,1,template_freq_list,TRANSFER_COUNT=template_N_freq
  	template_nuLnu_list = fltarr(template_N_Lbol,template_N_freq)
  	READU,1,template_nuLnu_list
  	CLOSE,1
  	template_Lbol_list = 7.0+0.1*findgen(template_N_Lbol)
end
