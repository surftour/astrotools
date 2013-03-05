;--------------------------------------------------------------
; Calculate BH luminosity/magnitude in different bands, given
;   an bolometric luminosity (in solar luminosities), 
;   based on the Marconi et al. (2004) bolometric corrections
;--------------------------------------------------------------
pro marconi_qso_band_lum, bolometric_lum, L_softX, L_hardX, L_B

  L0 = alog10(bolometric_lum)	; Log(L) for magnitudes, conversions, etc.
  L1 = L0 - 12.0 				; (script L of Marconi et al. 2004)

  ; Marconi et al. (2004) relations for bolometric & x-ray luminosity
  logL_soft_xray = L0 - (1.65 + 0.22*L1 + 0.012*L1*L1 - 0.0015*L1*L1*L1)	; 0.5-2 keV
  logL_hard_xray = L0 - (1.54 + 0.24*L1 + 0.012*L1*L1 - 0.0015*L1*L1*L1)	; 2-10 keV
  log_nuB_LB     = L0 - (0.80 - 0.067*L1 + 0.017*L1*L1 - 0.0023*L1*L1*L1)	; B-band
  nuB = 6.818d14		; assumes B-band at 0.44 microns
  L_B = 10^log_nuB_LB	; gives nu_B * L_nu_B


  ; Marconi alpha_OX (Zamorani 82, emperical Vignali, Brandt, & Scheider 03)
  		R = 2.605	; [defined as -log(nu(2500 Ang) / nu(2 keV)) ]
  logL_UV = (1.85*R + logL_soft_xray)						; 2500 Angstroms

  L_hardX = 10^(logL_hard_xray)
  L_softX = 10^(logL_soft_xray)

  return
end
