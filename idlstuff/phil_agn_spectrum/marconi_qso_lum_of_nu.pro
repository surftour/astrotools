;--------------------------------------------------------------------------
; Function to calculate the luminosity (L_solar per Hz) 
;   at a given frequency f [Hz], from a given bolometric luminosity
;   (Complete broken-power-law model intrinsic QSO spectrum up to 10keV)
;--------------------------------------------------------------------------
function marconi_qso_lum_of_nu, L_bol, f
  L0 = alog10(L_bol)			; Log(L) for magnitudes, conversions, etc.
  L1 = L0 - 12.0 				; (script L of Marconi et al. 2004)

  ; Marconi et al. (2004) relations for bolometric & x-ray luminosity
  logL_soft_xray = L0 - (1.65 + 0.22*L1 + 0.012*L1*L1 - 0.0015*L1*L1*L1)	; 0.5-2 keV
  logL_hard_xray = L0 - (1.54 + 0.24*L1 + 0.012*L1*L1 - 0.0015*L1*L1*L1)	;  2-10 keV

  ; Now assume a photon index (photons/keV ~ E^(-Gamma)) Gamma=1.9
  ;		(Marconi refs. George et al. 1998, Perola et al. 2002, Nandra & Pounds 1994)
  Gamma = 1.9
    integration_factor = (1.0/(2.0-Gamma)) * ((10.0/2.0)^(2.0-Gamma) - 1.0)
    L_hard_xray = 10.0^logL_hard_xray	; 2-10 keV band luminosity in L_solar
    f_2keV = 2.0 * (2.41775d17)			; frequency of 2keV in Hz
  dL_df_2keV = L_hard_xray / (integration_factor * f_2keV)

  ;print,logL_hard_xray,L_hard_xray

  ; Now can calculate other xray luminosities up to 10keV (at larger energies, the 
  ;   Compton reflection hump comes in, needs to be modeled)
  f_05keV = 0.25 * f_2keV		; frequency [Hz] at 0.5 keV
  f_10keV = 5.0 * f_2keV		; frequency [Hz] at 10 keV
  if ((f_05keV LE f) AND (f_10keV GE f)) then DL_DF = dL_df_2keV * ((f/f_2keV)^(1 - Gamma))
  if (f GT f_10keV) then DL_DF = 0.0  	; Sorry, not incorporated yet! (need Compton modeling)
  dL_df_05keV = dL_df_2keV * ((f_05keV/f_2keV)^(1 - Gamma))
  
  
  ; Using Marconi alpha_OX (Zamorani 82, emperical Vignali, Brandt, & Scheider 03)
  ;   which relates flux at 2500 Angstroms to that at 2keV  by:
  ;    a_OX = - LOG(L_2500 / L_2keV) / LOG(f_2500 / f_2keV)   and
  ;    a_OX = - 0.11*LOG(L_2500) + 1.85
  ;  (from Vignali, Brandt, & Schneider 2003, 
  ;     which needs everything in units of ergs/s/Hz
  R = -2.605	; [defined as log(nu(2500 Ang) / nu(2 keV)) ]
  logL_solar_ergs = 33.591 		; Log of solar luminosity in ergs/s
  log_dLdf2kev_ergs = alog10(dL_df_2kev) + logL_solar_ergs
  logdL_df_2500_ergs = (1.85*R + log_dLdf2kev_ergs)/(1.0 + 0.11*R)
  logdL_df_2500 = logdL_df_2500_ergs - logL_solar_ergs
  dL_df_2500 = 10.0^logdL_df_2500
  
  
  ; Now use the SDSS, etc. average optical power law L ~ f^(-0.44)
  ;   from 1micrometer to 1300 Angstroms  (Vanden Berk et al. 2001)
  f_micrometer = 2.998d14		; frequency in [Hz] at 1 micrometer
  f_2500ang    = 1.199d15		; frequency in [Hz] at 2500 Angstroms
  f_1300ang    = 2.306d15		; frequency in [Hz] at 1300 Angstroms
  alpha = -0.44					; spectral index in this range
  if ((f_micrometer LE f) AND (f_1300ang GE f)) then DL_DF = dL_df_2500 * ((f/f_2500ang)^(alpha))
  dL_df_micrometer = dL_df_2500 * ((f_micrometer/f_2500ang)^(alpha))
  dL_df_1300 = dL_df_2500 * ((f_1300ang/f_2500ang)^(alpha))
  
  
  ; For wavelengths larger than a micrometer, assume the spectrum cuts off with 
  ;   alpha = 2  (as in Rayleigh-Jeans tail of blackbody)
  alpha = 2.0
  if (f LT f_micrometer) then DL_DF = dL_df_micrometer * ((f/f_micrometer)^(alpha))

  
  ; Assume flat over small break 1300-1200 Angstroms
  f_1300ang    = 2.306d15		; frequency in [Hz] at 1300 Angstroms
  f_1200ang    = 2.498d15		; frequency in [Hz] at 1200 Angstroms
  if ((f_1300ang LE f) AND (f_1200ang GE f)) then DL_DF = dL_df_1300
  dL_df_1200 = dL_df_1300
  
  
  ; From 1200 - 500 Angstroms we have alpha = -1.76 
  ; 	(Telfer et al. 2002, Vanden Berk et al. 2001)
  f_500ang    = 5.996d15		; frequency in [Hz] at 500 Angstroms
  alpha = -1.76
  if ((f_1200ang LE f) AND (f_500ang GE f)) then DL_DF = dL_df_1200 * ((f/f_1200ang)^(alpha))
  dL_df_500 = dL_df_1200 * ((f_500ang/f_1200ang)^(alpha))
  
  
  ; Now, finally, we use the values of dL_df at 500 Angstroms and 0.5 keV
  ;   to create a power law connecting the two frequencies and project over it
  alpha = alog10(dL_df_05keV/dL_df_500) / alog10(f_05keV/f_500ang)
  if ((f_500ang LE f) AND (f_05keV GE f)) then DL_DF = dL_df_500 * ((f/f_500ang)^(alpha))

  ;print, alog10(f_micrometer*dL_df_micrometer)
  ;print, alog10(f_2500ang*dL_df_2500)
  ;print, alog10(f_1300ang*dL_df_1300)
  ;print, alog10(f_1200ang*dL_df_1200)
  ;print, alog10(f_500ang*dL_df_500)
  ;print, alog10(f_05keV*dL_df_05keV)
  ;print, alog10(f_2keV*dL_df_2keV)
  ;checked out, matches well with Marconi et al. (2004) model spectrum

  return, DL_DF
end
