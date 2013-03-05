;--------------------------------
;
;  This procedure converts between
; luminosity density and surface brightness,
; in this direction.
;
; This is problem 2.2 in Binney &
;    Merrifield.
;
;
;  band= observaed band, U,B,V, ...
;  lum_density= in L_solar/pc^2
;--------------------------------
pro astro_surfacebrightness, band, lum_density


	solarmag= astro_solarM(band)

	onearcsec_at_10pc= 20556.

	; -2.5 log (onearcsec_at_10pc ^2)

	SB_band_norm= solarmag + 5.0 * alog10(onearcsec_at_10pc)

	SB = SB_band_norm - 2.5 * alog10(lum_density)

	return, SB



end



