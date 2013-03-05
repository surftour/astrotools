;; use the Mandelbaum et al. recent stellar mass-halo mass correlation measurements
;;   to convert from stellar mass in Msun to halo mass in (Msun/h) 
;;
;; default to their early-type measurements
;;
function mhalo_of_mstellar, log_Mstellar, MEAN=MEAN, $
	LATE_TYPE=LATE_TYPE, EARLY_TYPE=EARLY_TYPE, ROBUST=ROBUST
	

	; default to early type 
	if ((keyword_set(EARLY_TYPE) EQ 0) AND (keyword_set(LATE_TYPE) EQ 0)) $
		then EARLY_TYPE=1
		
	if (keyword_set(EARLY_TYPE)) then begin
		Mstellar_obs = [ 0.76, 1.50, 3.00, 5.80, 11.20, 21.30, 39.60 ] * 1.0d10
		Mhalo_obs    = [ 3.18, 4.20, 4.90, 14.1, 34.00, 158.0, 716.0 ] * 1.0d11 ;/ h
		dMhalo_obs_p = [ 9.45, 6.63, 4.70,  5.6,  10.0,  37.0, 123.0 ] * 1.0d11 ;/ h
		dMhalo_obs_m = [ 3.14, 3.67, 3.20,  5.3,  9.00,  33.0, 190.0 ] * 1.0d11 ;/ h
	endif

	if (keyword_set(LATE_TYPE)) then begin
		Mstellar_obs = [ 0.73, 1.50, 2.90, 5.60, 10.80, 20.50, 40.40 ] * 1.0d10
		Mhalo_obs    = [ 0.02, 6.60, 6.10, 14.0, 13.00, 34.00, 180.0 ] * 1.0d11 ;/ h
		dMhalo_obs_p = [ 1.56, 5.10, 4.70,  8.0,  12.0,  33.0, 532.0 ] * 1.0d11 ;/ h
		dMhalo_obs_m = [0.018, 4.00, 3.40,  7.0,  9.00,  28.0, 173.0 ] * 1.0d11 ;/ h
	endif

	Mstar = alog10(Mstellar_obs)
	Mhalo = alog10(Mhalo_obs)
		M = INTERPOL(Mhalo,Mstar,log_Mstellar)
		
	if (keyword_set(ROBUST)) then $
		if (log_Mstellar LT Mstar[3]) then $
			M = log_Mstellar - (Mstar[3]-Mhalo[3])
		;;
		;; applies a flat efficiency below ~10^11 Msun, chosen at the lowest mass 
		;;  where the efficiencies are well-determined
		;;
		
	if (keyword_set(MEAN)) then M = log_Mstellar - alog10((0.046/0.27)*0.17/0.7)
		;; apply their mean efficiency (w.r.t. the baryon fraction) of 0.17

	return, M
end
