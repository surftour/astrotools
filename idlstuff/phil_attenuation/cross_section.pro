;--------------------------------------------------------------------------
; Function to calculate the cross section per H atom at a given 
;    frequency f (in Hz)
;--------------------------------------------------------------------------
function cross_section, f0, METALLICITY_OVER_SOLAR=METALLICITY_OVER_SOLAR

  SIGMA = 0.0
  sigma_temp = 0.0
  if (keyword_set(METALLICITY_OVER_SOLAR) EQ 0) then METALLICITY_OVER_SOLAR=1.0

  ; For 0.03 keV < E < 10 keV  
  ;   (7.2e15 < nu[Hz] < 2.4e18  or   1.2 < lambda[Angstroms] < 413)
  ;   we use the photoelectric absorption cross sections of 
  ;   Morrison & McCammon (1983)
  ;     NOTE: these assume solar abundances and no ionization, 
  ;             the appropriate number probably scales linearly with both
  ;   (this is all for the COMPTON THIN regime)
  
  f_003keV = 7.253d15
  f_H_edge = 1.362/3.0 * f_003keV
  f_10keV  = 2.418d18
  ;;if ((f_H_edge LE f0) AND (f_10keV GT f0)) then begin
  if (f_H_edge LE f0) then begin
  
	morrison_photoelec, f0, sigma_temp, ABUNDANCE=METALLICITY_OVER_SOLAR
	SIGMA = SIGMA + sigma_temp ;;/ 0.32 ;;(make sure this uses the total column, not just the 
							   ;;   neutral column, as ions still contribute)
	
  endif


  ; For optical-IR regions, we use the Pei numerical approximations below.
  ;
  ; xsi = tau(lambda)/tau(B) is the ratio of extinction at lambda to the 
  ;    extinction in the B-band. 
  ; k = 10^21 (tau_B / NH)   (NH in cm^2) gives the dimensionless gas-to-dust
  ;    ratio, with k=0.78 for MW, k=0.16 for LMC, k=0.08 for SMC.
  ;    k is INDEPENDENT of the grain properties, and seems to scale rougly
  ;    linearly with metallicity
  ; so, for now, assume solar metallicity and k = k_MW = 0.78. we can rescale later.
  ;
  ; tau_B = ( NH / (10^21 cm^-2) ) * k --> SIGMA_B = k*10^-21  cm^2
  ; tau_lambda = xsi * tau_B --> SIGMA = xsi * SIGMA_B
  ;
  ; k = 0.78 for the MW
  ; k = 0.08 for the SMC, approximately in line with the MW/LMC/SMC metallicity 
  ;  sequence, so we take a k_MW then scaled by the metallicity
  
  k = 0.78 * METALLICITY_OVER_SOLAR
  if (f0 LT f_003keV) then begin
	  lambda_microns = 3.0d14 / f0		; convert frequency [Hz] to wavelength [microns]
	  pei_dustparam, lambda_microns, xsi, SMC=1
	  SIGMA = SIGMA + xsi * k * 1.0d-21
  endif

  ; No double-counting, but I have checked it in detail, and there is really almost
  ;    absolutely no difference -- even up to NH~10^24-25, it's like a factor of 1.1
  ;    or so between including both these factors and not, and only in a very, very 
  ;    narrow frequency range


  ; + Thompson scattering cross sections (NR Compton scattering)
  ;    full compton scattering dsigma/dOmega = (1/2)*r0^2*(1+cos^2(theta)) 
  ;    ( sigma_thompson = (8pi/3)*r0^2 )

  sigma_thompson = 6.65d-25	;; cm^-2
  ;SIGMA = SIGMA + sigma_thompson * 0.5	;; b/c half scattered into viewing angle, 
  										;;      half out (*very* roughly)

  return, SIGMA
end
