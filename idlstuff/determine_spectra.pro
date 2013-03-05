;-------------------------------------------------------------------------
;
;
;From brobertson@cfa.harvard.edu Wed Jul 26 11:42:17 2006
;Date: Mon, 5 Jun 2006 16:25:53 -0400
;From: Brant Robertson <brobertson@cfa.harvard.edu>
;To: Sukanya Chakrabarti <schakrab@cfa.harvard.edu>
;Cc: Brant Robertson <brobertson@cfa.harvard.edu>,
;     Thomas J. Cox <tcox@cfa.harvard.edu>
;Subject: stellar spectra
;
;Hi Sukanya,
;
;	At the end of this email is an IDL script that will return a stellar  
;spectrum based on an input array of wavelengths in microns, age in  
;gyr, mass in Msun, and metallicity in solar units.  The returned  
;spectrum is log10(L_lambda) in units of L_sun / micron.
;
;	If there is a problem using this script please let me know.  I'm  
;fine with you using this for the far-IR calculation, but not for  
;anything w/ the PAH emission/near-IR.  If you want to use this script  
;for future projects, please ask, or of course you can develop  
;something independently.
;
;Thanks,
;B
;
;-------------------------------------------------------------------------

function stellar_spectrum,lambda,mass,age,zm

;RETURNS UNEXTINCTED STELLAR SPECTRUM
;in alog10(L_lambda) in L_sun / micron
;TAKE   LAMBDA IN microns
;TAKES  AGE    IN GYR (physical, no h)
;TAKES  MASS   IN MSUN (physical, no h)
;TAKES  ZM     IN SOLAR UNITS


nl = n_elements(lambda)
L_lambda = fltarr(nl)

n = n_elements(mass)
spectrum = {n_spec:nl,lambda:fltarr(nl),L_lambda:fltarr(nl)}

;S = CALL_EXTERNAL('/home/brant/code/idl/stellar_spectrum/ 
S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/stellar_spectrum/ 
stellar_spectrum', $
         'main', $
         n, $
         mass, $
         zm, $
         age, $
         nl, $
         lambda, $
         L_lambda)

         spectrum.lambda                 = lambda
         spectrum.L_lambda               = L_lambda

         return,spectrum
end 
  








;-------------------------------------------------------------------------
;From yxli@cfa.harvard.edu Wed Jul 26 11:42:42 2006
;Date: Mon, 17 Jul 2006 19:07:41 -0400
;From: Yuexing Li <yxli@cfa.harvard.edu>
;To: Sukanya Chakrabarti <schakrab@cfa.harvard.edu>
;Cc: Thomas J. Cox <tcox@cfa.harvard.edu>
;Subject: SB99.
;
;Hi Sukanya,
;
;So here's the script to calculate the stellar SED from a snapshot and 
;output an ascii file (lambda in Angstrom, Lum in Lsun). 
;/raid3/yxli/starburst99/idl/cal_spec.pro
;
;It calls a couple of other functions to interpolate the library SED 
;using SB99 (Padova STD metallicity + evolutionary tracks & Kroupa IMF) 
;in /raid3/yxli/starburst99/library/.
;
;A sample of the SED of the quasar at several redshifts is in: 
;/raid3/yxli/starburst99/data/*.dat
;It should be very easy to use. If you have any questions, please let me 
;know.
;
;Cheers,
;Yuexing




