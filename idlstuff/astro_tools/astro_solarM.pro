;------------------------------
;
;  Returns the magnitude of the
;  Sun in various bands.
;
;
;------------------------------

function astro_solarM, band

if not keyword_set(band) then band='B'

;VEGA system
;from www.ucolick.org/~cnaw/sun.html
;following Fukugita et al. 1995, PASP, 105, 945
s_UBVRIJHK= fltarr(8)
s_UBVRIJHK(0) = 5.56;  //U (BESSEL)
s_UBVRIJHK(1) = 5.45;  //B (BESSEL)
s_UBVRIJHK(2) = 4.80;  //V (BESSEL)
s_UBVRIJHK(3) = 4.46;  //R (KPNO)
s_UBVRIJHK(4) = 4.10;  //I (KPNO)
s_UBVRIJHK(5) = 3.66;  //J (BESSEL)
s_UBVRIJHK(6) = 3.32;  //H (BESSEL)
s_UBVRIJHK(7) = 3.28;  //K (BESSEL)

;AB  magnitudes -> SDSS unprimed
;from www.ucolick.org/~cnaw/sun.html
;following Fukugita et al. 1995, PASP, 105, 945
s_ugrizJHK = fltarr(5)
s_ugrizJHK(0) = 6.75; //u SDSS
s_ugrizJHK(1) = 5.33; //g SDSS
s_ugrizJHK(2) = 4.67; //r SDSS
s_ugrizJHK(3) = 4.48; //i SDSS
s_ugrizJHK(4) = 4.42; //z SDSS



    
case band of
   'U': this_solar_mag= s_UBVRIJHK(0)
   'B': this_solar_mag= s_UBVRIJHK(1)
   'V': this_solar_mag= s_UBVRIJHK(2)
   'R': this_solar_mag= s_UBVRIJHK(3)
   'I': this_solar_mag= s_UBVRIJHK(4)
   'J': this_solar_mag= s_UBVRIJHK(5)
   'H': this_solar_mag= s_UBVRIJHK(6)
   'K': this_solar_mag= s_UBVRIJHK(7)
   'u': this_solar_mag= s_ugrizJHK(0)
   'g': this_solar_mag= s_ugrizJHK(1)
   'r': this_solar_mag= s_ugrizJHK(2)
   'i': this_solar_mag= s_ugrizJHK(3)
   'z': this_solar_mag= s_ugrizJHK(4)
endcase


print, "Band= ", band
print, "Solar Magnitude= ", this_solar_mag


return, this_solar_mag

end



