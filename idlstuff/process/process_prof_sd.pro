; -----------------------------------------------------------------
;  compute a radial profile
;
;  Note that this procedure is completely general, all it
;  needs is x_var, y_var, the number of bins and the max
;  and min x range.  It then splits the x range into equal
;  bin
;
;     delta_x= (xmax - xmin) / bins
;
;  and then steps along these x bins to calculate the 
;  surface density.  Note that in order to get the area 
;  correct, we need to know if x is log or r^1/4.
;
;  We also allow the y variable to be weighted by a variable
;  that can be provided.
;
;  Returned:
;
;    radial_vals - the radial bins
;    surface_density - at each radial bin point
;    sd_1sig - error
;  
; -----------------------------------------------------------------
pro process_prof_sd, x_var, y_var, bins, xmax, xmin, $ 
                                        radial_vals, surface_density, $
					sd_1sig=sd_1sig, $
                                        x_is_log=x_is_log, $
                                        x_is_devac=x_is_devac, $
                                        y_weighting=y_weighting

	radial_vals= fltarr(bins)
	surface_density= fltarr(bins)
	sd_1sig= fltarr(bins)

        ;r_l= alog10(radius)    ; actually surface density is linear
        r_l= x_var

        binsize = float((xmax-xmin))/bins

        for i=1,bins do begin 
           lg_r = i*binsize + xmin
           sm_r = (i-1)*binsize + xmin
           radial_vals(i-1) = 0.5*(lg_r + sm_r)

           idx= where((r_l GE sm_r) AND (r_l LT lg_r))
           weight= 1.0
           if idx(0) lt 0 then begin
                m= [0.0]
           endif else begin
                m= y_var(idx)
                if keyword_set(y_weighting) then begin
                   m= y_var(idx)*y_weighting(idx)
                   weight= total(y_weighting(idx))/n_elements(y_weighting(idx))
                endif
           endelse
           ;thisarea= !PI*((10^(lg_r))^2 - (10^(sm_r))^2)       ; (kpc/h)^2
           thisarea= !PI*((lg_r)^2 - (sm_r)^2)       ; (kpc/h)^2
           if keyword_set(x_is_log) then thisarea= !PI*((10^(lg_r))^2 - (10^(sm_r))^2)
           if keyword_set(x_is_devac) then thisarea= !PI*((lg_r)^8 - (sm_r)^8)

           ;surface_density(i-1)= total(m)*1e10/thisarea                     ; h Msolar kpc-2
           surface_density(i-1)= total(m)/(thisarea*weight)          ; Gadget units - h 10^10 Msolar kpc-2

	   ;sd_1sig(i-1)= surface_density(i-1)/sqrt(n_elements(idx))
	   sd_1sig(i-1)= 2.0*surface_density(i-1)/sqrt(n_elements(idx))

	   ;print, i, n_elements(idx), surface_density(i-1), sd_1sig(i-1)
        endfor

end



