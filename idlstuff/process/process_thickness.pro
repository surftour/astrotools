; -------------------------------
;  compute average profile
;
;  This is a general procedure to
; average the y variables along the
; x range.
;
; -------------------------------
pro process_thickness, x_var, y_var, bins, xmax, xmin, $
			radial_vals, $
			thickness, $
                        y_weighting=y_weighting

	radial_vals= fltarr(bins)
	z_avg= fltarr(bins)
	zhalfmass= fltarr(bins)

        r_l= x_var

        binsize = float((xmax-xmin))/bins

        for i=1,bins do begin 
           lg_r = i*binsize + xmin
           sm_r = (i-1)*binsize + xmin
           radial_vals(i-1) = 0.5*(lg_r + sm_r)

           idx= where((r_l GE sm_r) AND (r_l LT lg_r))
           if idx(0) lt 0 then begin
                z_avg= [0.0]
                zhalfmass= 0.0
		weight= 1.0
           endif else begin
                z= y_var(idx)
		if keyword_set(y_weighting) then m= y_weighting(idx) else m= 0.0*z + 1.0
                weight= total(y_weighting(idx))
           endelse

	   ; check some simple quantities
	   ;print, "bin, N, Sigma= ", i, n_elements(idx), total(m)/(!PI*(lg_r*lg_r - sm_r*sm_r))

	   ; what is the average disk height
           z_avg(i-1)= total(z)/weight

	   ; now calculate the half-mass height
	   mhalf= 0.5 * weight
	   zabs=abs(z)

	   z_guess = 0.2
	   n_iterations= 0

	   repeat begin
	     idxx= where(zabs lt z_guess)
             if idxx(0) ne -1 then m_guess = total(m(idxx)) else m_guess= 0.0

             dm = m_guess-mhalf
             z_guess=z_guess * (1.0 - 0.5*dm/weight)
             n_iterations= n_iterations+1
             last_z_guess= z_guess
             second_to_last_z_guess= last_z_guess
             ;if second_to_last_z_guess eq z_guess then z_guess = z_guess + long((dm/abs(dm)) *  0.1 * z_guess)
             ;print, n_iterations, m_guess, mhalf, dm, z_guess
	   endrep until ((abs(dm/weight) lt 0.01) or (n_iterations gt 100))


	   zhalfmass(i-1)= z_guess
        endfor


	thickness= zhalfmass

end



