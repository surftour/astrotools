; -------------------------------
;  compute average profile
;
;  This is a general procedure to
; total y in each bin along x.
; -------------------------------
pro process_prof_tot, x_var, y_var, bins, xmax, xmin, $
			radial_vals, $
			tot_mass, $
			tot_1sig, $
                        y_weighting=y_weighting

	radial_vals= fltarr(bins)
	tot_mass= fltarr(bins)
	tot_1sig= fltarr(bins)

        r_l= x_var

        binsize = float((xmax-xmin))/bins

        for i=1,bins do begin 
           lg_r = i*binsize + xmin
           sm_r = (i-1)*binsize + xmin
           radial_vals(i-1) = 0.5*(lg_r + sm_r)

           idx= where((r_l GE sm_r) AND (r_l LT lg_r))
           if idx(0) lt 0 then begin
                m= [0.0]
                weight= 1.0
           endif else begin
                m= y_var(idx)
                weight= n_elements(idx)
                if keyword_set(y_weighting) then begin
                   m= y_var(idx)*y_weighting(idx)
                   weight= total(y_weighting(idx))
                endif
           endelse
           tot_mass(i-1)= total(m)
	   tot_1sig(i-1)= tot_mass(i-1)/sqrt(weight)
        endfor

end



