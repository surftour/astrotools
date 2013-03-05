; -------------------------------
;  compute average profile
;
;  This is a general procedure to
; average the y variables along the
; x range.
;
; -------------------------------
pro process_prof_avg, x_var, y_var, bins, xmax, xmin, $
			radial_vals, $
			avg_s, $
			avg_rootn, $
			avg_disp=avg_disp, $
                        y_weighting=y_weighting, $
			xaxis_min_sn=xaxis_min_sn

	radial_vals= fltarr(bins)
	avg_s= fltarr(bins)
	avg_rootn= fltarr(bins)
	avg_disp= fltarr(bins)

        binsize = float((xmax-xmin))/bins
	bini= 0

        for i=1,bins do begin 

	   ; min radius of bin
           sm_r = (i-1)*binsize + xmin
	   lg_r= sm_r

           add_a_bin_to_step:
	   ; max radius of bin
           lg_r = lg_r + binsize

	   
           idx= where((x_var GE sm_r) AND (x_var LT lg_r))
	   if keyword_set(xaxis_min_sn) and n_elements(idx) lt 300 then begin
	   ;if keyword_set(xaxis_min_sn) and n_elements(x_var)/n_elements(idx) gt (bins*200) then begin
		print, bini, i, sm_r, lg_r, n_elements(idx), n_elements(x_var)
		i= i + 1
		goto, add_a_bin_to_step
	   endif

           if idx(0) lt 0 then begin
                m= [0.0, 0.0, 0.0]
                weight= 3.0
           endif else begin
                m= y_var(idx)
                weight= n_elements(idx)
                if keyword_set(y_weighting) then begin
                   m= y_var(idx)*y_weighting(idx)
                   weight= total(y_weighting(idx))
                endif
           endelse

	   ; successfully computed bin, save info and move on
           radial_vals(bini) = 0.5*(lg_r + sm_r)
           avg_s(bini)= total(m)/weight
	   avg_rootn(bini)= avg_s(bini)/sqrt(n_elements(idx))
	   avg_disp(bini)= sqrt((moment(m))[1])
	   bini= bini + 1
        endfor

end



