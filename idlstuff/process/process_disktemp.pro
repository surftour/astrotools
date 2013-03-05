; -------------------------------
;  compute profile of disk sigmaz
;
;  This is a general procedure to
; average the y variables along the
; x range.
;
; -------------------------------
pro process_disktemp, x_var, y_var, bins, xmax, xmin, $
			radial_vals, $
			disktemp, $
                        y_weighting=y_weighting

	radial_vals= fltarr(bins)
	vz_avg= fltarr(bins)
	vz_disp= fltarr(bins)

        r_l= x_var

        binsize = float((xmax-xmin))/bins

        for i=1,bins do begin 
           lg_r = i*binsize + xmin
           sm_r = (i-1)*binsize + xmin
           radial_vals(i-1) = 0.5*(lg_r + sm_r)

           idx= where((r_l GE sm_r) AND (r_l LT lg_r))
           if idx(0) lt 0 then begin
                vz_avg= [0.0]
                vz_disp= 0.0
           endif else begin
		vz= y_var(idx)
		m= y_weighting(idx)
		;print, "<z>, sigma_z= ", mean(vz), sqrt(variance(vz))
		if keyword_set(y_weighting) then begin
			vz_avg(i-1)= total(vz*m)/total(m)
			vz_disp(i-1)= sqrt(total(((vz-vz_avg(i-1))*(vz-vz_avg(i-1)))*m)/total(m))
		endif else begin
			; mean vz
			vz_avg(i-1)= total(vz)/n_elements(vz)
			vz_disp(i-1)= sqrt(total((vz-vz_avg(i-1))*(vz-vz_avg(i-1)))/total(n_elements(vz)))
		endelse
           endelse
        endfor

	disktemp= vz_disp

end



