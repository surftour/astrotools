; ----------------------------
;  compute density profile
;   in logarithmis bins
;
;
; ----------------------------
pro process_prof_rho, radius, mass, bins, xmax, xmin,  r_s, rho, weight, $
			x_is_log=x_is_log


        ;sorta = sort(radius)
        ;r = radius(sorta)
	r = radius
        ;mass = mass(sorta)
        weight = 0.0*rho

        binsize = float((xmax-xmin))/bins

        for i=1,bins do begin
           lg_r = i*binsize + xmin
           sm_r = (i-1)*binsize + xmin
           r_s(i-1) = 0.5*(lg_r + sm_r)

           idx= where((r GE sm_r) AND (r LT lg_r))
           if idx(0) lt 0 then m= [0.0] else m= mass(idx)
           thisvolume= 4.0*(!PI)*( (lg_r)^3 - (sm_r)^3 )/3.0
           if keyword_set(x_is_log) then thisvolume= 4.0*(!PI)*( (10^(lg_r))^3 - (10^(sm_r))^3 )/3.0
           rho(i-1)= total(m)/thisvolume

           nwt= n_elements(idx)
           ;if nwt lt 2 then nwt= 2
           weight(i-1)= rho(i-1)/sqrt(nwt)
           ;if r_s(i-1) lt alog10(smoothlen) then weight(i-1)= rho(i-1)
;print, i, sm_r, lg_r, r_s(i-1), nwt, rho(i-1)
        endfor


        ; take out any zero's
        ; --------------------
        ;idx= where(rho GT 0)
        ;r_s= r_s(idx)
        ;rho= rho(idx)

end


