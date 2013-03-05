pro process_kenn_line, a, b, c, d, bins, nmass_sd, nsfr_sd


        sorta = sort(a)
        r = a(sorta)
        sfr = b(sorta)
        mass = c(sorta)
        mfs = d(sorta)

        sfr_sd = fltarr(bins)
        mass_sd = fltarr(bins)
        binsize = n_elements(r)/bins

        if binsize gt 0 then begin

        rs = fltarr(binsize)

        ; average annuli 
        ; ---------------
        for i=1,bins do begin
           rs(*) = r((i-1)*binsize:i*binsize-1)
           lg_r = max(rs)
           sm_r = min(rs)

           gm_inannulus = mass((i-1)*binsize:i*binsize-1)-mfs((i-1)*binsize:i*binsize-1)
           sfr_inannulus = sfr((i-1)*binsize:i*binsize-1)

           sfr_avg = total(sfr_inannulus)
           mass_tot = total(gm_inannulus)

           sd_kpc2 = (3.14159)*( (lg_r)^2 - (sm_r)^2 )
           sd_pc2 = sd_kpc2*(1e+6)
           mass_sd(i-1) = mass_tot*(1e10)/sd_pc2
           sfr_sd(i-1) = sfr_avg/sd_kpc2
        endfor

        endif



        ; check against non-star formation runs
        ; -------------------------------------
        idx = where(sfr_sd GT 0)
        if (idx(0) LT 0) then begin
           nmass_sd = [1e10]
           nsfr_sd = [1e10]
        endif else begin
           nmass_sd = mass_sd(idx)
           nsfr_sd = sfr_sd(idx)
        endelse

end




