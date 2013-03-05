
;--------------------
; Fit Sersic Profile
;--------------------
pro fitandoverplot_sersicprofile, radius, density, weight, $
                                        ylogaxis=ylogaxis, $
                                        ymagaxis=ymagaxis, $
                                        minfitradius=minfitradius, $
                                        sersic_n=sersic_n, $
                                        excesslight=excesslight, $
                                        r_e=r_e, $
                                        x_is_devac=x_is_devac, $
                                        x_is_log=x_is_log, $
                                        dontoverplot=dontoverplot

        if keyword_set(minfitradius) then begin
           idx=where(radius gt minfitradius)
           r_tofit= radius(idx)
           mass_tofit= density(idx)
           weight_tofit= weight(idx)
        endif else begin
           r_tofit= radius
           mass_tofit= density
           weight_tofit= weight
        endelse

        idx=where(mass_tofit le 0.0)
        if idx(0) ne -1 then begin
                idx= where(mass_tofit gt 0.0)
                r_tofit= r_tofit(idx)
                mass_tofit= mass_tofit(idx)
                weight_tofit= weight_tofit(idx)
        endif


        ; parameters
        ;---------------- 
        ; [0] - I_0             ; normalization
        ; [1] - R_e             ; effective radius
        ; [2] - n               ; sersic index n


        ; initial guess
        guess= [1.0d+5,1.0,2.0]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
        ; I_0 is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.0
        ; R_e is greater than 0.01 and less than 500
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.01
        pi[1].limited(1) = 1
        pi[1].limits(1) = 100.0
        ; n is greater than 0.5 and less than 10
        pi[2].limited(0) = 1
        pi[2].limits(0) = 0.5
        pi[2].limited(1) = 1
        pi[2].limits(1) = 10.0



        ; markwardt mpfit procedure
        sersic_result = MPFITFUN('func_sersic', r_tofit, mass_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        chi2_beta= strcompress(string(redchi2),/remove_all)
        chi2_beta= strmid(chi2_beta,0,4)


        I_0= sersic_result[0]
        print, "I_0= ", I_0
        Ilbl= strcompress(string(I_0),/remove_all)
        Ilbl= strmid(Ilbl,0,5)

        r_e= sersic_result[1]    ; dimensionless constant in sersic_func makes this R_e
        print, "R_e= ", r_e
        relbl= strcompress(string(r_e),/remove_all)
        relbl= strmid(relbl,0,5)
xyouts, 0.65, 0.65, relbl, /normal, charthick=3.0, size=1.2, color=0


        sersicn= sersic_result[2]
        sersic_n= sersicn
        print, "Sersic n= ", sersicn
        sersicnlbl= strcompress(string(sersicn),/remove_all)
        sersicnlbl= 'Sersic n= '+strmid(sersicnlbl,0,4)
xyouts, 0.65, 0.70, sersicnlbl, /normal, charthick=3.0, size=1.2, color=0


        ; draw fits
        x= radius

	x_toplot= x
        if keyword_set(x_is_devac) then x_toplot= x^(1/4.)
        if keyword_set(x_is_log) then x_toplot= alog10(x)

        y= func_sersic(x,sersic_result)
        if keyword_set(ylogaxis) then y=alog10(y)
        if keyword_set(ymagaxis) then y=-2.5*alog10(y)
        if not keyword_set(dontoverplot) then oplot, x_toplot, y, psym=-3, linestyle= 0, color= 0, thick=4.0



        ; calculate excess (or deficient) light
        Nr= n_elements(radius)

	dr_xaxis= (radius[Nr-1]-radius[0])/Nr
        if keyword_set(x_is_devac) then dr_xaxis= (radius[Nr-1])^(0.25) - (radius[0])^(0.25)
        if keyword_set(x_is_log) then dr_xaxis= alog10(radius[Nr-1]) - alog10(radius[0])
        dr_xaxis= dr_xaxis/Nr

        area= fltarr(Nr)

        for i=0,Nr-1 do begin
		thisr= radius[i]
                if keyword_set(x_is_devac) then thisr= (radius[i])^(0.25)
                if keyword_set(x_is_log) then thisr= alog10(radius[i])
                smr= thisr - 0.5*dr_xaxis
                lgr= thisr + 0.5*dr_xaxis

                if keyword_set(x_is_devac) then smr= smr^(4.)
                if keyword_set(x_is_log) then smr= 10^smr

                if keyword_set(x_is_devac) then lgr= lgr^(4.)
                if keyword_set(x_is_log) then lgr= 10^lgr

                area[i]= !PI * (lgr*lgr - smr*smr)
        endfor

        total_light= total(density*area)

        if keyword_set(ylogaxis) then y=10^y
        fit_light= total(y*area)

        ; 1.0e+6 comes from m_sun/pc2 conversion to kpc
        print, "*********** "
        print, "Total Light= ", total_light*1.0e+6
        print, "Fit Light= ", fit_light*1.0e+6, "   (",100.0*(total_light-fit_light)/total_light," % )"
        print, " "

        ; don't analyze outer profile
        ;idx=where(radius lt 5.0)
        idx=where(radius lt 2.0)
        total_2kpc_light= total(density(idx)*area(idx))
        fit_light= total(y(idx)*area(idx))

        print, "*********** "
        ;print, " just within 5 kpc/h"
        print, " just within 2 kpc/h"
        print, "Total2kpc_ Light= ", total_2kpc_light*1.0e+6
        print, "Fit Light= ", fit_light*1.0e+6
        print, " "
        excesslight= 100.0*(total_2kpc_light-fit_light)/total_light
        print, " excesslight= ", excesslight
        print, " "

end


