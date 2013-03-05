
;----------------------------
; Fit de Vacauluers profile
;----------------------------
pro fitandoverplot_devacprofile, radius, density, weight, $
                                        x_is_devac=x_is_devac, $
                                        x_is_log=x_is_log, $
                                        ylogaxis=ylogaxis

        r_tofit= radius
        mass_tofit= density
        weight_tofit= weight


        ; parameters
        ;---------------- 
        ; [0] - I_0             ; normalization
        ; [1] - R_e             ; core radius

        ; fitting to I(R), the projected
        ; surface density.  use formula (1)
        ; from Hernquist (1990)

        ; initial guess
        guess= [1.0d+5,3.0]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},2)
        ; rho_crit*rho_0 is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.0
        ;r_core is greater than 0.01 and less than 500
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.01
        pi[1].limited(1) = 1
        pi[1].limits(1) = 100.0



        ; markwardt mpfit procedure
        devac_result = MPFITFUN('func_devac', r_tofit, mass_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        chi2_beta= strcompress(string(redchi2),/remove_all)
        chi2_beta= strmid(chi2_beta,0,4)


        I_0= devac_result[0]
        print, "I_0= ", I_0
        Ilbl= strcompress(string(I_0),/remove_all)
        Ilbl= strmid(Ilbl,0,5)

        r_e= devac_result[1]
        print, "R_e= ", r_e
        relbl= strcompress(string(r_e),/remove_all)
        relbl= strmid(relbl,0,5)


        ; plot the fit
	; --------------
        if keyword_set(x_is_devac) then x= r_tofit^(1/4.)
        if keyword_set(x_is_log) then x= alog10(r_tofit)

        y= func_devac(r_tofit,devac_result)
        if keyword_set(ylogaxis) then y=alog10(y)

        oplot, x, y, psym=-3, linestyle= 0, color= 0, thick=3.0

        ;xyouts, 0.65, 0.65, '!7b!3-Model', /normal, size= 1.0, color=200
        ;xyouts, 0.65, 0.61, 'r!Dcore!N='+rclbl, /normal, size= 1.0, color=200
        ;xyouts, 0.79, 0.61, ', !7b!3='+blbl, /normal, size= 1.0, color=200

        ; reduced chi squared
        ;xyouts, 0.45, 0.22, '!7v!3!E2!N='+chi2_beta, /normal, size= 1.0, color=200


end


