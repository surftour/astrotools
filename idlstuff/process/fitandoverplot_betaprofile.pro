
;--------------------
; Fit Beta Profile
;--------------------
pro fitandoverplot_betaprofile, radius, density, weight, $
                                        ylogaxis=ylogaxis

        r_tofit= radius
        mass_tofit= density
        weight_tofit= weight


        ; parameters
        ;---------------- 
        ; [0] - rho_0             ; normalization
        ; [1] - r_core            ; core radius

; -> currently using option #1
        ; [2] - option #1 : 3.0 * beta - 0.5
        ; [2] - option #2 : beta


        ; initial guess
        guess= [1.0d+38,10.0,1.0]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
        ; rho_crit*rho_0 is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.0
        ;r_core is greater than 0.01 and less than 500
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.01
        pi[1].limited(1) = 1
        pi[1].limits(1) = 100.0
        ; beta is greater than 0
        pi[2].limited(0) = 1
        pi[2].limits(0) = 0.0



        ; markwardt mpfit procedure
        beta_result = MPFITFUN('func_beta', r_tofit, mass_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)


        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        chi2_beta= strcompress(string(redchi2),/remove_all)
        chi2_beta= strmid(chi2_beta,0,4)


        rho_0= beta_result[0]
        print, "Rho_0= ", rho_0
        rhoclbl= strcompress(string(rho_0),/remove_all)
        rhoclbl= strmid(rhoclbl,0,5)

        r_core= beta_result[1]
        print, "r_core= ", r_core
        rclbl= strcompress(string(r_core),/remove_all)
        rclbl= strmid(rclbl,0,5)


        betafit= (beta_result[2]+0.5)/3.0      ; option 1
        ;betafit= beta_result[2]               ; option 2
        print, "Beta_fit= ", betafit
        blbl= strcompress(string(betafit),/remove_all)
        blbl= strmid(blbl,0,4)


        ; draw fits
        x= findgen(50)/40.0
        x=alog10(min(r_tofit)) + (alog10(max(r_tofit))-alog10(min(r_tofit)))*x
        x_notlog= 10^(x)

        y= func_beta(x_notlog,beta_result)
        if keyword_set(ylogaxis) then y=alog10(y)
        oplot, x_notlog, y, psym=-3, linestyle= 0, color= 0, thick=3.0

        ;xyouts, 0.65, 0.65, '!7b!3-Model', /normal, size= 1.0, color=200
        ;xyouts, 0.65, 0.61, 'r!Dcore!N='+rclbl, /normal, size= 1.0, color=200
        ;xyouts, 0.79, 0.61, ', !7b!3='+blbl, /normal, size= 1.0, color=200

        ; reduced chi squared
        ;xyouts, 0.45, 0.22, '!7v!3!E2!N='+chi2_beta, /normal, size= 1.0, color=200


end


