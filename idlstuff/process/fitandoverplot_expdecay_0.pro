
;====================================================================


; Fit Exponential to Fading Starburst
;--------------------------------------
pro fitandoverplot_expdecay_0, time, sfr, $
			SFR_SB= SFR_SB, tau=tau, SBtime=SBtime

        x_tofit= time
        y_tofit= sfr
        weight_tofit= 1.0 + 0.0*sfr


        ;  parameters
        ;---------------- 
	; p[0] = Normalization, i.e., the SFR_max
        ; p[1] = Width of the dist.; sigma
        ; p[2] = mean of the dist., i.e., time of SFR_max

        ; initial guess
        guess= [100.0,1.0]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},2)
        ; SFR_SB is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ; tau is greater than 0
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.001


        ; markwardt mpfit procedure
        tau_result = MPFITFUN('func_exp', x_tofit, y_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        ;chi2= strcompress(string(redchi2),/remove_all)
        ;chi2= strmid(chi2,0,4)

        SFR_SB= tau_result[0]
	print, "SFR_SB= ", SFR_SB
        ;rho0lbl= strcompress(string(rho0lbl),/remove_all)
        ;rho0lbl= strmid(rho0lbl,0,5)

        tau= tau_result[1]
	print, "SB tau= ", tau," Gyr"
        ;rslbl= strcompress(string(rslbl),/remove_all)
        ;rslbl= strmid(rslbl,0,5)


        ; overplot fit
        ; ----------
        x=time
        y= func_exp(x,tau_result)
        x=x+SBtime
        oplot, x, y, psym=-3, linestyle= 0, color= 100

        ;xyouts, 0.65, 0.75, 'NFW', /normal, size= 1.0, color=150
        ;xyouts, 0.65, 0.71, 'r!Ds!N='+rslbl, /normal, size= 1.0, color=150
        ;xyouts, 0.77, 0.71, ', c='+clbl, /normal, size= 1.0, color=150


        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end





;====================================================================

