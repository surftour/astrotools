
;====================================================================




; Fit Exponential to Fading Starburst
;--------------------------------------
pro fitandoverplot_expdecay, time, sfr, $
                        SFR_SB= SFR_SB, tau=tau, SBtime=SBtime, decay_const=decay_const, $
                        quiet=quiet

        x_tofit= time
        y_tofit= sfr
        ;weight_tofit= 1.0 + 0.0*sfr
        weight_tofit= 0.1*sfr + 1.0


        ;  parameters
        ;---------------- 
        ; p[0] = Normalization, i.e., the SFR_max ( - p[2])
        ; p[1] = Exp. decay, i.e., tau
        ; p[2] = constant, so it decays to this, rather than 0

        ; initial guess
        guess= [100.0,1.0,0.1]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
        ; SFR_SB is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ; tau is greater than 0
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.001
        ; const is greater than 0
        pi[2].limited(0) = 1
        pi[2].limits(0) = 0.00001


        ; markwardt mpfit procedure
        tau_result = MPFITFUN('func_exp_sf', x_tofit, y_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)


        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        ;chi2= strcompress(string(redchi2),/remove_all)
        ;chi2= strmid(chi2,0,4)

        SFR_SB= tau_result[0]
        print, "SFR_SB= ", SFR_SB, " M_solar/Yr"
        ;rho0lbl= strcompress(string(rho0lbl),/remove_all)
        ;rho0lbl= strmid(rho0lbl,0,5)

        tau= tau_result[1]
        print, "SB tau= ", tau," Gyr"
        taulbl= strcompress(string(tau),/remove_all)
        taulbl= strmid(taulbl,0,4)

        decay_const= tau_result[2]
        print, "decay const= ", decay_const," M_solar/Yr"
        ;rslbl= strcompress(string(rslbl),/remove_all)
        ;rslbl= strmid(rslbl,0,5)


        ; overplot fit
        ; ----------
        x=time
        y= func_exp_sf(x,tau_result)
        x=x+SBtime
        if not keyword_set(quiet) then begin
                oplot, x, y, psym=-3, linestyle= 0, color= 200, thick= 6.0
        endif

        ;xyouts, 0.65, 0.75, 'NFW', /normal, size= 1.0, color=150
        ;xyouts, 0.65, 0.71, 'r!Ds!N='+rslbl, /normal, size= 1.0, color=150
        ;xyouts, 0.77, 0.71, ', c='+clbl, /normal, size= 1.0, color=150


        if not keyword_set(quiet) then begin
                ;xyouts, 0.65, 0.71, 'Tau Model', color=0, size=1.5, /normal
                ;xyouts, 0.65, 0.66, '!7s!6!DSF!N= '+taulbl+' Gyr', color=0, size=1.5, /normal
        endif

        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end





;====================================================================




