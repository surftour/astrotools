
;======================================================================

; Fit Gaussian to Starburst
;-----------------------------
pro fitandoverplot_gaussian, time, sfr, $
			SFR_max= SFR_max, width=width, SBtime=SBtime, $
			quiet=quiet

        x_tofit= time
        y_tofit= sfr
        weight_tofit= 1.0 + 0.0*sfr
        ;weight_tofit= 1.0 + 0.1*sfr


        ;  parameters
        ;---------------- 
	; p[0] = Normalization, i.e., the SFR_max
        ; p[1] = Width of the dist.; sigma
        ; p[2] = mean of the dist., i.e., time of SFR_max

        ; initial guess
        ;guess= [100.0,0.2,1.5]
	normz= max(sfr)
        idd= where(sfr eq normz)
        meanz= time(idd[0])
        wid= max(time)/3.0
        guess= [normz, wid, meanz]
print, guess

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
        ; Normalization is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ; width is greater than 0
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.001
        ; mean is greater than 0
        pi[2].limited(0) = 1
        pi[2].limits(0) = 0.001


        ; markwardt mpfit procedure
        sb_result = MPFITFUN('func_gaussian', x_tofit, y_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

        print, "--------------------------"
        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        chi2= strcompress(string(redchi2),/remove_all)
        chi2= strmid(chi2,0,4)

        SFR_max= sb_result[0]
        print, "SFR_max= ", SFR_max
        sfrmaxlbl= strcompress(string(SFR_max),/remove_all)
        sfrmaxlbl= strmid(sfrmaxlbl,0,5)

        width= sb_result[1]
        print, "Starburst Width= ", width," Gyr"
        wlbl= strcompress(string(width),/remove_all)
        wlbl= strmid(wlbl,0,4)

        SBtime= sb_result[2]
        print, "Time at SFR_max= ", SBtime," Gyr"
        sbtm= strcompress(string(SBtime),/remove_all)
        sbtm= strmid(sbtm,0,4)


        ; overplot fit
        ; ----------
        ;x=min(time) - 0.2 + (0.2+max(time)-min(time))*findgen(100)/100.0
        ;x= SBtime + 1.5*(findgen(100)/100.0 - 0.5)
        x= SBtime + 4.0*(findgen(100)/100.0 - 0.5)
        y= func_gaussian(x,sb_result)

	idx= where(y gt 0.3*min(sfr))
	y= y(idx)
	x= x(idx)

	if not keyword_set(quiet) then begin
        	oplot, x, y, psym=-3, linestyle= 2, color= 50, thick= 6.0
	endif


        if not keyword_set(quiet) then begin
                ;xyouts, 0.22, 0.75, 'Gaussian ('+chi2+')', size=1.2, color=0, /normal
                ;xyouts, 0.23, 0.71, 'SFRmax= '+sfrmaxlbl, size=1.0, color=0, /normal
                ;xyouts, 0.23, 0.67, 'sigma=  '+wlbl, size=1.0, color=0, /normal
                ;xyouts, 0.23, 0.63, 'sbtime= '+sbtm, size=1.0, color=0, /normal
                ;xyouts, 0.65, 0.84, 'Gaussian', size=1.6, color=0, /normal
                ;xyouts, 0.65, 0.80, '!7r!6= '+wlbl+' Gyr', size=1.5, color=0, /normal
        endif

        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end




;====================================================================


