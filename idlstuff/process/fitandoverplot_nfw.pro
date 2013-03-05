;--------------------
; Fit NFW Profile
;--------------------
pro fitandoverplot_nfw, radius, density, weight, xmin=xmin, xmax=xmax, $
				divide_by_rhocrit=divide_by_rhocrit

        r_tofit= radius
        mass_tofit= density
        weight_tofit= weight

        ; do we fit only a certain regions?
        ; -----------------------------------
        ;rmintofit= 100
        ;print, "fitting to r=",rmintofit," only"
        ;idx=where(r_tofit gt rmintofit)
        ;if idx(0) ne -1 then begin
        ;    r_tofit= r_tofit(idx)
        ;    mass_tofit= mass_tofit(idx)
        ;    weight_tofit= weight_tofit(idx)
        ;endif


        ;  parameters
        ;---------------- 
        ; [0] - rho_0 
        ; [1] - r_s

        ; initial guess
        guess= [10.0,20.0]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},2)
        ; rho_0 is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ;r_s is greater than 0 and less than 100
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.001
        pi[1].limited(1) = 1
        pi[1].limits(1) = 100.0


        ; markwardt mpfit procedure
        nfw_result = MPFITFUN('func_nfw', r_tofit, mass_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        chi2_nfw= strcompress(string(redchi2),/remove_all)
        chi2_nfw= strmid(chi2_nfw,0,4)


        ; determine concentration
        x= [12.0]
        rho_0= nfw_result[0]

	if keyword_set(divide_by_rhocrit) then begin
                rho_crit= fload_rho_crit(1)
                rho_0= rho_0 * rho_crit
        endif

        r_s= nfw_result[1]
        cc= NEWTON(x, 'func_c', /double)
        R_vir= cc*r_s
        M_vir= 4.0*!PI*rho_0*r_s*r_s*r_s*(alog(1+cc) - cc/(1+cc))

	print, "R_vir= ", R_vir
	print, "M_vir= ", M_vir

        rho0lbl= nfw_result[0]
        rho0lbl= strcompress(string(rho0lbl),/remove_all)
        rho0lbl= strmid(rho0lbl,0,5)

        rslbl= nfw_result[1]
        rslbl= strcompress(string(rslbl),/remove_all)
        rslbl= strmid(rslbl,0,5)

        clbl= cc
        clbl= strcompress(string(clbl),/remove_all)
        clbl= strmid(clbl,0,3)

        ; overplot fit
        ; ----------
        ;x=xmin + (xmax-xmin)*findgen(100)/100.0
        x=alog10(xmin) + (alog10(xmax)-alog10(xmin))*findgen(100)/100.0
        x_notlog= 10^(x)
	x= x_notlog

        y= func_nfw(x_notlog,nfw_result)
        ;if keyword_set(divide_by_rhocrit) then begin
        ;        rho_crit= fload_rho_crit(1)
        ;        y= y/rho_crit
        ;endif
        oplot, x, y, psym=-3, linestyle= 0, color= 150

        ;xyouts, 0.45, 0.38, 'NFW', /normal, size= 1.0, color=150
        ;xyouts, 0.45, 0.34, 'r!Ds!N='+rslbl, /normal, size= 1.0, color=150
        ;xyouts, 0.45, 0.30, 'c='+clbl, /normal, size= 1.0, color=150
        ;xyouts, 0.77, 0.71, ', c='+clbl, /normal, size= 1.0, color=150

        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end


