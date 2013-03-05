;==================================================================
; procedure to determine the effective radius
; given the x, y, and mass of all particles.  Note that this
; does the strict half-mass radius.
;
;
;
; --------------------------------------------
function process_2d_halfmass, x, y, mass, center=center, $
			do_quarter_mass=do_quarter_mass, $
			fixed_mass=fixed_mass, $
			FindMassFactor=FindMassFactor

        if not keyword_set(center) then c= [0,0] else c= center
	if not keyword_set(FindMassFactor) then FindMassFactor= 0.5

        rxy= sqrt((x-c[0])*(x-c[0]) + (y-c[1])*(y-c[1]))

        mtot = total(mass)
        if keyword_set(do_quarter_mass) then mtot= mtot*0.5
	if keyword_set(fixed_mass) then mtot= fixed_mass * 2.0    ; (corrects for the 0.5 below)
        if mtot le 0.0 then begin
                Reff= 0.0
                return, 0
        endif

        ;hm = 0.5*mtot
        hm = FindMassFactor*mtot
        nidx = n_elements(mass)
        print, " "
        print, " using FindMassFactor= ", FindMassFactor
        print, " total mass= ", mtot
        print, " hm= ", hm
        print, " "

        ; ok, preparations done
        sorta= sort(rxy)
        r= rxy(sorta)
        m= mass(sorta)

        ; find, effective radius
        n_guess = long(nidx/2.0)
	last_n_guess= 0L
	second_to_last_n_guess= 0L
        max_n_guess= long(1e6)
	min_n_guess= -100L
        n_iterations= 0

        if keyword_set(do_quarter_mass) then n_guess= long(n_guess/2.0)
        if keyword_set(do_quarter_mass) then nidx= nidx/2

        repeat begin
             m_guess = total(m[0:n_guess])
             dm = m_guess-hm
	     ;print, n_iterations, m_guess, mtot, hm, dm, n_guess, last_n_guess, second_to_last_n_guess, max_n_guess, min_n_guess

	     ; record guesses
	     second_to_last_n_guess= last_n_guess
	     last_n_guess= n_guess
	     if dm gt 0 and n_guess lt max_n_guess then max_n_guess= n_guess
	     if dm lt 0 and n_guess gt min_n_guess then min_n_guess= n_guess

	     ; now we set about finding the new guess
             n_guess=long(n_guess - nidx*0.5*dm/mtot)

	     ; do some checking to make sure we don't go negative and aren't in a crazy loop
	     if n_guess lt 1 then n_guess= 1L
	     ;if n_guess eq second_to_last_n_guess then n_guess= long(n_guess * 0.8)
	     ;if second_to_last_n_guess eq n_guess then n_guess = n_guess + long((dm/abs(dm)) *  0.1 * n_guess)
	     ;if (n_guess eq 1L) and (second_to_last_n_guess eq 1L) then n_guess= long(last_n_guess * 0.8)
	     if dm gt (0.5*hm) then n_guess= long(last_n_guess * 0.9)
	     if dm lt (-0.5*hm) then n_guess= max([long(last_n_guess * 1.1),last_n_guess+25L])
	     ;if n_guess ge max_n_guess then n_guess= max([long(max_n_guess * 0.95),max_n_guess-10L])
	     ;if n_guess le min_n_guess then n_guess= min([long(min_n_guess * 1.05),min_n_guess+10L])
	     if n_guess ge max_n_guess then n_guess= long(max_n_guess * 0.95)
	     if n_guess le min_n_guess then n_guess= long(min_n_guess * 1.05)

	     if n_guess eq last_n_guess then n_guess= long(0.5 * (max_n_guess + min_n_guess))

	     n_iterations= n_iterations+1
        endrep until ((abs(dm/mtot) lt 0.001) or (n_iterations gt 150))

        r_eff = r[n_guess]
        
        return, r_eff
        
end     




