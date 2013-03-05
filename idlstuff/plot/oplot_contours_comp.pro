;---------------------------------------------------------------
;
;
;
;---------------------------------------------------------------
pro oplot_contours_comp, comptype, x0, y0, x1, y1, $
			hmsb=hmsb, $
                        xlen=xlen, $
			pixels=pixels, $
			clr=clr, crude=crude, $
			center=center, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			HalfMassSB=HalfMassSB, $
                        fitellipse=fitellipse, $
			fitdiskyboxy=fitdiskyboxy, $
                        pa=pa, ellip=ellip, $
			fit_a=fit_a, fit_b=fit_b, $
			fit_ellip=fit_ellip, a4diva=a4diva


	if comptype eq 'gas' then begin
                ; default
                x= fload_gas_xyz('x',center=[0,0,0])
                y= fload_gas_xyz('y',center=[0,0,0])
                z= fload_gas_xyz('z',center=[0,0,0])
                m= fload_gas_mass(1)
                hsml= fload_gas_hsml(1)

		; sd, by default, is in 10^msun / kpc^2
		set_maxden= 1.0e-3
		set_dynrng= 1.0e+4

		; this is roughly 2.089 g / cm^2

		; or 1.25 x 10^24 protons cm^-2 

		set_maxden= 1.25e-2     ; or column density of 10^22 cm^-2
		set_dynrng= 1.0e+7
	endif


        if comptype eq 'dm' then begin
                ; default
                x= fload_halo_xyz('x',center=[0,0,0])
                y= fload_halo_xyz('y',center=[0,0,0])
                z= fload_halo_xyz('z',center=[0,0,0])
                m= fload_halo_mass(1)
                hsml= m*0.0 + 0.2
        endif

        if comptype eq 'special' then begin
                ; default
                x= fload_1gal_halo_xyz('x',200000L,500000L,center=[0,0,0])
                y= fload_1gal_halo_xyz('y',200000L,500000L,center=[0,0,0])
                z= fload_1gal_halo_xyz('z',200000L,500000L,center=[0,0,0])
                m= fload_1gal_halo_mass(200000L,500000L)
                hsml= m*0.0 + 4.0
        endif




        img_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, $
                        pixels=pixels, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        crude=crude, center=center, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage, HalfMassSB=HalfMassSB

        ContourMap= NxNImage

	oplot_contours, ContourMap, x0, y0, x1, y1, $
                        hmsb=hmsb, $
                        xlen=xlen, $
                        pixels=pixels, $
                        clr=clr, $
                        fitellipse=fitellipse, $
                        fitdiskyboxy=fitdiskyboxy, $
                        pa=pa, ellip=ellip, $
                        fit_a=fit_a, fit_b=fit_b, $
                        fit_ellip=fit_ellip, a4diva=a4diva


end












; =================================================================
; =================================================================

