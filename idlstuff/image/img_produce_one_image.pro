pro img_produce_one_image, frun, snapnum, xlen, sendto, $
			arepo=arepo, $
			center=center, $
			colorbar=colorbar, $
			crude=crude, $
			evensampling=evensampling, $
			filename=filename, $
			fitstoo=fitstoo, $
			imgcenter=imgcenter, $
			loadedsnap=loadedsnap, $
			msg=msg, $
			nolabels=nolabels, $
			particlesonly=particlesonly, $
			pixels=pixels, $
			pubstyle=pubstyle, $
			plotpts=plotpts, $
			plottype=plottype, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			showbhs=showbhs, $
			startid=startid, numpart=numpart, $
			thumbnail=thumbnail, $
			track_to_draw= track_to_draw, $
			xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
			xz=xz, yz=yz, $
			ReturnImage=ReturnImage, $
			ReturnTimeLbl= ReturnTimeLbl

if not keyword_set(xlen) then xlen=100.0
if not keyword_set(sendto) then sendto='ps'
if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(filename) then filename='contour.eps'
if not keyword_set(frun) then begin
   print, "  "
   print, " PROBLEM: img_produce_one_image"
   print, "  "
   print, "  "
   return
endif


center_orig= center
if not keyword_set(imgcenter) then imgcenter=''

if not keyword_set(plottype) then plottype="allstars"

; --------------------------------
;  Load variables for smoothing
; --------------------------------
    if not keyword_set(loadedsnap) then begin
      ;if (fload_snapshot(frun, snapnum)) then begin
      ;if (fload_snapshot_bh(frun, snapnum)) then begin
      ok= -1
      if keyword_set(arepo) then ok= fload_snapshot_bh(frun, snapnum, /nopot_in_snap, /arepo)
      if ok lt 0 then ok= fload_snapshot_bh(frun, snapnum, /nopot_in_snap)
      if ok eq 1 then begin
        print, "PROBLEM: opening file"
        return
      endif
    endif





;defaults

set_maxden= 10.0
;set_maxden= 1.0    ; 10^4 msolar/pc2
;set_dynrng= 1.0e+3
;set_dynrng= 1.0e+4
set_dynrng= 1.0e+5
;set_dynrng= 1.0e+6
;set_maxden= 0.5e-1
;set_dynrng= 1.0e+6





;---------------------------------
;  Gas Particles
;---------------------------------


if plottype eq "gas" then begin

    print, "OK, plotting gas"

    npart= fload_npart(0)
    N= long(npart)

    ; default
    ;x= (1./0.7)*fload_gas_xyz('x',center=[0,0,0])
    ;y= (1./0.7)*fload_gas_xyz('y',center=[0,0,0])
    ;z= (1./0.7)*fload_gas_xyz('z',center=[0,0,0])
    ;m= (1./0.7)*fload_gas_mass(1)
    x= fload_gas_xyz('x',center=[0,0,0])
    y= fload_gas_xyz('y',center=[0,0,0])
    z= fload_gas_xyz('z',center=[0,0,0])
    m= fload_gas_mass(1)
    hsml= fload_gas_hsml(1)
    ;hsml= 0.2 + m*0.0

    ;stop
    ;if keyword_set(arepo) then hsml=hsml/20000.
    ;if keyword_set(arepo) then hsml=hsml/5000.
    if keyword_set(arepo) then begin
	hsml= 2. * (fload_gas_volume(1)^(1./3.))
    endif

    ;; 1 galaxy
    ;if keyword_set(numpart) and keyword_set(startid) then begin
    ;	x= fload_1gal_gas_xyz('x',startid,numpart,center=[0,0,0])
    ;	y= fload_1gal_gas_xyz('y',startid,numpart,center=[0,0,0])
    ;	z= fload_1gal_gas_xyz('z',startid,numpart,center=[0,0,0])
    ;	m= fload_1gal_gas_mass(startid,numpart) - fload_1gal_gas_mfs(startid,numpart)
    ;endif

    ;xz= 1


    ; standard
    ; ---------
    ;set_maxden= 4.0
    set_maxden= 1.0
    ;set_maxden= 0.5
    ;set_maxden= 3.0e-3

    set_dynrng= 1.0e+7
    ;set_dynrng= 1.0e+6
    ;set_dynrng= 1.0e+4
    ;set_dynrng= 1.0e+3

    ; zoom into center
    ; -----------------
    ;set_maxden= 1.0e-1
    ;set_dynrng= 1.0e+3


endif





if plottype eq "gastemp" then begin

    ; not currently operational, check contour_gastemp (in old_crap), and will need to
    ; update both img_makepic and img_makepic_raw - but it shouldn't be too hard

    print, "OK, plotting gas"

    npart= fload_npart(0)
    N= long(npart)
    
    x= fload_gas_xyz('x',center=[0,0,0])
    y= fload_gas_xyz('y',center=[0,0,0])
    z= fload_gas_xyz('z',center=[0,0,0])
    m= fload_gas_mass(1)
    hsml= fload_gas_hsml(1)
    
    if keyword_set(arepo) then begin
        hsml= 2. * (fload_gas_volume(1)^(1./3.))
    endif

    ;temp= fload_gas_u(1)/0.012381322
    temp= fload_gas_temperature(1,/K)

    ; 


endif



if plottype eq "gasmetallicity" then  begin

	; check in old_crap/contour_metals.pro 

endif




if plottype eq "xrays" then begin

	; check in old_crap/contour_xrays.pro

endif



if plottype eq "sfr" then begin

    x= fload_gas_xyz('x',center=[0,0,0])
    y= fload_gas_xyz('y',center=[0,0,0])
    z= fload_gas_xyz('z',center=[0,0,0])
    m= fload_gas_sfr(1)
    hsml= fload_gas_hsml(1)

endif




if plottype eq "HI" then begin

    set_maxden= 1.0e-1
    set_dynrng= 1.0e+3

    x= fload_gas_xyz('x',center=[0,0,0])
    y= fload_gas_xyz('y',center=[0,0,0])
    z= fload_gas_xyz('z',center=[0,0,0])
    mtot= fload_gas_mass(1)
    print, "total gas mass= ", total(mtot)
    m= fload_gas_mass(1,/HI)
    print, "HI mass fraction= ", total(m)/total(mtot)
    hsml= fload_gas_hsml(1)

    ;
    ;
    ;   NOW, we did the HI cuts in fload_gas_mass (or at least
    ;        some proxy for it), but I'll now get rid of zero mass
    ;        particles, just to make things easier.
    ;       
    ;
    idx= where(m gt 0.0)
    x= x(idx)
    y= y(idx)
    z= z(idx)
    m= m(idx)
    hsml= hsml(idx)

endif





;--------------------------------
;  Dark Halo
;--------------------------------

if plottype eq "dm" or plottype eq "halo" then begin

    print, "OK, plotting dark halo"

    ndm= fload_npart(1)                         ; dm particles
    npart= ndm
    print, "Ntot= ",npart
    N= long(npart)

    x= fload_halo_xyz('x',center=[0,0,0])
    y= fload_halo_xyz('y',center=[0,0,0])
    z= fload_halo_xyz('z',center=[0,0,0])
    m= fload_halo_mass(1)

    hsml= x*0.0 + 2.0

endif





;--------------------------------
;  Stellar Particles
;--------------------------------

if plottype eq "allstars" then begin

    print, "OK, plotting stars"

    ndisk= fload_npart(2)                               ; disk particles
    nbulge= fload_npart(3)                              ; bulge particles
    nstars= fload_npart(4)
    npart= long(ndisk) + long(nbulge) + long(nstars)
    print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars
    
    x= fload_allstars_xyz('x',center=[0,0,0])
    y= fload_allstars_xyz('y',center=[0,0,0])
    z= fload_allstars_xyz('z',center=[0,0,0])
    m= fload_allstars_mass(1)
    ;hsml= fload_allstars_hsml(1)
    hsml= 0.2 + m*0.0

endif


if plottype eq "olstars" then begin

    ndisk= fload_npart(2)                               ; disk particles
    nbulge= fload_npart(3)                              ; bulge particles
    npart= ndisk + nbulge
    print, "Ntot,Ndisk,Nbulge= ",npart,ndisk,nbulge

    if ndisk gt 0 then begin
        x= fload_disk_xyz('x',center=[0,0,0])
        y= fload_disk_xyz('y',center=[0,0,0])
        z= fload_disk_xyz('z',center=[0,0,0])
        m= fload_disk_mass(1)
        N= long(ndisk)
    endif

    if nbulge gt 0 then begin
        if n_elements(x) gt 0 then begin
          x= [x,fload_bulge_xyz('x',center=[0,0,0])]
          y= [y,fload_bulge_xyz('y',center=[0,0,0])]
          z= [z,fload_bulge_xyz('z',center=[0,0,0])]
          m= [m,fload_bulge_mass(1)]
          N= long(ndisk+nbulge)
        endif else begin
          x= [fload_bulge_xyz('x',center=[0,0,0])]
          y= [fload_bulge_xyz('y',center=[0,0,0])]
          z= [fload_bulge_xyz('z',center=[0,0,0])]
          m= [fload_bulge_mass(1)]
          N= long(nbulge)
        endelse
    endif

    hsml= 0.2 + m*0.0

endif




;--------------------------------
;   New Stellar Particles
;--------------------------------

if plottype eq "newstars" then begin

	print, "plotting new stars only"

        nnewstars = fload_npart(4)                         ;newstars  particles
	npart= nnewstars

        print, "Ntot,Nnewstars= ",npart,nnewstars

        x= fload_newstars_xyz('x',center=[0,0,0])
        y= fload_newstars_xyz('y',center=[0,0,0])
        z= fload_newstars_xyz('z',center=[0,0,0])
        m= fload_newstars_mass(1)

        ;hsml= fload_newstars_hsml(1)
        hsml= 0.1 + 0.0*x

endif


if strmid(plottype,0,10) eq "youngstars" then begin

	;if youngstarsage ge 1.0 then youngstarsage= 0.1
	;youngstarsage= float(strmid(youngstars,12))
	youngstarsage= 0.1

        nnewstars = fload_npart(4)                         ;newstars  particles
        npart= nnewstars

        print, "Ntot,Nnewstars= ",npart,nnewstars

        x= fload_newstars_xyz('x',center=[0,0,0])
        y= fload_newstars_xyz('y',center=[0,0,0])
        z= fload_newstars_xyz('z',center=[0,0,0])
        m= fload_newstars_mass(1)

        ;hsml= fload_newstars_hsml(1)
        hsml= 0.1 + 0.0*x

	; what time is it?
        Ttime= float(fload_time(1))
        nsage= float(Ttime-fload_newstars_age(1))

	print, " "
	print, " selecting stars younger than: ", youngstarsage
	print, " "

        idx=where(nsage lt youngstarsage)
        if idx(0) ne -1 then begin
            N= long(n_elements(idx))
            x= x(idx)
            y= y(idx)
            z= z(idx)
            m= m(idx)
            hsml= hsml(idx)
        endif else begin
            print, " "
            print, " PROBLEM: no young stars"
            print, " "
            return
        endelse

endif






;===============================================================================




; determine_center
; -----------------

; old method
;if keyword_set(center) then goto, donewithcenters
;center= [0,0,0]
;center= fload_center_alreadycomp(1)

if total(abs(center_orig)) gt 0.0 then begin
	center=center_orig
	print, "center= ", center
	goto, donewithcenters
endif



if imgcenter eq "bh1" or imgcenter eq "bh2" then begin

	;center= [0,0,0]
	center_bh= fload_center_alreadycomp(1)

	; two black holes
	if fload_npart(5) gt 1 then begin
	        bhid= fload_blackhole_id(1)
	        bhid1= min(bhid)
	        bhid2= max(bhid)
	        center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
	        print, "Blackhole : ", bhid1, center1
	        center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
	        print, "Blackhole : ", bhid2, center2
	        idx= where(bhid eq min(bhid))
		if keyword_set(centerbh1) then center_bh= center1
		if keyword_set(centerbh2) then center_bh= center2
	endif

	; one black hole
	if fload_npart(5) eq 1 then begin
	        bhid= fload_blackhole_id(1)
	        center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
	        print, "Blackhole : ", bhid, center_bh
	endif

	center= center_bh
	print, "using this center= ", center

endif


if imgcenter eq "zero" then center= fload_center_alreadycomp(1)



donewithcenters:

center_used= center
print, "center_used= ", center_used



; ----------------------------------------------
;  Call our generic contour plotting procedure
; ----------------------------------------------



img_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
	hsml=hsml, $
	filename=filename, fitstoo=fitstoo, $
	xthickness=xthickness, ythickness=ythickness, $
	pixels=pixels, zthickness=zthickness, $
	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	crude=crude, center=center, $
	set_maxden= set_maxden, $
        set_dynrng= set_dynrng, $
	NxNImage=NxNImage, NxNTempImage=NxNTempImage

	ReturnImage= NxNImage




	ReturnTimeLbl= fload_timelbl(1,2,/noteq)
	;ReturnTimeLbl= fload_timelbl(1,3,/noteq)
	;ReturnTimeLbl= fload_timelbl(0.7,3,/noteq)
	;ReturnTimeLbl= fload_timelbl(0.7,2,/noteq)

	center= center_used
; -------------
;  Done
; -------------



end









