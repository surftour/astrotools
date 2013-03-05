;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Load a initial condition into IDL                    ;;
;; -takes as inputes snapshot directory and number      ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function fload_initial_conditions, frun, h=h

    COMMON GalaxyHeader, N,npart,massarr,time,redshift,flag_sfr,flag_feedbacktp,npartTotal, $
		flag_cooling,numfiles,cosmocrap, $
		flag_multiphase, $
		flag_stellarage, $
		flag_snaphaspot,$
		flag_metals, la, $
		flag_stargens, flag_energydetail, flag_parentid, flag_starorig
    COMMON HaloData, xhalo,yhalo,zhalo,vxhalo,vyhalo,vzhalo, mhalo
    COMMON DiskData, xdisk,ydisk,zdisk,vxdisk,vydisk,vzdisk, mdisk
    COMMON OtherData, id, mass
    COMMON GasData, xgas,ygas,zgas,vxgas,vygas,vzgas,u,rho,hsml, mgas
    COMMON OldStarData, diskage, diskmetals, bulgeage, bulgemetals
    COMMON SfrData, mfs, sfr, stellage, gmetals, smetals
    COMMON MultiphaseData, mclouds
    COMMON CoolingData, nume, numh
    COMMON FeedbackData, tpu
    COMMON BulgeData, xbulge,ybulge,zbulge,vxbulge,vybulge,vzbulge, mbulge
    COMMON NewStarData, xstars,ystars,zstars,vxstars,vystars,vzstars, mstars
    COMMON FileInfo, rundir, runnum
    COMMON PotData, pgas, phalo, pdisk, pbulge, pstars, pbh
    COMMON EnergyDetail, totradiated, totshocked, totfeedback
    COMMON ParentID, parentid
    COMMON StarOrig, origmass, orighsml
    COMMON Center, com
    COMMON BlackHoleData, xbh, ybh, zbh, vxbh, vybh, vzbh, mbh
;    COMMON Galaxy1, startid1, numpart1, gas_startid1, gas_numpart1
;    COMMON Galaxy2, startid2, numpart2, gas_startid2, gas_numpart2

    if not keyword_set(frun) then begin
	print, "  "
	print, "fload_initial_conditions, frun"
	return, 1
    endif

    ;spawn, "/usr/bin/ls "+frun, result
    ;if strlen(result) eq 0 then return, 1

    ;fname=result
    ;fname=strcompress(fname,/remove_all)
    ;rundir=strcompress(frun,/remove_all)
    ;runnum=num
    fname=strcompress(frun,/remove_all)
    rundir="junk"

    if keyword_set(h) then hubble= fload_cosmology('h')
    
    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedbacktp=0L
    npartTotal=lonarr(6)	
    flag_cooling=0L
    numfiles=0L
    cosmocrap=dblarr(4)
    flag_stellarage=0L
    flag_metals=0L
    bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4 - 2*4 - 4*8 - 2*4
    la=intarr(bytesleft/2)


    openr,1,fname,/f77_unformatted,ERROR=err
    if (err NE 0) then begin
	print, "  "
	print, !ERR_STRING
	print, "  "
	close, 1
	return, 1
    endif else begin
	print, "opening: ",fname
    endelse
    readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedbacktp,npartTotal,flag_cooling,numfiles,cosmocrap,flag_stellarage,flag_metals,la


    ; volker isn't using a flag, so do it manually
    ;flag_snaphaspot= 1L
    flag_snaphaspot= 0L
    flag_stargens= 1L
    flag_energydetail= 0L
    flag_parentid= 0L
    flag_starorig= 0L
    
   
    N=total(npart)
    pos=fltarr(3,N)
    vel=fltarr(3,N)
    id=lonarr(N)
 
    ind=where((npart gt 0) and (massarr eq 0)) 
    if ind(0) ne -1 then begin
	Nwithmass= total(npart(ind))
        mass=fltarr(Nwithmass)
    endif else begin	
        Nwithmass= 0
    endelse

    readu,1,pos
    readu,1,vel
    readu,1,id
    if Nwithmass gt 0 then begin
      readu,1,mass
      print, "Nwithmasses= ",Nwithmass
    endif

    Ngas=npart(0)
    Nhalo=npart(1)
    Ndisk=npart(2)
    Nbulge=npart(3)
    Nstars=npart(4)
    Nbh=npart(5)

    if Ngas gt 0 then begin
        u=fltarr(Ngas)
        readu,1,u
    endif

    if Ngas gt 0 then rho=fltarr(Ngas)

    if flag_cooling gt 0 and Ngas gt 0 then begin
	nume=fltarr(Ngas)
	numh=fltarr(Ngas)
    endif

    if Ngas gt 0 then hsml=fltarr(Ngas)
    if flag_sfr gt 0 and Ngas gt 0 then sfr=fltarr(Ngas)

	if flag_sfr gt 0 then begin
	    if flag_stellarage gt 0 and Nstars gt 0 then begin
		    stellage= fltarr(Nstars)
	    endif
	endif

    if flag_metals gt 0 then begin
	if (Ngas+Nstars) gt 0 then begin
		mets= fltarr(Ngas+Nstars)
		gmetals= fltarr(Ngas)
		if Nstars gt 0 then smetals= fltarr(Nstars)
		gmetals(*)= mets(0:Ngas-1)
		if Nstars gt 0 then smetals(*)= mets(Ngas:Ngas+Nstars-1)
	endif
    endif

    if flag_snaphaspot gt 0 then begin
        pot= fltarr(N)
    endif

    close,1


    if keyword_set(h) then begin
	hubble= fload_cosmology('h')
	print, "snapshot values corrected for h= ", hubble
	
	time= time/hubble
        pos= pos/hubble
        mass= mass/hubble
	massarr= massarr/hubble
        rho=rho*hubble*hubble
        mfs=mfs/hubble
        if flag_multiphase gt 0 then mclouds=mclouds/hubble
	hsml=hsml/hubble
	if flag_stellarage gt 0 and Nstars gt 0 then stellage= stellage/hubble
	if flag_metals gt 0 then gmetals=gmetals/hubble
	if flag_metals gt 0 and Nstars gt 0 then smetals=smetals/hubble
	if flag_starorig gt 0 and Nstars gt 0 then origmass=origmass/hubble
	if flag_starorig gt 0 and Nstars gt 0 then orighsml=orighsml/hubble
    endif


    ; ----------------------------------
    ;   now start parsing file
    ; ----------------------------------
    if Ngas gt 0 then begin
        xgas=fltarr(Ngas) &  ygas=fltarr(Ngas)  & zgas=fltarr(Ngas) & mgas=fltarr(Ngas)
        xgas(*)=pos(0,0:Ngas-1)
        ygas(*)=pos(1,0:Ngas-1)
        zgas(*)=pos(2,0:Ngas-1)

        vxgas=fltarr(Ngas) &  vygas=fltarr(Ngas)  & vzgas=fltarr(Ngas)
        vxgas(*)=vel(0,0:Ngas-1)
        vygas(*)=vel(1,0:Ngas-1)
        vzgas(*)=vel(2,0:Ngas-1)

        if massarr(0) eq 0 then begin
            mgas(*)=mass(0:Ngas-1)	
        endif else begin
            mgas(*)= massarr(0)
	endelse

	if flag_snaphaspot gt 0 then begin
	    pgas=fltarr(Ngas)
	    pgas(*)=pot(0:Ngas-1)
	endif
    endif

    if Nhalo gt 0 then begin
        xhalo=fltarr(Nhalo) &  yhalo=fltarr(Nhalo) & zhalo=fltarr(Nhalo) & mhalo=fltarr(Nhalo)
        xhalo(*)=pos(0,0+Ngas:Nhalo+Ngas-1)
        yhalo(*)=pos(1,0+Ngas:Nhalo+Ngas-1)
        zhalo(*)=pos(2,0+Ngas:Nhalo+Ngas-1)

        vxhalo=fltarr(Nhalo) &  vyhalo=fltarr(Nhalo) & vzhalo=fltarr(Nhalo)
        vxhalo(*)=vel(0,0+Ngas:Nhalo+Ngas-1)
        vyhalo(*)=vel(1,0+Ngas:Nhalo+Ngas-1)
        vzhalo(*)=vel(2,0+Ngas:Nhalo+Ngas-1)

        if massarr(1) eq 0 then begin
	    skip=0L
            for t=0,0 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mhalo(*)=mass(0+skip:Nhalo-1+skip)	
        endif else begin
            mhalo(*)= massarr(1)
	endelse

        if flag_snaphaspot gt 0 then begin
            phalo=fltarr(Nhalo)
            phalo(*)=pot(0+Ngas:Nhalo+Ngas-1)
        endif
    endif

    ; Disk
    ; ------
    if Ndisk gt 0 then begin
        xdisk=fltarr(Ndisk) &  ydisk=fltarr(Ndisk) &  zdisk=fltarr(Ndisk) & mdisk=fltarr(Ndisk)
        xdisk(*)=pos(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        ydisk(*)=pos(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        zdisk(*)=pos(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)

        vxdisk=fltarr(Ndisk) &  vydisk=fltarr(Ndisk) &  vzdisk=fltarr(Ndisk)
        vxdisk(*)=vel(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        vydisk(*)=vel(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        vzdisk(*)=vel(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)

        if massarr(2) eq 0 then begin
	    skip=0L
            for t=0,1 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mdisk(*)=mass(0+skip:Ndisk-1+skip)	
        endif else begin
            mdisk(*)= massarr(2)
	endelse

        if flag_snaphaspot gt 0 then begin
            pdisk=fltarr(Ndisk)
            pdisk(*)=pot(Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        endif

	diskage= -2.0*randomu(seed,Ndisk)
	diskmetals= 10^(1.5*randomu(seed,Ndisk) - 3.0)

    endif



    ; Bulge
    ; -------
    if Nbulge gt 0 then begin
        xbulge=fltarr(Nbulge) &  ybulge=fltarr(Nbulge) &  zbulge=fltarr(Nbulge) & mbulge=fltarr(Nbulge)
        xbulge(*)=pos(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        ybulge(*)=pos(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        zbulge(*)=pos(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)

        vxbulge=fltarr(Nbulge) &  vybulge=fltarr(Nbulge) &  vzbulge=fltarr(Nbulge)
        vxbulge(*)=vel(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        vybulge(*)=vel(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        vzbulge(*)=vel(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)

        if massarr(3) eq 0 then begin
	    skip=0L
            for t=0,2 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mbulge(*)=mass(0+skip:Nbulge-1+skip)	
        endif else begin
            mbulge(*)= massarr(3)
	endelse

        if flag_snaphaspot gt 0 then begin
            pbulge=fltarr(Nbulge)
            pbulge(*)=pot(Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        endif

	bulgeage= -2.0 + 0.0*mbulge
	bulgemetals= 0.001 + 0.0*mbulge
    endif


    if Nstars gt 0 then begin
        xstars=fltarr(Nstars) &  ystars=fltarr(Nstars)  & zstars=fltarr(Nstars)  & mstars=fltarr(Nstars)
        xstars(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        ystars(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        zstars(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)

        vxstars=fltarr(Nstars) &  vystars=fltarr(Nstars) &  vzstars=fltarr(Nstars)
        vxstars(*)=vel(0,Nhalo+Ndisk+Ngas+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        vystars(*)=vel(1,Nhalo+Ndisk+Ngas+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        vzstars(*)=vel(2,Nhalo+Ndisk+Ngas+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)

        if massarr(4) eq 0 then begin
	    skip=0L
            for t=0,3 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mstars(*)=mass(0+skip:Nstars-1+skip)	
        endif else begin
            mstars(*)= massarr(4)
	endelse

        if flag_snaphaspot gt 0 then begin
            pstars=fltarr(Nstars)
            pstars(*)=pot(Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        endif
    endif


    if Nbh gt 0 then begin
        xbh=fltarr(Nbh) &  ybh=fltarr(Nbh)  & zbh=fltarr(Nbh)  & mbh=fltarr(Nbh)
        xbh(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
        ybh(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
        zbh(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)

        vxbh=fltarr(Nbh) &  vybh=fltarr(Nbh) &  vzbh=fltarr(Nbh)
        vxbh(*)=vel(0,Nhalo+Ndisk+Ngas+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
        vybh(*)=vel(1,Nhalo+Ndisk+Ngas+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
        vzbh(*)=vel(2,Nhalo+Ndisk+Ngas+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)

        if massarr(5) eq 0 then begin
            skip=0L
            for t=0,4 do begin
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                        skip=skip + npart(t)
                endif
            endfor
            mbh(*)=mass(0+skip:Nbh-1+skip)
        endif else begin
            mbh(*)= massarr(5)
        endelse

        if flag_snaphaspot gt 0 then begin
            pbh=fltarr(Nbh)
            pbh(*)=pot(Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
        endif
    endif

    com=fload_center(1)

    return, 0

end


