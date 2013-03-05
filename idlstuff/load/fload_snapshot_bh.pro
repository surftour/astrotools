;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Load a snapshot into IDL                             ;;
;; -takes as inputes snapshot directory and number      ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                      ;;
;;                                                      ;;
;;                                                      ;;
;;  use as follows:                                     ;;
;;                                                      ;;
;;      ok= fload_snapshot_bh(frun,0,[....])            ;;
;;                                                      ;;
;;      frun= full path to data directory.              ;;
;;            it assumes that all snapshots are         ;;
;;            called snapshot_###                       ;;
;;                                                      ;;
;;      ##= snapshot number that you want to load       ;;
;;          (don't need preceeding 0's)                 ;;
;;                                                      ;;
;;      ....= options, see below, but likely need       ;;
;;            to ask TJ about these.                    ;;
;;                                                      ;;
;;                                                      ;;
;;                                                      ;;
;;                                                      ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro check_field, var, varlabel

	print, "checking: "+varlabel, "   ",size(var)
	idx=where(finite(var) eq 0, n_i)
	if idx(0) ne -1 then begin
	    print, varlabel+' has '+strcompress(string(n_i),/remove_all)+' bad number(s)'
	    ;if n_i lt 10 then begin
	    ;	print, idx
	    ;	var(idx)= 1.0e-9
	    ;endif
	endif

end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function fload_snapshot_bh, frun, num, $
			h0=h0, $
			ics=ics, $
			arepo=arepo, $
			nopot_in_snap=nopot_in_snap, $
			header_only=header_only, $
			show_header=show_header, $
			do_four=do_four, $
			skip_center=skip_center, $
			skip_extra_gas_info=skip_extra_gas_info

    COMMON GalaxyHeader, N,npart,massarr,time,redshift,flag_sfr,flag_feedbacktp,npartTotal, $
		flag_cooling,numfiles,cosmocrap, $
		flag_multiphase, $
		flag_stellarage, $
		flag_snaphaspot,$
		flag_metals, la, $
		flag_stargens, flag_energydetail, flag_parentid, flag_starorig
    COMMON HaloData, xhalo,yhalo,zhalo,vxhalo,vyhalo,vzhalo, mhalo,dmid
    COMMON DiskData, xdisk,ydisk,zdisk,vxdisk,vydisk,vzdisk, mdisk,did
    COMMON OtherData, id, mass
    COMMON GasData, xgas,ygas,zgas,vxgas,vygas,vzgas,u,rho,volume,hsml,mgas,gid
    COMMON OldStarData, diskage, diskmetals, bulgeage, bulgemetals
    COMMON SfrData, mfs, sfr, stellage, gmetals, smetals
    COMMON MultiphaseData, mclouds
    COMMON CoolingData, nume, numh
    COMMON FeedbackData, tpu
    COMMON BulgeData, xbulge,ybulge,zbulge,vxbulge,vybulge,vzbulge,mbulge,bib
    COMMON NewStarData, xstars,ystars,zstars,vxstars,vystars,vzstars,mstars,nsid
    COMMON FileInfo, rundir, runnum, fname, exts
    COMMON PotData, pgas, phalo, pdisk, pbulge, pstars, pbh
    COMMON EnergyDetail, totradiated, totshocked, totfeedback
    COMMON ParentID, parentid
    COMMON StarOrig, origmass, orighsml
    COMMON Center, com, alternate_com, b_com
    COMMON BlackHoleData, xbh, ybh, zbh, vxbh, vybh, vzbh, mbh, bhmass, bhaccrate, bhid
;    COMMON Galaxy1, startid1, numpart1, gas_startid1, gas_numpart1
;    COMMON Galaxy2, startid2, numpart2, gas_startid2, gas_numpart2
;    COMMON SnapInfo, frun

    if not keyword_set(num) then num=0
    if not keyword_set(frun) then begin
	print, "  "
	print, "fload_snapshot, frun, num"
	return, -1
    endif

    exts='0000'
    exts=exts+strcompress(string(num),/remove_all)

    if num ge 1000 then do_four= 1

    ; does the snapshot have four, or three,
    ; numbers after it
    ; ---------------------------------------
    ;do_four= 0
    ;do_four= 1
    ;if do_four eq 0 then begin
    if not keyword_set(do_four) then begin
	exts=strmid(exts,strlen(exts)-3,3)
    endif else begin
	exts=strmid(exts,strlen(exts)-4,4)
	print, "******************************"
	print, "******************************"
	print, "   fload_snapshot_bh "
	print, "          is  "
	print, "  momentarily changed"
	print, "******************************"
	print, "******************************"
    endelse

    if keyword_set(ics) then begin
	fname= frun
	result= frun
	if strlen(fname) gt 0 and strmid(fname,strlen(fname)-4,4) eq "hdf5" then goto, file_is_hdf5
	goto, foundfile
    endif

    basename='snap'

    ; HDF5 snapshot file?
    cmd= "/bin/ls "+frun+"/snapshot*_"+exts+".hdf5"
    spawn, cmd, result   &    print, "result= ", result
    if strlen(result) gt 0 and strmid(result,strlen(result)-4,4) eq "hdf5" then goto, file_is_hdf5

    ; regular snapshot file
    cmd= "/bin/ls "+frun+"/snapshot*_"+exts  &    print, "spawning cmd= ", cmd
    spawn, cmd, result   &    print, "result= ", result
    if strlen(result) eq 0 then return, -1

    fname=strcompress(result,/remove_all)
    rundir=strcompress(frun,/remove_all)
    runnum=num



;=========================================================================================
;=========================================================================================


foundfile:

    if keyword_set(h0) then hubble= 0.7   ; hubble= fload_cosmology('h')
    
    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedbacktp=0L
    npartTotal=lonarr(6)	
    flag_cooling=0L
    numfles=0L
    cosmocrap=dblarr(4)
    flag_stellarage=0L
    flag_metals=0L
    bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4 - 2*4 - 4*8 - 2*4
    la=intarr(bytesleft/2)


    if file_size(fname) le 0 then return, -1


    openr,1,fname,/f77_unformatted,ERROR=err
    if (err NE 0) then begin
	print, "  "
	print, !ERR_STRING
	print, "  "
	close, 1
	return, -1
    endif else begin
	print, "opening: ",fname
    endelse

    catch, error_status
    if error_status ne 0 then begin
	print, " "
	print, "Error: ", error_status
	print, "Error message: ", !ERROR_STATE.MSG
	close, 1
	catch, /cancel
	return, -1
    endif

    readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedbacktp,npartTotal,flag_cooling,numfiles,cosmocrap,flag_stellarage,flag_metals,la

    if keyword_set(show_header) then begin
	print, "** HEADER **"
	print, "npart= ", npart
	print, "massarr= ", massarr
	print, "time= ", time
	print, "flag_sfr", flag_sfr
	print, "flag_feedbacktp", flag_feedbacktp
	print, "npartTotal", npartTotal
	print, "flag_cooling", flag_cooling
	print, "numfiles", numfiles
	print, "cosmocrap", cosmocrap
	print, "flag_stellarage", flag_stellarage
	print, "flag_metals", flag_metals
	print, "strlen(la)", strlen(la)
	print, "** ******* **"
    endif

    if time lt 0.0 or time gt 1.0e4 then begin
	print, " "
	print, " PROBLEM with time read from header, returning .... "
	print, " "
	return, -1
    endif

    ; volker isn't using a flag, so do it manually
    if keyword_set(nopot_in_snap) then flag_snaphaspot= 0L else flag_snaphaspot= 1L
    ;flag_snaphaspot= 0L
    flag_stargens= 1L
    flag_energydetail= 0L
    flag_parentid= 0L
    flag_starorig= 0L
    

    if keyword_set(header_only) then begin
	close, 1
	return, 1
    endif
   
    N=long(total(npart,/double))    ; could also do   total(npart,/integer)
    ;N=long(npart(0)+npart(1)+npart(2)+npart(3)+npart(4)+npart(5))
    pos=fltarr(3,N)
    vel=fltarr(3,N)
    id=lonarr(N)
 
    ind=where((npart gt 0) and (massarr eq 0)) 
    if ind(0) ne -1 then begin
	Nwithmass= long(total(npart(ind),/double))
        mass=fltarr(Nwithmass)
        print, "Nwithmasses= ",Nwithmass
    endif else begin	
        Nwithmass= 0
    endelse

    readu,1,pos
    check_field, pos, 'pos'
    readu,1,vel
    check_field, vel, 'vel'
    readu,1,id
    check_field, id, 'id'
    if Nwithmass gt 0 then begin
      readu,1,mass
      check_field, mass, 'mass'
      ;print, "Nwithmasses= ",Nwithmass
    endif

    Ngas=npart(0)
    Nhalo=npart(1)
    Ndisk=npart(2)
    Nbulge=npart(3)
    Nstars=npart(4)
    Nbh=npart(5)

    N_baryons= Ngas + Ndisk + Nbulge + Nstars
    N_tot= total(npart)
    print, "Ntot= ", N_tot
    print, npart

    if Ngas gt 0 then begin
        u=fltarr(Ngas)
        readu,1,u
	check_field, u, 'u'
    endif

    if Ngas gt 0 then begin
        rho=fltarr(Ngas)
        if not keyword_set(ics) then readu,1,rho
        check_field, rho, 'rho'
    endif

    if (Ngas gt 0) and keyword_set(arepo) then begin
        volume=fltarr(Ngas)
        if not keyword_set(ics) then readu,1,volume
        check_field, volume, 'volume'
    endif

if keyword_set(skip_extra_gas_info) then goto, moveon

    if flag_cooling gt 0 and Ngas gt 0 then begin
	nume=fltarr(Ngas)
	numh=fltarr(Ngas)
	if not keyword_set(ics) then readu,1,nume
	check_field, nume, 'nume'
	if not keyword_set(ics) then readu,1,numh
	check_field, numh, 'numh'
    endif

    if Ngas gt 0 then begin
	hsml=fltarr(Ngas)
	if not keyword_set(ics) then readu,1,hsml
	check_field, hsml, 'hsml'
    endif

    if flag_sfr gt 0 and Ngas gt 0 then begin
        sfr=fltarr(Ngas)
        if not keyword_set(ics) then readu,1,sfr
	check_field, sfr, 'sfr'
	print, "sfr= ", total(sfr)
    endif

	if flag_sfr gt 0 then begin
	    if flag_stellarage gt 0 and Nstars gt 0 then begin
		    stellage= fltarr(Nstars)
		    if not keyword_set(ics) then readu,1,stellage
			check_field, stellage, 'stellage'
	    endif
	endif

    ;if eof(1) eq 1 then goto, moveon
    ;goto, moveon

    if flag_metals gt 0 then begin
	if (Ngas+Nstars) gt 0 then begin
		mets= fltarr(Ngas+Nstars)
		if not keyword_set(ics) then readu,1,mets
		check_field, mets, 'mets'
		gmetals= fltarr(Ngas)
		if Nstars gt 0 then smetals= fltarr(Nstars)
		gmetals(*)= mets(0:Ngas-1)
		if Nstars gt 0 then smetals(*)= mets(Ngas:Ngas+Nstars-1)
	endif
    endif

    if keyword_set(nopot_in_snap) then goto, moveon

    if flag_snaphaspot gt 0 then begin
        pot= fltarr(N)
        if not keyword_set(ics) then readu,1,pot
	check_field, pot, 'pot'
    endif

    if eof(1) eq 1 then goto, moveon

    if Nbh gt 0 then begin
        bhmass= fltarr(Nbh)
        readu, 1, bhmass
        check_field, bhmass, 'bhmass'
        bhaccrate= fltarr(Nbh)
        readu, 1, bhaccrate
        check_field, bhaccrate, 'bhaccrate'

	print, "bhmass= ", bhmass
	print, "bhaccrate= ", bhaccrate
    endif


moveon:
    close,1


    if keyword_set(h0) then begin
	;hubble= fload_cosmology('h0')
	hubble= 0.7
	print, "snapshot values corrected for h= ", hubble
	
	time= time/hubble
        pos= pos/hubble
        mass= mass/hubble
	massarr= massarr/hubble
        rho=rho*hubble*hubble
        volume=volume/hubble/hubble/hubble
        ;mfs=mfs/hubble
        ;if flag_multiphase gt 0 then mclouds=mclouds/hubble
	hsml=hsml/hubble
	if flag_stellarage gt 0 and Nstars gt 0 then stellage= stellage/hubble
	;if flag_metals gt 0 then gmetals=gmetals/hubble
	;if flag_metals gt 0 and Nstars gt 0 then smetals=smetals/hubble
	;if flag_starorig gt 0 and Nstars gt 0 then origmass=origmass/hubble
	;if flag_starorig gt 0 and Nstars gt 0 then orighsml=orighsml/hubble

	;
	; potential and velocity should be fine

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
    mdisk= 0
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

	diskage= -5.0*randomu(seed,Ndisk)
	;diskage= -1.0*randomu(seed,Ndisk)    ; changed this for Marijn
	diskmetals= 10^(1.5*randomu(seed,Ndisk) - 3.0)

    endif



    ; Bulge
    ; -------
    mbulge= 0
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

	bulgeage= -7.0 + 0.0*mbulge
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
	;id1= Nhalo+Ngas+Ndisk+Nbulge+Nstars
	;id2= Nhalo+Ngas+Ndisk+Nbulge+Nstars+Nbh-1
	;xbh(*)=pos(0,id1:id2)
	;ybh(*)=pos(1,id1:id2)
	;zbh(*)=pos(2,id1:id2)
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


    goto, calccenter




;=========================================================================================
;=========================================================================================

    file_is_hdf5:


    fname=strcompress(result,/remove_all)
    rundir=strcompress(frun,/remove_all)
    runnum=num

    if h5f_is_hdf5(fname) eq 0 then begin
	print, " "
	print, "Problem, not (or bad) HDF5 file: ", fname
	print, " "
	return, -1
    endif

    print, " "
    print, "Reading HDF5 file: ", fname
    print, " "

    ; open the file
    ; -------------------------
    ;print, h5_parse(fname)   ; prints a useful overview
    file_id= h5f_open(fname)


    ; open the header group
    ; -------------------------
    hdr_group_id= h5g_open(file_id,"Header")
    if hdr_group_id lt 0 then begin
	print, "Problem, bad header"
	return, -1
    endif

    npart= h5a_read(h5a_open_name(hdr_group_id,"NumPart_ThisFile"))
    massarr= h5a_read(h5a_open_name(hdr_group_id,"MassTable"))
    time= h5a_read(h5a_open_name(hdr_group_id,"Time"))
    redshift= h5a_read(h5a_open_name(hdr_group_id,"Redshift"))
    flag_sfr= h5a_read(h5a_open_name(hdr_group_id,"Flag_Sfr"))
    flag_feedbacktp= h5a_read(h5a_open_name(hdr_group_id,"Flag_Feedback"))
    npartTotal= h5a_read(h5a_open_name(hdr_group_id,"NumPart_Total"))
    flag_cooling= h5a_read(h5a_open_name(hdr_group_id,"Flag_Cooling"))
    numfles= h5a_read(h5a_open_name(hdr_group_id,"NumFilesPerSnapshot"))
    cosmocrap= fltarr(4)
    cosmocrap[0]= h5a_read(h5a_open_name(hdr_group_id,"BoxSize"))
    cosmocrap[1]= h5a_read(h5a_open_name(hdr_group_id,"Omega0"))
    cosmocrap[2]= h5a_read(h5a_open_name(hdr_group_id,"OmegaLambda"))
    cosmocrap[3]= h5a_read(h5a_open_name(hdr_group_id,"HubbleParam"))
    flag_stellarage= h5a_read(h5a_open_name(hdr_group_id,"Flag_StellarAge"))
    flag_metals= h5a_read(h5a_open_name(hdr_group_id,"Flag_Metals"))

    h5g_close, hdr_group_id


   ; process things a bit
   if time lt 0.0 or time gt 1.0e4 then begin
        print, " "
        print, " PROBLEM with time read from header, returning .... "
        print, " "
        return, -1
    endif

    ; volker isn't using a flag, so do these manually
    flag_stargens= 1L
    flag_energydetail= 0L
    flag_parentid= 0L
    flag_starorig= 0L


    if keyword_set(header_only) then begin
        close, 1
        return, 1
    endif

    N=long(total(npart,/double))    ; could also do   total(npart,/integer)

    ind=where((npart gt 0) and (massarr eq 0))
    if ind(0) ne -1 then begin
        Nwithmass= long(total(npart(ind),/double))
        mass=fltarr(Nwithmass)
        print, "Nwithmasses= ",Nwithmass
    endif else begin
        Nwithmass= 0
    endelse


    Ngas=npart(0)
    Nhalo=npart(1)
    Ndisk=npart(2)
    Nbulge=npart(3)
    Nstars=npart(4)
    Nbh=npart(5)

    N_baryons= Ngas + Ndisk + Nbulge + Nstars
    N_tot= total(npart)
    print, "Ntot= ", N_tot
    print, npart




    ; now get the particle data
    ; ------------------------------

    ;  Gas
    ; -----
    if npart(0) gt 0 then begin
	group_id= h5g_open(file_id,"PartType0")
	xyz= h5d_read(h5d_open(group_id,"Coordinates"))
	xgas= transpose(xyz[0,*])
	ygas= transpose(xyz[1,*])
	zgas= transpose(xyz[2,*])
	xyz= h5d_read(h5d_open(group_id,"Velocities"))
	vxgas= transpose(xyz[0,*])
	vygas= transpose(xyz[1,*])
	vzgas= transpose(xyz[2,*])
	if massarr[0] le 0 then mgas= h5d_read(h5d_open(group_id,"Masses")) else mgas= 0.0*xgas + massarr[0]
	u= h5d_read(h5d_open(group_id,"InternalEnergy"))
	rho= h5d_read(h5d_open(group_id,"Density"))
	if flag_metals gt 0 then gmetals= h5d_read(h5d_open(group_id,"Metallicity"))
	if flag_cooling gt 0 then nume= h5d_read(h5d_open(group_id,"ElectronAbundance"))
	if flag_cooling gt 0 then numh= h5d_read(h5d_open(group_id,"NeutralHydrogenAbundance"))
	if flag_sfr gt 0 then sfr= h5d_read(h5d_open(group_id,"StarFormationRate"))
	hsml= h5d_read(h5d_open(group_id,"SmoothingLength"))
	;pgas= h5d_read(h5d_open(group_id,"Potential"))
	gid= h5d_read(h5d_open(group_id,"ParticleIDs"))
	h5g_close, group_id
    endif


    ;  DM
    ; -----
    if npart(1) gt 0 then begin
	group_id= h5g_open(file_id,"PartType1")
	xyz= h5d_read(h5d_open(group_id,"Coordinates"))
	xhalo= transpose(xyz[0,*])
	yhalo= transpose(xyz[1,*])
	zhalo= transpose(xyz[2,*])
	xyz= h5d_read(h5d_open(group_id,"Velocities"))
	vxhalo= transpose(xyz[0,*])
	vyhalo= transpose(xyz[1,*])
	vzhalo= transpose(xyz[2,*])
	if massarr[1] le 0 then mhalo= h5d_read(h5d_open(group_id,"Masses")) else mhalo= 0.0*xhalo + massarr[1]
	;phalo= h5d_read(h5d_open(group_id,"Potential"))
	dmid= h5d_read(h5d_open(group_id,"ParticleIDs"))
	h5g_close, group_id
    endif


    ;  Disk
    ; -------
    if npart(2) gt 0 then begin
	group_id= h5g_open(file_id,"PartType2")
	xyz= h5d_read(h5d_open(group_id,"Coordinates"))
	xdisk= transpose(xyz[0,*])
	ydisk= transpose(xyz[1,*])
	zdisk= transpose(xyz[2,*])
	xyz= h5d_read(h5d_open(group_id,"Velocities"))
	vxdisk= transpose(xyz[0,*])
	vydisk= transpose(xyz[1,*])
	vzdisk= transpose(xyz[2,*])
	if massarr[2] le 0 then mdisk= h5d_read(h5d_open(group_id,"Masses")) else mdisk= 0.0*xdisk + massarr[2]
	;pdisk= h5d_read(h5d_open(group_id,"Potential"))
	did= h5d_read(h5d_open(group_id,"ParticleIDs"))

	diskage= -5.0*randomu(seed,Ndisk)
	;diskage= -1.0*randomu(seed,Ndisk)    ; changed this for Marijn
	diskmetals= 10^(1.5*randomu(seed,Ndisk) - 3.0)
	h5g_close, group_id
    endif


    ;  Bulge
    ; -------
    if npart(3) gt 0 then begin
	group_id= h5g_open(file_id,"PartType3")
	xyz= h5d_read(h5d_open(group_id,"Coordinates"))
	xbulge= transpose(xyz[0,*])
	ybulge= transpose(xyz[1,*])
	zbulge= transpose(xyz[2,*])
	xyz= h5d_read(h5d_open(group_id,"Velocities"))
	vxbulge= transpose(xyz[0,*])
	vybulge= transpose(xyz[1,*])
	vzbulge= transpose(xyz[2,*])
	if massarr[3] le 0 then mbulge= h5d_read(h5d_open(group_id,"Masses")) else mbulge= 0.0*xbulge + massarr[3]
	;pbulge= h5d_read(h5d_open(group_id,"Potential"))
	bid= h5d_read(h5d_open(group_id,"ParticleIDs"))

	bulgeage= -7.0 + 0.0*mbulge
	bulgemetals= 0.001 + 0.0*mbulge

	h5g_close, group_id
    endif


    ;  New Stars
    ; ------------
    if npart(4) gt 0 then begin
	group_id= h5g_open(file_id,"PartType4")
	xyz= h5d_read(h5d_open(group_id,"Coordinates"))
	xstars= transpose(xyz[0,*])
	ystars= transpose(xyz[1,*])
	zstars= transpose(xyz[2,*])
	xyz= h5d_read(h5d_open(group_id,"Velocities"))
	vxstars= transpose(xyz[0,*])
	vystars= transpose(xyz[1,*])
	vzstars= transpose(xyz[2,*])
	if massarr[4] le 0 then mstars= h5d_read(h5d_open(group_id,"Masses")) else mstars= 0.0*xstars + massarr[4]
	if flag_metals gt 0 then smetals= h5d_read(h5d_open(group_id,"Metallicity"))
	if flag_sfr gt 0 and flag_stellarage gt 0 then stellage= h5d_read(h5d_open(group_id,"StellarFormationTime"))
	;pstars= h5d_read(h5d_open(group_id,"Potential"))
	nsid= h5d_read(h5d_open(group_id,"ParticleIDs"))
	h5g_close, group_id
    endif


    ;  Black Holes
    ; --------------
    if npart(5) gt 0 then begin
        group_id= h5g_open(file_id,"PartType5")
        xyz= h5d_read(h5d_open(group_id,"Coordinates"))
        xbh= transpose(xyz[0,*])
        ybh= transpose(xyz[1,*])
        zbh= transpose(xyz[2,*])
        xyz= h5d_read(h5d_open(group_id,"Velocities"))
        vxbh= transpose(xyz[0,*])
        vybh= transpose(xyz[1,*])
        vzbh= transpose(xyz[2,*])
        if massarr[5] le 0 then mbh= h5d_read(h5d_open(group_id,"Masses")) else mbh= 0.0*xbh + massarr[5]
        ;pbh= h5d_read(h5d_open(group_id,"Potential"))
        bhid= h5d_read(h5d_open(group_id,"ParticleIDs"))
        bhmdot= h5d_read(h5d_open(group_id,"BH_Mdot"))
        bhmdot= h5d_read(h5d_open(group_id,"BH_Mass"))
        h5g_close, group_id
    endif



    ;
    ; a bit of data processing
    ;

    ; combine ids
    id= [-1]
    if n_elements(gid) gt 0 then id= [id, gid]
    if n_elements(dmid) gt 0 then id= [id, dmid]
    if n_elements(did) gt 0 then id= [id, did]
    if n_elements(bid) gt 0 then id= [id, bid]
    if n_elements(nsid) gt 0 then id= [id, nsid]
    if n_elements(bhid) gt 0 then id= [id, bhid]
    id= id[1:*]



    ; done, get out of here
    ; -------------------------
    h5f_close, file_id

    goto, calccenter

;=========================================================================================
;=========================================================================================


calccenter:

    if keyword_set(skip_center) then begin
        com=[0.0, 0.0, 0.0]
        alternate_com= com
        b_com= com
    endif else begin
        com=fload_center(1)
        alternate_com= fload_center_felix(1)
        if N_baryons gt 0 then b_com= fload_center_baryons(1) else b_com= com
    endelse

    return, 0




end


