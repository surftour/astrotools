Function fload_allstars_hsml, dummy, center=center

    COMMON GalaxyHeader
    COMMON FileInfo
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData
    COMMON OtherData
    COMMON Center

    if dummy NE '-1' then begin



	; do we already have a file of these
	; pre-calculated?
	;hsmlfile= rundir+'/allstars_hsml.txt'
	hsmlfile= rundir+'/allstars_hsml_'+exts+'.txt'
	openr, 1, hsmlfile, ERROR=err
	if (err EQ 0) then begin
		nstars= long(npart(2)+npart(3)+npart(4))
		hsml_data= fltarr(nstars)
		junk=''
		readf, 1, junk
		readf, 1, junk
		readf, 1, junk
		readf, 1, junk
		readf, 1, junk
		readf, 1, hsml_data
		close, 1
		Hsml= hsml_data
		return, Hsml
	endif else begin
		err= 0
		close, 1
	endelse
	


	if keyword_set(center) then begin
		c= center
	endif else begin
		c= com
	endelse

; -------------------------------------------------------------------
; now recalculate allxs, allys and allzs for returning correct x,y,z
; -------------------------------------------------------------------

        ; is there a disk?
        ; -----------------
        if npart(2) GT 0 then begin
	    allxs = [xdisk]
	    allys = [ydisk]
	    allzs = [zdisk]
	    allmass = [mdisk]
        endif

        ; is there a bulge?
        ; -----------------
        if npart(3) GT 0 then begin
	  if n_elements(allxs) gt 0 then begin
            allxs = [allxs, xbulge]
            allys = [allys, ybulge]
	    allzs = [allzs, zbulge]
	    allmass = [allmass, mbulge]
	  endif else begin
            allxs = [xbulge]
            allys = [ybulge]
            allzs = [zbulge]
            allmass = [mbulge]
	  endelse
        endif

        ; is there stars?
        ; -----------------
        if npart(4) GT 0 then begin
          if n_elements(allxs) gt 0 then begin
            allxs = [allxs, xstars]
            allys = [allys, ystars]
            allzs = [allzs, zstars]
            allmass = [allmass, mstars]
          endif else begin
            allxs = [xstars]
            allys = [ystars]
            allzs = [zstars]
            allmass = [mstars]
          endelse
        endif


	; Call C-program to compute
	; ---------------------------


	Nstars= long(n_elements(allxs))


	DesNgb=96L
	Hmax= 15.0

	if(Nstars gt 0) then begin

	        Coord=fltarr(3,Nstars)
	        Masses= fltarr(Nstars)

	        Coord(0,*)= allxs(*)
	        Coord(1,*)= allys(*)
	        Coord(2,*)= allzs(*)
	        Masses(*)= allmass(*)

	        Hsml= fltarr(Nstars)

		spawn, 'echo $TJHOME', result
		homedir= strcompress(result,/remove_all)
		libfile= homedir+'/Tools/C-Routines_for_IDL/StellarHsml/starhsml.so'
		S = CALL_EXTERNAL(libfile, $
	                'stellarhsml', $
	                Nstars, $
	                Coord, $
	                Masses, $
	                DesNgb, $
	                Hmax, $
	                Hsml, $
			/F_VALUE)

	endif

;stop

	; write a file of these values
	; -----------------------------
	openw, 1, hsmlfile, ERROR=err
	printf, 1, "#    allstars_hsml.txt, DesNgb= "+string(DesNgb)+"  Hmax= "+string(Hmax) 
	printf, 1, "# "
	printf, 1, "#  hsml (in gadget units kpc/h) for all stellar "
	printf, 1, "#  particles.  Should be in same order as snap."
	printf, 1, "# "
	for i=0L,Nstars-1 do printf, 1, FORMAT='(F9.5)', Hsml[i]
	close,1

	print, "Writing Hsml values to file:", hsmlfile


	return, Hsml

    endif

end


