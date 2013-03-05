;======================================================================
;======================================================================

;======================================================================

;======================================================================
;======================================================================
pro read_file_sfr_raw, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac

	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	    ; get sfr data (from file)
	    ; -------------------------
	    loopnum= 0
	    repeat begin
	    ;sfrfile= '/data/tcox/sfr/'+fload_getid(fruns[i])+'.sfr'        ; twopiter
	    ;sfrfile= '/home/tcox/data/sfr/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/raid4/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/raid2/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/data6/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard

                if strmid(frun,1,4) eq 'raid' then begin
                        spawn, '/bin/ls '+frun+'/sfr.txt ',result
                endif else begin
			case loopnum of
                          0: datadir='/n/home/tcox/data'
                          1: datadir='/n/circelfs/hernquist_lab'
                          2: datadir='/n/home/tcox'
                          3: datadir='/n/home'
                          4: datadir='/n/scratch/hernquist_lab/tcox'
                          5: datadir='/data7'
                          6: datadir=''
                          else: break
			endcase

			;sfrfile= datadir+'/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'

			nfrun= fload_getid(frun)
			spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/sfr.txt ',result
		endelse

		sfrfile=strcompress(result[0],/remove_all)



		get_lun, unit
		openr, unit, sfrfile, ERROR=err
		close, unit
		free_lun, unit

		if (err NE 0) then begin
		   ERR= 0
		   foundit= 0
		endif else begin
		   foundit= 1
		endelse

		loopnum= loopnum+1

	    endrep until (foundit eq 1)

	    if foundit eq 0 then begin
		msg= 'Problem: sfr file no found  file=[data]/tcox/'+fload_getid(frun)+'/sfr.txt'
	        print, "  "
	        print, msg
	        print, "  "
		return
	    endif

	    ; read to open
	    print, "opening: ",sfrfile
	    sfrdata= read_ascii(sfrfile)
	    sfrtime= sfrdata.field1[0,*]
	    sfrmfs= sfrdata.field1[1,*]
	    sfrsfr_gadu= sfrdata.field1[2,*]    ; amount of new stars
	    sfrsfr= sfrdata.field1[3,*]
	    sfrmfs_spawned= sfrdata.field1[4,*]
	    n_cols= n_elements(sfrdata.field1[*,0])
	    n_rows= n_elements(sfrdata.field1[0,*])

            maxsfr= max(sfrsfr)
            idx= where(sfrsfr eq maxsfr)
            maxsfrtime= sfrtime(idx)

            print, "-------------------------------------"
            print, "maximum sfr= ", maxsfr, " occurs at ", maxsfrtime
            print, "average sfr= ", total(sfrsfr)/n_rows
            print, "-------------------------------------"

	    N_sample= 5
            print, "-------------------------------------"
            print, "-------------------------------------"
	    print, "   RESAMPLING N_sample= ", N_sample
	    sfrsfr= transpose(sfrsfr)
	    sfrtime= transpose(sfrtime)
	    n= n_elements(sfrsfr)
	    n_sub= n/N_sample
	    n_extra= n-N_sample*n_sub
	    for i=0,(N_sample-n_extra-1) do begin
		sfrsfr= [sfrsfr, sfrsfr[n-1]]
		sfrtime= [sfrtime, sfrtime[n-1]]
	    endfor
	    sfrsfr= rebin(sfrsfr, n_sub+1)
	    sfrtime= rebin(sfrtime, n_sub+1)
            print, "-------------------------------------"
            print, "-------------------------------------"

end





;--------------------------------------------------------------------



