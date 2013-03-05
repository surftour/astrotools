;======================================================================
;
;
;
;
;======================================================================
pro minor_open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass

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

                ;if strmid(frun,1,4) eq 'raid' then begin
                ;        spawn, '/bin/ls '+frun+'/*.sfr ',result
                ;endif else begin
                ;        case loopnum of
                ;          0: datadir='/data/tcox/sfr'
                ;          1: datadir='/home'
                ;          2: datadir='/raid2'
                ;          3: datadir='/data'
                ;          4: datadir='/data6'
                ;          5: datadir='/data7'
                ;          else: break
                ;        endcase

                        ;sfrfile= datadir+'/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'

                ;        nfrun= fload_getid(frun)
                ;        spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/*.sfr ',result
                ;endelse
		;
                ;sfrfile=strcompress(result[0],/remove_all)
		;

		sfrfile= '/home2/tcox/upsand/sfr/'+fload_getid(frun)+'.sfr'



                openr, 1, sfrfile, ERROR=err
                close, 1

                if (err NE 0) then begin
                   ERR= 0
                   foundit= 0
                endif else begin
                   foundit= 1
                endelse

                loopnum= loopnum+1

            endrep until (foundit eq 1)

            if foundit eq 0 then begin
                msg= 'Problem: sfr file no found  file=[data]/tcox/'+fload_getid(frun)+'/'+fload_getid(frun)+'.sfr'
                print, "  "
                print, msg
                print, "  "
                return
            endif

            ; read to open
            print, "opening: ",sfrfile
            sfrdata= read_ascii(sfrfile)
            sfrtime= sfrdata.field1[0,*]
            sfrsfr= sfrdata.field1[1,*]
            sfrmfs= sfrdata.field1[3,*]    ; amount of new stars
            n_cols= n_elements(sfrdata.field1[*,0])
            n_rows= n_elements(sfrdata.field1[0,*])
            finalnewstarmass= sfrmfs[n_rows-1]
            ; gas mass
            if n_cols ge 6 then sfrgasmass= sfrdata.field1[6,*] else sfrgasmass=[20.0,20.0-finalnewstarmass]
            sfrmsg= ''


            t1per= 0
            idx=where(sfrtime ge 1.0)
            if idx(0) ne -1 and n_cols gt 5 then begin
                gm= sfrgasmass[0]-sfrgasmass(idx(0))
                t1per= 100.0*gm/sfrgasmass[0]
            endif

            t4per= 0
            idx=where(sfrtime ge 4.0)
            if idx(0) ne -1 and n_cols gt 5 then begin
                gm= sfrgasmass[0]-sfrgasmass(idx(0))
                t4per= 100.0*gm/sfrgasmass[0]
            endif

            maxsfr= max(sfrsfr)
            idx= where(sfrsfr eq maxsfr)
            maxsfrtime= sfrtime(idx)

            sfraftertmerg= 0.0
            smaftertmerg= 0.0
            taftertmerg= 0.0
            tmerger= 0.0
            if n_elements(tmerg) gt 1 then begin
                idx= where(sfrtime ge tmerg[i]+0.2)
                if idx(0) eq -1 then begin
                        sfraftertmerg= sfrsfr[n_rows-1]
                        smaftertmerg= sfrmfs[n_rows-1]
                        tmerger= sfrtime[n_rows-1]
                        taftertmerg= tmerger
                endif else begin
                        sfraftertmerg= sfrsfr[idx(0)]
                        smaftertmerg= sfrmfs[idx(0)]
                        tmerger= tmerg[i]
                        taftertmerg= sfrtime[idx(0)]
                endelse
            endif

            print, "-------------------------------------"
            print, "maximum sfr= ", maxsfr, " occurs at ", maxsfrtime
            print, "average sfr= ", total(sfrsfr)/n_rows
            print, "  "
            print, "      tmerg= ",tmerger,'  and 200 Myr later= ',taftertmerg
            print, "        sfr= ", sfraftertmerg,'  200 Myr after t_merger'
            print, "         sm= ", smaftertmerg,'  200 Myr after t_merger'
            n_mfs= n_elements(sfrgasmass)
            print, "-------------------------------------"
            print, "original gas mass   =", sfrgasmass[0]
            print, "remnant gas mass    =", sfrgasmass[n_mfs-1]
            print, "gas consumed        =", sfrgasmass[0]-sfrgasmass[n_mfs-1]
            print,"                      ", 100.0*(sfrgasmass[0]-sfrgasmass[n_mfs-1])/sfrgasmass[0],'  %'
            print,"                      ", t1per,' % after 1 Gyr'
            print,"                      ", t4per,' % after 4 Gyr'
            print, "-------------------------------------"

end






