;======================================================================

pro sfr_at_snaptimes, frun, h=h


if not keyword_set(frun) then begin
   print, "  "
   print, "sfr_at_snaptimes, frun, /h"
   print, "  "
   print, "  "
   return
endif


;-------------------------------------------
;-------------------------------------------


; determine the number of
; snapshots in frun directory
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snapshot* | wc ",result
nsnaps=long(result[0])




; ---------------------------------------------------


time= fltarr(nsnaps)

sfr_inst= fltarr(nsnaps)
sfr_avg= fltarr(nsnaps)
sfr_snap= fltarr(nsnaps)


ymin = 0


; physical units
if keyword_set(h) then begin
        h = fload_cosmology('h')
endif

; ---------------------------------------------------


        ;-------------------------------------------
        ;   Get SFR rate from txt - for each file
        ;-------------------------------------------
        read_file_sfr_processed, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac



            ; physical units
            if keyword_set(h) then sfrtime = sfrtime / h
            ;if i eq 1 then sfrtime = sfrtime/(0.7)    ; did this for comparing h and noh


; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"


        ; open snapshot
        ;ok=fload_snapshot(frun,snapnum)
        ;ok=fload_snapshot_bh(frun,i,/nopot_in_snap)
        ok=fload_snapshot_bh(frun,i,/skip_center)


        sfr_snap[i]= total(fload_gas_sfr(1))

        ; what time is it?
        time[i]= fload_time(1)

        idx_gtr_snaptime= where(sfrtime ge (time[i]-0.0001))
        idx_curr_snaptime= idx_gtr_snaptime[0]

        ; instantaneous sfr
        sfr_inst[i]= -1
        if idx_curr_snaptime(0) ne -1 then sfr_inst[i]= sfrsfr(idx_curr_snaptime)


        ; set time window
        dt= 0.01
        idx_window= where((sfrtime ge (time[i]-dt)) and (sfrtime le (time[i]+dt)))
        sfr_avg[i]= -1
        if idx_window(0) ne -1 then sfr_avg[i]= mean(sfrsfr(idx_window))

        print, "T= ", time[i], sfr_inst[i], sfr_avg[i]

endfor



;--------------------------------------
;--------------------------------------


; SFR(snaptimes) FILE
; --------------------
openw, 1, frun+'/sfr_snaptimes.txt', ERROR=err

printf, 1, "#   sfr_snaptimes.txt"
printf, 1, "# "
printf, 1, "#           inst.       avg     from snap   "
printf, 1, "# time       sfr        sfr        sfr      "
printf, 1, "# (Gyr)    (Mo /Yr)   (Mo /Yr)   (Mo /Yr)   "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"    ",3(F8.3,"   "))', $
                time[i], sfr_inst[i], sfr_avg[i], sfr_snap[i]
endfor
close, 1




; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            "
print, "---------------------------------"
print, "---------------------------------"





end




;====================================================================


