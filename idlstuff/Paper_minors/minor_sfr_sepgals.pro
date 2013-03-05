
;======================================================================








pro sfr_at_snaptimes, frun, $
		filename=filename, $
		h=h


if not keyword_set(frun) then begin
   print, "  "
   print, "sfr_at_snaptimes, frun, filename=filename, /h"
   print, "  "
   print, "  "
   return
endif


if not keyword_set(filename) then filename="sfr_snaptimes.txt"



;-------------------------------------------
;-------------------------------------------


; determine the number of
; snapshots in frun directory
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])




; ---------------------------------------------------


time= fltarr(nsnaps)

sfr= fltarr(nsnaps)
sfr_1= fltarr(nsnaps)
sfr_2= fltarr(nsnaps)




; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"


        ; open snapshot
        ok=fload_snapshot(frun,i)


        ; what time is it?
        time[i]= fload_time(1)

	; instantaneous sfr for galaxy 1 and 2
	;---------------------------------------
	; this is for G3G1
	;sfr[i]= total(fload_gas_sfr(1))
	;sfr_1[i]= total(fload_1gal_gas_sfr(1,1,240000))
	;sfr_2[i]= total(fload_1gal_gas_sfr(1,240001,95000))
	; this is for SbIm
	sfr[i]= total(fload_gas_sfr(1))
	sfr_1[i]= total(fload_1gal_gas_sfr(1,1,533332))
	sfr_2[i]= total(fload_1gal_gas_sfr(1,533333,56666))
	; this is for SbG1
	;sfr[i]= total(fload_gas_sfr(1))
	;sfr_1[i]= total(fload_1gal_gas_sfr(1,1,533332))
	;sfr_2[i]= total(fload_1gal_gas_sfr(1,533333,95000))


endfor



;--------------------------------------
;--------------------------------------


; SFR(snaptimes) FILE
; --------------------
openw, 1, frun+'/sfr_snaptimes.txt', ERROR=err

printf, 1, "#   sfr_snaptimes.txt"
printf, 1, "# "
printf, 1, "#           Total      Gal 1      Gal 2 "
printf, 1, "# time       sfr        sfr        sfr  "
printf, 1, "# (Gyr)    (Mo /Yr)   (Mo /Yr)   (Mo /Yr)   "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"    ",3(F8.3,"   "))', $
                time[i], sfr[i], sfr_1[i], sfr_2[i]
endfor
close, 1




; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            "
print, "---------------------------------"
print, "---------------------------------"





end








;======================================================================


; ----------------------------
;  Read the above file
; ----------------------------
pro read_sfratsnapt_file, frun, time, sfrtot, sfr1, sfr2

filename= frun+'/sfr_snaptimes.txt'

spawn, "wc "+filename,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then filedata= fltarr(4,lines)

openr, 1, filename

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, filedata
close, 1


time= filedata[0,*]
sfrtot= filedata[1,*]
sfr1= filedata[2,*]
sfr2= filedata[3,*]


end









