; ---------------------------------------------------------------------------
;
;  Compile some BH properties.
;
;
; ---------------------------------------------------------------------------


;-----------------------------
;  Generate bh_mass file
;-----------------------------

; We write a separate program to do this
; because the snapshot stores the particle
; mass which is mass of the blackhole AND
; the surrounding gas.

pro time_bh_mass, frun


if not keyword_set(frun) then begin
	print, " "
	print, " generate_bh_mass_file_frombh, frun"
	print, " "
	return
endif

;-------------------------------------
;   Get BH Info from fid.bh file
;-------------------------------------
; get bh data
;bhfile= '/home/tcox/data/bh/'+fload_getid(frun)+'_wBH.bh'
;bhfile= '/home/tcox/data/bh/'+fload_getid(frun)+'.bh'
bhfile= '/raid4/tcox/'+fload_getid(frun)+'/'+fload_getid(frun)+'.bh'    ; harvard

spawn, "wc "+bhfile,result
lines=long(result)
lines=lines(0)
if lines gt 0 then bhdata= fltarr(5,lines)

; open file
openr, 1, bhfile, ERROR=err

if (err NE 0) then begin
        print, "  "
        print, "Problem: ",!ERR_STRING
        print, "  "
        close, 1
        ERR= 0
        bhtime= [0]
        bhmsg= "No Black Hole"
endif else begin
        ;close, unit
        print, "opening: ",bhfile
        bhdata= read_ascii(bhfile)
        bhtime= bhdata.field1[0,*]
        bh_num= bhdata.field1[1,*]
        bh_mass= bhdata.field1[2,*]
        bh_mdot_gu= bhdata.field1[3,*]
        bh_mdot_sunyr= bhdata.field1[4,*]
        bh_totalmass= bhdata.field1[5,*]
        bh_mdot_edd= bhdata.field1[6,*]
        ;readf,1,bhdata
        ;bhtime= bhdata[0,*]
        ;sfrsfr= bhdata[1,*]
        close, 1
endelse

idx= where(bh_num eq 1)
if idx(0) ne -1 then time_merge= bhtime(idx(0))


; determine the number of
; snapshots in frun directory
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])
;nsnaps=31   ; don't have snaps in directory, do by hand

time= fltarr(nsnaps)
bmass= fltarr(nsnaps)


; ----------------------------------------
; This part loops through the snapshots
; and compiles BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

	;print, "--------------------------------------"


	; open snapshot
	;ok=fload_snapshot(frun,snapnum)
	ok=fload_snapshot_bh(frun,i)


	; what time is it?
	time[i]= fload_time(1)
	;time[i]= 0.05*i


	; get black hole mass
	idx= where(bhtime ge time[i])
	bh_mass[i]= 0
	if idx(0) ne -1 then bmass[i]=bh_mass(idx(0))
	if idx(0) eq 0 then bmass[i]=bh_mass(idx(1))


endfor 




; print bh information to file
; -------------------------------
openw, 1, frun+'/bh_mass.txt', ERROR=err

printf, 1, "#   bh_mass.txt"
printf, 1, "#   "
printf, 1, "#   "
printf, 1, "# time    bh mass"
printf, 1, "#(Gyr)     (gadu)"
for i=0,nsnaps-1 do begin
	printf, 1, FORMAT= '(F6.3,"  ",F9.7)', time[i], bmass[i]
endfor
close, 1

end





;=====================================================================




