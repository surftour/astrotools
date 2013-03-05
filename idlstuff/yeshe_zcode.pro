; ------------------------------------------------------------------------------------------
;
;
;  need to .run
;    1. time_hotgas
;    2. sfr_multi
;
;
pro zcode, junk


if not keyword_set(junk) then begin
	print, " "
	print, " zcode, junk"
	print, " "
	print, " "
	return
endif



;-------------------------------------------


frun="/raid4/tcox/vc3vc3e_2"


Ngb= 16


metaldir="nsmetals/"


;-------------------------------------------


; determine the number of
; snapshots in frun directory
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])



ztime= fltarr(nsnaps)

gas_sfr= fltarr(nsnaps)
gas_Ia= fltarr(nsnaps)
gas_II= fltarr(nsnaps)

newstar_time= fltarr(nsnaps)
newstar_mass= fltarr(nsnaps)
newstar_Ia= fltarr(nsnaps)
newstar_II= fltarr(nsnaps)
newstar_Ngbs= fltarr(nsnaps,Ngb)
newstar_metalstodistribute= fltarr(nsnaps)




; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

dt= 0.0

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"


        ; open snapshot
        ;ok=fload_snapshot(frun,snapnum)
        ;ok=fload_snapshot_bh(frun,i,/nopot_in_snap)
        ok=fload_snapshot_bh(frun,i)


        ; what time is it?
        Ttime= float(fload_time(1))
        ztime[i]= Ttime

	if i ne 0 then dt= Ttime-lastime
	lastime= Ttime


	gas_sfr= fload_gas_sfr(1)
	gas_II= fload_gas_z(1)
	Ngas= n_elements(gas_sfr)

	newstar_mass= 1.0e+10*fload_allstars_mass(1)/0.7
	newstar_time= float(Ttime-fload_newstars_age(1))
	newstar_II= fload_newstars_z(1)
	Nstars= n_elements(newstar_mass)


	Ia_rate= load_SN_rate(age)
	newstar_metalstodistribute= Ia_rate * newstar_mass


	; call a C-program to distribute the time-delayed metals
	; -------------------------------------------------------

        ; fields to pass
        Coord = fltarr(6,Ngas+Nstars)

        ; gas
        Coord(0,0:Ngas-1) = fload_gas_xyz('x', center=[0,0,0])
        Coord(1,0:Ngas-1) = fload_gas_xyz('y', center=[0,0,0])
        Coord(2,0:Ngas-1) = fload_gas_xyz('z', center=[0,0,0])

        Coord(3,Ngas:Ngas+Nstars-1) = fload_newstars_xyz('x', center=[0,0,0])
        Coord(4,Ngas:Ngas+Nstars-1) = fload_newstars_xyz('y', center=[0,0,0])
        Coord(5,Ngas:Ngas+Nstars-1) = fload_newstars_xyz('z', center=[0,0,0])

        Coord(6,Ngas:Ngas+Nstars-1) = newstar_metalstodistribute


        gas_new_Ia= fltarr(Ngas)

        print, " **** PASSING **** "
        print, " Ngas= ", Ngas
        print, " Nstars= ", Nstars
        help, Coord

        S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/DistributeMetals/distz.so', $
                'distz', $
                Ngas, $
                Nstars, $
                Coord, $
                gas_new_Ia)


	; we're currently screwing this
	; up, we need to keep track of
	; the correct id numbers
	;gas_Ia= gas_Ia + gas_new_Ia


endfor



;-------------------------------------------
;
; write new star metal file
;
;-------------------------------------------
openw, 1, frun+metaldir+'/nsmetals_'+snaplbl+'.txt', ERROR=err

printf, 1, "#   run: "+frun
printf, 1, "#   snap: "+snaplbl
printf, 1, "#         "
printf, 1, "#   time     New Star    Formation      Ia          II    "
printf, 1, "#   (Gyr)       ID      Time (Gyr)   Remnants    Remnants "
for i=0,nsnap-1 do begin
        printf, 1, FORMAT= '("  ",F6.3,"    ",4(F10.5,"  "))', $
                ztime[i], newstar_id[i], newstar_time[i], newstar_Ia[i], newstar_II[i]
endfor
close, 1

;-------------------------------------------

;spawn, 'gzip '+frun+metaldir+'/nsmetals_'+snaplbl+'.txt'

;-------------------------------------------


print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"




end



