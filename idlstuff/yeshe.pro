; ------------------------------------------------------------------------------------------
;
;
;  need to .run
;    1. time_hotgas
;    2. sfr_multi
;
;
pro prep_basicinfo, junk


if not keyword_set(junk) then begin
	print, " "
	print, " prep_basicinfo, junk"
	print, " "
	print, " "
	return
endif


;------------------------------
; first series (f=0.4, q=0.25)
;
;do_one_galaxy, "ds/d0e2_q"
;do_one_galaxy, "ds/d1e2_q"
;do_one_galaxy, "ds/d2e2_q"
;do_one_galaxy, "ds/d3e7"
;do_one_galaxy, "ds/d4e2_q"
;do_one_galaxy, "ds/d5e2_q"
;do_one_galaxy, "ds/d6e2_q"
;
;------------------------------
; second series (f=0.4, q=1.0)
;
;do_one_galaxy, 'ds/d0e', extra='50'
;do_one_galaxy, 'ds/d1e', extra='80'
;do_one_galaxy, 'ds/d2e', extra='115'
;do_one_galaxy, 'ds/d3e7', extra='160'
;do_one_galaxy, 'ds/d4e', extra='225'
;do_one_galaxy, 'ds/d5e', extra='320'
;do_one_galaxy, 'ds/d6e', extra='500'
;
;-------------------------
; third series (f=0.8, q=1.0)
;
;do_one_galaxy, 'bs/b0e', extra='50'
;do_one_galaxy, 'bs/b1e', extra='80'
;do_one_galaxy, 'bs/b2e', extra='115'
;do_one_galaxy, 'bs/b3e', extra='160'
;do_one_galaxy, 'bs/b4e', extra='225'
;do_one_galaxy, 'bs/b5e', extra='320'
;do_one_galaxy, 'bs/b6e', extra='500'
;
;-------------------------
; series #4 (f=0.2, q=0.25)
;
;do_one_galaxy, 'es/e0e', extra='50'
;do_one_galaxy, 'es/e1e', extra='80'
;do_one_galaxy, 'es/e2e', extra='115'
;do_one_galaxy, 'es/e3e', extra='160'
;do_one_galaxy, 'es/e4e', extra='225'
;do_one_galaxy, 'es/e5e', extra='320'
;do_one_galaxy, 'es/e6e', extra='500'
;
;-------------------------
; series #5 (f=0.05, q=0.25)
;
do_one_galaxy, 'zs/z0e', extra='50'
do_one_galaxy, 'zs/z1e', extra='80'
do_one_galaxy, 'zs/z2e', extra='115'
do_one_galaxy, 'zs/z3e', extra='160'
do_one_galaxy, 'zs/z4e', extra='225'
do_one_galaxy, 'zs/z5e', extra='320'
do_one_galaxy, 'zs/z6e', extra='500'
;
;-------------------------





end




;=========================================================================





pro do_one_galaxy, frun, extra=extra

;-------------------------------------------

fdata="/raid4/tcox/"

;-------------------------------------------


read_hotgas_file, fdata+frun, time, temp_keV_X, temp_keV, temp_K, entropy, $
				gas_tot, gas_hot, gas_cold, gas_sf


open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass


;-------------------------------------------

n_time= n_elements(time)
sfr_inst= fltarr(n_time)
sfr_avg= fltarr(n_time)

for i=0, n_time-1 do begin
	; get the *.sfr index
	idx_gtr_snaptime= where(sfrtime ge time[i])
	idx_curr_snaptime= idx_gtr_snaptime[0]


	; instantaneous sfr
	sfr_inst[i]= sfrsfr(idx_curr_snaptime)


	; set time window
	dt= 0.01
	idx_window= where((sfrtime ge (time[i]-dt)) and (sfrtime le (time[i]+dt)))
	sfr_dt= sfrsfr(idx_window)
	sfr_avg[i]= mean(sfr_dt)

	print, "T= ", time[i], sfr_inst[i], sfr_avg[i]
endfor


time= time/0.7
gas_tot= gas_tot/0.7
gas_sf= gas_sf/0.7
gas_cold= gas_cold/0.7
gas_hot= gas_hot/0.7


;-------------------------------------------
;
; write 
; Yeshe Gas Info File
;
;--------------------
openw, 1, fdata+frun+'/yinfo_'+extra+'.txt', ERROR=err

printf, 1, "#   yinfo for run:"+fdata+frun
printf, 1, "#   "
printf, 1, "#                       gas mass"
printf, 1, "# time          SFR        total         sf        cold         hot  "
printf, 1, "# (Gyr)  (Msolar/Yr)    (10^10 msolar ---->)   "
for i=0,n_time-1 do begin
        printf, 1, FORMAT= '(F6.3,"    ",F10.5,"  ",4(F10.5,"  "))', $
                time[i], sfr_avg[i], gas_tot[i], gas_sf[i], gas_cold[i], gas_hot[i]
endfor
close, 1

;-------------------------------------------

spawn, 'cp '+fdata+frun+'/yinfo_'+extra+'.txt  /home/tcox/'

;-------------------------------------------


end



