
;==========================================
;
;    Determine Time Evolution of Energy
;
;
;
;
;==========================================


pro time_energy, frun, fadd=fadd, $
			WindEfficiency=WindEfficiency, $
			WindSpeed=WindSpeed

if not keyword_set(frun) then begin
   print, "  "
   print, "time_energy, frun"
   print, "  "
   print, "  "
   return
endif


if not keyword_set(fadd) then fadd= '' else fadd='_'+fadd



; ---------------------------------
;
;  for the moment, we'll compute
; energies at the same frequency
; as the snapshots

;
; this assumes they are ordered
;  0 through [something]


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])


; set manual time interval
;dt= 0.1
;nsnaps= 31

time= fltarr(nsnaps)

; three types of fb
;
; std - multiphase model component based upon sfr
; wind - extra energy that goes into wind model
; bh - black hole feedback based upon bh accretion rate
;
fb_std_energy= fltarr(nsnaps)
fb_wind_energy= fltarr(nsnaps)
fb_bh_energy= fltarr(nsnaps)


; ----------------------------------------
; This part loops through the snapshots
; and compiles feedback energy information
; ----------------------------------------
snapi= 0

for si=0,nsnaps-1 do begin

        print, "--------------------------------------"

        ; open snapshot
	; ---------------
	thistrial= 0
        repeat begin
            ok=fload_snapshot_bh(frun,snapi)
            ;ok=fload_snapshot(frun,snapi)
	    ;ok= 0
            snapi= snapi+1
	    thistrial= thistrial+1
        endrep until ((ok eq 0) or (thistrial gt 20))

	if thistrial gt 20 then break

        ; what time is it?
        time[si]= fload_time(1)
	;time[si]= si*dt


	;-------------------------------------------
        ;   Get SFR rate from txt - for each file
	;
	;   (this needs sfr_multi.pro)
        ;-------------------------------------------
        open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass


	; get index of current time
        idx= where(sfrtime ge time[si])
	if idx(0) eq -1 then begin
		print, "PROBLEM: can not match times"
		print, "   N= ", n_elements(sfrtime)
		print, "   sfrtime(N-1)= ", sfrtime(n_elements(sfrtime)-1)
		print, "   si= ", si
		print, "   time[si]= ", time[si]
		print, "   "
		stop
	endif
	currentsfr= sfrsfr(idx(0))
	print, "current SFR= ", currentsfr
	if currentsfr le 0.0 then currentsfr= 1.0d-10


	; std fb energy - needed to maintain the multiphase model
	;
	;  for now we calculate this from std GADGET value
	;  as provided in the runtime output0.txt
	;  (if you redirect it to this, as I do)
	;
	;  Feedback energy (in ergs) per formed solar mass in stars
	FBEnergy_per_SN= 1.3955d+49     
	sec_per_year = 3.08568d+7
	fb_std_energy[si]= alog10(FBEnergy_per_SN * currentsfr / sec_per_year)



	; wind model
	;
	;  fb energy is calculated as 
	;
	convert_kms_cms= 1.0d+5
	solarmass_in_g = 1.989d+33
	;sec_per_year = 3.08568d+7     - set above
	if not keyword_set(WindSpeed) then WindSpeed= 837.0
	if not keyword_set(WindEfficiency) then WindEfficiency= 1.0
	WindVel= WindSpeed * convert_kms_cms
	fb_wind_energy[si]= alog10(WindEfficiency * WindVel * WindVel * currentsfr * solarmass_in_g / sec_per_year)



	; ---------------------------------
	;  Add BH Luminosity
	; ---------------------------------
	include_BH_lum= 1
	;include_BH_lum= 0
	if include_BH_lum eq 1 and fload_npart(5) gt 0 then begin

		; get BH lum
		; -------------
        	Nbh= fload_npart(5)
		bhids= fload_blackhole_id(1)
		bhaccrate= fload_blackhole_accrate(frun,time[si])

		if total(bhaccrate) gt 0.0 then begin

			; convert to energy
			cc= 2.9979d+10               ; speed of light in cm/sec
                	;convert_sunyr_gsec= 6.446d+27                  ; convert msun/yr -> g/sec  (wrong value - i think)
                	;convert_sunyr_gsec= 6.30428d+25                  ; convert msun/yr -> g/sec (used 365, or no)
                	;convert_sunyr_gsec= 6.446d+25                  ; convert msun/yr -> g/sec
                	convert_sunyr_gsec= solarmass_in_g / sec_per_year
			BHAccEfficiency= 0.1
			FBCouplingEfficiency= 0.05            ; std value
			bhfbenergy= FBCouplingEfficiency * BHAccEfficiency * cc * cc * bhaccrate * convert_sunyr_gsec

			fb_bh_energy[si]= alog10(total(bhfbenergy))
		endif else begin
			fb_bh_energy[si]= 0.0
		endelse

	endif



endfor


; ------------------------------
;  Done - now write text file
; ------------------------------

openw, 1, frun+'/fb_energy'+fadd+'.txt', ERROR=err

printf, 1, "#   fb_energy.txt (+fadd, maybe)"
printf, 1, "#   Compilation of various forms of FB energy "
printf, 1, "#        "
printf, 1, "# time        Std       Wind         BH "
printf, 1, "#(Gyr)    (erg/s)    (erg/s)    (erg/s)" 
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"  ",3(F9.4,"  "))', $
                time[i], fb_std_energy[i], fb_wind_energy[i], fb_bh_energy[i]
endfor
close, 1




end








; ================================================================================










; --------------------------------
;  Read FB Energy File
; ----------------------------------
pro read_fbenergy_file, frun, time, fb_std, fb_wind, fb_bh, $
				fbfile=fbfile

if keyword_set(fbfile) then fbfile= frun+'/'+fbfile else fbfile= frun+'/fb_energy.txt'

spawn, "wc "+fbfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then re_data= fltarr(4,lines)

openr, 1, fbfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, re_data
close, 1


time= re_data[0,*]
fb_std= re_data[1,*]
fb_wind= re_data[2,*]
fb_bh= re_data[3,*]


end






; ================================================================================







;-------------------------------------------------
;-------------------------------------------------
pro plot_fbenergy, junk


if not keyword_set(junk) then begin
        print, " "
        print, " plot_fbenergy, junk"
        print, " "
        print, " "
        return
endif

filename='fbenergy.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

;---------------------------

xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xmax = 3.0
xmin = 0.0

yaxistitle = "!6Log Energy (erg s!E-1!N)"
ymax = 46.0
ymin = 39.0


;---------------------------

!p.position= [0.16, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------


	;frun= "/raid4/tcox/vc3vc3e_2"
	;frun= "/raid4/tcox/vc3vc3e_sb8"
	;frun= "/raid4/tcox/vc3vc3e_sb1"
	;frun= "/raid4/tcox/vc3vc3e_sb10"
	frun= "/raid4/tcox/sbw/sb10BH"
	;frun= "/raid4/tcox/vc3vc3e_sb13"
	read_fbenergy_file, frun, time, fb_std, fb_wind, fb_bh

        oplot, time, fb_wind, thick=4.0, color= 50, psym=-2, linestyle=0
        oplot, time, fb_std, thick=4.0, color= 10, psym=-3
        oplot, time, fb_bh, thick=4.0, color= 150, psym=-5, linestyle=0

	xyouts, 0.68, 0.85, "Std", /normal, color= 10, charthick=3.0, charsize=1.3
	xyouts, 0.68, 0.80, "Wind", /normal, color= 50, charthick=3.0, charsize=1.3
	xyouts, 0.68, 0.75, "BH", /normal, color= 150, charthick=3.0, charsize=1.0
	

;---------------------------




; done
; ------
device, /close




end




; -----------------------------------------------------
; -----------------------------------------------------





pro plot_int_fbenergy, junk


if not keyword_set(junk) then begin
        print, " "
        print, " plot_int_fbenergy, junk"
        print, " "
        print, " "
        return
endif

filename='fbenergy_int.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

;---------------------------

xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xmax = 3.0
xmin = 0.0

yaxistitle = "!6Log Integrated FB Energy (erg)"
ymax = 60.0
ymin = 56.0


;---------------------------

!p.position= [0.16, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------


	;frun= "/raid4/tcox/vc3vc3e_2"
	;frun= "/raid4/tcox/vc3vc3e_sb8"
	;frun= "/raid4/tcox/vc3vc3e_sb1"
	;frun= "/raid4/tcox/vc3vc3e_sb10"
	frun= "/raid4/tcox/sbw/sb10BH"
	;frun= "/raid4/tcox/vc3vc3e_sb13"
	read_fbenergy_file, frun, time, fb_std, fb_wind, fb_bh

	gadunits_in_sec= 3.08568d+16

	ntime= n_elements(time)
	int_fb_wind= fltarr(ntime)
	int_fb_std= fltarr(ntime)
	int_fb_bh= fltarr(ntime)

	; put in more managable units
	unitsof40= 40.0
	fb_wind= fb_wind - unitsof40
	fb_std= fb_std - unitsof40
	fb_bh= fb_bh - unitsof40

	tot_wind_e= 0.0
	tot_std_e= 0.0
	tot_bh_e= 0.0

	for i=1,ntime-1 do begin

		;dt= (time[i]-time[i-1]) * gadunits_in_sec
		dt= (time[i]-time[i-1])

		avg_wind_e= dt * 0.5 * ((10^fb_wind[i-1]) + (10^fb_wind[i]))
		avg_std_e= dt * 0.5 * ((10^fb_std[i-1]) + (10^fb_std[i]))
		avg_bh_e= dt * 0.5 * ((10^fb_bh[i-1]) + (10^fb_bh[i]))

		if avg_wind_e gt 0.0 then tot_wind_e= tot_wind_e + avg_wind_e
		if avg_std_e gt 0.0 then tot_std_e= tot_std_e + avg_std_e
		if avg_bh_e gt 0.0 then tot_bh_e= tot_bh_e + avg_bh_e

		int_fb_wind[i]= unitsof40 + alog10(gadunits_in_sec) + alog10(tot_wind_e)
		int_fb_std[i]= unitsof40 + alog10(gadunits_in_sec) + alog10(tot_std_e)
		int_fb_bh[i]= unitsof40 + alog10(gadunits_in_sec) + alog10(tot_bh_e)

	endfor

        oplot, time, int_fb_wind, thick=4.0, color= 50, psym=-2, linestyle=0
        oplot, time, int_fb_std, thick=4.0, color= 10, psym=-3
        oplot, time, int_fb_bh, thick=4.0, color= 150, psym=-5, linestyle=0

	xyouts, 0.68, 0.85, "Std", /normal, color= 10, charthick=3.0, charsize=1.3
	xyouts, 0.68, 0.80, "Wind", /normal, color= 50, charthick=3.0, charsize=1.3
	xyouts, 0.68, 0.75, "BH", /normal, color= 150, charthick=3.0, charsize=1.0
	

;---------------------------




; done
; ------
device, /close




end







;-------------------------------------------------
;-------------------------------------------------



pro plot_fbepergas, junk


if not keyword_set(junk) then begin
        print, " "
        print, " plot_fbepergas, junk"
        print, " "
        print, " "
        return
endif

filename='fbepergas.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

;---------------------------

xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xmax = 3.0
xmin = 0.0

yaxistitle = "!6Log FB Energy per M!Dgas!N (erg s!E-1!N M!D!9n!6!N!E-1!N)"
ymax = 36.0
ymin = 29.0


;---------------------------

!p.position= [0.16, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------


	frun= "/raid4/tcox/vc3vc3e_2"
	;frun= "/raid4/tcox/vc3vc3e_sb8"
	;frun= "/raid4/tcox/vc3vc3e_sb1"
	;frun= "/raid4/tcox/vc3vc3e_sb10"
	;frun= "/raid4/tcox/vc3vc3e_sb13"
	read_fbenergy_file, frun, time, fb_std, fb_wind, fb_bh


	; need time_hotgas compiles for this
	read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf

	gas_tot= gas_tot * 1.0d+10

	fb_wind_gas= fb_wind - alog10(gas_tot)
	fb_std_gas= fb_std - alog10(gas_tot)
	fb_bh_gas= fb_bh - alog10(gas_tot)

        oplot, time, fb_wind_gas, thick=4.0, color= 50, psym=-2, linestyle=0
        oplot, time, fb_std_gas, thick=4.0, color= 10, psym=-3
        oplot, time, fb_bh_gas, thick=4.0, color= 150, psym=-5, linestyle=0

	xyouts, 0.68, 0.85, "Std", /normal, color= 10, charthick=3.0, charsize=1.3
	xyouts, 0.68, 0.80, "Wind", /normal, color= 50, charthick=3.0, charsize=1.3
	xyouts, 0.68, 0.75, "BH", /normal, color= 150, charthick=3.0, charsize=1.0
	

;---------------------------




; done
; ------
device, /close




end




; -----------------------------------------------------
; -----------------------------------------------------









;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------


