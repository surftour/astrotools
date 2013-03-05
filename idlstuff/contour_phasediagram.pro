pro contour_phasediagram, frun, snapnum, sendto, filename=filename, $
			arepo=arepo, $
			msg=msg, $
			q_eos=q_eos, $
			sfrhothresh=sfrhothresh, $
			drawindexn=drawindexn, indexn=indexn, $
			colortable=colortable, $
			loadedsnap=loadedsnap, nopot_in_snap=nopot_in_snap

if not keyword_set(sendto) then sendto='x'
if not keyword_set(snapnum) then snapnum=0
if not keyword_set(frun) then begin
   print, "  "
   print, "contour_phasediagram, frun, snapnum, sendto, filename=filename,"
   print, "                 msg=msg, sfrhothresh=sfrhothresh, indexn=indexn,"
   print, "                 loadedsnap=loadedsnap"
   print, "  "
   return
endif


        
;-------------------------------------
;  Setup Plot Stuff
;-------------------------------------

;initialize_plotinfo, 1

;setup_plot_stuff, sendto, filename=filename, colortable=colortable



;--------------------------------------
;  Set plotting max/min
;--------------------------------------

; this is for m_solar/pc^3
;xmax = 5
;xmin = -12
; this is for cm-3
xmax = 7
xmin = -6
ymax = 8
ymin = 2.5



;-------------------------------------
;  Load Snapshot
;-------------------------------------

if not keyword_set(loadedsnap) then begin
   ok= -1
   if keyword_set(arepo) then ok= fload_snapshot_bh(frun, snapnum, /nopot_in_snap, /arepo)
   if ok lt 0 then ok= fload_snapshot_bh(frun, snapnum, /nopot_in_snap)
   ;if (fload_snapshot(frun, snapnum)) then begin
   ;if (fload_snapshot_bh(frun, snapnum,nopot_in_snap=nopot_in_snap)) then begin
   ;if (fload_snapshot_bh(frun, snapnum,/nopot_in_snap)) then begin
   if ok eq 1 then begin
	print, "PROBLEM: opening file"
	return
   endif
endif



;temp= (fload_gas_u(1)+fload_gas_tpu(1))*fload_gas_mu(1)/0.012381322
;temp= (fload_gas_u(1)+fload_gas_tpu(1))/0.012381322
;temp= (fload_gas_u(1)*fload_gas_mu(1))/0.012381322
temp= fload_gas_temperature(1)


; convert this to Msolar pc-3
;rho= 10*fload_gas_rho(1)

; convert to cm-3
UnitDensity_in_cgs = 6.76991d-22
ProtonMass = 1.6726d-24
rho = fload_gas_rho(1) * UnitDensity_in_cgs / ProtonMass


if keyword_set(q_eos) then begin
	rhoc= 0.000854924 * 404.75368   ; critical density in cm-3
	idx= where(rho ge rhoc)
	temp(idx)= q_eos*temp(idx) + (1-q_eos)*1.0e+4   ; interpolate between 10^4
endif else begin
	q_eos= 1.0
endelse



temp= alog10(temp)
rho=alog10(rho)






;----------------------------------------
; Send it to plotting routines
;----------------------------------------

sendto= 'ps'


contour_makegeneralplot, rho, temp, xmax, xmin, ymax, ymin, sendto, filename=filename, $
                        msg=msg, $
                        ;xaxistitle='Log Density (M!D!9n!3!N/pc!E3!N)', $
                        xaxistitle='!6Log Density (cm!E-3!N)', $
			yaxistitle='!6Log T!Deff!N (K)', $
                        /phasediagram, $
			sfrhothresh=sfrhothresh, $
			drawindexn=drawindexn, indexn=indexn, $
			;/plotcoolingtime, $
                        colortable=colortable, $
			q_eos=q_eos




end


