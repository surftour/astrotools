pro methods_phase, junk


; warning, might need to check over the script here as things might
;  have changed a little, i.e. might need to turn off printing
;  of the fid.

   frun= "/data/tcox/Sbc201a-u4"
   sendto= 'ps'
   snapnum= 35

  ;contour_phasediagram, frun, snapnum, sendto, filename='phase2l.eps', msg='n2low', /drawindexn, indexn=2
  contour_phasediagram, frun, snapnum, sendto, filename='fig9c.eps', msg='n2low', /drawindexn, indexn=2

  ;frun= "/data/tcox/Sbc201a-u40"
  frun="execute/Sbc201a-u53"
  ;contour_phasediagram, frun, snapnum, sendto, filename='phase1l.eps', msg='n1low', /drawindexn, indexn=1
  contour_phasediagram, frun, snapnum, sendto, filename='fig9b.eps', msg='n1low', /drawindexn, indexn=1

  ;frun= "/data/tcox/Sbc201a-u8"
  frun="execute/Sbc201a-u50"
  ;frun= "/data/tcox/Sbc201a-u50"
  ;contour_phasediagram, frun, snapnum, sendto, filename='phase0l.eps', msg='n0low', /drawindexn, indexn=0, /gaslabels
  contour_phasediagram, frun, snapnum, sendto, filename='fig9a.eps', msg='n0low', /drawindexn, indexn=0, /gaslabels

  frun= "/data/tcox/Sbc201a-u43"
  ;contour_phasediagram, frun, snapnum, sendto, filename='phase2h.eps', msg='n2high', /drawindexn, indexn=2
  ;contour_phasediagram, frun, snapnum, sendto, filename='phase2h.eps', msg='n2high'
  ; later changed to drawing the line
  ; but with an increased amplitude
  ;contour_phasediagram, frun, snapnum, sendto, filename='phase2h.eps', msg='n2high', /drawindexn, indexn=2, /increaseamp
  contour_phasediagram, frun, snapnum, sendto, filename='fig9d.eps', msg='n2high', /drawindexn, indexn=2, /increaseamp

end





pro methods_iso_phase, junk


; warning, might need to check over the script here as things might
;  have changed a little, i.e. might need to turn off printing
;  of the fid.

   sendto= 'ps'
   snapnum= 10

  contour_phasediagram, "execute/Sbc11i4-u57", snapnum, sendto, filename='isop2l.eps', msg='n2low', /drawindexn, indexn=2
  contour_phasediagram, "/data/tcox/Sbc11i4-u4", snapnum, sendto, filename='isop2m.eps', msg='n2med', /drawindexn, indexn=2
  contour_phasediagram, "execute/Sbc11i4-u43", snapnum, sendto, filename='isop2h.eps', msg='n2high', /drawindexn, indexn=2

  contour_phasediagram, "execute/Sbc11i4-u56", snapnum, sendto, filename='isop1l.eps', msg='n1low', /drawindexn, indexn=1
  contour_phasediagram, "execute/Sbc11i4-u53", snapnum, sendto, filename='isop1m.eps', msg='n1med', /drawindexn, indexn=1
  contour_phasediagram, "execute/Sbc11i4-u52", snapnum, sendto, filename='isop1h.eps', msg='n1high', /drawindexn, indexn=1

  contour_phasediagram, "execute/Sbc11i4-u55", snapnum, sendto, filename='isop0l.eps', msg='n0low', /drawindexn, indexn=0
  contour_phasediagram, "execute/Sbc11i4-u50", snapnum, sendto, filename='isop0m.eps', msg='n0med', /drawindexn, indexn=0
  contour_phasediagram, "execute/Sbc11i4-u51", snapnum, sendto, filename='isop0h.eps', msg='n0high', /drawindexn, indexn=0

  ;contour_phasediagram, frun, snapnum, sendto, filename='phase0l.eps', msg='n0low', /drawindexn, indexn=0, /gaslabels
  ; later changed to drawing the line
  ; but with an increased amplitude
  ;contour_phasediagram, frun, snapnum, sendto, filename='phase2h.eps', msg='n2high', /drawindexn, indexn=2, /increaseamp

end





pro contour_phase_multi_9, junk



if not keyword_set(junk) then begin
   print, "  "
   print, "contour_phase_multi, junk"
   print, "  "
   return
endif

filename='phase_multi.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3, newxsize=16, newysize=14


;--------------------------------------
;--------------------------------------

;
;  WARNING: xmax, xmin, ymax, ymin, bins,
;   xaxistitle, yaxistitle, are all set
;   in process_one_phase_file
;

;  -------------------
;  |     |     |     |
;  |     |     |     |
;  |     |     |     |
;  -------------------
;  |     |     |     |
;  |     |     |     |
;  |     |     |     |
;  -------------------
;  |     |     |     |
;  |     |     |     |
;  |     |     |     |
;  -------------------
;


x0= 0.09
xsize= 0.30
x1= x0+xsize
x2= x0+xsize+xsize
x3= x0+xsize+xsize+xsize

y0= 0.11
ysize= 0.29
y1= y0+ysize
y2= y0+ysize+ysize
y3= y0+ysize+ysize+ysize

; ------------------------------------

;process_one_phase_file, 'execute/I1i-u8', x0, y2, x1, y3, /showyaxis, msg='!8n0low!6', teff= 2.5e4, indexn=0
;process_one_phase_file, 'execute/I1i-u2', x1, y2, x2, y3, msg='!8n0med!6', teff= 2.5e5, indexn= 0
;process_one_phase_file, 'execute/I1i-u9', x2, y2, x3, y3, msg='!8n0high!6', teff= 2.5e6, indexn= 0

;process_one_phase_file, 'execute/I1i-u16', x0, y1, x1, y2, /showyaxis, msg='!8n1low!6', teff= 2.5e4, indexn= 1
;process_one_phase_file, 'execute/I1i-u3',  x1, y1, x2, y2, msg='!8n1med!6', teff= 2.5e5, indexn= 1
;process_one_phase_file, 'execute/I1i-u17', x2, y1, x3, y2, msg='!8n1high!6', teff= 2.5e6, indexn= 1

;process_one_phase_file, 'execute/I1i-u5', x0, y0, x1, y1, /showyaxis, /showxaxis, msg='!8n2low!6', teff= 2.5e4, indexn= 2
;process_one_phase_file, 'execute/I1i-u1', x1, y0, x2, y1, /showxaxis, msg='!8n2med!6', teff= 2.5e5, indexn= 2
;process_one_phase_file, 'execute/I1i-u7', x2, y0, x3, y1, /showxaxis, msg='!8n2high!6', teff= 2.5e6, indexn= 2

; ------------------------------------

;process_one_phase_file, 'execute/Sbc11i4-u55', x0, y2, x1, y3, /showyaxis, msg='!8n0low!6', teff= 2.5e4, indexn=0
;process_one_phase_file, 'execute/Sbc11i4-u50', x1, y2, x2, y3, msg='!8n0med!6', teff= 1.6e5, indexn=0
;process_one_phase_file, 'execute/Sbc11i4-u51', x2, y2, x3, y3, msg='!8n0high!6', teff= 1.6e6, indexn=0, /xaxissolar

;process_one_phase_file, 'execute/Sbc11i4-u56', x0, y1, x1, y2, /showyaxis, msg='!8n1low!6', teff= 2.5e4, indexn=1
;process_one_phase_file, 'execute/Sbc11i4-u53', x1, y1, x2, y2, msg='!8n1med!6', teff= 1.6e5, indexn=1
;process_one_phase_file, 'execute/Sbc11i4-u52', x2, y1, x3, y2, msg='!8n1high!6', teff= 1.6e6, indexn=1, /xaxissolar

;process_one_phase_file, 'execute/Sbc11i4-u57',   x0, y0, x1, y1, /showyaxis, /showxaxis, msg='!8n2low!6', teff= 2.5e4, indexn=2
;;process_one_phase_file, '/data/tcox/Sbc11i4-u4', x1, y0, x2, y1, /showxaxis, msg='!8n2med!6', teff= 2.5e5, indexn=2
;process_one_phase_file, '/data/tcox/Sbc11i4-u4', x1, y0, x2, y1, msg='!8n2med!6', teff= 1.6e5, indexn=2
;process_one_phase_file, 'execute/Sbc11i4-u43',   x2, y0, x3, y1, /showxaxis, msg='!8n2high!6', teff= 1.6e6, indexn=2, /xaxissolar

; ------------------------------------

process_one_phase_file, 'execute/Sbc201a-u55', x0, y2, x1, y3, /showyaxis, msg='!8n0low!6', teff= 2.5e4, indexn=0
process_one_phase_file, 'execute/Sbc201a-u50', x1, y2, x2, y3, msg='!8n0med!6', teff= 1.6e5, indexn=0
process_one_phase_file, 'execute/Sbc201a-u51', x2, y2, x3, y3, msg='!8n0high!6', teff= 1.6e6, indexn=0, /xaxissolar

process_one_phase_file, 'execute/Sbc201a-u56', x0, y1, x1, y2, /showyaxis, msg='!8n1low!6', teff= 2.5e4, indexn=1
process_one_phase_file, 'execute/Sbc201a-u53', x1, y1, x2, y2, msg='!8n1med!6', teff= 1.6e5, indexn=1
process_one_phase_file, 'execute/Sbc201a-u52', x2, y1, x3, y2, msg='!8n1high!6', teff= 1.6e6, indexn=1, /xaxissolar

process_one_phase_file, 'execute/Sbc201a-u57',   x0, y0, x1, y1, /showyaxis, /showxaxis, msg='!8n2low!6', teff= 2.5e4, indexn=2
process_one_phase_file, '/data/tcox/Sbc201a-u4', x1, y0, x2, y1, msg='!8n2med!6', teff= 1.6e5, indexn=2
process_one_phase_file, '/data/tcox/Sbc201a-u43',   x2, y0, x3, y1, /showxaxis, msg='!8n2high!6', teff= 1.6e6, indexn=2, /xaxissolar


xyouts, 0.10, 0.892, 'HOT', size=1.0, color=0, /normal, charthick=3.0
xyouts, 0.10, 0.72, 'COLD', size=1.0, color=0, /normal, charthick=3.0
xyouts, 0.28, 0.75, 'STAR-', size=1.0, color=0, /normal, charthick=3.0
xyouts, 0.285, 0.72, 'FORMING', size=1.0, color=0, /normal, charthick=3.0


; ------------------------------------

device, /close


end




;
;
;===================================================================



;
;
;
; ---------------------------------------------------
pro process_one_phase_file, frun, x0, y0, x1, y1, $
                                showyaxis=showyaxis, showxaxis=showxaxis, $
                                msg=msg, $
                                secmsg=secmsg, $
                                showdelta=showdelta, $
				teff=teff, $
				indexn=indexn, xaxissolar=xaxissolar



xaxistitle='!6Log Density (cm!E-3!N)'
;xmax = 3.5    ; isolated disks
;xmin = -5.0
xmax = 4.5    ; Sbc major mergers
xmin = -7.0


yaxistitle='!6Log T!Deff!N (K)'
ymax = 7.5
ymin = 2.5



bins= 60



;snapnum= 5   ; isolated disks
snapnum= 35  ; Sbc majorm mergers

;-------------------------------------
;  Load Snapshot
;-------------------------------------

if not keyword_set(loadedsnap) then begin
   ;if (fload_snapshot(frun, snapnum)) then begin
   ;if (fload_snapshot_bh(frun, snapnum,nopot_in_snap=nopot_in_snap)) then begin
   if (fload_snapshot(frun, snapnum)) then begin
        print, "PROBLEM: opening file"
        return
   endif
endif


temp= (fload_gas_u(1)+fload_gas_tpu(1))*fload_gas_mu(1)/0.012381322
;temp= (fload_gas_u(1)+fload_gas_tpu(1))/0.012381322
;temp= (fload_gas_u(1)*fload_gas_mu(1))/0.012381322
temp= alog10(temp)

; convert this to Msolar pc-3
;rho= 10*fload_gas_rho(1)

; convert to cm-3
UnitDensity_in_cgs = 6.76991d-22
ProtonMass = 1.6726d-24
density_factor= UnitDensity_in_cgs / ProtonMass
print, "density factor= ", density_factor
rho = fload_gas_rho(1) * density_factor

if keyword_set(xaxissolar) then begin
	rho = 10.0 * rho / density_factor
	xmax= alog10(10^(xmax) / density_factor)     ; actually missing a factor of 10, but the graph looks good
	xmin= alog10(10^(xmin) / density_factor)
endif

rho=alog10(rho)



; ------------------
; compute histogram
; ------------------

contour_makegeneralpic, rho, temp, xmax, xmin, ymax, ymin, $
                                pixels= bins, $
                                NxNImage=NxNImage


;print, min(NxNImage), max(NxNImage)


x_size= (x1-x0)
y_size= (y1-y0)


tv, NxNImage, x0, y0, xsize=x_size, ysize=y_size, /normal




;----------------------------------------
; Generate plot
;----------------------------------------

!p.position=[x0, y0, x0+x_size,y0+y_size]

;tnames=['0',' ','0.2',' ','0.4',' ','0.6',' ']
;tnames=[' ','-2',' ','0',' ','2',' ','4',' ','6',' ','8',' ','10']

if keyword_set(showyaxis) and keyword_set(showxaxis) then begin
   plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $ 
        xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
        ;xticks= 7, xtickname=tnames, $
        xtitle=xaxistitle, ytitle=yaxistitle
endif else begin
   if keyword_set(showyaxis) then begin
        plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
        ;xticks= 7, $
        ytitle=yaxistitle, xtickformat='(a1)'
   endif

   if keyword_set(showxaxis) and not keyword_set(xaxissolar) then begin
        plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
        ;xticks= 7, xtickname=tnames, $
        xtitle=xaxistitle, ytickformat='(a1)'
   endif

   if keyword_set(showxaxis) and keyword_set(xaxissolar) then begin
	xaxistitle='!6Log Density (M!D!9n!6!N pc!E-3!N)'
	;xmax= alog10(10^(xmax) / density_factor)
	;xmin= alog10(10^(xmin) / density_factor)
        plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
        ;xticks= 7, xtickname=tnames, $
        xtitle=xaxistitle, ytickformat='(a1)'
   endif

   if ((not keyword_set(showyaxis)) and (not keyword_set(showxaxis))) then begin

	if keyword_set(xaxissolar) then begin
		;xmax= alog10(10^(xmax) / density_factor)
		;xmin= alog10(10^(xmin) / density_factor)
        	plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        	xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
        	;xticks= 7, $
        	xtickformat='(a1)', ytickformat='(a1)'
	endif else begin
                plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
                xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
                ;xticks= 7, $
                xtickformat='(a1)', ytickformat='(a1)'
	endelse
   endif

endelse



phasediagram= 1

        rhoc= 0.00171 * 404.75368   ; in cm-3

        if keyword_set(phasediagram) then begin
           if keyword_set(sfrhothresh) then begin
                y=[ymin,ymax]
                x= sfrhothresh + 0.0*y
                x= alog10(x)
                oplot, x, y, linestyle=0, color=0
                ok=oplot_coldtemp_cutoff_log(1,rhocut=sfrhothresh)
           endif else begin
                ;ok=oplot_rho_cutoff_y_log(1)
                ;ok=oplot_coldtemp_cutoff_log(1)

                ;ok=oplot_rho_cutoff_y_log(1,rhocut=0.00854924)
                ;ok=oplot_coldtemp_cutoff_log(1,rhocut=0.00854924)

                ;rhoc= 0.000854924 * 404.75368   ; in cm-3
                ;rhoc= 0.00171 * 404.75368   ; in cm-3

		if keyword_set(xaxissolar) then begin
			rhoc= 10.0 * rhoc / 404.75368   ; in solar mass / pc^3
		endif

                ok=oplot_rho_cutoff_y_log(1,rhocut=rhoc)
                ok=oplot_coldtemp_cutoff_log(1,rhocut=rhoc)
           endelse
        endif

        if keyword_set(plotcoolingtime) then begin
                ;ok=oplot_coolingtime_log(1.0,sendto=sendto)
                ;ok=oplot_coolingtime_log(0.001,sendto=sendto)
                ;ok=oplot_coolingtime_log(10.0,sendto=sendto)
                ;ok=oplot_coolingtime_log(100.0,sendto=sendto)
                ;xyouts, 0.37, 0.87, 't!Dcool!N=1 Gyr', size= 1.2, color= 0, /normal, charthick= 3.0
                ;xyouts, 0.47, 0.87, '1Gyr', size= 1.2, color= 0, /normal, charthick= 3.0
                ;xyouts, 0.62, 0.87, '1Myr', size= 1.2, color= 0, /normal, charthick= 3.0
        endif


        drawindexn= 1
        if keyword_set(drawindexn) then begin
                len= 1.0
                slope= indexn/2.0
                ;slope= indexn
                ;theta= atan(slope)
                ;x0n=0.0
		x0n= alog10(rhoc)
                ;y0n=4.8
		y0n= alog10(teff)
                ;x1= x0+(len*cos(theta))*(xmax-xmin)/(ymax-ymin)
                ;y1= y0+(len*sin(theta))
                x1n= xmax
                y1n= y0n+slope * (x1n-x0n)
                xn=[x0n,x1n]
                yn=[y0n,y1n]
                ;oplot, x, y, linestyle=0, color=0, thick=3.0
                oplot, xn, yn, linestyle=1, color=0, thick=3.0
        endif


        if keyword_set(msg) then begin
                xyouts, x0+0.02, y1-0.04, msg, size=1.5, color=0, /normal, charthick=3.0
                ;xyouts, x0+0.02, y1-0.08, fload_timelbl(1,2), size=1.33, color=0, /normal, charthick= 3.0
        endif else begin
                xyouts, 0.20, 0.90, fload_fid(1), size=1.33, color=0, /normal
                xyouts, 0.20, 0.85, fload_timelbl(1,2), size=1.33, color=0, /normal, charthick= 3.0
        endelse

end




