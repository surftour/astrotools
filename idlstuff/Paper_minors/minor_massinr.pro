pro minor_massinr, sendto, filename=filename

if not keyword_set(sendto) then begin
   print, "  "
   print, "  "
   print, "minor_massinr, sendto, filename=filename"
   print, "  "
   print, "  "
   return
endif


initialize_plotinfo, 1
;setup_plot_stuff, sendto, filename=filename
setup_plot_stuff, sendto, filename=filename, colortable= 4



normalized= 1
;normalized= 0

rlog= 1
;rlog= 0



;-------------
;  Plot axes 
;-------------

;if rlog eq 1 then xmax = alog10(50.0) else xmax = 25.0
if rlog eq 1 then xmax = 2.0 else xmax = 25.0
if rlog eq 1 then xmin = -1.25 else xmin = 0.0
if normalized eq 1 then ymax = 1.0 else ymax = 2.0
ymin = 0.0

if rlog eq 1 then xaxistitle='Log R (kpc)' else xaxistitle="R (kpc)"
yaxistitle="New Stellar Mass (r<R) (10!E10!N M!D!9n!3!N)"


logplot = 0
if logplot eq 1 then begin
        xmax = alog10(xmax)
        xmin = 0.1
        xaxistitle="Log R (kpc)"
endif


!p.position= [0.20, 0.15, 0.95, 0.95]

plot, [10.0], [10.0], psym=-3, linestyle=0, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata





;-----------------------
;  Put data on graph
;-----------------------


; Isolated
; ----------
snapnum= 500
frun= "execute/G3il-u1a"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-3, linestyle=3, color= 0, thick=1.0



; Mergers
; ---------
snapnum= 500
symsize= 1.5
;usersym,symsize*[-1,-1,1,1],symsize*[-1,1,1,-1],/fill
frun= "execute/G3G3b-u1"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-3, linestyle=0, color= 150, thick=5.0

snapnum= 510
frun= "execute/G3G2-u3"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-3, linestyle=2, color= 100, thick=5.0

snapnum= 500
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0   ;,/fill
frun= "execute/G3G1-u3"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-3, linestyle=1, color= 50, thick=5.0

snapnum= 500
;frun= "execute/G3G0e-u3"
;process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
;if normalized eq 1 then mass_inr= mass_inr/totalmass
;oplot, r_m, mass_inr, psym=-3, linestyle=1, color= 100, thick=1.0




; print arrow for smoothing length
smoothlen= 0.1
if rlog eq 1 then smoothlen= alog10(smoothlen)
toppos= 1.0
botpos= 0.5
if normalized eq 1 then begin
	toppos= 0.45
	botpos= 0.30
endif
xarrow= smoothlen
arrow, xarrow, toppos, xarrow, botpos, /data, COLOR=0, THICK=4.0, hthick=4.0




if normalized eq 1 then begin
	xyouts, 0.80, 0.83, 'G3G3', /normal, color= 150, charthick=2.0, size=1.2
        xyouts, 0.80, 0.54, 'G3G2', /normal, color= 100, charthick=2.0, size=1.2
        xyouts, 0.80, 0.35, 'G3G1', /normal, color= 50, charthick=2.0, size=1.2
        xyouts, 0.80, 0.18, 'G3G0', /normal, color= 0, charthick=2.0, size=1.2
endif else begin
	xyouts, 0.80, 0.83, 'G3G3', /normal, color= 150, charthick=2.0, size=1.2
	xyouts, 0.80, 0.54, 'G3G2', /normal, color= 100, charthick=2.0, size=1.2
	xyouts, 0.80, 0.35, 'G3G1', /normal, color= 50, charthick=2.0, size=1.2
	xyouts, 0.80, 0.18, 'G3G0', /normal, color= 0, charthick=2.0, size=1.2
endelse





;--------------------------------------
if (sendto EQ 'ps') then device, /close


end






;--------------------------------------
; Actually process the massinr 
;--------------------------------------
pro process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass

	bins= 20
	binsize= (xmax-xmin)/bins
	if xmin lt 0 then rlog= 1 else rlog= 0
	r_m= fltarr(bins+1)
	mass_inr= fltarr(bins+1)

	ok=fload_snapshot(frun,snapnum)

	;rxy = fload_luminousbaryons_xyz('rxy')
	;mass = fload_luminousbaryons_mass(1)
        r_s = fload_newstars_xyz('r')
        mass = fload_newstars_mass(1)

	totalmass= total(mass)

	r_m[0]= 0
	mass_inr[0]= 0

	if rlog eq 1 then begin
	    r_m[0]= xmin
	    r_s= alog10(r_s)
	endif

        for i=0,bins do begin
           lg_r = i*binsize + xmin
           ;sm_r = (i-1)*binsize + xmin
           ;r_m(i-1) = 0.5*(lg_r + sm_r)
	   r_m(i) = lg_r

           ;idx= where((r_s ge sm_r) and (r_s lt lg_r))
	   idx= where(r_s lt lg_r)
           if idx(0) lt 0 then m_in= [0.0] else m_in= mass(idx)

           mass_inr(i) = total(m_in)
        endfor

	print, "mass= ", mass_inr(bins), " at R= ", r_m(bins)

	;idx= where(mass_sd gt 0)
	;if idx(0) ne -1 then begin
	;	mass_sd= mass_sd(idx)
	;	r_m= r_m(idx)
	;endif

	;mass_sd= alog10(mass_sd) - 6.0     ; change to m_solar pc-2

end

