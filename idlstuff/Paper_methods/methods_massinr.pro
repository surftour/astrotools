pro methods_massinr, sendto, filename=filename

if not keyword_set(sendto) then begin
   print, "  "
   print, "  "
   print, "methods_massinr, sendto, filename=filename"
   print, "  "
   print, "  "
   return
endif


initialize_plotinfo, 1

;setup_plot_stuff, sendto, filename=filename
setup_plot_stuff, sendto, filename=filename, colortable= 4




snapnum= 60

;normalized= 1
normalized= 0


rlog= 1
;rlog= 0




;-------------
;  Plot axes 
;-------------

;if rlog eq 1 then xmax = alog10(50.0) else xmax = 25.0
if rlog eq 1 then xmax = 3.2 else xmax = 25.0
if rlog eq 1 then xmin = -1.5 else xmin = 0.0
if normalized eq 1 then ymax = 1.0 else ymax = 7.8
ymin = 0.0

if rlog eq 1 then xaxistitle='!6Log R (kpc)' else xaxistitle="R (kpc)"
yaxistitle="!6New Stellar Mass (r<R) (10!E10!N M!D!9n!6!N)"


logplot = 0
if logplot eq 1 then begin
        xmax = alog10(xmax)
        xmin = 0.1
        xaxistitle="Log R (kpc)"
endif


!p.position= [0.15, 0.14, 0.98, 0.98]

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



;  n=2
; -----
frun= "execute/Sbc201a-u57"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-2, linestyle=0, color= 150, thick=4.0, symsize=1.5

symsize= 1.5
;usersym,symsize*[-1,-1,1,1],symsize*[-1,1,1,-1],/fill
frun= "/data/tcox/Sbc201a-u4"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-3, linestyle=0, color= 150, thick=5.0

frun= "/data/tcox/Sbc201a-u43"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-5, linestyle=0, color= 150, thick=4.0, symsize=1.5



;  n=1
; -----
frun= "execute/Sbc201a-u56"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-2, linestyle=0, color= 100, thick=4.0, symsize=1.5

;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0   ;,/fill
frun= "/data/tcox/Sbc201a-u40"
frun= "execute/Sbc201a-u53"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-3, linestyle=1, color= 100, thick=5.0

frun= "/data/tcox/Sbc201a-u41"
frun= "execute/Sbc201a-u46"
frun= "execute/Sbc201a-u52"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-5, linestyle=1, color= 100, thick=4.0, symsize=1.5




;  n=0
; -----
frun= "execute/Sbc201a-u55"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-2, linestyle=0, color= 50, thick=4.0, symsize=1.5

;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=3.0, /fill
frun= "/data/tcox/Sbc201a-u8"
frun= "execute/Sbc201a-u50"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-3, linestyle=2, color= 50, thick=5.0

frun= "/data/tcox/Sbc201a-u42"
frun= "execute/Sbc201a-u45"
frun= "execute/Sbc201a-u51"
process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax, totalmass=totalmass
if normalized eq 1 then mass_inr= mass_inr/totalmass
oplot, r_m, mass_inr, psym=-5, linestyle=2, color= 50, thick=4.0, symsize=1.5




; low feedback
; ----------------------------
;  n=2
; -----
;frun= "/data/tcox/Sbc201a-u4"
;process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax
;oplot, r_m, mass_inr, psym=-3, linestyle=0, color= 150, thick=2.0
;  n=1
; -----
;frun= "/data/tcox/Sbc201a-u40"
;process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax
;oplot, r_m, mass_inr, psym=-3, linestyle=0, color= 100, thick=2.0
;  n=0
; -----
;frun= "/data/tcox/Sbc201a-u8"
;process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax
;oplot, r_m, mass_inr, psym=-3, linestyle=0, color= 50, thick=2.0






; high feedback
; ----------------------------
;  n=2
; -----
;frun= "/data/tcox/Sbc201a-u43"
;process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax
;oplot, r_m, mass_inr, psym=-3, linestyle=1, color= 150, thick=2.0
;  n=1
; -----
;frun= "/data/tcox/Sbc201a-u41"
;process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax
;oplot, r_m, mass_inr, psym=-3, linestyle=1, color= 100, thick=2.0
;  n=0
; -----
;frun= "/data/tcox/Sbc201a-u42"
;process_massinr, frun, snapnum, r_m, mass_inr, xmin, xmax
;oplot, r_m, mass_inr, psym=-3, linestyle=1, color= 50, thick=2.0




; print arrow for smoothing length
smoothlen= 0.1
if rlog eq 1 then smoothlen= alog10(smoothlen)
toppos= 2.8
botpos= 2.1
if normalized eq 1 then begin
	toppos= 0.45
	botpos= 0.30
endif
xarrow= smoothlen
arrow, xarrow, toppos, xarrow, botpos, /data, COLOR=0, THICK=4.0, hthick=4.0





display_legend= 1
;display_legend= 0
if display_legend eq 1 then begin
        oplot, [-1.1], [6.2], symsize= 1.5, psym=2, linestyle=0, color= 0, thick=2.0
        oplot, [-1.25, -0.9], [6.2,6.2], psym=-3, linestyle=0, color= 0, thick=4.0

        ;oplot, [-0.1], [-1.4], psym=3, symsize=1.5, linestyle=0, color=0, thick=2.0
        oplot, [-1.25, -0.9], [5.5,5.5], psym=-3, linestyle=0, color= 0, thick=4.0

        oplot, [-1.1], [4.8], psym=5, symsize=1.5, linestyle=0, color=0, thick=2.0
        oplot, [-1.25, -0.9], [4.8,4.8], psym=-3, linestyle=0, color= 0, thick=4.0

        xyouts, -0.7, 6.1, '!8low!6', color= 0, charthick=4.0, size=1.53

        xyouts, -0.7, 5.4, '!8med!6', color= 0, charthick=4.0, size=1.53

        xyouts, -0.7, 4.7, '!8high!6', color= 0, charthick=4.0, size=1.53
endif







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

