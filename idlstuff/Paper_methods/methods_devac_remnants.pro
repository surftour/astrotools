pro methods_devac_remnants, sendto, filename=filename

if not keyword_set(sendto) then begin
   print, "  "
   print, "  "
   print, "methods_devac_remnants, sendto"
   print, "  "
   print, "  "
   return
endif


initialize_plotinfo, 1

setup_plot_stuff, sendto, filename=filename
; should I try and send a different size for this?



;-------------
;  Plot axes 
;-------------

xmax = 2.4
xmin = 0.3
ymax = 7.0
ymin = 0.0

xaxistitle="R !E1/4 !N(kpc!E1/4!N)"
yaxistitle="Log !4R!3(r) !N(M!D!9n!3!N pc!E-2!N)"


logplot = 0
if logplot eq 1 then begin
        xmax = alog10(xmax)
        xmin = alog10(xmin)
        ymax = alog10(ymax)
        ymin = alog10(ymin)
        xaxistitle=""
        yaxistitle=""
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



snapnum= 60



;  n=2
; -----
symsize= 1.5
;usersym,symsize*[-1,-1,1,1],symsize*[-1,1,1,-1],/fill
frun= "/data/tcox/Sbc201a-u4"
process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax
oplot, r_devac, mass_sd, psym=-3, linestyle=0, color= 150, thick=2.0

frun= "/data/tcox/Sbc201a-u43"
process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax
oplot, r_devac, mass_sd, psym=-3, linestyle=0, color= 150, thick=2.0



;  n=1
; -----
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0   ;,/fill
frun= "/data/tcox/Sbc201a-u40"
process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax
oplot, r_devac, mass_sd, psym=-3, linestyle=0, color= 100, thick=2.0

frun= "/data/tcox/Sbc201a-u41"
process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax
oplot, r_devac, mass_sd, psym=-3, linestyle=0, color= 100, thick=2.0




;  n=0
; -----
;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=3.0, /fill
frun= "/data/tcox/Sbc201a-u8"
process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax
oplot, r_devac, mass_sd, psym=-3, linestyle=0, color= 50, thick=2.0

frun= "/data/tcox/Sbc201a-u42"
process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax
oplot, r_devac, mass_sd, psym=-3, linestyle=0, color= 50, thick=2.0






;--------------------------------------
if (sendto EQ 'ps') then device, /close


end






;--------------------------------------
; Actually process the devac point
;--------------------------------------
pro process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax

	bins= 100
	binsize= (xmax-xmin)/bins
	r_devac= fltarr(bins)
	mass_sd= fltarr(bins)

	ok=fload_snapshot(frun,snapnum)

	;rxy = fload_luminousbaryons_xyz('rxy')
	;mass = fload_luminousbaryons_mass(1)
        rxy = fload_stars_xyz('rxy')
        mass = fload_stars_mass(1)

	r_s= rxy^(0.25)

        for i=1,bins do begin
           lg_r = i*binsize + xmin
           sm_r = (i-1)*binsize + xmin
           r_devac(i-1) = 0.5*(lg_r + sm_r)

           idx= where((r_s ge sm_r) and (r_s lt lg_r))
           if idx(0) lt 0 then m_inannulus= [0.0] else m_inannulus= mass(idx)

           mass_tot = total(m_inannulus)
           sd_kpc2 = (3.14159)*( (lg_r)^8 - (sm_r)^8 )    ; from r^1/4 to r^2
           mass_sd(i-1) = mass_tot*(1e10)/sd_kpc2         ; m_solar/ kpc^2
        endfor


	idx= where(mass_sd gt 0)
	if idx(0) ne -1 then begin
		mass_sd= mass_sd(idx)
		r_devac= r_devac(idx)
	endif

	mass_sd= alog10(mass_sd) - 6.0     ; change to m_solar pc-2

end

