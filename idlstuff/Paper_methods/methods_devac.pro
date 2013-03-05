pro methods_devac_allofthem, junk


sendto='ps'


;  n=2
; -----
frun= "execute/Sbc201a-u57" & msg='!8n2low!6'
methods_devac_each, frun, msg, sendto, filename='fig10g.eps'

frun= "/data/tcox/Sbc201a-u4" & msg='!8n2med!6'
methods_devac_each, frun, msg, sendto, filename='fig10h.eps'

frun= "/data/tcox/Sbc201a-u43" & msg='!8n2high!6'
methods_devac_each, frun, msg, sendto, filename='fig10i.eps', /includecollisionless


;  n=1
; -----
frun= "execute/Sbc201a-u56" & msg='!8n1low!6'
methods_devac_each, frun, msg, sendto, filename='fig10d.eps'

frun= "/data/tcox/Sbc201a-u53" & msg='!8n1med!6'
methods_devac_each, frun, msg, sendto, filename='fig10e.eps'

;frun= "execute/Sbc201a-u52" & msg='n1high'
frun= "/data/tcox/Sbc201a-u52" & msg='!8n1high!6'
methods_devac_each, frun, msg, sendto, filename='fig10f.eps'


;  n=0
; -----
;frun= "execute/Sbc201a-u50" & msg='n0med'
frun= "execute/Sbc201a-u55" & msg='!8n0low!6'
methods_devac_each, frun, msg, sendto, filename='fig10a.eps'

frun= "/data/tcox/Sbc201a-u50" & msg='!8n0med!6'
methods_devac_each, frun, msg, sendto, filename='fig10b.eps'

;frun= "execute/Sbc201a-u51" & msg='n0high'
frun= "/data/tcox/Sbc201a-u51" & msg='!8n0high!6'
methods_devac_each, frun, msg, sendto, filename='fig10c.eps'






end



pro methods_devac_Sc, junk


sendto='ps'


;  n=2
; -----
frun= "execute/Sc201-u4" & msg='Sc'
methods_devac_each, frun, msg, sendto, filename='devac_Sc.eps'


end



pro methods_devac_timecourse, junk

sendto='ps'


;  n2low
; ---------
frun= "/data/tcox/Sbc201a-u4" & msg='Sbc'
methods_devac_each, frun, msg, sendto, filename='devac_1.eps', snapnum=40
methods_devac_each, frun, msg, sendto, filename='devac_2.eps', snapnum=45
methods_devac_each, frun, msg, sendto, filename='devac_3.eps', snapnum=50
methods_devac_each, frun, msg, sendto, filename='devac_4.eps', snapnum=55
methods_devac_each, frun, msg, sendto, filename='devac_5.eps', snapnum=60
        
        
end     







pro methods_devac_each, frun, msg, sendto, filename=filename, $
				includecollisionless=includecollisionless, $
				snapnum=snapnum

if not keyword_set(sendto) then begin
   print, "  "
   print, "  "
   print, "methods_devac_each, sendto"
   print, "  "
   print, "  "
   return
endif


initialize_plotinfo, 1

;setup_plot_stuff, sendto, filename=filename
setup_plot_stuff, sendto, filename=filename, colortable= 4
; should I try and send a different size for this?



;-------------
;  Plot axes 
;-------------

xmax = 2.4
xmin = 0.3
ymax = 7.0
ymin = 0.0

xaxistitle="!6R !E1/4 !N(kpc!E1/4!N)"
yaxistitle="!6Log !4R!6(r) !N(M!D!9n!6!N pc!E-2!N)"


logplot = 0
if logplot eq 1 then begin
        xmax = alog10(xmax)
        xmin = alog10(xmin)
        ymax = alog10(ymax)
        ymin = alog10(ymin)
        xaxistitle=""
        yaxistitle=""
endif


!p.position= [0.16, 0.14, 0.99, 0.98]

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



if not keyword_set(snapnum) then snapnum= 60



;  n=2
; -----
;frun= "/data/tcox/Sbc201a-u4" & msg='n2low'
;frun= "/data/tcox/Sbc201a-u43" & msg='n2high'

;  n=1
; -----
;frun= "/data/tcox/Sbc201a-u40" & msg='n1low'
;frun= "execute/Sbc201a-u46" & msg='n1high'

;  n=0
; -----
;frun= "/data/tcox/Sbc201a-u8" & msg='n0low'
;frun= "execute/Sbc201a-u45" & msg='n0high'



ok=fload_snapshot(frun,snapnum)


process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax 
oplot, r_devac, mass_sd, psym=-3, linestyle=0, color= 50, thick=6.0

process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax, newstars=1
oplot, r_devac, mass_sd, psym=-3, linestyle=1, color= 150, thick=4.0



; add on collisionless one?
if keyword_set(includecollisionless) then begin
	frun= "/data/tcox/Sbc201a2-u4"
	ok=fload_snapshot(frun,snapnum)
	process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax
	oplot, r_devac, mass_sd, psym=-3, linestyle=2, color= 0, thick=3.0
endif

; add high resolution case
frun= "execute/Sbc201a10x-u4"
ok=fload_snapshot(frun,snapnum)
process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax
oplot, r_devac, mass_sd, psym=-3, linestyle=3, color= 0, thick=8.0

; add N_g= 10 case
frun= "/data/tcox/Sbc201a-u4a"
ok=fload_snapshot(frun,snapnum)
process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax
oplot, r_devac, mass_sd, psym=-3, linestyle=2, color= 0, thick=3.0




xyouts, 0.65, 0.8, msg, /normal, charthick=3.0, size= 2.2, color= 0


; print arrow for smoothing length
smoothlen= 0.1
xarrow= smoothlen^(0.25)
arrow, xarrow, 2.0, xarrow, 1.0, /data, COLOR=0, THICK=4.0, hthick=4.0


;--------------------------------------
if (sendto EQ 'ps') then device, /close


end






;--------------------------------------
; Actually process the devac point
;--------------------------------------
pro process_devac, frun, snapnum, r_devac, mass_sd, xmin, xmax, newstars=newstars

	bins= 30
	binsize= (xmax-xmin)/bins
	r_devac= fltarr(bins)
	mass_sd= fltarr(bins)

	;ok=fload_snapshot(frun,snapnum)

	if keyword_set(newstars) then begin
	        rxy = fload_newstars_xyz('rxy')
	        mass = fload_newstars_mass(1)
	endif else begin
		;rxy = fload_luminousbaryons_xyz('rxy')
		;mass = fload_luminousbaryons_mass(1)
		rxy = fload_allstars_xyz('rxy')
		mass = fload_allstars_mass(1)
	endelse

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

