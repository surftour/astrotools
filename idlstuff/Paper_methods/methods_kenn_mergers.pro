pro methods_kenn_mergers, sendto, filename=filename

if not keyword_set(sendto) then begin
   print, "  "
   print, "  "
   print, "methods_kenn_mergers, sendto"
   print, "  "
   print, "  "
   return
endif


initialize_plotinfo, 1

setup_plot_stuff, sendto, filename=filename, colortable= 4



;-------------
;  Plot axes 
;-------------

xmax = 3.16e5      ; for direct comparison to kennicutt (or volker for that matter)
xmin = 0.316
ymax = 3.16e3
ymin = 3.16e-5


logplot = 1
if logplot eq 1 then begin
        xmax = alog10(xmax)
        xmin = alog10(xmin)
        ymax = alog10(ymax)
        ymin = alog10(ymin)
        xaxistitle="Log !4R !D!3gas !N(M!D!9n!3!N pc!E-2!N)"
        yaxistitle="Log !4R !D!3SFR !N(M!D!9n!3!N yr!E-1!N kpc!E-2!N)"
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


        if logplot eq 1 then ok=oplot_kenn_log(sendto,/linear) else ok = oplot_kenn_log(sendto)



;-----------------------
;  Put data on graph
;-----------------------



; all iso data
; --------------
oplot, [1.33262], [-1.59533], psym=7, linestyle=0, color= 0, thick=2.0         ; Sbc11i4-u4
oplot, [1.51443], [-1.84801], psym=7, linestyle=0, color= 0, thick=2.0         ; Sbc11i4-u43
oplot, [1.28718], [-1.59880], psym=7, linestyle=0, color= 0, thick=2.0         ; Sbc11i4-u53
oplot, [1.50097], [-1.80059], psym=7, linestyle=0, color= 0, thick=2.0         ; Sbc11i4-u52
oplot, [1.23101], [-1.65817], psym=7, linestyle=0, color= 0, thick=2.0         ; Sbc11i4-u50
oplot, [1.48304], [-1.52666], psym=7, linestyle=0, color= 0, thick=2.0         ; Sbc11i4-u51



;snapnum= 40
;snapnum= 38
;snapnum= 37
;snapnum= 36
snapnum= 35

r_kenn= 2.0


;  n=2
; -----
symsize= 1.5
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
frun= "/data/tcox/Sbc201a-u4"
process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 150, thick=2.0

usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
frun= "/data/tcox/Sbc201a-u43"
process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 150, thick=2.0



;  n=1
; -----
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;frun= "/data/tcox/Sbc201a-u40"
frun= "execute/Sbc201a-u53"
process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 100, thick=2.0

;frun= "/data/tcox/Sbc201a-u41"
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
frun= "execute/Sbc201a-u52"
process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 100, thick=2.0




;  n=0
; -----
usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
;frun= "/data/tcox/Sbc201a-u8"
frun= "execute/Sbc201a-u50"
process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 50, thick=2.0

usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
;frun= "/data/tcox/Sbc201a-u42"
frun= "execute/Sbc201a-u51"
process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 50, thick=2.0






;--------------------------------------
if (sendto EQ 'ps') then device, /close


end






;--------------------------------------
; Actually process the kennicutt point
;--------------------------------------
pro process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log

	ok=fload_snapshot(frun,snapnum)

	;c1= fload_1gal_center(1,170000)
	c1= fload_1gal_center(170001,170000)
	;c1= fload_center_alreadycomp(1)

	rxy = fload_gas_xyz('rxy',center=c1)
	;rxy = fload_gas_xyz('rxy',center=[0,0,0])
	;rxy = fload_gas_xyz('rxy')
	sfr = fload_gas_sfr(1)
	gas_mass= fload_gas_mass(1)
	gasmfs_mass= fload_gas_mfs(1)

	idx_within_r_kenn = where(rxy LE r_kenn)
        if idx_within_r_kenn(0) ne -1 then begin
		sfr_within_r_kenn = sfr(idx_within_r_kenn)
		mass_within_r_kenn = gas_mass(idx_within_r_kenn)-gasmfs_mass(idx_within_r_kenn)

		sd_kpc2 = !PI * r_kenn * r_kenn
		sd_pc2 = sd_kpc2 * 1e6
		gas_sd = total(mass_within_r_kenn)*(1e10)/sd_pc2    ; units Msolar/pc2
		sfr_sd = total(sfr_within_r_kenn)/sd_kpc2           ;   "   Msolar/kpc2 

		gas_sd_log= [alog10(gas_sd)]
		sfr_sd_log= [alog10(sfr_sd)]
	endif else begin
		gas_sd_log= [0]
		sfr_sd_log= [0]
	endelse


	print, "gas_sd, sfr_sd = ", gas_sd_log, sfr_sd_log


end










;===========================
;  Now for c_star mergers
;===========================

pro methods_kenn_mergers_cstar, sendto, filename=filename

if not keyword_set(sendto) then begin
   print, "  "
   print, "  "
   print, "methods_kenn_mergers_cstar, sendto"
   print, "  "
   print, "  "
   return
endif


initialize_plotinfo, 1

setup_plot_stuff, sendto, filename=filename, colortable= 4



;-------------
;  Plot axes 
;-------------

xmax = 3.16e5      ; for direct comparison to kennicutt (or volker for that matter)
xmin = 0.316
ymax = 3.16e3
ymin = 3.16e-5


logplot = 1
if logplot eq 1 then begin
        xmax = alog10(xmax)
        xmin = alog10(xmin)
        ymax = alog10(ymax)
        ymin = alog10(ymin)
        xaxistitle="!6Log !7R !D!6gas !N(M!D!9n!6!N pc!E-2!N)"
        yaxistitle="!6Log !7R !D!6SFR !N(M!D!9n!6!N yr!E-1!N kpc!E-2!N)"
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


        if logplot eq 1 then ok=oplot_kenn_log(sendto,/linear) else ok = oplot_kenn_log(sendto)



;-----------------------
;  Put data on graph
;-----------------------



; all iso data
; --------------
;oplot, [1.33262], [-1.59533], psym=7, linestyle=0, color= 0, thick=2.0         ; Sbc11i4-u4
;oplot, [1.51443], [-1.84801], psym=7, linestyle=0, color= 0, thick=2.0         ; Sbc11i4-u43
;oplot, [1.28718], [-1.59880], psym=7, linestyle=0, color= 0, thick=2.0         ; Sbc11i4-u53
;oplot, [1.50097], [-1.80059], psym=7, linestyle=0, color= 0, thick=2.0         ; Sbc11i4-u52
;oplot, [1.23101], [-1.65817], psym=7, linestyle=0, color= 0, thick=2.0         ; Sbc11i4-u50
;oplot, [1.48304], [-1.52666], psym=7, linestyle=0, color= 0, thick=2.0         ; Sbc11i4-u51



;snapnum= 40
snapnum= 38
;snapnum= 37
;snapnum= 36

r_kenn= 2.0


;  c_star= 0.3
; --------------
;symsize= 2.5
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
frun= "execute/Sbc201a-u74"
process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
;oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 100, thick=2.0
oplot, gas_sd_log, sfr_sd_log, psym=2, linestyle=0, color= 100, thick=2.0


;  c_star= 0.06
; --------------
;symsize= 2.00
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
frun= "execute/Sbc201a-u64"
process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
;oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 100, thick=2.0
oplot, gas_sd_log, sfr_sd_log, psym=2, linestyle=0, color= 100, thick=2.0


;  c_star= 0.03
; --------------
symsize= 1.5
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
frun= "/data/tcox/Sbc201a-u4"
process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 150, thick=2.0


;  c_star= 0.015
; --------------
;symsize= 1.25
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
frun= "execute/Sbc201a-u44"
process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
;oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 50, thick=2.0
oplot, gas_sd_log, sfr_sd_log, psym=7, linestyle=0, color= 50, thick=2.0


;  c_star= 0.004
; --------------
;symsize= 1.0
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
frun= "execute/Sbc201a-u54"
process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
;oplot, gas_sd_log, sfr_sd_log, psym=8, linestyle=0, color= 50, thick=2.0
oplot, gas_sd_log, sfr_sd_log, psym=7, linestyle=0, color= 50, thick=2.0





; Legend
; ---------
symsize= 1.0
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
oplot, [3.6], [-2.2], psym=8, linestyle=0, color= 50, thick=2.0  ; 0.3
oplot, [3.6], [-2.6], psym=7, linestyle=0, color= 100, thick=2.0  ; 0.06
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
oplot, [3.6], [-3.0], psym=8, linestyle=0, color= 150, thick=2.0   ; 0.03, the fiducial
oplot, [3.6], [-3.4], psym=2, linestyle=0, color= 200, thick=2.0  ; 0.015
usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
oplot, [3.6], [-3.8], psym=8, linestyle=0, color= 0, thick=2.0   ; 0.004

xyouts, 3.9, -1.85, '!6c!D*!N', color= 0, charthick=3.0, size=1.33
oplot, [3.9,4.3],[-1.95,-1.95], psym=-3, linestyle=0, color=0, thick=3.0
xyouts, 3.9, -2.3, '0.3', color= 50, charthick=3.0, size=1.33
xyouts, 3.9, -2.7, '0.15', color= 100, charthick=3.0, size=1.33
xyouts, 3.9, -3.1, '0.03', color= 150, charthick=3.0, size=1.33
xyouts, 3.9, -3.5, '0.015', color= 200, charthick=3.0, size=1.33
xyouts, 3.9, -3.9, '0.004', color= 0, charthick=3.0, size=1.33





;--------------------------------------
if (sendto EQ 'ps') then device, /close


end








;=========================
; Time evolution
;=========================

pro methods_kenn_mergers_time, sendto, filename=filename

if not keyword_set(sendto) then begin
   print, "  "
   print, "  "
   print, "methods_kenn_mergers_time, sendto, filename=filename"
   print, "  "
   print, "  "
   return
endif


initialize_plotinfo, 1

setup_plot_stuff, sendto, filename=filename, colortable= 4



;-------------
;  Plot axes 
;-------------

xmax = 3.16e5      ; for direct comparison to kennicutt (or volker for that matter)
xmin = 0.316
ymax = 3.16e3
ymin = 3.16e-5


logplot = 1
if logplot eq 1 then begin
        xmax = alog10(xmax)
        xmin = alog10(xmin)
        ymax = alog10(ymax)
        ymin = alog10(ymin)
        xaxistitle="Log !4R !D!3gas !N(M!D!9n!3!N pc!E-2!N)"
        yaxistitle="Log !4R !D!3SFR !N(M!D!9n!3!N yr!E-1!N kpc!E-2!N)"
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


        if logplot eq 1 then ok=oplot_kenn_log(sendto,/linear) else ok = oplot_kenn_log(sendto)



;-----------------------
;  Put data on graph
;-----------------------


r_kenn= 2.0

frun= "/data/tcox/Sbc201a-u4"

gas_sd= fltarr(30)
sfr_sd= fltarr(30)

for i=0,29 do begin

	snapnum= i*2 + 1

	process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
	gas_sd[i]= gas_sd_log
	sfr_sd[i]= sfr_sd_log


endfor



; now plot it up
; ----------------

; first prior to first passage
starti= 0
;endi= 10
endi= 5
;symsize= 1.5
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
;oplot, gas_sd[starti:endi], sfr_sd[starti:endi], psym=8, linestyle=0, color= 150, thick=2.0
oplot, gas_sd[starti:endi], sfr_sd[starti:endi], psym=7, linestyle=0, color= 0, thick=2.0


; between first passage and final merger
starti= 6   ; 11
endi= 17    ; 35
symsize= 1.5
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0, /fill
oplot, gas_sd[starti:endi], sfr_sd[starti:endi], psym=8, linestyle=0, color= 150, thick=2.0



; finally for the merger remnant
starti= 18   ; 36
endi= 29  ;60
symsize= 1.5
;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
;oplot, gas_sd[starti:endi], sfr_sd[starti:endi], psym=7, linestyle=0, color= 50, thick=2.0
oplot, gas_sd[starti:endi], sfr_sd[starti:endi], psym=8, linestyle=0, color= 50, thick=2.0







;--------------------------------------
if (sendto EQ 'ps') then device, /close


end













; this is for the second time sequence graph, fig 15
; in which we do time evol. for all c_star models
;
; --------------------------------------------------------
pro overplot_onesim_timekenn, frun, r_kenn, pointselection, $
			msg=msg, xm=xm, ym=ym


gas_sd= fltarr(20)
sfr_sd= fltarr(20)

; first 8
for i=0,7 do begin
	r_kenn= 2.0
	snapnum= i*4 + 1    ; 1, 5, 9, 13, 17, 21, 25, 29
	process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
	gas_sd[i]= gas_sd_log
	sfr_sd[i]= sfr_sd_log
endfor

; middle 8
for i=8,15 do begin
	r_kenn= 0.5
        snapnum= (i-8) + 32    ;  32, 33, 34 ,35, 36, 37, 38, 39
        process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
        gas_sd[i]= gas_sd_log
        sfr_sd[i]= sfr_sd_log
endfor

; last 4
for i=16,19 do begin
	r_kenn= 1.0
        snapnum= (i-16)*5 + 40    ; 40, 45, 50, 55
        process_kenn_info, frun, snapnum, r_kenn, gas_sd_log, sfr_sd_log
        gas_sd[i]= gas_sd_log
        sfr_sd[i]= sfr_sd_log
endfor







;  open (or filled) square
; --------------------------
if pointselection eq 1 then begin
	symsize= 1.0
	;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
	usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
	symsel= 8
	symcolor= 150
endif

;  x's
; -----
if pointselection eq 2 then begin
	symsel= 7
	symcolor= 180
endif

;  open circle
; --------------
if pointselection eq 3 then begin
	symsize= 1.0
	;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
	symsel= 8
	symcolor= 50
endif

;  open triangle
; ----------------
if pointselection eq 4 then begin
	symsize= 1.0
	;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
	usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
	symsel= 8
	symcolor= 90
endif

;  asterisk
; ----------
if pointselection eq 5 then begin
	symsel= 2
	symcolor= 10
endif


;  plus sign
; ----------
if pointselection eq 6 then begin
        symsel= 1
        symcolor= 30
endif


;  filled small box
; -------------------
if pointselection eq 7 then begin
	symsize= 0.5
	usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0, /fill
	symsel= 8
	symcolor= 110
endif


;  filled circle
; ----------
if pointselection eq 8 then begin
	symsize= 0.8
	usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	symsel= 8
	symcolor= 200
endif


;  diamond
; ----------
if pointselection eq 9 then begin
        symsel= 4
        symcolor= 100
endif


; now plot it up
; ----------------

oplot, gas_sd, sfr_sd, psym=symsel, color= symcolor, thick=2.0


if keyword_set(msg) then begin
	if not keyword_set(xm) then xm= 0.77
	if not keyword_set(ym) then ym= 0.50
	xyouts, xm, ym, msg, color= symcolor, charthick=3.0, size=1.2, /normal

	oplot, [3.6], [-0.4 - 10.0*(0.55-ym)], psym=symsel, color= symcolor, thick=2.0
endif


end

















;==================================
; Time evolution - of c_star runs
;==================================

pro methods_kenn_mergers_cstar_time, sendto, filename=filename

if not keyword_set(sendto) then begin
   print, "  "
   print, "  "
   print, "methods_kenn_mergers_cstar_time, sendto, filename=filename"
   print, "  "
   print, "  "
   return
endif


initialize_plotinfo, 1

setup_plot_stuff, sendto, filename=filename, colortable= 4



;-------------
;  Plot axes 
;-------------

xmax = 3.16e5      ; for direct comparison to kennicutt (or volker for that matter)
xmin = 0.316
ymax = 3.16e3
ymin = 3.16e-5


logplot = 1
if logplot eq 1 then begin
        xmax = alog10(xmax)
        xmin = alog10(xmin)
        ymax = alog10(ymax)
        ymin = alog10(ymin)
        xaxistitle="!6Log !7R !D!6gas !N(M!D!9n!6!N pc!E-2!N)"
        yaxistitle="!6Log !7R !D!6SFR !N(M!D!9n!6!N yr!E-1!N kpc!E-2!N)"
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


        if logplot eq 1 then ok=oplot_kenn_log(sendto,/linear) else ok = oplot_kenn_log(sendto)



;-----------------------
;  Put data on graph
;-----------------------


r_kenn= 1.0


; std, c_star= 0.03
frun= "/data/tcox/Sbc201a-u4"
overplot_onesim_timekenn, frun, r_kenn, 1


; lower sf normalizations
; -------------------------

; std, c_star= 0.015
frun= "execute/Sbc201a-u44"
overplot_onesim_timekenn, frun, r_kenn, 5

; std, c_star= 0.004
frun= "execute/Sbc201a-u54"
overplot_onesim_timekenn, frun, r_kenn, 4


; higher sf normalizations
; -------------------------

; std, c_star= 0.06
frun= "execute/Sbc201a-u64"
overplot_onesim_timekenn, frun, r_kenn, 2

; std, c_star= 0.3
frun= "execute/Sbc201a-u74"
overplot_onesim_timekenn, frun, r_kenn, 3



; Legend
; ---------
symsize= 1.0
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
oplot, [3.68], [-2.2], psym=8, linestyle=0, color= 50, thick=2.0  ; 0.3
oplot, [3.68], [-2.6], psym=7, linestyle=0, color= 100, thick=2.0  ; 0.06
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
oplot, [3.68], [-3.0], psym=8, linestyle=0, color= 150, thick=2.0   ; 0.03, the fiducial
oplot, [3.68], [-3.4], psym=2, linestyle=0, color= 200, thick=2.0  ; 0.015
usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
oplot, [3.68], [-3.8], psym=8, linestyle=0, color= 0, thick=2.0   ; 0.004

xyouts, 3.9, -1.85, '!6c!D*!N', color= 0, charthick=3.0, size=1.33
oplot, [3.9,4.3],[-1.95,-1.95], psym=-3, linestyle=0, color=0, thick=3.0
xyouts, 3.9, -2.3, '0.3', color= 50, charthick=3.0, size=1.33
xyouts, 3.9, -2.7, '0.15', color= 100, charthick=3.0, size=1.33
xyouts, 3.9, -3.1, '0.03', color= 150, charthick=3.0, size=1.33
xyouts, 3.9, -3.5, '0.015', color= 200, charthick=3.0, size=1.33
xyouts, 3.9, -3.9, '0.004', color= 0, charthick=3.0, size=1.33



;--------------------------------------
if (sendto EQ 'ps') then device, /close




end






;=======================================================================







;==================================
; Time evolution of all runs
;==================================

pro methods_kenn_mergers_sffb_time, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "  "
   print, "methods_kenn_mergers_sffb_time, junk"
   print, "  "
   print, "  "
   return
endif


filename='kenn_mergsffb.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4



;-------------
;  Plot axes 
;-------------

xmax = 3.16e5      ; for direct comparison to kennicutt (or volker for that matter)
xmin = 0.316
ymax = 3.16e3
ymin = 3.16e-5


logplot = 1
if logplot eq 1 then begin
        xmax = alog10(xmax)
        xmin = alog10(xmin)
        ymax = alog10(ymax)
        ymin = alog10(ymin)
        xaxistitle="!6Log !7R !D!6gas !N(M!D!9n!6!N pc!E-2!N)"
        yaxistitle="!6Log !7R !D!6SFR !N(M!D!9n!6!N yr!E-1!N kpc!E-2!N)"
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


        if logplot eq 1 then ok=oplot_kenn_log('ps',/linear) else ok = oplot_kenn_log('ps')



;-----------------------
;  Put data on graph
;-----------------------


r_kenn= 0.5


; ------------------------

; n0low
frun= "execute/Sbc201a-u55"
overplot_onesim_timekenn, frun, r_kenn, 1, msg='!8n0low', xm= 0.77, ym= 0.55

; n0med
frun= "/data/tcox/Sbc201a-u50"
overplot_onesim_timekenn, frun, r_kenn, 2, msg='!8n0med', xm= 0.77, ym= 0.52

; n0high
frun= "/data/tcox/Sbc201a-u51"
overplot_onesim_timekenn, frun, r_kenn, 8, msg='!8n0high', xm= 0.77, ym= 0.49

; ------------------------

; n1low
frun= "execute/Sbc201a-u56"
overplot_onesim_timekenn, frun, r_kenn, 3, msg='!8n1low', xm= 0.77, ym= 0.43

; n1med
frun= "/data/tcox/Sbc201a-u53"
overplot_onesim_timekenn, frun, r_kenn, 5, msg='!8n1med', xm= 0.77, ym= 0.40

; n1high
frun= "/data/tcox/Sbc201a-u52"
overplot_onesim_timekenn, frun, r_kenn, 6, msg='!8n1high', xm= 0.77, ym= 0.37

; ------------------------

; n2low
frun= "execute/Sbc201a-u57"
overplot_onesim_timekenn, frun, r_kenn, 4, msg='!8n2low!6', xm= 0.77, ym= 0.31

; std, n2med
frun= "/data/tcox/Sbc201a-u4"
overplot_onesim_timekenn, frun, r_kenn, 7, msg='!8n2med!6', xm= 0.77, ym= 0.28

; n2high
frun= "/data/tcox/Sbc201a-u43"
overplot_onesim_timekenn, frun, r_kenn, 9, msg='!8n2high!6', xm= 0.77, ym= 0.25


; Legend
; ---------
symsize= 1.0
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
;oplot, [3.68], [-2.2], psym=8, linestyle=0, color= 50, thick=2.0  ; 0.3
;oplot, [3.68], [-2.6], psym=7, linestyle=0, color= 100, thick=2.0  ; 0.06
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
;oplot, [3.68], [-3.0], psym=8, linestyle=0, color= 150, thick=2.0   ; 0.03, the fiducial
;oplot, [3.68], [-3.4], psym=2, linestyle=0, color= 200, thick=2.0  ; 0.015
;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
;oplot, [3.68], [-3.8], psym=8, linestyle=0, color= 0, thick=2.0   ; 0.004

;xyouts, 3.9, -1.85, 'c!D*!N', color= 0, charthick=3.0, size=1.33
;oplot, [3.9,4.3],[-1.95,-1.95], psym=-3, linestyle=0, color=0, thick=3.0
;xyouts, 3.9, -2.3, '0.3', color= 50, charthick=3.0, size=1.33
;xyouts, 3.9, -2.7, '0.15', color= 100, charthick=3.0, size=1.33
;xyouts, 3.9, -3.1, '0.03', color= 150, charthick=3.0, size=1.33
;xyouts, 3.9, -3.5, '0.015', color= 200, charthick=3.0, size=1.33
;xyouts, 3.9, -3.9, '0.004', color= 0, charthick=3.0, size=1.33



;--------------------------------------
device, /close




end




