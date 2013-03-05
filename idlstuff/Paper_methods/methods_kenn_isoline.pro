pro methods_kenn_isoline, junk


sendto= 'ps'
filename= 'kenni.eps'

initialize_plotinfo, 1
setup_plot_stuff, sendto, filename=filename, colortable=4


; must run prodir before this, because this routime calls
; process_kenn_line, which 



; ---------------------
;  Set up Parmaters
; ---------------------

xaxistitle="!6Log !7R !D!6gas !N(M!D!9n!6!N pc!E-2!N)"
yaxistitle="!6Log !7R !D!6SFR !N(M!D!9n!6!N yr!E-1!N kpc!E-2!N)"


;xmax = 3.16e5      ; for direct comparison to kennicutt (or volker for that matter)
;xmin = 0.316
;ymax = 3.16e3
;ymin = 3.16e-5

; we'll make the plot log
;xmax = 3.5
;xmax = 2.5
xmax = 2.3
;xmin = -0.5
;xmin = 0.0
xmin = -0.4

;ymax = 1.0
;ymax = 1.5
ymax = 0.0
ymin = -6.0




; ------------------
; Plot this up
; ------------------
!p.position= [0.20, 0.15, 0.95, 0.95]
;!p.font= 0

plot, [10.0], [10.0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




snapnum= 10
;load_and_plot, "execute/I1i-u1", snapnum, plotcolor=150, plotthick=4.0
;load_and_plot, "execute/I1i-u4", snapnum, plotcolor=50, plotthick=4.0
;load_and_plot, "execute/I1i-u5", snapnum, plotcolor=100, plotthick=4.0

;snapnum= 7
;load_and_plot, "execute/I1gf1i-u1", snapnum, plotcolor=150, plotthick=4.0
;load_and_plot, "execute/I1gf1i-u2", snapnum, plotcolor=150, plotthick=4.0
;load_and_plot, "execute/I1gf1i-u4", snapnum, plotcolor=50, plotthick=4.0
;load_and_plot, "execute/I1gf1i-u5", snapnum, plotcolor=100, plotthick=4.0

;load_and_plot, "execute/I1gf2i-u1", snapnum, plotcolor=150, plotthick=4.0
;load_and_plot, "execute/I1gf2i-u4", snapnum, plotcolor=50, plotthick=4.0


; show n (sfr n) dependence
;snapnum= 1
;load_and_plot, "execute/I1i-u5", snapnum, plotcolor=0, plotthick=4.0
;load_and_plot, "execute/I1i-u10", snapnum, plotcolor=140, plotthick=4.0
;load_and_plot, "execute/I1i-u11", snapnum, plotcolor=80, plotthick=4.0
;load_and_plot, "execute/I1i-u12", snapnum, plotcolor=50, plotthick=4.0
;load_and_plot, "execute/I1i-u13", snapnum, plotcolor=150, plotthick=4.0
;load_and_plot, "execute/I1i-u14", snapnum, plotcolor=100, plotthick=4.0
;load_and_plot, "execute/I1i-u15", snapnum, plotcolor=50, plotthick=4.0



; time sequence
;frun= "execute/I1i-u1"
;frun= "execute/I1gf1i-u1"
;frun= "execute/I1gf2i-u1"
;load_and_plot, frun, 1, plotcolor=20, plotthick=4.0
;load_and_plot, frun, 2, plotcolor=40, plotthick=4.0
;load_and_plot, frun, 3, plotcolor=60, plotthick=4.0
;load_and_plot, frun, 4, plotcolor=80, plotthick=4.0
;load_and_plot, frun, 5, plotcolor=100, plotthick=4.0
;load_and_plot, frun, 6, plotcolor=120, plotthick=4.0
;load_and_plot, frun, 7, plotcolor=140, plotthick=4.0
;load_and_plot, frun, 8, plotcolor=160, plotthick=4.0
;load_and_plot, frun, 9, plotcolor=180, plotthick=4.0
;load_and_plot, frun, 10, plotcolor=200, plotthick=4.0



; ---------------------------


do_Is= 1
;do_Is= 0

if do_Is eq 1 then begin

	snapnum= 2


	; without cut-off

	; n= 2
	;load_and_plot, 'execute/I1i-u7a', snapnum, plotcolor=150, plotthick=4.0;, pstyle= 5
	load_and_plot, 'execute/I1i-u1a', snapnum, plotcolor=150, plotthick=4.0;, pstyle= 5
	;load_and_plot, 'execute/I1i-u5a', snapnum, plotcolor=150, plotthick=4.0;, pstyle= 5
	;load_and_plot, 'execute/I1i-u6a', snapnum, plotcolor=150, plotthick=4.0;, pstyle= 5

	; n= 1
	;load_and_plot, 'execute/I1i-u17a', snapnum, plotcolor=100, plotthick=4.0;, pstyle= 4
	load_and_plot, 'execute/I1i-u3a', snapnum, plotcolor=100, plotthick=4.0;, pstyle= 4
	;load_and_plot, 'execute/I1i-u16a', snapnum, plotcolor=100, plotthick=4.0;, pstyle= 4

	; n= 0
	;load_and_plot, 'execute/I1i-u9a', snapnum, plotcolor=50, plotthick=4.0;, pstyle= 2
	load_and_plot, 'execute/I1i-u2a', snapnum, plotcolor=50, plotthick=4.0;, pstyle= 2
	;load_and_plot, 'execute/I1i-u8a', snapnum, plotcolor=50, plotthick=4.0;, pstyle= 2


	; ---------------------------


	; with cut-off

	; n= 2
	load_and_plot, 'execute/I1i-u5', snapnum, plotcolor=150, plotthick=3.0, pstyle= 2
	load_and_plot, 'execute/I1i-u6', snapnum, plotcolor=150, plotthick=3.0, pstyle= 2
	load_and_plot, 'execute/I1i-u1', snapnum, plotcolor=150, plotthick=3.0;, pstyle= 5
	load_and_plot, 'execute/I1i-u7', snapnum, plotcolor=150, plotthick=3.0, pstyle= 5

	; n= 1
	load_and_plot, 'execute/I1i-u16', snapnum, plotcolor=100, plotthick=2.0, pstyle= 2
	load_and_plot, 'execute/I1i-u3', snapnum, plotcolor=100, plotthick=2.0;, pstyle= 4
	load_and_plot, 'execute/I1i-u17', snapnum, plotcolor=100, plotthick=2.0, pstyle= 5

	; n= 0
	load_and_plot, 'execute/I1i-u8', snapnum, plotcolor=50, plotthick=2.0, pstyle= 2
	load_and_plot, 'execute/I1i-u2', snapnum, plotcolor=50, plotthick=2.0;, pstyle= 2
	load_and_plot, 'execute/I1i-u9', snapnum, plotcolor=50, plotthick=2.0, pstyle= 5

endif


; ---------------------------


;load_and_plot, 'execute/Sbc11i4-u70', snapnum, plotcolor=150, plotthick=2.0
;load_and_plot, 'execute/Sbc11i4-u71', snapnum, plotcolor=120, plotthick=2.0
;load_and_plot, 'execute/Sbc11i4-u72', snapnum, plotcolor=100, plotthick=2.0
;load_and_plot, 'execute/Sbc11i4-u73', snapnum, plotcolor=50, plotthick=2.0

;load_and_plot, 'execute/Sbc11i4-u4aa', snapnum, plotcolor=200, plotthick=4.0
;load_and_plot, 'execute/Sbc11i4-u43a', snapnum, plotcolor=180, plotthick=4.0


; ---------------------------


;do_AltNs= 1
do_AltNs= 0

if do_AltNs eq 1 then begin

	load_and_plot, 'execute/I1i-u1a', snapnum, plotcolor=150, plotthick=3.0, pstyle= 10
	load_and_plot, 'execute/I1i-u5a', snapnum, plotcolor=150, plotthick=3.0, pstyle= 10
	load_and_plot, 'execute/I1i-u6a', snapnum, plotcolor=150, plotthick=3.0, pstyle= 10

	load_and_plot, 'execute/I1i-u70', snapnum, plotcolor=110, plotthick=3.0, pstyle= 2
	load_and_plot, 'execute/I1i-u72', snapnum, plotcolor=100, plotthick=3.0, pstyle= 2

	load_and_plot, 'execute/I1i-u71', snapnum, plotcolor=50, plotthick=3.0, pstyle= 5
	load_and_plot, 'execute/I1i-u73', snapnum, plotcolor=50, plotthick=3.0, pstyle= 5

endif


; ---------------------------




;do_Sbcs= 1
do_Sbcs= 0

if do_Sbcs eq 1 then begin

        snapnum= 5


        ; without cut-off

        ; n= 2
        load_and_plot, 'execute/Sbc11i4-u57a', snapnum, plotcolor=150, plotthick=4.0;, pstyle= 5
        load_and_plot, 'execute/Sbc11i4-u4aa', snapnum, plotcolor=150, plotthick=4.0;, pstyle= 5
        load_and_plot, 'execute/Sbc11i4-u43a', snapnum, plotcolor=150, plotthick=4.0;, pstyle= 5

        ; n= 1
        load_and_plot, 'execute/Sbc11i4-u56a', snapnum, plotcolor=100, plotthick=3.0;, pstyle= 4
        load_and_plot, 'execute/Sbc11i4-u53a', snapnum, plotcolor=100, plotthick=3.0;, pstyle= 4
	load_and_plot, 'execute/Sbc11i4-u52a', snapnum, plotcolor=100, plotthick=3.0;, pstyle= 4

        ; n= 0
        load_and_plot, 'execute/Sbc11i4-u55a', snapnum, plotcolor=50, plotthick=3.0;, pstyle= 2
        load_and_plot, 'execute/Sbc11i4-u50a', snapnum, plotcolor=50, plotthick=3.0;, pstyle= 2
        load_and_plot, 'execute/Sbc11i4-u51a', snapnum, plotcolor=50, plotthick=3.0;, pstyle= 2


        ; ---------------------------


        ; with cut-off

        ; n= 2
        load_and_plot, 'execute/Sbc11i4-u57', snapnum, plotcolor=150, plotthick=3.0;, pstyle= 5
        load_and_plot, '/data/tcox/Sbc11i4-u4', snapnum, plotcolor=150, plotthick=3.0;, pstyle= 5
        load_and_plot, 'execute/Sbc11i4-u43', snapnum, plotcolor=150, plotthick=3.0, pstyle= 5

        ; n= 1
        load_and_plot, 'execute/Sbc11i4-u56', snapnum, plotcolor=100, plotthick=2.0;, pstyle= 4
        load_and_plot, 'execute/Sbc11i4-u53', snapnum, plotcolor=100, plotthick=2.0;, pstyle= 4
        load_and_plot, 'execute/Sbc11i4-u52', snapnum, plotcolor=100, plotthick=2.0, pstyle= 4

        ; n= 0
        load_and_plot, 'execute/Sbc11i4-u55', snapnum, plotcolor=50, plotthick=2.0;, pstyle= 2
        load_and_plot, 'execute/Sbc11i4-u50', snapnum, plotcolor=50, plotthick=2.0;, pstyle= 2
        load_and_plot, 'execute/Sbc11i4-u51', snapnum, plotcolor=50, plotthick=2.0, pstyle= 2

endif







;ok=oplot_kenn_log(sendto,/linear)
ok=oplot_kenn_log(sendto,/linear,/ditch_boundary)


;xyouts, 0.25, 0.85, fload_timelbl(1,2), /normal, size=1, color= 0, charthick= 2.5



; 
; Legend
; ---------
;display_sfrlaw_legend= 1
display_sfrlaw_legend= 0
if display_sfrlaw_legend eq 1 then begin
	symsize= 1.5
	oplot, [1.0], [-4.4], symsize= 1.5, psym=2, linestyle=0, color= 100, thick=2.0

	usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0, /fill
	oplot, [1.0], [-4.8], psym=8, linestyle=0, color= 150, thick=2.0

	oplot, [1.0], [-5.2], psym=5, symsize=1.5, linestyle=0, color=50, thick=2.0
	
	xyouts, 0.8, -4.0, '!6SFR !9? !7q!6!EN!N', color= 0, charthick=4.0, size=2.5

	xyouts, 1.2, -4.5, '!6N= 1.0', color= 100, charthick=4.0, size=1.53

	xyouts, 1.2, -4.9, '!6N= 1.5', color= 150, charthick=4.0, size=1.53

	xyouts, 1.2, -5.3, '!6N= 2.0', color= 50, charthick=4.0, size=1.53
endif




display_iso_legend= 1
;display_iso_legend= 0
if display_iso_legend eq 1 then begin
        oplot, [-0.1], [-1.0], symsize= 1.5, psym=2, linestyle=0, color= 0, thick=2.0
	oplot, [-0.25, 0.05], [-1.0,-1.0], psym=-3, linestyle=0, color= 0, thick=4.0

        ;oplot, [-0.1], [-1.4], psym=3, symsize=1.5, linestyle=0, color=0, thick=2.0
	oplot, [-0.25, 0.05], [-1.4,-1.4], psym=-3, linestyle=0, color= 0, thick=4.0

        oplot, [-0.1], [-1.8], psym=5, symsize=1.5, linestyle=0, color=0, thick=2.0
	oplot, [-0.25, 0.05], [-1.8,-1.8], psym=-3, linestyle=0, color= 0, thick=4.0

        xyouts, 0.1, -1.1, '!8low!6', color= 0, charthick=4.0, size=1.53

        xyouts, 0.1, -1.5, '!8med!6', color= 0, charthick=4.0, size=1.53

        xyouts, 0.1, -1.9, '!8high!6', color= 0, charthick=4.0, size=1.53
endif






if sendto eq 'ps' then device, /close




end









pro load_and_plot, frun, snapnum, plotcolor=plotcolor, $
				plotlinest=plotlinest, $
				plotthick=plotthick, $
				pstyle=pstyle

	bins = 40

	plotlinest= 0    ; didn't work that will with other linestyles

	ok=fload_snapshot(frun,snapnum)
	a= fload_gas_xyz('rxy')
	b= fload_gas_sfr(1)
	c= fload_gas_mass(1)
	d= 0.0 + 0.0*c

	process_kenn_line, a, b, c, d, bins, mass_sd, sfr_sd

	mass_sd = alog10(mass_sd)
	sfr_sd = alog10(sfr_sd)

print, "SFR SD", max(sfr_sd), min(sfr_sd)
print, "Mass SD", max(mass_sd), min(mass_sd)

	ppsym= -3
	psymsize= 1.5
	if not keyword_set(pstyle) then pstyle= 0 else ppsym= -pstyle

	if pstyle eq 10 then begin
		ppsym= -8
		psymsize= 2.0
		usersym,[-1,-1,1,1,-1],[-1,1,1,-1,-1],thick=4.0,/fill
	endif


	oplot, mass_sd, sfr_sd, psym=ppsym, linestyle= plotlinest, $
				color= plotcolor, thick= plotthick, $
				symsize= psymsize



end





