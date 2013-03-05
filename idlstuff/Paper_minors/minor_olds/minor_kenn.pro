pro minor_kenn, junk


sendto= 'ps'
filename= 'minkenn.eps'

initialize_plotinfo, 1
setup_plot_stuff, sendto, filename=filename, colortable=4


; must run prodir before this, because this routime calls
; process_kenn_line, which 



; ---------------------
;  Set up Parmaters
; ---------------------

xaxistitle="Log !4R !D!3gas !N(M!D!9n!3!N pc!E-2!N)"
yaxistitle="Log !4R !D!3sfr !N(M!D!9n!3!N yr!E-1!N kpc!E-2!N)"


;xmax = 3.16e5      ; for direct comparison to kennicutt (or volker for that matter)
;xmin = 0.316
;ymax = 3.16e3
;ymin = 3.16e-5

; we'll make the plot log
;xmax = 4.0
xmax = 2.5
;xmin = -1.0
xmin = 0.0
;ymax = 1.0
ymax = 0.0
ymin = -6.0




; ------------------
; Plot this up
; ------------------
!p.position= [0.20, 0.15, 0.95, 0.95]

plot, [10.0], [10.0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata






; ------------------
;  Process - 1
; ------------------

frun= "execute/G3il-u1a"
snapnum= 150
load_and_plot, frun, snapnum, plotcolor=150, plotthick=2.0


; ------------------
;  Process - 2
; ------------------

frun= "execute/G2im-u1a"
load_and_plot, frun, snapnum, plotcolor=100, plotthick=4.0


; ------------------
;  Process - 3
; ------------------

frun= "execute/G1i-u1a"
load_and_plot, frun, snapnum, plotcolor=50, plotthick=7.0


; ------------------
;  Process - 4
; ------------------

frun= "execute/G0i-u1a"
load_and_plot, frun, snapnum, plotcolor=0, plotthick=10.0




ok=oplot_kenn_log(sendto,/linear)
;xyouts, 0.25, 0.85, fload_timelbl(1,2), /normal, size=1, color= 0, charthick= 2.5


if sendto eq 'ps' then device, /close




end









pro load_and_plot, frun, snapnum, plotcolor=plotcolor, plotlinest=plotlinest, plotthick=plotthick

	bins = 50

	plotlinest= 0    ; didn't work that will with other linestyles

	ok=fload_snapshot(frun,snapnum)
	a= fload_gas_xyz('rxy')
	b= fload_gas_sfr(1)
	c= fload_gas_mass(1)
	d= 0.0 + 0.0*c

	process_kenn_line, a, b, c, d, bins, mass_sd, sfr_sd

	mass_sd = alog10(mass_sd)
	sfr_sd = alog10(sfr_sd)


	oplot, mass_sd, sfr_sd, psym=-3, linestyle= plotlinest, color= plotcolor, thick= plotthick



end





