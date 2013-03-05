pro sfr_nofb, junk, $
		filename=filename, $
		cumulative=cumulative, $
		gasmass=gasmass, $
		h=h


if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_nofb, junk, filename=filename, /h"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12.0, newysize=20.0


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"

xaxistitle = "!6Time (Gyr)"

xmax = 4.25
xmin = 0

ymax = 300
ymin= 0.001


;---------------------------

x0= 0.18
x1= 0.98

y0= 0.08
y1= 0.53
y2= 0.98

;---------------------------

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	ytickformat='exp_label', $
	xtickformat='(a1)', $
	;xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata





;-------------------------------------------

;process_one_sfr, "sbw/sb10", lcolor=200, ctab= 1, lthick= 2.0, msg='no BH', x0= 0.75, y0= 0.88, h=0.7
process_one_sfr, "sbw/sb10BH", lcolor=50, lthick= 1.0, msg=' ', x0= 0.75, y0= 0.84, h=0.7
process_one_sfr, "sbw/sb10BHnofb", lcolor=150, lthick= 3.0, lstyle= 1, msg=' ', x0= 0.75, y0= 0.80, h=0.7


xyouts, 0.66, 0.90, 'std BH fb', color= 50, size=2.0, charthick=3.4, /normal
xyouts, 0.66, 0.85, 'no BH fb', color= 150, size=2.0, charthick=3.4, /normal

; ---------------------------------------------------



yaxistitle="!6Black Hole Mass (M!D!9n!6!N)"

ymax= 7.0e+10
ymin= 7.0e+4


!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, /noerase, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;open_blackholes_file, "sbw/sb10", bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass
;bhtime= bhtime/0.7
;bh_mass= 1.0e+10*bh_mass / 0.7
;oplot, bhtime, bh_mass, psym=-3, color= zcolor, linestyle= 0, thick= 12.0

open_blackholes_file, "sbw/sb10BH", bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass
bhtime= bhtime/0.7
bh_mass= 1.0e+10*bh_mass / 0.7
oplot, bhtime, bh_mass, psym=-3, color= 50, linestyle= 0, thick= 4.0

open_blackholes_file, "sbw/sb10BHnofb", bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass
bhtime= bhtime/0.7
bh_mass= 1.0e+10*bh_mass / 0.7
oplot, bhtime, bh_mass, psym=-3, color= 150, linestyle= 1, thick= 12.0



;   done
;  --------

device, /close


print, "done"



end


; =================================================================================


