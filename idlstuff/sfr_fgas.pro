;======================================================================
;
;
;
;      |----------------------|
;      |                      |
;      |                      |
;      |                      |
;      |      SFR             |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;      |                      |
;      |                      |
;      |                      |
;      |      Gas Mass        |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;
;
;
;
;======================================================================


pro sfr_fgas, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_fgas, junk"
   print, "  "
   print, "  "
   return
endif


filename='sfr_fgas.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=10, newysize=16




;--------------------------------------
;--------------------------------------

;xaxistitle = "Time (Gyr h!E-1!D)"
xaxistitle = "!6Time (Gyr)"
xmax = 3.25
xmin = 0


x0= 0.20
x1= 0.97

y0= 0.10
y1= 0.545
y2= 0.99




;--------------------------------------
; Top Plot
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
;ymax = 980.0
ymax = 80.0
ymin = 4.0e-1
;ymin = 2.0e-3

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, $
        xtickformat='(a1)', ytitle=yaxistitle, /ylog, ytickformat='exp_label'



        ;-------------------------------------------
        ;
	;  Show 3 SFR's
        ;
        ;-------------------------------------------
        ;process_and_plot_sfr_history, 'bs/b3e', lcolor=50,  lthick= 2.0, msg='80%', x0= 0.82, y0= 0.92, h=0.7
        ;process_and_plot_sfr_history, 'ds/d3e7', lcolor=100,  lthick= 2.0, msg='40%', x0= 0.82, y0= 0.88, h=0.7
        ;process_and_plot_sfr_history, 'es/e3e', lcolor=150,  lthick= 2.0, msg='20%', x0= 0.82, y0= 0.84, h=0.7

        ;process_and_plot_sfr_history, 'ds/vc3vc3f', lcolor=50,  lthick= 2.0, msg='w/ BH', x0= 0.32, y0= 0.74, h=0.7
        ;process_and_plot_sfr_history, 'ds/vc3vc3f_no', lcolor=150,  lthick= 2.0, msg='w/o BH', x0= 0.32, y0= 0.68, h=0.7
        process_and_plot_sfr_history, 'ds/vc3vc3e', lcolor=50,  lthick= 2.0, msg='w/ BH', x0= 0.32, y0= 0.74, h=0.7
        process_and_plot_sfr_history, 'ds/vc3vc3e_no', lcolor=150,  lthick= 2.0, msg='w/o BH', x0= 0.32, y0= 0.68, h=0.7



;--------------------------------------
; Bottom Plot
;--------------------------------------

;yaxistitle="!6Gas Fraction (M!Dgas!N/ M!Dgas+stellar disk!N)"
yaxistitle="!6M!Dgas!N (t) / M!Dgas+stellar disk!N (t=0)"
;yaxistitle="!6Gas Fraction"
ymax = 0.95
ymin = 0

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle


        ;-------------------------------------------
        ;
	;  Show 3 Gas Fractions
        ;
        ;-------------------------------------------
        ;process_and_plot_sfr_history, 'bs/b3e', lcolor=50,  lthick= 2.0, msg=' ', x0= 0.25, y0= 0.90, gasfraction= 0.8, h=0.7
        ;process_and_plot_sfr_history, 'ds/d3e7', lcolor=100,  lthick= 2.0, msg=' ', x0= 0.25, y0= 0.86, gasfraction= 0.4, h=0.7
        ;process_and_plot_sfr_history, 'es/e3e', lcolor=150,  lthick= 2.0, msg=' ', x0= 0.25, y0= 0.82, gasfraction= 0.2, h=0.7

        ;process_and_plot_sfr_history, 'ds/vc3vc3f', lcolor=50,  lthick= 2.0, msg=' ', x0= 0.25, y0= 0.90, gasfraction= 0.8, h=0.7
        ;process_and_plot_sfr_history, 'ds/vc3vc3f_no', lcolor=150,  lthick= 2.0, msg=' ', x0= 0.25, y0= 0.90, gasfraction= 0.8, h=0.7
        ;process_and_plot_sfr_history, 'ds/vc3vc3e', lcolor=50,  lthick= 2.0, msg=' ', x0= 0.25, y0= 0.90, gasfraction= 0.8, h=0.7
        ;process_and_plot_sfr_history, 'ds/vc3vc3e_no', lcolor=150,  lthick= 2.0, msg=' ', x0= 0.25, y0= 0.90, gasfraction= 0.8, h=0.7



;--------------------------------------
;--------------------------------------




; done
; ----------
device, /close

end









