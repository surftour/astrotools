pro plot_sfr, galinfo, xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, $
		filename=filename, $
		ynotlog=ynotlog, $
		h=h


if not keyword_set(galinfo) then begin
   print, "  "
   print, "  "
   print, "  WARNING: check script!  "
   print, "  "
   print, "  default filename: sfr.eps"
   print, "  "
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"

;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7


if not keyword_set(xmax) then begin
;xmax = 14.0
;xmax = 12.0
;xmax = 10.0
;xmax = 7.5
;xmax = 6.0
;xmax = 5.0
;xmax = 4.25
;xmax = 4.0
;xmax = 3.5
;xmax = 3.0
;xmax = 2.8
;xmax = 2.5
;xmax = 2.4
xmax = 2.2
;xmax = 2.0
;xmax = 1.5
;xmax = 1.4
;xmax = 1.3
;xmax = 1.0
;xmax = 0.5
endif 


if not keyword_set(xmin) then xmin = 0

if not keyword_set(ymax) then begin
;ymax= 8e+5
;ymax= 8e+4
;ymax= 2e+4
;ymax = 8000
;ymax = 4000
;ymax = 2500
;ymax = 1750
;ymax = 1500
ymax = 1000
;ymax = 800
;ymax = 600
;ymax = 400
;ymax = 300
;ymax = 250.0
;ymax = 202
;ymax = 180
;ymax = 120
;ymax = 100
;ymax = 90
;ymax = 80
;ymax = 60
;ymax = 50.0
;ymax = 40.0
;ymax = 30.0
;ymax = 20.0
;ymax = 16.0
;ymax = 15
;ymax = 13
;ymax = 12
;ymax = 10.0
;ymax = 8.0
;ymax = 6.0
;ymax = 5.0
;ymax = 3.5
;ymax = 2.0
;ymax = 1.5
;ymax = 1.3
;ymax = 1.2
;ymax = 0.75
;ymax = 0.6
;ymax = 0.5
endif

if not keyword_set(ymin) then begin
;ymin= 8.0
;ymin= 1.0
;ymin = 0  & ynotlog= 1
;ymin = 0.8
;ymin = 0.4
;ymin = 0.1
;ymin= 0.07
ymin= 0.01
;ymin= 0.001
;ymin= 0.0001
;ymin= 0.00001
endif

if keyword_set(ynotlog) then ymin = 0 

; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	;h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;---------------------------

!p.position= [0.16, 0.14, 0.98, 0.98]

if keyword_set(ynotlog) then begin
    plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	;ytickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata
endif else begin
    plot, [1.0],[1.0], psym=-3, $
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
endelse




;-------------------------------------------
;   Load the runs to display
;-------------------------------------------


;n_runs= (size(galinfo))[2]
;for i=0, n_runs-1 do begin
;	process_and_plot_sfr_history, string(galinfo[0,i]), $
;			lcolor=long(galinfo[1,i]), $
;			lthick=float(galinfo[2,i]), $
;			msg=string(galinfo[3,i]), $
;			x0= float(galinfo[4,i]), $
;			y0= float(galinfo[5,i]), h=h


n_runs= (size(galinfo))[1]
for i=0, n_runs-1 do begin

	print, " "
	print, "i= ", i
	print, " "
	process_and_plot_sfr_history, galinfo[i].frun, $
			lcolor=galinfo[i].color, $
			lthick=galinfo[i].lthick, $
			msg=galinfo[i].msg, $
			x0= galinfo[i].x0, $
			y0= galinfo[i].y0, h=h


print, "+++++++++++++++++++++++++++++++++++++++"
print, " "
print, " WARNING: h^-1 conversion is occuring"
print, " "
print, "+++++++++++++++++++++++++++++++++++++++"


endfor



;--------------------------------------
;--------------------------------------

device, /close


end





;=====================================================================
;=====================================================================
;=====================================================================
;=====================================================================


