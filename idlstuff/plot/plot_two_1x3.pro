; =================================================
;
;   Generic plotting program using the following
;   format:
;
;
;
;       ----------------       -----------------
;       |              |       |               |
;       |              |       |               |
;       |              |       |               |
;       |    1         |       |    4          |
;       |              |       |               |
;       |              |       |               |
;       |              |       |               |
;       ----------------       -----------------
;       |              |       |               |
;       |              |       |               |
;       |    2         |       |    5          |
;       |              |       |               |
;       |              |       |               |
;       |              |       |               |
;       |              |       |               |
;       ----------------       -----------------
;       |              |       |               |
;       |              |       |               |
;       |    3         |       |    6          |
;       |              |       |               |
;       |              |       |               |
;       |              |       |               |
;       |              |       |               |
;       ----------------       -----------------
;
;
;
; ==================================================



pro plot_two_1x3, filename=filename, xtitle, xmax, xmin, xlog=xlog, $
        ytitle_1, x_1, y_1, ymax_1, ymin_1, ylog_1=ylog_1



if not keyword_set(filename) then begin
	print, " "
	print, " PROBLEM: plot_two_1x3  "
	print, " "
	return
endif



;--------------------------------------------


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=22, newysize=18

x0= 0.10 & xsize= 0.38
x1= x0+xsize

; put the same x0 inbetween the two columns
x2= 2.*x0+xsize
x3= 2.*x0+xsize+xsize


y0= 0.09 & ysize=0.29
y1= y0+ysize
y2= y0+ysize+ysize
y3= y0+ysize+ysize+ysize




; --------------------------------------------------------------------------------------


overplot_panel, x_1, y_1, lx_rem_5, lx_rem_err_5, $
                                x0, y2, x1, y3, $
                                xmax, xmin, ymax, ymin, $
				yaxistitle=yaxistitle, /plot_xray_lum, /plabel

overplot_panel, xvar1, lx_max_20, lx_max_err_20, lx_rem_20, lx_rem_err_20, $
                                x1, y2, x2, y3, $ 
                                xmax, xmin, ymax, ymin, /plot_xray_lum

overplot_panel, xvar1, lx_max_40, lx_max_err_40, lx_rem_40, lx_rem_err_40, $
                                x2, y2, x3, y3, $ 
                                xmax, xmin, ymax, ymin, /plot_xray_lum

overplot_panel, xvar1, lx_max_100, lx_max_err_100, lx_rem_100, lx_rem_err_100, $
                                x3, y2, x4, y3, $ 
                                xmax, xmin, ymax, ymin, /plot_xray_lum




; --------------------------------------------------------------------------------------



;  Mass fractions
; ---------------------------
yaxistitle='!6mass fraction'
ymax = 0.44
;ymax = 0.25
ymin= 0.0


overplot_panel, xvar1, xraygm_5, xraygm_err_5, windgm_5, windgm_err_5, $
				x0, y1, x1, y2, $
				xmax, xmin, ymax, ymin, $
				yaxistitle=yaxistitle, /showtotal

overplot_panel, xvar1, xraygm_20, xraygm_err_20, windgm_20, windgm_err_20, $
				x1, y1, x2, y2, $
				xmax, xmin, ymax, ymin, /showtotal

overplot_panel, xvar1, xraygm_40, xraygm_err_40, windgm_40, windgm_err_40, $
				x2, y1, x3, y2, $
				xmax, xmin, ymax, ymin, /showtotal, /plabel

overplot_panel, xvar1, xraygm_100, xraygm_err_100, windgm_100, windgm_err_100, $
				x3, y1, x4, y2, $
				xmax, xmin, ymax, ymin, /showtotal






; --------------------------------------------------------------------------------------



;  Metallicity
; -------------------
yaxistitle='!6Metallicity (Z!D!9n!6!N)'
ymax = 3.0
ymin= 0.05


overplot_panel, xvar1, xraygasz_5, xraygasz_err_5, windz_5, windz_err_5, $
				x0, y0, x1, y1, $
				xmax, xmin, ymax, ymin, /ylog, $
				yaxistitle=yaxistitle, $
				xaxistitle=xaxistitle

overplot_panel, xvar1, xraygasz_20, xraygasz_err_20, windz_20, windz_err_20, $
				x1, y0, x2, y1, $
				xmax, xmin, ymax, ymin, /ylog, $
				xaxistitle=xaxistitle

overplot_panel, xvar1, xraygasz_40, xraygasz_err_40, windz_40, windz_err_40, $
				x2, y0, x3, y1, $
				xmax, xmin, ymax, ymin, /ylog, $
				xaxistitle=xaxistitle

overplot_panel, xvar1, xraygasz_100, xraygasz_err_100, windz_100, windz_err_100, $
				x3, y0, x4, y1, $
				xmax, xmin, ymax, ymin, /ylog, $
				xaxistitle=xaxistitle






; --------------------------------------------------------------------------------------


;---------
if not keyword_set(dontclose) then device, /close

end




; --------------------------------------------------------------------------------------


