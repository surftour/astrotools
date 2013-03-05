; =================================================
;
;   Generic plotting program to insert one
;   panel in a multiplot figure.
;
;
; ==================================================
;
pro plot_panel, x, y, y_err=y_err, $
		x0, y0, x1, y1, $
		xmax, xmin, ymax, ymin, $
		yaxistitle=yaxistitle, $
		xaxistitle=xaxistitle, $
		xlog=xlog, $
		ylog=ylog, $
		pointt=pointt


mode= 0
if keyword_set(xlog) then mode= mode + 1
if keyword_set(ylog) then mode= mode + 2
if keyword_set(xaxistitle) then mode= mode + 4
if keyword_set(yaxistitle) then mode= mode + 8


!p.position= [x0, y0, x1, y1]

case mode of

15: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	ytickformat='exp_label', $
	xtitle=xaxistitle, ytitle=yaxistitle, /xlog, /ylog
    end

14: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	ytickformat='exp_label', $
	xtitle=xaxistitle, ytitle=yaxistitle, /ylog
    end

13: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle, /xlog
    end

12: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle
    end

11: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	ytickformat='exp_label', $
	xtickformat='(a1)', ytitle=yaxistitle, /xlog, /ylog
    end

10: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	ytickformat='exp_label', $
	xtickformat='(a1)', ytitle=yaxistitle, /ylog
    end

9: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtickformat='(a1)', ytitle=yaxistitle, /xlog
    end

8: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtickformat='(a1)', ytitle=yaxistitle
    end

7: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	ytickformat='exp_label', $
	xtitle=xaxistitle, ytickformat='(a1)', /xlog, /ylog
    end

6: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	ytickformat='exp_label', $
	xtitle=xaxistitle, ytickformat='(a1)', /ylog
    end

5: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)', /xlog
    end

4: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytickformat='(a1)'
    end

3: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	ytickformat='exp_label', $
	xtickformat='(a1)', ytickformat='(a1)', /xlog, /ylog
    end

2: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	ytickformat='exp_label', $
	xtickformat='(a1)', ytickformat='(a1)', /ylog
    end

1: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtickformat='(a1)', ytickformat='(a1)', /xlog
    end

0: begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtickformat='(a1)', ytickformat='(a1)'
    end

endcase



; now, actually put the data down

select_thispoint, pointt, thispsym, thiscolor

oplot, x, y, psym=-thispsym, color=thiscolor, thick= 7.0, symsize= 1.5
if keyword_set(y_err) then begin
	oploterror, x, y, y_err, psym=-thispsym, errcolor=thiscolor, color=thiscolor, thick= 3.0, errthick= 3.0
endif



end



;================================================================================================
