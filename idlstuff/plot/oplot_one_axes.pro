;----------------------------------
;
;
;
;----------------------------------
pro oplot_one_axes, xx0, yy0, xx1, yy1, thetimelbl, xlen, $
			extramsg=extramsg, add_j_vector=add_j_vector, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			xtitle=xtitle, ytitle=ytitle

        !p.position=[xx0,yy0,xx1,yy1]
	!p.ticklen=0.03

	xmin= -xlen
	xmax=  xlen
	ymin= -xlen
	ymax=  xlen


	if keyword_set(xtitle) and keyword_set(ytitle) then goto, bothtitles
	if not keyword_set(xtitle) and keyword_set(ytitle) then goto, ytitleonly
	if keyword_set(xtitle) and not keyword_set(ytitle) then goto, xtitleonly

	goto, notitle


xtitleonly:
        	plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
                      xthick=3.0, ythick=3.0, ytickformat='(a1)', xtitle=xtitle
		goto, moveon


ytitleonly:
        	plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
                      xthick=3.0, ythick=3.0, xtickformat='(a1)', ytitle=ytitle
		goto, moveon


bothtitles:
        	plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
                      xthick=3.0, ythick=3.0, xtitle=xtitle, ytitle=ytitle
		goto, moveon


notitle:
        	plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xstyle=1, ystyle=1, /normal, /nodata, $
                      xthick=3.0, ythick=3.0, xtickformat='(a1)', ytickformat='(a1)'
		goto, moveon

moveon:

        ; white box
        ;bx0=xmin+3.0
        ;bx1=xmin+21.5
        ;by0=ymax-12.0
        ;by1=ymax-4.0
        ;polyfill, [bx0,bx0,bx1,bx1], [by0,by1,by1,by0], color= 1
        ;oplot, [bx0,bx0,bx1,bx1,bx0], [by0,by1,by1,by0,by0] , color= 0, thick=3.0

	if keyword_set(thetimelbl) and thetimelbl ne '-1' then begin
        	xyouts, xx0+0.03, yy1-0.065, '!6'+thetimelbl, /normal, size= 1.2, charthick=3.0, color= 0
        	;xyouts, xx0+0.02, yy1-0.03, '!6'+thetimelbl, /normal, size= 2.0, charthick=3.0, color= 0
	endif

	if keyword_set(extramsg) then begin
		xyouts, xx0+0.02, yy1-0.10, extramsg, /normal, size= 1.2, charthick=3.0, color= 0
	endif

	xlenlbl= strcompress(string(2.0*xlen),/remove_all)
	if (2.0*xlen) ge 1.0 then digs= 1
	if (2.0*xlen) ge 10.0 then digs= 2
	if (2.0*xlen) ge 100.0 then digs= 3
        xlenlbl = strmid(xlenlbl,0,digs)        ; T=0.x (4+digits after decimal)
	xlenlbl= '!94!6'+xlenlbl+'!96!6'
	;xyouts, xx1-0.10, yy0+0.015, xlenlbl, /normal, size= 1.2, charthick=3.0, color= 0

end




