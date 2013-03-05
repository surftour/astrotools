pro contour_makegeneralplot, x, y, xmax, xmin, ymax, ymin, sendto, filename=filename, $
			msg=msg, $
			q_eos=q_eos, $
			xaxistitle=xaxistitle, yaxistitle=yaxistitle, $
			phasediagram=phasediagram, sfrhothresh=sfrhothresh, $
			drawindexn=drawindexn, indexn=indexn, $
			colortable=colortable, plotcoolingtime=plotcoolingtime, $
			loadedsnap=loadedsnap

if not keyword_set(x) then begin
   print, "  "
   print, "PROBLEM: contour_makegeneralplot"
   print, "  "
   return
endif



;if not keyword_set(xmax) then xmax= 1.0
;if not keyword_set(xmin) then xmin= -1.0
;if not keyword_set(ymax) then ymax= 1.0
;if not keyword_set(ymin) then ymin= -1.0
        
;-------------------------------------
;  Setup Plot Stuff
;-------------------------------------

;initialize_plotinfo, 1

;setup_plot_stuff, sendto, filename=filename, colortable=3
;setup_plot_stuff, sendto, filename=filename, colortable=4



contour_makegeneralpic, x, y, xmax, xmin, ymax, ymin, NxNImage=NxNImage




;---------------------------------------
;  Print Tv image
;-----------------------------------------

x0= 0.15
y0= 0.15
x_size= 0.82
y_size= 0.82


!p.position=[x0, y0, x0+x_size,y0+y_size]

tv, NxNImage, x0, y0, xsize=x_size, ysize=y_size, /normal


;----------------------------------------
; Generate plot
;----------------------------------------

!p.position=[x0, y0, x0+x_size,y0+y_size]

plot, [1,0], [1,0], psym=3,xstyle=1,ystyle=1, $
	xrange=[xmin,xmax],$
	yrange=[ymin,ymax],$
	color=0, $
	xcharsize=1.5, ycharsize=1.5, $
	xthick=4.0, ythick=4.0, $
	charthick=3.0, $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/noerase, $
	/nodata


;---------------------------------------
;  Put contours on top of it
;---------------------------------------
docontour= 0
if docontour EQ 1 then begin
        levels = 6
        step = (Max(Pic) - Min(Pic)) / levels
        userLevels = IndGen(levels) * step + Min(Pic)

        ; load contours into variables (this doesn't print)
        contour, Pic, path_xy=xy_contours, path_info=info

        ; now print contours
        for i=0, (n_elements(info)-1) do begin
             cnt_index= [indgen(info(i).N),0]
             plots, xy_contours(*,info(i).offset+cnt_index), /normal, color= 0
        endfor
endif


;-----------------------------------------
; Print lines or words if desired
;-----------------------------------------


	if keyword_set(phasediagram) then begin
	   if keyword_set(sfrhothresh) then begin
		y=[ymin,ymax]
		x= sfrhothresh + 0.0*y
		x= alog10(x)
		oplot, x, y, linestyle=0, color=0
		ok=oplot_coldtemp_cutoff_log(1,rhocut=sfrhothresh)
	   endif else begin
		;ok=oplot_rho_cutoff_y_log(1)
		;ok=oplot_coldtemp_cutoff_log(1)

		;ok=oplot_rho_cutoff_y_log(1,rhocut=0.00854924)
		;ok=oplot_coldtemp_cutoff_log(1,rhocut=0.00854924)

		rhoc= 0.000854924 * 404.75368   ; in cm-3
		ok=oplot_rho_cutoff_y_log(1,rhocut=rhoc)
		ok=oplot_coldtemp_cutoff_log(1,rhocut=rhoc)
	   endelse
	endif

	if keyword_set(plotcoolingtime) then begin
        	;ok=oplot_coolingtime_log(1.0,sendto=sendto)
        	;ok=oplot_coolingtime_log(0.001,sendto=sendto)
		;ok=oplot_coolingtime_log(10.0,sendto=sendto)
		;ok=oplot_coolingtime_log(100.0,sendto=sendto)
        	;xyouts, 0.37, 0.87, 't!Dcool!N=1 Gyr', size= 1.2, color= 0, /normal, charthick= 3.0
		;xyouts, 0.47, 0.87, '1Gyr', size= 1.2, color= 0, /normal, charthick= 3.0
        	;xyouts, 0.62, 0.87, '1Myr', size= 1.2, color= 0, /normal, charthick= 3.0
	endif


	;indexn= 0
	if keyword_set(drawindexn) then begin
		len= 1.0
		slope= indexn/2.0
		theta= atan(slope)
		x0=0.0
		y0=4.8
		x1= x0+(len*cos(theta))*(xmax-xmin)/(ymax-ymin)
		y1= y0+(len*sin(theta))
		x=[x0,x1]
		y=[y0,y1]
		oplot, x, y, linestyle=0, color=0, thick=3.0
	endif


	if keyword_set(msg) then begin
		xyouts, 0.20, 0.88, msg, size=1.5, color=0, /normal, charthick=3.0
		xyouts, 0.20, 0.84, fload_timelbl(1,2), size=1.33, color=0, /normal, charthick= 3.0
	endif else begin
		xyouts, 0.20, 0.90, fload_fid(1), size=1.33, color=0, /normal
		xyouts, 0.20, 0.85, fload_timelbl(1,2), size=1.33, color=0, /normal, charthick= 3.0
	endelse


	if keyword_set(q_eos) then begin
        	qeoslbl = '!6q!DEQS!N='+strcompress(string(q_eos),/remove_all)
        	qeoslbl = strmid(qeoslbl,0,15)
        	xyouts, 0.20, 0.80, qeoslbl, size=1.33, color=0, /normal, charthick= 3.0
	endif



;if (sendto EQ 'ps') then device, /close



end


