
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Seperation of MW-M31 vs. time
;     -------------------------------------------
;  plot the distance between the two interacting galaxies as a function
;  of time
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------











;----------------------------------
;
;
;
;----------------------------------
pro showcen, junk


if not keyword_set(junk) then begin
	print, "  "
	print, "  show_centers_time, junk"
	print, "  "
	print, "  needs: .run time_centers"
	print, "  "
	return
endif


filename='lgcenters.eps'
;filename='lgcenters_test.eps'

timemin= 0.0
timemax= 12.5
;timemax= 15.0
;timemax= 25.0
;timemax= 35.0

ymax= 1.4


;  mark today's position
;
t_sep= 0.780 ; in Mpc
t_sep_err= 0.035




; ----------------------------------
;   Try this plot thingy
; ----------------------------------

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4

!p.position= [0.19, 0.15, 0.97, 0.98]


; -------------------------------
; plot axes
;
plot, [1.0], [1.0], psym=-3, xrange=[timemin,timemax], yrange=[0,ymax], xstyle=1, ystyle=1, $
        color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=2.0, /nodata, $
        xtitle="!6Time (Gyr)", ytitle="!6MW/M31 Separation (Mpc)"


; -------------------------------
; hatched region
;
do_hatched= 0
;do_hatched= 1
if do_hatched eq 1 then begin
	y= [t_sep-t_sep_err, t_sep-t_sep_err, t_sep+t_sep_err, t_sep+t_sep_err, t_sep-t_sep_err]
	x= [timemin,timemax,timemax,timemin,timemin]

	polyfill, x, y, /data, color= 220, /fill
	xyouts, 0.57, 0.61, 'separation today', charthick=2.0, size=1.2, color= 0, /normal

	; replot axes
	plot, [1.0], [1.0], psym=-3, xrange=[timemin,timemax], yrange=[0,ymax], xstyle=1, ystyle=1, $
	        color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=2.0, /nodata, $
	        xtitle="!6Time (Gyr)", ytitle="!6MW/M31 Separation (Mpc)", /noerase
endif


; -------------------------------
; go through files
; 

x0= 0.65
y0= 0.90

;process_onesim_centers, "/raid4/tcox/localgroup/v4", 2, x0, y0
;process_onesim_centers, "/raid4/tcox/localgroup/v5", 1, x0, y0-0.03, /nolabel
;process_onesim_centers, "/raid4/tcox/localgroup/v6", 3, x0, y0-0.06
;process_onesim_centers, "/raid4/tcox/localgroup/v7", 3, x0, y0-0.09
;process_onesim_centers, "/raid4/tcox/localgroup/v8", 3, x0, y0-0.12
;process_onesim_centers, "/raid4/tcox/localgroup/bhires", 4, x0, y0-0.09
;process_onesim_centers, "data/v4", 2, x0, y0
;process_onesim_centers, "data/v5", 1, x0, y0-0.03
;process_onesim_centers, "/raid4/tcox/localgroup/v5_noigm", 1, x0, y0-0.03, /nolabel
;process_onesim_centers, "data/v6", 3, x0, y0-0.06
;process_onesim_centers, "data/bhires", 4, x0, y0-0.09

;process_onesim_centers, "/raid4/tcox/localgroup/v1", 11, x0, y0-0.03, /nolabel
;process_onesim_centers, "/raid4/tcox/localgroup/v2", 11, x0, y0-0.03, /nolabel
;process_onesim_centers, "/raid4/tcox/localgroup/v3", 11, x0, y0-0.03, /nolabel

process_onesim_centers, "/raid4/tcox/localgroup/v4", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v5", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v6", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v7", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v8", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v9", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v10", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v11", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v12", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v13", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v14", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v15", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v16", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v17", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v18", 11, x0, y0-0.03, /nolabel
process_onesim_centers, "/raid4/tcox/localgroup/v19", 11, x0, y0-0.03, /nolabel


;process_onesim_centers, "/raid4/tcox/localgroup/v5", 1, x0, y0-0.03, /nolabel
;process_onesim_centers, "/raid4/tcox/localgroup/v5_noigm", 3, x0, y0-0.07, /nolabel
;xyouts, 0.64, 0.85, 'fiducial', charthick=2.0, size=1.3, color= 150, /normal
;xyouts, 0.64, 0.80, 'no IGM', charthick=2.0, size=1.3, color= 50, /normal




; -------------------------------
;  Draw Circle
;
;draw_re_circle= 1
draw_re_circle= 0
if draw_re_circle eq 1 then begin
	;timetoday= 5.0
	timetoday= 4.7
	cfac= 0.05
   
	phi = dindgen(101)*2D*!dpi/100
	; data coord
	re_x = (timemax-timemin) * cfac * cos(phi)
	re_y = (ymax-0.0) * cfac * sin(phi)
   
	re_x= re_x + timetoday
	re_y= re_y + t_sep

	polyfill, re_x, re_y, /data, color= 50, /line_fill, linestyle=0, $
                                thick=1.0, orientation=45.0

	oplot, re_x, re_y, color= 50, thick=4.0, linestyle=0, psym=-3
	xyouts, 0.54, 0.66, 'today', charthick=2.0, size=1.3, color= 50, /normal

endif


; -------------------------------
;  dotted line
;
x=[timemin,timemax]
y=[t_sep, t_sep]
oplot, x, y, psym=-3, linestyle=1, thick=3.0, color= 0
;xyouts, 0.57, 0.60, 'separation today', charthick=2.0, size=1.2, color= 0, /normal


; -------------------------------

device, /close


end



;=================================================================================




;----------------------------------
;
;
;
;----------------------------------
pro showcen_norm, junk


if not keyword_set(junk) then begin
	print, "  "
	print, "  showcen_norm, junk"
	print, "  "
	print, "  needs: .run time_centers"
	print, "  "
	return
endif


filename='lgcenters_norm.eps'

timemin= -6.5
timemax= 8.5
;timemax= 15.0
;timemax= 25.0
;timemax= 35.0

ymax= 1.4


;  mark today's position
;
t_sep= 0.780 ; in Mpc
t_sep_err= 0.035




; ----------------------------------
;   Try this plot thingy
; ----------------------------------

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4

!p.position= [0.19, 0.15, 0.97, 0.98]


; -------------------------------
; plot axes
;
plot, [1.0], [1.0], psym=-3, xrange=[timemin,timemax], yrange=[0,ymax], xstyle=1, ystyle=1, $
        color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=2.0, /nodata, $
        xtitle="!6Time Relative To Now (Gyr)", ytitle="!6MW/M31 Separation (Mpc)"


; -------------------------------
; go through files
; 

x0= 0.65
y0= 0.90

;process_onesim_centers, "/raid4/tcox/localgroup/v1", 11, x0, y0-0.03, /nolabel, /timenorm
;process_onesim_centers, "/raid4/tcox/localgroup/v2", 11, x0, y0-0.03, /nolabel, /timenorm
;process_onesim_centers, "/raid4/tcox/localgroup/v3", 11, x0, y0-0.03, /nolabel, /timenorm

process_onesim_centers, "/raid4/tcox/localgroup/v4", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v5", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v6", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v7", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v8", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v9", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v10", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v11", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v12", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v13", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v14", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v15", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v16", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v17", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v18", 11, x0, y0-0.03, /nolabel, /timenorm
process_onesim_centers, "/raid4/tcox/localgroup/v19", 11, x0, y0-0.03, /nolabel, /timenorm

; list of first passage times (relative to now)
; t_firstpassage= [2.38, 2.26, 2.31, 2.60, 2.67, 2.41, 2.45, 2.41, 2.58, 2.68, 3.13, 3.01, 3.40, 3.41, 3.80, 3.78]
;
; list of final merger times (relative to now)
; f_mergertimes= [6.38, 4.98, 5.88, 5.02, 5.38, 5.70, 6.02, 5.84, 5.29, 5.11, 4.98, 5.01, 5.11, 5.12, 5.37, 5.50]
;



; -------------------------------
;  dotted line
;
x=[timemin,timemax]
y=[t_sep, t_sep]
oplot, x, y, psym=-3, linestyle=1, thick=3.0, color= 0
;xyouts, 0.57, 0.60, 'separation today', charthick=2.0, size=1.2, color= 0, /normal


; -------------------------------

device, /close


end



;=================================================================================
;=================================================================================





; do the dirty work
; --------------------
pro process_onesim_centers, frun, pointselection, x0, y0, $
				nolabel=nolabel, timenorm=timenorm

	print, "processing: ", frun

	; read centers
	; -------------
	read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=frun+"/centers.txt"

	cp_time= transpose(cp_time)
	cp_cen1= transpose(cp_cen1)
	cp_cen2= transpose(cp_cen2)

	x = cp_cen1[*,0] - cp_cen2[*,0]
	y = cp_cen1[*,1] - cp_cen2[*,1]
	z = cp_cen1[*,2] - cp_cen2[*,2]

	rdiff = sqrt(x*x + y*y + z*z)

	;rdiff= rdiff/0.7
	rdiff= rdiff/ 1000.0 / 0.7   ; transfor from kpc/h to Mpc
	cp_time= cp_time/0.7


;idx= where(rdiff lt 0.820)
;print, "rdiff= ", rdiff(idx(0)), cp_time(idx(0))
;print, "rdiff= ", rdiff(idx(1)), cp_time(idx(1))
;print, "rdiff= ", rdiff(idx(2)), cp_time(idx(2))

	now_separation= 0.78

	if keyword_set(timenorm) then begin
		idx= where(rdiff le now_separation) 
		print, "now separation= ", now_separation

		i1= idx(0)
		i0= i1-1

		dt= cp_time[i1]-cp_time[i0]

		m= (rdiff[i1]-rdiff[i0])/dt
		b= rdiff[i1] - m*cp_time[i1]

		t_now= (0.7800 - b)/m 
		print, "t_now= ", t_now, idx(0)

		cp_time= cp_time - t_now

	endif

	for i=0,n_elements(rdiff)-2 do begin
		if rdiff(i+1) gt rdiff(i) then begin
		   print, "t_firstp= ", cp_time(i+1), i+1
		   break
		endif
	endfor

	idx=where(rdiff le 0)
	if idx(0) ne -1 then print, "t_merge= ", cp_time(idx(0)), idx(0)

	; get line color and point type
	select_thispoint, pointselection, thispsym, thiscolor

	oplot, cp_time, rdiff, psym=-thispsym, color= thiscolor, thick=3.0

	if not keyword_set(nolabel) then begin
		xyouts, x0, y0, fload_getid(frun), /normal, size= 1.1, color= thiscolor, charthick= 3.0
	endif

end






;======================================================================





pro veli,  junk, filename=filename


if not keyword_set(junk) then begin
        print, "  "
        print, "  velocity_info, junk"
        print, "  "
	print, "  needs: .run time_centers"
	print, "  "
        return
endif


; ----------------------------------
; ----------------------------------

timemin= 0.0
timemax= 12.0

ymax= 550.0

if not keyword_set(filename) then filename='lgvelocity_info.eps'

; ----------------------------------
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4

!p.position= [0.18, 0.13, 0.97, 0.99]


; ----------------------------------
;  plot axes
;
plot, [1.0], [1.0], psym=-3, xrange=[timemin,timemax], yrange=[0,ymax], xstyle=1, ystyle=1, $
        color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=2.0, /nodata, $
        xtitle="!6Time (Gyr)", ytitle="!6Velocity (km s!E-1!N)"


; -------------------------------
; hatched region
;
;do_hatched= 0
do_hatched= 1
if do_hatched eq 1 then begin
	velt= 120.0
	velt_err= 13.0

	todayt= 4.7
	tterr= 0.35

        y= [velt-velt_err, velt-velt_err, velt+velt_err, velt+velt_err, velt-velt_err]
        x= [todayt-tterr,todayt+tterr,todayt+tterr,todayt-tterr,todayt-tterr]

        polyfill, x, y, /data, color= 50, /line_fill, linestyle=0, $
                                thick=1.0, orientation=45.0

        xyouts, 0.45, 0.24, 'today', charthick=2.0, size=1.2, color= 50, /normal
	oplot, x, y, color= 50, thick=4.0, linestyle=0, psym=-3

        ; replot axes
	;plot, [1.0], [1.0], psym=-3, xrange=[timemin,timemax], yrange=[0,ymax], xstyle=1, ystyle=1, $
	;        color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=2.0, /nodata, $
	;        xtitle="!6Time (Gyr)", ytitle="!6Velocity (km sec!E-1!N)"
endif


; ----------------------------------
; go through files

x0= 0.65
y0= 0.90

if not keyword_set(junk) then frun= "/raid4/tcox/localgroup/v8"
process_onesim_centervelocity, junk, 3, x0, y0-0.06


;process_onesim_centervelocity, "/raid4/tcox/localgroup/v4", 2, x0, y0
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v5", 1, x0, y0-0.03
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v6", 3, x0, y0-0.06
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v7", 3, x0, y0-0.06
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v8", 3, x0, y0-0.06
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v9", 3, x0, y0-0.06
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v10", 3, x0, y0-0.06
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v11", 3, x0, y0-0.06
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v12", 3, x0, y0-0.06
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v13", 3, x0, y0-0.06
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v14", 3, x0, y0-0.06
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v15", 3, x0, y0-0.06
;process_onesim_centervelocity, "/raid4/tcox/localgroup/bhires", 4, x0, y0-0.09
;process_onesim_centervelocity, "data/v4", 2, x0, y0
;process_onesim_centervelocity, "data/v5", 1, x0, y0-0.03
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v5", 1, x0, y0-0.03, /nolabel
;process_onesim_centervelocity, "data/v6", 3, x0, y0-0.06
;process_onesim_centervelocity, "data/bhires", 4, x0, y0-0.09

; ----------------------------------

;process_onesim_centervelocity, "/raid4/tcox/localgroup/v5", 1, x0, y0-0.03, /nolabel
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v5_noigm", 2, x0, y0-0.08, /nolabel
oplot, [0.5, 1.4], [520.0, 520.0], color= 0, psym= -3, linestyle= 0, thick= 5.0
xyouts, 0.30, 0.93, 'total relative', /normal, size= 1.2, color= 0, charthick= 3.0
xyouts, 0.33, 0.89, 'velocity', /normal, size= 1.2, color= 0, charthick= 3.0
oplot, [0.5, 1.4], [460.0, 460.0], color= 150, psym= -3, linestyle= 1, thick= 5.0
xyouts, 0.30, 0.84, 'radial', /normal, size= 1.2, color= 150, charthick= 3.0
oplot, [0.5, 1.4], [430.0, 430.0], color= 0, psym= -3, linestyle= 2, thick= 5.0
xyouts, 0.30, 0.79, 'tangential', /normal, size= 1.2, color= 50, charthick= 3.0


; ----------------------------------
; today's estimated velocity?? 
;
;x= [5.4,5.4]
;y=[0,300.0]
;oplot, x, y, psym=-3, color= 0, thick=3.0, linestyle= 1


; -------------------------------

device, /close



end


;======================================================================





pro showvel,  junk, filename=filename


if not keyword_set(junk) then begin
        print, "  "
        print, "  showvel, junk"
        print, "  "
	print, "  needs: .run time_centers"
	print, "  "
        return
endif


; ----------------------------------
; ----------------------------------

timemin= 0.0
timemax= 12.0

ymax= 550.0

if not keyword_set(filename) then filename='lgvelocity_info.eps'

; ----------------------------------
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4

!p.position= [0.18, 0.13, 0.97, 0.99]


; ----------------------------------
;  plot axes
;
plot, [1.0], [1.0], psym=-3, xrange=[timemin,timemax], yrange=[0,ymax], xstyle=1, ystyle=1, $
        color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=2.0, /nodata, $
        xtitle="!6Time (Gyr)", ytitle="!6Velocity (km s!E-1!N)"


; -------------------------------
; hatched region
;
do_hatched= 0
;do_hatched= 1
if do_hatched eq 1 then begin
	velt= 120.0
	velt_err= 13.0

	todayt= 4.7
	tterr= 0.35

        y= [velt-velt_err, velt-velt_err, velt+velt_err, velt+velt_err, velt-velt_err]
        x= [todayt-tterr,todayt+tterr,todayt+tterr,todayt-tterr,todayt-tterr]

        polyfill, x, y, /data, color= 50, /line_fill, linestyle=0, $
                                thick=1.0, orientation=45.0

        xyouts, 0.45, 0.24, 'today', charthick=2.0, size=1.2, color= 50, /normal
	oplot, x, y, color= 50, thick=4.0, linestyle=0, psym=-3

        ; replot axes
	;plot, [1.0], [1.0], psym=-3, xrange=[timemin,timemax], yrange=[0,ymax], xstyle=1, ystyle=1, $
	;        color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=2.0, /nodata, $
	;        xtitle="!6Time (Gyr)", ytitle="!6Velocity (km sec!E-1!N)"
endif


; ----------------------------------
; go through files

x0= 0.65
y0= 0.90

if not keyword_set(junk) then frun= "/raid4/tcox/localgroup/v8"

;process_onesim_centervelocity, "/raid4/tcox/localgroup/v1", 11, x0, y0, /nolabel
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v2", 11, x0, y0, /nolabel
;process_onesim_centervelocity, "/raid4/tcox/localgroup/v3", 11, x0, y0, /nolabel

process_onesim_centervelocity, "/raid4/tcox/localgroup/v4", 11, x0, y0, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v5", 11, x0, y0-0.03, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v6", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v7", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v8", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v9", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v10", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v11", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v12", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v13", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v14", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v15", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v16", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v17", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v18", 11, x0, y0-0.06, /nolabel
process_onesim_centervelocity, "/raid4/tcox/localgroup/v19", 11, x0, y0-0.06, /nolabel


;process_onesim_centervelocity, "/raid4/tcox/localgroup/bhires", 11, x0, y0-0.09, /nolabel


;xyouts, 0.26, 0.82, 'radial', /normal, size= 1.8, color= 0, charthick= 3.0
xyouts, 0.26, 0.82, 'transverse', /normal, size= 1.8, color= 0, charthick= 3.0
xyouts, 0.26, 0.75, 'velocity', /normal, size= 1.8, color= 0, charthick= 3.0


; -------------------------------
;  dotted line
;

;x=[timemin,timemax]
;y=[120.0, 120.0]
;oplot, x, y, psym=-3, linestyle=1, thick=3.0, color= 150
;xyouts, 0.57, 0.60, 'separation today', charthick=2.0, size=1.2, color= 0, /normal


; ----------------------
device, /close



end







pro showvel_both,  junk, filename=filename


if not keyword_set(junk) then begin
        print, "  "
        print, "  showvel, junk"
        print, "  "
	print, "  needs: .run time_centers"
	print, "  "
        return
endif


; ----------------------------------
; ----------------------------------

timemin= 0.0
timemax= 12.0

ymax= 550.0


;frun1= "/raid4/tcox/localgroup/v5"
;frun2= "/raid4/tcox/localgroup/v14"
;frun3= "/raid4/tcox/localgroup/v16"

;frun1= "/raid4/tcox/localgroup/v6"
;frun2= "/raid4/tcox/localgroup/v4"
;frun3= "/raid4/tcox/localgroup/v12"

frun1= "/raid4/tcox/localgroup/v7"
frun2= "/raid4/tcox/localgroup/v10"
frun3= "/raid4/tcox/localgroup/v18"

;frun1= "/raid4/tcox/localgroup/v8"
;frun2= "/raid4/tcox/localgroup/v17"
;frun3= "/raid4/tcox/localgroup/v19"


; ----------------------------------


filename='lgv_tan.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4

!p.position= [0.18, 0.13, 0.97, 0.99]

plot, [1.0], [1.0], psym=-3, xrange=[timemin,timemax], yrange=[0,ymax], xstyle=1, ystyle=1, $
        color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=2.0, /nodata, $
        xtitle="!6Time (Gyr)", ytitle="!6Velocity (km s!E-1!N)"


x0= 0.65
y0= 0.90


process_onesim_centervelocity, frun1, 5, x0, y0, /do_tan
process_onesim_centervelocity, frun2, 3, x0, y0-0.03, /do_tan
process_onesim_centervelocity, frun3, 2, x0, y0-0.06, /do_tan


xyouts, 0.26, 0.82, 'transverse', /normal, size= 1.8, color= 0, charthick= 3.0
xyouts, 0.26, 0.75, 'velocity', /normal, size= 1.8, color= 0, charthick= 3.0


device, /close


; -------------------------------
; -------------------------------



filename='lgv_rad.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4

!p.position= [0.18, 0.13, 0.97, 0.99]

plot, [1.0], [1.0], psym=-3, xrange=[timemin,timemax], yrange=[0,ymax], xstyle=1, ystyle=1, $
        color= 0, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=2.0, /nodata, $
        xtitle="!6Time (Gyr)", ytitle="!6Velocity (km s!E-1!N)"


x0= 0.65
y0= 0.90


process_onesim_centervelocity, frun1, 5, x0, y0, /do_rad
process_onesim_centervelocity, frun2, 3, x0, y0-0.03, /do_rad
process_onesim_centervelocity, frun3, 2, x0, y0-0.06, /do_rad


xyouts, 0.26, 0.82, 'radial', /normal, size= 1.8, color= 0, charthick= 3.0
xyouts, 0.26, 0.75, 'velocity', /normal, size= 1.8, color= 0, charthick= 3.0


device, /close


; -------------------------------
; -------------------------------









end









;======================================================================================








; do the dirty work
; --------------------
pro process_onesim_centervelocity, frun, pointselection, x0, y0, $
				nolabel=nolabel, do_tan=do_tan, do_rad=do_rad

	if (not keyword_set(do_tan)) and (not keyword_set(do_rad)) then do_tan= 1

	print, "processing: ", frun

        ; read centers
        ; -------------
        read_centerpositions, center_time, center_1, center_2, filename=frun+"/centers.txt"

        center_time= transpose(center_time)
        center_1= transpose(center_1)
        center_2= transpose(center_2)

        center_time= center_time/0.7
        center_1= center_1/0.7
        center_2= center_2/0.7

        x = center_1[*,0] - center_2[*,0]
        y = center_1[*,1] - center_2[*,1]
        z = center_1[*,2] - center_2[*,2]

        rdiff = sqrt(x*x + y*y + z*z)


        ; read com velocity
        ; ------------------
        ;read_centerpositions, comvel_time, comvel_1, comvel_2, filename=frun+"/comvel.txt"
        read_centerpositions, comvel_time, comvel_1, comvel_2, filename=frun+"/cenvel.txt"

        comvel_time= transpose(comvel_time)
        comvel_1= transpose(comvel_1)
        comvel_2= transpose(comvel_2)


        vx= comvel_1[*,0] - comvel_2[*,0] 
        vy= comvel_1[*,1] - comvel_2[*,1]
        vz= comvel_1[*,2] - comvel_2[*,2]

;print, "multiplying by factor 0.88"
;	vx= vx*0.88
;	vy= vy*0.88
;	vz= vz*0.88

        velrel= sqrt(vx*vx + vy*vy + vz*vz)

	idx=where(rdiff gt 0.0)
        vrad_x= vx(idx)*x(idx)/rdiff(idx) 
        vrad_y= vy(idx)*y(idx)/rdiff(idx)
        vrad_z= vz(idx)*z(idx)/rdiff(idx)

        vrad= sqrt(vrad_x*vrad_x + vrad_y*vrad_y + vrad_z*vrad_z)
        vtan = sqrt(velrel*velrel - vrad*vrad)

        ;for i=0,n_elements(center_time)-1 do begin
        ;        print, center_time[i], rdiff[i], velrel[i], vrad[i], vtan[i]
        ;endfor
	fiddle= 0
	;fiddle= 1
	if fiddle eq 1 then begin
	   vrad[0]= 0      &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[1]= 1.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[2]= 2.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[3]= 3.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[4]= 4.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[5]= 6.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[6]= 8.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[7]= 10.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[8]= 12.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[9]= 14.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[10]= 17.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[11]= 20.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[12]= 24.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[13]= 28.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[14]= 32.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[15]= 36.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[16]= 40.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[17]= 44.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[18]= 48.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[19]= 52.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[20]= 57.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[21]= 62.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[22]= 67.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[23]= 72.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[24]= 77.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[25]= 83.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[26]= 90.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[27]= 98.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[28]= 106.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[29]= 113.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[30]= 120.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])
	   vrad[31]= 127.0    &  velrel[0]= sqrt(vrad[0]*vrad[0] + vtan[0]*vtan[0])

	   for i=0, 31 do velrel[i]= sqrt(vrad[i]*vrad[i] + vtan[i]*vtan[i])
        endif

print, comvel_time(0), rdiff(0), vrad(0), vtan(0)

;idx= where(rdiff le 820.0)
;print, "TODAY"
;print, idx
;print, "Time= ", comvel_time(idx[0:5])
;print, "Separation= ", rdiff(idx[0:5])
;print, "Radial Velocity= ", vrad(idx[0:5])
;print, "Tangential Velocity= ", vtan(idx[0:5])

;nct= n_elements(comvel_time)-1
;print, "Time= ", comvel_time(nct-5:nct)
;print, "Separation= ", rdiff(nct-5:nct)
;print, "Radial Velocity= ", vrad(nct-5:nct)
;print, "Tangential Velocity= ", vtan(nct-5:nct)
;stop

	;
        ; get line color and point type
        select_thispoint, pointselection, thispsym, thiscolor


	; ====  1  =====
	; version one: plot one for several runs
	;
	;oplot, center_time, velrel, psym=-thispsym, color= thiscolor, thick=5.0
	if keyword_set(do_rad) then oplot, center_time, vrad, psym=-thispsym, color= thiscolor, thick=5.0, linestyle= 0
	if keyword_set(do_tan) then oplot, center_time, vtan, psym=-thispsym, color= thiscolor, thick=5.0, linestyle= 0


	; ====  2  =====
	; version two: plot three components for one run
	;
	;oplot, center_time, velrel, psym=-3, color= 0, thick=5.0
	;oplot, center_time, vrad, psym=-3, color= 150, thick=5.0, linestyle= 1
	;oplot, center_time, vtan, psym=-3, color= 50, thick=5.0, linestyle= 2
	;xyouts, 0.22, 0.93, 'total relative', /normal, size= 1.2, color= 0, charthick= 3.0
	;xyouts, 0.25, 0.89, 'velocity', /normal, size= 1.2, color= 0, charthick= 3.0
	;xyouts, 0.22, 0.84, 'radial', /normal, size= 1.2, color= 150, charthick= 3.0
	;xyouts, 0.22, 0.79, 'tangential', /normal, size= 1.2, color= 50, charthick= 3.0


	; always do this
	if not keyword_set(nolabel) then begin
		xyouts, x0, y0, fload_getid(frun), /normal, size= 1.1, color= thiscolor, charthick= 3.0
	endif

end





;===============================================================





; grab veli
; --------------------
pro grab_veli, frun, t_now, v_tot, v_rad, v_tan


	print, "processing: ", frun

        ; read centers
        ; -------------
        read_centerpositions, center_time, center_1, center_2, filename=frun+"/centers.txt"

        center_time= transpose(center_time)
        center_1= transpose(center_1)
        center_2= transpose(center_2)

        center_time= center_time/0.7
        center_1= center_1/0.7
        center_2= center_2/0.7

        x = center_1[*,0] - center_2[*,0]
        y = center_1[*,1] - center_2[*,1]
        z = center_1[*,2] - center_2[*,2]

        rdiff = sqrt(x*x + y*y + z*z)


        ; read com velocity
        ; ------------------
        ;read_centerpositions, comvel_time, comvel_1, comvel_2, filename=frun+"/comvel.txt"
        read_centerpositions, comvel_time, comvel_1, comvel_2, filename=frun+"/cenvel.txt"

        comvel_time= transpose(comvel_time)
        comvel_1= transpose(comvel_1)
        comvel_2= transpose(comvel_2)


        vx= comvel_1[*,0] - comvel_2[*,0] 
        vy= comvel_1[*,1] - comvel_2[*,1]
        vz= comvel_1[*,2] - comvel_2[*,2]

;print, "multiplying by factor 0.88"
;	vx= vx*0.88
;	vy= vy*0.88
;	vz= vz*0.88

        velrel= sqrt(vx*vx + vy*vy + vz*vz)

	idx=where(rdiff gt 0.0)
        vrad_x= vx(idx)*x(idx)/rdiff(idx) 
        vrad_y= vy(idx)*y(idx)/rdiff(idx)
        vrad_z= vz(idx)*z(idx)/rdiff(idx)

        vrad= sqrt(vrad_x*vrad_x + vrad_y*vrad_y + vrad_z*vrad_z)
        vtan = sqrt(velrel*velrel - vrad*vrad)


idx= where(rdiff le 780.0)
print, "now separation= ", 780.0

i1= idx(0)
i0= i1-1

dt= comvel_time[i1]-comvel_time[i0]

m= (rdiff[i1]-rdiff[i0])/dt
b= rdiff[i1] - m*comvel_time[i1]

t_now= (780.0 - b)/m 
print, "t_now= ", t_now

m= (vrad[i1]-vrad[i0])/dt
b = vrad[i1] - m*comvel_time[i1]
v_rad= m * t_now + b
print, "v_rad= ", v_rad

m= (vtan[i1]-vtan[i0])/dt
b = vtan[i1] - m*comvel_time[i1]
v_tan= m * t_now + b
print, "v_tan= ", v_tan


v_tot= sqrt(v_rad*v_rad + v_tan*v_tan)
print, "v_tot= ", v_tot

end




