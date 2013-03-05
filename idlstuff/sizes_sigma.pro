;==================================================================
;==================================================================
;
;
;
;==================================================================
;==================================================================


pro sigma_time, junk


if not keyword_set(junk) then begin
	print, " "
	print, " sigma_time, junk"
	print, " "
	print, " "
	return
endif

filename='sigma.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4


yaxistitle = "!7r!6(km/sec)"
xaxistitle = "!6Time (Gyr)"
xmax= 6.0
;xmax= 4.4
;xmax= 3.8
xmin = 0
;ymax = 15.0
ymax = 400.0
;ymax = 6.0
ymin = 0


;---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
	color= 0, $
	xstyle=1, $
	ystyle=1, $
	;ystyle=8, $     ; this will suppress the second y axis
	xcharsize=1.50, ycharsize=1.50, $
	xthick=4.0, ythick=4.0, $
	charthick=3.0, $
	;ytickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata

;---------------------------


do_one_sig, "data1/remergers/b3eb3e2", pt= 0
xyouts, 0.22, 0.94, "b3eb3e2", /normal, color=0, size=1.2

do_one_sig, "data1/remergers/b3ev2", pt= 200
xyouts, 0.22, 0.78, "b3ev2", /normal, color=200, size=1.2

do_one_sig, "data1/remergers/b3ev3", pt= 100
xyouts, 0.22, 0.82, "b3ev3", /normal, color=100, size=1.2

do_one_sig, "data1/remergers/b3ev4", pt= 150
xyouts, 0.22, 0.86, "b3ev4", /normal, color=150, size=1.2

do_one_sig, "data1/remergers/b3esph", pt= 50
xyouts, 0.22, 0.90, "b3esph", /normal, color=50, size=1.2


;---------------------------

; set the colortable to 1
; use xmax= 4.4 and ymax= 6.0

;xyouts, 0.52, 0.94, "ErrTolIntAccuracy", /normal, color= 0, size=1.2, charthick=2.0

;do_one_sig, "data/bs/b3e_i2", pt= 240
;xyouts, 0.55, 0.90, '1.0', /normal, color=240, size=1.2, charthick=2.0

;do_one_sig, "data/bs/b3e_i1", pt= 200
;xyouts, 0.55, 0.86, '0.1', /normal, color=200, size=1.2, charthick=2.0

;do_one_sig, "data/bs/b3e", pt= 160
;xyouts, 0.55, 0.82, '0.025', /normal, color=160, size=1.2, charthick=2.0

;do_one_sig, "data/bs/b3e_i4", pt= 120
;xyouts, 0.55, 0.78, '0.025 (diff Gad)', /normal, color=120, size=1.2, charthick=2.0

;do_one_sig, "data/bs/b3e_i3", pt= 80
;xyouts, 0.55, 0.74, '0.0025', /normal, color=80, size=1.2, charthick=2.0

;do_one_sig, "data/bs/b3e_i5", pt= 40
;xyouts, 0.55, 0.70, '0.00125', /normal, color=40, size=1.2, charthick=2.0

;do_one_sig, "data/bs/b3e_i6", pt= 0
;xyouts, 0.55, 0.66, '0.00025', /normal, color=0, size=1.2, charthick=2.0



;---------------------------

; turn off printing of error range

;xyouts, 0.52, 0.94, "Quiescent Galaxies", /normal, color= 0, size=1.2, charthick=2.0

;do_one_sig, "data1/remergers/iso_b3e", pt= 0
;xyouts, 0.75, 0.24, 'b3e,', /normal, color=0, size=1.2, charthick=2.0

;do_one_sig, "data1/remergers/iso_b3e_2", pt= 0
;xyouts, 0.82, 0.24, 'b3e_2', /normal, color=0, size=1.2, charthick=2.0

;do_one_sig, "data1/remergers/iso_sph", pt= 240
;xyouts, 0.75, 0.80, 'sph', /normal, color=240, size=1.2, charthick=2.0

;do_one_sig, "data1/remergers/iso_v2v2", pt= 150
;xyouts, 0.75, 0.32, 'v2v2', /normal, color=150, size=1.2, charthick=2.0

;do_one_sig, "data1/remergers/iso_v3x2", pt= 200
;xyouts, 0.75, 0.43, 'v3x2', /normal, color=200, size=1.2, charthick=2.0

;do_one_sig, "data1/remergers/iso_v4x3", pt= 100
;xyouts, 0.75, 0.58, 'v4x3', /normal, color=100, size=1.2, charthick=2.0

;do_one_sig, "data1/remergers/iso_v4x4", pt= 70
;xyouts, 0.65, 0.66, 'v4x4,', /normal, color=70, size=1.2, charthick=2.0

;do_one_sig, "data1/remergers/iso_v4x5", pt= 30
;xyouts, 0.83, 0.66, 'v4x5', /normal, color=30, size=1.2, charthick=2.0

;do_one_sig, "data1/remergers/iso_v4x6", pt= 50
;xyouts, 0.74, 0.66, 'v4x6,', /normal, color=50, size=1.2, charthick=2.0



;---------------------------


; done
; ------
device, /close


end






;==================================================================
;==================================================================



pro do_one_sig, frun, pt=pt

spawn, "/bin/ls "+frun+"/sigma_*.txt",result

for i=0, n_elements(result)-1 do begin

	sigfile= result[i]
	read_file_sigma, sigfile, snapext, time, sigxy, sigxz, sigyz, sigavg, sigerr



	time= time/0.7
	sigavg= sigavg
	sigerr= sigerr

	;select_thispoint, pt, thispsym, thiscolor
	thispsym= 3
	thiscolor= pt

	oplot, time, sigavg, thick=4.0, psym=-thispsym, color=thiscolor
	;oplot, time, sigavg+sigerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
	;oplot, time, sigavg-sigerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
endfor


end



;==================================================================
;==================================================================



