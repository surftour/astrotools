
; =================================================
;
;  Wind mass versus time
;
;
;
;
; ==================================================


pro plot_wind_comparison, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_wind_comparison, junk"
	print, " "
	print, " "
	return
endif

filename='wind.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;yaxistitle = "Mass (10!E10!N M!D!9n!E!N)"
yaxistitle = "!6Log Wind Mass (M!D!9n!6!N)"
;yaxistitle = "Log Mass (M)"
;xaxistitle = "Time after merger (Gyr)"
xaxistitle = "!6Time (Gyr)"
;xmax = max(time)
xmax = 4.25
;xmax = 3.0
;xmax = 2.0
xmin = 0
;ymax = 8.0
;ymin = 0.0
;ymin = 5.0e-4    ; 10^7 solar masses
;ymax = 12.0
ymax = 11.5
;ymax = 10.0
ymin = 8.0
;ymin = 6.0


;---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




x0= 0.65
y0= 0.35

; Galaxy 1
;----------
;frun="pool/vc3vc3"
;frun="pool/vc3bvc3b"
;frun="/raid4/tcox/vc3avc3a_3"
;frun="/raid4/tcox/vc3ba_2_3"
;frun="/raid4/tcox/vc1vc1"
;frun="/raid4/tcox/As/A1"
;frun="/raid4/tcox/vc3vc3e_2"
;frun="/raid4/tcox/vc3vc3e"
;frun="/raid4/tcox/vc3vc3h"

;process_one_wind, "/raid4/tcox/vc3vc3h", lcolor=50, msg='std', x0= 0.65, y0= 0.35



; Galaxy 2
;----------
;frun="pool/vc3vc3_wBH"
;frun="pool/vc3bvc3b_wBH"
;frun="/raid4/tcox/vc3avc3a_2"
;frun="/raid4/tcox/vc3ba_2_2"
;frun="/raid4/tcox/vc2vc2"
;frun="/raid4/tcox/As/A2"
;frun="/raid4/tcox/vc3vc3e_no"
;frun="/raid4/tcox/vc3vc3e"

;process_one_wind, "/raid4/tcox/vc3vc3e", lcolor=150, msg='tilted', x0= 0.65, y0= 0.31




; Galaxy 3
;----------
;frun="/raid4/tcox/vc3avc3a_1"
;frun="/raid4/tcox/vc3ba_2_1"
;frun="/raid4/tcox/vc3bvc3b"
;frun="/raid4/tcox/vc3vc3"
;frun="/raid4/tcox/As/A3"


; Galaxy 4
;----------
;frun="/raid4/tcox/vc3avc3a_4"
;frun="/raid4/tcox/vc3ba_2_4"
;frun="/raid4/tcox/vc4vc4a"
;frun="/raid4/tcox/As/A4"


; Galaxy 5
;----------
;frun="/raid4/tcox/vc3avc3a_6"      ; vc3avc3a_5 isn't there, don't know why
;frun="/raid4/tcox/vc3ba_2_5"
;frun="/raid4/tcox/vc5vc5a"
;frun="/raid4/tcox/As/A5"


; Galaxy 6
;----------
;frun="/raid4/tcox/vc6vc6a"


;----------
;xyouts, 0.70, 0.90, 'bulge', /normal, charthick=3, size=1.33, color=0


;process_one_wind, "/raid4/tcox/vc3vc3e_2", lcolor=150, lpsym= 5, msg='std', x0= 0.54, y0= 0.34

;process_one_wind, "/raid4/tcox/sbw/sb8", lcolor=200, msg='sb8 - small sfr', x0= 0.54, y0= 0.30
;process_one_wind, "/raid4/tcox/sbw/sb10", lcolor=200, msg='sb10', x0= 0.54, y0= 0.30
;process_one_wind, "/raid4/tcox/sbw/sb15", lcolor=100, msg='sb15', x0= 0.54, y0= 0.26
;process_one_wind, "/raid4/tcox/sbw/sb13", lcolor=50, msg='sb13 - osc. sfr', x0= 0.54, y0= 0.22
;process_one_wind, "/raid4/tcox/sbw/sb16", lcolor=50, msg='sb16', x0= 0.54, y0= 0.22
;process_one_wind, "/raid4/tcox/sbw/sb1", lcolor=0, msg='sb1 - no change', x0= 0.54, y0= 0.18
;process_one_wind, "/raid4/tcox/sbw/sb17", lcolor=0, msg='sb17', x0= 0.54, y0= 0.18

; eta= 0.05
;process_one_wind, "/raid4/tcox/sbw/sb9", lcolor=200, msg='sb9', x0= 0.54, y0= 0.30
;process_one_wind, "/raid4/tcox/sbw/sb18", lcolor=100, msg='sb18', x0= 0.54, y0= 0.26
;process_one_wind, "/raid4/tcox/sbw/sb19", lcolor=50, msg='sb19', x0= 0.54, y0= 0.22
;process_one_wind, "/raid4/tcox/sbw/sb7", lcolor=0, msg='sb7', x0= 0.54, y0= 0.18

; eta= 0.5
;process_one_wind, "/raid4/tcox/sbw/sb10", lcolor=200, msg='sb10', x0= 0.54, y0= 0.30
;process_one_wind, "/raid4/tcox/sbw/sb15", lcolor=100, msg='sb15', x0= 0.54, y0= 0.26
;process_one_wind, "/raid4/tcox/sbw/sb16", lcolor=50, msg='sb16', x0= 0.54, y0= 0.22
;process_one_wind, "/raid4/tcox/sbw/sb17", lcolor=0, msg='sb17', x0= 0.54, y0= 0.18

;process_one_wind, "/raid4/tcox/sbw/sb17", lcolor=0, msg='sb17', x0= 0.54, y0= 0.18
;process_one_wind, "/raid4/tcox/sbw/sb17", lcolor=0, msg='sb17BH_2', x0= 0.54, y0= 0.18


; eta= 2.0
;process_one_wind, "/raid4/tcox/sbw/sb8", lcolor=200, msg='sb8', x0= 0.54, y0= 0.30
;process_one_wind, "/raid4/tcox/sbw/sb12", lcolor=100, msg='sb12', x0= 0.54, y0= 0.26
;process_one_wind, "/raid4/tcox/sbw/sb13", lcolor=50, msg='sb13', x0= 0.54, y0= 0.22
;process_one_wind, "/raid4/tcox/sbw/sb14", lcolor=0, msg='sb14', x0= 0.54, y0= 0.18


process_one_wind, "/raid4/tcox/bs/b5e", lcolor=0, msg='no wind', x0= 0.54, y0= 0.30
process_one_wind, "/raid4/tcox/sb10_mass/b5e", lcolor=150, msg='sb10', x0= 0.54, y0= 0.26
process_one_wind, "/raid4/tcox/sb8_mass/b5e", lcolor=100, msg='sb8', x0= 0.54, y0= 0.22
process_one_wind, "/raid4/tcox/sb13_mass/b5e", lcolor=50, msg='sb13', x0= 0.54, y0= 0.18




device, /close




end






pro process_one_wind, frun, lcolor=lcolor, lpsym=lpsym, $
				msg=msg, x0=x0, y0=y0

	if not keyword_set(lcolor) then lcolor= 0
	if not keyword_set(lpsym) then lpsym= 3


	read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
	idx= where(mass_egt le 0.0)
	if idx(0) ne -1 then mass_egt(idx)= 1e-6
	mass_egt= 10.0 + alog10(mass_egt/0.7)
	;mass_gtx= 10.0 + alog10(mass_gtx/0.7)
	time=time/0.7

	oplot, time, mass_egt, thick=5.0, psym=-lpsym, color= lcolor
	;oplot, time, mass_gtx, thick=5.0, psym=-lpsym, color= lcolor, linestyle=1
	;oplot, time, mass_gtx, thick=5.0, psym=-lpsym, color= lcolor

	if not keyword_set(x0) then x0= 0.75
	if not keyword_set(y0) then y0= 0.85

	if keyword_set(msg) then begin
		xyouts, x0, y0, msg, /normal, charthick=3, size=1.33, color=lcolor
	endif else begin
		xyouts, x0, y0, fload_getid(frun), /normal, charthick=1, size=1.33, color=150
	endelse



end








;==========================================================================




