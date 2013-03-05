pro tully_fisher, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "tully_fisher, junk"
   return
endif

sendto='ps'

initialize_plotinfo, 1




;=========================================================
;  tully fisher relation
; -----------------------
;
; we'll use the baryonic TF from
;=========================================================

setup_plot_stuff, sendto, filename='tf.eps', colortable=4



;--------------------------------------
;  Print the Shit
;--------------------------------------


; baryonic mass
; 
ymax = 10.0e+12
ymin = 1.0e+8

xmax = 10.0e+2
xmin = 30.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
	ystyle=1, $
	;ystyle=8, $                     ; suppress right y axis
	/ylog, $
	/xlog, $
        xcharsize=1.5, ycharsize=1.5, $
        charthick=2.0, $
        xthick=4.0, ythick=4.0, $
	;ytickformat='noexp_label', $
	ytickformat='exp_label', $
        xtitle="V!Dc!N (km/sec)", $
        ytitle="Log M!Dbaryon!N (M!D!9n!3!N)", $
        /nodata



if sendto NE 'ps' then begin
       print, 'we are only equiped to print to a ps file' 
endif else begin

        ; vc0
        ; ---
        mbaryon= [0.168e+10]
        vc= [72.9]
        xyouts, vc, mbaryon, 'vc0', charthick=5.0, color= 40, size=1.5

        ; vc1
        ; ---
        mbaryon= [0.488e+10]
        vc= [104.0]
        xyouts, vc, mbaryon, 'vc1', charthick=5.0, color= 50, size=1.5

        ; vc2
        ; ---
        mbaryon= [1.45e+10]
        vc= [147.0]
        xyouts, vc, mbaryon, 'vc2', charthick=5.0, color= 60, size=1.5

        ; vc3
        ; ---
        mbaryon= [3.90e+10]
        vc= [208.0]
        xyouts, vc, mbaryon, 'vc3', charthick=5.0, color= 70, size=1.5

	; bulge version
        ;mbaryon= 1.33*[3.90e+10]
        ;vc= [225.0]
        ;xyouts, vc, mbaryon, '(b)', charthick=5.0, color= 70, size=1.5
        ;oplot, vc, mbaryon, psym=2, thick=5.0, color= 70, symsize=1.5

	; brant's high c bulge version
        ;mbaryon= 1.33*[3.90e+10]
        ;vc= [246.0]
        ;xyouts, vc, mbaryon, '(c)', charthick=5.0, color= 70, size=1.5
        ;oplot, vc, mbaryon, psym=5, thick=5.0, color= 70, symsize=1.5

	; vc4
	; ----
	mbaryon= [11.00e+10]
	vc= [310.0]
        xyouts, vc, mbaryon, 'vc4', charthick=5.0, color= 90, size=1.5

	; vc5
	; ---
        mbaryon= [31.24e+10]
        vc= [418.0]
	xyouts, vc, mbaryon, 'vc5', charthick=5.0, color= 110, size=1.5

        ; vc6
        ; ---
        mbaryon= [119.167e+10]
        vc= [617.0]
        xyouts, vc, mbaryon, 'vc6', charthick=5.0, color= 130, size=1.5



endelse


;--------------------------------------
; plot Bell & de Jong baryonic TF
; 
alpha= 3.51
M_100= 10^(9.79-(2*alpha))

x=[xmin,xmax]
y= M_100 * (x^alpha)
oplot, x, y, psym=-3, thick=3.0, linestyle=3, color= 0


;--------------------------------------
;--------------------------------------

if (sendto EQ 'ps') then device, /close








end







