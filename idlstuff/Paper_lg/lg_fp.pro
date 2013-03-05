;------------------------------------------------------------------------
;
;     Plot FP of LG Remnant
;   ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro lgfp, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "lgfp, junk"
   print, "  "
   return 
endif

filename='lgfp.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


;--------------------------------------
;--------------------------------------


xaxistitle= '!6Log !7r!6 (km/s) + 0.205 <!7m!6!DK!N>!Deff!N'
xmax= 7.0
xmin= 4.5
;xaxistitle= '!6Log !7r!6 (km/s) + 0.516 Log!4R!6!Deff!N'
;xmax= 8.4
;xmin= 5.5

yaxistitle= '!6Log R!Deff!N (kpc)'
ymax = 1.5
ymin = -1.0


bins= 40



re= 3.44 / 0.7    ; in kpc
sigma= 175.3         ; in km/sec
Mstar= 9.53e+10 / 0.7  ; in m_sun




;----------------------------------------
; Generate plot 
;----------------------------------------

!p.position= [0.2, 0.15, 0.98, 0.98]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, $
        xrange=[xmin,xmax],$
        yrange=[ymin,ymax],$
        color=0, $
        xcharsize=1.5, ycharsize=1.5, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;xticks= 14, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /noerase, $
        /nodata




; --------------------

re= 3.44 / 0.7    ; in kpc
sigma= 175.3         ; in km/sec
Mstar= 9.53e+10 / 0.7  ; in m_sun
y= [alog10(re)]
;x= [alog10(sigma) + 0.516*alog10(Mstar/re/re/2/!PI)]
x= [5.75]

oplot, x, y, psym=5, symsize=3.0, color=150, thick=8.0

; --------------------

; near-IR FP from Pahre et al. (1998)
;re=[0.01,0.1,1.0,10.0,100.0]
x= findgen(10)
normalization= 8.3    ; if x is surface brightness
;normalization= 10.1    ; if x is surface density
y= 1.53 * x - normalization
oplot, x, y, psym=-3, color=0, thick=3.0

xyouts, 0.47, 0.35, /normal, '!6R!Deff!N !9A!6 !7r!6!E1.53!N !4R!6!Deff!N!E0.79!N', color= 0, size=1.7, charthick=4.0
xyouts, 0.5, 0.29, /normal, '(Pahre et al. 1998)', color= 0, size=1.1, charthick=3.0



xyouts, 0.8, 0.9, '(f)', /normal, size=1.8, color= 0, charthick=3.0


; --------------------

device, /close



end


