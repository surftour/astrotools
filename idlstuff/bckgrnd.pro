pro background, junk


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename='background.eps', colortable= 0, newxsize=14, newysize=14


;loadct, 0
;loadct, 1
;loadct, 4
;tvlct,r,g,b,/get
;v1=[0,255]
;v2=[0,255]
;v3=[0,255]
;tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)

;tvlct,r1,g1,b1,/get     ; reget new color table with the black and white



; png file
; ---------
;read_png, 'background/large.png', pngimage, r1, g1, b1



; fits file
; ----------
fitsimg= readfits('background/NGC_1637_1.fits')

;newfitsimg= 256.0 * fitsimg / max(fitsimg)
newfitsimg= 256.0 * (fitsimg-min(fitsimg)) / (max(fitsimg)-min(fitsimg))
newfitsimg= 256-newfitsimg

; cut out galaxy
newfitsimg[425:620,455:620]= newfitsimg[225:420,255:420]

tv, newfitsimg, 0.0, 0.0, xsize=1.0, ysize=1.0, /normal


device, /close



end


