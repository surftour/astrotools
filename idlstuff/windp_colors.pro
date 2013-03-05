




;-------------------------------------------------
;  g-u color evolution during a merger
;   (motivated by figure 3. of Springel
;     di Matteo, Hernquist red gals)
;-------------------------------------------------



pro read_colors_file, frun, time, mags, $
			funky=funky

cfile= frun+'/colors.txt'
;if keyword_set(funky) then cfile= frun+'/colors.txt.disk=1.0'
if keyword_set(funky) then cfile= frun+'/colors_dt.txt'

spawn, "wc "+cfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then c_data= fltarr(15,lines)

openr, 1, cfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, c_data
close, 1


time= c_data[0,*]
mags= c_data[1:14,*]

end






;==================================================================================


pro read_outputlist, OListTimes

outputlistfile= '/raid4/tcox/vc3vc3e_2/outputlist.txt'

spawn, "wc "+outputlistfile,result
lines=long(result)
lines=lines(0)
if lines GT 0 then c_data= fltarr(lines)

openr, 1, outputlistfile
readf, 1, c_data
close, 1


OListTimes= c_data[*]

end






;==================================================================================





; do the dirty work
; --------------------
pro process_onesim_gmu, frun, pointselection, x0, y0, $
                                funky=funky, msg=msg

read_colors_file, frun, time, mags, funky=funky

;-----------------
if keyword_set(funky) then begin
    read_outputlist, OListTimes
    nlist= n_elements(OListTimes)
    oldtime= time
    ncolors= n_elements(oldtime)
    oldmags= mags
    time= fltarr(nlist)
    mags= fltarr(14,nlist)

    for i=0, nlist-1 do begin
	idx= where(oldtime ge OListTimes[i])
	if OListTimes[i] gt oldtime[ncolors-1] then idx= [ncolors-1]
	time[i]= oldtime(idx(0))
	mags[*,i]= oldmags[*,idx(0)]
    endfor

endif
;-----------------

time= time/0.7

symthick= 3.0

msg=' '


;bolo= mags[0,*]
;U= mags[1,*]
;Bmag= mags[2,*]
;V= mags[3,*]
;R= mags[4,*]
;I= mags[5,*]
;J= mags[6,*]
;H= mags[7,*]
;Kmag= mags[8,*]
u_sdss= mags[9,*]
;g_sdss= mags[10,*]
r_sdss= mags[11,*]
;i_sdss= mags[12,*]
;z_sdss= mags[13,*]


select_thispoint, pointselection, thispsym, thiscolor


if thispsym eq 3 then thisthick= 4.0 else thisthick= 1.0

oplot, time, u_sdss-r_sdss, thick=thisthick, psym=-thispsym, color= thiscolor

if not keyword_set(msg) then begin
        xyouts, x0, y0, fload_getid(frun), /normal, charthick=1, size=1.33, color=thiscolor
endif else begin
        xyouts, x0, y0, msg, /normal, charthick=1, size=1.33, color=thiscolor
endelse


end






;==================================================================================






pro c1, junk


if not keyword_set(junk) then begin
        print, " "
        print, " color1, junk"
        print, " "
        print, " "
        return
endif

;filename=frun+'/color.eps'
filename='color.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;-----------------------------------------

x0= 0.18
x1= 0.98

y0= 0.15
y1= 0.98



;---------------------------

yaxistitle = "!6u-r"
;yaxistitle = "U-V"
;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"

xmax = 4.25
;xmax = 3.0
xmin = 0.0

ymax = 2.2
;ymax = 1.5
ymin = 0.0
;ymin = -0.5

;---------------------------


!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /noerase


;---------------------------

;process_onesim_gmu, "/raid4/tcox/vc3vc3e_2", 12, 0.3, 0.68, msg=' '

;oplot, [2.25,2.55], [0.68,0.68], psym=-3, color=150, thick= 8.0, symsize=1.0
;xyouts, x1-0.30, y0+0.25, '!6BH', /normal, charthick=3, size=1.2, color=150


;---------------------------

process_onesim_gmu, "/raid4/tcox/vc3vc3e_no", 11, 0.3, 0.64, msg=' '

oplot, [2.25,2.55], [0.55,0.55], psym=-3, color=0, thick= 8.0
xyouts, x1-0.30, y0+0.20, '!6no SB winds', /normal, charthick=3, size=1.2, color=0


;---------------------------

process_onesim_gmu, "/raid4/tcox/sbw/sb10", 1, 0.6, 0.68, msg=' '

usersym,1.0*[-1,-1,1,1,-1],1.0*[-1,1,1,-1,-1],thick=4.0
oplot, [2.3,2.5], [0.42,0.42], psym=-8, color=150, thick= 8.0
xyouts, x1-0.30, y0+0.15, '!7g!6=0.5, !8v!6!DW!N=837', /normal, charthick=3, size=1.2, color=150


;---------------------------

process_onesim_gmu, "/raid4/tcox/sbw/sb8", 2, 0.3, 0.28, msg=' '

oplot, [2.3,2.5], [0.29,0.29], psym=-7, color=100, thick= 8.0
xyouts, x1-0.30, y0+0.10, '!7g!6=2.0, !8v!6!DW!N=837', /normal, charthick=3, size=1.2, color=100


;---------------------------

process_onesim_gmu, "/raid4/tcox/sbw/sb13", 3, 0.6, 0.28, msg=' '

usersym,1.0*cos(findgen(49)/49*2*!pi),1.0*sin(findgen(49)/49*2*!pi), thick=4.0
oplot, [2.3,2.5], [0.15,0.15], psym=-8, color=50, thick= 8.0
xyouts, x1-0.30, y0+0.05, '!7g!6=2.0, !8v!6!DW!N=209', /normal, charthick=3, size=1.2, color=50





;  Done
; ------
device, /close




end





;==================================================================================




pro c2, junk


if not keyword_set(junk) then begin
        print, " "
        print, " c2, junk"
        print, " "
        print, " "
        return
endif

;filename=frun+'/color.eps'
filename='color.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize= 12.0, newysize= 16.0


;-----------------------------------------

x0= 0.16
x1= 0.98

y0= 0.10
ys= (0.99-y0)/2.0
y1= y0+ys
y2= y0+ys+ys

;---------------------------

yaxistitle = "!6u-r"
xaxistitle = "!6Time (Gyr)"

xmax = 4.25
xmin = 0.0
ymax = 2.2
ymin = 0.0

;---------------------------


!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /noerase


;---------------------------
process_onesim_gmu, "/raid4/tcox/vc3vc3e_2", 12, 0.3, 0.68, msg=' '
oplot, [2.25,2.55], [0.78,0.78], psym=-3, color=220, thick= 8.0, symsize=1.0
xyouts, x1-0.30, y1+0.15, '!6BH', /normal, charthick=3, size=1.2, color=220
;---------------------------
process_onesim_gmu, "/raid4/tcox/vc3vc3e_no", 11, 0.3, 0.64, msg=' '
oplot, [2.25,2.55], [0.61,0.61], psym=-3, color=0, thick= 8.0
xyouts, x1-0.30, y1+0.12, '!6no SB winds', /normal, charthick=3, size=1.2, color=0
;---------------------------
process_onesim_gmu, "/raid4/tcox/sbw/sb10", 1, 0.6, 0.68, msg=' '
usersym,1.0*[-1,-1,1,1,-1],1.0*[-1,1,1,-1,-1],thick=4.0
oplot, [2.3,2.5], [0.48,0.48], psym=-8, color=150, thick= 8.0
xyouts, x1-0.30, y1+0.09, '!7g!6=0.5, !8v!6!DW!N=837', /normal, charthick=3, size=1.2, color=150
;---------------------------
process_onesim_gmu, "/raid4/tcox/sbw/sb8", 2, 0.3, 0.28, msg=' '
oplot, [2.3,2.5], [0.31,0.31], psym=-7, color=100, thick= 8.0
xyouts, x1-0.30, y1+0.06, '!7g!6=2.0, !8v!6!DW!N=837', /normal, charthick=3, size=1.2, color=100
;---------------------------
process_onesim_gmu, "/raid4/tcox/sbw/sb13", 3, 0.6, 0.28, msg=' '
usersym,1.0*cos(findgen(49)/49*2*!pi),1.0*sin(findgen(49)/49*2*!pi), thick=4.0
oplot, [2.3,2.5], [0.16,0.16], psym=-8, color=50, thick= 8.0
xyouts, x1-0.30, y1+0.03, '!7g!6=2.0, !8v!6!DW!N=209', /normal, charthick=3, size=1.2, color=50


xyouts, x0+0.05, y2-0.05, '!8f!6!Dg!N=0.4', /normal, charthick=3, size=1.8, color=0


;---------------------------


!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /noerase


;---------------------------
process_onesim_gmu, "/raid4/tcox/bs/b3e", 12, 0.3, 0.68, msg=' ', /funky
;oplot, [2.25,2.55], [0.68,0.68], psym=-3, color=150, thick= 8.0, symsize=1.0
;xyouts, x1-0.30, y0+0.25, '!6BH', /normal, charthick=3, size=1.2, color=150
;---------------------------
process_onesim_gmu, "/raid4/tcox/bs/b3e_no", 11, 0.3, 0.64, msg=' ', /funky
;oplot, [2.25,2.55], [0.55,0.55], psym=-3, color=0, thick= 8.0
;xyouts, x1-0.30, y0+0.20, '!6no SB winds', /normal, charthick=3, size=1.2, color=0
;---------------------------
process_onesim_gmu, "/raid4/tcox/sb10_mass/b3e_no", 1, 0.6, 0.68, msg=' ', /funky
;usersym,1.0*[-1,-1,1,1,-1],1.0*[-1,1,1,-1,-1],thick=4.0
;oplot, [2.3,2.5], [0.42,0.42], psym=-8, color=150, thick= 8.0
;xyouts, x1-0.30, y0+0.15, '!7g!6=0.5, !8v!6!DW!N=837', /normal, charthick=3, size=1.2, color=150
;---------------------------
process_onesim_gmu, "/raid4/tcox/sb8_mass/b3e_no", 2, 0.3, 0.28, msg=' ', /funky
;oplot, [2.3,2.5], [0.29,0.29], psym=-7, color=100, thick= 8.0
;xyouts, x1-0.30, y0+0.10, '!7g!6=2.0, !8v!6!DW!N=837', /normal, charthick=3, size=1.2, color=100
;---------------------------
process_onesim_gmu, "/raid4/tcox/sb13_mass/b3e_no", 3, 0.6, 0.28, msg=' ', /funky
;usersym,1.0*cos(findgen(49)/49*2*!pi),1.0*sin(findgen(49)/49*2*!pi), thick=4.0
;oplot, [2.3,2.5], [0.15,0.15], psym=-8, color=50, thick= 8.0
;xyouts, x1-0.30, y0+0.05, '!7g!6=2.0, !8v!6!DW!N=209', /normal, charthick=3, size=1.2, color=50




xyouts, x0+0.05, y1-0.05, '!8f!6!Dg!N=0.8', /normal, charthick=3, size=1.8, color=0




;  Done
; ------
device, /close




end





;==================================================================================





pro c4, junk


if not keyword_set(junk) then begin
        print, " "
        print, " color4, junk"
        print, " "
        print, " "
        return
endif

;filename=frun+'/color.eps'
filename='color.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=18.0, newysize=16.0



;-----------------------------------------
;
;   Load the panels
;
;-----------------------------------------
;
;     |---------------------|
;     |          |          |
;     |          |          |
;     |    1     |    2     |
;     |          |          |
;     |          |          |
;     |---------------------|
;     |          |          |
;     |          |          |
;     |    3     |    4     |
;     |          |          |
;     |          |          |
;     |---------------------|
;
;

x0= 0.12
xs= 0.5*(0.98-x0)
x1= x0+xs
x2= x0+xs+xs

y0= 0.10
ys= 0.5*(0.98-y0)
y1= y0+ys
y2= y0+ys+ys



;---------------------------

yaxistitle = "!6u-r"
;yaxistitle = "U-V"
;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"

xmax = 4.25
;xmax = 3.0
xmin = 0.0
ymax = 2.2
;ymax = 1.5
ymin = 0.0
;ymin = -0.5




;   1
; -----

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle


process_onesim_gmu, "/raid4/tcox/vc3vc3e_2", 12, 0.3, 0.68, msg='std (with BH)'
process_onesim_gmu, "/raid4/tcox/vc3vc3e_no", 11, 0.3, 0.64, msg='std (no BH)'

xyouts, x0+0.03, y2-0.06, '!6no SB winds', /normal, charthick=3, size=1.33, color=0

fload_newcolortable, 4
oplot, [2.3,2.7], [0.6,0.6], psym=-3, color=200, thick= 8.0, symsize=1.0
xyouts, 2.9, 0.58, 'no BH', /data, charthick=3, size=1.33, color=200

oplot, [2.3,2.7], [0.4,0.4], psym=-4, color=150, thick= 4.0
xyouts, 2.9, 0.38, 'BH', /data, charthick=3, size=1.2, color=150



;   2
; -----

!p.position= [x1, y1, x2, y2]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


fload_newcolortable, 1
process_onesim_gmu, "/raid4/tcox/sbw/sb10", 10, 0.6, 0.68, msg='no BH'
process_onesim_gmu, "/raid4/tcox/sbw/sb10BH", 7, 0.6, 0.64, msg='BH'


xyouts, x2-0.15, y2-0.21, '!7g!6=0.5', /normal, charthick=3, size=1.33, color=0
xyouts, x2-0.15, y2-0.26, '!8v!6!DW!N=837', /normal, charthick=3, size=1.33, color=0

fload_newcolortable, 1
oplot, [2.3,2.7], [0.6,0.6], psym=-3, color=20, thick= 8.0
xyouts, 2.9, 0.58, 'no BH', /data, charthick=3, size=1.2, color=20

oplot, [2.3,2.7], [0.4,0.4], psym=-4, color=200, thick= 4.0
xyouts, 2.9, 0.38, 'BH', /data, charthick=3, size=1.2, color=200




;   3
; -----

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /noerase


process_onesim_gmu, "/raid4/tcox/sbw/sb8", 10, 0.3, 0.28, msg='no BH'
process_onesim_gmu, "/raid4/tcox/sbw/sb8BH", 7, 0.3, 0.24, msg='BH'

xyouts, x1-0.15, y1-0.21, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0
xyouts, x1-0.15, y1-0.26, '!8v!6!DW!N=837', /normal, charthick=3, size=1.33, color=0

fload_newcolortable, 1
oplot, [2.3,2.7], [0.6,0.6], psym=-3, color=20, thick= 8.0
xyouts, 2.9, 0.58, 'no BH', /data, charthick=3, size=1.2, color=20

oplot, [2.3,2.7], [0.4,0.4], psym=-4, color=200, thick= 4.0
xyouts, 2.9, 0.38, 'BH', /data, charthick=3, size=1.2, color=200






;   4
; -----

!p.position= [x1, y0, x2, y1]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtitle=xaxistitle, ytickformat='(a1)', /noerase


process_onesim_gmu, "/raid4/tcox/sbw/sb13", 10, 0.6, 0.28, msg='no BH'
process_onesim_gmu, "/raid4/tcox/sbw/sb13BH", 7, 0.6, 0.24, msg='BH'


xyouts, x2-0.15, y1-0.21, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0
xyouts, x2-0.15, y1-0.26, '!8v!6!DW!N=209', /normal, charthick=3, size=1.33, color=0

fload_newcolortable, 1
oplot, [2.3,2.7], [0.6,0.6], psym=-3, color=20, thick= 8.0
xyouts, 2.9, 0.58, 'no BH', /data, charthick=3, size=1.2, color=20

oplot, [2.3,2.7], [0.4,0.4], psym=-4, color=200, thick= 4.0
xyouts, 2.9, 0.38, 'BH', /data, charthick=3, size=1.2, color=200





;  Done
; ------
device, /close




end






;==============================================================================================




pro c6, junk


if not keyword_set(junk) then begin
        print, " "
        print, " c2, junk"
        print, " "
        print, " "
        return
endif

filename='color.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=20.0, newysize=16.0



;-----------------------------------------
;
;   Load the panels
;
;-----------------------------------------

;     |--------------------------------|
;     |          |          |          |
;     |          |          |          |
;     |    1     |    2     |    3     |
;     |          |          |          |
;     |          |          |          |
;     |---------------------|----------|
;     |          |          |          |
;     |          |          |          |
;     |    4     |    5     |    6     |
;     |          |          |          |
;     |          |          |          |
;     |---------------------|----------|
;
;

x0= 0.10
xs= (0.99-x0) * (1./3.)
x1= x0+xs
x2= x0+xs+xs
x3= x0+xs+xs+xs

y0= 0.09
ys= (0.99-y0) * (1./3.)
y1= y0+ys
y2= y0+ys+ys
y3= y0+ys+ys+ys



;---------------------------

yaxistitle = "!6u-r"
;yaxistitle = "U-V"
;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"

xmax = 4.25
;xmax = 3.0
xmin = 0.0
ymax = 2.2
;ymax = 1.5
ymin = 0.0
;ymin = -0.5

;---------------------------




; -----------------------------------------------------------------------
;  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
; -----------------------------------------------------------------------

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle


;---------------------------
;process_onesim_gmu, "/raid4/tcox/vc3vc3e_2", 12, 0.3, 0.68, msg=' '
;oplot, [2.25,2.55], [0.68,0.68], psym=-3, color=150, thick= 8.0, symsize=1.0
;xyouts, x1-0.30, y0+0.25, '!6BH', /normal, charthick=3, size=1.2, color=150
;---------------------------
process_onesim_gmu, "/raid4/tcox/vc3vc3e_no", 11, 0.3, 0.64, msg=' '
oplot, [2.25,2.55], [0.55,0.55], psym=-3, color=0, thick= 8.0
xyouts, 2.68, 0.53, '!6no SB winds', /data, charthick=3, size=1.2, color=0
;---------------------------
process_onesim_gmu, "/raid4/tcox/sbw/sb10", 1, 0.6, 0.68, msg=' '
usersym,1.0*[-1,-1,1,1,-1],1.0*[-1,1,1,-1,-1],thick=4.0
oplot, [2.3,2.5], [0.42,0.42], psym=-8, color=150, thick= 8.0
xyouts, 2.63, 0.40, '!7g!6=0.5, !8v!6!DW!N=837', /data, charthick=3, size=1.2, color=150
;---------------------------
process_onesim_gmu, "/raid4/tcox/sbw/sb8", 2, 0.3, 0.28, msg=' '
oplot, [2.3,2.5], [0.29,0.29], psym=-7, color=100, thick= 8.0
xyouts, 2.63, 0.27, '!7g!6=2.0, !8v!6!DW!N=837', /data, charthick=3, size=1.2, color=100
;---------------------------
process_onesim_gmu, "/raid4/tcox/sbw/sb13", 3, 0.6, 0.28, msg=' '
usersym,1.0*cos(findgen(49)/49*2*!pi),1.0*sin(findgen(49)/49*2*!pi), thick=4.0
oplot, [2.3,2.5], [0.15,0.15], psym=-8, color=50, thick= 8.0
xyouts, 2.63, 0.13, '!7g!6=2.0, !8v!6!DW!N=209', /data, charthick=3, size=1.2, color=50





; -----------------------------------------------------------------------
;  2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2
; -----------------------------------------------------------------------

!p.position= [x1, y1, x2, y2]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtickformat='(a1)', ytickformat='(a1)', /noerase


fload_newcolortable, 1
process_onesim_gmu, "/raid4/tcox/sbw/sb10", 10, 0.6, 0.68, msg='no BH'
process_onesim_gmu, "/raid4/tcox/sbw/sb10BH", 7, 0.6, 0.64, msg='BH'


xyouts, x2-0.15, y2-0.21, '!7g!6=0.5', /normal, charthick=3, size=1.33, color=0
xyouts, x2-0.15, y2-0.26, '!8v!6!DW!N=837', /normal, charthick=3, size=1.33, color=0

fload_newcolortable, 1
oplot, [2.3,2.7], [0.6,0.6], psym=-3, color=20, thick= 8.0
xyouts, 2.9, 0.58, 'no BH', /data, charthick=3, size=1.2, color=20

oplot, [2.3,2.7], [0.4,0.4], psym=-4, color=200, thick= 4.0
xyouts, 2.9, 0.38, 'BH', /data, charthick=3, size=1.2, color=200




; -----------------------------------------------------------------------
;  3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3
; -----------------------------------------------------------------------

!p.position= [x2, y1, x3, y2]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /noerase


process_onesim_gmu, "/raid4/tcox/sbw/sb8", 10, 0.3, 0.28, msg='no BH'
process_onesim_gmu, "/raid4/tcox/sbw/sb8BH", 7, 0.3, 0.24, msg='BH'

xyouts, x1-0.15, y1-0.21, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0
xyouts, x1-0.15, y1-0.26, '!8v!6!DW!N=837', /normal, charthick=3, size=1.33, color=0

fload_newcolortable, 1
oplot, [2.3,2.7], [0.6,0.6], psym=-3, color=20, thick= 8.0
xyouts, 2.9, 0.58, 'no BH', /data, charthick=3, size=1.2, color=20

oplot, [2.3,2.7], [0.4,0.4], psym=-4, color=200, thick= 4.0
xyouts, 2.9, 0.38, 'BH', /data, charthick=3, size=1.2, color=200






; -----------------------------------------------------------------------
;  4   4   4   4   4   4   4   4   4   4   4   4   4   4   4   4   4   4
; -----------------------------------------------------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtitle=xaxistitle, ytickformat='(a1)', /noerase


process_onesim_gmu, "/raid4/tcox/sbw/sb13", 10, 0.6, 0.28, msg='no BH'
process_onesim_gmu, "/raid4/tcox/sbw/sb13BH", 7, 0.6, 0.24, msg='BH'


xyouts, x2-0.15, y1-0.21, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0
xyouts, x2-0.15, y1-0.26, '!8v!6!DW!N=209', /normal, charthick=3, size=1.33, color=0

fload_newcolortable, 1
oplot, [2.3,2.7], [0.6,0.6], psym=-3, color=20, thick= 8.0
xyouts, 2.9, 0.58, 'no BH', /data, charthick=3, size=1.2, color=20

oplot, [2.3,2.7], [0.4,0.4], psym=-4, color=200, thick= 4.0
xyouts, 2.9, 0.38, 'BH', /data, charthick=3, size=1.2, color=200




; -----------------------------------------------------------------------
;  5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5
; -----------------------------------------------------------------------

!p.position= [x1, y0, x2, y1]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtitle=xaxistitle, ytickformat='(a1)', /noerase





; -----------------------------------------------------------------------
;  6   6   6   6   6   6   6   6   6   6   6   6   6   6   6   6   6   6
; -----------------------------------------------------------------------

!p.position= [x2, y0, x3, y1]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtitle=xaxistitle, ytickformat='(a1)', /noerase





;  Done
; ------
device, /close




end



