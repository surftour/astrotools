




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
Bmag= mags[2,*]
V= mags[3,*]
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


;select_thispoint, pointselection, thispsym, thiscolor

if pointselection eq 11 then thiscolor= 150
if pointselection eq 12 then thiscolor= 50

thispsym= 3

if thispsym eq 3 then thisthick= 4.0 else thisthick= 1.0

;oplot, time, u_sdss-r_sdss, thick=thisthick, psym=-thispsym, color= thiscolor
oplot, time, Bmag-V, thick=thisthick, psym=-thispsym, color= thiscolor

if not keyword_set(msg) then begin
        xyouts, x0, y0, fload_getid(frun), /normal, charthick=1, size=1.33, color=thiscolor
endif else begin
        xyouts, x0, y0, msg, /normal, charthick=1, size=1.33, color=thiscolor
endelse


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

filename='color.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4, newxsize=9.0, newysize=16.0

;-----------------------------------------

;frun1= "/raid4/tcox/bs/b3e"
;frun2= "/raid4/tcox/bs/b3e_no"
;cmt= "80% gas"

frun1= "/raid4/tcox/vc3vc3e_2"
frun2= "/raid4/tcox/vc3vc3e_no"
cmt= '40% gas'

;-----------------------------------------
;
;   Load the panels
;
;-----------------------------------------
;
;     |----------|
;     |          |
;     |          |
;     |   SFR    |
;     |          |
;     |          |
;     |----------|
;     |          |
;     |          |
;     |   B-V    |
;     |          |
;     |          |
;     |----------|
;
;

x0= 0.20
xs= 1.0*(0.98-x0)
x1= x0+xs

y0= 0.10
ys= 0.5*(0.98-y0)
y1= y0+ys
y2= y0+ys+ys



;---------------------------

;yaxistitle = "!6u-r"
;yaxistitle = "U-V"
yaxistitle = "B-V"
;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"

xmax = 4.25
xmin = 0.0
ymax = 2.2
ymin = 0.0




;   1
; -----

!p.position= [x0, y1, x1, y2]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle


process_onesim_gmu, frun1, 12, 0.3, 0.68, msg='with BH'
process_onesim_gmu, frun2, 11, 0.3, 0.64, msg='no BH'

xyouts, x0+0.03, y2-0.06, cmt, /normal, charthick=3, size=1.33, color=0

fload_newcolortable, 4
oplot, [2.3,2.7], [0.6,0.6], psym=-3, color=150, thick= 8.0, symsize=1.0
xyouts, 2.9, 0.58, 'no BH', /data, charthick=3, size=1.33, color=150

oplot, [2.3,2.7], [0.4,0.4], psym=-3, color=50, thick= 4.0
xyouts, 2.9, 0.38, 'BH', /data, charthick=3, size=1.2, color=50





;   3
; -----

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
ymax = 300.0
ymin = 0.001

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, color= 0, xstyle=1, ystyle=1, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], charthick=3.0, /ylog, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
	ytickformat='exp_label', $
        xtitle=xaxistitle, ytitle=yaxistitle, /noerase


process_one_sfr, frun1, lcolor=50, lthick= 2.0, msg=' ', x0= 0.25, y0= 0.28, h=0.7
process_one_sfr, frun2, lcolor=150, lthick= 2.0, msg=' ', x0= 0.25, y0= 0.28, h=0.7

;xyouts, x1-0.15, y1-0.21, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0
;xyouts, x1-0.15, y1-0.26, '!8v!6!DW!N=837', /normal, charthick=3, size=1.33, color=0

;fload_newcolortable, 1
;oplot, [2.3,2.7], [0.6,0.6], psym=-3, color=20, thick= 8.0
;xyouts, 2.9, 0.58, 'no BH', /data, charthick=3, size=1.2, color=20

;oplot, [2.3,2.7], [0.4,0.4], psym=-4, color=200, thick= 4.0
;xyouts, 2.9, 0.38, 'BH', /data, charthick=3, size=1.2, color=200





;  Done
; ------
device, /close




end






;==============================================================================================







