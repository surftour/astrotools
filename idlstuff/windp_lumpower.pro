
;==========================================
;
;    Luminosity Contribution
;
;
;
;
;==========================================






; --------------------------------
;  Read FB Energy File
; ----------------------------------
pro read_fbenergy_file, frun, time, fb_std, fb_wind, fb_bh, $
				fbfile=fbfile

if keyword_set(fbfile) then fbfile= frun+'/'+fbfile else fbfile= frun+'/fb_energy.txt'

spawn, "wc "+fbfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then re_data= fltarr(4,lines)

openr, 1, fbfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, re_data
close, 1


time= re_data[0,*]
fb_std= re_data[1,*]
fb_wind= re_data[2,*]
fb_bh= re_data[3,*]


end




;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------






; --------------------------------
;  Read Bolometric Lum File
; ----------------------------------
pro read_bololum_file, frun, time, sbololum, bhbololum, $
                bolofile=bolofile

;if not keyword_set(bolofile) then bolofile= frun+'/lum_bolo.txt'
if not keyword_set(bolofile) then bolofile= frun+'/luminosity.txt'

spawn, "wc "+bolofile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then re_data= fltarr(3,lines)

openr, 1, bolofile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, re_data
close, 1


time= re_data[0,*]
sbololum= re_data[1,*]
bhbololum= re_data[2,*]


end




;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------






;pro get_luminfo, frun, subdir=subdir, bololum, powerratio
pro get_luminfo, frun, subdir=subdir, bololum, sblum, bhlum, time=time


        thisfrun= '/raid4/tcox/'+subdir+'/'+frun
        print, "processing: "+ thisfrun
        read_bololum_file, thisfrun, time, sbololum, bhbololum, bolofile=bolofile

        this_bololum= alog10(10^(sbololum) + 10^(bhbololum))
        ;this_powerratio= 10^(bhbololum - sbololum)
	this_sblum= sbololum
	this_bhlum= bhbololum

        ; get time range for active phase
        ; --------------------------------
        ;thisfrun= '/raid4/tcox/'+subdir+'/'+frun
        t_merg= fload_bh_mergertime(thisfrun)
        ;t_merg= t_merg/0.7
        ;print, "t_merg= ", t_merg
        t_sb= t_merg-0.3
        t_postsb= t_merg-0.3+1.0
        ;t_postsb= t_merg-0.3+0.7
        ;print, "beginning of SB time= ", t_sb


        time= time / 0.7
        idx= where((time ge t_sb) and (time le t_postsb))


        bololum= this_bololum(idx)
        ;powerratio= this_powerratio(idx)
	sblum= this_sblum(idx)
	bhlum= this_bhlum(idx)
	time= time(idx)

end









;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------



pro get_luminfo_fruns, fruns, subdir=subdir, bololum, sblum, bhlum


nfruns= n_elements(fruns)
bololum= fltarr(nfruns)
sblum= fltarr(nfruns)
bhlum= fltarr(nfruns)

for i= 0, nfruns-1 do begin

    get_luminfo, fruns[i], subdir=subdir, ibololum, isblum, ibhlum

    bololum[i]= max(ibololum)
    sblum[i]= max(isblum)
    bhlum[i]= max(ibhlum)
print, bololum[i], sblum[i], bhlum[i]

endfor



end




;----------------------------------------------------------------


   
pro get_tdinfo_fruns, fruns, subdir=subdir, td


nfruns= n_elements(fruns)
bololum= fltarr(nfruns)
sblum= fltarr(nfruns)
bhlum= fltarr(nfruns)
td= fltarr(nfruns)

for i= 0, nfruns-1 do begin

    get_luminfo, fruns[i], subdir=subdir, ibololum, isblum, ibhlum, time=time

    bololum[i]= max(ibololum)

    sblum[i]= max(isblum)
    sbidx= where(isblum eq sblum[i])
    sbtime= time(sbidx)

    bhlum[i]= max(ibhlum)
    bhidx= where(ibhlum eq bhlum[i])
    bhtime= time(bhidx)

    td[i]= 1000.0*(bhtime - sbtime)

print, bololum[i], sblum[i], bhlum[i], td[i]

endfor



end




;---------------------------------------------------------------




pro avg_two, bl_1, sb_1, bh_1, bl_2, sb_2, bh_2, $
	bl_avg, bl_1sig, sb_avg, sb_1sig, bh_avg, bh_1sig


ni= n_elements(bl_1)
bl_avg= fltarr(ni)
bl_1sig= fltarr(ni)
sb_avg= fltarr(ni)
sb_1sig= fltarr(ni)
bh_avg= fltarr(ni)
bh_1sig= fltarr(ni)
for i=0, ni-1 do begin
   
   x= [bl_1[i], bl_2[i]]
   bl_avg[i]= mean(x)
   bl_1sig[i]= sqrt(variance(x))

   x= [sb_1[i], sb_2[i]]
   sb_avg[i]= mean(x)
   sb_1sig[i]= sqrt(variance(x))

   x= [bh_1[i], bh_2[i]]
   bh_avg[i]= mean(x)
   bh_1sig[i]= sqrt(variance(x))
endfor



end




;---------------------------------------------------------------




pro process_and_plot_lums, totmass, fruns1, subdir1, fruns2, subdir2, pointt, $
				msg, xpt1, xpt2, xpt3, ylbl


get_luminfo_fruns, fruns1, subdir=subdir1, bololum_80e, sblum_80e, bhlum_80e
get_luminfo_fruns, fruns2, subdir=subdir2, bololum_40e, sblum_40e, bhlum_40e


avg_two, bololum_80e, sblum_80e, bhlum_80e, bololum_40e, sblum_40e, bhlum_40e, $
        bl_avg, bl_1sig, sb_avg, sb_1sig, bh_avg, bh_1sig


select_thispoint, pointt, thispsym, thiscolor
;oplot, totmass, bl_avg, psym=-thispsym, color=thiscolor, thick=3.0
;oploterror, totmass, bl_avg, bl_1sig, color=thiscolor, errcolor=thiscolor, errthick=3.0,  thick=3.0
oplot, totmass, sb_avg, psym=-thispsym, color=thiscolor, thick=3.0
oploterror, totmass, sb_avg, sb_1sig, color=thiscolor, errcolor=thiscolor, errthick=3.0,  thick=3.0

oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, msg, size=1.2, color=thiscolor, /data, charthick=4.0


end






;---------------------------------------------------------------




pro process_and_plot_td, totmass, fruns1, subdir1, fruns2, subdir2, pointt, $
                                msg, xpt1, xpt2, xpt3, ylbl


get_tdinfo_fruns, fruns1, subdir=subdir1, td1
get_tdinfo_fruns, fruns2, subdir=subdir2, td2


select_thispoint, pointt, thispsym, thiscolor
oplot, totmass, td1, psym=-thispsym, color=thiscolor, thick=3.0
oplot, totmass, td2, psym=-thispsym, color=thiscolor, thick=3.0

oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, msg, size=1.2, color=thiscolor, /data, charthick=4.0


end









;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------





;-------------------------------------------------
;-------------------------------------------------
pro lumpower, junk


if not keyword_set(junk) then begin
        print, " "
        print, " lumpower, junk"
        print, " "
        print, " "
        return
endif

filename='lumpower.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

;---------------------------



xaxistitle = "!6Log Bolometric Luminosity (L!D!9n!6!N)"
xmax = 13.9
xmin = 8.7

yaxistitle = "!6L!DBH!N/L!DStarburst!N "
ymax = 40.0
ymin = 0.0004


;---------------------------

x0= 0.18
x1= 0.97
y0= 0.15
y1= 0.98


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
	;/xlog, $
	/ylog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        ;xtickformat='exp_label', $
        ;xtickformat='(a1)', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------

xpt1= 9.0
xpt2= 9.0
xpt3= 9.25



;---------------------------
fruns= ['d0e2_q', 'd1e2_q', 'd2e2_q', 'd3e7', 'd4e2_q', 'd5e2_q', 'd6e2_q'] & nfruns= n_elements(fruns)
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='ds', bololum, powerratio
    oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
endfor
ylbl= 9.4
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40e_q', size=1.2, color=thiscolor, /data, charthick=4.0


;---------------------------
fruns= ['d0h2_q', 'd1h2_q', 'd2h2_q', 'd3h7', 'd4h2_q', 'd5h2_q', 'd6h2_q'] & nfruns= n_elements(fruns)
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='ds', bololum, powerratio
    oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
endfor
ylbl= 7.1
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40h_q', size=1.2, color=thiscolor, /data, charthick=4.0


;---------------------------
fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e'] & nfruns= n_elements(fruns)
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='bs', bololum, powerratio
    oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
endfor
ylbl= 20.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80e', size=1.2, color=thiscolor, /data, charthick=4.0


;---------------------------
fruns=['b0h', 'b1h', 'b2h', 'b3h', 'b4h', 'b5h', 'b6h'] & nfruns= n_elements(fruns)
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='bs', bololum, powerratio
    oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
endfor
ylbl= 12.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80h', size=1.2, color=thiscolor, /data, charthick=4.0


;---------------------------
fruns= ['e0e', 'e1e', 'e2e', 'e3e', 'e4e', 'e5e', 'e6e'] & nfruns= n_elements(fruns)
pointt= 6
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='es', bololum, powerratio
    oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
endfor
ylbl= 5.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '20e', size=1.2, color=thiscolor, /data, charthick=4.0


;---------------------------
fruns= ['z0e', 'z1e', 'z2e', 'z3e', 'z4e', 'z5e', 'z6e'] & nfruns= n_elements(fruns)
pointt= 7
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='zs', bololum, powerratio
    oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
endfor
ylbl= 3.8
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '5e', size=1.2, color=thiscolor, /data, charthick=4.0


;---------------------------
fruns= ['vc3vc3e_2'] & nfruns= n_elements(fruns)
pointt= 8
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='', bololum, powerratio
    oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
endfor
ylbl= 2.8
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, 'vc3vc3e_2', size=1.2, color=thiscolor, /data, charthick=4.0


;---------------------------
fruns= ['sb10','sb10BH'] & nfruns= n_elements(fruns)
pointt= 9
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='sbw', bololum, powerratio
    oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
endfor
ylbl= 1.8
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, 'sbw', size=1.2, color=thiscolor, /data, charthick=4.0



;---------------------------




; done
; ------
device, /close




end






;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------








;------------------------------------------------------------------------------------
; ===================================================================================
;------------------------------------------------------------------------------------




pro lumpower_panel, junk


if not keyword_set(junk) then begin
        print, " "
        print, " lumpower, junk"
        print, " "
        print, " "
        return
endif 

filename='lumpower.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=22, newysize=18



x0= 0.10 & xsize= 0.178
x1= x0+xsize
x2= x0+xsize+xsize
x3= x0+xsize+xsize+xsize
x4= x0+xsize+xsize+xsize+xsize
x5= x0+xsize+xsize+xsize+xsize+xsize


y0= 0.09 & ysize=0.225
y1= y0+ysize
y2= y0+ysize+ysize
y3= y0+ysize+ysize+ysize
y4= y0+ysize+ysize+ysize+ysize


;---------------------------


    
xaxistitle = "!6Log Bolometric Luminosity (L!D!9n!6!N)"
xmax = 13.9
xmin = 8.7

yaxistitle = "!6L!DBH!N/L!DStarburst!N "
ymax = 40.0
ymin = 0.0004


xvar=[0]
yvar=[0]




;---------------------------
;   top row:  5% gas
;---------------------------

fruns= ['z1e', 'z2e', 'z3e', 'z4e', 'z5e']
fruns2= ['z1h', 'z2h', 'z3h', 'z4h', 'z5h']

pointt= 7
select_thispoint, pointt, thispsym, thiscolor

overplot_panel, xvar, yvar, x0, y3, x1, y4, xmax, xmin, ymax, ymin, yaxistitle=yaxistitle, /ylog
get_luminfo, fruns[0], subdir='zs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[0], subdir='zs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x1, y3, x2, y4, xmax, xmin, ymax, ymin, /ylog
get_luminfo, fruns[1], subdir='zs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[1], subdir='zs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x2, y3, x3, y4, xmax, xmin, ymax, ymin, /ylog
get_luminfo, fruns[2], subdir='zs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[2], subdir='zs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x3, y3, x4, y4, xmax, xmin, ymax, ymin, /ylog
get_luminfo, fruns[3], subdir='zs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[3], subdir='zs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x4, y3, x5, y4, xmax, xmin, ymax, ymin, /ylog
get_luminfo, fruns[4], subdir='zs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[4], subdir='zs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0


xyouts, 0.12, 0.95, '5%', size=1.2, color=thiscolor, /normal, charthick=4.0




;---------------------------
;   top row:  20% gas
;---------------------------

fruns= ['e1e', 'e2e', 'e3e', 'e4e', 'e5e']
fruns2= ['e1h', 'e2h', 'e3h', 'e4h', 'e5h']

pointt= 5
select_thispoint, pointt, thispsym, thiscolor

overplot_panel, xvar, yvar, x0, y2, x1, y3, xmax, xmin, ymax, ymin, yaxistitle=yaxistitle, /ylog
get_luminfo, fruns[0], subdir='es', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[0], subdir='es', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x1, y2, x2, y3, xmax, xmin, ymax, ymin, /ylog
get_luminfo, fruns[1], subdir='es', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[1], subdir='es', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x2, y2, x3, y3, xmax, xmin, ymax, ymin, /ylog
get_luminfo, fruns[2], subdir='es', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[2], subdir='es', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x3, y2, x4, y3, xmax, xmin, ymax, ymin, /ylog
get_luminfo, fruns[3], subdir='es', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[3], subdir='es', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
   
overplot_panel, xvar, yvar, x4, y2, x5, y3, xmax, xmin, ymax, ymin, /ylog
get_luminfo, fruns[4], subdir='es', bololum, powerratio 
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[4], subdir='es', bololum, powerratio 
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0


xyouts, 0.12, 0.72, '20%', size=1.2, color=thiscolor, /normal, charthick=4.0



;---------------------------
;   top row:  40% gas
;---------------------------

fruns= ['d1e2_q', 'd2e2_q', 'd3e7', 'd4e2_q', 'd5e2_q']
fruns2= ['d1h2_q', 'd2h2_q', 'd3h7', 'd4h2_q', 'd5h2_q']

pointt= 3
select_thispoint, pointt, thispsym, thiscolor

overplot_panel, xvar, yvar, x0, y1, x1, y2, xmax, xmin, ymax, ymin, yaxistitle=yaxistitle, /ylog
get_luminfo, fruns[0], subdir='ds', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[0], subdir='ds', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x1, y1, x2, y2, xmax, xmin, ymax, ymin, /ylog
get_luminfo, fruns[1], subdir='ds', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[1], subdir='ds', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x2, y1, x3, y2, xmax, xmin, ymax, ymin, /ylog
get_luminfo, fruns[2], subdir='ds', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[2], subdir='ds', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x3, y1, x4, y2, xmax, xmin, ymax, ymin, /ylog
get_luminfo, fruns[3], subdir='ds', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[3], subdir='ds', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
   
overplot_panel, xvar, yvar, x4, y1, x5, y2, xmax, xmin, ymax, ymin, /ylog
get_luminfo, fruns[4], subdir='ds', bololum, powerratio 
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[4], subdir='ds', bololum, powerratio 
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0


xyouts, 0.12, 0.50, '40%', size=1.2, color=thiscolor, /normal, charthick=4.0



;---------------------------
;   top row:  80% gas
;---------------------------

fruns= ['b1e', 'b2e', 'b3e', 'b4e', 'b5e']
fruns2= ['b1h', 'b2h', 'b3h', 'b4h', 'b5h']

pointt= 6
select_thispoint, pointt, thispsym, thiscolor
;thispsym=-3

overplot_panel, xvar, yvar, x0, y0, x1, y1, xmax, xmin, ymax, ymin, yaxistitle=yaxistitle, /ylog, xaxistitle=' '
get_luminfo, fruns[0], subdir='bs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[0], subdir='bs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x1, y0, x2, y1, xmax, xmin, ymax, ymin, /ylog, xaxistitle=' '
get_luminfo, fruns[1], subdir='bs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[1], subdir='bs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x2, y0, x3, y1, xmax, xmin, ymax, ymin, /ylog, xaxistitle=xaxistitle
get_luminfo, fruns[2], subdir='bs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[2], subdir='bs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0

overplot_panel, xvar, yvar, x3, y0, x4, y1, xmax, xmin, ymax, ymin, /ylog, xaxistitle=' '
get_luminfo, fruns[3], subdir='bs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[3], subdir='bs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
   
overplot_panel, xvar, yvar, x4, y0, x5, y1, xmax, xmin, ymax, ymin, /ylog, xaxistitle=' '
get_luminfo, fruns[4], subdir='bs', bololum, powerratio 
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0
get_luminfo, fruns2[4], subdir='bs', bololum, powerratio
oplot, bololum, powerratio, psym=thispsym, color=thiscolor, thick=3.0


xyouts, 0.12, 0.28, '80%', size=1.2, color=thiscolor, /normal, charthick=4.0







;---------------------------




; done
; ------
device, /close




end





; --------------------------------------------------------------------------------------



pro overplot_panel, xvar1, yvar1, $
				x0, y0, x1, y1, $
				xmax, xmin, ymax, ymin, $
				yaxistitle=yaxistitle, $
				xaxistitle=xaxistitle, ylog=ylog, $
				showtotal=showtotal, $
				plot_xray_lum=plot_xray_lum, $
				plabel=plabel


f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill

!p.position= [x0, y0, x1, y1]

if keyword_set(xaxistitle) and keyword_set(yaxistitle) then begin
   if keyword_set(ylog) then begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
                xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, $
                xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
		ytickformat='exp_label', $
                xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase
   endif else begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
	        xrange=[xmin,xmax], yrange=[ymin,ymax], $
	        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase
   endelse
endif

if keyword_set(xaxistitle) and (not keyword_set(yaxistitle)) then begin
   if keyword_set(ylog) then begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
                xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, $
                xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
                xtitle=xaxistitle, ytickformat='(a1)', /nodata, /noerase
   endif else begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
	        xrange=[xmin,xmax], yrange=[ymin,ymax], $
	        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	        xtitle=xaxistitle, ytickformat='(a1)', /nodata, /noerase
   endelse
endif

if not keyword_set(xaxistitle) and keyword_set(yaxistitle) then begin
   if keyword_set(ylog) then begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
                xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, $
                xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
		ytickformat='exp_label', $
                xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
   endif else begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
                xrange=[xmin,xmax], yrange=[ymin,ymax], $
                xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
                xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
   endelse
endif

;if not keyword_set(xaxistitle) and keyword_set(yaxistitle) then begin
;	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
;	        xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, $
;	        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
;	        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
;endif

if not keyword_set(xaxistitle) and (not keyword_set(yaxistitle)) then begin
   if keyword_set(ylog) then begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
                xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, $
                xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
                xtickformat='(a1)', ytickformat='(a1)', /nodata, /noerase
   endif else begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
                xrange=[xmin,xmax], yrange=[ymin,ymax], $
                xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
                xtickformat='(a1)', ytickformat='(a1)', /nodata, /noerase
   endelse
endif


x= [11.0, 11.0]
y= [ymin, ymax]
oplot, x, y, psym=-3, linestyle= 1, color= 0

x= [12.0, 12.0]
y= [ymin, ymax]
oplot, x, y, psym=-3, linestyle= 0, color= 0, thick= 3.0

x= [13.0, 13.0]
y= [ymin, ymax]
oplot, x, y, psym=-3, linestyle= 1, color= 0


;  bh lum= sf lum
; ------------------
x= [xmin, xmax]
y= [1.0, 1.0]
oplot, x, y, psym=-3, linestyle= 2, color= 0


end






;================================================================================================





;-------------------------------------------------
;-------------------------------------------------
pro peaklum, junk


if not keyword_set(junk) then begin
        print, " "
        print, " peaklum, junk"
        print, " "
        print, " "
        return
endif

filename='peaklum.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

;---------------------------


yaxistitle = "!6Log Peak Bolometric Luminosity (L!D!9n!6!N)"
ymax = 13.9
ymin = 8.7


xaxistitle = "!6Stellar Mass (M!D!9n!6!N)"
xmax = 7.0e+12
xmin = 2.0e+9


;---------------------------

x0= 0.18
x1= 0.97
y0= 0.15
y1= 0.98


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /xlog, $
        ;/ylog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtickformat='exp_label', $
        ;xtickformat='(a1)', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------

xpt1= 3.0e+11
xpt2= 3.0e+11
xpt3= 7.0e+11
    


;---------------------------


fruns= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e'] & nfruns= n_elements(fruns)
stmass= [0.291, 0.876, 2.53, 6.53, 18.0, 47.7, 212.9] * 1.0d+10 / 0.7
bololum= stmass
pointt= 1
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='bs', ibololum, powerratio
    bololum[i]= max(ibololum)
endfor
oplot, stmass, bololum, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 10.6
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80e', size=1.2, color=thiscolor, /data, charthick=4.0



;---------------------------


fruns=['b0h', 'b1h', 'b2h', 'b3h', 'b4h', 'b5h', 'b6h'] & nfruns= n_elements(fruns)
stmass= [0.302, 0.814, 2.43, 6.63, 18.7, 53.8, 219.1] * 1.0d+10 / 0.7
bololum= stmass
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='bs', ibololum, powerratio
    bololum[i]= max(ibololum)
endfor
oplot, stmass, bololum, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 10.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '80h', size=1.2, color=thiscolor, /data, charthick=4.0


;---------------------------


fruns= ['d0e2_q', 'd1e2_q', 'd2e2_q', 'd3e7', 'd4e2_q', 'd5e2_q', 'd6e2_q'] & nfruns= n_elements(fruns)
stmass= [0.310, 0.899, 2.66, 7.22, 18.9, 57.8, 227.4] * 1.0d+10 / 0.7
bololum= stmass
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='ds', ibololum, powerratio
    bololum[i]= max(ibololum)
endfor
oplot, stmass, bololum, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 10.0
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40e', size=1.2, color=thiscolor, /data, charthick=4.0


;---------------------------


fruns= ['d0h2_q', 'd1h2_q', 'd2h2_q', 'd3h7', 'd4h2_q', 'd5h2_q', 'd6h2_q'] & nfruns= n_elements(fruns)
stmass= [0.314, 0.871, 2.68, 7.37, 19.5, 57.8, 226.9] * 1.0d+10 / 0.7
bololum= stmass
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='ds', ibololum, powerratio
    bololum[i]= max(ibololum)
endfor
oplot, stmass, bololum, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 9.7
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '40h', size=1.2, color=thiscolor, /data, charthick=4.0


;---------------------------


fruns= ['e0e', 'e1e', 'e2e', 'e3e', 'e4e', 'e5e', 'e6e'] & nfruns= n_elements(fruns)
stmass= [0.321, 0.919, 2.69, 7.38, 19.1, 59.4, 228.4] * 1.0d+10 / 0.7
bololum= stmass
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='es', ibololum, powerratio
    bololum[i]= max(ibololum)
endfor
oplot, stmass, bololum, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 9.4
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '20e', size=1.2, color=thiscolor, /data, charthick=4.0


;---------------------------


fruns= ['z0e', 'z1e', 'z2e', 'z3e', 'z4e', 'z5e', 'z6e'] & nfruns= n_elements(fruns)
stmass= [0.328, 0.963, 2.85, 7.64, 21.3, 61.3, 234.0] * 1.0d+10 / 0.7
bololum= stmass
pointt= 6
select_thispoint, pointt, thispsym, thiscolor
for i= 0, nfruns-1 do begin
    get_luminfo, fruns[i], subdir='zs', ibololum, powerratio
    bololum[i]= max(ibololum)
endfor
oplot, stmass, bololum, psym=-thispsym, color=thiscolor, thick=3.0
ylbl= 9.1
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl, '5e', size=1.2, color=thiscolor, /data, charthick=4.0


;---------------------------


x= [xmin, xmax]
y= [12.0, 12.0]
oplot, x, y, psym=-3, linestyle= 2, thick= 3.0, color= 0




; done
; ------
device, /close




end




;================================================================================================





;-------------------------------------------------
;-------------------------------------------------
pro peaklum2, junk


if not keyword_set(junk) then begin
        print, " "
        print, " peaklum2, junk"
        print, " "
        print, " "
        return
endif

filename='peaklum.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

;---------------------------


;yaxistitle = "!6Log Peak Bolometric Luminosity (L!D!9n!6!N)"
yaxistitle = "!6Log Peak Starburst Luminosity (L!D!9n!6!N)"
ymax = 13.9
ymin = 9.5


;xaxistitle = "!6Stellar Mass (M!D!9n!6!N)"
;xmax = 7.0e+12
;xmin = 2.0e+9
xaxistitle = "!6Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10




;---------------------------

x0= 0.18
x1= 0.97
y0= 0.15
y1= 0.98


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /xlog, $
        ;/ylog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtickformat='exp_label', $
        ;xtickformat='(a1)', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------

xpt1= 3.0e+12
xpt2= 3.0e+12
xpt3= 7.0e+12
    


;---------------------------


totmass= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass

;---------------------------


fruns1= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
subdir1= 'bs'

fruns2= ['d0e2_q', 'd1e2_q', 'd2e2_q', 'd3e7', 'd4e2_q', 'd5e2_q', 'd6e2_q']
subdir2= 'ds'

process_and_plot_lums, totmass, fruns1, subdir1, fruns2, subdir2, 1, $
				'std - w/ BH', xpt1, xpt2, xpt3, 10.6

;---------------------------


fruns1= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
subdir1= 'sb10_mass'

fruns2= ['d0e_no', 'd1e_no', 'd2e_no', 'd3e_no', 'd4e_no', 'd5e_no', 'd6e_no']
subdir2= subdir1


process_and_plot_lums, totmass, fruns1, subdir1, fruns2, subdir2, 3, $
				'sb10 - w/o BH', xpt1, xpt2, xpt3, 10.3

;---------------------------


fruns1= ['b0e_no', 'b1e_no', 'b2e_no', 'b3e_no', 'b4e_no', 'b5e_no', 'b6e_no']
subdir1= 'sb8_mass'

fruns2= ['d0e_no', 'd1e_no', 'd2e_no', 'd3e_no', 'd4e_no', 'd5e_no', 'd6e_no'] 
subdir2= subdir1


process_and_plot_lums, totmass, fruns1, subdir1, fruns2, subdir2, 5, $
                                'sb8 - w/o BH', xpt1, xpt2, xpt3, 10.0




;---------------------------


x= [xmin, xmax]
y= [12.0, 12.0]
oplot, x, y, psym=-3, linestyle= 2, thick= 3.0, color= 0




; done
; ------
device, /close




end




;-------------------------------------------------
;-------------------------------------------------
pro peaktd, junk


if not keyword_set(junk) then begin
        print, " "
        print, " peaktd, junk"
        print, " "
        print, " "
        return
endif

filename='peaktd.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

;---------------------------


;yaxistitle = "!6Log Peak Bolometric Luminosity (L!D!9n!6!N)"
yaxistitle = "!7D!6T (BH - Sb Peak) (Myr)"
ymax = 300.0
ymin = -300.0


;xaxistitle = "!6Stellar Mass (M!D!9n!6!N)"
;xmax = 7.0e+12
;xmin = 2.0e+9
xaxistitle = "!6Mass (M!D!9n!6!N)"
xmax = 1.6e+14
xmin = 6.0e+10




;---------------------------

x0= 0.18
x1= 0.97
y0= 0.15
y1= 0.98


;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /xlog, $
        ;/ylog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtickformat='exp_label', $
        ;xtickformat='(a1)', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;---------------------------

xpt1= 3.0e+12
xpt2= 3.0e+12
xpt3= 7.0e+12
    


;---------------------------


totmass= [8.1668, 23.81, 70.727, 190.48, 495.174, 1523.84, 5812.99] * 1.0d+10 / 0.7      ; total mass

;---------------------------


fruns1= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
subdir1= 'bs'

fruns2= ['d0e2_q', 'd1e2_q', 'd2e2_q', 'd3e7', 'd4e2_q', 'd5e2_q', 'd6e2_q']
subdir2= 'ds'

process_and_plot_td, totmass, fruns1, subdir1, fruns2, subdir2, 1, $
				'std - w/ BH', xpt1, xpt2, xpt3, 10.6

;---------------------------


fruns1= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
subdir1= 'sb10_mass'

fruns2= ['d0e', 'd1e', 'd2e', 'd3e', 'd4e', 'd5e', 'd6e']
subdir2= subdir1


process_and_plot_td, totmass, fruns1, subdir1, fruns2, subdir2, 3, $
				'sb10 - w/o BH', xpt1, xpt2, xpt3, 10.3

;---------------------------


fruns1= ['b0e', 'b1e', 'b2e', 'b3e', 'b4e', 'b5e', 'b6e']
subdir1= 'sb8_mass'

fruns2= ['d0e', 'd1e', 'd2e', 'd3e', 'd4e', 'd5e', 'd6e'] 
subdir2= subdir1


process_and_plot_td, totmass, fruns1, subdir1, fruns2, subdir2, 5, $
                                'sb8 - w/o BH', xpt1, xpt2, xpt3, 10.0




;---------------------------


x= [xmin, xmax]
y= [0.0, 0.0]
oplot, x, y, psym=-3, linestyle= 2, thick= 3.0, color= 0




; done
; ------
device, /close




end


