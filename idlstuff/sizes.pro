;==================================================================
;==================================================================
;
;
;
;==================================================================
;==================================================================


pro re_time, junk


if not keyword_set(junk) then begin
	print, " "
	print, " re_time, junk"
	print, " "
	print, " "
	return
endif

filename='re.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4


yaxistitle = "!6R!De!N(kpc)"
xaxistitle = "!6Time (Gyr)"
;xmax= 6.0
;xmax= 4.4
xmax= 3.8
xmin = 0
;ymax = 15.0
ymax = 13.0
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


;do_one_re, "data1/remergers/b3eb3e2", pt= 0
;xyouts, 0.22, 0.94, "b3eb3e2", /normal, color=0, size=1.2

;do_one_re, "data1/remergers/b3ev2", pt= 200
;xyouts, 0.22, 0.78, "b3ev2", /normal, color=200, size=1.2

;do_one_re, "data1/remergers/b3ev3", pt= 100
;xyouts, 0.22, 0.82, "b3ev3", /normal, color=100, size=1.2

;do_one_re, "data1/remergers/b3ev4", pt= 150
;xyouts, 0.22, 0.86, "b3ev4", /normal, color=150, size=1.2

;do_one_re, "data1/remergers/b3esph", pt= 50
;xyouts, 0.22, 0.90, "b3esph", /normal, color=50, size=1.2


;---------------------------

; set the colortable to 1
; use xmax= 4.4 and ymax= 6.0

;xyouts, 0.52, 0.94, "ErrTolIntAccuracy", /normal, color= 0, size=1.2, charthick=2.0

;do_one_re, "data/bs/b3e_i2", pt= 240
;xyouts, 0.55, 0.90, '1.0', /normal, color=240, size=1.2, charthick=2.0

;do_one_re, "data/bs/b3e_i1", pt= 200
;xyouts, 0.55, 0.86, '0.1', /normal, color=200, size=1.2, charthick=2.0

;do_one_re, "data/bs/b3e", pt= 160
;xyouts, 0.55, 0.82, '0.025', /normal, color=160, size=1.2, charthick=2.0

;do_one_re, "data/bs/b3e_i4", pt= 120
;xyouts, 0.55, 0.78, '0.025 (diff Gad)', /normal, color=120, size=1.2, charthick=2.0

;do_one_re, "data/bs/b3e_i3", pt= 80
;xyouts, 0.55, 0.74, '0.0025', /normal, color=80, size=1.2, charthick=2.0

;do_one_re, "data/bs/b3e_i5", pt= 40
;xyouts, 0.55, 0.70, '0.00125', /normal, color=40, size=1.2, charthick=2.0

;do_one_re, "data/bs/b3e_i6", pt= 0
;xyouts, 0.55, 0.66, '0.00025', /normal, color=0, size=1.2, charthick=2.0



;---------------------------

; turn off printing of error range

xyouts, 0.52, 0.94, "Quiescent Galaxies", /normal, color= 0, size=1.2, charthick=2.0

do_one_re, "data1/remergers/iso_b3e", pt= 0
xyouts, 0.75, 0.24, 'b3e,', /normal, color=0, size=1.2, charthick=2.0

do_one_re, "data1/remergers/iso_b3e_2", pt= 0
xyouts, 0.82, 0.24, 'b3e_2', /normal, color=0, size=1.2, charthick=2.0

do_one_re, "data1/remergers/iso_sph", pt= 240
xyouts, 0.75, 0.80, 'sph', /normal, color=240, size=1.2, charthick=2.0

do_one_re, "data1/remergers/iso_v2v2", pt= 150
xyouts, 0.75, 0.32, 'v2v2', /normal, color=150, size=1.2, charthick=2.0

do_one_re, "data1/remergers/iso_v3x2", pt= 200
xyouts, 0.75, 0.43, 'v3x2', /normal, color=200, size=1.2, charthick=2.0

do_one_re, "data1/remergers/iso_v4x3", pt= 100
xyouts, 0.75, 0.58, 'v4x3', /normal, color=100, size=1.2, charthick=2.0

do_one_re, "data1/remergers/iso_v4x4", pt= 70
xyouts, 0.65, 0.66, 'v4x4,', /normal, color=70, size=1.2, charthick=2.0

do_one_re, "data1/remergers/iso_v4x5", pt= 30
xyouts, 0.83, 0.66, 'v4x5', /normal, color=30, size=1.2, charthick=2.0

do_one_re, "data1/remergers/iso_v4x6", pt= 50
xyouts, 0.74, 0.66, 'v4x6,', /normal, color=50, size=1.2, charthick=2.0



;---------------------------


; done
; ------
device, /close


end






;==================================================================
;==================================================================



pro do_one_re, frun, pt=pt

spawn, "/bin/ls "+frun+"/R_e_*.txt",result

for i=0, n_elements(result)-1 do begin

	Refile= result[i]
	read_file_r_e, Refile, snapext, time, re, reerr

	time= time/0.7
	re= re/0.7
	reerr= reerr/0.7

	;select_thispoint, pt, thispsym, thiscolor
	thispsym= 3
	thiscolor= pt

	oplot, time, re, thick=4.0, psym=-thispsym, color=thiscolor
	;oplot, time, re+reerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
	;oplot, time, re-reerr, thick=1.0, psym=-3, color=thiscolor, linestyle=1
endfor


end



;==================================================================
;==================================================================




pro size_mstar, junk


if not keyword_set(junk) then begin
	print, " "
	print, " re_time, junk"
	print, " "
	print, " "
	return
endif

filename='remstar.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


yaxistitle = "!6R!Deff!N(kpc)"
xaxistitle = "!6Stellar Mass (M!D!9n!6!N)"
ymax= 14.5
ymin= 0.27
xmax = 2.1e+12
xmin = 0.8e+10


;---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
	color= 0, $
	xstyle=1, $
	ystyle=1, $
	/ylog, /xlog, $
	;ystyle=8, $     ; this will suppress the second y axis
	xcharsize=1.50, ycharsize=1.50, $
	xthick=4.0, ythick=4.0, $
	charthick=3.0, $
	;ytickformat='exp_label', $
	xtickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata

;---------------------------

; Add NYU Value Added Catalog

;---------------------------

; Add Bezanson et al. local E galaxies

;---------------------------

; hand-made z=0 spheroid size relation

fload_newcolortable, 0
x= [xmin, xmin, 3.0e+11, xmax, xmax, xmin]
y= [0.5, 3.0, ymax, ymax, 9.0, 0.5]
polyfill, x, y, /data, color= 200, /fill
fload_newcolortable, 4

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, xstyle=1, ystyle=1, /ylog, /xlog, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, $
        charthick=3.0, xtickformat='exp_label', $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase


;---------------------------

; Add the z~2.5 compact galaxies from
; van Dokkum et al. 2008

sizes=  [0.47, 0.49, 0.76, 0.78, 0.92,  0.93, 1.42, 1.89, 2.38]     ; in kpc
mstars= [6.0,  8.8, 32.02, 23.7, 10.3, 20.01, 17.5, 15.0, 20.7]    ; in 10^10 Msolar
mstars= mstars * 1.0d+10        ; in m_solar

symsize=1.2
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
oplot, mstars, sizes, psym=8, color= 0

;---------------------------

; Add Shen et al. 2003 scaling

x= [xmin,xmin*2.,xmax/2.,xmax]
y= 2.88e-6 * x^(0.56)
oplot, x, y, psym=-3, color= 0, linestyle=1

;---------------------------

; just do this by hand
;grab_re_mstar, "data1/remergers/iso_b3e", 10, re=re, mstar=mstar, pt=150
re_0= 0.62 / 0.7     ; in kpc (b3e remnant as determined in the 001 snapshot for each merger)
mstar_0= 1.0e+10 * (1.56 + 4.826) / 0.7   ; in m_solar

;---------------------------

grab_and_plot_re_mstar, "data1/remergers/b3eb3e2", re_0, mstar_0, pt= 150

grab_and_plot_re_mstar, "data1/remergers/b3esph", re_0, mstar_0, pt= 200

grab_and_plot_re_mstar, "data1/remergers/b3ev2", re_0, mstar_0, pt= 25

grab_and_plot_re_mstar, "data1/remergers/b3ev3", re_0, mstar_0, pt= 50

grab_and_plot_re_mstar, "data1/remergers/b3ev4", re_0, mstar_0, pt= 75


;---------------------------

symsize=2.0
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
oplot, [mstar_0], [re_0], thick=4.0, psym=8, color=150

;---------------------------


; done
; ------
device, /close


end




;==================================================================
;==================================================================


pro grab_and_plot_re_mstar, frun, re_0, mstar_0, pt=pt

nsnaps= fload_frun_nsnaps(frun)

ok= fload_snapshot_bh(frun,nsnaps-1,/nopot_in_snap,/skip_center)
bhids= fload_blackhole_id(1,n_bh=n_bh)

if n_bh gt 1 then begin
	print, " "
	print, " run not merged yet: ", frun
	print, " "
	return
endif

mstar= 1.0e+10 * total(fload_allstars_mass(1)) / 0.7
ttime= fload_time(1)-1.0e-3

Refile= frun+"/R_e_"+strcompress(string(bhids),/remove_all)+".txt"
read_file_r_e, Refile, snapext, time, ret, reerr

idx=where(time ge ttime)
re= ret[idx(0)]/0.7


thiscolor= pt
oplot, [mstar], [re], thick=4.0, psym=4, color=thiscolor
oplot, [mstar_0, mstar], [re_0, re], thick=4.0, psym=-3, color=thiscolor

end





;==================================================================
;==================================================================


