;-------------------------------------------------------
;-------------------------------------------------------
;
;
;
;
;-------------------------------------------------------
;-------------------------------------------------------
pro toomre_multi, frun, filename=filename, $
			snaparr=snaparr, $
			xlen=xlen, $
			interaction=interaction, $
			msg=msg


if not keyword_set(frun) then begin
   print, "  "
   print, "toomre_multi, frun, filename=filename, "
   print, "    snaparr=snaparr, xlen=xlen "
   print, "  "
   print, "  "
   return
endif

; 2 x 3
npanels= 6
newxsize= 30
newysize= 20
x0= 0.008
xs= 0.328   ; assumes 3 panels
y0= 0.010
ys= 0.490   ; assumes 4 panels

; 4 x 3
;npanels= 12
;newxsize= 18
;newysize= 24
;x0= 0.008
;xs= 0.328   ; assumes 3 panels
;y0= 0.006
;ys= 0.247   ; assumes 4 panels


if not keyword_set(filename) then filename='panelfig.eps'
if not keyword_set(xlen) then xlen= 30.0
xlenarr= fltarr(npanels)
xlenarr(*)= xlen

if not keyword_set(snaparr) then begin
	print, "  "
	print, "  WARNING: snaparr not set - using /evensampling"
	print, "  "
	;return
	evensampling= 1
endif else begin
	evensampling= 0
endelse

if keyword_set(evensampling) then begin

	; determine the number of
	; snapshots in frun directory
	spawn, "/bin/ls "+frun+"/snapshot* | wc ",result
	nsnaps=long(result[0])

	snaparr= (indgen(12) + 1) * long(nsnaps/12)
	idx=where(snaparr ge nsnaps)
	if idx(0) ne -1 then snaparr(idx)= nsnaps-1
	snaparr[11]= nsnaps-1

	if not keyword_set(filename) then filename=frun+"/panel_gas_sample.eps"
	if filename eq "panelfig.eps" then filename=frun+"/panel_gas_sample.eps"

endif





;--------------------------------------
;  Now plot this mess (to postscript)
;--------------------------------------


        initialize_plotinfo, 1

	;
	;  need to match xloadct at line 528
	;
	;setup_plot_stuff, 'ps', filename=filename, colortable= 4
        ;setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=newxsize, newysize=newysize
        setup_plot_stuff, 'ps', filename=filename, colortable= 3, newxsize=newxsize, newysize=newysize
        ;setup_plot_stuff, 'ps', filename=filename, colortable= 1, newxsize=newxsize, newysize=newysize
        ;setup_plot_stuff, 'ps', filename=filename, colortable= 0, newxsize=newxsize, newysize=newysize


	x1= x0+xs
	x2= x0+xs+xs
	x3= x0+xs+xs+xs

        y1= y0+ys
	y2= y0+ys+ys
	y3= y0+ys+ys+ys
	y4= y0+ys+ys+ys+ys



	;  place the image(s) down
	; ----------------------

	if not keyword_set(interaction) then begin
		; good for the toomredisk_dv_pro*
		toomre_do_one_panel, frun, snaparr[0], xlenarr[0], x0, y1, xs, ys, /pt1center
		toomre_do_one_panel, frun, snaparr[1], xlenarr[1], x1, y1, xs, ys, /pt1center
		toomre_do_one_panel, frun, snaparr[2], xlenarr[2], x2, y1, xs, ys, /pt1center
		toomre_do_one_panel, frun, snaparr[3], xlenarr[3], x0, y0, xs, ys, /pt1center
		toomre_do_one_panel, frun, snaparr[4], xlenarr[4], x1, y0, xs, ys, /pt1center
		toomre_do_one_panel, frun, snaparr[5], xlenarr[5], x2, y0, xs, ys, /pt1center
	endif else begin
		; good for the toomrepro*
		;toomre_do_one_panel, frun, snaparr[0], xlenarr[0], x0, y1, xs, ys, /pt2center, /remove_traj2, /remove_disk1
		;toomre_do_one_panel, frun, snaparr[1], xlenarr[1], x1, y1, xs, ys, /pt2center, /remove_traj2, /remove_disk1
		;toomre_do_one_panel, frun, snaparr[2], xlenarr[2], x2, y1, xs, ys, /pt2center, /remove_traj2, /remove_disk1
		;toomre_do_one_panel, frun, snaparr[3], xlenarr[3], x0, y0, xs, ys, /pt2center, /remove_traj2, /remove_disk1
		;toomre_do_one_panel, frun, snaparr[4], xlenarr[4], x1, y0, xs, ys, /pt2center, /remove_traj2, /remove_disk1
		;toomre_do_one_panel, frun, snaparr[5], xlenarr[5], x2, y0, xs, ys, /pt2center, /remove_traj2, /remove_disk1
		toomre_do_one_panel, frun, snaparr[0], xlenarr[0], x0, y1, xs, ys, /pt2center, /remove_disk1
		toomre_do_one_panel, frun, snaparr[1], xlenarr[1], x1, y1, xs, ys, /pt2center, /remove_disk1
		toomre_do_one_panel, frun, snaparr[2], xlenarr[2], x2, y1, xs, ys, /pt2center, /remove_disk1
		toomre_do_one_panel, frun, snaparr[3], xlenarr[3], x0, y0, xs, ys, /pt2center, /remove_disk1
		toomre_do_one_panel, frun, snaparr[4], xlenarr[4], x1, y0, xs, ys, /pt2center, /remove_disk1
		toomre_do_one_panel, frun, snaparr[5], xlenarr[5], x2, y0, xs, ys, /pt2center, /remove_disk1
	endelse



; --------------------
;msg= ' '
if keyword_set(msg) then xyouts, x0+0.04, y1+0.05, msg, /normal, size= 1.8, charthick=3.0, color= 0   ; 0=black, 1=white


; put time units on panel
;xyouts, x0+0.09, y1-0.03, '!6Gyr/h', /normal, size= 1.2, charthick=3.0, color= 0
xyouts, x0+0.08, y2-0.06, '!6Gyr', /normal, size= 1.8, charthick=3.0, color= 0


; put side length on panel
xlenlbl= strcompress(string(2.0*xlen / 0.7),/remove_all)
if (2.0*xlen/0.7) ge 1.0 then digs= 1
if (2.0*xlen/0.7) ge 10.0 then digs= 2
if (2.0*xlen/0.7) ge 100.0 then digs= 3
xlenlbl = strmid(xlenlbl,0,digs)        ; T=0.x (4+digits after decimal)
xlenlbl= '!94!6'+xlenlbl+'!96!6'
xyouts, x0+0.025, y2-0.14, xlenlbl, /normal, size= 1.8, charthick=3.0, color= 0






        ; done, close this up
        ; --------------------
	device, /close






; -------------
;  Done
; -------------



end






pro doit, junk

;xlen= 28.0
xlen= 42.0
;xlen= 56.0

;snaparr= [0, 100, 225, 325, 425, 550]
snaparr= [0, 50, 100, 150, 215, 320]

toomre_multi, "data1/tides/toomredisk_dv_pro1", filename="tddvp1.eps", snaparr=snaparr, xlen=(xlen+14.0), msg="impulse"
toomre_multi, "data1/tides/toomredisk_dv_pro2", filename="tddvp2.eps", snaparr=snaparr, xlen=xlen, msg="impulse"
toomre_multi, "data1/tides/toomredisk_dv_pro3", filename="tddvp3.eps", snaparr=snaparr, xlen=xlen, msg="impulse"
toomre_multi, "data1/tides/toomredisk_dv_pro4", filename="tddvp4.eps", snaparr=snaparr, xlen=xlen, msg="impulse"


; at r_peri in snap=1271
toomre_multi, "data1/tides/toomrepro1", filename="tp1.eps", snaparr=1271+snaparr, xlen=xlen+14.0, /interaction, msg='3-body'
; at r_peri in snap=1169
toomre_multi, "data1/tides/toomrepro2", filename="tp2.eps", snaparr=1169+snaparr, xlen=xlen, /interaction, msg='3-body'
; at r_peri in snap=460
toomre_multi, "data1/tides/toomrepro3", filename="tp3.eps", snaparr=460+snaparr, xlen=xlen, /interaction, msg='3-body'
; at r_peri in snap=413
toomre_multi, "data1/tides/toomrepro4", filename="tp4.eps", snaparr=413+snaparr, xlen=xlen, /interaction, msg='3-body'



end




