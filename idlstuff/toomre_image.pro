;-------------------------------------------------------
;-------------------------------------------------------
;
;
;
;
;-------------------------------------------------------
;-------------------------------------------------------
pro toomre_image, frun, snapnum, filename=filename, xlen=xlen, $
			pt1center=pt1center, $
			pt2center=pt2center, $
        		remove_disk1=remove_disk1, $
        		remove_disk2=remove_disk2, $
        		remove_pt1=remove_pt1, $
        		remove_pt2=remove_pt2, $
			remove_traj1=remove_traj1, $
			remove_traj2=remove_traj2


if not keyword_set(frun) then begin
   print, "  "
   print, "toomre_image, frun, filename=filename, xlen=xlen "
   print, "  "
   print, "  "
   return
endif


if not keyword_set(filename) then filename='toomreimage.eps'
if not keyword_set(snapnum) then snapnum=0
if not keyword_set(xlen) then xlen= 70.0

center=[0.0,0.0,0.0]



;--------------------------------------
;  Now plot this mess (to postscript)
;--------------------------------------


initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize= 25.0, newysize= 25.0
;setup_plot_stuff, 'ps', filename=filename, colortable= 3
;setup_plot_stuff, 'ps', filename=filename, colortable= 1
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


x0= 0.01
x1= 0.99

y0= 0.01
y1= 0.99


print, x0, y0, x1, y1

toomre_do_one_panel, frun, snapnum, xlen, x0, y0, (x1-x0), (y1-y0), center=center, $
			pt1center=pt1center, $
			pt2center=pt2center, $
        		remove_disk1=remove_disk1, $
        		remove_disk2=remove_disk2, $
        		remove_pt1=remove_pt1, $
        		remove_pt2=remove_pt2, $
			remove_traj1=remove_traj1, $
			remove_traj2=remove_traj2




; --------------------
msg= ' '
xyouts, 0.02, 0.94, msg, /normal, size= 1.2, charthick=3.0, color= 0   ; 0=black, 1=white

print, x0, y0, x1, y1

; put time units on panel
;xyouts, x0+0.09, y1-0.03, '!6Gyr/h', /normal, size= 1.2, charthick=3.0, color= 0
xyouts, x0+0.14, y1-0.07, '!6Gyr', /normal, size= 1.8, charthick=3.0, color= 0


; put side length on panel
xlenlbl= strcompress(string(2.0*xlen / 0.7),/remove_all)
if (2.0*xlen) ge 1.0 then digs= 1
if (2.0*xlen) ge 10.0 then digs= 2
if (2.0*xlen) ge 100.0 then digs= 3
xlenlbl = strmid(xlenlbl,0,digs)        ; T=0.x (4+digits after decimal)
xlenlbl= '!94!6'+xlenlbl+'!96!6'
xyouts, x1-0.14, y0+0.05, xlenlbl, /normal, size= 1.8, charthick=3.0, color= 0


; done, close this up
; --------------------
device, /close


end






