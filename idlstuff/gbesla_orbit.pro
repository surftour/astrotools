
pro doit, junk


  show_centers_time, junk
  velocity_info, junk

end



; ------------------------------------------------------------------
; ------------------------------------------------------------------
; ------------------------------------------------------------------
; ------------------------------------------------------------------
; ------------------------------------------------------------------
; ------------------------------------------------------------------


pro show_centers_time, junk


if not keyword_set(junk) then begin
	print, "  "
	print, "  show_centers_time, junk"
	print, "  "
endif


filename='center.eps'

timemin= 0.0
;timemin= 1.0
;timemax= 2.0
;timemax= 4.0
;timemax= 4.25
;timemax= 4.6
;timemax= 6.0
timemax= 21.5
;timemax= 25.0

;ymax= 1000.0
;ymax= 300.0
ymax= 200.0
;ymax= 175.0
;ymax= 120.0
;ymax= 100.0
;ymax= 1.4
;ymin= 0.1
ymin= 0


; ----------------------------------
;   Try this plot thingy
; ----------------------------------

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename=filename, colortable= 4



!p.position= [0.19, 0.15, 0.97, 0.98]

plot, [1.0], [1.0], psym=-3, $
        xrange=[timemin,timemax], $
	yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
	;ytickformat='exp_label', /ylog, $
        color= 0, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle="!6Time (Gyr)", $
        ;ytitle="Center Separation (kpc)", $
        ytitle="Center Separation (!8h!6!E-1!Nkpc)", $
        ;ytitle="MW/M31 Separation (kpc)", $
        ;ytitle="!6MW/M31 Separation (Mpc)", $
        /nodata



; go through files
; ------------------

x0= 0.22
y0= 0.98


; -------------------------------

;process_onesim_centers, "data1/lg_m33/Sbm33_1", 1, x0, y0-0.03, /m31m33
;process_onesim_centers, "data1/lg_m33/Sbm33_2", 2, x0, y0-0.06, /m31m33
process_onesim_centers, "data1/lg_m33/Sbm33_3", 3, x0, y0-0.09, /m31m33

;process_onesim_centers, "data1/lg_m33/SbSbm33_1", 6, x0, y0-0.15, /mwm31m33
;process_onesim_centers, "data1/lg_m33/SbSbm33_2", 6, x0, y0-0.18, /mwm31m33
;process_onesim_centers, "data1/lg_m33/SbSbm33_3", 6, x0, y0-0.21, /mwm31m33
;process_onesim_centers, "data1/lg_m33/SbSbm33_4", 6, x0, y0-0.24, /mwm31m33
;process_onesim_centers, "data1/lg_m33/SbSbm33_5", 6, x0, y0-0.27, /mwm31m33
;process_onesim_centers, "data1/lg_m33/SbSbm33_6", 6, x0, y0-0.30, /mwm31m33
;process_onesim_centers, "data1/lg_m33/SbSbm33_7", 6, x0, y0-0.33, /mwm31m33

; -------------------------------


xyouts, 0.76, 0.90, 'centers', /normal, size= 1.2, color= 0, charthick= 6.0
xyouts, 0.76, 0.85, 'COM', /normal, size= 1.2, color= 0, charthick= 1.0


; -------------------------------


device, /close


end








; do the dirty work
; --------------------
pro process_onesim_centers, frun, pointselection, x0, y0, ylog=ylog, msg=msg, $
						m31m33=m31m33, mwm31m33=mwm31m33

;--------------------------------------
print, "processing: ", frun

; read centers
; -------------
if keyword_set(mwm31m33) then begin
	read_center, cp_time, cp_cen1, filename=frun+"/centers_den_2600002.txt"        ; m31
	read_center, cp_time, cp_cen2, filename=frun+"/centers_den_2760003.txt"        ; m33
endif

if keyword_set(m31m33) then begin
	;read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=frun+"/centers.txt"
	read_center, cp_time, cp_cen1, filename=frun+"/centers_den_1300001.txt"        ; m31
	read_center, cp_time, cp_cen2, filename=frun+"/centers_den_1460002.txt"        ; m33
	;read_center, cp_time, cp_cen1, filename=frun+"/centers_bh_1300001.txt"        ; m31
	;read_center, cp_time, cp_cen2, filename=frun+"/centers_bh_1460002.txt"        ; m33
endif

cp_time= transpose(cp_time)
cp_cen1= transpose(cp_cen1)
cp_cen2= transpose(cp_cen2)


x = cp_cen1[*,0] - cp_cen2[*,0]
y = cp_cen1[*,1] - cp_cen2[*,1]
z = cp_cen1[*,2] - cp_cen2[*,2]

rdiff = sqrt(x*x + y*y + z*z)

;rdiff= rdiff/0.7
;rdiff= rdiff/1000.0
cp_time= cp_time/0.7

;--------------------------------------

if keyword_set(ylog) then begin
idx= where(rdiff lt 1.0e-6)
if idx(0) ne -1 then rdiff(idx)= 1.0e-6
endif

select_thispoint, pointselection, thispsym, thiscolor

;oplot, cp_time, rdiff, psym=-thispsym, color= thiscolor, thick=1.0
;oplot, cp_time, rdiff, psym=-3, color= thiscolor, thick=6.0
oplot, cp_time, rdiff, psym=-3, color= 0, thick=6.0

if keyword_set(msg) then begin
	xyouts, x0, y0, msg, /normal, size= 1.1, color= thiscolor, charthick= 3.0
endif else begin
	;xyouts, x0, y0, fload_getid(frun), /normal, size= 1.1, color= thiscolor, charthick= 3.0
	xyouts, x0, y0, fload_getid(frun), /normal, size= 1.1, color= thiscolor, charthick= 3.0
endelse

;=======================================================


; read centers
; -------------
if keyword_set(mwm31m33) then begin
        read_center, cp_time, cp_cen1, filename=frun+"/com_2600002.txt"        ; m31
        read_center, cp_time, cp_cen2, filename=frun+"/com_2760003.txt"        ; m33
endif

if keyword_set(m31m33) then begin
        ;read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=frun+"/centers.txt"
        read_center, cp_time, cp_cen1, filename=frun+"/com_1300001.txt"        ; m31
        read_center, cp_time, cp_cen2, filename=frun+"/com_1460002.txt"        ; m33
        ;read_center, cp_time, cp_cen1, filename=frun+"/centers_bh_1300001.txt"        ; m31
        ;read_center, cp_time, cp_cen2, filename=frun+"/centers_bh_1460002.txt"        ; m33
endif

cp_time= transpose(cp_time)
cp_cen1= transpose(cp_cen1)
cp_cen2= transpose(cp_cen2)


x = cp_cen1[*,0] - cp_cen2[*,0]
y = cp_cen1[*,1] - cp_cen2[*,1]
z = cp_cen1[*,2] - cp_cen2[*,2]

rdiff = sqrt(x*x + y*y + z*z)

;rdiff= rdiff/0.7
;rdiff= rdiff/1000.0
cp_time= cp_time/0.7

;--------------------------------------

if keyword_set(ylog) then begin
idx= where(rdiff lt 1.0e-6)
if idx(0) ne -1 then rdiff(idx)= 1.0e-6
endif

select_thispoint, pointselection, thispsym, thiscolor

;oplot, cp_time, rdiff, psym=-thispsym, color= thiscolor, thick=1.0
;oplot, cp_time, rdiff, psym=-3, color= thiscolor, thick=1.0, linestyle= 2
oplot, cp_time, rdiff, psym=-3, color= 0, thick=1.0, linestyle= 2


end





; ------------------------------------------------------------------
; ------------------------------------------------------------------
; ------------------------------------------------------------------
; ------------------------------------------------------------------
; ------------------------------------------------------------------





pro velocity_info, frun

	if not keyword_set(frun) then begin
		print, " "
		print, " velocity_info, frun  (need centers.txt and cenvel.txt) "
		print, " "
		return
	endif



; ----------------------------------
;   Try this plot thingy
; ----------------------------------

;timemin= 0.0
;timemax= max(center_time)
timemin= 1.0
;timemax= 2.0
timemax= 21.5

;ymax= 500.0
ymax= 250.0
;ymax= 120.0

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename='velocity_info.eps', colortable= 4



!p.position= [0.18, 0.15, 0.98, 0.99]

plot, [1.0], [1.0], psym=-3, $
        xrange=[timemin,timemax], $
	yrange=[0,ymax], $
        xstyle=1, ystyle=1, $
        color= 0, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        ;xtitle="!6Time (h!E-1!N Gyr)", $
        xtitle="!6Time (Gyr)", $
        ytitle="!6Relative Velocity (km sec!E-1!N)", $
        /nodata


; -------------------------------

x0= 0.22
y0= 0.98


;process_onesim_velocity, "data1/lg_m33/Sbm33_1", 1, x0, y0-0.03, /m31m33
;process_onesim_velocity, "data1/lg_m33/Sbm33_2", 2, x0, y0-0.06, /m31m33
process_onesim_velocity, "data1/lg_m33/Sbm33_3", 3, x0, y0-0.09, /m31m33

;process_onesim_velocity, "data1/lg_m33/SbSbm33_1", 5, x0, y0-0.03, /mwm31m33
;process_onesim_velocity, "data1/lg_m33/SbSbm33_2", 5, x0, y0-0.06, /mwm31m33
;process_onesim_velocity, "data1/lg_m33/SbSbm33_3", 5, x0, y0-0.09, /mwm31m33
;process_onesim_velocity, "data1/lg_m33/SbSbm33_4", 5, x0, y0-0.12, /mwm31m33
;process_onesim_velocity, "data1/lg_m33/SbSbm33_5", 5, x0, y0-0.15, /mwm31m33
;process_onesim_velocity, "data1/lg_m33/SbSbm33_6", 5, x0, y0-0.18, /mwm31m33
;process_onesim_velocity, "data1/lg_m33/SbSbm33_7", 5, x0, y0-0.21, /mwm31m33


; -------------------------------

xyouts, 0.76, 0.90, 'centers', /normal, size= 1.2, color= 0, charthick= 6.0
xyouts, 0.76, 0.85, 'COM', /normal, size= 1.2, color= 0, charthick= 1.0

;xyouts, 0.62, 0.93, 'total relative', /normal, size= 1.2, color= 0, charthick= 3.0
;xyouts, 0.65, 0.89, 'velocity', /normal, size= 1.2, color= 0, charthick= 3.0
;xyouts, 0.62, 0.84, 'radial', /normal, size= 1.2, color= 150, charthick= 3.0
;xyouts, 0.62, 0.79, 'tangential', /normal, size= 1.2, color= 50, charthick= 3.0


;x= [5.4,5.4]
;y=[0,300.0]
;oplot, x, y, psym=-3, color= 0, thick=3.0, linestyle= 1

; -------------------------------

device, /close


end



;===============================================================





; do the dirty work
; --------------------
pro process_onesim_velocity, frun, pointselection, x0, y0, ylog=ylog, msg=msg, $
						m31m33=m31m33, mwm31m33=mwm31m33

;--------------------------------------

; read centers
; -------------
if keyword_set(mwm31m33) then begin
	read_center, center_time, center_1, filename=frun+"/centers_den_2600002.txt"        ; m31
	read_center, center_time, center_2, filename=frun+"/centers_den_2760003.txt"        ; m33
endif

if keyword_set(m31m33) then begin
	read_center, center_time, center_1, filename=frun+"/centers_den_1300001.txt"        ; m31
	read_center, center_time, center_2, filename=frun+"/centers_den_1460002.txt"        ; m33
	;read_center, center_time, center_1, filename=frun+"/centers_bh_1300001.txt"        ; m31
	;read_center, center_time, center_2, filename=frun+"/centers_bh_1460002.txt"        ; m33
endif

        center_time= transpose(center_time)
        center_1= transpose(center_1)
        center_2= transpose(center_2)

        ;center_time= center_time/0.7
        center_1= center_1/0.7
        center_2= center_2/0.7

        x = center_2[*,0] - center_1[*,0]
        y = center_2[*,1] - center_1[*,1]
        z = center_2[*,2] - center_1[*,2]

        rdiff = sqrt(x*x + y*y + z*z)



; read com velocity
; ------------------
if keyword_set(mwm31m33) then begin
	read_center, cenvel_time, cenvel_1, filename=frun+"/cenvel_den_2600002.txt"        ; m31
	read_center, center_time, cenvel_2, filename=frun+"/cenvel_den_2760003.txt"        ; m33
endif

if keyword_set(m31m33) then begin
	read_center, cenvel_time, cenvel_1, filename=frun+"/cenvel_den_1300001.txt"        ; m31
	read_center, center_time, cenvel_2, filename=frun+"/cenvel_den_1460002.txt"        ; m33
	;read_center, center_time, cenvel_1, filename=frun+"/cenvel_bh_1300001.txt"        ; m31
	;read_center, center_time, cenvel_2, filename=frun+"/cenvel_bh_1460002.txt"        ; m33
endif

        cenvel_time= transpose(cenvel_time)
        cenvel_1= transpose(cenvel_1)
        cenvel_2= transpose(cenvel_2)


        vx= cenvel_2[*,0] - cenvel_1[*,0]
        vy= cenvel_2[*,1] - cenvel_1[*,1]
        vz= cenvel_2[*,2] - cenvel_1[*,2]

        velrel= sqrt(vx*vx + vy*vy + vz*vz)

        vrad_x= vx*x/rdiff
        vrad_y= vy*y/rdiff
        vrad_z= vz*z/rdiff

        vrad= sqrt(vrad_x*vrad_x + vrad_y*vrad_y + vrad_z*vrad_z)
        vtan = sqrt(velrel*velrel - vrad*vrad)

        ;for i=0,n_elements(center_time)-1 do begin
        ;       print, center_time[i], rdiff[i], velrel[i], vrad[i], vtan[i]
        ;endfor


;--------------------------------------

center_time= center_time/0.7

if keyword_set(ylog) then begin
idx= where(rdiff lt 1.0e-6)
if idx(0) ne -1 then rdiff(idx)= 1.0e-6
endif

select_thispoint, pointselection, thispsym, thiscolor

oplot, center_time, velrel, psym=-3, color= 0, thick=6.0
;oplot, center_time, vrad, psym=-3, color= 150, thick=5.0, linestyle= 1
;oplot, center_time, vtan, psym=-3, color= 50, thick=5.0, linestyle= 2


if keyword_set(msg) then begin
	xyouts, x0, y0, msg, /normal, size= 1.1, color= thiscolor, charthick= 3.0
endif else begin
	;xyouts, x0, y0, fload_getid(frun), /normal, size= 1.1, color= thiscolor, charthick= 3.0
	xyouts, x0, 0.22, fload_getid(frun), /normal, size= 1.1, color= thiscolor, charthick= 3.0
endelse



;====================================

if keyword_set(m31m33) then begin 
        read_center, cenvel_time, cenvel_1, filename=frun+"/comvel_den_1300001.txt"        ; m31
        read_center, center_time, cenvel_2, filename=frun+"/comvel_den_1460002.txt"        ; m33
        ;read_center, center_time, cenvel_1, filename=frun+"/cenvel_bh_1300001.txt"        ; m31
        ;read_center, center_time, cenvel_2, filename=frun+"/cenvel_bh_1460002.txt"        ; m33
endif

        cenvel_time= transpose(cenvel_time)
        cenvel_1= transpose(cenvel_1)
        cenvel_2= transpose(cenvel_2)


        vx= cenvel_2[*,0] - cenvel_1[*,0]
        vy= cenvel_2[*,1] - cenvel_1[*,1]
        vz= cenvel_2[*,2] - cenvel_1[*,2]

        velrel= sqrt(vx*vx + vy*vy + vz*vz)

        vrad_x= vx*x/rdiff
        vrad_y= vy*y/rdiff
        vrad_z= vz*z/rdiff

        vrad= sqrt(vrad_x*vrad_x + vrad_y*vrad_y + vrad_z*vrad_z)
        vtan = sqrt(velrel*velrel - vrad*vrad)

        ;for i=0,n_elements(center_time)-1 do begin
        ;       print, center_time[i], rdiff[i], velrel[i], vrad[i], vtan[i]
        ;endfor


;--------------------------------------

center_time= center_time/0.7


if keyword_set(ylog) then begin
idx= where(rdiff lt 1.0e-6)
if idx(0) ne -1 then rdiff(idx)= 1.0e-6
endif

select_thispoint, pointselection, thispsym, thiscolor

oplot, center_time, velrel, psym=-3, color= 0, thick=1.0, linestyle= 2
;oplot, center_time, vrad, psym=-3, color= 150, thick=5.0, linestyle= 1
;oplot, center_time, vtan, psym=-3, color= 50, thick=5.0, linestyle= 2



;====================================





end
