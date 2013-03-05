;------------------------------------------------------------------------
;
;    Random Procedures related to Angular Momentum
;
;
;
;
;------------------------------------------------------------------------









; ---------------------------------
;
;  Radial dependence of J_tot
;
; ---------------------------------
pro radial_j, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "radial_j, junk"
   print, "  "
   return
endif

filename='jradial.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12, newysize=13




; ----------------------------- 
; set up constants      
; -----------------------------


;frun="/raid4/tcox/vc3vc3f"
frun="/raid4/tcox/cvc3vc3f"
snapnum=30
          
bins = 30
          
xmax = 50
xmin = 0.1

;do_z= 1
do_z= 0

xaxistitle= 'radius (kpc)'
yaxistitle= 'Normalized Specific J!Dtotal!N'
if do_z eq 1 then yaxistitle= 'Normalized Specific J!Dz!N'
ymax = 1.2
ymin = -1.2

x0= 0.18
x1= 0.98
y0= 0.15
y1= 0.98

!p.position=[x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, $ 
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase



; ----------------------
;  Now Do the Work
; ----------------------

ok=fload_snapshot_bh(frun,snapnum)




;  New Stars
; -----------

r_ns= fload_newstars_xyz('r')

;specific j's
ns_jx= fload_newstars_j(21)
ns_jy= fload_newstars_j(22)
ns_jz= fload_newstars_j(23)
;print, "J_tot= ", total(ns_jx), total(ns_jy), total(ns_jz)


;--------------------
; process the cumulative profile
if fload_npart(4) gt 0 then begin
	rs = fltarr(bins)
	j_tot= fltarr(bins)
	j_z= fltarr(bins)

	tempxmin= alog10(xmin)
	tempxmax= alog10(xmax)

	binsize = float((tempxmax-tempxmin))/bins

	for i=0, bins-1 do begin
	    sm_r_log= tempxmin + i*binsize
	    lg_r_log= tempxmin + (i+1)*binsize
	    sm_r= 10^(sm_r_log)
	    lg_r= 10^(lg_r_log)
	    rs[i]= 10^(sm_r_log + 0.5*binsize)

	    idx= where((r_ns gt sm_r) and (r_ns le lg_r))
	    if idx(0) ne -1 then begin
		j_x= total(ns_jx(idx))
		j_y= total(ns_jy(idx))
		j_z[i]= total(ns_jz(idx))
		j_tot[i]= sqrt(j_x*j_x + j_y*j_y + j_z[i]*j_z[i])
	    endif
	endfor

	j_z= j_z/1.0e+7
	j_tot= j_tot/1.0e+7

	if do_z eq 1 then begin
		oplot, rs, j_z, psym=-2, linestyle=0, color= 150, thick=4.0
	endif else begin
		oplot, rs, j_tot, psym=-2, linestyle=0, color= 150, thick=4.0
	endelse
endif


;  Old Stars
; -----------

r_old= fload_oldstars_xyz('r')

;specific j's
old_jx= fload_oldstars_j(21)
old_jy= fload_oldstars_j(22)
old_jz= fload_oldstars_j(23)
;print, "J_tot= ", total(ns2_jx), total(ns2_jy), total(ns2_jz)

;--------------------
; process the cumulative profile

rs = fltarr(bins)
j_tot= fltarr(bins)
j_z= fltarr(bins)

tempxmin= alog10(xmin)
tempxmax= alog10(xmax)

binsize = float((tempxmax-tempxmin))/bins

for i=0, bins-1 do begin
    sm_r_log= tempxmin + i*binsize
    lg_r_log= tempxmin + (i+1)*binsize
    sm_r= 10^(sm_r_log)
    lg_r= 10^(lg_r_log)
    rs[i]= 10^(sm_r_log + 0.5*binsize)

    idx= where((r_old gt sm_r) and (r_old le lg_r))
    if idx(0) ne -1 then begin
	j_x= total(old_jx(idx))
	j_y= total(old_jy(idx))
	j_z[i]= total(old_jz(idx))
	j_tot[i]= sqrt(j_x*j_x + j_y*j_y + j_z[i]*j_z[i])
    endif
endfor

j_z= j_z/1.0e+7
j_tot= j_tot/1.0e+7

if do_z eq 1 then begin
	oplot, rs, j_z, psym=-5, linestyle=0, color= 50, thick=4.0
endif else begin
	oplot, rs, j_tot, psym=-5, linestyle=0, color= 50, thick=4.0
endelse




; --------------------

; Normalize this
;normalizeJ= max([sage1,sage2])
;sage1=sage1/normalizeJ
;sage2=sage2/normalizeJ




; -----------------
;  Plot Extras
; -----------------


xyouts, 0.3, 0.85, 'New Stars', /normal, charthick=3.5, size= 1.4, color= 150
xyouts, 0.3, 0.80, 'Old Stars', /normal, charthick=3.5, size= 1.4, color= 50


;--------------------------------------
;--------------------------------------

device, /close





end





; ====================================================================================






; ---------------------------------
;
;  Radial dependence of J_z
;
; ---------------------------------
pro radial_total_j, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "radial_total_j, junk"
   print, "  "
   return
endif

filename='jradialtot.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12, newysize=13




; ----------------------------- 
; set up constants      
; -----------------------------


;frun="/raid4/tcox/vc3vc3f"
;frun="/raid4/tcox/cvc3vc3f"
;frun="/raid4/tcox/vc3vc3h"
;frun="/raid4/tcox/cvc3vc3h"
frun="/raid4/tcox/vc3vc3e"
;frun="/raid4/tcox/cvc3vc3h"
snapnum=30
          
bins = 30
          
xmax = 50
xmin = 0.1


xaxistitle= 'radius (kpc)'
yaxistitle= 'Normalized Specific J!Dz!N'
ymax = 1.2
ymin = -1.2

x0= 0.18
x1= 0.98
y0= 0.15
y1= 0.98

!p.position=[x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, $ 
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase



; ----------------------
;  Now Do the Work
; ----------------------

ok=fload_snapshot_bh(frun,snapnum)




;  New Stars
; -----------

r_ns= fload_newstars_xyz('r')

;specific j's
ns_jx= fload_newstars_j(21)
ns_jy= fload_newstars_j(22)
ns_jz= fload_newstars_j(23)
;print, "J_tot= ", total(ns_jx), total(ns_jy), total(ns_jz)


;--------------------
; process the cumulative profile
if fload_npart(4) gt 0 then begin
	rs = fltarr(bins)
	j_z_plus= fltarr(bins)
	j_z_min= fltarr(bins)

	tempxmin= alog10(xmin)
	tempxmax= alog10(xmax)

	binsize = float((tempxmax-tempxmin))/bins

	for i=0, bins-1 do begin
	    ; radius
	    sm_r_log= tempxmin + i*binsize
	    lg_r_log= tempxmin + (i+1)*binsize
	    sm_r= 10^(sm_r_log)
	    lg_r= 10^(lg_r_log)
	    rs[i]= 10^(sm_r_log + 0.5*binsize)

	    ; j
	    idx= where((r_ns gt sm_r) and (r_ns le lg_r))
	    if idx(0) ne -1 then begin
		thisjz= ns_jz(idx)
		idx= where(thisjz gt 0)
		if idx(0) ne -1 then j_z_plus[i]= total(thisjz(idx))
		idx= where(thisjz le 0)
		if idx(0) ne -1 then j_z_min[i]= total(thisjz(idx))
	    endif
	endfor

	j_z_plus= j_z_plus/1.0e+7
	j_z_min= j_z_min/1.0e+7

	oplot, rs, j_z_plus, psym=-2, linestyle=0, color= 150, thick=4.0
	oplot, rs, j_z_min, psym=-2, linestyle=0, color= 150, thick=4.0
endif


;  Old Stars
; -----------

r_old= fload_oldstars_xyz('r')

;specific j's
old_jx= fload_oldstars_j(21)
old_jy= fload_oldstars_j(22)
old_jz= fload_oldstars_j(23)
;print, "J_tot= ", total(ns2_jx), total(ns2_jy), total(ns2_jz)

;--------------------
; process the cumulative profile

rs = fltarr(bins)
j_z_plus= fltarr(bins)
j_z_min= fltarr(bins)

tempxmin= alog10(xmin)
tempxmax= alog10(xmax)

binsize = float((tempxmax-tempxmin))/bins

for i=0, bins-1 do begin
    sm_r_log= tempxmin + i*binsize
    lg_r_log= tempxmin + (i+1)*binsize
    sm_r= 10^(sm_r_log)
    lg_r= 10^(lg_r_log)
    rs[i]= 10^(sm_r_log + 0.5*binsize)

    idx= where((r_old gt sm_r) and (r_old le lg_r))
    if idx(0) ne -1 then begin
	thisjz= old_jz(idx)
	idx= where(thisjz gt 0)
	if idx(0) ne -1 then j_z_plus[i]= total(thisjz(idx))
	idx= where(thisjz le 0)
	if idx(0) ne -1 then j_z_min[i]= total(thisjz(idx))
    endif
endfor

j_z_plus= j_z_plus/1.0e+7
j_z_min= j_z_min/1.0e+7


oplot, rs, j_z_plus, psym=-5, linestyle=0, color= 50, thick=4.0
oplot, rs, j_z_min, psym=-5, linestyle=0, color= 50, thick=4.0




; --------------------

; Normalize this
;normalizeJ= max([sage1,sage2])
;sage1=sage1/normalizeJ
;sage2=sage2/normalizeJ




; -----------------
;  Plot Extras
; -----------------


xyouts, 0.3, 0.85, 'New Stars', /normal, charthick=3.5, size= 1.4, color= 150
xyouts, 0.3, 0.80, 'Old Stars', /normal, charthick=3.5, size= 1.4, color= 50


;--------------------------------------
;--------------------------------------

device, /close





end




;===================================================================================
;===================================================================================
;===================================================================================
;===================================================================================







pro stellar_age_hist, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "stellar_age_hist, junk"
   print, "  "
   return
endif

filename='agehist.eps' 

;frun1="/raid4/tcox/vc3vc3f"
;frun2="/raid4/tcox/vc3vc3f_no"
;frun1="/raid4/tcox/vc3vc3e"
;frun2="/raid4/tcox/vc3vc3e_no"
frun1="/raid4/tcox/vc3vc3h"
frun2="/raid4/tcox/vc3vc3h_no"
snapnum= 30

msg1='40% gas'
msg2='40% gas, no BH'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;--------------------------------------
;--------------------------------------

xaxistitle= 'Stellar Age (Gyr)'
xticklbls=['0','0.5','1.0','1.5','2.0','2.5','3.0']
xticknum= n_elements(xticklbls)-1
xmax= 3.0
xmin= 0.0

ymax = 1.2
ymin = 0.0

;----------------------------------------
; Generate plot
;----------------------------------------

x0= 0.02
y0= 0.15
x_size= 0.96
y_size= 0.83

!p.position=[x0, y0, x0+x_size,y0+y_size]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, /nodata, $
        color= 0, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, charthick=3.0, $
        xticks= xticknum, yticks= 1, ytickformat='(a1)', xtitle=xaxistitle, ytitle=yaxistitle, $
	xtickname=xticklbls


; --------------------

ok= fload_snapshot_bh(frun1,snapnum)
stellarage= fload_newstars_age(1)

sage1=stellarage


; --------------------

ok= fload_snapshot_bh(frun2,snapnum)
stellarage= fload_newstars_age(1)

sage2=stellarage


; --------------------



levels= 50.0


xyouts, 0.67, 0.85, msg1, /normal, charthick=3.0, size=1.5, color=150
xyouts, 0.67, 0.79, msg2, /normal, charthick=3.0, size=1.5, color=50



; --------------------


step= (xmax-xmin)/(levels)
bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
;bins=[bins,xmax]
bins=bins+(step*0.5)


; std histogram
hist_sage1= histogram(sage1, binsize=step, max=xmax, min=xmin)

; std histogram
hist_sage2= histogram(sage2, binsize=step, max=xmax, min=xmin)


; make sure it extends to 0
bins= [0,bins]
hist_sage1= [hist_sage1(0),hist_sage1]
hist_sage2= [hist_sage2(0),hist_sage2]


normalization= float(max([hist_sage1,hist_sage2]))
print, "histogram normalization= ", normalization
hist_sage1= hist_sage1/normalization
hist_sage2= hist_sage2/normalization

oplot, bins, hist_sage1, psym=10, color=150, thick=4.0
oplot, bins, hist_sage2, psym=10, color=50, thick=8.0

; fill in histogram
; ------------------
nbins= bins+(step*0.5)          ; make x coord
nbins[0]= 0.0
nbins=[nbins,nbins]
nbins= nbins(sort(nbins)) 

ntts= fltarr(2.*levels + 2) 
nctts= fltarr(2.*levels + 2) 
for i=1,levels do begin
   ntts[2.*i-1]= hist_sage1[i]
   ntts[2.*i]= hist_sage1[i]
   nctts[2.*i-1]= hist_sage2[i]
   nctts[2.*i]= hist_sage2[i]
endfor

; collisionless
;
;polyfill, nbins, nctts, /data, color= 50, /line_fill, linestyle=0, $
;				thick=3.0
;polyfill, nbins, nctts, /data, color= 50, /fill, linestyle=0, $
;				thick=3.0
;polyfill, nbins, nctts, /data, color= 220, /fill, linestyle=0, $
;				thick=3.0

; 40% gas
;polyfill, nbins, ntts, /data, color= 150, /line_fill, linestyle=0, $
;				thick=3.0, orientation=90.0
polyfill, nbins, ntts, /data, color= 150, /line_fill, linestyle=0, $
				thick=3.0, orientation=45.0



device, /close


end









;===================================================================================
;===================================================================================
;===================================================================================
;===================================================================================







pro jz_hist, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "jz_hist, junk"
   print, "  "
   return
endif

filename='jzhist.eps' 

;frun1="/raid4/tcox/vc3vc3f"
;frun2="/raid4/tcox/vc3vc3f_no"
;frun1="/raid4/tcox/vc3vc3e"
;frun2="/raid4/tcox/vc3vc3e_no"
frun1="/raid4/tcox/vc3vc3h"
frun2="/raid4/tcox/vc3vc3h_no"
snapnum= 30

msg1='40% gas'
msg2='40% gas, no BH'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;--------------------------------------
;--------------------------------------

;xaxistitle= 'j!Dz!N (km sec!E-1!N)'
xaxistitle= 'Normalized Specific j!Dtot!N'
xticklbls=['0','0.2','0.4','0.6','0.8','1']
xticknum= n_elements(xticklbls)-1
xmax= 1.0
xmin= 0.0

ymax = 1.2
ymin = 0.0

;----------------------------------------
; Generate plot
;----------------------------------------

x0= 0.02
y0= 0.15
x_size= 0.96
y_size= 0.83

!p.position=[x0, y0, x0+x_size,y0+y_size]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, /nodata, $
        color= 0, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, charthick=3.0, $
        xticks= xticknum, yticks= 1, ytickformat='(a1)', xtitle=xaxistitle, ytitle=yaxistitle, $
	xtickname=xticklbls


; --------------------

ok= fload_snapshot_bh(frun1,snapnum)

; --------------------

;  New Stars
; -----------

; stellar age
;sage1= fload_newstars_age(1)

;specific j's
ns1_jx= fload_newstars_j(21)
ns1_jy= fload_newstars_j(22)
ns1_jz= fload_newstars_j(23)

print, "J_tot= ", total(ns1_jx), total(ns1_jy), total(ns1_jz)

ns1_jtot= sqrt(ns1_jx*ns1_jx + ns1_jy*ns1_jy + ns1_jz*ns1_jz)

sage1= ns1_jtot

print,"New J min/max=", min(ns1_jtot), max(ns1_jtot)


; --------------------

;  Old Stars
; -----------

;sage2=stellarage

;specific j's
ns2_jx= fload_oldstars_j(21)
ns2_jy= fload_oldstars_j(22)
ns2_jz= fload_oldstars_j(23)

print, "J_tot= ", total(ns2_jx), total(ns2_jy), total(ns2_jz)

ns2_jtot= sqrt(ns2_jx*ns2_jx + ns2_jy*ns2_jy + ns2_jz*ns2_jz)

sage2= ns2_jtot

print,"Old J min/max=", min(ns2_jtot), max(ns2_jtot)


; --------------------

; Normalize this
normalizeJ= max([sage1,sage2])
sage1=sage1/normalizeJ
sage2=sage2/normalizeJ


; --------------------



levels= 50.0


xyouts, 0.67, 0.85, 'New Stars', /normal, charthick=3.0, size=1.5, color=150
xyouts, 0.67, 0.79, 'Old Stars', /normal, charthick=3.0, size=1.5, color=50



; --------------------


step= (xmax-xmin)/(levels)
bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
;bins=[bins,xmax]
bins=bins+(step*0.5)


; std histogram
hist_sage1= histogram(sage1, binsize=step, max=xmax, min=xmin)

; second histogram
hist_sage2= histogram(sage2, binsize=step, max=xmax, min=xmin)


; make sure it extends to 0
bins= [0,bins]
hist_sage1= [hist_sage1(0),hist_sage1]
hist_sage2= [hist_sage2(0),hist_sage2]


; make y log
hist_sage1= float(hist_sage1)
hist_sage2= float(hist_sage2)
idx=where(hist_sage1 gt 0)
hist_sage1(idx)= alog10(hist_sage1(idx))
idx=where(hist_sage2 gt 0)
hist_sage2(idx)= alog10(hist_sage2(idx))


normalization= float(max([hist_sage1]))
normalization= float(max([hist_sage1,hist_sage2]))
print, "histogram normalization= ", normalization
hist_sage1= hist_sage1/normalization
hist_sage2= hist_sage2/normalization

oplot, bins, hist_sage1, psym=10, color=150, thick=4.0
oplot, bins, hist_sage2, psym=10, color=50, thick=8.0

; fill in histogram
; ------------------
nbins= bins+(step*0.5)          ; make x coord
nbins[0]= 0.0
nbins=[nbins,nbins]
nbins= nbins(sort(nbins)) 

ntts= fltarr(2.*levels + 2) 
nctts= fltarr(2.*levels + 2) 
for i=1,levels do begin
   ntts[2.*i-1]= hist_sage1[i]
   ntts[2.*i]= hist_sage1[i]
   ;nctts[2.*i-1]= hist_sage2[i]
   ;nctts[2.*i]= hist_sage2[i]
endfor

; collisionless
;
;polyfill, nbins, nctts, /data, color= 50, /line_fill, linestyle=0, $
;				thick=3.0
;polyfill, nbins, nctts, /data, color= 50, /fill, linestyle=0, $
;				thick=3.0
;polyfill, nbins, nctts, /data, color= 220, /fill, linestyle=0, $
;				thick=3.0

; 40% gas
;polyfill, nbins, ntts, /data, color= 150, /line_fill, linestyle=0, $
;				thick=3.0, orientation=90.0
polyfill, nbins, ntts, /data, color= 150, /line_fill, linestyle=0, $
				thick=3.0, orientation=45.0



device, /close


end









;===================================================================================
;===================================================================================
;===================================================================================
;===================================================================================







pro j_age, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "j_age, junk"
   print, "  "
   return
endif

filename='jage.eps' 

frun1="/raid4/tcox/vc3vc3f"
;frun2="/raid4/tcox/vc3vc3f_no"
;frun1="/raid4/tcox/vc3vc3e"
;frun2="/raid4/tcox/vc3vc3e_no"
;frun1="/raid4/tcox/vc3vc3h"
;frun2="/raid4/tcox/vc3vc3h_no"
snapnum= 30

;msg1='40% gas'
;msg2='40% gas, no BH'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=4
setup_plot_stuff, 'ps', filename=filename, colortable=3


;--------------------------------------
;--------------------------------------


bins= 40

yaxistitle= 'Stellar Age (Gyr)'
xmax= 3.0
xmin= 0.0

yaxistitle= 'j!Dz!N (km sec!E-1!N)'
ymax = 1.2
ymin = 0.0


; --------------------

ok= fload_snapshot_bh(frun1,snapnum)

; stellar age
sage1= fload_newstars_age(1)

;specific j's
ns1_jx= fload_newstars_j(21)
ns1_jy= fload_newstars_j(22)
ns1_jz= fload_newstars_j(23)


ns1_jtot= sqrt(ns1_jx*ns1_jx + ns1_jy*ns1_jy + ns1_jz*ns1_jz)

ns1_jtot= ns1_jtot/max(ns1_jtot)
print, 'N= ',n_elements(ns1_jtot), max(ns1_jtot), min(ns1_jtot)
; --------------------



;----------------------------------------
; Generate plot
;----------------------------------------


contour_makegeneralpic, sage1, ns1_jtot, xmax, xmin, ymax, ymin, $
				pixels= bins, $
				NxNImage=NxNImage



;---------------------------------------
;  Print Tv image
;-----------------------------------------

x0= 0.15
y0= 0.15
x_size= 0.80
y_size= 0.80
           
                
tv, NxNImage, x0, y0, xsize=x_size, ysize=y_size, /normal



;------------------------
; Now the rest of it
;------------------------

!p.position=[x0, y0, x0+x_size,y0+y_size]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, /nodata, $
        color= 0, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, charthick=3.0, $
        ;xticks= xticknum, yticks= 1, ytickformat='(a1)', xtitle=xaxistitle, ytitle=yaxistitle, $
        xtitle=xaxistitle, ytitle=yaxistitle



; --------------------


device, /close


end









;=================================================================================






pro compile_js, junk


determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3b", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3c", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3d", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3e", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3f", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3g", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3h", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3i", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3j", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3k", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3l", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3m", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3n", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3o", 30
determine_one_sim_js, "/raid4/tcox/collisionless/cvc3vc3p", 30

end





;
;
;-------------------------------------------------
pro determine_one_sim_js, frun, snapnum

; ----------------------
;  Now Do the Work
; ----------------------

ok=fload_snapshot_bh(frun,snapnum)



; open angular momentum file
snaplbl= strcompress(string(snapnum),/remove_all)
openw, 1, frun+'/j_'+snaplbl+'.txt', ERROR=err
printf, 1, "#   "
printf, 1, "#  Angular Momentum  (Gadget Units) "
printf, 1, "#   "
printf, 1, "#              J_tot        J_x        J_y        J_z      theta       phi  "
printf, 1, "#               (GU)       (GU)       (GU)       (GU)     (deg.)     (deg.) "



;  All Stars
; -----------
print, "---------"
print, "all stars"
print, " "

r_as= fload_allstars_xyz('r')

; total j's
jx= fload_allstars_j(1)   ; j_x
jy= fload_allstars_j(2)   ; j_y
jz= fload_allstars_j(3)   ; j_z
jtot= sqrt(jx*jx + jy*jy + jz*jz)
jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
jphi= atan(jy,jx) * 180.0 / !PI
print, "tot= ", jtot, "  (",jx,jy,jz," )"
print, "theta= ", jtheta
print, "phi= ", jphi

printf, 1, FORMAT= '("AllStars   ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi



;  New Stars
; -----------
print, "---------"
print, "new stars"
print, " "

r_ns= fload_newstars_xyz('r')

; total j's
jx= fload_newstars_j(1)   ; j_x
jy= fload_newstars_j(2)   ; j_y
jz= fload_newstars_j(3)   ; j_z
jtot= sqrt(jx*jx + jy*jy + jz*jz)
jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
jphi= atan(jy,jx) * 180.0 / !PI
print, "tot= ", jtot, "  (",jx,jy,jz," )"
print, "theta= ", jtheta
print, "phi= ", jphi

printf, 1, FORMAT= '("NewStars   ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

; galaxy 1
startid= 1L
numpart= 200000L
jx= fload_1gal_newstars_j(1,startid=startid,numpart=numpart)   ; j_x
jy= fload_1gal_newstars_j(2,startid=startid,numpart=numpart)   ; j_y
jz= fload_1gal_newstars_j(3,startid=startid,numpart=numpart)   ; j_z
jtot= sqrt(jx*jx + jy*jy + jz*jz)
jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
jphi= atan(jy,jx) * 180.0 / !PI
print, "tot= ", jtot, "  (",jx,jy,jz," )"
;print, "theta= ", jtheta
;print, "phi= ", jphi

printf, 1, FORMAT= '("NewStars#1 ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

; galaxy 2
startid= 200001L
numpart= 200000L
jx= fload_1gal_newstars_j(1,startid=startid,numpart=numpart)   ; j_x
jy= fload_1gal_newstars_j(2,startid=startid,numpart=numpart)   ; j_y
jz= fload_1gal_newstars_j(3,startid=startid,numpart=numpart)   ; j_z
jtot= sqrt(jx*jx + jy*jy + jz*jz)
jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
jphi= atan(jy,jx) * 180.0 / !PI
print, "tot= ", jtot, "  (",jx,jy,jz," )"
;print, "theta= ", jtheta
;print, "phi= ", jphi

printf, 1, FORMAT= '("NewStars#2 ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi



;  Old Stars
; -----------
print, "---------"
print, "old stars"
print, " "
r_old= fload_oldstars_xyz('r')

;specific j's
jx= fload_oldstars_j(1)
jy= fload_oldstars_j(2)
jz= fload_oldstars_j(3)
jtot= sqrt(jx*jx + jy*jy + jz*jz)
jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
jphi= atan(jy,jx) * 180.0 / !PI
print, "tot= ", jtot, "  (",jx,jy,jz," )"
print, "theta= ", jtheta
print, "phi= ", jphi

printf, 1, FORMAT= '("OldStars   ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

; galaxy 1
startid= 1L
numpart= 200000L
jx= fload_1gal_oldstars_j(1,startid=startid,numpart=numpart)   ; j_x
jy= fload_1gal_oldstars_j(2,startid=startid,numpart=numpart)   ; j_y
jz= fload_1gal_oldstars_j(3,startid=startid,numpart=numpart)   ; j_z
jtot= sqrt(jx*jx + jy*jy + jz*jz)
jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
jphi= atan(jy,jx) * 180.0 / !PI
print, "tot= ", jtot, "  (",jx,jy,jz," )"
;print, "theta= ", jtheta
;print, "phi= ", jphi

printf, 1, FORMAT= '("OldStars#1 ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi

; galaxy 2
startid= 200001L
numpart= 200000L
jx= fload_1gal_oldstars_j(1,startid=startid,numpart=numpart)   ; j_x
jy= fload_1gal_oldstars_j(2,startid=startid,numpart=numpart)   ; j_y
jz= fload_1gal_oldstars_j(3,startid=startid,numpart=numpart)   ; j_z
jtot= sqrt(jx*jx + jy*jy + jz*jz)
jtheta= atan(sqrt(jx*jx + jy*jy), jz) * 180.0 / !PI
jphi= atan(jy,jx) * 180.0 / !PI
print, "tot= ", jtot, "  (",jx,jy,jz," )"
;print, "theta= ", jtheta
;print, "phi= ", jphi

printf, 1, FORMAT= '("OldStars#2 ", 6(F9.2,"  "))', $
		jtot, jx, jy, jz, jtheta, jphi




close, 1


;--------------------


end



