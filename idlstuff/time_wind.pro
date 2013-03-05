
; =================================================
;
;  This section analyzes the properties of the
;  hot, outflowing gas.
;
;
;
;
; ==================================================
pro time_windproperties, frun


if not keyword_set(frun) then begin
        print, " "
        print, " time_windproperties, frun"
        print, " "
        print, " "
        print, " "
        print, " "
        return
endif

; this assumes they are ordered
;  0 through


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])



; --------------
; Parameters
; --------------

r_wind= 100.0      ; used in wind_gtx


; ---------------
; Arrays
; ---------------
time= fltarr(nsnaps)

mass_gtx=  fltarr(nsnaps)
mass_egt0=  fltarr(nsnaps)

mass_wind= fltarr(nsnaps)
wind_ke= fltarr(nsnaps)
wind_radial_ke=   fltarr(nsnaps)
wind_z=fltarr(nsnaps)


;t_merg= fload_mergertime(frun)


; ----------------------------------------
; This part loops through the snapshots
; and compiles various information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"


        ; open snapshot
        ;ok=fload_snapshot(frun,snapnum)
        ok=fload_snapshot_bh(frun,i)


        ; what time is it?
        time[i]= fload_time(1)
        TTime= float(time[i])


	; gas properties
	; -----------------
	rad= fload_gas_xyz('r')
	gas_mass= fload_gas_mass(1)

	;e= fload_gas_energy()
	ke= fload_gas_energy(1,/kinetic)
	ke_radial= fload_gas_energy(1,/radial_kinetic)
	pe= fload_gas_energy(1,/potential)
	the= fload_gas_energy(1,/thermal)
	e= ke + pe + the
	temp= fload_gas_temperature(1)
	keV= fload_gas_temperature(1,/keV)
	entropy= fload_gas_entropy(1)

	; not really used at current, and slightly
	; outdated
	;xr= fload_gas_xray_luminosity(1)

	metallicity= fload_gas_metallicity(1)


	; define wind
	; --------------

	; greater than r_wind
	gtr_idx= where(rad gt r_wind)
	if gtr_idx(0) ne -1 then mass_gtx[i]= total(gas_mass(gtr_idx)) else mass_gtx[i]= 0

	; positive energy
	egt0_idx= where(e gt 0.0)
	if egt0_idx(0) ne -1 then mass_egt0[i]= total(gas_mass(egt0_idx)) else mass_egt0[i]= 0

	; wind index
	;wind_idx= gtr_idx
	wind_idx= egt0_idx


	; wind properties
	; -----------------
	if n_elements(wind_idx) gt 1 then begin
		mass_wind[i]= total(gas_mass(wind_idx))

		z_moment= moment(metallicity(wind_idx))
		print, "Metallicity= ",z_moment[0]," +/- ", sqrt(z_moment[1])
		wind_z[i]= z_moment(0)/0.02

		wind_ke[i]= alog10(total(ke(wind_idx)))

		wind_radial_ke[i]= alog10(total(ke_radial(wind_idx)))
	endif else begin
		mass_wind[i]= 0
		wind_z[i]= 0
		wind_ke[i]= 0
		wind_radial_ke[i]= 0
	endelse


endfor

; ----------------------------------------
; write to file

openw, 1, frun+'/wind.txt', ERROR=err
        
printf, 1, "#   wind.txt"     
printf, 1, "#       --> masses in 10^10 msolar,energy in log ers"
printf, 1, "#               "
printf, 1, "# time     gtx      egt0      w_e      wrke        z  "               
printf, 1, "# (Gyr)    (m)       (m)      (e)       (e)    (solar) "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3," ",F8.5,"  ",F8.5," ", 2(F8.3," "),F8.5)', $
                time[i], mass_gtx[i], mass_egt0[i], $
                        wind_ke[i], wind_radial_ke[i], wind_z[i]
endfor      
close, 1    
            
            
        

; ----------------------------------------
        
print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"
        




end












; ----------------------------
;  Read wind.txt file
; ----------------------------
pro read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z

filename= frun+'/wind.txt'

spawn, "wc "+filename,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then filedata= fltarr(6,lines)

openr, 1, filename

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, filedata
close, 1


time= filedata[0,*]
mass_gtx= filedata[1,*]
mass_egt= filedata[2,*]
wind_ke= filedata[3,*]
wind_rke= filedata[4,*]
wind_z= filedata[5,*]


end






;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro plot_wind_comparison, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_wind_comparison, junk"
	print, " "
	print, " "
	return
endif

filename='wind.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;yaxistitle = "Mass (10!E10!N M!D!9n!E!N)"
yaxistitle = "!6Log Wind Mass (M!D!9n!6!N)"
;yaxistitle = "Log Mass (M)"
;xaxistitle = "Time after merger (Gyr)"
xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
;xmax = max(time)
;xmax = 4.5
xmax = 3.0
;xmax = 2.0
xmin = 0
;ymax = 8.0
;ymin = 0.0
;ymin = 5.0e-4    ; 10^7 solar masses
;ymax = 12.0
ymax = 11.0
;ymax = 10.0
ymin = 7.0
;ymin = 6.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




x0= 0.65
y0= 0.35

; Galaxy 1
;----------
;frun="pool/vc3vc3"
;frun="pool/vc3bvc3b"
;frun="/raid4/tcox/vc3avc3a_3"
;frun="/raid4/tcox/vc3ba_2_3"
;frun="/raid4/tcox/vc1vc1"
;frun="/raid4/tcox/As/A1"
;frun="/raid4/tcox/vc3vc3e_2"
;frun="/raid4/tcox/vc3vc3e"
;frun="/raid4/tcox/vc3vc3h"

;process_one_wind, "/raid4/tcox/vc3vc3h", lcolor=50, msg='std', x0= 0.65, y0= 0.35



; Galaxy 2
;----------
;frun="pool/vc3vc3_wBH"
;frun="pool/vc3bvc3b_wBH"
;frun="/raid4/tcox/vc3avc3a_2"
;frun="/raid4/tcox/vc3ba_2_2"
;frun="/raid4/tcox/vc2vc2"
;frun="/raid4/tcox/As/A2"
;frun="/raid4/tcox/vc3vc3e_no"
;frun="/raid4/tcox/vc3vc3e"

;process_one_wind, "/raid4/tcox/vc3vc3e", lcolor=150, msg='tilted', x0= 0.65, y0= 0.31




; Galaxy 3
;----------
;frun="/raid4/tcox/vc3avc3a_1"
;frun="/raid4/tcox/vc3ba_2_1"
;frun="/raid4/tcox/vc3bvc3b"
;frun="/raid4/tcox/vc3vc3"
;frun="/raid4/tcox/As/A3"


; Galaxy 4
;----------
;frun="/raid4/tcox/vc3avc3a_4"
;frun="/raid4/tcox/vc3ba_2_4"
;frun="/raid4/tcox/vc4vc4a"
;frun="/raid4/tcox/As/A4"


; Galaxy 5
;----------
;frun="/raid4/tcox/vc3avc3a_6"      ; vc3avc3a_5 isn't there, don't know why
;frun="/raid4/tcox/vc3ba_2_5"
;frun="/raid4/tcox/vc5vc5a"
;frun="/raid4/tcox/As/A5"


; Galaxy 6
;----------
;frun="/raid4/tcox/vc6vc6a"


;----------
;xyouts, 0.70, 0.90, 'bulge', /normal, charthick=3, size=1.33, color=0


process_one_wind, "/raid4/tcox/vc3vc3e_2", lcolor=150, lpsym= 5, msg='std', x0= 0.54, y0= 0.34

process_one_wind, "/raid4/tcox/vc3vc3e_sb8", lcolor=200, msg='sb8 - small sfr', x0= 0.54, y0= 0.30
process_one_wind, "/raid4/tcox/vc3vc3e_sb10", lcolor=100, msg='sb10 - mod. sfr', x0= 0.54, y0= 0.26
process_one_wind, "/raid4/tcox/vc3vc3e_sb13", lcolor=50, msg='sb13 - osc. sfr', x0= 0.54, y0= 0.22
process_one_wind, "/raid4/tcox/vc3vc3e_sb1", lcolor=0, msg='sb1 - no change', x0= 0.54, y0= 0.18




device, /close




end






pro process_one_wind, frun, lcolor=lcolor, lpsym=lpsym, $
				msg=msg, x0=x0, y0=y0

	if not keyword_set(lcolor) then lcolor= 0
	if not keyword_set(lpsym) then lpsym= 3


	read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
	mass_egt= 10.0 + alog10(mass_egt/0.7)
	mass_gtx= 10.0 + alog10(mass_gtx/0.7)
	;time=time/0.7

	;oplot, time, mass_egt, thick=5.0, psym=-lpsym, color= lcolor
	;oplot, time, mass_gtx, thick=5.0, psym=-lpsym, color= lcolor, linestyle=1
	oplot, time, mass_gtx, thick=5.0, psym=-lpsym, color= lcolor

	if not keyword_set(x0) then x0= 0.75
	if not keyword_set(y0) then y0= 0.85

	if keyword_set(msg) then begin
		xyouts, x0, y0, msg, /normal, charthick=3, size=1.33, color=lcolor
	endif else begin
		xyouts, x0, y0, fload_getid(frun), /normal, charthick=1, size=1.33, color=150
	endelse



end








;==========================================================================







;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro plot_wind_z_comparison, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_wind_z_comparison, junk"
	print, " "
	print, " "
	return
endif

filename='windz.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



yaxistitle = "Metallicity (Z!D!9n!3!N)"
xaxistitle = "Time (Gyr)"
;xmax = max(time)
;xmax = 4.5
xmax = 3.0
;xmax = 2.0
xmin = 0
ymax = 7.0
;ymin = 0.0
;ymin = 1.0e-3    ; in solar units
ymin = 1.0e-3    ; in solar units


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


label_x0= 0.70
label_y0= 0.40

; Galaxy 1
;----------
;frun="pool/vc3vc3"
;frun="pool/vc3bvc3b"
;frun="/raid4/tcox/vc3avc3a_3"
;frun="/raid4/tcox/vc3ba_2_3"
;frun="/raid4/tcox/vc1vc1"
;frun="/raid4/tcox/As/A1"
frun="/raid4/tcox/vc3vc3e_2"
read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
time= time/0.7
oplot, time, wind_z, thick=3.0, psym=-2, color= 50

xyouts, label_x0, label_y0, fload_getid(frun), /normal, charthick=1, size=1.33, color=50


; Galaxy 2
;----------
;frun="pool/vc3vc3_wBH"
;frun="pool/vc3bvc3b_wBH"
;frun="/raid4/tcox/vc3avc3a_2"
;frun="/raid4/tcox/vc3ba_2_2"
;frun="/raid4/tcox/vc2vc2"
;frun="/raid4/tcox/As/A2"
frun="/raid4/tcox/vc3vc3e_no"
read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
time= time/0.7
oplot, time, wind_z, thick=3.0, psym=-6, color= 75

xyouts, label_x0, label_y0-0.04, fload_getid(frun), /normal, charthick=1, size=1.33, color=75


; Galaxy 3
;----------
;frun="/raid4/tcox/vc3avc3a_1"
;frun="/raid4/tcox/vc3bvc3b"
;frun="/raid4/tcox/vc3ba_2_1"
;frun="/raid4/tcox/vc3vc3"
;frun="/raid4/tcox/As/A3"
;read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
;time=time/0.7
;oplot, time, wind_z, thick=3.0, psym=-6, color= 100
;xyouts, label_x0, label_y0-0.08, fload_getid(frun), /normal, charthick=1, size=1.33, color=100


; Galaxy 4
;----------
;frun="/raid4/tcox/vc3avc3a_4"
;frun="/raid4/tcox/vc3ba_2_4"
;frun="/raid4/tcox/vc4vc4a"
;frun="/raid4/tcox/As/A4"
;read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
;time=time/0.7
;oplot, time, wind_z, thick=3.0, psym=-6, color= 125
;xyouts, label_x0, label_y0-0.12, fload_getid(frun), /normal, charthick=1, size=1.33, color=125


; Galaxy 5
;----------
;frun="/raid4/tcox/vc3avc3a_6"      ; vc3avc3a_5 isn't there, don't know why
;frun="/raid4/tcox/vc3ba_2_5"
;frun="/raid4/tcox/vc5vc5a"
;frun="/raid4/tcox/As/A5"
;read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
;time=time/0.7
;oplot, time, wind_z, thick=3.0, psym=-6, color= 150
;xyouts, label_x0, label_y0-0.16, fload_getid(frun), /normal, charthick=1, size=1.33, color=150


; Galaxy 6
;----------
;frun="/raid4/tcox/vc6vc6a"
;read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
;time=time/0.7
;oplot, time, wind_z, thick=3.0, psym=-6, color= 175
;xyouts, label_x0, label_y0-0.20, fload_getid(frun), /normal, charthick=1, size=1.33, color=175


;----------
;xyouts, 0.70, 0.92, 'bulge', /normal, charthick=1, size=1.33, color= 0


device, /close




end





;--------------------------------------------------------------------------------


;===================================
;
;  Wind Relations
;
;===================================


; wind metallicity versus
; black hole seed mass
;
pro plot_wind_bhmet, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_wind_bhm, junk"
	print, " "
	print, " "
	return
endif

;do_bhmass= 1
do_bhmass= 0
do_galsizes= 1
;do_galsizes= 0

filename='wind_info.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



yaxistitle = "Wind Metallicity (Z!D!9n!3!N)"
if do_bhmass eq 1 then xaxistitle = "Log BH Seed Mass (M!D!9n!3!N)"
if do_galsizes eq 1 then xaxistitle = "Log Halo Mass (M!D!9n!3!N)"
;xmax = max(time)
if do_bhmass eq 1 then xmax = 7.2
if do_bhmass eq 1 then xmin = 2.8
if do_galsizes eq 1 then xmax = 13.6
if do_galsizes eq 1 then xmin = 10.6
ymax = 7.0
ymin = 1.0e-2    ; 10^-3 solar


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; bulge bhmass vc3 runs
; ----------------------
if do_bhmass eq 1 then begin
   frun="/raid4/tcox/vc3ba_2_3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc3ba_2_2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc3ba_2_1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc3ba_2_4"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc3ba_2_5"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]


   xpoints= [3.0,4.0,5.0,6.0,7.0]
   oplot, xpoints, wind_zs, thick=3.0, psym=-2, color= 50
   oplot, [3.0,3.2], [5.0,5.0], thick=3.0, psym=-2, color= 50
   xyouts, 0.28, 0.90, 'bulge', /normal,charthick=3.0, size=1.33, color=50
endif


; bulge-less bhmass vc3 runs
; ----------------------------
if do_bhmass eq 1 then begin
   frun="/raid4/tcox/vc3avc3a_3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc3avc3a_2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc3avc3a_1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc3avc3a_4"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc3avc3a_6"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]


   xpoints= [3.0,4.0,5.0,6.0,7.0]
   oplot, xpoints, wind_zs, thick=3.0, psym=-5, color= 150
   oplot, [3.0,3.2], [3.5,3.5], thick=3.0, psym=-5, color= 150
   xyouts, 0.28, 0.85, 'no bulge', /normal, charthick=1, size=1.33, color=150
endif


; Galaxy Masses
;---------------
if do_galsizes eq 1 then begin
   frun="/raid4/tcox/vc1vc1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc2vc2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc3vc3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc4vc4a"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc5vc5a"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/vc6vc6a"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]

   xpoints= [8.2e10,23.8e10,70.8e10,190.4e10,536.8e10,1523.8e10,5813.0e10]
   xpoints= alog10(xpoints)
   oplot, xpoints, wind_zs, thick=3.0, psym=-2, color= 50
   oplot, [10.8,11.0], [5.0,5.0], thick=3.0, psym=-2, color= 50
   xyouts, 0.33, 0.90, '20% gas', /normal,charthick=3.0, size=1.33, color=50
endif


; Galaxy Masses
;---------------
doggytown=1
if do_galsizes eq 1 and doggytown eq 1 then begin
   frun="/raid4/tcox/As/A1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/As/A2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/As/A3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/As/A4"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]
   frun="/raid4/tcox/As/A5"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   wind_zs= [wind_zs,wind_z(n_elements(time)-1)]


   xpoints= [23.8e10,67.1e10,190.5e10,536.8e10,1523.8e10]
   xpoints= alog10(xpoints)
   oplot, xpoints, wind_zs, thick=3.0, psym=-5, color= 150
   oplot, [10.8,11.0], [3.3,3.3], thick=3.0, psym=-5, color= 150
   xyouts, 0.31, 0.85, '100% gas', /normal, charthick=1, size=1.33, color=150
endif


;----------
device, /close

end





; wind mass versus
; black hole seed mass
;
pro plot_wind_bhmass, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_wind_bhmass, junk"
	print, " "
	print, " "
	return
endif

;do_bhmass= 1
do_bhmass= 0
do_galsizes= 1
;do_galsizes= 0

filename='wind_info.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



yaxistitle = "Log Wind Mass (M!D!9n!3!N)"
if do_bhmass eq 1 then xaxistitle = "Log BH Seed Mass (M!D!9n!3!N)"
if do_galsizes eq 1 then xaxistitle = "Log Halo Mass (M!D!9n!3!N)"
;xmax = max(time)
if do_bhmass eq 1 then xmax = 7.2
if do_bhmass eq 1 then xmin = 2.8
if do_galsizes eq 1 then xmax = 13.6
if do_galsizes eq 1 then xmin = 10.6
ymax = 11.3
ymin = 7.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; bulge bhmass vc3 runs
; ----------------------
if do_bhmass eq 1 then begin
   frun="/raid4/tcox/vc3ba_2_3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc3ba_2_2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc3ba_2_1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc3ba_2_4"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc3ba_2_5"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]


   xpoints= [3.0,4.0,5.0,6.0,7.0]
   oplot, xpoints, wind_wm, thick=3.0, psym=-2, color= 50
   oplot, [3.0,3.2], [10.12,10.12], thick=3.0, psym=-2, color= 50
   xyouts, 0.28, 0.90, 'bulge', /normal,charthick=3.0, size=1.33, color=50
endif


; bulge-less bhmass vc3 runs
; ----------------------------
if do_bhmass eq 1 then begin
   frun="/raid4/tcox/vc3avc3a_3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc3avc3a_2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc3avc3a_1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc3avc3a_4"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc3avc3a_6"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]


   xpoints= [3.0,4.0,5.0,6.0,7.0]
   oplot, xpoints, wind_wm, thick=3.0, psym=-5, color= 150
   oplot, [3.0,3.2], [9.9,9.9], thick=3.0, psym=-5, color= 150
   xyouts, 0.28, 0.85, 'no bulge', /normal, charthick=1, size=1.33, color=150
endif


; Galaxy Masses
;---------------
if do_galsizes eq 1 then begin
   frun="/raid4/tcox/vc1vc1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc2vc2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc3vc3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc4vc4a"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc5vc5a"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/vc6vc6a"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]

   xpoints= [8.2e10,23.8e10,70.8e10,190.4e10,536.8e10,1523.8e10,5813.0e10]
   xpoints= alog10(xpoints)
   oplot, xpoints, wind_wm, thick=3.0, psym=-2, color= 50
   oplot, [10.8,11.0], [11.10,11.10], thick=3.0, psym=-2, color= 50
   xyouts, 0.33, 0.90, '20% gas', /normal, charthick=1, size=1.33, color=50
endif

;frun="pool/vc3vc3"
;frun="pool/vc3bvc3b"
;frun="/raid4/tcox/vc4avc4a"
;frun="/raid4/tcox/vc5avc5a"
;frun="/raid4/tcox/vc6avc6a"


; Galaxy Masses
;---------------
doggytown= 1
if do_galsizes eq 1 and doggytown eq 1 then begin
   frun="/raid4/tcox/As/A1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/As/A2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/As/A3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/As/A4"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]
   frun="/raid4/tcox/As/A5"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   mass_egt= 10.0 + alog10(mass_egt)
   wind_wm= [wind_wm,mass_egt(n_elements(time)-1)]

   xpoints= [23.8e10,67.1e10,190.5e10,536.8e10,1523.8e10]
   xpoints= alog10(xpoints)
   oplot, xpoints, wind_wm, thick=3.0, psym=-5, color= 150
   oplot, [10.8,11.0], [10.9,10.9], thick=3.0, psym=-5, color= 150
   xyouts, 0.31, 0.85, '100% gas', /normal, charthick=1, size=1.33, color=150
endif



;----------
device, /close

end




; wind mass fraction versus
; black hole seed mass
;
pro plot_wind_bhmf, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_wind_bhmf, junk"
	print, " "
	print, " "
	return
endif

;do_bhmass= 1
do_bhmass= 0
do_galsizes= 1
;do_galsizes= 0

filename='wind_info.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



yaxistitle = "Wind Mass Fraction"
if do_bhmass eq 1 then xaxistitle = "Log BH Seed Mass (M!D!9n!3!N)"
if do_galsizes eq 1 then xaxistitle = "Log Halo Mass (M!D!9n!3!N)"
;xmax = max(time)
if do_bhmass eq 1 then xmax = 7.2
if do_bhmass eq 1 then xmin = 2.8
if do_galsizes eq 1 then xmax = 13.6
if do_galsizes eq 1 then xmin = 10.6
ymax = 1.0
ymin = 0.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; bulge bhmass vc3 runs
; ----------------------
if do_bhmass eq 1 then begin
   frun="/raid4/tcox/vc3ba_2_3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3ba_2_2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3ba_2_1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3ba_2_4"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3ba_2_5"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]


   xpoints= [3.0,4.0,5.0,6.0,7.0]
   oplot, xpoints, wind_mf, thick=3.0, psym=-2, color= 50
   oplot, [3.0,3.2], [0.95,0.95], thick=3.0, psym=-2, color= 50
   xyouts, 0.28, 0.90, 'bulge', /normal,charthick=3.0, size=1.33, color=50
endif


; bulge-less bhmass vc3 runs
; ----------------------------
if do_bhmass eq 1 then begin
   frun="/raid4/tcox/vc3avc3a_3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3avc3a_2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3avc3a_1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3avc3a_4"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3avc3a_6"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]


   xpoints= [3.0,4.0,5.0,6.0,7.0]
   oplot, xpoints, wind_mf, thick=3.0, psym=-5, color= 150
   oplot, [3.0,3.2], [0.87,0.87], thick=3.0, psym=-5, color= 150
   xyouts, 0.28, 0.85, 'no bulge', /normal, charthick=1, size=1.33, color=150
endif


; Galaxy Masses
;----------------
if do_galsizes eq 1 then begin
   frun="/raid4/tcox/vc1vc1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc2vc2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3vc3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc4vc4a"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc5vc5a"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc6vc6a"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]


   xpoints= [8.2e10,23.8e10,70.8e10,190.4e10,536.8e10,1523.8e10,5813.0e10]
   xpoints= alog10(xpoints)
   oplot, xpoints, wind_mf, thick=3.0, psym=-2, color= 50
   oplot, [10.8,11.0], [0.94,0.94], thick=3.0, psym=-2, color= 50
   xyouts, 0.33, 0.90, '20% gas', /normal, charthick=3, size=1.33, color=50
endif



; Galaxy Masses
;----------------
doggytown= 1
if do_galsizes eq 1 and doggytown eq 1 then begin
   frun="/raid4/tcox/As/A1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/As/A2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/As/A3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/As/A4"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/As/A5"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]


   xpoints= [23.8e10,70.8e10,190.4e10,536.8e10,1523.8e10]
   xpoints= alog10(xpoints)
   oplot, xpoints, wind_mf, thick=3.0, psym=-5, color= 150
   oplot, [10.8,11.0], [0.89,0.89], thick=3.0, psym=-5, color= 150
   xyouts, 0.31, 0.85, '100% gas', /normal, charthick=3, size=1.33, color=150
endif


;----------
device, /close

end








; wind energy versus
; black hole seed mass
; or progenitor halo
pro plot_wind_energy, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_wind_energy, junk"
	print, " "
	print, " "
	return
endif

;do_bhmass= 1
do_bhmass= 0
do_galsizes= 1
;do_galsizes= 0

filename='wind_info.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



yaxistitle = "Log Wind Energy (ergs)"
if do_bhmass eq 1 then xaxistitle = "Log BH Seed Mass (M!D!9n!3!N)"
if do_galsizes eq 1 then xaxistitle = "Log Halo Mass (M!D!9n!3!N)"
;xmax = max(time)
if do_bhmass eq 1 then xmax = 7.2
if do_bhmass eq 1 then xmin = 2.8
if do_galsizes eq 1 then xmax = 13.6
if do_galsizes eq 1 then xmin = 10.6
ymax = 62.0
ymin = 55.0

Gadget_unit_energy= 1.989d+53

;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; bulge bhmass vc3 runs
; ----------------------
if do_bhmass eq 1 then begin
   frun="/raid4/tcox/vc3ba_2_3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3ba_2_2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3ba_2_1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3ba_2_4"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3ba_2_5"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]


   xpoints= [3.0,4.0,5.0,6.0,7.0]
   oplot, xpoints, wind_mf, thick=3.0, psym=-2, color= 50
   oplot, [3.0,3.2], [0.95,0.95], thick=3.0, psym=-2, color= 50
   xyouts, 0.28, 0.90, 'bulge', /normal,charthick=3.0, size=1.33, color=50
endif


; bulge-less bhmass vc3 runs
; ----------------------------
if do_bhmass eq 1 then begin
   frun="/raid4/tcox/vc3avc3a_3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3avc3a_2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3avc3a_1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3avc3a_4"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]
   frun="/raid4/tcox/vc3avc3a_6"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_mf= [wind_mf,mass_egt(n_elements(time)-1)/gas_tot(0)]


   xpoints= [3.0,4.0,5.0,6.0,7.0]
   oplot, xpoints, wind_mf, thick=3.0, psym=-5, color= 150
   oplot, [3.0,3.2], [0.87,0.87], thick=3.0, psym=-5, color= 150
   xyouts, 0.28, 0.85, 'no bulge', /normal, charthick=1, size=1.33, color=150
endif


; Galaxy Masses
;----------------
if do_galsizes eq 1 then begin
   frun="/raid4/tcox/vc1vc1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_e= [wind_ke(n_elements(time)-1)]
   frun="/raid4/tcox/vc2vc2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_e= [wind_e,wind_ke(n_elements(time)-1)]
   frun="/raid4/tcox/vc3vc3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_e= [wind_e,wind_ke(n_elements(time)-1)]
   frun="/raid4/tcox/vc4vc4a"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_e= [wind_e,wind_ke(n_elements(time)-1)]
   frun="/raid4/tcox/vc5vc5a"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_e= [wind_e,wind_ke(n_elements(time)-1)]
   frun="/raid4/tcox/vc6vc6a"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_e= [wind_e,wind_ke(n_elements(time)-1)]


   xpoints= [8.2e10,23.8e10,70.8e10,190.4e10,536.8e10,1523.8e10,5813.0e10]
   xpoints= alog10(xpoints)
   wind_e= wind_e + alog10(Gadget_unit_energy)
   oplot, xpoints, wind_e, thick=3.0, psym=-2, color= 50
   oplot, [10.8,11.0], [0.94,0.94], thick=3.0, psym=-2, color= 50
   xyouts, 0.33, 0.90, '20% gas', /normal, charthick=3, size=1.33, color=50
endif



; Galaxy Masses
;----------------
doggytown= 1
if do_galsizes eq 1 and doggytown eq 1 then begin
   frun="/raid4/tcox/As/A1"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_e= [wind_ke(n_elements(time)-1)]
   frun="/raid4/tcox/As/A2"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_e= [wind_e,wind_ke(n_elements(time)-1)]
   frun="/raid4/tcox/As/A3"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_e= [wind_e,wind_ke(n_elements(time)-1)]
   frun="/raid4/tcox/As/A4"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_e= [wind_e,wind_ke(n_elements(time)-1)]
   frun="/raid4/tcox/As/A5"
   read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
   ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, xray, xray_sf, gas_tot, gas_hot, gas_cold, gas_sf
   wind_e= [wind_e,wind_ke(n_elements(time)-1)]


   xpoints= [23.8e10,70.8e10,190.4e10,536.8e10,1523.8e10]
   xpoints= alog10(xpoints)
   wind_e= wind_e + alog10(Gadget_unit_energy)
   oplot, xpoints, wind_e, thick=3.0, psym=-5, color= 150
   oplot, [10.8,11.0], [0.89,0.89], thick=3.0, psym=-5, color= 150
   xyouts, 0.31, 0.85, '100% gas', /normal, charthick=3, size=1.33, color=150
endif


;----------
device, /close

end





