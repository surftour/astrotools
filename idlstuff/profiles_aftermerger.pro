pro aftermerger_masses, frun


if not keyword_set(frun) then begin
        print, " "
        print, " aftermerger_masses, frun"
        print, " "
        print, " needs: *.bh file"
	print, "        and"
	print, "       .run bh.pro"
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


x_rad= 5.0   ; follow mass inside of this radius

time= fltarr(nsnaps)

mass_gas=  fltarr(9,nsnaps)
mass_disk= fltarr(9,nsnaps)
mass_bulge=   fltarr(9,nsnaps)
mass_newstars= fltarr(9,nsnaps)
mass_darkmatter= fltarr(9,nsnaps)


;t_merg= fload_mergertime(1)
t_merg= bh_mergertime(frun)    ; this needs bh.pro file to be compiled

;t_merg= 1.15684  & print, "WARNING: manually setting merger time"


; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"


        ; open snapshot
        ;ok=fload_snapshot(frun,snapnum)
        ok=fload_snapshot_bh(frun,i)


        ; what time is it?
        time[i]= fload_time(1)
        TTime= float(time[i])

	if time[i] ge t_merg then begin

		; gas
		r= fload_gas_xyz('r')
		gmass= fload_gas_mass(1)
		mass_gas[*,i]= process_mass_within_r(r, gmass)

		; disk
		r= fload_disk_xyz('r')
		dmass= fload_disk_mass(1)
		mass_disk[*,i]= process_mass_within_r(r, dmass)

		; bulge
		r= fload_bulge_xyz('r')
		bmass= fload_bulge_mass(1)
		mass_bulge[*,i]= process_mass_within_r(r, bmass)

		; new stars
		r= fload_newstars_xyz('r')
		nsmass= fload_newstars_mass(1)
		mass_newstars[*,i]= process_mass_within_r(r, nsmass)

		; dark matter
		r= fload_halo_xyz('r')
		hmass= fload_halo_mass(1)
		mass_darkmatter[*,i]= process_mass_within_r(r, hmass)

	endif

endfor

; ----------------------------------------
; write to files
;
;
;   0 kpc  - essentially this will simply store the total mass
;  --------------------------------------------------------------
openw, 1, frun+'/massin_0_aftermerger.txt', ERROR=err

printf, 1, "#   massin_0_aftermerger.txt"
printf, 1, "#      ---> total mass of each component, in gadget units (10^10 msolar/h) "
printf, 1, "#        "
printf, 1, "# time   "
printf, 1, "#(Gyr)      gas          disk        bulge     newstars           dm  "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3, 5(F11.6,"  "))', $
                time[i], mass_gas[0,i], mass_disk[0,i], mass_bulge[0,i], mass_newstars[0,i], mass_darkmatter[0,i]
endfor
close, 1
;
;
;   1 kpc
;  --------
openw, 1, frun+'/massin_1_aftermerger.txt', ERROR=err
        
printf, 1, "#   massin_1_aftermerger.txt"     
printf, 1, "#      ---> mass of each component, in gadget units (10^10 msolar/h), within 1 kpc"
printf, 1, "#               of the galactic center as measured after the bh merger "
printf, 1, "# time   "               
printf, 1, "#(Gyr)      gas          disk        bulge     newstars           dm  "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3, 5(F11.6,"  "))', $
                time[i], mass_gas[1,i], mass_disk[1,i], mass_bulge[1,i], mass_newstars[1,i], mass_darkmatter[1,i]
endfor      
close, 1
;
;
;   2 kpc
;  --------
openw, 1, frun+'/massin_2_aftermerger.txt', ERROR=err

printf, 1, "#   massin_2_aftermerger.txt"
printf, 1, "#      ---> mass of each component, in gadget units (10^10 msolar/h), within 2 kpc"
printf, 1, "#               of the galactic center as measured after the bh merger "
printf, 1, "# time   "
printf, 1, "#(Gyr)      gas          disk        bulge     newstars           dm  "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3, 5(F11.6,"  "))', $
                time[i], mass_gas[2,i], mass_disk[2,i], mass_bulge[2,i], mass_newstars[2,i], mass_darkmatter[2,i]
endfor  
close, 1
;
;
;   3 kpc
;  --------
openw, 1, frun+'/massin_3_aftermerger.txt', ERROR=err

printf, 1, "#   massin_3_aftermerger.txt"
printf, 1, "#      ---> mass of each component, in gadget units (10^10 msolar/h), within 3 kpc"
printf, 1, "#               of the galactic center as measured after the bh merger "
printf, 1, "# time   "
printf, 1, "#(Gyr)      gas          disk        bulge     newstars           dm  "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3, 5(F11.6,"  "))', $
                time[i], mass_gas[3,i], mass_disk[3,i], mass_bulge[3,i], mass_newstars[3,i], mass_darkmatter[3,i]
endfor  
close, 1
;
;
;   4 kpc
;  --------
openw, 1, frun+'/massin_4_aftermerger.txt', ERROR=err

printf, 1, "#   massin_4_aftermerger.txt"
printf, 1, "#      ---> mass of each component, in gadget units (10^10 msolar/h), within 4 kpc"
printf, 1, "#               of the galactic center as measured after the bh merger "
printf, 1, "# time   "
printf, 1, "#(Gyr)      gas          disk        bulge     newstars           dm  "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3, 5(F11.6,"  "))', $
                time[i], mass_gas[4,i], mass_disk[4,i], mass_bulge[4,i], mass_newstars[4,i], mass_darkmatter[4,i]
endfor  
close, 1
;
;
;   5 kpc
;  --------
openw, 1, frun+'/massin_5_aftermerger.txt', ERROR=err

printf, 1, "#   massin_5_aftermerger.txt"
printf, 1, "#      ---> mass of each component, in gadget units (10^10 msolar/h), within 5 kpc"
printf, 1, "#               of the galactic center as measured after the bh merger "
printf, 1, "# time   "
printf, 1, "#(Gyr)      gas          disk        bulge     newstars           dm  "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3, 5(F11.6,"  "))', $
                time[i], mass_gas[5,i], mass_disk[5,i], mass_bulge[5,i], mass_newstars[5,i], mass_darkmatter[5,i]
endfor
close, 1
;
;
;   10 kpc
;  --------
openw, 1, frun+'/massin_10_aftermerger.txt', ERROR=err

printf, 1, "#   massin_10_aftermerger.txt"
printf, 1, "#      ---> mass of each component, in gadget units (10^10 msolar/h), within 10 kpc"
printf, 1, "#               of the galactic center as measured after the bh merger "
printf, 1, "# time   "
printf, 1, "#(Gyr)      gas          disk        bulge     newstars           dm  "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3, 5(F11.6,"  "))', $
                time[i], mass_gas[6,i], mass_disk[6,i], mass_bulge[6,i], mass_newstars[6,i], mass_darkmatter[6,i]
endfor
close, 1
;
;
;   15 kpc
;  --------
openw, 1, frun+'/massin_15_aftermerger.txt', ERROR=err

printf, 1, "#   massin_15_aftermerger.txt"
printf, 1, "#      ---> mass of each component, in gadget units (10^10 msolar/h), within 15 kpc"
printf, 1, "#               of the galactic center as measured after the bh merger "
printf, 1, "# time   "
printf, 1, "#(Gyr)      gas          disk        bulge     newstars           dm  "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3, 5(F11.6,"  "))', $
                time[i], mass_gas[7,i], mass_disk[7,i], mass_bulge[7,i], mass_newstars[7,i], mass_darkmatter[7,i]
endfor
close, 1
;
;
;   20 kpc
;  --------
openw, 1, frun+'/massin_20_aftermerger.txt', ERROR=err

printf, 1, "#   massin_20_aftermerger.txt"
printf, 1, "#      ---> mass of each component, in gadget units (10^10 msolar/h), within 20 kpc"
printf, 1, "#               of the galactic center as measured after the bh merger "
printf, 1, "# time   "
printf, 1, "#(Gyr)      gas          disk        bulge     newstars           dm  "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3, 5(F11.6,"  "))', $
                time[i], mass_gas[8,i], mass_disk[8,i], mass_bulge[8,i], mass_newstars[8,i], mass_darkmatter[8,i]
endfor
close, 1



            
            
        

; ----------------------------------------
        
print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"
        




end








;-----------------------------
; Process the amount of mass
;  within r
;-----------------------------
function process_mass_within_r, radius, mass


	mass_array= fltarr(9)

	mass_array[0]= total(mass)

	idx= where(radius le 1.0)
	if idx(0) ne -1 then mass_array[1]= total(mass(idx))

	idx= where(radius le 2.0)
	if idx(0) ne -1 then mass_array[2]= total(mass(idx))

	idx= where(radius le 3.0)
	if idx(0) ne -1 then mass_array[3]= total(mass(idx))

	idx= where(radius le 4.0)
	if idx(0) ne -1 then mass_array[4]= total(mass(idx))

	idx= where(radius le 5.0)
	if idx(0) ne -1 then mass_array[5]= total(mass(idx))

	idx= where(radius le 10.0)
	if idx(0) ne -1 then mass_array[6]= total(mass(idx))

	idx= where(radius le 15.0)
	if idx(0) ne -1 then mass_array[7]= total(mass(idx))

	idx= where(radius le 20.0)
	if idx(0) ne -1 then mass_array[8]= total(mass(idx))


	return, mass_array

end






; ----------------------------
;  Read hotgas.txt file
; ----------------------------
pro read_massinr_file, filename, time, mass_g, mass_d, mass_b, mass_ns, mass_dm

thisfile= filename

spawn, "wc "+thisfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then hgas_data= fltarr(6,lines)

openr, 1, thisfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, hgas_data
close, 1


time= hgas_data[0,*]
mass_g= hgas_data[1,*]
mass_d= hgas_data[2,*]
mass_b= hgas_data[3,*]
mass_ns= hgas_data[4,*]
mass_dm= hgas_data[5,*]


end






;===================================================================================





;--------------------------------------
;  Plot the mass within R
;----------------------------------------
pro plot_mass_within_r, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_mass_within_r, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='massinr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


yaxistitle = "Mass (h!E-1!N10!E10!N)"
xaxistitle = "Time After Merger (h!E-1!NGyr)"


xmax = 2.0
xmin = 0



ymax = 25.0
;ymax = 5.0
;ymin = 0.0
ymin= 0.1
;ymin= 0.001


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




; Data Set
;-----------
;frun="/raid4/tcox/vc3vc3e"
frun="/raid4/tcox/vc3vc3e_no"
;frun="/raid4/tcox/vc3vc3e_2"

; ---------------
read_massinr_file, frun+'/massin_0_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
nsnaps= n_elements(time)
mass_gas= fltarr(9,nsnaps)
mass_disk= fltarr(9,nsnaps)
mass_bulge= fltarr(9,nsnaps)
mass_newstars= fltarr(9,nsnaps)
mass_darkm= fltarr(9,nsnaps)
ti= 0
mass_gas[ti,*]= mass_g & mass_disk[ti,*]=mass_d & mass_bulge[ti,*]= mass_b & mass_newstars[ti,*]= mass_ns & mass_darkm[ti,*]= mass_dm
; ---------------
read_massinr_file, frun+'/massin_1_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
ti= 1
mass_gas[ti,*]= mass_g & mass_disk[ti,*]=mass_d & mass_bulge[ti,*]= mass_b & mass_newstars[ti,*]= mass_ns & mass_darkm[ti,*]= mass_dm
; ---------------
read_massinr_file, frun+'/massin_2_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
ti= 2
mass_gas[ti,*]= mass_g & mass_disk[ti,*]=mass_d & mass_bulge[ti,*]= mass_b & mass_newstars[ti,*]= mass_ns & mass_darkm[ti,*]= mass_dm
; ---------------
read_massinr_file, frun+'/massin_3_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
ti= 3
mass_gas[ti,*]= mass_g & mass_disk[ti,*]=mass_d & mass_bulge[ti,*]= mass_b & mass_newstars[ti,*]= mass_ns & mass_darkm[ti,*]= mass_dm
; ---------------
read_massinr_file, frun+'/massin_4_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
ti= 4
mass_gas[ti,*]= mass_g & mass_disk[ti,*]=mass_d & mass_bulge[ti,*]= mass_b & mass_newstars[ti,*]= mass_ns & mass_darkm[ti,*]= mass_dm
; ---------------
read_massinr_file, frun+'/massin_5_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
ti= 5
mass_gas[ti,*]= mass_g & mass_disk[ti,*]=mass_d & mass_bulge[ti,*]= mass_b & mass_newstars[ti,*]= mass_ns & mass_darkm[ti,*]= mass_dm
; ---------------
read_massinr_file, frun+'/massin_10_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
ti= 6
mass_gas[ti,*]= mass_g & mass_disk[ti,*]=mass_d & mass_bulge[ti,*]= mass_b & mass_newstars[ti,*]= mass_ns & mass_darkm[ti,*]= mass_dm
; ---------------
read_massinr_file, frun+'/massin_15_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
ti= 7
mass_gas[ti,*]= mass_g & mass_disk[ti,*]=mass_d & mass_bulge[ti,*]= mass_b & mass_newstars[ti,*]= mass_ns & mass_darkm[ti,*]= mass_dm
; ---------------
read_massinr_file, frun+'/massin_20_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
ti= 8
mass_gas[ti,*]= mass_g & mass_disk[ti,*]=mass_d & mass_bulge[ti,*]= mass_b & mass_newstars[ti,*]= mass_ns & mass_darkm[ti,*]= mass_dm

;----------------------------------------

;massam= mass_gas
;massam= mass_disk
;massam= mass_bulge
massam= mass_darkm
;massam= mass_newstars


idx= where(mass_darkm[0,*] gt 0.0)
mergidx= idx[0]

;t_merg= bh_mergertime(frun)    ; this needs bh.pro file to be compiled
t_merg= 1.15684  &  print, "WARNING: manually setting the merger time"
print, "T_merg= ",t_merg

time= time-t_merg


oplot, time[mergidx:*], massam[1,mergidx:*], thick=3.0, psym=-3, color= 200
oplot, time[mergidx:*], massam[2,mergidx:*], thick=3.0, psym=-3, color= 180
oplot, time[mergidx:*], massam[3,mergidx:*], thick=3.0, psym=-3, color= 160
oplot, time[mergidx:*], massam[4,mergidx:*], thick=3.0, psym=-3, color= 140
oplot, time[mergidx:*], massam[5,mergidx:*], thick=3.0, psym=-3, color= 120
oplot, time[mergidx:*], massam[6,mergidx:*], thick=3.0, psym=-3, color= 100
oplot, time[mergidx:*], massam[7,mergidx:*], thick=3.0, psym=-3, color= 80
oplot, time[mergidx:*], massam[8,mergidx:*], thick=3.0, psym=-3, color= 60


;  plot extras
; -------------
xyouts, 0.7, 0.85, fload_getid(frun), /normal, charthick=1, size=1.33, color=150

device, /close




end










;===================================================================================





;--------------------------------------
;  Plot the mass diff within R
;----------------------------------------
pro plot_massdiff_within_r, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_massdiff_within_r, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='dmassinr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


yaxistitle = "!4D!3Mass (h!E-1!N10!E10!N)"
xaxistitle = "Time After Merger (h!E-1!NGyr)"


xmax = 2.0
xmin = 0



ymax = 1.0
ymin= -1.0


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




; Data Set
;-----------
;frun="/raid4/tcox/vc3vc3e"
;frun="/raid4/tcox/vc3vc3e_no"
frun="/raid4/tcox/vc3vc3e_2"

; ---------------
;read_massinr_file, frun+'/massin_0_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
;read_massinr_file, frun+'/massin_1_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
read_massinr_file, frun+'/massin_2_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
;read_massinr_file, frun+'/massin_3_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
;read_massinr_file, frun+'/massin_4_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
;read_massinr_file, frun+'/massin_5_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
;read_massinr_file, frun+'/massin_10_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
;read_massinr_file, frun+'/massin_15_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm
;read_massinr_file, frun+'/massin_20_aftermerger.txt', time, mass_g, mass_d, mass_b, mass_ns, mass_dm

;----------------------------------------

; get non-zero values
idx= where(mass_dm gt 0.0)
mergidx= idx[0]

dtime= time[mergidx:*]


dmass_g= mass_g[mergidx:*]
dmass_d= mass_d[mergidx:*]
;dmass_b= mass_b[mergidx:*]
dmass_ns= mass_ns[mergidx:*]
dmass_dm= mass_dm[mergidx:*]

dmass_gns= mass_g[mergidx:*]+mass_ns[mergidx:*]


;subtract initial value
dmass_g= dmass_g - mass_g[mergidx]
dmass_d= dmass_d - mass_d[mergidx]
;dmass_b= dmass_b - mass_b[mergidx]
dmass_ns= dmass_ns - mass_ns[mergidx]
dmass_dm= dmass_dm - mass_dm[mergidx]
dmass_gns= dmass_gns - (mass_g[mergidx]+mass_ns[mergidx])


;----------------------------------------

;t_merg= bh_mergertime(frun)    ; this needs bh.pro file to be compiled
t_merg= 1.15684  &  print, "WARNING: manually setting the merger time"
print, "T_merg= ",t_merg

dtime= dtime-t_merg


oplot, dtime, dmass_g, thick=4.0, psym=-3, color= 150
oplot, dtime, dmass_d, thick=3.0, psym=-3, color= 50
;oplot, dtime, dmass_b, thick=3.0, psym=-3, color= 160
oplot, dtime, dmass_ns, thick=4.0, psym=-3, color= 100
oplot, dtime, dmass_dm, thick=3.0, psym=-3, color= 0
oplot, dtime, dmass_gns, thick=3.0, psym=-4, color= 220, linestyle=2


;  plot extras
; -------------
xyouts, 0.7, 0.85, fload_getid(frun), /normal, charthick=1, size=1.33, color=150

device, /close




end





