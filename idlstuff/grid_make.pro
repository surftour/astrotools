;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   Calculated grid quantities
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------





pro do_desika, junk

   ;frun="/raid4/tcox/As/A3nopot"
   ;starti= 90
   ;endi= 507

   ;frun="/raid4/tcox/vc3vc3e_2"
   ;frun="/raid4/tcox/vc3vc3e_no"
   frun="/raid4/tcox/vc3vc3h_2"
   ;frun="/raid4/tcox/sbw/sb10"
   ;frun="/raid4/tcox/sbw/sb10BH"
   ;frun="/raid4/tcox/bs/b3e_2"
   ;frun="/raid4/tcox/ds/d5e2_2"

   starti= 22
   ;starti= 32
   ;starti= 39
   ;starti= 40
   ;starti= 61
   endi= 31
   ;endi= 72
   ;endi= 99
   ;endi= 129

   ;desika_grid, frun, 22, /use_calc_center, /stellar_ages
   ;desika_grid, frun, 82, /use_calc_center, /stellar_ages
   ;desika_grid, frun, 87, /use_calc_center, /stellar_ages
   ;desika_grid, frun, 97, /use_calc_center, /stellar_ages
   for i=starti,endi do begin
	desika_grid, frun, i, /use_calc_center, /stellar_ages
   endfor

end




; =========================
;
;      Desika Grid
;
; =========================
pro desika_grid, frun, snapnum, center=center, $
				use_calc_center=use_calc_center, $
				stellar_ages=stellar_ages


if not keyword_set(frun) then begin
	print, "  "
	print, " desika_grid, frun, snapnum"
	print, "  "
	return
endif

if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(stellar_ages) then stellar_ages= 0
if not keyword_set(center) then center=[0,0,0]


;-------------------------------------------
;
;   Define Grid
;
;n_side= 25L
n_side= 50L
;n_side= 75L
;n_side= 100L
;n_side= 500L
gridlbl= strcompress(string(n_side),/remove_all)
griddim= gridlbl+'x'+gridlbl+'x'+gridlbl
;side_len= 70.0
;side_len= 12.0
side_len= 8.0
sidellbl= strcompress(string(side_len),/remove_all)
sidellbl= ' box len='+sidellbl
element_len= side_len/n_side
cellVolume= element_len*element_len*element_len

print, "Total grid volume= ", side_len*side_len*side_len
print, "cellVolume= ", cellVolume

Ngrid= n_side*n_side*n_side

Grid= fltarr(3,Ngrid)

n_i= long(0)
for i=0,n_side-1 do begin
  for j=0,n_side-1 do begin
    for k=0,n_side-1 do begin
	Grid[0,n_i]= (-1.0*side_len/2.0) + 0.5*element_len + i*element_len
	Grid[1,n_i]= (-1.0*side_len/2.0) + 0.5*element_len + j*element_len
	Grid[2,n_i]= (-1.0*side_len/2.0) + 0.5*element_len + k*element_len
	n_i= n_i + 1
    endfor
  endfor
endfor



;-------------------------------------------
;
;   Open Snapshot
;

;spawn, "/bin/ls "+frun+"/snap* | wc ",result                                    
;lastsnapshot=long(result[0])-1                                                  
;ok=fload_snapshot_bh(frun,lastsnapshot)         
;ok=fload_snapshot_bh(frun,snapnum)         
ok=fload_snapshot_bh(frun,snapnum, /nopot_in_snap)         


if keyword_set(use_calc_center) then center=fload_center_alreadycomp(1)

orig_center= center

;-------------------------------------------
;
;   Density
;

center= orig_center
x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
m= fload_gas_mass(1) / 0.7
ids= fload_gas_id(1)

N= long(n_elements(x))

print, "Total Gas Mass= ", total(m)
idx= where((abs(x) lt 0.5*side_len) and (abs(y) lt 0.5*side_len) and (abs(z) lt 0.5*side_len))
if idx(0) ne -1 then print, "Gas mass within cube= ", total(m(idx)) else print, "No mass in cube"


DesNgb=96L
Hmax= 100.0

if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= m(*)

        Density= fltarr(Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
                'density', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Density)

	; Make the Density Log
	; ---------------------
	print, " "
	print, "Cube Gas Mass= ", total(cellVolume*Density)

	; log and cm-3
	Density= alog10(674.59*Density)     ; convert to cm-3  (is mu right?)

	; this is probably (?) better
	;UnitDensity_in_cgs = 6.76991d-22
	;PROTONMASS = 1.6726d-24
	;factor= UnitDensity_in_cgs / PROTONMASS
	;print, "Density factor= ", factor
	;Density_TotalGas= alog10(factor*Density)     ; convert to cm-3  (is mu right?)


endif else begin
        print,'No stars, no problem.'
        return
endelse

help, Density



;-------------------------------------------
;
;   Temperature
;

temp_k= fload_gas_temperature_multi(1,/cold)

N= long(n_elements(x))

DesNgb=12L
Hmax= 100.0

if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	ColdTemps= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	ColdTemps(0:N-1)= temp_k(*)

        TempAvg= fltarr(Ngrid)

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
                'density', $
                N, $
                Ngrid, $
                Coord, $
                ColdTemps, $
                DesNgb, $
                Hmax, $
                TempAvg)

        ; Make the Density Log
        ; ---------------------
	print, " "
        Temp= alog10(TempAvg)


endif else begin
        print,'No gas, no problem.'
        return
endelse

help, Temp



;-------------------------------------------
;
;   Velocity
;

center= orig_center
comvel= fload_all_comvel(1,center=center)
vx= fload_gas_v('x',comvel=comvel,center=center)
vy= fload_gas_v('y',comvel=comvel,center=center)
vz= fload_gas_v('z',comvel=comvel,center=center)

N= long(n_elements(vx))


DesNgb=96L
Hmax= 100.0

if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(10,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= vx(*)
	Coord(4,0:N-1)= vy(*)
	Coord(5,0:N-1)= vz(*)
	Coord(6,0:N-1)= ids(*)
	Coord(7,0:Ngrid-1)= Grid(0,*)
	Coord(8,0:Ngrid-1)= Grid(1,*)
	Coord(9,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= m(*)

        Velocity= fltarr(3,Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Velocity/velocity.so', $
                'velocity', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Velocity)
	print, " "


endif else begin
        print,'No stars, no problem.'
        return
endelse

help, Velocity



;-------------------------------------------
;
;   Stellar Mass
;

center= orig_center
x= fload_allstars_xyz('x',center=center) / 0.7
y= fload_allstars_xyz('y',center=center) / 0.7
z= fload_allstars_xyz('z',center=center) / 0.7
m= fload_allstars_mass(1) / 0.7
ids= fload_allstars_ids(1)
age= (fload_time(1) - fload_allstars_age(1)) / 0.7

N= long(n_elements(x))

print, "Total Stellar Mass= ", total(m)
idx= where((abs(x) lt 0.5*side_len) and (abs(y) lt 0.5*side_len) and (abs(z) lt 0.5*side_len))
if idx(0) ne -1 then print, "Stellar mass within cube= ", total(m(idx)) else print, "No mass in cube"

DesNgb=96L
Hmax= 100.0


if(N gt 0) then begin

	if keyword_set(stellar_ages) then begin

		; Stellar Age 1
                ;-------------------
                idx= where(age le 0.010)           ; below 10 Myr
		help, idx
		;if idx(0) ne -1 then begin
		if n_elements(idx) gt 1 then begin
                   N= long(n_elements(idx))
                   Ncoord= max([N,Ngrid]) & Coord=fltarr(7,Ncoord) & Masses= fltarr(Ncoord)

                   Coord(0,0:N-1)= x(idx) & Coord(1,0:N-1)= y(idx) & Coord(2,0:N-1)= z(idx) & Coord(3,0:N-1)= ids(idx)
                   Coord(4,0:Ngrid-1)= Grid(0,*) & Coord(5,0:Ngrid-1)= Grid(1,*) & Coord(6,0:Ngrid-1)= Grid(2,*)
                   Masses(0:N-1)= m(idx)

                   StellarDensity= fltarr(Ngrid)      ;  full velocity field

                   S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
                        'density', $
                        N, $
                        Ngrid, $
                        Coord, $
                        Masses, $
                        DesNgb, $
                        Hmax, $
                        StellarDensity)

		   print, " "
                   print, "Total Stellar mass (according to grid)= ", total(cellVolume*StellarDensity)
                   StellarMass_1= 10.0+alog10(cellVolume*StellarDensity)
		   help, StellarMass_1
		endif else begin
		   StellarMass_1= fltarr(Ngrid)
		   StellarMass_1(*)= 0.0
		endelse


		; Stellar Age 2
                ;-------------------
                idx= where((age gt 0.010) and (age le 0.100))
		help, idx
		;if idx(0) ne -1 then begin
		if n_elements(idx) gt 1 then begin
                   N= long(n_elements(idx))
                   Ncoord= max([N,Ngrid]) & Coord=fltarr(7,Ncoord) & Masses= fltarr(Ncoord)

                   Coord(0,0:N-1)= x(idx) & Coord(1,0:N-1)= y(idx) & Coord(2,0:N-1)= z(idx) & Coord(3,0:N-1)= ids(idx)
                   Coord(4,0:Ngrid-1)= Grid(0,*) & Coord(5,0:Ngrid-1)= Grid(1,*) & Coord(6,0:Ngrid-1)= Grid(2,*)
                   Masses(0:N-1)= m(idx)

                   StellarDensity= fltarr(Ngrid)      ;  full velocity field

                   S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
                        'density', $
                        N, $
                        Ngrid, $
                        Coord, $
                        Masses, $
                        DesNgb, $
                        Hmax, $
                        StellarDensity)

		   print, " "
                   print, "Total Stellar mass (according to grid)= ", total(cellVolume*StellarDensity)
                   StellarMass_2= 10.0+alog10(cellVolume*StellarDensity)
		   help, StellarMass_2
		endif else begin
		   StellarMass_2= fltarr(Ngrid)
		   StellarMass_2(*)= 0.0
		endelse


		; Stellar Age 3
		;-------------------
		idx= where((age gt 0.100) and (age le 1.0))
		help, idx
		;if idx(0) ne -1 then begin
		if n_elements(idx) gt 1 then begin
		   N= long(n_elements(idx))
                   Ncoord= max([N,Ngrid]) & Coord=fltarr(7,Ncoord) & Masses= fltarr(Ncoord)

                   Coord(0,0:N-1)= x(idx) & Coord(1,0:N-1)= y(idx) & Coord(2,0:N-1)= z(idx) & Coord(3,0:N-1)= ids(idx)
                   Coord(4,0:Ngrid-1)= Grid(0,*) & Coord(5,0:Ngrid-1)= Grid(1,*) & Coord(6,0:Ngrid-1)= Grid(2,*)
                   Masses(0:N-1)= m(idx)

                   StellarDensity= fltarr(Ngrid)      ;  full velocity field

                   S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
                        'density', $
                        N, $
                        Ngrid, $
                        Coord, $
                        Masses, $
                        DesNgb, $
                        Hmax, $
                        StellarDensity)

		   print, " "
                   print, "Total Stellar mass (according to grid)= ", total(cellVolume*StellarDensity)
                   StellarMass_3= 10.0+alog10(cellVolume*StellarDensity)
		   help, StellarMass_3
		endif else begin
		   StellarMass_3= fltarr(Ngrid)
		   StellarMass_3(*)= 0.0
		endelse


                ; Stellar Age 4
                ;-------------------
                idx= where(age gt 1.0)
		help, idx
		;if idx(0) ne -1 then begin
		if n_elements(idx) gt 1 then begin
                   N= long(n_elements(idx))
                   Ncoord= max([N,Ngrid]) & Coord=fltarr(7,Ncoord) & Masses= fltarr(Ncoord)

                   Coord(0,0:N-1)= x(idx) & Coord(1,0:N-1)= y(idx) & Coord(2,0:N-1)= z(idx) & Coord(3,0:N-1)= ids(idx)
                   Coord(4,0:Ngrid-1)= Grid(0,*) & Coord(5,0:Ngrid-1)= Grid(1,*) & Coord(6,0:Ngrid-1)= Grid(2,*)
                   Masses(0:N-1)= m(idx)

                   StellarDensity= fltarr(Ngrid)      ;  full velocity field

                   S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
                        'density', $
                        N, $
                        Ngrid, $
                        Coord, $
                        Masses, $
                        DesNgb, $
                        Hmax, $
                        StellarDensity)

		   print, " "
                   print, "Total Stellar mass (according to grid)= ", total(cellVolume*StellarDensity)
                   StellarMass_4= 10.0+alog10(cellVolume*StellarDensity)
		   help, StellarMass_4
		endif else begin
		   StellarMass_4= fltarr(Ngrid)
		   StellarMass_4(*)= 0.0
		endelse

	endif else begin
		Ncoord= max([N,Ngrid]) & Coord=fltarr(7,Ncoord) & Masses= fltarr(Ncoord)

		Coord(0,0:N-1)= x(*) & Coord(1,0:N-1)= y(*) & Coord(2,0:N-1)= z(*) & Coord(3,0:N-1)= ids(*)
		Coord(4,0:Ngrid-1)= Grid(0,*) & Coord(5,0:Ngrid-1)= Grid(1,*) & Coord(6,0:Ngrid-1)= Grid(2,*) 
		Masses(0:N-1)= m(*)

	        StellarDensity= fltarr(Ngrid)      ;  full velocity field

	        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
	                'density', $
			N, $
			Ngrid, $
	                Coord, $
	                Masses, $
			DesNgb, $
			Hmax, $
	                StellarDensity)

		; Make the Density Log
		; ---------------------
		print, " "
		print, "Total Stellar mass (according to grid)= ", total(cellVolume*StellarDensity)
		StellarMass= 10.0+alog10(cellVolume*StellarDensity)
		help, StellarMass
	endelse

endif else begin
        print,'No stars, no problem.'
        return
endelse





idx=where(finite(Grid) eq 0, n_i)
if idx(0) ne -1 then print, 'Grid has '+strcompress(string(n_i),/remove_all)+' bad numbers'

idx=where(finite(Velocity) eq 0, n_i)
if idx(0) ne -1 then print, 'Velocity has '+strcompress(string(n_i),/remove_all)+' bad numbers'

idx=where(finite(Density) eq 0, n_i)
if idx(0) ne -1 then print, 'Density has '+strcompress(string(n_i),/remove_all)+' bad numbers'

idx=where(finite(Temp) eq 0, n_i)
if idx(0) ne -1 then print, 'Temp has '+strcompress(string(n_i),/remove_all)+' bad numbers'

idx=where(finite(StellarMass_1) eq 0, n_i)
if idx(0) ne -1 then print, 'StellarMass_1 has '+strcompress(string(n_i),/remove_all)+' bad numbers'
idx=where(finite(StellarMass_2) eq 0, n_i)
if idx(0) ne -1 then print, 'StellarMass_2 has '+strcompress(string(n_i),/remove_all)+' bad numbers'
idx=where(finite(StellarMass_3) eq 0, n_i)
if idx(0) ne -1 then print, 'StellarMass_3 has '+strcompress(string(n_i),/remove_all)+' bad numbers'
idx=where(finite(StellarMass_4) eq 0, n_i)
if idx(0) ne -1 then print, 'StellarMass_4 has '+strcompress(string(n_i),/remove_all)+' bad numbers'


;=========================================================================================
;----------------------------------------------
;
;   Write to Grid file
;

;gridfile=frun+'/desika_grids/griddata_'+strcompress(string(snapnum),/remove_all)+'.txt'
gridfile=frun+'/desikagrids/griddata_'+strcompress(string(snapnum),/remove_all)+'.txt'
openw, 1, gridfile, ERROR=err

; actually print the rows
if keyword_set(stellar_ages) then begin
	printf, 1, '# '+griddim+sidellbl+' frun='+frun+'  snapnum=',strcompress(string(snapnum),/remove_all)
	printf, 1, "# "
	printf, 1, "#                                                 Cold Gas Cold Gas              Stellar Mass (Msolar)"
	printf, 1, "#   x     y      z       Velx      Vely      Velz  Density     Temp      Age   >10Myr  >100Myr"
	printf, 1, "#(kpc) (kpc)  (kpc)    (km/s)    (km/s)    (km/s)   (cm-3)      (K)  0-10Myr <=100Myr   <=1Gyr    >1Gyr"

	for i=0L,Ngrid-1 do begin
		printf, 1, FORMAT='(3(F5.2,"  "),3(F8.2,"  "),6(F7.3,"  "))', $
			Grid[0,i],Grid[1,i],Grid[2,i], $
			Velocity[0,i],Velocity[1,i],Velocity[2,i], $
			Density[i], $
			Temp[i], $
			StellarMass_1[i], $
			StellarMass_2[i], $
			StellarMass_3[i], $
			StellarMass_4[i]
	endfor
endif else begin
	printf, 1, '# '+griddim+sidellbl+' frun='+frun+'  snapnum=',strcompress(string(snapnum),/remove_all)
	printf, 1, "# "
	printf, 1, "#                                                 Cold Gas Cold Gas Stellar"
	printf, 1, "#   x     y      z       Velx      Vely      Velz  Density     Temp    Mass"
	printf, 1, "#(kpc) (kpc)  (kpc)    (km/s)    (km/s)    (km/s)   (cm-3)      (K) (Msolar)"

	for i=0L,Ngrid-1 do begin
		printf, 1, FORMAT='(3(F5.2,"  "),3(F8.2,"  "),3(F7.3,"  "))', $
			Grid[0,i],Grid[1,i],Grid[2,i], $
			Velocity[0,i],Velocity[1,i],Velocity[2,i], $
			Density[i], $
			Temp[i], $
			StellarMass[i]
	endfor
endelse
close, 1

;----------------------------------------------

end











;==========================================================================================





; =========================
;
;      Simple Grid
;
; =========================
pro simple_grid, frun, snapnum, center=center, $
				use_calc_center=use_calc_center, $
				stellar_ages=stellar_ages


if not keyword_set(frun) then begin
	print, "  "
	print, " simple_grid, frun, snapnum"
	print, "  "
	return
endif

if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(stellar_ages) then stellar_ages= 0
if not keyword_set(center) then center=[0,0,0]


;-------------------------------------------
;
;   Define Grid
;
;n_side= 25L
;n_side= 50L
;n_side= 75L
;n_side= 100L
n_side= 200L
gridlbl= strcompress(string(n_side),/remove_all)
griddim= gridlbl+'x'+gridlbl+'x'+gridlbl
side_len= 70.0
;side_len= 12.0
;side_len= 8.0
sidellbl= strcompress(string(side_len),/remove_all)
sidellbl= ' box len='+sidellbl
element_len= side_len/n_side
cellVolume= element_len*element_len*element_len

print, "Total grid volume= ", side_len*side_len*side_len
print, "cellVolume= ", cellVolume

Ngrid= n_side*n_side*n_side

Grid= fltarr(3,Ngrid)

n_i= long(0)
for i=0,n_side-1 do begin
  for j=0,n_side-1 do begin
    for k=0,n_side-1 do begin
	Grid[0,n_i]= (-1.0*side_len/2.0) + 0.5*element_len + i*element_len
	Grid[1,n_i]= (-1.0*side_len/2.0) + 0.5*element_len + j*element_len
	Grid[2,n_i]= (-1.0*side_len/2.0) + 0.5*element_len + k*element_len
	n_i= n_i + 1
    endfor
  endfor
endfor



;-------------------------------------------
;
;   Open Snapshot
;

;spawn, "/bin/ls "+frun+"/snap* | wc ",result                                    
;lastsnapshot=long(result[0])-1                                                  
;ok=fload_snapshot_bh(frun,lastsnapshot)         
;ok=fload_snapshot_bh(frun,snapnum)         
ok=fload_snapshot_bh(frun,snapnum, /nopot_in_snap)         


if keyword_set(use_calc_center) then center=fload_center_alreadycomp(1)

orig_center= center

;-------------------------------------------
;
;   Density
;

center= orig_center
x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
m= fload_gas_mass(1) / 0.7
ids= fload_gas_id(1)

N= long(n_elements(x))

print, "Total Gas Mass= ", total(m)
idx= where((abs(x) lt 0.5*side_len) and (abs(y) lt 0.5*side_len) and (abs(z) lt 0.5*side_len))
if idx(0) ne -1 then print, "Gas mass within cube= ", total(m(idx)) else print, "No mass in cube"


DesNgb=96L
Hmax= 100.0

if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= m(*)

        Density= fltarr(Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
                'density', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Density)

	; Make the Density Log
	; ---------------------
	print, " "
	print, "Cube Gas Mass= ", total(cellVolume*Density)

	; log and cm-3
	Density= alog10(674.59*Density)     ; convert to cm-3  (is mu right?)

	; this is probably (?) better
	;UnitDensity_in_cgs = 6.76991d-22
	;PROTONMASS = 1.6726d-24
	;factor= UnitDensity_in_cgs / PROTONMASS
	;print, "Density factor= ", factor
	;Density_TotalGas= alog10(factor*Density)     ; convert to cm-3  (is mu right?)


endif else begin
        print,'No gas, no problem.'
        return
endelse

help, Density




idx=where(finite(Grid) eq 0, n_i)
if idx(0) ne -1 then print, 'Grid has '+strcompress(string(n_i),/remove_all)+' bad numbers'

idx=where(finite(Density) eq 0, n_i)
if idx(0) ne -1 then print, 'Density has '+strcompress(string(n_i),/remove_all)+' bad numbers'



;=========================================================================================
;----------------------------------------------
;
;   Write to Grid file
;

gridfile=frun+'/simplegriddata_'+strcompress(string(snapnum),/remove_all)+'.txt'
openw, 1, gridfile, ERROR=err

; actually print the rows
	printf, 1, '# '+griddim+sidellbl+' frun='+frun+'  snapnum=',strcompress(string(snapnum),/remove_all)
	printf, 1, "# "
	printf, 1, "# "
	printf, 1, "#   x     y      z    Density"
	printf, 1, "#(kpc) (kpc)  (kpc)    (cm-3)"

	for i=0L,Ngrid-1 do begin
		printf, 1, FORMAT='(3(F5.2,"  "),1(F8.2,"  "))', $
			Grid[0,i],Grid[1,i],Grid[2,i], $
			Density[i]
	endfor

close, 1

;----------------------------------------------

end









;==========================================================================================








; =========================
;
;      Grid 2
;
; =========================




;  grid provided to sukanya
; ------------------------------------
pro sukanya_grid, frun, snapnum, extra=extra


if not keyword_set(frun) then begin
	print, "  "
	print, " sukanya_grid, frun, snapnum"
	print, "  "
	return
endif

;frun= "/raid4/tcox/vc3vc3e"
if not keyword_set(snapnum) then snapnum= 12
;snapnum= 31
;ok=fload_snapshot_bh(frun,snapnum)

if not keyword_set(extra) then extra=''



;-------------------------------------------
;
;   Define Grid
;  --------------

;  spherical
; ------------
spherical_grid= 1
;spherical_grid= 0
if spherical_grid eq 1 then begin

	; need the following order

	; grid_*
	; ----------
	N_r= 200L     &  grid_r= fltarr(N_r)
	N_theta= 59L  &  grid_theta= fltarr(N_theta)
	N_phi= 118L   &  grid_phi= fltarr(N_phi)
	; cent_*
	; ---------
	N_r= 199L     &  grid_cen_r= fltarr(N_r)
	N_theta= 58L  &  grid_cen_theta= fltarr(N_theta)
	N_phi= 117L   &  grid_cen_phi= fltarr(N_phi)
	; test_*
	; ---------
	;N_r= 23L     &  r_data= fltarr(N_r)
	;N_theta= 9L  &  theta_data= fltarr(N_theta)
	;N_phi= 25L   &  phi_data= fltarr(N_phi)


	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_r.dat"
	readf, 1, grid_r
	close, 1
	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_r.dat"
	readf, 1, grid_cen_r
	close, 1

	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_theta.dat"
	readf, 1, grid_theta
	close, 1
	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_theta.dat"
	readf, 1, grid_cen_theta
	close, 1

	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_phi.dat"
	readf, 1, grid_phi
	close, 1
	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_phi.dat"
	readf, 1, grid_cen_phi
	close, 1


	r_len= max(grid_cen_r)
	total_Volume= (4./3.)*!PI*r_len*r_len*r_len
	print, "total_Volume= ", total_Volume

	Ngrid= N_r*N_theta*N_phi

	Grid= fltarr(3,Ngrid)
	Grid_rthetaphi= fltarr(3,Ngrid)
	Grid_Volume= fltarr(Ngrid)

	n_i= long(0)
	for i=0,N_r-1 do begin
	  for j=0,N_theta-1 do begin
	    for k=0,N_phi-1 do begin

		R= grid_cen_r[i]
		theta= grid_cen_theta[j]
		phi= grid_cen_phi[k]
		Grid_rthetaphi[0,n_i]= R
		Grid_rthetaphi[1,n_i]= theta
		Grid_rthetaphi[2,n_i]= phi

		x= R*sin(theta)*cos(phi)
		y= R*sin(theta)*sin(phi)
		z= R*cos(theta)
		Grid[0,n_i]= x
		Grid[1,n_i]= y
		Grid[2,n_i]= z

		;Grid_Volume[n_i]= R * 10^(d_log_r) * d_phi * d_cos_theta
		;Grid_Volume[n_i]= R * 10^(d_log_r) * d_phi * sin(theta) * d_theta
		;Grid_Volume[n_i]= R * R * d_log_r * d_phi * sin(theta) * d_theta    ;  dr= r * d_logr
		delta_r= grid_r[i+1] - grid_r[i]
		delta_theta= grid_theta[j+1] - grid_theta[j]
		delta_phi= grid_phi[k+1] - grid_phi[k]
		Grid_Volume[n_i]= R * delta_r * delta_phi * sin(theta) * delta_theta

		n_i= n_i + 1
	    endfor
	  endfor
	endfor
endif



;  Volume check
; --------------
sphere_area= 0.0
for j=0,N_theta-1 do begin
  for k=0,N_phi-1 do begin
        delta_theta= grid_theta[j+1] - grid_theta[j]
        delta_phi= grid_phi[k+1] - grid_phi[k]
        sphere_area= sphere_area + delta_phi * sin(grid_theta[j]) * delta_theta
  endfor
endfor
print, "sphere area= ", sphere_area
print, '4*PI= ', 4.0*!PI




cube_grid= 0
;cube_grid= 1
if cube_grid eq 1 then begin
	n_side= 50L
	side_len= 12.0
	element_len= side_len/n_side
	cellVolume= element_len*element_len*element_len

	Ngrid= n_side*n_side*n_side

	Grid= fltarr(3,Ngrid)
	Grid_rthetaphi= fltarr(3,Ngrid)

	n_i= long(0)
	for i=0,n_side-1 do begin
	  for j=0,n_side-1 do begin
	    for k=0,n_side-1 do begin
		Grid[0,n_i]= (-1.0*side_len/2.0) + 0.5*element_len + i*element_len
		Grid[1,n_i]= (-1.0*side_len/2.0) + 0.5*element_len + j*element_len
		Grid[2,n_i]= (-1.0*side_len/2.0) + 0.5*element_len + k*element_len
		n_i= n_i + 1
	    endfor
	  endfor
	endfor
endif


;-------------------------------------------
;
;   Open Snapshot
;

;spawn, "/bin/ls "+frun+"/snap* | wc ",result                                    
;lastsnapshot=long(result[0])-1                                                  
;ok=fload_snapshot_bh(frun,lastsnapshot)         
ok=fload_snapshot_bh(frun,snapnum)         


;-------------------------------------------
;
;   Density - Total Gas
;

print, " ---------------------"
print, "  Total Gas Density  "
print, " ---------------------"

; determine_center
; -----------------
center= [0,0,0]
center_bh= fload_center_alreadycomp(1)

; two black holes
if fload_npart(5) gt 1 then begin
	bhid= fload_blackhole_id(1)
	bhid1= bhid[0]
	bhid2= bhid[1]
	print, "Blackhole ID's: ", bhid1, bhid2
	center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
	center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
	center_bh= center1
	;center_bh= center2
endif

; one black hole
if fload_npart(5) eq 1 then begin
	bhid= fload_blackhole_id(1)
	center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
endif


;  grab variables
; ----------------
print, "using center= ", center_bh
x= fload_gas_xyz('x',center=center_bh)
y= fload_gas_xyz('y',center=center_bh)
z= fload_gas_xyz('z',center=center_bh)
m= fload_gas_mass(1)
ids= fload_gas_id(1)

print, "Total Gas Mass= ", total(m)
rr= sqrt(x*x + y*y + z*z)
if spherical_grid eq 1 then idx= where(rr lt r_len)
if cube_grid eq 1 then begin
	idx= where((abs(x) lt 0.5*side_len) and (abs(y) lt 0.5*side_len) and (abs(z) lt 0.5*side_len))
endif
if idx(0) ne -1 then print, "Gas mass within Sphere = ", total(m(idx)) else print, "No Mass in Sphere"

; ------
rho= fload_gas_rho(1)
idx= where(rho gt 0.000854924)
print, "SF Gas Mass= ", total(m(idx))


; -----------------------------------------------
print, "================================="
print, "      (Total) Gas Density "
print, "================================="
print, "*** std density routine ***"

N= long(n_elements(x))


DesNgb=32L
Hmax= 100.0

if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= m(*)

        Density_TotalGas= fltarr(Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
                'density', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Density_TotalGas)

	; Make the Density Log
	; ---------------------
	print, "Cube Gas Mass= ", total(Grid_Volume*Density_TotalGas)

	UnitDensity_in_cgs = 6.76991d-22
	PROTONMASS = 1.6726d-24
	factor= UnitDensity_in_cgs / PROTONMASS
	print, "Density factor= ", factor
	Density_TotalGas= alog10(factor*Density_TotalGas)     ; convert to cm-3  (is mu right?)
	;Density_TotalGas= alog10(674.59*Density_TotalGas)     ; convert to cm-3  (is mu right?)


endif else begin
        print,'No stars, no problem.'
        return
endelse

help, Density_TotalGas


;-------------------------------------------
;
;  test Average Routine



	print, "================================="
	print, "      (Total) Gas Density "
	print, "================================="
	print, "*** test the average routine ***"
	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)
	RhoGas= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= m(*)
	RhoGas(0:N-1)= rho(*)

        Density_Gas= fltarr(Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Average/average.so', $
                'average', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Density_Gas, $
		RhoGas)

	print, "Cube Gas Mass (via Average routine)= ", total(Grid_Volume*Density_Gas)
	print, "Density_Gas: max/min= ", max(Density_Gas), min(Density_Gas)




;-------------------------------------------
;
;   Density - Cold Gas
;
;    (via the density-weighted average method.
;     does that make sense?)
;

print, " ---------------------"
print, "  Dense Cold Gas  "
print, " ---------------------"

x= fload_gas_xyz('x',center=center_bh)
y= fload_gas_xyz('y',center=center_bh)
z= fload_gas_xyz('z',center=center_bh)
m= fload_gas_mass(1)
ids= fload_gas_id(1)
rho= fload_gas_rho(1)
print, "rho min/max= ", min(rho), max(rho)
rho_cold= fload_gas_rho(1,/cold)
print, "rho_cold max/min: ", max(rho_cold), min(rho_cold)
idx= where(rho_cold gt 0.0)
print, "    n (greater than zero)= ", n_elements(idx)
print, "SF (or Multiphase) Gas Mass= ", total(m(idx))

coldf= fload_gas_coldfraction(1)
print, "coldf max/min: ", max(coldf), min(coldf)

ColdGasMass= m*coldf
print, "Cold Gas Mass= ", total(ColdGasMass)

rr= sqrt(x*x + y*y + z*z)
if spherical_grid eq 1 then idx= where(rr lt r_len)
if cube_grid eq 1 then begin
	idx= where((abs(x) lt 0.5*side_len) and (abs(y) lt 0.5*side_len) and (abs(z) lt 0.5*side_len))
endif
if idx(0) ne -1 then print, "Cold gas mass within Cube= ", total(m(idx)*coldf(idx)) else print, "No Mass in Cube"

N= long(n_elements(x))


DesNgb=32L
Hmax= 100.0

if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)
	RhoColdGas= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= ColdGasMass(*)
	RhoColdGas(0:N-1)= rho_cold(*)

        Density_ColdGas= fltarr(Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Average/average.so', $
                'average', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Density_ColdGas, $
		RhoColdGas)

	; Make the Density Log
	; ---------------------
	print, "Cube Cold Gas Mass= ", total(Grid_Volume*Density_ColdGas)

	print, "Density_ColdGas: max/min= ", max(Density_ColdGas), min(Density_ColdGas)
	idx=where(Density_ColdGas le 0.0)
	if idx(0) ne -1 then Density_ColdGas(idx)= 1.0e-8
	Density_ColdGas= alog10(factor*Density_ColdGas)     ; convert to cm-3  (is mu right?)
	;Density_ColdGas= alog10(674.59*Density_ColdGas)     ; convert to cm-3  (is mu right?)
	if idx(0) ne -1 then Density_ColdGas(idx)= 0.0


endif else begin
        print,'No stars, no problem.'
        return
endelse





;
;  test alternate version
; -------------------------------------------
print, "---------------------------------------------"
print, "*** test aternate version of cold density ***"
print, "---------------------------------------------"
;DesNgb=5L
;Hmax= 0.050
DesNgb=2L
Hmax= 0.015

if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= ColdGasMass(*)

        Density_ColdGas2= fltarr(Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
                'density', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Density_ColdGas2)

	; Make the Density Log
	; ---------------------
	print, "Cube Gas Mass (II)= ", total(Grid_Volume*Density_ColdGas2)
	idx=where(Density_ColdGas2 le 0.0)
	if idx(0) ne -1 then Density_ColdGas2(idx)= 1.0e-8
	Density_ColdGas2= alog10(factor*Density_ColdGas2)     ; convert to cm-3  (is mu right?)
	if idx(0) ne -1 then Density_ColdGas2(idx)= 0.0


	; use this one!
	print, " ************************** "
	print, "   HEY: using the density   "
	print, "        version not the " 
	print, "        average one    "
	print, " ************************** "
	Density_ColdGas= Density_ColdGas2
	print, "Density_ColdGas= Density_ColdGas2"
	print, " "

	;extra= 'test2_'


endif else begin
        print,'No gas, no problem.'
        return
endelse




;-------------------------------------------
;
;   Cold Gas Filling Factor
;
;    (via the density-weighted average method.
;     does that make sense?)
;

print, " -------------------------"
print, "  Cold Gas Filling Factor "
print, " -------------------------"

coldgasff= fload_gas_coldgasfillingfactor(1)
print, "coldgasff max/min: ", max(coldgasff), min(coldgasff)


N= long(n_elements(x))


DesNgb=12L
Hmax= 1.0

if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)
	ColdGasFillingFac= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= m(*)
	ColdGasFillingFac(0:N-1)= coldgasff(*)

        ColdFillingFac= fltarr(Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Average/average.so', $
                'average', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                ColdFillingFac, $
		ColdGasFillingFac)

	; Make the Density Log
	; ---------------------
	;print, "cellVolume= ", cellVolume
	print, "Cube Cold Gas Mass (without filling factor considered)= ", total(Grid_Volume*Density_ColdGas)
	print, "Cube Cold Gas Mass (with filling factor considered)= ", total(Grid_Volume*Density_ColdGas*ColdFillingFac)

	print, "ColdFillingFac: max/min= ", max(ColdFillingFac), min(ColdFillingFac)
	idx=where(ColdFillingFac le 0.0)
	if idx(0) ne -1 then ColdFillingFac(idx)= 1.0e-8
	ColdFillingFac= alog10(ColdFillingFac)     ; convert to cm-3  (is mu right?)
	if idx(0) ne -1 then ColdFillingFac(idx)= 0.0


endif else begin
        print,'No stars, no problem.'
        return
endelse



;-------------------------------------------
;
;   Temperature
;

print, " ---------------------"
print, "  Gas Temperature  " 
print, " ---------------------"

temp_k= fload_gas_temperature_multi(1,/cold)
print, "Cold Temp max/min= ", max(temp_k), min(temp_k)

N= long(n_elements(x))

DesNgb=12L
Hmax= 1.0

if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)
	ColdTemps= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= m(*)
	ColdTemps(0:N-1)= temp_k(*)

        TempAvg= fltarr(Ngrid)

        ;S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
        ;        'density', $
        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Average/average.so', $
                'average', $
                N, $
                Ngrid, $
                Coord, $
                Masses, $
                DesNgb, $
                Hmax, $
                TempAvg, $
		ColdTemps)

        ; Make the Density Log
        ; ---------------------
        Temp= alog10(TempAvg)


endif else begin
        print,'No gas, no problem.'
        return
endelse

help, Temp



;-------------------------------------------
;
;   Velocity
;

;print, " ---------------------"
;print, "  Gas Velocity  "
;print, " ---------------------"

;comvel= fload_all_comvel(1)
;vx= fload_gas_v('x',comvel=comvel)
;vy= fload_gas_v('y',comvel=comvel)
;vz= fload_gas_v('z',comvel=comvel)

;N= long(n_elements(x))


;DesNgb=96L
;Hmax= 100.0

;if(N gt 0) then begin

	;Ncoord= max([N,Ngrid])

	;Coord=fltarr(10,Ncoord)
	;Masses= fltarr(Ncoord)

	;Coord(0,0:N-1)= x(*)
	;Coord(1,0:N-1)= y(*)
	;Coord(2,0:N-1)= z(*)
	;Coord(3,0:N-1)= vx(*)
	;Coord(4,0:N-1)= vy(*)
	;Coord(5,0:N-1)= vz(*)
	;Coord(6,0:N-1)= ids(*)
	;Coord(7,0:Ngrid-1)= Grid(0,*)
	;Coord(8,0:Ngrid-1)= Grid(1,*)
	;Coord(9,0:Ngrid-1)= Grid(2,*)
	;Masses(0:N-1)= m(*)

        ;Velocity= fltarr(3,Ngrid)      ;  full velocity field

        ;S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/Cube/Velocity/velocity.so', $
        ;        'velocity', $
	;	N, $
	;	Ngrid, $
        ;        Coord, $
        ;        Masses, $
	;	DesNgb, $
	;	Hmax, $
        ;        Velocity)


;endif else begin
;        print,'No stars, no problem.'
;        return
;endelse

;help, Velocity



;-------------------------------------------
;
;   Stellar Mass
;

print, " ---------------------"
print, "  Stellar Mass  "
print, " ---------------------"

x= fload_allstars_xyz('x',center=center_bh)
y= fload_allstars_xyz('y',center=center_bh)
z= fload_allstars_xyz('z',center=center_bh)
m= fload_allstars_mass(1)
ids= fload_allstars_ids(1)

print, "Total Stellar Mass= ", total(m)
rr= sqrt(x*x + y*y + z*z)
if spherical_grid eq 1 then idx= where(rr lt r_len)
if cube_grid eq 1 then begin
	idx= where((abs(x) lt 0.5*side_len) and (abs(y) lt 0.5*side_len) and (abs(z) lt 0.5*side_len))
endif
if idx(0) ne -1 then print, "Stellar Mass within Cube= ", total(m(idx)) else print, "No Mass in Cube"

N= long(n_elements(x))

DesNgb=96L
Hmax= 100.0


if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= m(*)

        StellarDensity= fltarr(Ngrid)      ;  full velocity field

        ;S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
        ;        'density', $
;		N, $
;		Ngrid, $
;                Coord, $
;                Masses, $
;		DesNgb, $
;		Hmax, $
;                StellarDensity)
	StellarDensity(*)= 1.0

	; Make the Density Log
	; ---------------------
	print, "Cube Stellar Mass= ", total(Grid_Volume*StellarDensity)
	StellarMass= 10.0+alog10(Grid_Volume*StellarDensity)



endif else begin
        print,'No stars, no problem.'
        return
endelse

help, StellarMass





;-------------------------------------------
;
;   Stellar Luminosity
;

print, " ---------------------"
print, "  Stellar Luminosity  "
print, " ---------------------"

    ndisk= fload_npart(2)                               ; disk particles
    nbulge= fload_npart(3)                              ; bulge particles
    nstars= fload_npart(4)
    npart= long(ndisk) + long(nbulge) + long(nstars)
    print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars
    TTime= float(fload_time(1))

    N= npart
    m= 1.0e+10*fload_allstars_mass(1)
    age=fload_allstars_age(1)
    age=float(TTime-age)
    zmets=fload_allstars_z(1)


    ; get the luminosities
    ;  - in units of solar luminosities
    print, "load luminosities"
    load_all_stellar_luminosities, N, TTime, m, age, zmets, $
                                Lum_Bol, Lum_U, Lum_B, Lum_V, Lum_R, Lum_I, Lum_J, Lum_H, Lum_K, $
                                Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z   ;, /notsolar

    ; trap for any NaN's
    idx=where(finite(Lum_Bol) eq 0)
    if idx(0) ne -1 then begin
        Lum_Bol(idx)= 100.0
        print, "Fixed ",n_elements(idx), " luminosities when trapping for NaN values."
    endif

    print, "Total Bolo Lum= ", total(Lum_Bol)


x= fload_allstars_xyz('x',center=center_bh)
y= fload_allstars_xyz('y',center=center_bh)
z= fload_allstars_xyz('z',center=center_bh)
lum= Lum_Bol(*)
ids= fload_allstars_ids(1)

rr= sqrt(x*x + y*y + z*z)
if spherical_grid eq 1 then idx= where(rr lt r_len)
if cube_grid eq 1 then begin
	idx= where((abs(x) lt 0.5*side_len) and (abs(y) lt 0.5*side_len) and (abs(z) lt 0.5*side_len))
endif
if idx(0) ne -1 then print, "Stellar Lum within Cube= ", total(m(idx)) else print, "No Mass in Cube"

N= long(n_elements(x))

DesNgb=96L
Hmax= 100.0


if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= lum(*)

        StellarLumDensity= fltarr(Ngrid)      ;  full velocity field

        ;S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/Cube/Density/density.so', $
        ;        'density', $
	;	N, $
	;	Ngrid, $
        ;        Coord, $
        ;        Masses, $
	;	DesNgb, $
	;	Hmax, $
        ;        StellarLumDensity)
	StellarLumDensity(*)= 1.0

	; Make the Density Log
	; ---------------------
	print, "Cube Stellar Lum= ", total(Grid_Volume*StellarLumDensity)
	StellarLum= alog10(Grid_Volume*StellarLumDensity)



endif else begin
        print,'No stars, no problem.'
        return
endelse

help, StellarLum



;----------------------------------------------

;----------------------------------------------
;
;   Write to Grid file
;




gridfile=extra+'griddata_'+strcompress(string(snapnum),/remove_all)+'.txt'
openw, 1, gridfile, ERROR=err
printf, 1, '# 50x50x50 ',frun,'  snapnum=',strcompress(string(snapnum),/remove_all)
printf, 1, "#                                                           Log     Log      Log      Log       Log       Log"
printf, 1, "#                                                        Tot Gas Cold Gas Cold Gas Cold Gas  Stellar  Stellar"
printf, 1, "#      x        y         z         r    theta      phi  Density  Density     Temp  Filling     Mass      Lum"
printf, 1, "#   (kpc)    (kpc)     (kpc)    (kpc)    (rad)    (rad)   (cm-3)   (cm-3)      (K)   Factor (Msolar) (Lsolar)"

for i=0L,Ngrid-1 do begin
	printf, 1, FORMAT='(3(F8.5,"  "),3(F7.3,"  "),6(F7.3,"  "))', $
		Grid[0,i],Grid[1,i],Grid[2,i], $
		Grid_rthetaphi[0,i],Grid_rthetaphi[1,i],Grid_rthetaphi[2,i], $
		Density_TotalGas[i], $
		Density_ColdGas[i], $
		Temp[i], $
		ColdFillingFac[i], $
		StellarMass[i], StellarLum[i]
endfor


close, 1

;----------------------------------------------
;
;   Write ColdGas Info to file
;


;coldgridfile='grid_coldgas.txt'
coldgridfile=extra+'griddata_'+strcompress(string(snapnum),/remove_all)+'_cold.txt'
openw, 1, coldgridfile, ERROR=err

for i=0L,Ngrid-1 do begin
        printf, 1, FORMAT='(F7.3)', Density_ColdGas[i]
endfor


close, 1

;----------------------------------------------
;
;   Write StellarLum to file
;


;slgridfile='grid_stellarlum.txt'
slgridfile=extra+'griddata_'+strcompress(string(snapnum),/remove_all)+'_slum.txt'
openw, 1, slgridfile, ERROR=err

for i=0L,Ngrid-1 do begin
        printf, 1, FORMAT='(F7.3)', StellarLum[i]
endfor


close, 1

;----------------------------------------------




end






; =============================================================================





; =============================================================================






; SHORT version of sukanya grid




;  grid provided to sukanya
; ------------------------------------
pro sukanya_grid_coldonly_original, frun, snapnum, $
		extra=extra, $
		grid_radius=grid_radius


if not keyword_set(frun) then begin
	print, "  "
	print, " sukanya_grid_coldonly, frun, snapnum"
	print, "  "
	return
endif

;frun= "/raid4/tcox/vc3vc3e"
if not keyword_set(snapnum) then snapnum= 12
;snapnum= 31
;ok=fload_snapshot_bh(frun,snapnum)

if not keyword_set(extra) then extra=''



;-------------------------------------------
;
;   Define Grid
;  --------------

;  spherical
; ------------
spherical_grid= 1
;spherical_grid= 0
if spherical_grid eq 1 then begin

	; need the following order

	; grid_*
	; ----------
	N_r= 200L     &  grid_r= fltarr(N_r)
	N_theta= 59L  &  grid_theta= fltarr(N_theta)
	N_phi= 118L   &  grid_phi= fltarr(N_phi)
	; cent_*
	; ---------
	N_r= 199L     &  grid_cen_r= fltarr(N_r)
	N_theta= 58L  &  grid_cen_theta= fltarr(N_theta)
	N_phi= 117L   &  grid_cen_phi= fltarr(N_phi)
	; test_*
	; ---------
	;N_r= 23L     &  r_data= fltarr(N_r)
	;N_theta= 9L  &  theta_data= fltarr(N_theta)
	;N_phi= 25L   &  phi_data= fltarr(N_phi)


	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_r.dat"
	readf, 1, grid_r
	close, 1
	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_r.dat"
	readf, 1, grid_cen_r
	close, 1

	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_theta.dat"
	readf, 1, grid_theta
	close, 1
	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_theta.dat"
	readf, 1, grid_cen_theta
	close, 1

	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_phi.dat"
	readf, 1, grid_phi
	close, 1
	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_phi.dat"
	readf, 1, grid_cen_phi
	close, 1

	r_len= max(grid_r)
	if keyword_set(grid_radius) then begin
		print, " ******  MANUALLY setting grid max radius to ", grid_radius, "  ******** "
		grid_r= grid_r * ( grid_radius / r_len)
		grid_cen_r= grid_cen_r * ( grid_radius / r_len)
		r_len= max(grid_r)
	endif

	print, "max(grid_r)= ", max(grid_r)
	print, "max(grid_cen_r)= ", max(grid_cen_r)
	print, "r_len= ", r_len
	total_Volume= (4./3.)*!PI*r_len*r_len*r_len
	print, "total_Volume= ", total_Volume

	Ngrid= N_r*N_theta*N_phi

	Grid= fltarr(3,Ngrid)
	Grid_rthetaphi= fltarr(3,Ngrid)
	Grid_Volume= fltarr(Ngrid)

	n_i= long(0)
	for i=0,N_r-1 do begin
	  for j=0,N_theta-1 do begin
	    for k=0,N_phi-1 do begin

		R= grid_cen_r[i]
		theta= grid_cen_theta[j]
		phi= grid_cen_phi[k]
		Grid_rthetaphi[0,n_i]= R
		Grid_rthetaphi[1,n_i]= theta
		Grid_rthetaphi[2,n_i]= phi

		x= R*sin(theta)*cos(phi)
		y= R*sin(theta)*sin(phi)
		z= R*cos(theta)
		Grid[0,n_i]= x
		Grid[1,n_i]= y
		Grid[2,n_i]= z

		;Grid_Volume[n_i]= R * 10^(d_log_r) * d_phi * d_cos_theta
		;Grid_Volume[n_i]= R * 10^(d_log_r) * d_phi * sin(theta) * d_theta
		;Grid_Volume[n_i]= R * R * d_log_r * d_phi * sin(theta) * d_theta    ;  dr= r * d_logr
		delta_r= grid_r[i+1] - grid_r[i]
		delta_theta= grid_theta[j+1] - grid_theta[j]
		delta_phi= grid_phi[k+1] - grid_phi[k]
		Grid_Volume[n_i]= R * R * delta_r * delta_phi * sin(theta) * delta_theta

		n_i= n_i + 1
	    endfor
	  endfor
	endfor
endif


print, "Grid_Volume max/min= ", max(Grid_Volume), min(Grid_Volume)
print, "Volume check= ", total(Grid_Volume)


;  Volume check
; --------------
sphere_area= 0.0
for j=0,N_theta-1 do begin
  for k=0,N_phi-1 do begin
        delta_theta= grid_theta[j+1] - grid_theta[j]
        delta_phi= grid_phi[k+1] - grid_phi[k]
        sphere_area= sphere_area + delta_phi * sin(grid_theta[j]) * delta_theta
  endfor
endfor
print, "sphere area= ", sphere_area
print, '4*PI= ', 4.0*!PI




;-------------------------------------------
;
;   Open Snapshot
;

;spawn, "/bin/ls "+frun+"/snap* | wc ",result                                    
;lastsnapshot=long(result[0])-1                                                  
;ok=fload_snapshot_bh(frun,lastsnapshot)         
ok=fload_snapshot_bh(frun,snapnum)         


;-------------------------------------------
;
;   Density - Total Gas
;

print, " ---------------------"
print, "  DETERMINE: center   "
print, " ---------------------"

; determine_center
; -----------------
center= [0,0,0]
center_bh= fload_center_alreadycomp(1)

; two black holes
if fload_npart(5) gt 1 then begin
	bhid= fload_blackhole_id(1)
	bhid1= bhid[0]
	bhid2= bhid[1]
	print, "Blackhole ID's: ", bhid1, bhid2
	center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
	center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
	center_bh= center1
	;center_bh= center2
endif

; one black hole
if fload_npart(5) eq 1 then begin
	bhid= fload_blackhole_id(1)
	center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
endif


;  grab variables
; ----------------
print, "using center= ", center_bh
x= fload_gas_xyz('x',center=center_bh)
y= fload_gas_xyz('y',center=center_bh)
z= fload_gas_xyz('z',center=center_bh)
m= fload_gas_mass(1)
ids= fload_gas_id(1)
hsml= fload_gas_hsml(1)

print, "Total Gas Mass= ", total(m)
rr= sqrt(x*x + y*y + z*z)
if spherical_grid eq 1 then idx= where(rr lt r_len)
if idx(0) ne -1 then print, "Gas mass within sphere = ", total(m(idx)) else print, "No Mass in Sphere"

; ------
rho= fload_gas_rho(1)
idx= where(rho gt 0.000854924)
print, "SF gas mass= ", total(m(idx))
idx= where((rho gt 0.000854924) and (rr lt r_len))
print, "SF gas mass within sphere= ", total(m(idx))

rho= fload_gas_rho(1)
print, "rho min/max= ", min(rho), max(rho)
rho_cold= fload_gas_rho(1,/cold)
print, "rho_cold max/min: ", max(rho_cold), min(rho_cold)
idx= where(rho_cold gt 0.0)
print, "    n (greater than zero)= ", n_elements(idx)
print, "SF (or Multiphase) Gas Mass= ", total(m(idx))

coldf= fload_gas_coldfraction(1)
print, "coldf max/min: ", max(coldf), min(coldf)

ColdGasMass= m*coldf
print, "Cold Gas Mass= ", total(ColdGasMass)
if spherical_grid eq 1 then idx= where(rr lt r_len)
if idx(0) ne -1 then print, "Cold gas mass within Sphere = ", total(ColdGasMass(idx)) else print, "No Cold Gas in Sphere"


; -----------------------------------------------
print, "================================="
print, "      (Total) Gas Density "
print, "================================="
print, "*** std density routine ***"

N= long(n_elements(x))


DesNgb=256L
;DesNgb=32L

Hmax= 0.1

if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= m(*)

        Density_TotalGas= fltarr(Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density2/density.so', $
                'density', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Density_TotalGas)


	; do some checks
	; ----------------
	print, "Density_TotalGas   max/min= ", max(Density_TotalGas), min(Density_TotalGas)
	print, "         N= ", n_elements(Density_TotalGas)
	idx=where(Density_TotalGas gt 0.0)
	print, " N_nonzero= ", n_elements(idx)
	idx=where(Density_TotalGas lt 1.0e-10)
	Density_TotalGas(idx)= 1.0d-12

	print, "Cube Gas Mass= ", total(Grid_Volume*Density_TotalGas)


	; Make the Density Log
	; ---------------------
	UnitDensity_in_cgs = 6.76991d-22
	PROTONMASS = 1.6726d-24
	factor= UnitDensity_in_cgs / PROTONMASS
	print, "Density factor= ", factor
	Density_TotalGas= alog10(factor*Density_TotalGas)     ; convert to cm-3  (is mu right?)
	;Density_TotalGas= alog10(674.59*Density_TotalGas)     ; convert to cm-3  (is mu right?)
	


endif else begin
        print,'No stars, no problem.'
        return
endelse

;help, Density_TotalGas

;----------------------------------------------
;
;   Write Total Gas Info to file
;


;coldgridfile='grid_coldgas.txt'
totalgridfile=extra+'griddata_'+strcompress(string(snapnum),/remove_all)+'_total.txt'
openw, 1, totalgridfile, ERROR=err

for i=0L,Ngrid-1 do begin
        printf, 1, FORMAT='(F9.4)', Density_TotalGas[i]
endfor


close, 1

;----------------------------------------------




print, "======================================="
print, "      (Total) Gas Density (method 2)  "
print, "======================================="
print, "      *** std density routine ***      "

N= long(n_elements(x))


DesNgb=256L
;DesNgb=32L

Hmax= 0.1

if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(8,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Coord(7,0:N-1)= hsml(*)
	Masses(0:N-1)= m(*)

        Mass_TotalGas= fltarr(Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Mass/mass.so', $
                'mass', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Mass_TotalGas)


	; do some checks
	; ----------------
	print, "Mass_TotalGas   max/min= ", max(Mass_TotalGas), min(Mass_TotalGas)
	print, "         N= ", n_elements(Mass_TotalGas)

	idx=where(Mass_TotalGas gt 0.0)
	print, " N_nonzero= ", n_elements(idx)


	idx=where(Mass_TotalGas lt 1.0e-10)
	Mass_TotalGas(idx)= 1.0d-12

	Density_TotalGas= Mass_TotalGas/Grid_Volume

	print, "Cube Gas Mass= ", total(Mass_TotalGas)


	; Make the Density Log
	; ---------------------
	UnitDensity_in_cgs = 6.76991d-22
	PROTONMASS = 1.6726d-24
	factor= UnitDensity_in_cgs / PROTONMASS
	print, "Density factor= ", factor
	Density_TotalGas= alog10(factor*Density_TotalGas)     ; convert to cm-3  (is mu right?)
	;Density_TotalGas= alog10(674.59*Density_TotalGas)     ; convert to cm-3  (is mu right?)
	


endif else begin
        print,'No stars, no problem.'
        return
endelse

;help, Density_TotalGas

;----------------------------------------------
;
;   Write Total Gas Info to file
;


;coldgridfile='grid_coldgas.txt'
totalgridfile=extra+'griddata_'+strcompress(string(snapnum),/remove_all)+'_mass.txt'
openw, 1, totalgridfile, ERROR=err

for i=0L,Ngrid-1 do begin
        printf, 1, FORMAT='(F9.4)', Density_TotalGas[i]
endfor


close, 1







;-------------------------------------------
;
;    COLD GAS
;  
; -------------------------------------------
; -----------------------------------------------
print, "================================="
print, "      Cold Gas Density "
print, "================================="
print, "*** std density routine ***"

;DesNgb=5L
;Hmax= 0.050

DesNgb=96L
;DesNgb=16L

;Hmax= 0.030
Hmax= 0.015
;Hmax= 0.0075

;DesNgb=2L
;Hmax= 0.10

if(N gt 0) then begin

	Ncoord= max([N,Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Masses(0:N-1)= ColdGasMass(*)

	print, "N_masses= ", n_elements(Masses(0:N-1))
	idx= where(ColdGasMass gt 0.0)
	print, "N_masses with non-zero mass= ", n_elements(ColdGasMass(idx))

        Density_ColdGas2= fltarr(Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density2/density.so', $
                'density', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Density_ColdGas2)


	; do some checks
	; ----------------
	print, "Density_ColdGas2   max/min= ", max(Density_ColdGas2), min(Density_ColdGas2)
	print, "         N= ", n_elements(Density_ColdGas2)
	idx=where(Density_ColdGas2 gt 0.0)
	print, " N_nonzero= ", n_elements(idx)

	print, "Cube Gas Mass (II)= ", total(Grid_Volume*Density_ColdGas2)

	; Make the Density Log
	; ---------------------
	idx=where(Density_ColdGas2 le 0.0)
	if idx(0) ne -1 then Density_ColdGas2(idx)= 1.0e-12
	Density_ColdGas2= alog10(factor*Density_ColdGas2)     ; convert to cm-3  (is mu right?)
	if idx(0) ne -1 then Density_ColdGas2(idx)= 0.0


	; use this one!
	print, " ************************** "
	print, "   HEY: using the density   "
	print, "        version not the " 
	print, "        average one    "
	print, " ************************** "
	Density_ColdGas= Density_ColdGas2
	print, "Density_ColdGas= Density_ColdGas2"
	print, " "

	;extra= 'test2_'


endif else begin
        print,'No gas, no problem.'
        return
endelse


;----------------------------------------------
;
;   Write ColdGas Info to file
;


;coldgridfile='grid_coldgas.txt'
coldgridfile=extra+'griddata_'+strcompress(string(snapnum),/remove_all)+'_cold.txt'
openw, 1, coldgridfile, ERROR=err

for i=0L,Ngrid-1 do begin
        printf, 1, FORMAT='(F9.4)', Density_ColdGas[i]
endfor


close, 1

;----------------------------------------------

end






; =============================================================================
; =============================================================================
; =============================================================================
; =============================================================================







; SHORT version of sukanya grid
;
;  this is the newest, latest, and
; greatest version
;


; ============================================
;
;  WARNING: this is NOT new, the following
;           procedure, sukanya_sph_info, is
;           more recent
;
; ============================================



;  grid provided to sukanya
; ------------------------------------
pro sukanya_grid_new, frun, snapnum, $
		extra=extra, $
		grid_radius=grid_radius


if not keyword_set(frun) then begin
	print, "  "
	print, " sukanya_grid_new, frun, snapnum"
	print, "  "
	return
endif

if not keyword_set(snapnum) then snapnum= 12

if not keyword_set(extra) then extra=''



;-------------------------------------------
;
;   Define Grid
;
;-------------------------------------------

;  spherical
; ------------
spherical_grid= 1
;spherical_grid= 0
if spherical_grid eq 1 then begin

	; need the following order

	; grid_*
	; ----------
	N_r= 200L     &  grid_r= fltarr(N_r)
	N_theta= 59L  &  grid_theta= fltarr(N_theta)
	N_phi= 118L   &  grid_phi= fltarr(N_phi)
	; cent_*
	; ---------
	N_r= 199L     &  grid_cen_r= fltarr(N_r)
	N_theta= 58L  &  grid_cen_theta= fltarr(N_theta)
	N_phi= 117L   &  grid_cen_phi= fltarr(N_phi)
	; test_*
	; ---------
	;N_r= 23L     &  r_data= fltarr(N_r)
	;N_theta= 9L  &  theta_data= fltarr(N_theta)
	;N_phi= 25L   &  phi_data= fltarr(N_phi)


	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_r.dat"
	readf, 1, grid_r
	close, 1
	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_r.dat"
	readf, 1, grid_cen_r
	close, 1

	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_theta.dat"
	readf, 1, grid_theta
	close, 1
	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_theta.dat"
	readf, 1, grid_cen_theta
	close, 1

	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_phi.dat"
	readf, 1, grid_phi
	close, 1
	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_phi.dat"
	readf, 1, grid_cen_phi
	close, 1

	r_len= max(grid_r)
	if keyword_set(grid_radius) then begin
		print, " ******  MANUALLY setting grid max radius to ", grid_radius, "  ******** "
		grid_r= grid_r * ( grid_radius / r_len)
		grid_cen_r= grid_cen_r * ( grid_radius / r_len)
		r_len= max(grid_r)
	endif

	print, "max(grid_r)= ", max(grid_r)
	print, "max(grid_cen_r)= ", max(grid_cen_r)
	print, "r_len= ", r_len
	total_Volume= (4./3.)*!PI*r_len*r_len*r_len
	print, "total_Volume= ", total_Volume

	Ngrid= N_r*N_theta*N_phi

	Grid= fltarr(3,Ngrid)
	Grid_rthetaphi= fltarr(3,Ngrid)
	Grid_Volume= fltarr(Ngrid)

	n_i= long(0)
	for i=0,N_r-1 do begin
	  for j=0,N_theta-1 do begin
	    for k=0,N_phi-1 do begin

		R= grid_cen_r[i]
		theta= grid_cen_theta[j]
		phi= grid_cen_phi[k]
		Grid_rthetaphi[0,n_i]= R
		Grid_rthetaphi[1,n_i]= theta
		Grid_rthetaphi[2,n_i]= phi

		x= R*sin(theta)*cos(phi)
		y= R*sin(theta)*sin(phi)
		z= R*cos(theta)
		Grid[0,n_i]= x
		Grid[1,n_i]= y
		Grid[2,n_i]= z

		;Grid_Volume[n_i]= R * 10^(d_log_r) * d_phi * d_cos_theta
		;Grid_Volume[n_i]= R * 10^(d_log_r) * d_phi * sin(theta) * d_theta
		;Grid_Volume[n_i]= R * R * d_log_r * d_phi * sin(theta) * d_theta    ;  dr= r * d_logr
		delta_r= grid_r[i+1] - grid_r[i]
		delta_theta= grid_theta[j+1] - grid_theta[j]
		delta_phi= grid_phi[k+1] - grid_phi[k]
		Grid_Volume[n_i]= R * R * delta_r * delta_phi * sin(theta) * delta_theta

		n_i= n_i + 1
	    endfor
	  endfor
	endfor
endif


print, "Grid_Volume max/min= ", max(Grid_Volume), min(Grid_Volume)
print, "Volume check= ", total(Grid_Volume)


;  Volume check
; --------------
sphere_area= 0.0
for j=0,N_theta-1 do begin
  for k=0,N_phi-1 do begin
        delta_theta= grid_theta[j+1] - grid_theta[j]
        delta_phi= grid_phi[k+1] - grid_phi[k]
        sphere_area= sphere_area + delta_phi * sin(grid_theta[j]) * delta_theta
  endfor
endfor
print, "sphere area= ", sphere_area
print, '4*PI= ', 4.0*!PI




;-------------------------------------------
;
;   Open Snapshot
;

;spawn, "/bin/ls "+frun+"/snap* | wc ",result                                    
;lastsnapshot=long(result[0])-1                                                  
;ok=fload_snapshot_bh(frun,lastsnapshot)         
ok=fload_snapshot_bh(frun,snapnum)         





; ----------------------------------------------------------------------
;
;   Center
;
; ----------------------------------------------------------------------

print, " ---------------------"
print, "  DETERMINE: center   "
print, " ---------------------"

; determine_center
; -----------------
center= [0,0,0]
center_bh= fload_center_alreadycomp(1)

; two black holes
if fload_npart(5) gt 1 then begin
	bhid= fload_blackhole_id(1)
	bhid1= bhid[0]
	bhid2= bhid[1]
	print, "Blackhole ID's: ", bhid1, bhid2
	center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
	center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
	center_bh= center1
	;center_bh= center2
endif

; one black hole
if fload_npart(5) eq 1 then begin
	bhid= fload_blackhole_id(1)
	center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
endif


; ----------------------------------------------------------------------


;  grab variables
; ----------------
print, "using center= ", center_bh
x= fload_gas_xyz('x',center=center_bh)
y= fload_gas_xyz('y',center=center_bh)
z= fload_gas_xyz('z',center=center_bh)
m= fload_gas_mass(1)
ids= fload_gas_id(1)
hsml= fload_gas_hsml(1)


; ----------------------------------------------------------------------


print, " ================================ "
print, " "
print, " N_Gas= ", n_elements(m)
print, " Total Gas Mass= ", total(m)
print, " "
rr= sqrt(x*x + y*y + z*z)
if spherical_grid eq 1 then idx= where(rr lt r_len)
if idx(0) ne -1 then print, " Gas mass within sphere = ", total(m(idx)) else print, "No Mass in Sphere"


; -------------------
; density information
rho= fload_gas_rho(1)
print, " density: rho min/max= ", min(rho), max(rho)


; -------------------
; -------------------

  ; define the cut-off density
  ; for multi-phase model

  rho_crit= 0.000854924     ; std GADGET value

; -------------------
; -------------------


idx= where(rho gt rho_crit)
print, "   Star-forming Gas "
print, "      N_sf= ", n_elements(idx)
print, "      Mass= ", total(m(idx))
idx= where((rho gt rho_crit) and (rr lt r_len))
print, "         (within sphere ", total(m(idx)), " )"
print, " "


; --------------------
;  This is
;
;     TurbulentFactor = 1 + (P_turb / P_thermal)
; 
;if not keyword_set(TurbulentFactor) then begin
        ;TurbulentFactor= 1.0
        ;TurbulentFactor= 6.0
        ;TurbulentFactor= 11.0
        TurbulentFactor= 101.0
;endif
; --------------------
h= 0.7
; --------------------


rho_cold= fload_gas_rho(1,/cold,TurbulentFactor=TurbulentFactor)
idx= where(rho_cold gt 0.0)
print, "   Cold Gas (should be same as Star-forming) "
print, "      N_cold (rho_cold > 0)= ", n_elements(idx)
;print, "      Mass= ", total(m(idx))
print, "      rho_cold max/min: ", max(rho_cold), min(rho_cold(idx)), min(rho_cold)

coldf= fload_gas_coldfraction(1)
print, "      cold gas fraction max/min: ", max(coldf), min(coldf)

ColdGasMass= m*coldf
print, "      Mass= ", total(ColdGasMass)

; volume_cold / volume_hot
coldgasff= fload_gas_coldgasfillingfactor(1,TurbulentFactor=TurbulentFactor)
print, "      coldgasff max/min: ", max(coldgasff), min(coldgasff)

; volume_cold / volume_total
coldgasvf= fload_gas_coldgasvolumefraction(1,TurbulentFactor=TurbulentFactor)
print, "      coldgasvf max/min: ", max(coldgasvf), min(coldgasvf)
print, " "


; BAD clumpsize calculation.  The hsml is not a good 
; indicator of the size of the SPH particle
;
;clumpsize= hsml(idx) * (coldgasvf(idx)^(1./3.))
;print, "      clumpsize max/min: ", max(clumpsize), min(clumpsize)
;print, "                   mean: ", mean(clumpsize)
;

; cold clump size
newclumpsize= ( (4.*!PI/ 3.) * rho_cold(idx) / ColdGasMass(idx))^(-1./3.)
clumpsize= newclumpsize
print, "(new) clumpsize max/min: ", max(clumpsize), min(clumpsize)
print, "                   mean: ", mean(clumpsize)

print, " "
print, " "



;if spherical_grid eq 1 then idx= where(rr lt r_len)
;if idx(0) ne -1 then print, " Cold gas mass within Sphere = ", total(ColdGasMass(idx)) else print, "No Cold Gas in Sphere"



; -----------------------------------------------
;
;    Set the variables that go into file

        N_ColdGas= long(n_elements(idx))
        ColdGas_x= x(idx) / h
        ColdGas_y= y(idx) / h
        ColdGas_z= z(idx) / h
	ColdGas_Mass= ColdGasMass(idx) / h
        ColdGas_Density= rho_cold(idx) * h * h
        ColdGas_Size= clumpsize / h
	ColdGas_Ids= ids(idx)
        SphGas_Hsml= hsml(idx) / h


; -----------------------------------------------
;
;    Add any lower density gas

add_lowdensity_gas= 1
;add_lowdensity_gas= 0

if add_lowdensity_gas eq 1 then begin

        ; low density gas
        gastemp= fload_gas_temperature(1)   ; in K
        idx= where((rho lt rho_crit) and (gastemp lt 4.0e+4))

        print, "-----------------------------"
        print, "Adding Cold Low-density gas"
        print, " "
        print, "  N= ", n_elements(idx)
        N_ColdGas= N_ColdGAs + n_elements(idx)
        r= sqrt(x(idx)*x(idx) + y(idx)*y(idx) + z(idx)*z(idx))
        print, "  radius max/min: ", max(r), min(r)
        ColdGas_x= [ColdGas_x, x(idx) / h]
        ColdGas_y= [ColdGas_y, y(idx) / h]
        ColdGas_z= [ColdGas_z, z(idx) / h] 
        ColdGas_Density= [ColdGas_Density, rho(idx) * h * h]     ; use regular rho
        print, "  max/min density: ", max(rho(idx)), min(rho(idx))
        clumpsize= ( (4.*!PI/ 3.) * rho(idx) / m(idx))^(-1./3.)
        print, "  clumpsize max/min: ", max(clumpsize), min(clumpsize)
        ColdGas_Size= [ColdGas_Size, clumpsize / h]
	ColdGas_Mass= [ColdGas_Mass, m(idx)]
	ColdGas_Ids= [ColdGas_Ids, ids(idx)]
        SphGas_Hsml= [SphGas_Hsml, hsml(idx) / h]
        print, "-----------------------------"
endif


rr= sqrt(ColdGas_x*ColdGas_x + ColdGas_y*ColdGas_y + ColdGas_z*ColdGas_z)
if spherical_grid eq 1 then idx= where(rr lt r_len)
if idx(0) ne -1 then print, " Cold gas mass within sphere = ", total(ColdGas_Mass(idx)) else print, "No Mass in Sphere"
print, "-----------------------------"




; -----------------------------------------------
;
;  Determine Total Gas Density on the Grid
;
;      version #1
;
; -----------------------------------------------

do_std1= 1
;do_std1= 0

if do_std1 eq 1 then begin
	print, "==================================================================="
	print, "                Total Gas Density   (version #1)                   "
	print, "==================================================================="
	print, "using: /home2/tcox/Tools/C-Routines_for_IDL/Cube/Density2/density.so"
	print, " "

	N= long(n_elements(x))


	DesNgb=256L
	;DesNgb=32L

	Hmax= 0.1

	if(N gt 0) then begin

		Ncoord= max([N,Ngrid])

		Coord=fltarr(7,Ncoord)
		Masses= fltarr(Ncoord)

		Coord(0,0:N-1)= x(*)
		Coord(1,0:N-1)= y(*)
		Coord(2,0:N-1)= z(*)
		Coord(3,0:N-1)= ids(*)
		Coord(4,0:Ngrid-1)= Grid(0,*)
		Coord(5,0:Ngrid-1)= Grid(1,*)
		Coord(6,0:Ngrid-1)= Grid(2,*)
		Masses(0:N-1)= m(*)

	        Density_TotalGas= fltarr(Ngrid)      ;  full velocity field

	        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density2/density.so', $
	                'density', $
			N, $
			Ngrid, $
	                Coord, $
	                Masses, $
			DesNgb, $
			Hmax, $
	                Density_TotalGas)


		; do some checks
		; ----------------
		print, "Density_TotalGas   max/min= ", max(Density_TotalGas), min(Density_TotalGas)
		print, "         N= ", n_elements(Density_TotalGas)
		idx=where(Density_TotalGas gt 0.0)
		print, " N_nonzero= ", n_elements(idx)
		idx=where(Density_TotalGas lt 1.0e-10)
		Density_TotalGas(idx)= 1.0d-12

		print, "Cube Gas Mass= ", total(Grid_Volume*Density_TotalGas)


		; Make the Density Log
		; ---------------------
		UnitDensity_in_cgs = 6.76991d-22
		PROTONMASS = 1.6726d-24
		factor= UnitDensity_in_cgs / PROTONMASS
		print, "Density factor= ", factor
		Density_TotalGas= alog10(factor*Density_TotalGas)     ; convert to cm-3  (is mu right?)
		;Density_TotalGas= alog10(674.59*Density_TotalGas)     ; convert to cm-3  (is mu right?)
	


	endif else begin
	        print,'No stars, no problem.'
	        return
	endelse

	;help, Density_TotalGas

	;----------------------------------------------
	;
	;   Write Total Gas Info to file
	;


	;coldgridfile='grid_coldgas.txt'
	totalgridfile=extra+'totalgrid1_'+strcompress(string(snapnum),/remove_all)+'.txt'
	openw, 1, totalgridfile, ERROR=err

	for i=0L,Ngrid-1 do begin
	        printf, 1, FORMAT='(F9.4)', Density_TotalGas[i]
	endfor


	close, 1
endif



; -----------------------------------------------
;
;  Determine Total Gas Density on the Grid
;
;      version #2
;
; -----------------------------------------------

do_std2= 1
;do_std2= 0

if do_std2 eq 1 then begin
	print, "==================================================================="
	print, "                Total Gas Density   (version #2)                   "
	print, "==================================================================="
	print, "using: /home2/tcox/Tools/C-Routines_for_IDL/Cube/Mass/mass.so"
	print, " "

	N= long(n_elements(x))

        sphsize= ( (4.*!PI/ 3.) * rho / m)^(-1./3.)

	; not currently being used by the Mass routine
	DesNgb=256L
	;DesNgb=32L
	Hmax= 1.0

	if(N gt 0) then begin

		Ncoord= max([N,Ngrid])

		Coord=fltarr(8,Ncoord)
		Masses= fltarr(Ncoord)

		Coord(0,0:N-1)= x(*)
		Coord(1,0:N-1)= y(*)
		Coord(2,0:N-1)= z(*)
		Coord(3,0:N-1)= ids(*)
		Coord(4,0:Ngrid-1)= Grid(0,*)
		Coord(5,0:Ngrid-1)= Grid(1,*)
		Coord(6,0:Ngrid-1)= Grid(2,*)
		Coord(7,0:N-1)= sphsize(*)
		Masses(0:N-1)= m(*)

	        Mass_TotalGas= fltarr(Ngrid)      ;  full velocity field

        	S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Mass/mass.so', $
                	'mass', $
			N, $
			Ngrid, $
                	Coord, $
                	Masses, $
			DesNgb, $
			Hmax, $
                	Mass_TotalGas)



	        ; do some checks
	        ; ----------------
	        print, "Mass_TotalGas   max/min= ", max(Mass_TotalGas), min(Mass_TotalGas)
	        print, "         N= ", n_elements(Mass_TotalGas)

	        idx=where(Mass_TotalGas gt 0.0)
	        print, " N_nonzero= ", n_elements(idx)


        	idx=where(Mass_TotalGas lt 1.0e-10)
        	Mass_TotalGas(idx)= 1.0d-12

        	Density_TotalGas= Mass_TotalGas/Grid_Volume
		print, "Density_TotalGas   max/min= ", max(Density_TotalGas), min(Density_TotalGas)

        	print, " Total gas mass in grid= ", total(Mass_TotalGas)


		; how many non-zero
		idx=where(Density_TotalGas gt 0.0)
		print, " N_nonzero= ", n_elements(idx)
		idx=where(Density_TotalGas lt 1.0e-10)
		if idx(0) ne -1 then Density_TotalGas(idx)= 1.0d-12

		print, " Total gas mass in grid #2= ", total(Grid_Volume*Density_TotalGas)


		; Make the Density Log
		; ---------------------
		UnitDensity_in_cgs = 6.76991d-22
		PROTONMASS = 1.6726d-24
		factor= UnitDensity_in_cgs / PROTONMASS
		print, "Density factor= ", factor
		Density_TotalGas= alog10(factor*Density_TotalGas)     ; convert to cm-3  (is mu right?)
		;Density_TotalGas= alog10(674.59*Density_TotalGas)     ; convert to cm-3  (is mu right?)
	


	endif else begin
	        print,'No stars, no problem.'
	        return
	endelse

	;help, Density_TotalGas

	;----------------------------------------------
	;
	;   Write Total Gas Info to file
	;


	;coldgridfile='grid_coldgas.txt'
	totalgridfile=extra+'totalgrid2_'+strcompress(string(snapnum),/remove_all)+'.txt'
	openw, 1, totalgridfile, ERROR=err

	for i=0L,Ngrid-1 do begin
	        printf, 1, FORMAT='(F9.4)', Density_TotalGas[i]
	endfor


	close, 1
endif


;----------------------------------------------



print, "==================================================================="
print, "                        Cold Gas Density                           "
print, "==================================================================="
print, "using: /home2/tcox/Tools/C-Routines_for_IDL/Cube/Mass/mass.so"
print, " "



; neither of these are really used in the new routine
DesNgb=256L
;DesNgb=32L
Hmax= 1.0

if(N_ColdGas gt 0) then begin

	N= long(N_ColdGas)
	Ncoord= max([N,Ngrid])

	Coord=fltarr(8,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= ColdGas_x(*)
	Coord(1,0:N-1)= ColdGas_y(*)
	Coord(2,0:N-1)= ColdGas_z(*)
	Coord(3,0:N-1)= ColdGas_Ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Coord(7,0:N-1)= ColdGas_Size(*)
	Masses(0:N-1)= ColdGas_Mass(*)
	print, "      Total(Masses)= ", total(Masses)

        Mass_ColdGas= fltarr(Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Mass/mass.so', $
                'mass', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Mass_ColdGas)


	; do some checks
	; ----------------
	print, "Mass_ColdGas   max/min= ", max(Mass_ColdGas), min(Mass_ColdGas)
	print, "      Total(Mass)= ", total(Mass_ColdGas)
	print, "         N= ", n_elements(Mass_ColdGas)

	idx=where(Mass_ColdGas gt 0.0)
	print, " N_nonzero= ", n_elements(idx)


	idx=where(Mass_ColdGas lt 1.0e-10)
	Mass_ColdGas(idx)= 1.0d-12

	Density_ColdGas= Mass_ColdGas/Grid_Volume

	print, " Total cold gas mass in grid= ", total(Mass_ColdGas)


	; Make the Density Log
	; ---------------------
	UnitDensity_in_cgs = 6.76991d-22
	PROTONMASS = 1.6726d-24
	factor= UnitDensity_in_cgs / PROTONMASS
	print, "Density factor= ", factor
	Density_ColdGas= alog10(factor*Density_ColdGas)     ; convert to cm-3  (is mu right?)
	;Density_ColdGas= alog10(674.59*Density_ColdGas)     ; convert to cm-3  (is mu right?)
	



	; -------------------------------------------
	;
	;   Write Total Gas Info to file
	;


	gridfile=extra+'coldgrid_'+strcompress(string(snapnum),/remove_all)+'.txt'
	openw, 1, gridfile, ERROR=err

	for i=0L,Ngrid-1 do begin
	        printf, 1, FORMAT='(F9.4)', Density_ColdGas[i]
	endfor

	close, 1


endif else begin
        print,'No stars, no problem.'
        return
endelse




; -------------------------------------------

end






; =============================================================================
; =============================================================================
; =============================================================================
; =============================================================================









; ============================================
;
;   Version as of 2/16/06
;
;  has the new finer fixed grid, and it
; the newest, most up-to-date procedure.
;
; some of the grid stuff needs some help.  
;
; ============================================

pro sukanya_grid_xxx, frun, snapnum, $
		extra=extra, $
		grid_radius=grid_radius


if not keyword_set(frun) then begin
	print, "  "
	print, " sukanya_grid_xxx, frun, snapnum"
	print, "  "
	return
endif

if not keyword_set(snapnum) then snapnum= 12

if not keyword_set(extra) then extra=''



;-------------------------------------------
;
;   Define Grid
;
;-------------------------------------------


get_grid_info, Ngrid, Grid, Grid_Volume, Grid_MaxRadius




;-------------------------------------------
;
;   Open Snapshot
;

;spawn, "/bin/ls "+frun+"/snap* | wc ",result                                    
;lastsnapshot=long(result[0])-1                                                  
;ok=fload_snapshot_bh(frun,lastsnapshot)         
ok=fload_snapshot_bh(frun,snapnum)         





; ----------------------------------------------------------------------
;
;   Center
;
; ----------------------------------------------------------------------

print, " ---------------------"
print, "  DETERMINE: center   "
print, " ---------------------"

; determine_center
; -----------------
center= [0,0,0]
center_bh= fload_center_alreadycomp(1)

; two black holes
if fload_npart(5) gt 1 then begin
	bhid= fload_blackhole_id(1)
	bhid1= bhid[0]
	bhid2= bhid[1]
	print, "Blackhole ID's: ", bhid1, bhid2
	center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
	center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
	center_bh= center1
	;center_bh= center2
endif

; one black hole
if fload_npart(5) eq 1 then begin
	bhid= fload_blackhole_id(1)
	center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
endif


; ----------------------------------------------------------------------


;  grab variables
; ----------------
print, "using center= ", center_bh
x= fload_gas_xyz('x',center=center_bh)
y= fload_gas_xyz('y',center=center_bh)
z= fload_gas_xyz('z',center=center_bh)
m= fload_gas_mass(1)
ids= fload_gas_id(1)
hsml= fload_gas_hsml(1)


; ----------------------------------------------------------------------


print, " ================================ "
print, " "
print, " N_Gas= ", n_elements(m)
print, " Total Gas Mass= ", total(m)
print, " "
rr= sqrt(x*x + y*y + z*z)
idx= where(rr lt Grid_MaxRadius)
if idx(0) ne -1 then print, " Gas mass within sphere = ", total(m(idx)) else print, "No Mass in Sphere"


; -------------------
; density information
rho= fload_gas_rho(1)
print, " density: rho min/max= ", min(rho), max(rho)


; -------------------
; -------------------

  ; define the cut-off density
  ; for multi-phase model

  rho_crit= 0.000854924     ; std GADGET value

; -------------------
; -------------------


idx= where(rho gt rho_crit)
if idx(0) eq -1 then begin
	print, " "
	print, " "
	print, " "
	print, " HEY: there is no cold gas in this snapshot!!"
	print, " "
	print, " "
	print, " "
	return
endif
print, "   Star-forming Gas "
print, "      N_sf= ", n_elements(idx)
print, "      Mass= ", total(m(idx))
idx= where((rho gt rho_crit) and (rr lt Grid_MaxRadius))
print, "         (within sphere ", total(m(idx)), " )"
print, " "


; --------------------
;  This is
;
;     TurbulentFactor = 1 + (P_turb / P_thermal)
; 
;if not keyword_set(TurbulentFactor) then begin
        ;TurbulentFactor= 1.0
        ;TurbulentFactor= 6.0
        ;TurbulentFactor= 11.0
        TurbulentFactor= 101.0
;endif
; --------------------
h= 0.7
; --------------------


rho_cold= fload_gas_rho(1,/cold,TurbulentFactor=TurbulentFactor)
idx= where(rho_cold gt 0.0)
print, "   Cold Gas (should be same as Star-forming) "
print, "      N_cold (rho_cold > 0)= ", n_elements(idx)
;print, "      Mass= ", total(m(idx))
print, "      rho_cold max/min: ", max(rho_cold), min(rho_cold(idx)), min(rho_cold)

coldf= fload_gas_coldfraction(1)
print, "      cold gas fraction max/min: ", max(coldf), min(coldf)

ColdGasMass= m*coldf
print, "      Mass= ", total(ColdGasMass)

; volume_cold / volume_hot
coldgasff= fload_gas_coldgasfillingfactor(1,TurbulentFactor=TurbulentFactor)
print, "      coldgasff max/min: ", max(coldgasff), min(coldgasff)

; volume_cold / volume_total
coldgasvf= fload_gas_coldgasvolumefraction(1,TurbulentFactor=TurbulentFactor)
print, "      coldgasvf max/min: ", max(coldgasvf), min(coldgasvf)
print, " "


; BAD clumpsize calculation.  The hsml is not a good 
; indicator of the size of the SPH particle
;
;clumpsize= hsml(idx) * (coldgasvf(idx)^(1./3.))
;print, "      clumpsize max/min: ", max(clumpsize), min(clumpsize)
;print, "                   mean: ", mean(clumpsize)
;

; cold clump size
newclumpsize= ( (4.*!PI/ 3.) * rho_cold(idx) / ColdGasMass(idx))^(-1./3.)
clumpsize= newclumpsize
print, "(new) clumpsize max/min: ", max(clumpsize), min(clumpsize)
print, "                   mean: ", mean(clumpsize)

print, " "
print, " "



;if spherical_grid eq 1 then idx= where(rr lt r_len)
;if idx(0) ne -1 then print, " Cold gas mass within Sphere = ", total(ColdGasMass(idx)) else print, "No Cold Gas in Sphere"



; -----------------------------------------------
;
;    Set the variables that go into file

        N_ColdGas= long(n_elements(idx))
        ColdGas_x= x(idx) / h
        ColdGas_y= y(idx) / h
        ColdGas_z= z(idx) / h
	ColdGas_Mass= ColdGasMass(idx) / h
        ColdGas_Density= rho_cold(idx) * h * h
        ColdGas_Size= clumpsize / h
	ColdGas_Ids= ids(idx)
        SphGas_Hsml= hsml(idx) / h


; -----------------------------------------------
;
;    Add any lower density gas

;add_lowdensity_gas= 1
add_lowdensity_gas= 0

if add_lowdensity_gas eq 1 then begin

        ; low density gas
        gastemp= fload_gas_temperature(1)   ; in K
        idx= where((rho lt rho_crit) and (gastemp lt 4.0e+4))

        print, "-----------------------------"
        print, "Adding Cold Low-density gas"
        print, " "
        print, "  N= ", n_elements(idx)
        N_ColdGas= N_ColdGAs + n_elements(idx)
        r= sqrt(x(idx)*x(idx) + y(idx)*y(idx) + z(idx)*z(idx))
        print, "  radius max/min: ", max(r), min(r)
        ColdGas_x= [ColdGas_x, x(idx) / h]
        ColdGas_y= [ColdGas_y, y(idx) / h]
        ColdGas_z= [ColdGas_z, z(idx) / h] 
        ColdGas_Density= [ColdGas_Density, rho(idx) * h * h]     ; use regular rho
        print, "  max/min density: ", max(rho(idx)), min(rho(idx))
        clumpsize= ( (4.*!PI/ 3.) * rho(idx) / m(idx))^(-1./3.)
        print, "  clumpsize max/min: ", max(clumpsize), min(clumpsize)
        ColdGas_Size= [ColdGas_Size, clumpsize / h]
	ColdGas_Mass= [ColdGas_Mass, m(idx)]
	ColdGas_Ids= [ColdGas_Ids, ids(idx)]
        SphGas_Hsml= [SphGas_Hsml, hsml(idx) / h]
        print, "-----------------------------"
endif


rr= sqrt(ColdGas_x*ColdGas_x + ColdGas_y*ColdGas_y + ColdGas_z*ColdGas_z)
idx= where(rr lt Grid_MaxRadius)
if idx(0) ne -1 then print, " Cold gas mass within sphere = ", total(ColdGas_Mass(idx)) else print, "No Mass in Sphere"
print, "-----------------------------"




; -----------------------------------------------
;
;  Determine Total Gas Density on the Grid
;
;      version #1
;
; -----------------------------------------------

;do_std1= 1
do_std1= 0

if do_std1 eq 1 then begin
	print, "==================================================================="
	print, "                Total Gas Density   (version #1)                   "
	print, "==================================================================="
	print, "using: /home2/tcox/Tools/C-Routines_for_IDL/Cube/Density2/density.so"
	print, " "

	N= long(n_elements(x))


	DesNgb=256L
	;DesNgb=32L

	Hmax= 0.1

	if(N gt 0) then begin

		Ncoord= max([N,Ngrid])

		Coord=fltarr(7,Ncoord)
		Masses= fltarr(Ncoord)

		Coord(0,0:N-1)= x(*)
		Coord(1,0:N-1)= y(*)
		Coord(2,0:N-1)= z(*)
		Coord(3,0:N-1)= ids(*)
		Coord(4,0:Ngrid-1)= Grid(0,*)
		Coord(5,0:Ngrid-1)= Grid(1,*)
		Coord(6,0:Ngrid-1)= Grid(2,*)
		Masses(0:N-1)= m(*)

	        Density_TotalGas= fltarr(Ngrid)      ;  full velocity field

	        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Density2/density.so', $
	                'density', $
			N, $
			Ngrid, $
	                Coord, $
	                Masses, $
			DesNgb, $
			Hmax, $
	                Density_TotalGas)


		; do some checks
		; ----------------
		print, "Density_TotalGas   max/min= ", max(Density_TotalGas), min(Density_TotalGas)
		print, "         N= ", n_elements(Density_TotalGas)
		idx=where(Density_TotalGas gt 0.0)
		print, " N_nonzero= ", n_elements(idx)
		idx=where(Density_TotalGas lt 1.0e-10)
		Density_TotalGas(idx)= 1.0d-12

		print, "Cube Gas Mass= ", total(Grid_Volume*Density_TotalGas)


		; Make the Density Log
		; ---------------------
		UnitDensity_in_cgs = 6.76991d-22
		PROTONMASS = 1.6726d-24
		factor= UnitDensity_in_cgs / PROTONMASS
		print, "Density factor= ", factor
		Density_TotalGas= alog10(factor*Density_TotalGas)     ; convert to cm-3  (is mu right?)
		;Density_TotalGas= alog10(674.59*Density_TotalGas)     ; convert to cm-3  (is mu right?)
	


	endif else begin
	        print,'No stars, no problem.'
	        return
	endelse

	;help, Density_TotalGas

	;----------------------------------------------
	;
	;   Write Total Gas Info to file
	;


	;coldgridfile='grid_coldgas.txt'
	totalgridfile=extra+'totalgrid1_'+strcompress(string(snapnum),/remove_all)+'.txt'
	openw, 1, totalgridfile, ERROR=err

	for i=0L,Ngrid-1 do begin
	        printf, 1, FORMAT='(F9.4)', Density_TotalGas[i]
	endfor


	close, 1
endif



; -----------------------------------------------
;
;  Determine Total Gas Density on the Grid
;
;      version #2
;
; -----------------------------------------------

;do_std2= 1
do_std2= 0

if do_std2 eq 1 then begin
	print, "==================================================================="
	print, "                Total Gas Density   (version #2)                   "
	print, "==================================================================="
	print, "using: /home2/tcox/Tools/C-Routines_for_IDL/Cube/Mass/mass.so"
	print, " "

	N= long(n_elements(x))

        sphsize= ( (4.*!PI/ 3.) * rho / m)^(-1./3.)

	; not currently being used by the Mass routine
	DesNgb=256L
	;DesNgb=32L
	Hmax= 1.0

	if(N gt 0) then begin

		Ncoord= max([N,Ngrid])

		Coord=fltarr(8,Ncoord)
		Masses= fltarr(Ncoord)

		Coord(0,0:N-1)= x(*)
		Coord(1,0:N-1)= y(*)
		Coord(2,0:N-1)= z(*)
		Coord(3,0:N-1)= ids(*)
		Coord(4,0:Ngrid-1)= Grid(0,*)
		Coord(5,0:Ngrid-1)= Grid(1,*)
		Coord(6,0:Ngrid-1)= Grid(2,*)
		Coord(7,0:N-1)= sphsize(*)
		Masses(0:N-1)= m(*)

	        Mass_TotalGas= fltarr(Ngrid)      ;  full velocity field

        	S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Mass/mass.so', $
                	'mass', $
			N, $
			Ngrid, $
                	Coord, $
                	Masses, $
			DesNgb, $
			Hmax, $
                	Mass_TotalGas)



	        ; do some checks
	        ; ----------------
	        print, "Mass_TotalGas   max/min= ", max(Mass_TotalGas), min(Mass_TotalGas)
	        print, "         N= ", n_elements(Mass_TotalGas)

	        idx=where(Mass_TotalGas gt 0.0)
	        print, " N_nonzero= ", n_elements(idx)


        	idx=where(Mass_TotalGas lt 1.0e-10)
        	Mass_TotalGas(idx)= 1.0d-12

        	Density_TotalGas= Mass_TotalGas/Grid_Volume
		print, "Density_TotalGas   max/min= ", max(Density_TotalGas), min(Density_TotalGas)

        	print, " Total gas mass in grid= ", total(Mass_TotalGas)


		; how many non-zero
		idx=where(Density_TotalGas gt 0.0)
		print, " N_nonzero= ", n_elements(idx)
		idx=where(Density_TotalGas lt 1.0e-10)
		if idx(0) ne -1 then Density_TotalGas(idx)= 1.0d-12

		print, " Total gas mass in grid #2= ", total(Grid_Volume*Density_TotalGas)


		; Make the Density Log
		; ---------------------
		UnitDensity_in_cgs = 6.76991d-22
		PROTONMASS = 1.6726d-24
		factor= UnitDensity_in_cgs / PROTONMASS
		print, "Density factor= ", factor
		Density_TotalGas= alog10(factor*Density_TotalGas)     ; convert to cm-3  (is mu right?)
		;Density_TotalGas= alog10(674.59*Density_TotalGas)     ; convert to cm-3  (is mu right?)
	


	endif else begin
	        print,'No stars, no problem.'
	        return
	endelse

	;help, Density_TotalGas

	;----------------------------------------------
	;
	;   Write Total Gas Info to file
	;


	;coldgridfile='grid_coldgas.txt'
	totalgridfile=extra+'totalgrid2_'+strcompress(string(snapnum),/remove_all)+'.txt'
	openw, 1, totalgridfile, ERROR=err

	for i=0L,Ngrid-1 do begin
	        printf, 1, FORMAT='(F9.4)', Density_TotalGas[i]
	endfor


	close, 1
endif


;----------------------------------------------



; neither of these are really used in the new routine
DesNgb=256L
;DesNgb=32L
Hmax= 1.0

skip_this= 1
if (skip_this ne 1) then begin
;if(N_ColdGas gt 0) then begin

	print, "==================================================================="
	print, "                        Cold Gas Density                           "
	print, "==================================================================="
	print, "using: /home2/tcox/Tools/C-Routines_for_IDL/Cube/Mass/mass.so"
	print, " "


	N= long(N_ColdGas)
	Ncoord= max([N,Ngrid])

	Coord=fltarr(8,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= ColdGas_x(*)
	Coord(1,0:N-1)= ColdGas_y(*)
	Coord(2,0:N-1)= ColdGas_z(*)
	Coord(3,0:N-1)= ColdGas_Ids(*)
	Coord(4,0:Ngrid-1)= Grid(0,*)
	Coord(5,0:Ngrid-1)= Grid(1,*)
	Coord(6,0:Ngrid-1)= Grid(2,*)
	Coord(7,0:N-1)= ColdGas_Size(*)
	Masses(0:N-1)= ColdGas_Mass(*)
	print, "      Total(Masses)= ", total(Masses)

        Mass_ColdGas= fltarr(Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Mass/mass.so', $
                'mass', $
		N, $
		Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Mass_ColdGas)


	; do some checks
	; ----------------
	print, "Mass_ColdGas   max/min= ", max(Mass_ColdGas), min(Mass_ColdGas)
	print, "      Total(Mass)= ", total(Mass_ColdGas)
	print, "         N= ", n_elements(Mass_ColdGas)

	idx=where(Mass_ColdGas gt 0.0)
	print, " N_nonzero= ", n_elements(idx)


	idx=where(Mass_ColdGas lt 1.0e-10)
	Mass_ColdGas(idx)= 1.0d-12

	Density_ColdGas= Mass_ColdGas/Grid_Volume

	print, " Total cold gas mass in grid= ", total(Mass_ColdGas)


	; Make the Density Log
	; ---------------------
	UnitDensity_in_cgs = 6.76991d-22
	PROTONMASS = 1.6726d-24
	factor= UnitDensity_in_cgs / PROTONMASS
	print, "Density factor= ", factor
	Density_ColdGas= alog10(factor*Density_ColdGas)     ; convert to cm-3  (is mu right?)
	;Density_ColdGas= alog10(674.59*Density_ColdGas)     ; convert to cm-3  (is mu right?)
	



	; -------------------------------------------
	;
	;   Write Total Gas Info to file
	;


	gridfile=extra+'coldgrid_'+strcompress(string(snapnum),/remove_all)+'.txt'
	openw, 1, gridfile, ERROR=err

	for i=0L,Ngrid-1 do begin
	        printf, 1, FORMAT='(F9.4)', Density_ColdGas[i]
	endfor

	close, 1


endif


;----------------------------------------------



if(N_ColdGas gt 0) then begin

	print, "==================================================================="
	print, "                        Cold Gas Density                           "
	print, "==================================================================="
	print, " using homemade IDL method"
	print, " "

	N= long(N_ColdGas)
	;Ncoord= max([N,Ngrid])

	;Coord=fltarr(8,Ncoord)
	;Masses= fltarr(Ncoord)

	;Coord(0,0:N-1)= ColdGas_x(*)
	;Coord(1,0:N-1)= ColdGas_y(*)
	;Coord(2,0:N-1)= ColdGas_z(*)
	;Coord(3,0:N-1)= ColdGas_Ids(*)
	;Coord(4,0:Ngrid-1)= Grid(0,*)
	;Coord(5,0:Ngrid-1)= Grid(1,*)
	;Coord(6,0:Ngrid-1)= Grid(2,*)
	;Coord(7,0:N-1)= ColdGas_Size(*)
	;Masses(0:N-1)= ColdGas_Mass(*)
	;print, "      Total(Masses)= ", total(Masses)
	print, "      Total(ColdGas_Mass)= ", total(ColdGas_Mass)

        Mass_ColdGas= fltarr(Ngrid)

	; go through particles and distribute cold gas
	; ---------------------------------------------
	print, "N_ColdGas= ", N_ColdGas
	for i= 0L, N_ColdGas-1 do begin

		if (i mod 100) eq 0 then print, "i= ", i

		xi= ColdGas_x(i)
		yi= ColdGas_y(i)
		zi= ColdGas_z(i)
		grid_rel_point_x= Grid(0,*) - xi
		grid_rel_point_y= Grid(1,*) - yi
		grid_rel_point_z= Grid(2,*) - zi

		grid_rel_r= sqrt(grid_rel_point_x*grid_rel_point_x + $
					grid_rel_point_y*grid_rel_point_y + $
					grid_rel_point_z*grid_rel_point_z)

		idx= where(grid_rel_r le ColdGas_Size(i))

		if idx(0) eq -1 then begin
			;print, " "
			;print, " no close grid center"
			;print, " i=", i
			;print, " x|y|z=", xi, yi, zi
			;print, " size=", ColdGas_Size(i)
			;print, " min dist to grid center=", min(grid_rel_r)
			;print, " "

			idx= where(grid_rel_r eq min(grid_rel_r))
			Mass_ColdGas(idx)= Mass_ColdGas(idx) + ColdGas_Mass(i)
		endif else begin
			mi= ColdGas_Mass(i) / n_elements(idx)
			Mass_ColdGas(idx)= Mass_ColdGas(idx) + mi
		endelse

	endfor


	; do some checks
	; ----------------
	print, "Mass_ColdGas   max/min= ", max(Mass_ColdGas), min(Mass_ColdGas)
	print, "      Total(Mass)= ", total(Mass_ColdGas)
	print, "         N= ", n_elements(Mass_ColdGas)

	idx=where(Mass_ColdGas gt 0.0)
	print, " N_nonzero= ", n_elements(idx)


	idx=where(Mass_ColdGas lt 1.0e-10)
	Mass_ColdGas(idx)= 1.0d-12

	Density_ColdGas= Mass_ColdGas/Grid_Volume

	print, " Total cold gas mass in grid= ", total(Mass_ColdGas)


	; Make the Density Log
	; ---------------------
	UnitDensity_in_cgs = 6.76991d-22
	PROTONMASS = 1.6726d-24
	factor= UnitDensity_in_cgs / PROTONMASS
	print, "Density factor= ", factor
	Density_ColdGas= alog10(factor*Density_ColdGas)     ; convert to cm-3  (is mu right?)
	;Density_ColdGas= alog10(674.59*Density_ColdGas)     ; convert to cm-3  (is mu right?)
	



	; -------------------------------------------
	;
	;   Write Total Gas Info to file
	;


	gridfile=extra+'coldgrid_'+strcompress(string(snapnum),/remove_all)+'.txt'
	openw, 1, gridfile, ERROR=err

	for i=0L,Ngrid-1 do begin
	        printf, 1, FORMAT='(F9.4)', Density_ColdGas[i]
	endfor

	close, 1


endif






; -------------------------------------------

end








; =============================================================================
; =============================================================================
; =============================================================================
; =============================================================================




pro loop_over_dir, junk


frun= "/raid4/tcox/vc3vc3e_2"
;frun= "/raid4/tcox/vc3vc3h_2"
;frun= "/raid4/tcox/sbw/sb1"

for i=44, 57 do begin
;for i=1, 107 do begin
;for i=94, 107 do begin
	sukanya_sph_info, frun, i, TurbulentFactor= 101.0
endfor


end





;   SPH Particle Version



; sukanya grid
;
;  this (no, seriously, this is REALLY it!) 
; is the newest, latest, and
; greatest version
;




;  grid provided to sukanya
; ------------------------------------
pro sukanya_sph_info, frun, snapnum, $
		extra=extra, $
		TurbulentFactor=TurbulentFactor, $
		grid_radius=grid_radius


if not keyword_set(frun) then begin
	print, "  "
	print, " sukanya_sph_info, frun, snapnum"
	print, "  "
	return
endif

if not keyword_set(snapnum) then snapnum= 12

if not keyword_set(extra) then extra=''



;-------------------------------------------
;
;   Open Snapshot
;

;spawn, "/bin/ls "+frun+"/snap* | wc ",result                                    
;lastsnapshot=long(result[0])-1                                                  
;ok=fload_snapshot_bh(frun,lastsnapshot)         
ok=fload_snapshot_bh(frun,snapnum)
;ok=fload_snapshot_bh(frun,snapnum,/nopot_in_snap)




; ----------------------------------------------------------------------
;
;   Center
;
; ----------------------------------------------------------------------

print, " ---------------------"
print, "  DETERMINE: center   "
print, " ---------------------"

; determine_center
; -----------------
center= [0,0,0]
center_bh= fload_center_alreadycomp(1)

; two black holes
if fload_npart(5) gt 1 then begin
	bhid= fload_blackhole_id(1)
	bhid1= bhid[0]
	bhid2= bhid[1]
	print, "Blackhole ID's: ", bhid1, bhid2
	center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
	center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
	center_bh= center1
	;center_bh= center2
endif

; one black hole
if fload_npart(5) eq 1 then begin
	bhid= fload_blackhole_id(1)
	center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
endif




; ----------------------------------------------------------------------
;
;  grab variables
;
; ----------------------------------------------------------------------

print, "using center= ", center_bh
x= fload_gas_xyz('x',center=center_bh)
y= fload_gas_xyz('y',center=center_bh)
z= fload_gas_xyz('z',center=center_bh)
m= fload_gas_mass(1)
ids= fload_gas_id(1)
hsml= fload_gas_hsml(1)
zmets= fload_gas_metallicity(1)




; ----------------------------------------------------------------------
;
;  determine cold gas
;
; ----------------------------------------------------------------------


print, " ================================ "
print, " "
print, " N_Gas= ", n_elements(m)
print, " Total Gas Mass= ", total(m)
print, " "

; -------------------
; density information
rho= fload_gas_rho(1)
print, " density: rho min/max= ", min(rho), max(rho)


; -------------------
; -------------------

  ; define the cut-off density
  ; for multi-phase model

  rho_crit= 0.000854924     ; std GADGET value

; -------------------
; -------------------


; --------------------
;  This is
;
;     TurbulentFactor = 1 + (P_turb / P_thermal)
; 
if not keyword_set(TurbulentFactor) then begin
        TurbulentFactor= 1.0
        ;TurbulentFactor= 6.0
        ;TurbulentFactor= 11.0
        ;TurbulentFactor= 101.0
endif
print, "xxxxxxxxxxxxxxxxxxxxxx"
print, " TurbulentFactor= ", TurbulentFactor
print, "xxxxxxxxxxxxxxxxxxxxxx"
; --------------------
h= 0.7
; --------------------




idx= where(rho gt rho_crit)
if idx(0) eq -1 then begin
	N_ColdGas= 0
	snaplbl= strcompress(string(snapnum),/remove_all)
	gridfile=extra+'griddata_'+snaplbl+'_sph.txt'
	openw, 1, gridfile, ERROR=err
	printf, 1, '#  SPH cold clumps for: '+frun+', snap= '+snaplbl+' (T= '+fload_timelbl(h,3,/noteq)+' Gyr)'
	printf, 1, "#             N_clumps: "+string(N_ColdGas)+"    TurbFactor: "+string(TurbulentFactor)
	printf, 1, "#               center: ", center_bh
	printf, 1, "#     x          y          z   log(density) log(size) "
	printf, 1, "#   (kpc)      (kpc)      (kpc)     (cm-3)      (kpc)  "
	printf, 1, " "
	printf, 1, " NO COLD GAS"
	printf, 1, " "
	close, 1
	return
endif
print, "   Star-forming Gas "
print, "      N_sf= ", n_elements(idx)
print, "      Mass= ", total(m(idx))
print, " "




rho_cold= fload_gas_rho(1,/cold,TurbulentFactor=TurbulentFactor)
idx= where(rho_cold gt 0.0)
print, "   Cold Gas (should be same as Star-forming) "
print, "      N_cold (rho_cold > 0)= ", n_elements(idx)
;print, "      Mass= ", total(m(idx))
print, "      rho_cold max/min: ", max(rho_cold), min(rho_cold(idx)), min(rho_cold)

coldf= fload_gas_coldfraction(1)
print, "      cold gas fraction max/min: ", max(coldf), min(coldf)

ColdGasMass= m*coldf
print, "      Mass= ", total(ColdGasMass)

; volume_cold / volume_hot
coldgasff= fload_gas_coldgasfillingfactor(1,TurbulentFactor=TurbulentFactor)
print, "      coldgasff max/min: ", max(coldgasff), min(coldgasff)

; volume_cold / volume_total
coldgasvf= fload_gas_coldgasvolumefraction(1,TurbulentFactor=TurbulentFactor)
print, "      coldgasvf max/min: ", max(coldgasvf), min(coldgasvf)
print, " "


; BAD clumpsize calculation.  The hsml is not a good 
; indicator of the size of the SPH particle
;
;clumpsize= hsml(idx) * (coldgasvf(idx)^(1./3.))
;print, "      clumpsize max/min: ", max(clumpsize), min(clumpsize)
;print, "                   mean: ", mean(clumpsize)
;

; cold clump size
newclumpsize= ( (4.*!PI/ 3.) * rho_cold(idx) / ColdGasMass(idx))^(-1./3.)
clumpsize= newclumpsize
print, "(new) clumpsize max/min: ", max(clumpsize), min(clumpsize)
print, "                   mean: ", mean(clumpsize)

print, " "
print, " "




; -----------------------------------------------
;
;    Set the variables that go into file

	N_ColdGas= n_elements(idx)
	ColdGas_x= x(idx) / h
	ColdGas_y= y(idx) / h
	ColdGas_z= z(idx) / h
	ColdGas_Density= rho_cold(idx) * h * h
	ColdGas_Size= clumpsize / h
	SphGas_Hsml= hsml(idx) / h
	SphGas_Z= zmets(idx) / 0.02


; -----------------------------------------------

print, "xxxxxxxxxxxxxxxxxx"
print, "  Check SPH info  "
print, "  "
cgmass= 4.*!PI/3. * ColdGas_Size * ColdGas_Size * ColdGas_Size * ColdGas_Density
print, "total mass= ", total(cgmass)
print, "  "
print, "xxxxxxxxxxxxxxxxxx"

; -----------------------------------------------
;
;    Add any lower density gas

;add_lowdensity_gas= 1
add_lowdensity_gas= 0

if add_lowdensity_gas eq 1 then begin

	; low density gas
	gastemp= fload_gas_temperature(1)   ; in K
	idx= where((rho lt rho_crit) and (gastemp lt 4.0e+4))

	print, "-----------------------------"
	print, "Adding Cold Low-density gas"
	print, " "
	print, "  N= ", n_elements(idx)
        N_ColdGas= N_ColdGas + n_elements(idx)
	r= sqrt(x(idx)*x(idx) + y(idx)*y(idx) + z(idx)*z(idx))
	print, "  radius max/min: ", max(r), min(r)
        ColdGas_x= [ColdGas_x, x(idx) / h]
        ColdGas_y= [ColdGas_y, y(idx) / h]
        ColdGas_z= [ColdGas_z, z(idx) / h]
        ColdGas_Density= [ColdGas_Density, rho(idx) * h * h]     ; use regular rho
	print, "  max/min density: ", max(rho(idx)), min(rho(idx))
	clumpsize= ( (4.*!PI/ 3.) * rho(idx) / m(idx))^(-1./3.)
	print, "  clumpsize max/min: ", max(clumpsize), min(clumpsize)
        ColdGas_Size= [ColdGas_Size, clumpsize / h]
        SphGas_Hsml= [SphGas_Hsml, hsml(idx) / h]
	print, "-----------------------------"


	print, "xxxxxxxxxxxxxxxxxx"
	print, "  Check SPH info  "
	print, "  "
	cgmass= 4.*!PI/3. * ColdGas_Size * ColdGas_Size * ColdGas_Size * ColdGas_Density
	print, "total mass= ", total(cgmass)
	print, "  "
	print, "xxxxxxxxxxxxxxxxxx"

endif

; -----------------------------------------------

;    Order these by radius

	radius= sqrt(ColdGas_x * ColdGas_x + $
			ColdGas_y * ColdGas_y + $
			ColdGas_z * ColdGas_z)
	sortidx= sort(radius)
	ColdGas_x= ColdGas_x(sortidx)
	ColdGas_y= ColdGas_y(sortidx)
	ColdGas_z= ColdGas_z(sortidx)
	ColdGas_Density= ColdGas_Density(sortidx)
	ColdGas_Size= ColdGas_Size(sortidx)
	SphGas_Hsml= SphGas_Hsml(sortidx)
	SphGas_Z= SphGas_Z(sortidx)


	print, "min radius= ", min(radius)
	print, "max radius= ", max(radius)

; -----------------------------------------------


; do some checks
; ----------------
;print, "Density_TotalGas   max/min= ", max(Density_TotalGas), min(Density_TotalGas)
;print, "         N= ", n_elements(Density_TotalGas)
;idx=where(Density_TotalGas gt 0.0)
;print, " N_nonzero= ", n_elements(idx)
;idx=where(Density_TotalGas lt 1.0e-10)
;Density_TotalGas(idx)= 1.0d-12

;print, "Cube Gas Mass= ", total(Grid_Volume*Density_TotalGas)


; Make the Density Log
; ---------------------
UnitDensity_in_cgs = 6.76991d-22
PROTONMASS = 1.6726d-24
factor= UnitDensity_in_cgs / PROTONMASS
print, "Density factor= ", factor
ColdGas_Density= alog10(factor*ColdGas_Density)     ; convert to cm-3  (is mu right?)
	


; Make the Clump Size Log
; ------------------------
ColdGas_Size= alog10(ColdGas_Size)



; -------------------------------------------
;
;   Write SPH Gas Info to file
;

snaplbl= strcompress(string(snapnum),/remove_all)
gridfile=extra+'griddata_'+snaplbl+'_sph.txt'
openw, 1, gridfile, ERROR=err

printf, 1, '#  SPH cold clumps for: '+frun+', snap= '+snaplbl+' (T= '+fload_timelbl(h,3,/noteq)+' Gyr)'
printf, 1, "#             N_clumps: "+string(N_ColdGas)+"    TurbFactor: "+string(TurbulentFactor)
printf, 1, "#               center: ", center_bh
printf, 1, "#     x          y          z   log(density) log(size) "
printf, 1, "#   (kpc)      (kpc)      (kpc)     (cm-3)      (kpc)  "

for i=0L,N_ColdGas-1 do begin
        printf, 1, FORMAT='(6(F9.4,"  "))', $
		ColdGas_x[i], $
		ColdGas_y[i], $
		ColdGas_z[i], $
		ColdGas_Density[i], $
		ColdGas_Size[i], $
		SphGas_Z[i]
		;SphGas_Hsml[i]
endfor

close, 1




; -------------------------------------------

end





; =============================================================================
; =============================================================================
; =============================================================================
; =============================================================================






pro sukanya_distributed_sources_grid, junk




;-------------------------------------------
;
;   Define Grid
;  --------------

;  spherical
; ------------
spherical_grid= 1
;spherical_grid= 0
if spherical_grid eq 1 then begin

	; need the following order

	; grid_*
	; ----------
	N_r= 200L     &  grid_r= fltarr(N_r)
	N_theta= 59L  &  grid_theta= fltarr(N_theta)
	N_phi= 118L   &  grid_phi= fltarr(N_phi)
	; cent_*
	; ---------
	N_r= 199L     &  grid_cen_r= fltarr(N_r)
	N_theta= 58L  &  grid_cen_theta= fltarr(N_theta)
	N_phi= 117L   &  grid_cen_phi= fltarr(N_phi)
	; test_*
	; ---------
	;N_r= 23L     &  r_data= fltarr(N_r)
	;N_theta= 9L  &  theta_data= fltarr(N_theta)
	;N_phi= 25L   &  phi_data= fltarr(N_phi)


	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_r.dat"
	readf, 1, grid_r
	close, 1
	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_r.dat"
	readf, 1, grid_cen_r
	close, 1

	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_theta.dat"
	readf, 1, grid_theta
	close, 1
	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_theta.dat"
	readf, 1, grid_cen_theta
	close, 1

	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_phi.dat"
	readf, 1, grid_phi
	close, 1
	openr, 1, "/home2/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_phi.dat"
	readf, 1, grid_cen_phi
	close, 1


	r_len= max(grid_r)
	total_Volume= (4./3.)*!PI*r_len*r_len*r_len
	print, "total_Volume= ", total_Volume

	Ngrid= N_r*N_theta*N_phi

	Grid= fltarr(3,Ngrid)
	Grid_rthetaphi= fltarr(3,Ngrid)
	Grid_Volume= fltarr(Ngrid)

	n_i= long(0)
	for i=0,N_r-1 do begin
	  for j=0,N_theta-1 do begin
	    for k=0,N_phi-1 do begin

		R= grid_cen_r[i]
		theta= grid_cen_theta[j]
		phi= grid_cen_phi[k]
		Grid_rthetaphi[0,n_i]= R
		Grid_rthetaphi[1,n_i]= theta
		Grid_rthetaphi[2,n_i]= phi

		x= R*sin(theta)*cos(phi)
		y= R*sin(theta)*sin(phi)
		z= R*cos(theta)
		Grid[0,n_i]= x
		Grid[1,n_i]= y
		Grid[2,n_i]= z

		;Grid_Volume[n_i]= R * 10^(d_log_r) * d_phi * d_cos_theta
		;Grid_Volume[n_i]= R * 10^(d_log_r) * d_phi * sin(theta) * d_theta
		;Grid_Volume[n_i]= R * R * d_log_r * d_phi * sin(theta) * d_theta    ;  dr= r * d_logr
		delta_r= grid_r[i+1] - grid_r[i]
		delta_theta= grid_theta[j+1] - grid_theta[j]
		delta_phi= grid_phi[k+1] - grid_phi[k]
		Grid_Volume[n_i]= R * R * delta_r * delta_phi * sin(theta) * delta_theta

		n_i= n_i + 1
	    endfor
	  endfor
	endfor
endif



;  Volume check
; --------------
sphere_area= 0.0
for j=0,N_theta-1 do begin
  for k=0,N_phi-1 do begin
        delta_theta= grid_theta[j+1] - grid_theta[j]
        delta_phi= grid_phi[k+1] - grid_phi[k]
        sphere_area= sphere_area + delta_phi * sin(grid_theta[j]) * delta_theta
  endfor
endfor
print, "sphere area= ", sphere_area
print, '4*PI= ', 4.0*!PI

print, "Total Grid_Volume= ", total(Grid_Volume)


;----------------------------------------------

; define lums

;bololum= 10^(11.5226)   & snapnum= 10    ; vc3vc3e_2/snapshot_010
;bololum= 10^(10.9360)   & snapnum= 18    ; vc3vc3e_2/snapshot_018
;bololum= 10^(11.2716)   & snapnum= 32    ; vc3vc3e_2/snapshot_032
bololum= 10^(11.4220)   & snapnum= 52    ; vc3vc3e_2/snapshot_052
;bololum= 10^(11.3048)   & snapnum= 54    ; vc3vc3e_2/snapshot_054
;bololum= 10^(11.1312)   & snapnum= 72    ; vc3vc3e_2/snapshot_072



; number of sources (only 1 or 6 right now)
;source_num= 1
source_num= 6


;-------------------------------------------
;

StellarLum= fltarr(Ngrid)


if source_num eq 1 then begin
	lum_radius= 0.050956954
	;lum_radius= 0.10181827
	;lum_radius= 0.20344542
	lum_theta= 1.5963744
	lum_phi= 0.026851219

	idx=where((Grid_rthetaphi[0,*] gt (lum_radius-0.001)) and (Grid_rthetaphi[0,*] lt (lum_radius+0.001)) and $
			(Grid_rthetaphi[1,*] gt (lum_theta-0.01)) and (Grid_rthetaphi[1,*] lt (lum_theta+0.01)) and $
			(Grid_rthetaphi[2,*] gt (lum_phi-0.01)) and (Grid_rthetaphi[2,*] lt (lum_phi+0.01)))

	if idx(0) ne -1 then begin
		print, "found cell number", idx," to luminate"
		StellarLum(idx)= alog10(bololum)
	endif else begin
		print, "PROBLEM: can't find cell"
	endelse
endif


if source_num eq 6 then begin

	;lum_radius= 0.050956954
	;lum_radius= 0.10181827
	lum_radius= 0.20344542
        lum_theta= 0.069211335
        lum_phi= 0.026851219

        idx=where((Grid_rthetaphi[0,*] gt (lum_radius-0.001)) and (Grid_rthetaphi[0,*] lt (lum_radius+0.001)) and $
                        (Grid_rthetaphi[1,*] gt (lum_theta-0.01)) and (Grid_rthetaphi[1,*] lt (lum_theta+0.01)) and $
                        (Grid_rthetaphi[2,*] gt (lum_phi-0.01)) and (Grid_rthetaphi[2,*] lt (lum_phi+0.01)))

        if idx(0) ne -1 then StellarLum(idx)= alog10(bololum / source_num)
	print, '1: found cell number ', idx,'  to luminate'

        lum_theta= 1.5963744

        idx=where((Grid_rthetaphi[0,*] gt (lum_radius-0.001)) and (Grid_rthetaphi[0,*] lt (lum_radius+0.001)) and $
                        (Grid_rthetaphi[1,*] gt (lum_theta-0.01)) and (Grid_rthetaphi[1,*] lt (lum_theta+0.01)) and $
                        (Grid_rthetaphi[2,*] gt (lum_phi-0.01)) and (Grid_rthetaphi[2,*] lt (lum_phi+0.01)))

        if idx(0) ne -1 then StellarLum(idx)= alog10(bololum / source_num)
	print, '2: found cell number ', idx,'  to luminate'

        lum_phi= 1.5842219

        idx=where((Grid_rthetaphi[0,*] gt (lum_radius-0.001)) and (Grid_rthetaphi[0,*] lt (lum_radius+0.001)) and $
                        (Grid_rthetaphi[1,*] gt (lum_theta-0.01)) and (Grid_rthetaphi[1,*] lt (lum_theta+0.01)) and $
                        (Grid_rthetaphi[2,*] gt (lum_phi-0.01)) and (Grid_rthetaphi[2,*] lt (lum_phi+0.01)))

        if idx(0) ne -1 then StellarLum(idx)= alog10(bololum / source_num)
	print, '3: found cell number ', idx,'  to luminate'

        lum_phi= 3.1415927

        idx=where((Grid_rthetaphi[0,*] gt (lum_radius-0.001)) and (Grid_rthetaphi[0,*] lt (lum_radius+0.001)) and $
                        (Grid_rthetaphi[1,*] gt (lum_theta-0.01)) and (Grid_rthetaphi[1,*] lt (lum_theta+0.01)) and $
                        (Grid_rthetaphi[2,*] gt (lum_phi-0.01)) and (Grid_rthetaphi[2,*] lt (lum_phi+0.01)))

        if idx(0) ne -1 then StellarLum(idx)= alog10(bololum / source_num)
	print, '4: found cell number ', idx,'  to luminate'

        lum_phi= 4.6989633

        idx=where((Grid_rthetaphi[0,*] gt (lum_radius-0.001)) and (Grid_rthetaphi[0,*] lt (lum_radius+0.001)) and $
                        (Grid_rthetaphi[1,*] gt (lum_theta-0.01)) and (Grid_rthetaphi[1,*] lt (lum_theta+0.01)) and $
                        (Grid_rthetaphi[2,*] gt (lum_phi-0.01)) and (Grid_rthetaphi[2,*] lt (lum_phi+0.01)))

        if idx(0) ne -1 then StellarLum(idx)= alog10(bololum / source_num)
	print, '5: found cell number ', idx,'  to luminate'

        lum_theta= 3.0723813
        lum_phi= 0.026851219

        idx=where((Grid_rthetaphi[0,*] gt (lum_radius-0.001)) and (Grid_rthetaphi[0,*] lt (lum_radius+0.001)) and $
                        (Grid_rthetaphi[1,*] gt (lum_theta-0.01)) and (Grid_rthetaphi[1,*] lt (lum_theta+0.01)) and $
                        (Grid_rthetaphi[2,*] gt (lum_phi-0.01)) and (Grid_rthetaphi[2,*] lt (lum_phi+0.01)))

        if idx(0) ne -1 then StellarLum(idx)= alog10(bololum / source_num)
	print, '6: found cell number ', idx,'  to luminate'

endif



;----------------------------------------------
;
;   Write StellarLum to file
;


;slgridfile='grid_stellarlum.txt'
slgridfile='griddata_'+strcompress(string(snapnum),/remove_all)+'_slum.txt'
openw, 1, slgridfile, ERROR=err

for i=0L,Ngrid-1 do begin
        printf, 1, FORMAT='(F7.3)', StellarLum[i]
endfor


close, 1

;----------------------------------------------





end








; =============================================================================






; =============================================================================






