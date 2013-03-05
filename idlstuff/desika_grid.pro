;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   Calculated grid quantities
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------




;
;
;
;-------------------------------
pro do_desika, frun

   ;frun="/raid4/tcox/ds/vc3vc3e_2"
   ;frun="/raid4/tcox/vc3vc3e_no"
   ;frun="/raid4/tcox/vc3vc3h_2"
   ;frun="/raid4/tcox/sbw/sb10"
   ;frun="/raid4/tcox/sbw/sb10BH"
   ;frun="/raid4/tcox/bs/b3e_2"
   ;frun="/raid4/tcox/ds/d5e2_2"


   if not keyword_set(frun) then begin
	print, "  "
	print, " do_desika, frun"
	print, "  "
	return
   endif

   starti= 0L
   ;starti= 4L
   ;starti= 14L
   ;starti= 39
   ;starti= 40
   ;starti=64L
   ;starti= 61
   ;starti= 64

   ;endi= 31
   ;endi= 52
   ;endi= 72
   ;endi= 99
   ;endi= 107

   snapinterval= 1
   ;snapinterval= 10

   spawn, "/bin/ls "+frun+"/snap* | wc ",result                                    
   endi=long(result[0])-1                                                  

   spawn, "mkdir "+frun+"/desikagrids_new"

   print, "frun= ", frun
   print, "starti= ", starti
   print, "endi= ", endi

   for i=starti,endi do begin

     thisi= i
     if ((thisi mod snapinterval) eq 0) then begin
	print, "do_one_snapshot, "+frun+", "+strcompress(string(thisi),/remove_all)
	do_one_snapshot, frun, thisi
     endif

   endfor

end



;
;
;
;-------------------------------
pro do_one_snapshot, frun, snapshotnum

     ;do_bh_center= 0
     do_bh_center= 1

     if do_bh_center eq 1 then begin
        ok=fload_snapshot_bh(frun, snapshotnum, /nopot_in_snap)
	if ok ne 0 then return   ; make sure we found a snapshot

        bhid= fload_blackhole_id(1)
        ;bhid= bhid[0]
        ;bhid= bhid[1]
        ;bhid= 200001L
        ;bhid= 280002L   ; used for z3/b4e
        ;bhid= 400002L   ; used for ds/vc3vc3e_2
        bhid= long(max(bhid))   ; uses minimum bhid - then you don't need to know the BH ID
        center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)

        print, "Blackhole ID: ", bhid
        print, "Blackhole center: ", center_bh

        grid_v1, frun, snapshotnum, center= center_bh, /loadedsnap
     endif else begin
        grid_v1, frun, snapshotnum, /use_calc_center
     endelse


end










; =========================
;
;      Desika Grid (v0)
;
;
;
;   OLD OLD OLD OLD   - do not use!!
;
; =========================
pro grid_v0, frun, snapnum, center=center, $
				use_calc_center=use_calc_center, $
				stellar_ages=stellar_ages

spawn, "date", result
print, "date=", result

if not keyword_set(frun) then begin
	print, "  "
	print, " grid_v0, frun, snapnum"
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
print, " ======================="
print, "    Density  "
print, " ======================="


center= orig_center
print, "using center: ", center
x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
m= fload_gas_mass(1) / 0.7
ids= fload_gas_id(1)

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

	print, "-----------------------------"
	print, " "
	print, "Total according to particles = ", total(m)
	idx= where((abs(x) lt 0.5*side_len) and (abs(y) lt 0.5*side_len) and (abs(z) lt 0.5*side_len))
	if idx(0) ne -1 then mass_inside_pts= total(m(idx)) else mass_inside_pts= 0.0
	print, "Total in box according to particles = ", mass_inside_pts



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
        print, "Total according to Grid ", total(cellVolume*Density)
        print, " "
        perlbl= strmid(string(100.0 * (total(cellVolume*Density)-mass_inside_pts)/mass_inside_pts),5)
        print, "   error: "+perlbl, " %"
        print, " "

	print, "gas_density min/max= ", min(Density), max(Density)
	idx= where(Density le 0.0)
	print, "idx info"
	help, idx

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
print, "using center: ", center
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
print, "using center: ", center
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

	StellarMass_1= fltarr(Ngrid)
	StellarMass_2= fltarr(Ngrid)
	StellarMass_3= fltarr(Ngrid)
	StellarMass_4= fltarr(Ngrid)

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

		print, "-----------------------------"
		print, " "
		minside= m(idx)
		print, "Total according to particles = ", total(minside)
		inside_idx= where((abs(x(idx)) lt 0.5*side_len) and (abs(y(idx)) lt 0.5*side_len) and (abs(z(idx)) lt 0.5*side_len))
		if inside_idx(0) ne -1 then mass_inside_pts= total(minside(inside_idx)) else mass_inside_pts= 0.0
		print, "Total in box according to particles = ", mass_inside_pts


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

	           print, "Total according to Grid ", total(cellVolume*StellarDensity)
	           print, " "
	           perlbl= strmid(string(100.0 * (total(cellVolume*StellarDensity)-mass_inside_pts)/mass_inside_pts),5)
	           print, "   error: "+perlbl, " %"
	           print, " "

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


		print, "-----------------------------"
		print, " "
		minside= m(idx)
		print, "Total according to particles = ", total(minside)
		inside_idx= where((abs(x(idx)) lt 0.5*side_len) and (abs(y(idx)) lt 0.5*side_len) and (abs(z(idx)) lt 0.5*side_len))
		if inside_idx(0) ne -1 then mass_inside_pts= total(minside(inside_idx)) else mass_inside_pts= 0.0
		print, "Total in box according to particles = ", mass_inside_pts


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

	           print, "Total according to Grid ", total(cellVolume*StellarDensity)
	           print, " "
	           perlbl= strmid(string(100.0 * (total(cellVolume*StellarDensity)-mass_inside_pts)/mass_inside_pts),5)
	           print, "   error: "+perlbl, " %"
	           print, " "

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

		print, "-----------------------------"
		print, " "
		minside= m(idx)
		print, "Total according to particles = ", total(minside)
		inside_idx= where((abs(x(idx)) lt 0.5*side_len) and (abs(y(idx)) lt 0.5*side_len) and (abs(z(idx)) lt 0.5*side_len))
		if inside_idx(0) ne -1 then mass_inside_pts= total(minside(inside_idx)) else mass_inside_pts= 0.0
		print, "Total in box according to particles = ", mass_inside_pts


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

	           print, "Total according to Grid ", total(cellVolume*StellarDensity)
	           print, " "
	           perlbl= strmid(string(100.0 * (total(cellVolume*StellarDensity)-mass_inside_pts)/mass_inside_pts),5)
	           print, "   error: "+perlbl, " %"
	           print, " "

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

		   print, "-----------------------------"
		   print, " "
		   minside= m(idx)
		   print, "Total according to particles = ", total(minside)
		   inside_idx= where((abs(x(idx)) lt 0.5*side_len) and (abs(y(idx)) lt 0.5*side_len) and (abs(z(idx)) lt 0.5*side_len))
		   if inside_idx(0) ne -1 then mass_inside_pts= total(minside(inside_idx)) else mass_inside_pts= 0.0
		   print, "Total in box according to particles = ", mass_inside_pts


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

	           print, "Total according to Grid ", total(cellVolume*StellarDensity)
	           print, " "
	           perlbl= strmid(string(100.0 * (total(cellVolume*StellarDensity)-mass_inside_pts)/mass_inside_pts),5)
	           print, "   error: "+perlbl, " %"
	           print, " "

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

		print, "-----------------------------"
		print, " "
		minside= m(idx)
		print, "Total according to particles = ", total(minside)
		inside_idx= where((abs(x(idx)) lt 0.5*side_len) and (abs(y(idx)) lt 0.5*side_len) and (abs(z(idx)) lt 0.5*side_len))
		if inside_idx(0) ne -1 then mass_inside_pts= total(minside(inside_idx)) else mass_inside_pts= 0.0
		print, "Total in box according to particles = ", mass_inside_pts


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
	        print, "Total according to Grid ", total(cellVolume*StellarDensity)
	        print, " "
	        perlbl= strmid(string(100.0 * (total(cellVolume*StellarDensity)-mass_inside_pts)/mass_inside_pts),5)
	        print, "   error: "+perlbl, " %"
	        print, " "

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



centerlbl= string(center[0])+','+string(center[1])+','+string(center[2])



;=========================================================================================
;----------------------------------------------
;
;   Write to Grid file
;

;gridfile=frun+'/desika_grids/griddata_'+strcompress(string(snapnum),/remove_all)+'.txt'
gridfile=frun+'/desikagrids/griddata_'+strcompress(string(snapnum),/remove_all)+'.txt'
;gridfile=frun+'/desikagrids/griddata_test_'+strcompress(string(snapnum),/remove_all)+'.txt'
;gridfile=frun+'/griddata_'+strcompress(string(snapnum),/remove_all)+'.txt'
openw, 1, gridfile, ERROR=err

; actually print the rows
if keyword_set(stellar_ages) then begin
	printf, 1, '# '+griddim+sidellbl+' frun='+frun+'  snapnum=',strcompress(string(snapnum),/remove_all)+'  center='+centerlbl
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

spawn, "date", result
print, "date=", result

;----------------------------------------------

print, "done"

end










;==========================================================================================
;
;
;      =========================
;     
;           Desika Grid (v1)
;     
;      =========================
;
;
;==========================================================================================
pro grid_v1, frun, snapnum, center=center, $
			use_calc_center=use_calc_center, $
			loadedsnap=loadedsnap


spawn, "date", result
print, "date=", result

if not keyword_set(frun) then begin
	print, "  "
	print, " grid_v1, frun, snapnum, center=center"
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

print, " "
print, "n_side= ", n_side
print, "side_len= ", side_len
print, "Total grid volume= ", side_len*side_len*side_len
print, "cellVolume= ", cellVolume
print, " "

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

if not keyword_set(loadedsnap) then begin
	;ok=fload_snapshot_bh(frun,lastsnapshot)         
	;ok=fload_snapshot_bh(frun,snapnum)         
	ok=fload_snapshot_bh(frun,snapnum, /nopot_in_snap)         
endif


if keyword_set(use_calc_center) then center=fload_center_alreadycomp(1)

orig_center= center




;-------------------------------------------
;
;   Density
;

mass= fload_gas_mass(1) / 0.7

; option 1
;print, " ======================="
;print, "    Density (trial 1)   "
;print, " ======================="
;m= mass & ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
;gas_masses= fload_total_ongrid(m, ng, g, cv, sl, c)
;gas_masses= fload_total_ongrid_1(m, ng, g, cv, sl, c)
;Density= alog10(gas_masses * (1.0/cellVolume) * 674.59)     ; factor converts to cm-3
;
; these were bad
;

; option 2
print, " ======================="
print, "    Density (trial 2)   "
print, " ======================="
center= orig_center
print, "using center: ", center
m= mass & ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
gas_density= fload_density_ongrid(m, ng, g, cv, sl, c)   ; this one produces smoother dist. which Desika prefers, this is also the same as the "old" version
;gas_density= fload_density_ongrid_1(m, ng, g, cv, sl, c)
print, "gas_density min/max= ", min(gas_density), max(gas_density)
idx= where(gas_density le 0.0)
print, "idx info"
help, idx
if idx(0) ne -1 then gas_density(idx) = min(gas_density(where(gas_density gt 0.0))) / 100000.0
Density= alog10(gas_density * 674.59)     ; factor converts to cm-3




;-------------------------------------------
;
;   Temperature
;

;temp_k= fload_gas_temperature_multi(1,/cold)
temp_k= fload_gas_temperature(1)

print, " ======================="
print, "       Temperature   "
print, " ======================="
;m= temp_k & ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
;gas_temperatures= fload_density_ongrid(m, ng, g, cv, sl, c)
;print, "gas_temperatures min/max= ", min(gas_temperatures), max(gas_temperatures)
;Temp= alog10(gas_temperatures)


center= orig_center
print, "using center: ", center
m= temp_k & ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
gas_temperatures= fload_average_ongrid(m, ng, g, cv, sl, c)
print, "gas_temperatures min/max= ", min(gas_temperatures), max(gas_temperatures)
;Temp= alog10(gas_temperatures) ; do this later so that we can properly convert to degrees





;-------------------------------------------
;
;  Compute Multiphase Quantities
;

simtime = float(fload_time(1))

simrho= gas_density
simne= gas_density*0.0 + 1.0

; need these to match the run
; ----------------------------
hubbleparam= 0.7
tstar = 4.5
tSN   = float(3.0e+8)
fEVP  = 3000.0
fSN   = 0.1
comoving = 0.0

print, "t_star=    ",tstar
print, "tSN=       ",tSN
print, "fEVP=      ",fEVP
print, "fSN=       ",fSN
;print, "coming=    ",comoving


cold_mass_fraction   = fltarr(Ngrid)     ;mass fraction
cold_volume_fraction = fltarr(Ngrid)     ;volume fraction
hot_temperature      = fltarr(Ngrid)     ;hot phase temp


; Brant's external code to calculate the cooling and
; return the calculated mass fraction
; ----------------------------------------------------
;S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/ComputeMultiphase/multiphase.so', $
;S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/ComputeMultiphase/mphase', $
S = CALL_EXTERNAL('/n/home/tcox/C-Routines_for_IDL/ComputeMultiphase/mphase', $
	'multiphase_info', $
	Ngrid, $
	tstar, $
	tSN, $
	fEVP, $
	fSN, $
	hubbleparam, $
	simtime, $
	comoving, $
	simrho, $
	simne, $
	cold_mass_fraction, $
	cold_volume_fraction, $
	hot_temperature )




; this is the x from SH03  (eq. 17, rho_cold/rho)
x= cold_mass_fraction
;print, "x (mass fraction)= ", x[100:103]
ColdFraction= x

        ;  Contants
        ; ---------------------------
        HYDROGEN_MASSFRAC= 0.76
        GAMMA_MINUS1= (5./3.) - 1.0
        BOLTZMANN=   1.3806d-16
        PROTONMASS=  1.6726d-24
        UnitMass_in_g= 1.989d+43
        UnitDensity_in_cgs = 6.76991d-22
        UnitEnergy_in_cgs = 1.989d+53

        meanweight = 4. / (1. + 3. * HYDROGEN_MASSFRAC)
	temp_factor= 1. / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * UnitMass_in_g / UnitEnergy_in_cgs

	TurbulentFactor= 1.0


        ;  Cold Phase Temperature
        ; ---------------------------
        TempClouds    = 1000.0           ; in K
        u_cold = TempClouds * temp_factor 


        ;  calculate grid "u" and "rho"
        ; -------------------------------
	u= gas_temperatures * temp_factor
	rho= gas_density


        ; fix cold "temperature"
        u_cold= u*0.0 + u_cold*TurbulentFactor
        idx= where(u_cold gt 0.75*u)
        if idx(0) ne -1 then u_cold(idx)= 0.75*u(idx)

        ; we assume that the total energy (u) is a 
        ; mass weighted average, so
        ;  rho*u =  rho_hot*u_hot + rho_cold*u_cold
        ;or
        ; u = (rho_hot/rho)u_hot  +  (rho_cold/rho)u_cold
        ; u = (1-x)u_hot  +  (x)u_cold
        ;u_hot = (u - x*u_cold) / (1.0-x)
        u_hot = (u - x*u_cold) / (1.0-x)
	;u_hot = hot_temperature

        ; ********
        ; not sure if this is correct, but
        ; add this to trap for all cold particles
        idx= where(u_hot lt 1.1*u_cold) 
        if idx(0) ne -1 then u_hot(idx) = 1.1*u_cold(idx)
        ; ********

        density_ratio = u_cold / u_hot
        ;print, "rho_hot / rho_cold = ", density_ratio[10318:10320]

        ; this is volume_cold / volume_hot 
        volume_filling_factor= (x/(1.0-x)) * density_ratio
        ;print, "cold volume_filling_factor= ", volume_filling_factor[10318:10320]

        ; what's original rho?
        ;print, "rho= ", rho[10318:10320]

        ; ***  tj's method ***
        ;rho_cold= rho * u / (2.0 * u_cold)
        rho_cold= rho * u /  u_cold
        rho_hot = rho_cold * density_ratio
        ;print, "---- TJ's methods ----"
        ;print, "rho_cold= ", rho_cold[10318:10320]
        ;print, "rho_hot= ", rho_hot[10318:10320]
        ;print, " "

        ; ***  phil's method  ***
        xc0= (1.0 + volume_filling_factor) * x / (1.e-10 + volume_filling_factor)
        xh0= xc0 * density_ratio
        rho_cold= rho *  xc0
        rho_hot= rho * xh0
        ;print, "---- Phil's methods ----"
        ;print, "rho_cold= ", rho_cold[10318:10320]
        ;print, "rho_hot= ", rho_hot[10318:10320]
        ;print, " "

        ;print, "Pressures:"
        ;print, "tot:   ",  rho[10318:10320]*u[10318:10320]
        ;print, "cold:  ",  rho_cold[10318:10320]*u_cold
        ;print, "hot:   ",  rho_hot[10318:10320]*u_hot[10318:10320]





;xxxxxxxx

Density_Cold = rho_cold
Temp_Cold = fltarr(Ngrid)


; these default to straight particle values
Density_Hot = gas_density
Temp_Hot= gas_temperatures

;
; these are actually set above.
;
; save the "totals"
;Density = gas_density
Temp= gas_temperatures
;
;idx= where(Density gt 0.0)
;if idx(0) ne -1 then Density(idx)= alog10(Density(idx) * 674.59)
idx= where(Temp gt 0.0)
if idx(0) ne -1 then Temp(idx)= alog10(Temp(idx) / temp_factor)



idx= where(rho_cold gt 0.0)
;idx= [-1]
if idx(0) ne -1 then begin

	Density_Cold(idx)= rho_cold(idx)
	Temp_Cold(idx)= 3.0                 ; immediately put this in Log K

	Density_Hot(idx)= rho_hot(idx)
	Temp_Hot(idx)= hot_temperature(idx)

endif

idx= where(Density_Cold gt 0.0)
if idx(0) ne -1 then Density_Cold(idx)= alog10(Density_Cold(idx) * 674.59)     ; factor converts to cm-3

idx= where(Density_Hot gt 0.0)
if idx(0) ne -1 then Density_Hot(idx)= alog10(Density_Hot(idx) * 674.59)     ; factor converts to cm-3

idx= where(Temp_Hot gt 0.0)
Temp_Hot(idx)= alog10(Temp_Hot(idx) / temp_factor)           ;  now in K





;-------------------------------------------
;
;   Velocity
;

print, " ======================="
print, "       Velocity   "
print, " ======================="
center= orig_center
print, "using center: ", center
ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
Velocity= fload_velocity_ongrid(ng, g, cv, sl, c, comvel=comvel)





;-------------------------------------------
;
;   SFR
;

thissfr= fload_gas_sfr(1)
print, " ======================="
print, "      SFR   "
print, "      Total (Msolar/Yr) = ", total(thissfr)
print, " ======================="

center= orig_center
print, "using center: ", center
ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
;SFR= fload_total_ongrid(thissfr, ng, g, cv, sl, c)
;SFR= fload_total_ongrid_1(thissfr, ng, g, cv, sl, c)
SFR_density= fload_density_ongrid(thissfr, ng, g, cv, sl, c, /with_mass_only)
;SFR_density= fload_density_ongrid_1(thissfr, ng, g, cv, sl, c)
SFR= cellVolume * SFR_density

idx= where(SFR gt 0.0)
SFR(idx) = alog10(SFR(idx))



;-------------------------------------------
;
;   Stellar Mass
;


DesNgb=96L
Hmax= 100.0

StellarMass_0= fltarr(Ngrid)
StellarMass_1= fltarr(Ngrid)
StellarMass_2= fltarr(Ngrid)
StellarMass_3= fltarr(Ngrid)
StellarMass_4= fltarr(Ngrid)

N= fload_npart(2) + fload_npart(3) + fload_npart(4)

if(n_elements(N) gt 0) then begin


	; option 2
	mass= fload_allstars_mass(1) / 0.7
	age= (fload_time(1) - fload_allstars_age(1)) / 0.7


	; All Stellar Ages
	;-------------------
	print, " ======================="
	print, "   Total Stellar Mass   "
	print, " ======================="
	this_mass= mass
	center= orig_center
	print, "using center: ", center
	ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
	;
	; trial 0 
	;print, "stellar mass: fload_total_ongrid"
	;stellar_mass= fload_total_ongrid(this_mass, ng, g, cv, sl, c, /stars)
	;
	; trial 1
	;stellar_mass= fload_total_ongrid_1(this_mass, ng, g, cv, sl, c, /stars)
	;
	; trial 2
	stellar_density= fload_density_ongrid(this_mass, ng, g, cv, sl, c, /stars, /with_mass_only)
	;stellar_density= fload_density_ongrid_1(this_mass, ng, g, cv, sl, c, /stars)
	;stellar_density= fload_density_ongrid_1(this_mass, ng, g, cv, sl, c, /stars, /with_mass_only)
	stellar_mass= cellVolume * stellar_density
	;
	StellarMass_0= stellar_mass
	idx= where(StellarMass_0 lt 0.0)
	if idx(0) ne -1 then print, "WARNING: there are negative masses in StellarMass_0"
	idx= where(StellarMass_0 gt 0.0)
	if idx(0) ne -1 then StellarMass_0(idx)= 10.0 + alog10(StellarMass_0(idx))



	; Stellar Age 1
        ;-------------------
	print, " ======================="
	print, "   New Stellar Mass   "
	print, " ======================="
	this_mass= mass
	center= orig_center
	print, "using center: ", center
	idx= where(age ge 0.010)           ; set non- below 10 Myr to 0
	if n_elements(idx) gt 1 then this_mass(idx)= 0.0
	ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
	;stellar_mass= fload_total_ongrid(this_mass, ng, g, cv, sl, c, /stars)
	;stellar_mass= fload_total_ongrid_1(this_mass, ng, g, cv, sl, c, /stars)
	stellar_density= fload_density_ongrid(this_mass, ng, g, cv, sl, c, /stars, /with_mass_only)
	;stellar_density= fload_density_ongrid_1(this_mass, ng, g, cv, sl, c, /stars)
	;stellar_density= fload_density_ongrid_1(this_mass, ng, g, cv, sl, c, /stars, /with_mass_only)
	stellar_mass= cellVolume * stellar_density
	StellarMass_1= stellar_mass
	idx= where(StellarMass_1 lt 0.0)
	if idx(0) ne -1 then print, "WARNING: there are negative masses in StellarMass_1"
	idx= where(StellarMass_1 gt 0.0)
	if idx(0) ne -1 then StellarMass_1(idx)= 10.0 + alog10(StellarMass_1(idx))



	; Stellar Age 2
        ;-------------------
	print, " ======================="
	print, "   Young Stellar Mass   "
	print, " ======================="
	this_mass= mass
	center= orig_center
	print, "using center: ", center
        idx= where((age lt 0.010) or (age ge 0.100))
	if n_elements(idx) gt 1 then this_mass(idx)= 0.0
	ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
	;stellar_mass= fload_total_ongrid(this_mass, ng, g, cv, sl, c, /stars)
	;stellar_mass= fload_total_ongrid_1(this_mass, ng, g, cv, sl, c, /stars)
	;stellar_density= fload_density_ongrid(this_mass, ng, g, cv, sl, c, /stars)
	stellar_density= fload_density_ongrid(this_mass, ng, g, cv, sl, c, /stars, /with_mass_only)
	;stellar_density= fload_density_ongrid_1(this_mass, ng, g, cv, sl, c, /stars)
	;stellar_density= fload_density_ongrid_1(this_mass, ng, g, cv, sl, c, /stars, /with_mass_only)
	stellar_mass= cellVolume * stellar_density
	StellarMass_2= stellar_mass
	idx= where(StellarMass_2 lt 0.0)
	if idx(0) ne -1 then print, "WARNING: there are negative masses in StellarMass_2"
	idx= where(StellarMass_2 gt 0.0)
	if idx(0) ne -1 then StellarMass_2(idx)= 10.0 + alog10(StellarMass_2(idx))




	; Stellar Age 3
	;-------------------
	print, " ======================="
	print, "   Mod Stellar Mass   "
	print, " ======================="
	this_mass= mass
	center= orig_center
	print, "using center: ", center
	idx= where((age lt 0.100) or (age ge 1.0))
	if n_elements(idx) gt 1 then this_mass(idx)= 0.0
	ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
	;stellar_mass= fload_total_ongrid(this_mass, ng, g, cv, sl, c, /stars)
	;stellar_mass= fload_total_ongrid_1(this_mass, ng, g, cv, sl, c, /stars)
	stellar_density= fload_density_ongrid(this_mass, ng, g, cv, sl, c, /stars, /with_mass_only)
	;stellar_density= fload_density_ongrid_1(this_mass, ng, g, cv, sl, c, /stars)
	;stellar_density= fload_density_ongrid_1(this_mass, ng, g, cv, sl, c, /stars, /with_mass_only)
	stellar_mass= cellVolume * stellar_density
	StellarMass_3= stellar_mass
	idx= where(StellarMass_3 lt 0.0)
	if idx(0) ne -1 then print, "WARNING: there are negative masses in StellarMass_3"
	idx= where(StellarMass_3 gt 0.0)
	if idx(0) ne -1 then StellarMass_3(idx)= 10.0 + alog10(StellarMass_3(idx))




	; Stellar Age 4
	;-------------------
	print, " ======================="
	print, "   Old Stellar Mass   "
	print, " ======================="
	this_mass= mass
	center= orig_center
	print, "using center: ", center
        idx= where(age lt 1.0)
	if n_elements(idx) gt 1 then this_mass(idx)= 0.0
	ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
	;stellar_mass= fload_total_ongrid(this_mass, ng, g, cv, sl, c, /stars)
	;stellar_mass= fload_total_ongrid_1(this_mass, ng, g, cv, sl, c, /stars)
	stellar_density= fload_density_ongrid(this_mass, ng, g, cv, sl, c, /stars, /with_mass_only)
	;stellar_density= fload_density_ongrid_1(this_mass, ng, g, cv, sl, c, /stars)
	;stellar_density= fload_density_ongrid_1(this_mass, ng, g, cv, sl, c, /stars, /with_mass_only)
	stellar_mass= cellVolume * stellar_density
	StellarMass_4= stellar_mass
	idx= where(StellarMass_4 lt 0.0)
	if idx(0) ne -1 then print, "WARNING: there are negative masses in StellarMass_4"
	idx= where(StellarMass_4 gt 0.0)
	if idx(0) ne -1 then StellarMass_4(idx)= 10.0 + alog10(StellarMass_4(idx))


endif




;=========================================================================================
;



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




centerlbl= string(center[0])+','+string(center[1])+','+string(center[2])

comvellbl= string(comvel[0])+','+string(comvel[1])+','+string(comvel[2])


;=========================================================================================
;
;   Write to Grid file
;
;



;gridfile=frun+'/desikagrids/griddata_new_'+strcompress(string(snapnum),/remove_all)+'.txt'
gridfile=frun+'/desikagrids_new/griddata_'+strcompress(string(snapnum),/remove_all)+'.txt'
;gridfile=frun+'/desikagrids_new/griddata_test_'+strcompress(string(snapnum),/remove_all)+'.txt'
print, "opening: "+gridfile
openw, 1, gridfile, ERROR=err

; actually print the rows
printf, 1, '# '+griddim+sidellbl+' frun='+frun+'  snapnum='+strcompress(string(snapnum),/remove_all)+'  center='+centerlbl
printf, 1, '#                                                               comvel='+comvellbl
printf, 1, "#                                                    lg         lg   lg Cold  lg Cold            lg Hot   lg Hot     lg             Stellar Mass (log Msolar)"
printf, 1, "#   x     y      z       Velx      Vely      Velz  Density     Temp  Density     Temp    Cold   Density     Temp    SFR               Age     >10Myr  >100Myr"
printf, 1, "#(kpc) (kpc)  (kpc)    (km/s)    (km/s)    (km/s)   (cm-3)      (K)   (cm-3)      (K) Fraction   (cm-3)      (K)(Msun/Yr)    Total  0-10Myr <=100Myr   <=1Gyr    >1Gyr"
;printf, 1, "#                                                  lg Cold  lg Cold            lg Hot   lg Hot                    Stellar Mass (log Msolar)"
;printf, 1, "#   x     y      z       Velx      Vely      Velz  Density     Temp    Cold   Density     Temp    SFR               Age     >10Myr  >100Myr"
;printf, 1, "#(kpc) (kpc)  (kpc)    (km/s)    (km/s)    (km/s)   (cm-3)      (K) Fraction   (cm-3)      (K)(Msun/Yr)    Total  0-10Myr <=100Myr   <=1Gyr    >1Gyr"
;printf, 1, "#                                                     Cold     Cold      Hot      Hot         Stellar Mass (Msolar)"
;printf, 1, "#   x     y      z       Velx      Vely      Velz  Density     Temp  Density     Temp      Age   >10Myr  >100Myr"
;printf, 1, "#(kpc) (kpc)  (kpc)    (km/s)    (km/s)    (km/s)   (cm-3)      (K)   (cm-3)      (K)  0-10Myr <=100Myr   <=1Gyr    >1Gyr"

for i=0L,Ngrid-1 do begin
	printf, 1, FORMAT='(3(F5.2,"  "),3(F8.2,"  "),13(F7.3,"  "))', $
		Grid[0,i],Grid[1,i],Grid[2,i], $
		Velocity[0,i],Velocity[1,i],Velocity[2,i], $
		Density[i], $
		Temp[i], $
		Density_Cold[i], $
		Temp_Cold[i], $
		ColdFraction[i], $
		Density_Hot[i], $
		Temp_Hot[i], $
		SFR[i], $
		StellarMass_0[i], $
		StellarMass_1[i], $
		StellarMass_2[i], $
		StellarMass_3[i], $
		StellarMass_4[i]
endfor
close, 1

;----------------------------------------------

spawn, "date", result
print, "date=", result
print, "done"

;----------------------------------------------

end




;==========================================================================================




pro grid_v1_info, junk



   read_new_desika_grid, "/raid4/tcox/ds/vc3vc3e_2/griddata_new_52.txt", Grid, Velocity, $
			Density, Temp, Density_Cold, Temp_Cold, ColdFraction, Density_Hot, Temp_Hot, SFR, StellarMass

   n_3= (size(Grid))[2]
   n= long(n_3^(1/3.))

   dg= abs(Grid[2,0]-Grid[2,1])
   cellVolume= dg*dg*dg


   print, " "
   print, "Total SFR (Msolar/Yr)= ", total(SFR)


   print, " "
   print, "Cold Gas Mass (Msolar)= ", total(10^(Density_Cold)) / 674.59
   print, "Hot  Gas Mass (Msolar)= ", total(10^(Density_Hot)) / 674.59

   print, " "
   print, "Stellar Mass Total (Msolar)= ", total(10^(StellarMass[0,*]))
   print, "Stellar Mass Age0 (Msolar)=  ", total(10^(StellarMass[1,*]))
   print, "Stellar Mass Age1 (Msolar)=  ", total(10^(StellarMass[2,*]))
   print, "Stellar Mass Age2 (Msolar)=  ", total(10^(StellarMass[3,*]))
   print, "Stellar Mass Age3 (Msolar)=  ", total(10^(StellarMass[4,*]))

   print, " "

end





;==========================================================================================





; =========================
;
;    Test Density Calc.
;
; =========================
pro grid_test_density, frun, snapnum, center=center, $
				use_calc_center=use_calc_center, $
				stellar_ages=stellar_ages, $
				loadedsnap=loadedsnap


if not keyword_set(frun) then begin
	print, "  "
	print, " grid_test_density, frun, snapnum"
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
;n_side= 200L
gridlbl= strcompress(string(n_side),/remove_all)
griddim= gridlbl+'x'+gridlbl+'x'+gridlbl
;side_len= 70.0
;side_len= 12.0
side_len= 8.0
sidellbl= strcompress(string(side_len),/remove_all)
sidellbl= ' box len='+sidellbl
element_len= side_len/n_side
cellVolume= element_len*element_len*element_len

print, " "
print, "n_side= ", n_side
print, "side_len= ", side_len
print, "Total grid volume= ", side_len*side_len*side_len
print, "cellVolume= ", cellVolume
print, " "

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
if not keyword_set(loadedsnap) then ok=fload_snapshot_bh(frun,snapnum, /nopot_in_snap)         


if keyword_set(use_calc_center) then center=fload_center_alreadycomp(1)

orig_center= center

;-------------------------------------------
;
;   Density
;


;
; original method to calculate density along
; the grid
;
;
;method_one= 1
method_one= 0
if method_one eq 1 then begin
	center= orig_center
	x= fload_gas_xyz('x',center=center) / 0.7
	y= fload_gas_xyz('y',center=center) / 0.7
	z= fload_gas_xyz('z',center=center) / 0.7
	m= fload_gas_mass(1) / 0.7
	ids= fload_gas_id(1)

	N= long(n_elements(x))

	print, "Total Gas Mass= ", total(m)
	idx= where((abs(x) lt 0.5*side_len) and (abs(y) lt 0.5*side_len) and (abs(z) lt 0.5*side_len))
	if idx(0) ne -1 then mass_inside_pts= total(m(idx)) else mass_inside_pts= 0.0


	DesNgb=96L
	;DesNgb=26L
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
		print, "Within Cube"
		print, "==========="
		print, "Grid Gas Mass= ", total(cellVolume*Density)
		print, "Particles Gas Mass within cube= ", mass_inside_pts
		print, "   error: ", 100.0 * (total(cellVolume*Density)-mass_inside_pts)/mass_inside_pts, " %"
		print, " "

	

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

endif



;
;  new method which is basically identical to the above only it
;  is wrapped inside another function
;

method_two= 1
;method_two= 0
if method_two eq 1 then begin

	mass= fload_gas_mass(1)
	;mass= fload_gas_sfr(1)

	m= mass & ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
	gas_density= fload_density_ongrid(m, ng, g, cv, sl, c)
	Density= alog10(gas_density * 674.59)     ; factor converts to cm-3

endif



;
;  alternate method that calculated totals, and then divides by volume
;

;method_three= 1
method_three= 0
if method_three eq 1 then begin

	mass= fload_gas_mass(1)
	;mass= fload_gas_sfr(1)

	;m= mass & ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
	;gas_masses= fload_total_ongrid(m, ng, g, cv, sl, c)
	;Density= alog10(gas_masses * (1.0/cellVolume) * 674.59)     ; factor converts to cm-3

        m= mass & ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
        gas_density= fload_density_ongrid_1(m, ng, g, cv, sl, c)
        Density= alog10(gas_density * 674.59)     ; factor converts to cm-3

endif


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

gridfile=frun+'/griddata_testdensity_'+strcompress(string(snapnum),/remove_all)+'.txt'
openw, 1, gridfile, ERROR=err

; actually print the rows
	printf, 1, '# '+griddim+sidellbl+' frun='+frun+'  snapnum=',strcompress(string(snapnum),/remove_all)
	printf, 1, "# "
	printf, 1, "# "
	printf, 1, "#   x     y      z    Density"
	printf, 1, "#(kpc) (kpc)  (kpc)    (cm-3)"

	for i=0L,Ngrid-1 do begin
		printf, 1, FORMAT='(3(F5.2,"  "),1(F8.3,"  "))', $
			Grid[0,i],Grid[1,i],Grid[2,i], $
			Density[i]
	endfor

close, 1

;----------------------------------------------

end









;==========================================================================================







; =========================
;
;    Test Mass Calc.
;
; =========================
pro grid_test_mass, frun, snapnum, center=center, $
				use_calc_center=use_calc_center, $
				stellar_ages=stellar_ages, $
				loadedsnap=loadedsnap

frun="/raid4/tcox/ds/vc3vc3e_2"
if not keyword_set(snapnum) then snapnum= 52

if not keyword_set(frun) then begin
	print, "  "
	print, " grid_test_mass, frun, snapnum"
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
;n_side= 200L
gridlbl= strcompress(string(n_side),/remove_all)
griddim= gridlbl+'x'+gridlbl+'x'+gridlbl
;side_len= 70.0
;side_len= 12.0
side_len= 8.0
sidellbl= strcompress(string(side_len),/remove_all)
sidellbl= ' box len='+sidellbl
element_len= side_len/n_side
cellVolume= element_len*element_len*element_len

print, " "
print, "n_side= ", n_side
print, "side_len= ", side_len
print, "Total grid volume= ", side_len*side_len*side_len
print, "cellVolume= ", cellVolume
print, " "

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
if not keyword_set(loadedsnap) then ok=fload_snapshot_bh(frun,snapnum, /nopot_in_snap)         


if keyword_set(use_calc_center) then center=fload_center_alreadycomp(1)

orig_center= center



;-------------------------------------------
;
;   Density
;


;
; original method to calculate density along
; the grid
;
;
;method_one= 1
method_one= 0
if method_one eq 1 then begin
	center= orig_center
	x= fload_gas_xyz('x',center=center) / 0.7
	y= fload_gas_xyz('y',center=center) / 0.7
	z= fload_gas_xyz('z',center=center) / 0.7
	hsml= fload_gas_hsml(1)/0.7
	m= fload_gas_mass(1) / 0.7
	ids= fload_gas_id(1)

	N= long(n_elements(x))

	print, "Total Gas Mass= ", total(m)
	idx= where((abs(x) lt 0.5*side_len) and (abs(y) lt 0.5*side_len) and (abs(z) lt 0.5*side_len))
	if idx(0) ne -1 then mass_inside_pts= total(m(idx)) else mass_inside_pts= 0.0


	DesNgb=96L
	Hmax= 100.0

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

        	Mass_onGrid= fltarr(Ngrid)      ;  full velocity field

        	S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Cube/Mass/mass.so', $
                	'mass', $
			N, $
			Ngrid, $
                	Coord, $
                	Masses, $
			DesNgb, $
			Hmax, $
                	Mass_onGrid)

		; Make the Density Log
		; ---------------------
		print, " "
		print, "Within Cube"
		print, "==========="
		print, "Grid Gas Mass= ", total(Mass_onGrid)
		print, "Particles Gas Mass within cube= ", mass_inside_pts
		print, "   error: ", 100.0 * (total(Mass_onGrid)-mass_inside_pts)/mass_inside_pts, " %"
		print, " "

	
		; check for zeros
		idx= where(Mass_onGrid le 0.0)
		;if idx(0) ne -1 then Mass_onGrid(idx)= min(Mass_onGrid(where(Mass_onGrid gt 0.0)))

		; log and cm-3
		Density= alog10(674.59*Mass_onGrid/cellVolume)     ; convert to cm-3  (is mu right?)

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
endif




;
;  new method which is basically identical to the above only it
;  is wrapped inside another function
;

method_two= 1
;method_two= 0
if method_two eq 1 then begin

        mass= fload_gas_mass(1)
        ;mass= fload_gas_sfr(1)

        m= mass & ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
        gas_mass= fload_total_ongrid(m, ng, g, cv, sl, c)
        Density= alog10(674.59 * gas_mass / cellVolume)     ; factor converts to cm-3

endif



;
;  alternate method that calculated totals, and then divides by volume
;

method_three= 1
;method_three= 0
if method_three eq 1 then begin

        mass= fload_gas_mass(1)
        ;mass= fload_gas_sfr(1)

        m= mass & ng= Ngrid & g= Grid & cv= cellVolume & sl= side_len & c= center
        gas_masses= fload_total_ongrid_1(m, ng, g, cv, sl, c)
        Density= alog10(gas_masses * (1.0/cellVolume) * 674.59)     ; factor converts to cm-3

endif



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

gridfile=frun+'/griddata_testmass_'+strcompress(string(snapnum),/remove_all)+'.txt'
openw, 1, gridfile, ERROR=err

; actually print the rows
	printf, 1, '# '+griddim+sidellbl+' frun='+frun+'  snapnum=',strcompress(string(snapnum),/remove_all)
	printf, 1, "# "
	printf, 1, "# "
	printf, 1, "#   x     y      z    Density"
	printf, 1, "#(kpc) (kpc)  (kpc)    (cm-3)"

	for i=0L,Ngrid-1 do begin
		printf, 1, FORMAT='(3(F5.2,"  "),1(F8.3,"  "))', $
			Grid[0,i],Grid[1,i],Grid[2,i], $
			Density[i]
	endfor

close, 1

;----------------------------------------------

end





;============================================================================






; =========================
;
;    Test Velocity Calc.
;
; =========================
pro grid_test_velocity, frun, snapnum, center=center, $
				use_calc_center=use_calc_center, $
				loadedsnap=loadedsnap

;frun="/raid4/tcox/z3/b6h"
if not keyword_set(snapnum) then snapnum= 39

if not keyword_set(frun) then begin
	print, "  "
	print, " grid_test_mass, frun, snapnum"
	print, "  "
	return
endif

if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(center) then center=[0,0,0]



;-------------------------------------------
;
;   Define Grid
;
;n_side= 25L
n_side= 50L
;n_side= 75L
;n_side= 100L
;n_side= 200L
gridlbl= strcompress(string(n_side),/remove_all)
griddim= gridlbl+'x'+gridlbl+'x'+gridlbl
;side_len= 70.0
;side_len= 12.0
side_len= 8.0
sidellbl= strcompress(string(side_len),/remove_all)
sidellbl= ' box len='+sidellbl
element_len= side_len/n_side
cellVolume= element_len*element_len*element_len

print, " "
print, "n_side= ", n_side
print, "side_len= ", side_len
print, "Total grid volume= ", side_len*side_len*side_len
print, "cellVolume= ", cellVolume
print, " "

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
if not keyword_set(loadedsnap) then ok=fload_snapshot_bh(frun,snapnum, /nopot_in_snap)         


if keyword_set(use_calc_center) then center=fload_center_alreadycomp(1)


	bhid= fload_blackhole_id(1)
        ;bhid= bhid[0]
        ;bhid= bhid[1]
        ;bhid= 200001L
        ;bhid= 280002L   ; used for z3/b4e
        ;bhid= 400002L   ; used for ds/vc3vc3e_2
        bhid= long(max(bhid))   ; uses minimum bhid - then you don't need to know the BH ID
        center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)

        print, "Blackhole ID: ", bhid
        print, "Blackhole center: ", center_bh

	center= center_bh

	; use black hole velocity
	comvel= fload_blackhole_v('xyz', center=center, idtofollow=bhid)
	print, "#1= using comvel: ", comvel

	comvel= fload_all_comvel(1,center=center)
	print, "#2= using comvel: ", comvel

	comvel= fload_all_comvel(1,center=center,justcenter=2.0)
	print, "#3= using comvel: ", comvel

	comvel= fload_comvel(1,center=center,rfact=0.001)
	print, "#4= using comvel: ", comvel


orig_center= center




;-------------------------------------------
;
;   Velocity
;


center= orig_center
print, "using center: ", center
x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
m= fload_gas_mass(1) / 0.7
ids= fload_gas_id(1)

center= orig_center
print, "using center: ", center
print, "using comvel: ", comvel
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





idx=where(finite(Grid) eq 0, n_i)
if idx(0) ne -1 then print, 'Grid has '+strcompress(string(n_i),/remove_all)+' bad numbers'

idx=where(finite(Velocity) eq 0, n_i)
if idx(0) ne -1 then print, 'Velocity has '+strcompress(string(n_i),/remove_all)+' bad numbers'



;==============================================
;----------------------------------------------
;
;   Write to Grid file
;

;gridfile=frun+'/griddata_testvelocity_'+strcompress(string(snapnum),/remove_all)+'.txt'
gridfile=frun+'/desikagrids_new/griddata_testvelocity_'+strcompress(string(snapnum),/remove_all)+'.txt'
print, "opening: "+gridfile
openw, 1, gridfile, ERROR=err

; actually print the rows
	printf, 1, '# '+griddim+sidellbl+' frun='+frun+'  snapnum=',strcompress(string(snapnum),/remove_all)
	printf, 1, "# "
	printf, 1, "# "
	printf, 1, "#   x     y      z    Velocity"
	printf, 1, "#(kpc) (kpc)  (kpc)   (km s-1)"

	for i=0L,Ngrid-1 do begin
		printf, 1, FORMAT='(3(F5.2,"  "),3(F8.3,"  "))', $
			Grid[0,i],Grid[1,i],Grid[2,i], $
			Velocity[0,i],Velocity[1,i],Velocity[2,i]
	endfor

close, 1

;----------------------------------------------

end







;==========================================================================================
;
;
;   load the new grid
;
;

pro read_new_desika_grid, gridfile, Grid, Velocity, $
		Density, Temp, Density_Cold, Temp_Cold, ColdFraction, Density_Hot, Temp_Hot, SFR, StellarMass

if not keyword_set(gridfile) then begin
        print, " "
        print, " PROBLEM: gridfile not set "
        print, " "
        return
endif else begin
        print, "opening: "+gridfile
endelse

spawn, "wc "+gridfile, result
lines=long(result)
datalines=lines(0)-5

Grid= fltarr(3,datalines)
Velocity= fltarr(3,datalines)
Density= fltarr(datalines)
Temp= fltarr(datalines)
Density_Cold= fltarr(datalines)
Temp_Cold= fltarr(datalines)
ColdFraction= fltarr(datalines)
Density_Hot= fltarr(datalines)
Temp_Hot= fltarr(datalines)
SFR= fltarr(datalines)
StellarMass= fltarr(5,datalines)

openr, 1, gridfile
junk=''

; read header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read data
for i=0L,lines(0)-6 do begin
        readf, 1, junk
        ;print, junk
        tempjunk= strsplit(junk,/extract,count=count)
        Grid(0,i)= float(tempjunk(0))
        Grid(1,i)= float(tempjunk(1))
        Grid(2,i)= float(tempjunk(2))
        Velocity(0,i)= float(tempjunk(3))
        Velocity(1,i)= float(tempjunk(4))
        Velocity(2,i)= float(tempjunk(5))
        Density(i)= float(tempjunk(6))
        Temp(i)= float(tempjunk(7))
        Density_Cold(i)= float(tempjunk(8))
        Temp_Cold(i)= float(tempjunk(9))
        ColdFraction(i)= float(tempjunk(10))
        Density_Hot(i)= float(tempjunk(11))
        Temp_Hot(i)= float(tempjunk(12))
        SFR(i)= float(tempjunk(13))
        StellarMass(0,i)= float(tempjunk(14))
        StellarMass(1,i)= float(tempjunk(15))
        StellarMass(2,i)= float(tempjunk(16))
        StellarMass(3,i)= float(tempjunk(17))
        StellarMass(4,i)= float(tempjunk(18))
endfor

close, 1


end






;==========================================================================================
;
;
;   load the density of something on
;   the grid.
;
;   two versions of this:
;        _0, i.e. nothing    goes through the grid and finds the particles nearby
;                            (quicker, but has a larger error)
;        _1      goes through the particles and puts them on the grid
;


;
;   _0
;   
function fload_density_ongrid, input_x, input_Ngrid, input_Grid, input_cellVolume, input_side_len, input_center, $
						stars=stars, input_idx=input_idx, with_mass_only=with_mass_only

;
m= input_x
center= input_center

x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
hsml= fload_gas_hsml(1)
ids= fload_gas_id(1)


if keyword_set(stars) then begin
        x= fload_allstars_xyz('x',center=center) / 0.7
        y= fload_allstars_xyz('y',center=center) / 0.7
        z= fload_allstars_xyz('z',center=center) / 0.7
        hsml= x*0.0 + 0.2
        ids= fload_allstars_ids(1)
endif


; make a cut
;
if keyword_set(input_idx) then begin
	m= m(input_idx)
	x= x(input_idx)
	y= y(input_idx)
	z= z(input_idx)
	hsml= hsml(input_idx)
	idx= idx(input_idx)
endif


; select only particles with mass
;
if keyword_set(with_mass_only) then begin
	midx= where(m gt 0.0)
	if midx(0) ne -1 then begin
		m= m(midx)
		x= x(midx)
		y= y(midx)
		z= z(midx)
		hsml= hsml(midx)
		ids= ids(midx)
	endif else begin
		x= [-1]
	endelse
endif



N= long(n_elements(x))

print, "-----------------------------"
print, " "
print, "Total according to particles = ", total(m)
idx= where((abs(x) lt 0.5*input_side_len) and (abs(y) lt 0.5*input_side_len) and (abs(z) lt 0.5*input_side_len))
if idx(0) ne -1 then mass_inside_pts= total(m(idx)) else mass_inside_pts= 0.0
print, "Total in box according to particles = ", mass_inside_pts


DesNgb=96L
;Hmax= input_side_len/3.0
Hmax= input_side_len

if(N gt 1) then begin

	Ncoord= max([N,input_Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:input_Ngrid-1)= input_Grid(0,*)
	Coord(5,0:input_Ngrid-1)= input_Grid(1,*)
	Coord(6,0:input_Ngrid-1)= input_Grid(2,*)
	Masses(0:N-1)= m(*)

        Density= fltarr(input_Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/n/home/tcox/C-Routines_for_IDL/Cube/Density/density.so', $
                'density', $
		N, $
		input_Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Density)

	; 
	; ---------------------
	print, "Total according to Grid ", total(input_cellVolume*Density)
	print, " "
	perlbl= strmid(string(100.0 * (total(input_cellVolume*Density)-mass_inside_pts)/mass_inside_pts),5)
	print, "   error: "+perlbl, " %"
	print, " "

endif else begin
        print,' Cannot find more then 1 particle'
	print,'N= ', N
        Density= fltarr(input_Ngrid)
endelse

print, " "
print, "-----------------------------"

return, Density


end






;
;   _1
;   
function fload_density_ongrid_1, input_x, input_Ngrid, input_Grid, input_cellVolume, input_side_len, input_center, $
						stars=stars, input_idx=input_idx, with_mass_only=with_mass_only

;
m= input_x
center= input_center

x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
hsml= fload_gas_hsml(1)
ids= fload_gas_id(1)

if keyword_set(stars) then begin
        x= fload_allstars_xyz('x',center=center) / 0.7
        y= fload_allstars_xyz('y',center=center) / 0.7
        z= fload_allstars_xyz('z',center=center) / 0.7
        hsml= x*0.0 + 0.2
        ids= fload_allstars_ids(1)
endif


; make a cut
;
if keyword_set(input_idx) then begin
	m= m(input_idx)
	x= x(input_idx)
	y= y(input_idx)
	z= z(input_idx)
	hsml= hsml(input_idx)
	ids= ids(input_idx)
endif


; select only particles with mass
;
if keyword_set(with_mass_only) then begin
	midx= where(m gt 0.0)
	m= m(midx)
	x= x(midx)
	y= y(midx)
	z= z(midx)
	hsml= hsml(midx)
	ids= ids(midx)
endif



N= long(n_elements(x))

print, "-----------------------------"
print, " "
print, "Total according to particles = ", total(m)
idx= where((abs(x) lt 0.5*input_side_len) and (abs(y) lt 0.5*input_side_len) and (abs(z) lt 0.5*input_side_len))
if idx(0) ne -1 then mass_inside_pts= total(m(idx)) else mass_inside_pts= 0.0
print, "Total in box according to particles = ", mass_inside_pts
print, "center= ", input_center


DesNgb=96L
Hmax= 100.0

if(N gt 1) then begin

	Ncoord= max([N,input_Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:input_Ngrid-1)= input_Grid(0,*)
	Coord(5,0:input_Ngrid-1)= input_Grid(1,*)
	Coord(6,0:input_Ngrid-1)= input_Grid(2,*)
	Masses(0:N-1)= m(*)

        Density= fltarr(input_Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/n/home/tcox/C-Routines_for_IDL/Cube/Density2/density.so', $
                'density', $
		N, $
		input_Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Density)

	; 
	; ---------------------
	print, "Total according to Grid ", total(input_cellVolume*Density)
	print, " "
	perlbl= strmid(string(100.0 * (total(input_cellVolume*Density)-mass_inside_pts)/mass_inside_pts),5)
	print, "   error: "+perlbl, " %"
	print, " "

endif else begin
        print,' Cannot find anything'
	print,'N= ', N
        Density= fltarr(input_Ngrid)
endelse

print, " "
print, "-----------------------------"

return, Density


end







;==========================================================================================
;
;
;   load the total amount of something on
;   the grid.
;
;
;   two versions of this:
;        _0  does a smoothing (max r_smooth is hsml)
;        _1  each particle should be put on only one grid cell (doesn't use hsml)
;
;     (both go through the particle and put on the grid)


;
;   _0
;   

function fload_total_ongrid, input_x, input_Ngrid, input_Grid, input_cellVolume, input_side_len, input_center, $
				stars=stars

;
m= input_x
center= input_center

x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
hsml= fload_gas_hsml(1)
ids= fload_gas_id(1)

if keyword_set(stars) then begin
	x= fload_allstars_xyz('x',center=center) / 0.7
	y= fload_allstars_xyz('y',center=center) / 0.7
	z= fload_allstars_xyz('z',center=center) / 0.7
	hsml= x*0.0 + 0.2
	ids= fload_allstars_ids(1)
endif


N= long(n_elements(x))

print, "-----------------------------"
print, " "
print, "Total according to particles = ", total(m)
idx= where((abs(x) lt 0.5*input_side_len) and (abs(y) lt 0.5*input_side_len) and (abs(z) lt 0.5*input_side_len))
if idx(0) ne -1 then mass_inside_pts= total(m(idx)) else mass_inside_pts= 0.0
print, "Total in box according to particles = ", mass_inside_pts
print, "center= ", input_center


DesNgb=96L
Hmax= 100.0

if(N gt 1) then begin

	Ncoord= max([N,input_Ngrid])

	Coord=fltarr(8,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:input_Ngrid-1)= input_Grid(0,*)
	Coord(5,0:input_Ngrid-1)= input_Grid(1,*)
	Coord(6,0:input_Ngrid-1)= input_Grid(2,*)
	Coord(7,0:N-1)= hsml(*)
	Masses(0:N-1)= m(*)

        TotalOnGrid= fltarr(input_Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/n/home/tcox/C-Routines_for_IDL/Cube/Mass/mass.so', $
                'mass', $
		N, $
		input_Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                TotalOnGrid)

	; Make the Density Log
	; ---------------------
	print, "Total according to Grid ", total(TotalOnGrid)
	print, " "
	perlbl= strmid(string(100.0 * (total(TotalOnGrid)-mass_inside_pts)/mass_inside_pts),5)
	print, "   error: "+perlbl, " %"
	print, " "

endif else begin
        print,' Cannot find anything'
	print,'N= ', N
        TotalOnGrid= fltarr(input_Ngrid)
endelse

print, " "
print, "-----------------------------"

return, TotalOnGrid


end







;
;   _1
;   

function fload_total_ongrid_1, input_x, input_Ngrid, input_Grid, input_cellVolume, input_side_len, input_center, $
				stars=stars

;
m= input_x
center= input_center

x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
hsml= fload_gas_hsml(1)
ids= fload_gas_id(1)

if keyword_set(stars) then begin
	x= fload_allstars_xyz('x',center=center) / 0.7
	y= fload_allstars_xyz('y',center=center) / 0.7
	z= fload_allstars_xyz('z',center=center) / 0.7
	hsml= x*0.0 + 0.2
	ids= fload_allstars_ids(1)
endif


N= long(n_elements(x))

print, "-----------------------------"
print, " "
print, "Total according to particles = ", total(m)
idx= where((abs(x) lt 0.5*input_side_len) and (abs(y) lt 0.5*input_side_len) and (abs(z) lt 0.5*input_side_len))
if idx(0) ne -1 then mass_inside_pts= total(m(idx)) else mass_inside_pts= 0.0
print, "Total in box according to particles = ", mass_inside_pts


DesNgb=5L
;Hmax= 1.0
;Hmax= 1.0
ds= abs(input_Grid[2,0]-input_Grid[2,1])
Hmax= 1.2 * sqrt((ds/2.)*(ds/2.)*3.)      ; all particles should be within this distance

if(N gt 1) then begin

	Ncoord= max([N,input_Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:input_Ngrid-1)= input_Grid(0,*)
	Coord(5,0:input_Ngrid-1)= input_Grid(1,*)
	Coord(6,0:input_Ngrid-1)= input_Grid(2,*)
	Masses(0:N-1)= m(*)

        TotalOnGrid= fltarr(input_Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/n/home/tcox/C-Routines_for_IDL/Cube/Mass2/mass.so', $
                'mass', $
		N, $
		input_Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                TotalOnGrid)

	; Make the Density Log
	; ---------------------
	print, " "
	print, "Total according to Grid ", total(TotalOnGrid)
	print, " "
	perlbl= strmid(string(100.0 * (total(TotalOnGrid)-mass_inside_pts)/mass_inside_pts),5)
	print, "   error: "+perlbl, " %"
	print, " "

endif else begin
        print,' Cannot find anything - need more then one particle'
	print,'N= ', N
        TotalOnGrid= fltarr(input_Ngrid)
endelse

print, " "
print, "-----------------------------"

return, TotalOnGrid


end




;
;   _2
;   

function fload_total_ongrid_2, input_x, input_Ngrid, input_Grid, input_cellVolume, input_side_len, input_center, $
				stars=stars, input_idx=input_idx

;
m= input_x
center= input_center

x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
hsml= fload_gas_hsml(1)
ids= fload_gas_id(1)

if keyword_set(stars) then begin
	x= fload_allstars_xyz('x',center=center) / 0.7
	y= fload_allstars_xyz('y',center=center) / 0.7
	z= fload_allstars_xyz('z',center=center) / 0.7
	hsml= x*0.0 + 0.2
	ids= fload_allstars_ids(1)
endif


N= long(n_elements(x))

print, "-----------------------------"
print, " "
print, "Total according to particles = ", total(m)
idx= where((abs(x) lt 0.5*input_side_len) and (abs(y) lt 0.5*input_side_len) and (abs(z) lt 0.5*input_side_len))
if idx(0) ne -1 then mass_inside_pts= total(m(idx)) else mass_inside_pts= 0.0
print, "Total in box according to particles = ", mass_inside_pts


DesNgb=5L
;Hmax= 1.0
;Hmax= 1.0
ds= abs(input_Grid[2,0]-input_Grid[2,1])
Hmax= 1.2 * sqrt((ds/2.)*(ds/2.)*3.)      ; all particles should be within this distance

if(N gt 1) then begin

	Ncoord= max([N,input_Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:input_Ngrid-1)= input_Grid(0,*)
	Coord(5,0:input_Ngrid-1)= input_Grid(1,*)
	Coord(6,0:input_Ngrid-1)= input_Grid(2,*)
	Masses(0:N-1)= m(*)

        TotalOnGrid= fltarr(input_Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/n/home/tcox/C-Routines_for_IDL/Cube/Mass2/mass.so', $
                'mass', $
		N, $
		input_Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                TotalOnGrid)

	; Make the Density Log
	; ---------------------
	print, " "
	print, "Total according to Grid ", total(TotalOnGrid)
	print, " "
	perlbl= strmid(string(100.0 * (total(TotalOnGrid)-mass_inside_pts)/mass_inside_pts),5)
	print, "   error: "+perlbl, " %"
	print, " "

endif else begin
        print,' Cannot find anything'
	print,' N= ',N
        TotalOnGrid= fltarr(input_Ngrid)
endelse

print, " "
print, "-----------------------------"

return, TotalOnGrid


end






;==========================================================================================
;
;
;   load the velocity on
;   the grid.
;
;

function fload_velocity_ongrid, input_Ngrid, input_Grid, input_cellVolume, input_side_len, input_center, $
                                stars=stars, comvel=comvel

print, "=============================="
print, "entering fload_velocity_ongrid"

;
center= input_center
print, "using center: ", center
;
;comvel= fload_all_comvel(1,center=center)
;print, "comvel #1: ", comvel
comvel= fload_all_comvel(1,center=center,justcenter=2.0)
print, "comvel #2: ", comvel
;comvel= fload_comvel(1,center=center,rfact=0.001)
;print, "comvel #3: ", comvel





x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
m= fload_gas_mass(1)
hsml= fload_gas_hsml(1)
ids= fload_gas_id(1)
;comvel= fload_all_comvel(1,center=center)
vx= fload_gas_v('x',comvel=comvel,center=center)
vy= fload_gas_v('y',comvel=comvel,center=center)
vz= fload_gas_v('z',comvel=comvel,center=center)


if keyword_set(stars) then begin
        x= fload_allstars_xyz('x',center=center) / 0.7
        y= fload_allstars_xyz('y',center=center) / 0.7
        z= fload_allstars_xyz('z',center=center) / 0.7
	m= fload_allstars_mass(1)
        hsml= x*0.0 + 0.2
        ids= fload_allstars_ids(1)
	vx= fload_allstars_v('x',comvel=comvel,center=center)
	vy= fload_allstars_v('y',comvel=comvel,center=center)
	vz= fload_allstars_v('z',comvel=comvel,center=center)

endif


N= long(n_elements(vx))


DesNgb=96L
Hmax= 100.0

if(N gt 1) then begin

	Ncoord= max([N,input_Ngrid])

	Coord=fltarr(10,Ncoord)
	Masses= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= vx(*)
	Coord(4,0:N-1)= vy(*)
	Coord(5,0:N-1)= vz(*)
	Coord(6,0:N-1)= ids(*)
	Coord(7,0:input_Ngrid-1)= input_Grid(0,*)
	Coord(8,0:input_Ngrid-1)= input_Grid(1,*)
	Coord(9,0:input_Ngrid-1)= input_Grid(2,*)
	Masses(0:N-1)= m(*)

        Velocity= fltarr(3,input_Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/n/home/tcox/C-Routines_for_IDL/Cube/Velocity/velocity.so', $
                'velocity', $
		N, $
		input_Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Velocity)
	print, " "


endif else begin
        print,'No stars, no problem.'
	print, 'N= ', N
        Velocity= fltarr(3, input_Ngrid)
endelse


print, " "
print, "-----------------------------"

return, Velocity


end




;==========================================================================================
;
;
;   load average quantities on
;   the grid.
;
;


function fload_average_ongrid, input_x, input_Ngrid, input_Grid, input_cellVolume, input_side_len, input_center, $
						star=star

;
qtoavg= input_x
center= input_center

x= fload_gas_xyz('x',center=center) / 0.7
y= fload_gas_xyz('y',center=center) / 0.7
z= fload_gas_xyz('z',center=center) / 0.7
hsml= fload_gas_hsml(1)
ids= fload_gas_id(1)
rho= fload_gas_rho(1)
m= fload_gas_mass(1)

if keyword_set(stars) then begin
        x= fload_allstars_xyz('x',center=center) / 0.7
        y= fload_allstars_xyz('y',center=center) / 0.7
        z= fload_allstars_xyz('z',center=center) / 0.7
        hsml= x*0.0 + 0.2
        ids= fload_allstars_ids(1)
endif


N= long(n_elements(x))

print, "-----------------------------"
print, " "
print, "particles: mean, max, min = ", mean(qtoavg), max(qtoavg), min(qtoavg)


DesNgb=96L
Hmax= 100.0

if(N gt 1) then begin

	Ncoord= max([N,input_Ngrid])

	Coord=fltarr(7,Ncoord)
	Masses= fltarr(Ncoord)
	Quant_to_Average= fltarr(Ncoord)

	Coord(0,0:N-1)= x(*)
	Coord(1,0:N-1)= y(*)
	Coord(2,0:N-1)= z(*)
	Coord(3,0:N-1)= ids(*)
	Coord(4,0:input_Ngrid-1)= input_Grid(0,*)
	Coord(5,0:input_Ngrid-1)= input_Grid(1,*)
	Coord(6,0:input_Ngrid-1)= input_Grid(2,*)
	Masses(0:N-1)= m(*)
	Quant_to_Average(0:N-1)= qtoavg(*)

        Average= fltarr(input_Ngrid)      ;  full velocity field

        S = CALL_EXTERNAL('/n/home/tcox/C-Routines_for_IDL/Cube/Average/average.so', $
                'average', $
		N, $
		input_Ngrid, $
                Coord, $
                Masses, $
		DesNgb, $
		Hmax, $
                Average, $
		Quant_to_Average)

	; 
	; ---------------------
	print, "Grid: mean, max, min = ", mean(Average), max(Average), min(Average)
	print, " "

endif else begin
        print,' Cannot find anything'
        Average= [-1]
endelse

print, " "
print, "-----------------------------"

return, Average


end





; =============================================================================
; =============================================================================
;
;
;
;
;
;       MAKING IMAGES
;
;
;
;
;
;
;
;
;
; =============================================================================
; =============================================================================


;
;  Used by plotting routine to scale the 
;  images.
; ------------------------------------------------
function ScalePic, Pic, maxvalue, minvalue, $
			var_already_inlog=var_already_inlog

    Map= Pic
    ;ma= maxvalue+maxvalue*0.1
    ;mi= minvalue-minvalue*0.1
    ma= maxvalue
    mi= minvalue

    ind=where(Map lt mi) 
    if ind(0) ne -1 then Map(ind)=mi
    ind=where(Map gt ma) 
    if ind(0) ne -1 then Map(ind)=ma


    Cols=255 ; number of colors

    if keyword_set(var_already_inlog) then begin
        Pic=(Map-mi)/(ma-mi) * (cols-3) + 2
    endif else begin
	print, "Taking log of Map"
	if mi le 0.0 then begin
	    print, "PROBLEM: mi is less than or equal to zero"
	    print, "         setting to ma/1.0e+5"
	    mi= ma/1.0d-5
	endif
        Pic=(alog10(Map)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2
    endelse

    ind=where(Pic ge 256)
    if ind(0) ne -1 then Pic(ind)=255

    ind=where(Pic le 0)
    if ind(0) ne -1 then Pic(ind)=1


    invertcolors= 1
    if invertcolors eq 1 then begin
        Pic=256-Pic                          ; invert color table
        idx= where(Pic EQ 254)               ; set background to white
    endif else begin
        idx= where(Pic EQ 2)
    endelse

    if idx(0) ne -1 then Pic(idx)= 1       ; this is setting it to white
    ;if idx(0) ne -1 then Pic(idx)= 0       ; this sets background to black

   NewScaledPic= Pic

   return, NewScaledPic

end




;------------------------------------------------------------------------


pro fload_all_colortables, junk


COMMON BLUECOLORTABLE, blue0, blue1, blue2
COMMON REDCOLORTABLE, red0, red1, red2
COMMON REGCOLORTABLE, c0, c1, c2


; get plotting stuff ready
;SET_PLOT, 'ps'
set_plot, 'z'
dxx= 120 & dyy= 40
device, set_resolution= [dxx, dyy]

LOADCT, 1
v1=[0,255]
v2=[0,255]
v3=[0,255]
tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)
TVLCT, blue0, blue1, blue2, /GET

LOADCT, 3
v1=[0,255]
v2=[0,255]
v3=[0,255]
tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)
TVLCT, red0, red1, red2, /GET

LOADCT, 4
;LOADCT, 13
v1=[0,255]
v2=[0,255]
v3=[0,255]
tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)
TVLCT, c0, c1, c2, /GET


end





; =============================================================================
;
;   Image of the grid
;
;
;  -----------------------
;  |          |          |
;  | Gas      | Gas      |
;  | Density  | Velocity |
;  |          |          |
;  |          |          |
;  -----------------------
;  |          |          |
;  | Gas      |          |
;  | Temp     |  SFR     |
;  |          |          |
;  |          |          |
;  -----------------------
;  |          |          |
;  | Cold     | Stellar  |
;  | Fraction | Mass     |
;  |          |          |
;  |          |          |
;  -----------------------
;
;



pro grid_to_image, gridfile


COMMON BLUECOLORTABLE
COMMON REDCOLORTABLE
COMMON REGCOLORTABLE


read_new_desika_grid, gridfile, Grid, Velocity, $
	Density, Temp, Density_Cold, Temp_Cold, ColdFraction, Density_Hot, Temp_Hot, SFR, StellarMass


; -------------------------------------------------------


;
Ncubed= n_elements(Grid)/3.0
N= long(Ncubed^(1/3.0))
Pic= fltarr(N,N)


zvalue= min(abs(Grid(2,*)))
;------------------------------
idx=where(Grid(2,*) eq zvalue)
Velocity=Velocity[2,*] & Velocity=Velocity(idx)
Density=Density(idx)
Temp=Temp(idx)
ColdFraction=ColdFraction(idx)
Density_Cold=Density_Cold(idx)
Temp_Cold=Temp_Cold(idx)
Density_Hot=Density_Hot(idx)
Temp_Hot=Temp_Hot(idx)
SFR=SFR(idx)
StellarMass_0=StellarMass[0,*] & StellarMass_0=StellarMass_0(idx)
StellarMass_1=StellarMass[1,*] & StellarMass_1=StellarMass_1(idx)
StellarMass_2=StellarMass[2,*] & StellarMass_2=StellarMass_2(idx)
StellarMass_3=StellarMass[3,*] & StellarMass_3=StellarMass_3(idx)
StellarMass_4=StellarMass[4,*] & StellarMass_4=StellarMass_4(idx)
;------------------------------
;zs= Grid(2,*)
;idx=where(zs gt zvalue)
;zvalue= min(zs(idx))
;idx=where(zs eq zvalue)
;Density=Density(idx)
;Temp=Temp(idx)
;StellarMass=StellarMass(idx)


;
; The following are already in LOG!!!!
;

; Density
; ------------
for i=0,N-1 do Pic(i,*)= Density(i*N:(i+1)*N-1)
maxtemp=max(Pic)
mintemp=min(Pic)
print, "Density: max/min ", maxtemp, mintemp
maxtemp=2.0
mintemp=-4.0
print, "setting pic scale to: max/min ", maxtemp, mintemp
Density= ScalePic(Pic,maxtemp,mintemp,/var_already_inlog)


; Velocity
; ------------
for i=0,N-1 do Pic(i,*)= Velocity(i*N:(i+1)*N-1)
maxtemp=max(Pic)
mintemp=min(Pic)
print, "Velocity: max/min ", maxtemp, mintemp
maxtemp=300.0
mintemp=-300.0
print, "setting pic scale to: max/min ", maxtemp, mintemp
Velocity= ScalePic(Pic,maxtemp,mintemp,/var_already_inlog)


; Temp
; ------------
for i=0,N-1 do Pic(i,*)= Temp(i*N:(i+1)*N-1)
maxtemp=max(Pic)
mintemp=min(Pic)
print, "Temp: max/min ", maxtemp, mintemp
maxtemp=9.0
mintemp=2.0
print, "setting pic scale to: max/min ", maxtemp, mintemp
Temp= ScalePic(Pic,maxtemp,mintemp,/var_already_inlog)


; SFR
; ------------
for i=0,N-1 do Pic(i,*)= SFR(i*N:(i+1)*N-1)
maxtemp=max(Pic)
mintemp=min(Pic)
print, "SFR: max/min, total ", maxtemp, mintemp, total(Pic)
maxtemp= -3.0
mintemp= -10.0
print, "setting pic scale to: max/min ", maxtemp, mintemp
SFR= ScalePic(Pic,maxtemp,mintemp)


; ColdFraction
; ------------
for i=0,N-1 do Pic(i,*)= ColdFraction(i*N:(i+1)*N-1)
maxtemp=max(Pic)
mintemp=min(Pic)
print, "ColdFraction: max/min ", maxtemp, mintemp
maxtemp=1.0
mintemp=0.8
print, "setting pic scale to: max/min ", maxtemp, mintemp
ColdFraction= ScalePic(Pic,maxtemp,mintemp,/var_already_inlog)

for i=0,N-1 do Pic(i,*)= Density_Cold(i*N:(i+1)*N-1)
print, "Density_Cold: max/min ", max(Pic), min(Pic)
for i=0,N-1 do Pic(i,*)= Temp_Cold(i*N:(i+1)*N-1)
print, "Temp_Cold: max/min ", max(Pic), min(Pic)
for i=0,N-1 do Pic(i,*)= Density_Hot(i*N:(i+1)*N-1)
print, "Density_Hot: max/min ", max(Pic), min(Pic)
for i=0,N-1 do Pic(i,*)= Temp_Hot(i*N:(i+1)*N-1)
print, "Temp_Hot: max/min ", max(Pic), min(Pic)


; StellarMass
; ------------
for i=0,N-1 do Pic(i,*)= StellarMass_0(i*N:(i+1)*N-1)
maxtemp=max(Pic)
mintemp=min(Pic)
print, "StellarMass: max/min ", maxtemp, mintemp
maxtemp=9.0
mintemp=0.0
print, "setting pic scale to: max/min ", maxtemp, mintemp
StellarMass= ScalePic(Pic,maxtemp,mintemp,/var_already_inlog)



;---------------------------------------------


;   write jpg files
; --------------------
pixels= (size(Density))[1]
zvalue= idx(0)
image3d = BYTARR(3, 2*pixels, 3*pixels)

; row 3

image3d(0, *, 0:pixels-1)= [blue0(ColdFraction), c0(StellarMass)]
image3d(1, *, 0:pixels-1)= [blue1(ColdFraction), c1(StellarMass)]
image3d(2, *, 0:pixels-1)= [blue2(ColdFraction), c2(StellarMass)]


; row 2

image3d(0, *, pixels:2*pixels-1)= [red0(Temp), c0(SFR)]
image3d(1, *, pixels:2*pixels-1)= [red1(Temp), c1(SFR)]
image3d(2, *, pixels:2*pixels-1)= [red2(Temp), c2(SFR)]

; row 1

image3d(0, *, 2*pixels:3*pixels-1)= [blue0(Density), c0(Velocity)]
image3d(1, *, 2*pixels:3*pixels-1)= [blue1(Density), c1(Velocity)]
image3d(2, *, 2*pixels:3*pixels-1)= [blue2(Density), c2(Velocity)]


;---------------------------------------------

; resize the array

image3d_big = REBIN(image3d, 3, 600, 900, /SAMPLE)


;---------------------------------------------


; assume filename is *.txt
imagefilename= strmid(gridfile,0,strlen(gridfile)-4)+'.jpg'
WRITE_JPEG, imagefilename, image3d_big, TRUE=1, quality= 95



;fout= imagedir+'/'+imagename+'_'+ilbl+'.jpg'
;timelbl=fload_timelbl(0.7,2,/noteq)

;cmd = 'convert -font arial -antialias -pointsize 30 -fill black -draw '
;cmd = cmd + "'"
;cmd = cmd + 'text 20,46 '
;cmd = cmd + """
;cmd = cmd + timelbl
;cmd = cmd + """
;cmd = cmd + "'"
;cmd = cmd + ' -quality 95 '
;cmd = cmd + filename + " " + fout

;print, cmd
;spawn, cmd


end






; =============================================================================





; =============================================================================



pro grid_movie_images, junk

    fload_all_colortables, 1

    ;------
    frun='/raid4/tcox/ds/vc3vc3e_2'
    ;------
    gridbase= 'griddata'
    griddir= frun+'/desikagrids_new/'
    ;------

    ;spawn, "/bin/ls "+griddir+gridbase+"* | wc ",result                                    
    ;ngrds=long(result[0])-1                                                  
    ;ngrds=10
    ngrds=500


    ;for i=0,ngrds do begin
    for i=11L,ngrds do begin
    ;for i=63L,ngrds do begin
   
	gridfile= griddir+gridbase+'_'+strcompress(string(i),/remove_all)+'.txt'
	openr, 1, gridfile, ERROR=err
	close,1
	if err eq 0 then grid_to_image, gridfile

    endfor

end







