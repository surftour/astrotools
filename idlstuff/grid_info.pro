;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   Calculated grid properties
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------




;---------------------
;
;   Read Grid Basics
;
;---------------------
;

pro get_grid_info, Grid_Volume, Grid_rthetaphi=Grid_rthetaphi, $
				Grid_dr=Grid_dr


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


	openr, 1, "/home/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_r.dat"
	readf, 1, grid_r
	close, 1
	openr, 1, "/home/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_r.dat"
	readf, 1, grid_cen_r
	close, 1

	openr, 1, "/home/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_theta.dat"
	readf, 1, grid_theta
	close, 1
	openr, 1, "/home/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_theta.dat"
	readf, 1, grid_cen_theta
	close, 1

	openr, 1, "/home/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/grid_phi.dat"
	readf, 1, grid_phi
	close, 1
	openr, 1, "/home/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/cent_phi.dat"
	readf, 1, grid_cen_phi
	close, 1


	r_len= max(grid_r)
	total_Volume= (4./3.)*!PI*r_len*r_len*r_len
	print, "total_Volume= ", total_Volume

	Ngrid= N_r*N_theta*N_phi

	Grid= fltarr(3,Ngrid)
	Grid_rthetaphi= fltarr(3,Ngrid)
	Grid_Volume= fltarr(Ngrid)
	Grid_dr= fltarr(Ngrid)

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
		Grid_dr[n_i]= delta_r

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




end



; =============================================================================




;  spherical
; ------------

pro get_grid_info2, Ngrid, Grid, Grid_Volume, Grid_MaxRadius
;pro get_grid_info, grid_r, grid_theta, grid_phi, $
;               grid_cen_r, grid_cen_theta, grid_cen_phi, $


        ; need the following order

        ; rarr, phiarr, tharr
        ; --------------------

        ; old, r^-2 grids
        ;N_r= 400L     &  grid_r= fltarr(N_r)
        ;N_theta= 119L  &  grid_theta= fltarr(N_theta)
        ;N_phi= 238L   &  grid_phi= fltarr(N_phi)

        ; log r grids
        N_r= 600L     &  grid_r= fltarr(N_r)
        N_theta= 131L  &  grid_theta= fltarr(N_theta)
        N_phi= 262L   &  grid_phi= fltarr(N_phi)

        ; centers
        ; ---------
        grid_cen_r= fltarr(N_r-1)
        grid_cen_theta= fltarr(N_theta-1)
        grid_cen_phi= fltarr(N_phi-1)


        ; radius
        ; --------
        openr, 1, "/home/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/rarr.dat"
        junk=''
        readf, 1, junk
        print, junk
        readf, 1, junk
        r_data=fltarr(3,N_r)
        readf, 1, r_data
        close, 1
        grid_r= r_data[1,*]

        for i=0,N_r-2 do begin
                grid_cen_r[i]= 0.5*(grid_r[i+1]+grid_r[i])
        endfor

        ; theta
        ; --------
        openr, 1, "/home/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/tharr.dat"
        junk=''
        readf, 1, junk
        print, junk
        readf, 1, junk
        th_data=fltarr(5,N_theta)
        readf, 1, th_data
        grid_theta= th_data[1,*]
        close, 1

        for i=0,N_theta-2 do begin
                grid_cen_theta[i]= 0.5*(grid_theta[i+1]+grid_theta[i])
        endfor


        ; phi
        ; --------
        openr, 1, "/home/tcox/Tools/C-Routines_for_IDL/Cube/SphericalGridPoint/phiarr.dat"
        junk=1.0
        readf, 1, junk
        readf, 1, grid_phi
        close, 1

        ; it's in degrees
        maxphi= max(grid_phi)
        grid_phi= 2*!PI * grid_phi / maxphi

        for i=0,N_phi-2 do begin
                grid_cen_phi[i]= 0.5*(grid_phi[i+1]+grid_phi[i])
        endfor

        r_len= max(grid_r)

        grid_radius= 5.0
        if keyword_set(grid_radius) then begin
                print, " ******  MANUALLY setting grid max radius to ", grid_radius, "  ******** "
                grid_r= grid_r * ( grid_radius / r_len)
                grid_cen_r= grid_cen_r * ( grid_radius / r_len)
                r_len= max(grid_r)
        endif

        Grid_MaxRadius= r_len

        print, "max(grid_r)= ", max(grid_r)
        print, "max(grid_cen_r)= ", max(grid_cen_r)
        print, "r_len= ", r_len
        total_Volume= (4./3.)*!PI*r_len*r_len*r_len
        print, "total_Volume= ", total_Volume

        Ngrid= (N_r-1)*(N_theta-1)*(N_phi-1)

        Grid= fltarr(3,Ngrid)
        Grid_rthetaphi= fltarr(3,Ngrid)
        Grid_Volume= fltarr(Ngrid)



        n_i= long(0)
        for i=0,N_r-2 do begin
          for j=0,N_theta-2 do begin
            for k=0,N_phi-2 do begin

                if (n_i mod 500000) eq 0 then print, "n_i= ", n_i

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


print, "Grid_Volume max/min= ", max(Grid_Volume), min(Grid_Volume)
print, "Volume check= ", total(Grid_Volume)


;  Volume check
; --------------
sphere_area= 0.0
for j=0,N_theta-2 do begin
  for k=0,N_phi-2 do begin
        delta_theta= grid_theta[j+1] - grid_theta[j]
        delta_phi= grid_phi[k+1] - grid_phi[k]
        sphere_area= sphere_area + delta_phi * sin(grid_theta[j]) * delta_theta
  endfor
endfor
print, "sphere area= ", sphere_area
print, '4*PI= ', 4.0*!PI



end



; =============================================================================







; =============================================================================




pro read_sukanya_grid_full, Grid, Velocity, Density, Density_Cold, Temp, ColdFF, StellarMass


gridfile= '/home/tcox/grid.txt'

spawn, "wc "+gridfile, result
lines=long(result)
datalines=lines(0)-5

Grid= fltarr(3,lines)
SphericalGrid= fltarr(3,lines)
Density= fltarr(lines)
Density_Cold= fltarr(lines)
Temp= fltarr(lines)
ColdFF= fltarr(lines)
StellarMass= fltarr(lines)

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
        SphericalGrid(0,i)= float(tempjunk(3))
        SphericalGrid(1,i)= float(tempjunk(4))
        SphericalGrid(2,i)= float(tempjunk(5))
        Density(i)= float(tempjunk(6))
        Density_Cold(i)= float(tempjunk(7))
        Temp(i)= float(tempjunk(8))
        ColdFF(i)= float(tempjunk(9))
        StellarMass(i)= float(tempjunk(10))
endfor

close, 1


end





; =============================================================================





pro read_sukanya_grid_coldgasonly, Density_Cold, gridfile=gridfile


if not keyword_set(gridfile) then gridfile= '/home/tcox/grid_cold.txt'

spawn, "wc "+gridfile, result
lines=long(result)
;datalines=lines(0)-5
datalines=lines(0)

Density_Cold= fltarr(datalines)

openr, 1, gridfile
;junk=''

readf, 1, Density_Cold

; read data
;for i=0L,lines(0)-1 do begin
;        readf, 1, junk
;        ;print, junk
;        tempjunk= strsplit(junk,/extract,count=count)
;	if count ne 1 then print, "what's up with line: ", i
;        Density_Cold(i)= float(tempjunk(0))
;endfor

close, 1


end





; =============================================================================




pro read_desika_grid, gridfile, Grid, Velocity, Density, Temp, StellarMass


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

Grid= fltarr(3,lines)
Velocity= fltarr(3,lines)
Density= fltarr(lines)
Temp= fltarr(lines)
StellarMass= fltarr(4,lines)

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
        StellarMass(0,i)= float(tempjunk(8))
        StellarMass(1,i)= float(tempjunk(9))
        StellarMass(2,i)= float(tempjunk(10))
        StellarMass(3,i)= float(tempjunk(11))
endfor

close, 1


end






; =============================================================================







;
;  Used by plotting routine to scale the 
;  images.
; ------------------------------------------------
function ScalePic, Pic, maxvalue, minvalue, $
			var_already_inlog=var_already_inlog

    Map= Pic
    ma= maxvalue+maxvalue*0.1
    mi= minvalue-minvalue*0.1

    ind=where(Map lt mi) 
    if ind(0) ne -1 then Map(ind)=mi
    ind=where(Map gt ma) 
    if ind(0) ne -1 then Map(ind)=ma


    Cols=255 ; number of colors

    if keyword_set(var_already_inlog) then begin
        Pic=(Map-mi)/(ma-mi) * (cols-3) + 2
    endif else begin
	print, "Taking log of Map"
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








; =============================================================================


;   Image the grid we just created


;  -------------------
;  |        |        |
;  | Vel(1) |Density |
;  |        |        |
;  |        |        |
;  --------------------
;  |        |        |
;  | Gas    |Stellar |
;  | Temp   |  Mass  |
;  |        |        |
;  --------------------



pro look_at_grid, junk

;frun= '/raid4/tcox/As/A5e'
;------
;gridf= 'griddata_12.txt'
;gridfile= frun+'/desikagrids/'+gridf
;filename='gridpic.eps'
;------
;gridf= 'griddata_12.txt'
;gridfile= frun+'/desikagrids_bad/'+gridf
;filename='gridpicbad.eps'
;------
;frun= 'yuexing'
;gridf= 'griddata_8.txt'
;gridfile= frun+'/50_8/'+gridf
;filename='yuexing/50_8/griddata_8.eps'
;gridfile= frun+'/desikagrids/'+gridf
;filename='yuexing/desikagrids/griddata_8.eps'
;gridfile= frun+'/desikagrids_orig/'+gridf
;filename='yuexing/desikagrids_orig/griddata_8.eps'
;------
frun='/raid4/tcox/vc3vc3e_2'
;gridf= 'griddata_42.txt'
;gridfile= frun+'/desikagrids_v1/'+gridf
;filename=frun+'/desikagrids_v1/griddata_42.eps'
;gridf= 'griddata_42.txt.12.orig'
;gridfile= frun+'/desikagrids/'+gridf
;filename=frun+'/desikagrids/griddata_42orig.eps'
;gridf= 'griddata_42.txt.8.orig'
;gridfile= frun+'/desikagrids/'+gridf
;filename=frun+'/desikagrids/griddata_42orig2.eps'
gridf= 'griddata_42.txt'
gridfile= frun+'/desikagrids/'+gridf
filename=frun+'/desikagrids/griddata_42.eps'
;------

read_desika_grid, gridfile, Grid, Velocity, Density, Temp, StellarMass
;read_sukanya_file, Grid, Velocity, Density, Density_Cold, Temp, ColdFF, StellarMass


; -------------------------------------------------------


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newysize=16, newxsize=16

;
Ncubed= n_elements(Grid)/3.0
N= long(Ncubed^(1/3.0))
Pic= fltarr(N,N)

;zvalue= 2.04
zvalue= -0.12
;zvalue= -0.08
;------------------------------
idx=where(Grid(2,*) eq zvalue)
Density=Density(idx)
Temp=Temp(idx)
StellarMass=StellarMass(idx)
;------------------------------
;stop
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
Pic= ScalePic(Pic,maxtemp,mintemp,/var_already_inlog)
tv, Pic, 0.5, 0.5, xsize=0.5, ysize=0.5, /normal
xyouts, 0.53, 0.95, "Gas Density", /normal, size=1.2, color= 0, charthick=4.0


; Temp
; ------------
for i=0,N-1 do Pic(i,*)= Temp(i*N:(i+1)*N-1)
maxtemp=max(Pic)
mintemp=min(Pic)
print, "Temp: max/min ", maxtemp, mintemp
maxtemp=6.5
mintemp=0.0
print, "setting pic scale to: max/min ", maxtemp, mintemp
Pic= ScalePic(Pic,maxtemp,mintemp,/var_already_inlog)
tv, Pic, 0.0, 0.0, xsize=0.5, ysize=0.5, /normal
xyouts, 0.03, 0.45, "Gas Temp", /normal, size=1.2, color= 0, charthick=4.0



; StellarMass
; ------------
for i=0,N-1 do Pic(i,*)= StellarMass(i*N:(i+1)*N-1)
maxtemp=max(Pic)
mintemp=min(Pic)
print, "StellarMass: max/min ", maxtemp, mintemp
Pic= ScalePic(Pic,maxtemp,mintemp,/var_already_inlog)
tv, Pic, 0.5, 0.0, xsize=0.5, ysize=0.5, /normal
xyouts, 0.53, 0.45, "Stellar Mass", /normal, size=1.2, color= 0, charthick=4.0


xyouts, 0.13, 0.85, frun, /normal, size=1.0, color= 0, charthick=2.0
xyouts, 0.13, 0.81, gridf, /normal, size=1.0, color= 0, charthick=2.0
xyouts, 0.13, 0.70, "Slice through z="+strcompress(string(zvalue),/remove_all), /normal, size=1.2, color= 0, charthick=4.0


device, /close



end





; =============================================================================


;   Some Random Grid Info

pro grid_info1, gridfile

if not keyword_set(gridfile) then begin
        print, "  "
        print, " grid_info, gridfile"
        print, "  "
        return
endif


;frun='/raid4/tcox/vc3vc3e_2'
;gridf= 'griddata_42.txt'
;gridfile= frun+'/desikagrids_v1/'+gridf 
;gridf= 'griddata_42.txt.12.orig'
;gridfile= frun+'/desikagrids/'+gridf
;gridf= 'griddata_42.txt.8.orig'
;gridfile= frun+'/desikagrids/'+gridf
;gridf= 'griddata_42.txt'
;gridfile= frun+'/desikagrids/'+gridf
;------

read_desika_grid, gridfile, Grid, Velocity, Density, Temp, StellarMass
;read_sukanya_file, Grid, Velocity, Density, Density_Cold, Temp, ColdFF, StellarMass


print, "Density:     max/min ", max(Density), min(Density)
print, "Velocity:    max/min ", max(Velocity), min(Velocity)
print, "Temp:        max/min ", max(Temp), min(Temp)
print, "StellarMass: max/min ", max(StellarMass), min(StellarMass)


;-------------------------------------------
;
;   Define Grid
;
n_side= 50L
;n_side= 75L
gridlbl= strcompress(string(n_side),/remove_all)
griddim= gridlbl+'x'+gridlbl+'x'+gridlbl
side_len= 12.0
;side_len= 8.0
sidellbl= strcompress(string(side_len),/remove_all)
sidellbl= ' box len='+sidellbl
element_len= side_len/n_side
cellVolume= element_len*element_len*element_len

print, "Total grid volume= ", side_len*side_len*side_len
print, "cellVolume= ", cellVolume

;-------------------------------------------

Density= 10^(Density) / 674.59

print, " "
print, "Cube Gas Mass= ", total(cellVolume*Density)

StellarMass= 10^(StellarMass - 10.0)

print, " "
print, "Cube Stellar Mass= ", total(StellarMass)




end







; =============================================================================


;   Some Random Grid Info

pro grid_info, gridfile

if not keyword_set(gridfile) then begin
	print, "  "
	print, " grid_info, gridfile"
	print, "  "
	return
endif


; full grid test
do_full_gridtest= 0
if do_full_gridtest eq 1 then begin
	read_sukanya_grid_full, Grid, Velocity, Density, Density_Cold, Temp, ColdFF, StellarMass

	print, "Value            Max         Min "
	print, "-------          ------      -----"
	print, "Density      ", max(Density), min(Density)
	print, "Density_Cold ", max(Density_Cold), min(Density_Cold)
	print, "Temp", max(Temp), min(Temp)
	print, "ColdFF", max(ColdFF), min(ColdFF)
	print, "StellarMass", max(StellarMass), min(StellarMass)
endif


; just cold test
do_cold_test= 1
if do_cold_test eq 1 then begin
        read_coldgas_grid, Density_Cold

        print, " "
        print, " "
        print, " Reading: /home/tcox/grid_cold.txt"
        print, " "
        print, " Ngrid= ", n_elements(Density_Cold)
        print, " "
        print, " Value            Max         Min "
        print, " -------          ------      -----"
        print, " Density_Cold ", max(Density_Cold), min(Density_Cold)
        print, " "
        print, " "
endif


;----------------------------------------------
;
;


;get_oldgrid_info, Ngrid, Grid, Grid_Volume, Grid_MaxRadius
get_grid_info, Ngrid, Grid, Grid_Volume, Grid_MaxRadius

print, " "
print, " "
print, " Ngrid = ", Ngrid
print, " Grid Volume= ", total(Grid_Volume)
print, " Grid_MaxRadius= ", Grid_MaxRadius
print, " "
print, " "


;----------------------------------------------
;
;



;get_grid_info, Grid_Volume, $
;		Grid_rthetaphi=Grid_rthetaphi, $
;		Grid_dr=Grid_dr



;-------------------------------------------
;
;
;   Does the mass work out?
;
;

print, "Grid Volume= ", total(Grid_Volume)

if do_full_gridtest eq 1 then begin
	TotalGasMass= Grid_Volume * 10^(Density) / 404.75368
	print, "TotalGasMass= ", total(TotalGasMass)
endif


ColdGasMass= Grid_Volume * 10^(Density_Cold) / 404.75368
print, "ColdGasMass= ", total(ColdGasMass)




end






; =============================================================================






; ========================================
;
;      Histogram of LOS column densities
;
; ========================================



pro los_cd_hist, junk


filename='columnhist.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

xaxistitle= '!6Log Column Density, N!DH!N (cm!E-2!N)'
xticklbls=['21','22','23','24','25','26']
xticknum= n_elements(xticklbls)-1
;xticknum= 10
xmax= 26.0
xmin= 21.0

ymax = 1.2
ymin = 0.0

x0= 0.02
y0= 0.15
x_size= 0.96
y_size= 0.83

!p.position=[x0, y0, x0+x_size,y0+y_size]
!p.ticklen= 0.02

plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, /nodata, $
        color= 0, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, charthick=3.0, $
        xticks= xticknum, yticks= 1, ytickformat='(a1)', xtitle=xaxistitle, ytitle=yaxistitle, $
        xtickname=xticklbls



; --------------------

; do the heavy lifting

;process_los_hist, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_52_cold.txt", xmin=xmin, xmax=xmax

;process_los_hist, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_80_cold.txt", xmin=xmin, xmax=xmax, $
										hcolor= 50


; --------------------

;process_los_hist, "/home/tcox/griddata_52_total_orig.txt", xmin=xmin, xmax=xmax, hcolor= 50, /fill
;process_los_hist, "/home/tcox/griddata_52_total.txt", xmin=xmin, xmax=xmax, hcolor= 50

;process_los_hist, "/home/tcox/griddata_52_cold_15_ngb2.txt", xmin=xmin, xmax=xmax, /fill
;process_los_hist, "/home/tcox/griddata_52_cold_075.txt", xmin=xmin, xmax=xmax, hcolor= 200
;process_los_hist, "/home/tcox/griddata_52_cold_15.txt", xmin=xmin, xmax=xmax
;process_los_hist, "/home/tcox/griddata_52_cold_30.txt", xmin=xmin, xmax=xmax, hcolor= 100

; --------------------

;process_los_hist, "/home/tcox/griddata_52_total.txt", xmin=xmin, xmax=xmax, hcolor= 50
;process_los_hist, "/home/tcox/griddata_52_cold.txt", xmin=xmin, xmax=xmax
;process_los_hist, "/home/tcox/griddata_52_mass.txt", xmin=xmin, xmax=xmax, hcolor= 100

; --------------------

;process_los_hist, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_44_new.txt", xmin=xmin, xmax=xmax, hcolor= 20
;process_los_hist, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_52_new.txt", xmin=xmin, xmax=xmax, hcolor= 40
;process_los_hist, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_54_new.txt", xmin=xmin, xmax=xmax, hcolor= 80
;process_los_hist, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_64_new.txt", xmin=xmin, xmax=xmax, hcolor= 120
;process_los_hist, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_72_new.txt", xmin=xmin, xmax=xmax, hcolor= 160
;process_los_hist, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_84_new.txt", xmin=xmin, xmax=xmax, hcolor= 200

; --------------------

;process_los_hist, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_52_new.txt", xmin=xmin, xmax=xmax, hcolor= 50, normalization=450.0
;process_los_hist, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_52_cold.txt", xmin=xmin, xmax=xmax, hcolor= 150, normalization=450.0
;xyouts, 0.18, 0.92, 'New Grid', /normal, size=1.2, charthick=3.0, color=50
;xyouts, 0.18, 0.87, 'Old Grid', /normal, size=1.2, charthick=3.0, color=150

; --------------------

process_los_hist, "coldgrid_52.txt", xmin=xmin, xmax=xmax, hcolor= 50, normalization=450.0
process_los_hist, "totalgrid1_52.txt", xmin=xmin, xmax=xmax, hcolor= 100, normalization=450.0
process_los_hist, "totalgrid2_52.txt", xmin=xmin, xmax=xmax, hcolor= 150, normalization=450.0
xyouts, 0.10, 0.92, 'Cold Gas', /normal, size=1.2, charthick=3.0, color=50
xyouts, 0.10, 0.87, 'Total Gas (#1)', /normal, size=1.2, charthick=3.0, color=100
xyouts, 0.10, 0.82, 'Total Gas (#2)', /normal, size=1.2, charthick=3.0, color=150

; --------------------

device, /close


end



; =============================================================================




pro process_los_hist, frun, xmin=xmin, xmax=xmax, $
				normalization=normalization, $
				hcolor=hcolor, $
				fill=fill


cm_per_kpc= 3.085678d+21


get_grid_info, Grid_Volume, $
		Grid_rthetaphi=Grid_rthetaphi, $
		Grid_dr=Grid_dr



n_angles= long(117 * 58)    ; 6786

columndens= fltarr(n_angles)

read_sukanya_grid_coldgasonly, Density_Cold, gridfile=frun

print, "Value            Max         Min "
print, "-------          ------      -----"
print, "Density     ", max(Density_Cold), min(Density_Cold)


idx= n_angles*indgen(199)
print, "Column length= ", total(Grid_dr(idx))

idx= where(abs(Density_Cold) lt 1.0e-6)
print, "Number with zero Density= ", n_elements(idx), "  (", n_elements(idx)*100./n_angles," %)"
if idx(0) ne -1 then Density_Cold(idx)= -15.0

; determine histogram
for n_i= 0, n_angles-1 do begin

	idx= n_angles*indgen(199) + n_i

	columndens[n_i]= cm_per_kpc * total(Grid_dr(idx) * 10^(Density_Cold(idx)))

endfor

idx= where(columndens gt 0.0)
if idx(0) ne -1 then columndens(idx)= alog10(columndens(idx))

process_and_print_hist, columndens, xmin=xmin, xmax=xmax, $
				normalization=normalization, $
				hcolor=hcolor, fill=fill


end



; =============================================================================



;
; Compute Histogram
;------------------------------
pro process_and_print_hist, allvars, xmin=xmin, xmax=xmax, $
			normalization=normalization, $
			fill=fill, $
			hcolor=hcolor

if not keyword_set(hcolor) then hcolor= 150


levels= 50.0

step= (xmax-xmin)/(levels)
bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
bins=bins+(step*0.5)


; std histogram
hist_allvars= histogram(allvars, binsize=step, max=xmax, min=xmin)


; make sure it extends to 0
bins= [0,bins,xmax]
hist_allvars= [hist_allvars(0),hist_allvars,hist_allvars(levels-1)]


if not keyword_set(normalization) then normalization= float(max(hist_allvars))
print, "histogram normalization= ", normalization
hist_allvars= hist_allvars/normalization


; print it
; ----------
oplot, bins, hist_allvars, psym=10, color=hcolor, thick=6.0


; fill in histogram
; ------------------
if keyword_set(fill) then begin
	nbins= bins+(step*0.5)          ; make x coord
	;nbins[0]= 0.0
	nbins[0]= xmin
	nbins=[nbins,nbins]
	nbins= nbins(sort(nbins))
	;nbins= nbins[nbins, xmax]

	filledhist= fltarr(2.*levels + 2)
	for i=1,levels do begin
	   filledhist[2.*i-1]= hist_allvars[i]
	   filledhist[2.*i]= hist_allvars[i]
	endfor


	;                               thick=3.0, orientation=90.0
	polyfill, nbins, filledhist, /data, color= hcolor, /line_fill, linestyle=0, $
                                thick=3.0, orientation=45.0
endif


end




; =============================================================================




; ========================================
;
;      Quantify the Clumpiness of 
;            the grid density
;
; ========================================






pro grid_clumpiness, junk




filename='clumpiness.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

xaxistitle= '!6Log Density, n!DH!N (cm!E-3!N)'
;xticklbls=['21','22','23','24','25','26']
;xticknum= n_elements(xticklbls)-1
;xticknum= 10
xmax= 6.0
xmin= -5.0

yaxistitle= '!6Volume Fraction'
ymax = 1.0
ymin = 1.0e-8

x0= 0.17
y0= 0.15
x_size= 0.81
y_size= 0.83

!p.position=[x0, y0, x0+x_size,y0+y_size]
!p.ticklen= 0.02

plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, /nodata, $
        color= 0, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, charthick=3.0, /ylog, $
        ;xticks= xticknum, yticks= 1, $
	ytickformat='exp_label', $
	xtitle=xaxistitle, ytitle=yaxistitle, $
        xtickname=xticklbls



; --------------------

; do the heavy lifting

;process_clumpiness, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_52_cold.txt", xmin=xmin, xmax=xmax

;process_clumpiness, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_80_cold.txt", xmin=xmin, xmax=xmax, $
;                                                                                hcolor= 50

; --------------------

;process_clumpiness, "/home/tcox/griddata_52_cold_075.txt", xmin=xmin, xmax=xmax, hcolor= 200
;process_clumpiness, "/home/tcox/griddata_52_cold_15.txt", xmin=xmin, xmax=xmax
;process_clumpiness, "/home/tcox/griddata_52_cold_15_ngb2.txt", xmin=xmin, xmax=xmax, hpsym=2
;process_clumpiness, "/home/tcox/griddata_52_cold_30.txt", xmin=xmin, xmax=xmax

;process_clumpiness, "/home/tcox/griddata_52_total_orig.txt", xmin=xmin, xmax=xmax, hcolor= 50
;process_clumpiness, "/home/tcox/griddata_52_total.txt", xmin=xmin, xmax=xmax, hcolor= 50, hpsym=5

; --------------------

;process_clumpiness, "/home/tcox/griddata_52_mass.txt", xmin=xmin, xmax=xmax, hcolor= 100
;process_clumpiness, "/home/tcox/griddata_52_cold.txt", xmin=xmin, xmax=xmax
;process_clumpiness, "/home/tcox/griddata_52_total.txt", xmin=xmin, xmax=xmax, hcolor= 50

; --------------------

;process_clumpiness, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_44_new.txt", xmin=xmin, xmax=xmax, hcolor= 20
;process_clumpiness, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_52_new.txt", xmin=xmin, xmax=xmax, hcolor= 40
;process_clumpiness, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_54_new.txt", xmin=xmin, xmax=xmax, hcolor= 80
;process_clumpiness, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_64_new.txt", xmin=xmin, xmax=xmax, hcolor= 120
;process_clumpiness, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_72_new.txt", xmin=xmin, xmax=xmax, hcolor= 160
;process_clumpiness, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_84_new.txt", xmin=xmin, xmax=xmax, hcolor= 200

; --------------------

;process_clumpiness, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_52_new.txt", xmin=xmin, xmax=xmax, hcolor= 50
;process_clumpiness, "/raid4/tcox/vc3vc3e_2/sukanya_grids/griddata_52_cold.txt", xmin=xmin, xmax=xmax, hcolor= 150
;xyouts, 0.18, 0.92, 'New Grid', /normal, size=1.2, charthick=3.0, color=50
;xyouts, 0.18, 0.87, 'Old Grid', /normal, size=1.2, charthick=3.0, color=150

; --------------------

process_clumpiness, "coldgrid_52.txt", xmin=xmin, xmax=xmax, hcolor= 50
process_clumpiness, "totalgrid1_52.txt", xmin=xmin, xmax=xmax, hcolor= 100
process_clumpiness, "totalgrid2_52.txt", xmin=xmin, xmax=xmax, hcolor= 150
xyouts, 0.68, 0.92, 'Cold Gas', /normal, size=1.2, charthick=3.0, color=50
xyouts, 0.68, 0.87, 'Total Gas (#1)', /normal, size=1.2, charthick=3.0, color=100
xyouts, 0.68, 0.82, 'Total Gas (#2)', /normal, size=1.2, charthick=3.0, color=150

; --------------------

device, /close


end






; =============================================================================




pro process_clumpiness, frun, xmin=xmin, xmax=xmax, $
                                hcolor=hcolor, $
				hpsym=hpsym



if not keyword_set(hcolor) then hcolor= 150
if not keyword_set(hpsym) then hpsym= 3

get_grid_info, Grid_Volume, $
                Grid_rthetaphi=Grid_rthetaphi, $
                Grid_dr=Grid_dr



nbins= 50

; we assume density is in log here
d_log_den= (xmax-xmin)/nbins


volumef= fltarr(nbins)
dens= fltarr(nbins)

read_sukanya_grid_coldgasonly, Density_Cold, gridfile=frun

print, "Value            Max         Min "
print, "-------          ------      -----"
print, "Density_Cold ", max(Density_Cold), min(Density_Cold)


idx= where(abs(Density_Cold) lt 1.0e-6)
if idx(0) ne -1 then Density_Cold(idx)= -15.0
;print, "Volume Fraction with Cold Gas= ", (total(Grid_Volume)-total(Grid_Volume(idx)))/total(Grid_Volume)


; determine histogram
for n_i= 0, nbins-1 do begin

	dens_min= xmin + n_i * d_log_den
	dens[n_i]= dens_min

        idx= where(Density_Cold gt dens_min)

	if idx(0) ne -1 then begin
		volumef[n_i]= total(Grid_Volume(idx))/total(Grid_Volume)
	endif else begin
		volumef[n_i]= 0
	endelse

endfor


; print it
; ----------
oplot, dens, volumef, psym=-hpsym, color=hcolor, thick=4.0




end



; =============================================================================




pro grid_slice, junk



;
;  Get Grid Info
; ----------------
get_grid_info, Grid_Volume, $
                Grid_rthetaphi=Grid_rthetaphi, $
                Grid_dr=Grid_dr


read_sukanya_grid_coldgasonly, Density_Cold, gridfile=frun


idx= where(abs(Density_Cold) lt 1.0e-6)
if idx(0) ne -1 then Density_Cold(idx)= -15.0
;print, "Volume Fraction with Cold Gas= ", (total(Grid_Volume)-total(Grid_Volume(idx)))/total(Grid_Volume)





;
;  Make Plot
; ------------


filename='gridslice.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

;x0= 0.17
;y0= 0.15
;x_size= 0.81
;y_size= 0.83

;!p.position=[x0, y0, x0+x_size,y0+y_size]
;!p.ticklen= 0.02

;plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, /nodata, $
;        color= 0, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, charthick=3.0, /ylog, $
;        ;xticks= xticknum, yticks= 1, $
;        ytickformat='exp_label', $
;        xtitle=xaxistitle, ytitle=yaxistitle, $
;        xtickname=xticklbls


;oplot, dens, volumef, psym=-hpsym, color=hcolor, thick=4.0


   ; Establish plot coordinates.

Plot, radius, angle, /Polar, XStyle=5, YStyle=5, $
   /NoData, Background=1, Position=Aspect(1.0)
   
   ; Draw axis through center.
   
Axis, /XAxis, 0, 0, Color=0
Axis, /YAxis, 0, 0, Color=0

   ; Plot data.
   
OPlot, radius, angle, PSym=2, /Polar, Color=3


end



