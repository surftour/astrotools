; =================================================
;
;  When is gas not rotationally supported?
;
;
; ==================================================
pro time_grot, frun


if not keyword_set(frun) then begin
        print, " "
        print, " time_grot, frun"
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



; ---------------
; Arrays
; ---------------
time= fltarr(nsnaps)
rotfrac=  fltarr(nsnaps)
radfrac=  fltarr(nsnaps)
rsupfrac= fltarr(nsnaps)


; std vc3c mergers
startid1= 1L
numpart1= 200001L
startid2= 200002L
numpart2= 200001L




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


        ; instead use BH positions
        ; ------------------------
	use_bh_positions= 1
        if use_bh_positions eq 1 then begin

           ; get blackhole id's
           if i eq 0 then begin
                bhid= fload_blackhole_id(1)
                bhid1= bhid[0]
                bhid2= bhid[1]
                print, "Blackhole ID's: ", bhid1, bhid2
           endif 

           if fload_npart(5) gt 1 then begin
                center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
                center_2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
                ;center_1= fload_blackhole_xyz('xyz',idtofollow=bhid1)
                ;center_2= fload_blackhole_xyz('xyz',idtofollow=bhid2)
           endif else begin
                bhid= fload_blackhole_id(1)
                center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
                ;center_1= fload_blackhole_xyz('xyz',idtofollow=bhid)
                center_2= center_1
           endelse
        endif


        ; center-of-mass velocity, and central density velocity
        ;   both of these need the particle numbers to be set
        ;   correctly.
	; -------------------------------------------------------
        comvel_1= fload_1gal_comvel(startid1,numpart1,center=center_1,rfact=1.0)
        comvel_2= fload_1gal_comvel(startid2,numpart2,center=center_2,rfact=1.0)

        cenvel_1= fload_1gal_comvel(startid1,numpart1,center=center_1,rfact=0.05)
        cenvel_2= fload_1gal_comvel(startid2,numpart2,center=center_2,rfact=0.05)



	; gas properties
	; -----------------
	gas_r1= fload_1gal_gas_xyz('r', startid1, numpart1, center=center_1)
	gas_r2= fload_1gal_gas_xyz('r', startid2, numpart2, center=center_2)
	gas_mass1= fload_1gal_gas_mass(1, startid1, numpart1)
	gas_mass2= fload_1gal_gas_mass(1, startid2, numpart2)

	;rad= [fload_1gal_gas_xyz('r',startid1,numpart1), fload_1gal_gas_xyz('r',startid2,numpart2)]
	;x= [fload_1gal_gas_xyz('x',startid1,numpart1), fload_1gal_gas_xyz('x',startid2,numpart2)]
	;y= [fload_1gal_gas_xyz('y',startid1,numpart1), fload_1gal_gas_xyz('y',startid2,numpart2)]
	;z= [fload_1gal_gas_xyz('z',startid1,numpart1), fload_1gal_gas_xyz('z',startid2,numpart2)]
	;gas_mass= [fload_1gal_gas_mass(1,startid1,numpart1), fload_1gal_gas_mass(1,startid2,numpart2)]

	v_rad1= fload_1gal_gas_v('r', startid1, numpart1, center=center_1, com_velocity=cenvel_1)
	;v_rot1= fload_1gal_gas_v('theta', startid1, numpart1, center=center_1, com_velocity=cenvel_1)
	;v_rot1= fload_1gal_gas_v('tan', startid1, numpart1, center=center_1, com_velocity=cenvel_1)
	v_tot1= fload_1gal_gas_v('tot', startid1, numpart1, center=center_1, com_velocity=cenvel_1)
	v_rot1= sqrt(v_tot1*v_tot1 - v_rad1*v_rad1)

	v_rad2= fload_1gal_gas_v('r', startid2, numpart2, center=center_2, com_velocity=cenvel_2)
	;v_rot2= fload_1gal_gas_v('theta', startid2, numpart2, center=center_2, com_velocity=cenvel_2)
	;v_rot2= fload_1gal_gas_v('tan', startid2, numpart2, center=center_2, com_velocity=cenvel_2)
	v_tot2= fload_1gal_gas_v('tot', startid2, numpart2, center=center_2, com_velocity=cenvel_2)
	v_rot2= sqrt(v_tot2*v_tot2 - v_rad2*v_rad2)

	v_rotsup1= v_rot1
	v_rotsup2= v_rot2
	v_rotsup1(*)= 1.0e+6    ; large value guarantees notrotsupport
	v_rotsup2(*)= 1.0e+6


	; ------------------------------

	minr= 0.0
	maxr= 50.0
	dr= 0.1
	nr= long(maxr/dr)

	minr1= fltarr(nr)
	minr2= fltarr(nr)

	allr1= fload_all_xyz('r', center= center_1)
	allr2= fload_all_xyz('r', center= center_2)
	all_mass= fload_all_mass(1)

	for gi= 0L, (nr-1) do begin

	  radiusx= dr*(gi+1)
	  radiusx0= dr*gi

	  idx1= where(allr1 lt radiusx)
	  idx2= where(allr2 lt radiusx)

	  minr1= total(all_mass(idx1))
	  minr2= total(all_mass(idx2))

	  circvel1= sqrt(43007.1 * minr1 / radiusx)
	  circvel2= sqrt(43007.1 * minr2 / radiusx)

	  idx1= where((gas_r1 gt radiusx0) and (gas_r1 le radiusx))
	  idx2= where((gas_r2 gt radiusx0) and (gas_r2 le radiusx))
	  
	  if idx1(0) ne -1 then v_rotsup1(idx1)= circvel1
	  if idx2(0) ne -1 then v_rotsup2(idx2)= circvel2

	endfor



	;mass_in_r= gas_mass
	;xyz= fltarr(3)
	;nloop= n_elements(gas_mass)
	;for gi= 0L, (nloop-1) do begin
	;  xyz[0]= x[gi]
	;  xyz[1]= y[gi]
	;  xyz[2]= z[gi]
	;  all_r= fload_all_xyz('r', center= xyz)
	;  idx= where(all_r lt rad[i])
	;  mass_in_r[gi]= total(all_mass(idx))
        ;  if (gi mod 2000) eq 0 then begin
	;	print, "xxxxxx"
	;	print, "   gi= ",gi, "  xyz= ", xyz
	;	print, "    r= ", rad[gi], ' m= ', mass_in_r[gi],' v= ', sqrt(43007.1*mass_in_r[gi]/rad[gi])
	;  endif
	;endfor

	;rotsupported_v= sqrt(43007.1 * mass_in_r / rad)

	; ------------------------------


	v_rad= [v_rad1, v_rad2]
	v_rot= [v_rot1, v_rot2]
	v_tot= [v_tot1, v_tot2]
	radf= v_rad/v_tot
	rotf= v_rot/v_tot
	rotfrac[i]=  mean(rotf)
	radfrac[i]=  mean(radf)

	; ------------------------------

	v_rotsup= [v_rotsup1, v_rotsup2]
	gas_mass= [gas_mass1, gas_mass2]
	idx= where(v_rot gt 0.8*v_rotsup)
	if idx(0) ne -1 then rsuppmass= gas_mass(idx) else rsuppmass= 0.0
	rsupfrac[i]= total(rsuppmass)/total(gas_mass)


endfor



; ----------------------------------------
; write to file

openw, 1, frun+'/rotrad.txt', ERROR=err
        
printf, 1, "#   rotrad.txt"     
printf, 1, "#             "
printf, 1, "#               "
printf, 1, "#    time     v_rot/      v_rad/   rot support"
printf, 1, "#  (Gyr/h)     v_tot       v_tot       M_gas"
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F9.4,"    ",F9.5,"    ",F9.5,"    ",F9.5)', $
                time[i], rotfrac[i], radfrac[i], rsupfrac[i]
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
;  Read rotrad.txt file
; ----------------------------
pro read_rotrad_file, frun, time, rotfrac, radfrac, rotsuppfrac

filename= frun+'/rotrad.txt'

spawn, "wc "+filename,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then filedata= fltarr(4,lines)

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
rotfrac= filedata[1,*]
radfrac= filedata[2,*]
rotsuppfrac= filedata[3,*]


end





;===================================================================================





;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro plot_rotrad, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_rotrad, junk"
	print, " "
	print, " "
	return
endif

filename='rotrad.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



yaxistitle = "!6Rotationally Supported Fraction"
xaxistitle = "!6Time (Gyr)"
xmax = 4.25
xmin = 0
ymax = 0.95
;ymin = 0.01
ymin = 0.0


;---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]

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

;----------


;process_one_rotrad, "/raid4/tcox/vc3vc3e_2", lcolor=150, lpsym= 5, msg='No Wind', x0= 0.54, y0= 0.84
process_one_rotrad, "/raid4/tcox/vc3vc3e_no", lcolor=150, lpsym= 5, msg='No Wind', x0= 0.64, y0= 0.87

; eta=0.5 (v_w= 837, 105)
; ------------------------
;process_one_rotrad, "/raid4/tcox/sbw/sb10", lcolor=100, msg='sb10', x0= 0.54, y0= 0.26
;process_one_rotrad, "/raid4/tcox/sbw/sb17", lcolor=50, msg='sb17', x0= 0.54, y0= 0.22


; eta=2.0 (v_w= 837, 105)
; ------------------------
;process_one_rotrad, "/raid4/tcox/sbw/sb8", lcolor=100, msg='sb8', x0= 0.54, y0= 0.26
;process_one_rotrad, "/raid4/tcox/sbw/sb14", lcolor=50, msg='sb14', x0= 0.54, y0= 0.22



process_one_rotrad, "/raid4/tcox/sbw/sb17", lcolor=200, msg='!7g!6= 0.05', x0= 0.64, y0= 0.82
process_one_rotrad, "/raid4/tcox/sbw/sb14", lcolor=100, msg='!7g!6= 2.0', x0= 0.64, y0= 0.77
process_one_rotrad, "/raid4/tcox/sbw/sb22", lcolor=50, msg='!7g!6= 5.0', x0= 0.64, y0= 0.72


device, /close




end






pro process_one_rotrad, frun, lcolor=lcolor, lpsym=lpsym, $
				msg=msg, x0=x0, y0=y0

	if not keyword_set(lcolor) then lcolor= 0
	if not keyword_set(lpsym) then lpsym= 3


	read_rotrad_file, frun, time, rotfrac, radfrac, rotsuppfrac
	time=time/0.7

	;oplot, time, radfrac, thick=5.0, psym=-lpsym, color= lcolor, linestyle=1
	;oplot, time, rotfrac, thick=5.0, psym=-lpsym, color= lcolor
	oplot, time, rotsuppfrac, thick=5.0, psym=-lpsym, color= lcolor

	if not keyword_set(x0) then x0= 0.75
	if not keyword_set(y0) then y0= 0.85

	if keyword_set(msg) then begin
		xyouts, x0, y0, msg, /normal, charthick=3, size=1.5, color=lcolor
	endif else begin
		xyouts, x0, y0, fload_getid(frun), /normal, charthick=1, size=1.33, color=150
	endelse



end








;==========================================================================





