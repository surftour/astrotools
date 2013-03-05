pro time_mags, frun



;------------------------------------
;
;  MAKE SURE TO
;
;  .run determine_lums
;
;------------------------------------


if not keyword_set(frun) then begin
	print, " "
	print, " time_mags, frun"
	print, " "
	print, " needs: .run determine_lums"
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



time= fltarr(nsnaps)

bolo=     fltarr(nsnaps)
U_mag=    fltarr(nsnaps)
B_mag=    fltarr(nsnaps)
V_mag=    fltarr(nsnaps)
R_mag=    fltarr(nsnaps)
I_mag=    fltarr(nsnaps)
J_mag=    fltarr(nsnaps)
H_mag=    fltarr(nsnaps)
K_mag=    fltarr(nsnaps)
u_sdss_mag=    fltarr(nsnaps)
g_sdss_mag=    fltarr(nsnaps)
r_sdss_mag=    fltarr(nsnaps)
i_sdss_mag=    fltarr(nsnaps)
z_sdss_mag=    fltarr(nsnaps)




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

	; load all the stellar info
        ; --------------------------
        N= fload_npart(2)+fload_npart(3)+fload_npart(4)
        m= 1.0e+10*fload_allstars_mass(1)
        age=fload_allstars_age(1)
        age=float(TTime-age)
        zmets=fload_allstars_z(1)


	; get the luminosities
	print, "load luminosities"
	load_all_stellar_luminosities, N, TTime, m, age, zmets, $
				Lum_Bol, Lum_U, Lum_B, Lum_V, Lum_R, Lum_I, Lum_J, Lum_H, Lum_K, $
				Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z, /notsolar

	; trap for any NaN's
	idx=where(finite(Lum_B) eq 0)
	if idx(0) ne -1 then begin
	    Lum_Bol(idx)= 100.0 & Lum_U(idx)= 100.0 & Lum_B(idx)= 100.0 & Lum_V(idx)= 100.0 & Lum_R(idx)= 100.0
	    Lum_I(idx)= 100.0 & Lum_J(idx)= 100.0 & Lum_H(idx)= 100.0 & Lum_K(idx)= 100.0
	    Lum_sdss_u(idx)= 100.0 & Lum_sdss_g(idx)= 100.0 & Lum_sdss_r(idx)= 100.0
	    Lum_sdss_i(idx)= 100.0 & Lum_sdss_z(idx)= 100.0
	    print, "Fixed ",n_elements(idx), " luminosities when trapping for NaN values."
	endif


	; convert to a total magnitude
	bolo[i]= -2.5*alog10(total(Lum_Bol)) 
	U_mag[i]= -2.5*alog10(total(Lum_U)) 
	B_mag[i]= -2.5*alog10(total(Lum_B)) 
	V_mag[i]= -2.5*alog10(total(Lum_V))
	R_mag[i]= -2.5*alog10(total(Lum_R)) 
	I_mag[i]= -2.5*alog10(total(Lum_I)) 
	J_mag[i]= -2.5*alog10(total(Lum_J)) 
	H_mag[i]= -2.5*alog10(total(Lum_H)) 
	K_mag[i]= -2.5*alog10(total(Lum_K)) 
	u_sdss_mag[i]= -2.5*alog10(total(Lum_sdss_u)) 
	g_sdss_mag[i]= -2.5*alog10(total(Lum_sdss_g)) 
	r_sdss_mag[i]= -2.5*alog10(total(Lum_sdss_r)) 
	i_sdss_mag[i]= -2.5*alog10(total(Lum_sdss_i)) 
	z_sdss_mag[i]= -2.5*alog10(total(Lum_sdss_z))



endfor
        
; ----------------------------------------
; write to file

openw, 1, frun+'/colors.txt', ERROR=err

printf, 1, "#   colors.txt"
printf, 1, "#   file contains: time, and total unobscured magnitudes in the following bands:"
printf, 1, "#                        bolometric, U, B, V, R, I, J, H, K, u, g, r, i ,z "
printf, 1, "#   "
printf, 1, "#   "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3," ", 14(F8.3," "))', $
                time[i], bolo[i], U_mag[i], B_mag[i], V_mag[i], R_mag[i], $
			I_mag[i], J_mag[i], H_mag[i], K_mag[i], u_sdss_mag[i], $
			g_sdss_mag[i], r_sdss_mag[i], i_sdss_mag[i], z_sdss_mag[i]
endfor
close, 1




; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end






;==================================================================================



pro time_mags_smalldt, frun


if not keyword_set(frun) then begin
	print, " "
	print, " time_mags, frun"
	print, " "
	print, " needs: .run determine_lums"
	print, "        .run sfr_multi"
	print, " "
	print, " "
	print, " "
	return
endif




;  get sfr info
;------------------
open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass


Ni= n_elements(sfrtime)-1

Time_Colors= fltarr(Ni)

bolo=     fltarr(Ni)
U_mag=    fltarr(Ni)
B_mag=    fltarr(Ni)
V_mag=    fltarr(Ni)
R_mag=    fltarr(Ni)
I_mag=    fltarr(Ni)
J_mag=    fltarr(Ni)
H_mag=    fltarr(Ni)
K_mag=    fltarr(Ni)
u_sdss_mag=    fltarr(Ni)
g_sdss_mag=    fltarr(Ni)
r_sdss_mag=    fltarr(Ni)
i_sdss_mag=    fltarr(Ni)
z_sdss_mag=    fltarr(Ni)


Stellar_Age= fltarr(Ni+1)
Stellar_Z= fltarr(Ni+1)
Stellar_Mass= fltarr(Ni+1)






; grab initial stellar mass
; ----------------------------
; open snapshot 0
;ok=fload_snapshot_bh(frun,snapi, /nopot_in_snap)
ok=fload_snapshot_bh(frun,0)
Stellar_Mass[0]= total(1.0e+10*fload_allstars_mass(1) / 0.7)
Stellar_Age[0]= -1.0
Stellar_Z[0]= 1.0/3.0     ; one-third solar





;-------------------------------------------

for i= 1, Ni do begin

        Ttime= 0.5*(sfrtime[i] + sfrtime[i-1])   ; gyr/h
        Time_Colors[i-1]= Ttime



        ;  stars
        ; =======
        dt= (sfrtime[i]-sfrtime[i-1]) * 1.0d+9 / 0.7     ; put in yr (take out h)
        Stellar_Mass[i]= 0.5*(sfrsfr[i]+sfrsfr[i-1]) * dt
        Stellar_Age[i]= Ttime - 1.0e-4
        Stellar_Z[i]= 0.5

        idx= where(Stellar_Mass gt 0.0)     ; don't waste time if no stars
        iN= n_elements(idx)
        iMass= Stellar_Mass(idx)
        iAge= Stellar_Age(idx)
        iZ= Stellar_Z(idx)

        itime= 1.0    ; never actually uses this
        iAge= (Ttime - iAge) / 0.7


        ; get the luminosities
        ;  - in units of solar luminosities
        print, "load luminosities"
        load_all_stellar_luminosities, iN, itime, iMass, iAge, iZ, $
                                Lum_Bol, Lum_U, Lum_B, Lum_V, Lum_R, Lum_I, Lum_J, Lum_H, Lum_K, $
                                Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z, /notsolar


	; trap for any NaN's
	idx=where(finite(Lum_B) eq 0)
	if idx(0) ne -1 then begin
	    Lum_Bol(idx)= 100.0 & Lum_U(idx)= 100.0 & Lum_B(idx)= 100.0 & Lum_V(idx)= 100.0 & Lum_R(idx)= 100.0
	    Lum_I(idx)= 100.0 & Lum_J(idx)= 100.0 & Lum_H(idx)= 100.0 & Lum_K(idx)= 100.0
	    Lum_sdss_u(idx)= 100.0 & Lum_sdss_g(idx)= 100.0 & Lum_sdss_r(idx)= 100.0
	    Lum_sdss_i(idx)= 100.0 & Lum_sdss_z(idx)= 100.0
	    print, "Fixed ",n_elements(idx), " luminosities when trapping for NaN values."
	endif


	; convert to a total magnitude
	bolo[i-1]= -2.5*alog10(total(Lum_Bol)) 
	U_mag[i-1]= -2.5*alog10(total(Lum_U)) 
	B_mag[i-1]= -2.5*alog10(total(Lum_B)) 
	V_mag[i-1]= -2.5*alog10(total(Lum_V))
	R_mag[i-1]= -2.5*alog10(total(Lum_R)) 
	I_mag[i-1]= -2.5*alog10(total(Lum_I)) 
	J_mag[i-1]= -2.5*alog10(total(Lum_J)) 
	H_mag[i-1]= -2.5*alog10(total(Lum_H)) 
	K_mag[i-1]= -2.5*alog10(total(Lum_K)) 
	u_sdss_mag[i-1]= -2.5*alog10(total(Lum_sdss_u)) 
	g_sdss_mag[i-1]= -2.5*alog10(total(Lum_sdss_g)) 
	r_sdss_mag[i-1]= -2.5*alog10(total(Lum_sdss_r)) 
	i_sdss_mag[i-1]= -2.5*alog10(total(Lum_sdss_i)) 
	z_sdss_mag[i-1]= -2.5*alog10(total(Lum_sdss_z))


endfor


        
; ----------------------------------------
; write to file

openw, 1, frun+'/colors_dt.txt', ERROR=err

printf, 1, "#   colors_dt.txt"
printf, 1, "#   file contains: time, and total unobscured magnitudes in the following bands:"
printf, 1, "#                        bolometric, U, B, V, R, I, J, H, K, u, g, r, i ,z "
printf, 1, "#   "
printf, 1, "#   "
for i=0,Ni-1 do begin
        printf, 1, FORMAT= '(F6.3," ", 14(F8.3," "))', $
                Time_Colors[i], bolo[i], U_mag[i], B_mag[i], V_mag[i], R_mag[i], $
			I_mag[i], J_mag[i], H_mag[i], K_mag[i], u_sdss_mag[i], $
			g_sdss_mag[i], r_sdss_mag[i], i_sdss_mag[i], z_sdss_mag[i]
endfor
close, 1




; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end









;==================================================================================






pro read_colors_file, frun, time, mags, $
			funky=funky

cfile= frun+'/colors.txt'
if keyword_set(funky) then cfile= frun+'/colors.txt.disk=1.0'

spawn, "wc "+cfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then c_data= fltarr(15,lines)

openr, 1, cfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, c_data
close, 1


time= c_data[0,*]
mags= c_data[1:14,*]

end




;==================================================================================




;-------------------------------------------------
;
;-------------------------------------------------
pro plot_color_color, frun


if not keyword_set(frun) then begin
        print, " "
        print, " plot_color_color, frun"
        print, " "
        print, " "
        return
endif

filename=frun+'/cc_diagram.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


read_colors_file, frun, time, mags

bolo= mags[0,*]
U= mags[1,*]
B= mags[2,*]
V= mags[3,*]
R= mags[4,*]
I= mags[5,*]
J= mags[6,*]
H= mags[7,*]
K= mags[8,*]
u_sdss= mags[9,*]
g_sdss= mags[10,*]
r_sdss= mags[11,*]
i_sdss= mags[12,*]
z_sdss= mags[13,*]


b_minus_v= B-V
u_minus_b= U-B

xaxistitle = "B-V"
yaxistitle = "U-B"
xmax = 1.2
xmin = -0.1
;ymax = -1.0
ymax = 0.0
;ymin = 0.8
ymin = 1.2

;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        ;/xlog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


merger_time= 0.91

; pre-merger
idx=where(time le (merger_time+0.1))
b_m_v_PM= b_minus_v(idx)
u_m_b_PM= u_minus_b(idx)
oplot, b_m_v_PM, u_m_b_PM, thick=3.0, psym=-2, color= 150

; post-merger
idx=where(time ge (merger_time-0.1))
b_m_v_PM= b_minus_v(idx)
u_m_b_PM= u_minus_b(idx)
oplot, b_m_v_PM, u_m_b_PM, thick=3.0, psym=-7, color= 50


; whole thing
;oplot, b_minus_v, u_minus_b, thick=3.0, psym=-2, color= 150


xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0


;  Done
; ------
device, /close




end







;==================================================================================



;-------------------------------------------------
;
;-------------------------------------------------
pro plot_Bmag_time, frun


if not keyword_set(frun) then begin
        print, " "
        print, " plot_Bmag_time, frun"
        print, " "
        print, " "
        return
endif

filename=frun+'/Bmag.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


read_colors_file, frun, time, mags

;bolo= mags[0,*]
;U= mags[1,*]
Bmag= mags[2,*]
;V= mags[3,*]
;R= mags[4,*]
;I= mags[5,*]
;J= mags[6,*]
;H= mags[7,*]
Kmag= mags[8,*]
;u_sdss= mags[9,*]
;g_sdss= mags[10,*]
;r_sdss= mags[11,*]
;i_sdss= mags[12,*]
;z_sdss= mags[13,*]



yaxistitle = "B,K"
xaxistitle = "Time (Gyr h!E-1!N)"
xmax = 3.0
xmin = 0.0
ymax = -25.0
ymin = -20.0

;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        ;/xlog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


;oplot, time, bolo, thick=3.0, psym=-3, color= 0
oplot, time, Bmag, thick=3.0, psym=-3, color= 50
oplot, time, Kmag, thick=3.0, psym=-3, color= 150

xyouts, 0.6, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0


;  Done
; ------
device, /close




end

;==================================================================================




;-------------------------------------------------
;  g-u color evolution during a merger
;   (motivated by figure 3. of Springel
;     di Matteo, Hernquist red gals)
;-------------------------------------------------
pro plot_color_time, junk


if not keyword_set(junk) then begin
        print, " "
        print, " plot_color_time, junk"
        print, " "
        print, " "
        return
endif

;filename=frun+'/color.eps'
filename='color.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;---------------------------

yaxistitle = "!6u-r"
;yaxistitle = "U-V"
;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"
xmax = 7.0
xmin = 0.0
ymax = 2.5
;ymax = 1.5
ymin = 0.3
;ymin = -0.5

;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        ;/xlog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


x0= 0.60

;process_onesim_gmu, "/raid4/tcox/vc3vc3", 4, x0, 0.48
;process_onesim_gmu, "/raid4/tcox/vc3rem_vc3", 6, x0, 0.44
;process_onesim_gmu, "/raid4/tcox/vc3rem_vc3rem", 3, x0, 0.40
;process_onesim_gmu, "/raid4/tcox/vc3rem_vc2", 2, x0, 0.36
;process_onesim_gmu, "/raid4/tcox/vc3rem_vc1", 7, x0, 0.32
;process_onesim_gmu, "/raid4/tcox/vc3rem_vc1a", 1, x0, 0.28
;process_onesim_gmu, "/raid4/tcox/vc3rem_vc1b", 5, x0, 0.24
;process_onesim_gmu, "/raid4/tcox/vc3rem_vc1c", 8, x0, 0.20

;process_onesim_gmu, "/raid4/tcox/vc3vc3e", 4, x0, 0.48
;process_onesim_gmu, "/raid4/tcox/vc3vc3e", 6, x0, 0.44, /funky


process_onesim_gmu, "/raid4/tcox/ds/vc3vc3e_2", 6, x0, 0.48, msg='std (with BH)'
process_onesim_gmu, "/raid4/tcox/ds/vc3vc3e_no", 7, x0, 0.44, msg='std (no BH)'

;process_onesim_gmu, "/raid4/tcox/sbw/sb8", 1, x0, 0.40, msg='!7g!6=2.0, !8v!6= 837'
;process_onesim_gmu, "/raid4/tcox/sbw/sb8", 1, x0, 0.40

;process_onesim_gmu, "/raid4/tcox/sbw/sb10", 2, x0, 0.36, msg='!7g!6=0.5, !8v!6= 837'
;;process_onesim_gmu, "/raid4/tcox/sbw/sb10", 2, x0, 0.36

;process_onesim_gmu, "/raid4/tcox/sbw/sb13", 3, x0, 0.32, msg='!7g!6=2.0, !8v!6= 209'
;process_onesim_gmu, "/raid4/tcox/sbw/sb13", 3, x0, 0.32

;process_onesim_gmu, "/raid4/tcox/sbw/sb9", 4, x0, 0.28, msg='!7g!6=0.05, !8v!6= 837'
;process_onesim_gmu, "/raid4/tcox/sbw/sb11", 5, x0, 0.24, msg='!7g!6=5.0, !8v!6= 837'



;  Done
; ------
device, /close




end


;==================================================================================


; do the dirty work
; --------------------
pro process_onesim_gmu, frun, pointselection, x0, y0, $
				funky=funky, msg=msg

read_colors_file, frun, time, mags, funky=funky


;bolo= mags[0,*]
;U= mags[1,*]
;Bmag= mags[2,*]
;V= mags[3,*]
;R= mags[4,*]
;I= mags[5,*]
;J= mags[6,*]
;H= mags[7,*]
;Kmag= mags[8,*]
u_sdss= mags[9,*]
;g_sdss= mags[10,*]
r_sdss= mags[11,*]
;i_sdss= mags[12,*]
;z_sdss= mags[13,*]


select_thispoint, pointselection, thispsym, thiscolor


time=time/0.7

u_sdss= u_sdss+0.2 + 0.02*time



oplot, time, u_sdss-r_sdss, thick=3.0, psym=-thispsym, color= thiscolor
;oplot, time, U-V, thick=3.0, psym=-symsel, color= symcolor

if not keyword_set(msg) then begin
	xyouts, x0, y0, fload_getid(frun), /normal, charthick=1, size=1.33, color=thiscolor
endif else begin
	xyouts, x0, y0, msg, /normal, charthick=1, size=1.33, color=thiscolor
endelse

;if not keyword_set(funky) then begin
;	xyouts, x0, y0, 'old=[0,5.0]', /normal, charthick=1, size=1.33, color=symcolor
;endif else begin
;	xyouts, x0, y0, 'old=[0,1.0]', /normal, charthick=1, size=1.33, color=symcolor
;endelse


end
