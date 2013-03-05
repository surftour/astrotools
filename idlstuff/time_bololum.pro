




pro time_bololum, frun, rotate_theta= rotate_theta, rotate_phi=rotate_phi, $
			fadd=fadd, snapnum=snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "time_bololum, frun, rotate_theta= rotate_theta, rotate_phi=rotate_phi"
   print, "  "
   print, "  "
   return
endif


if not keyword_set(fadd) then fadd= '' else fadd='_'+fadd
if not keyword_set(rotate_theta) then rotate_theta=0.0
if not keyword_set(rotate_phi) then rotate_phi= 0.0


	; ---------------------------------
	;  Fix the Band Information
	; ---------------------------------

	Bolometric_Solar_Luminosity= 3.826d+33    ; in ergs/sec

	; Wavelength(s) in Angstroms
	; --------------------------

	; Number of Bands
	NBands= 13

	Band= strarr(NBands)
	Band_Lambda= fltarr(NBands)
	Band_width= fltarr(NBands)
	Band_Solar_Luminosity= dblarr(NBands)

	Band[0]= "U"  & Band_Lambda[0] = 3650   &  Band_width[0]= 700   & Band_Solar_Luminosity[0]= 2.0d+32       ; U
	Band[1]= "B"  & Band_Lambda[1] = 4400   &  Band_width[1]= 980   & Band_Solar_Luminosity[1]= 5.2d+32       ; B
	Band[2]= "V"  & Band_Lambda[2] = 5500   &  Band_width[2]= 980   & Band_Solar_Luminosity[2]= 5.2d+32       ; V
	Band[3]= "R"  & Band_Lambda[3] = 6500   &  Band_width[3]= 1180  & Band_Solar_Luminosity[3]= 7.7d+32       ; R
	Band[4]= "I"  & Band_Lambda[4] = 8000   &  Band_width[4]= 1400  & Band_Solar_Luminosity[4]= 5.2d+32       ; I
	Band[5]= "J"  & Band_Lambda[5] = 12150  &  Band_width[5]= 2600  & Band_Solar_Luminosity[5]= 2.8d+32       ; J
	Band[6]= "H"  & Band_Lambda[6] = 16540  &  Band_width[6]= 2900  & Band_Solar_Luminosity[6]= 1.8d+32       ; H
	Band[7]= "K"  & Band_Lambda[7] = 21790  &  Band_width[7]= 4100  & Band_Solar_Luminosity[7]= 0.8d+32       ; K
	Band[8]= "u"  & Band_Lambda[8] = 3543   &  Band_width[8]= 650   & Band_Solar_Luminosity[8]= 2.0d+32       ; u
	Band[9]= "g"  & Band_Lambda[9] = 4770   &  Band_width[9]= 1480  & Band_Solar_Luminosity[9]= 5.2d+32       ; g
	Band[10]= "r" & Band_Lambda[10] = 6231  &  Band_width[10]= 1385 & Band_Solar_Luminosity[10]= 7.5d+32      ; r
	Band[11]= "i" & Band_Lambda[11] = 7625  &  Band_width[11]= 1565 & Band_Solar_Luminosity[11]= 5.2d+32      ; i
	Band[12]= "z" & Band_Lambda[12] = 9134  &  Band_width[12]= 1130 & Band_Solar_Luminosity[12]= 4.0d+32      ; z
	print, "----------------------------------------"
	print, "----------------------------------------"

min_Lambda= min(Band_Lambda - Band_width) - 100.0
max_Lambda= max(Band_Lambda - Band_width) + 100.0


; ---------------------------------

; this assumes they are ordered
;  0 through [something]


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snapshot* | wc ",result
nsnaps=long(result[0])



time= fltarr(nsnaps)

B_Lum= fltarr(nsnaps)
S_B_Lum= fltarr(nsnaps)
S_Bolo_Lum= fltarr(nsnaps)
S_AttBLum= fltarr(nsnaps)
BH_B_Lum= fltarr(nsnaps)
BH_Bolo_Lum= fltarr(nsnaps)
BH_AttBLum= fltarr(nsnaps)

Total_Band_Lums= fltarr(nsnaps,NBands)
Atten_Band_Lums= fltarr(nsnaps,NBands)

IR_Lum= fltarr(nsnaps)


; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------
snapi= 0

if keyword_set(snapnum) then begin
	snapi= snapnum
	nsnaps= 1
	fadd='_snap'+strcompress(string(snapnum),/remove_all)
endif

for si=0,nsnaps-1 do begin

        print, "--------------------------------------"


        ; open snapshot
	thistrial= 0
        repeat begin
            ;ok=fload_snapshot_bh(frun,snapi)
            ok=fload_snapshot_bh(frun,snapi, /nopot_in_snap)
            ;ok=fload_snapshot(frun,snapi)
            snapi= snapi+1
	    thistrial= thistrial+1
        endrep until ((ok eq 0) or (thistrial gt 200))

	if thistrial gt 200 then break

        ; what time is it?
        time[si]= fload_time(1)



	; grab luminosities
	; ----------------------
	ndisk= fload_npart(2)				; disk particles
	nbulge= fload_npart(3)				; bulge particles
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
                                Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z   ;, $
                                ;/notsolar

	; trap for any NaN's
	;idx=where(finite(Lum_B) eq 0)
	idx=where(finite(Lum_K) eq 0)
	if idx(0) ne -1 then begin
		Lum_Bol(idx)= 100.0 & Lum_U(idx)= 100.0 & Lum_B(idx)= 100.0 & Lum_V(idx)= 100.0 & Lum_R(idx)= 100.0
		Lum_I(idx)= 100.0 & Lum_J(idx)= 100.0 & Lum_H(idx)= 100.0 & Lum_K(idx)= 100.0
		Lum_sdss_u(idx)= 100.0 & Lum_sdss_g(idx)= 100.0 & Lum_sdss_r(idx)= 100.0
		Lum_sdss_i(idx)= 100.0 & Lum_sdss_z(idx)= 100.0
	print, "Fixed ",n_elements(idx), " luminosities when trapping for NaN values."
	endif

	print, "Bolometric luminosity= ", total(Lum_Bol)
	print, "Total K band luminosity= ", total(Lum_K)
	print, "Total B band luminosity= ", total(Lum_B)

	S_Bolo_Lum[si]= alog10(total(Lum_Bol))

	Band_Lum= fltarr(NBands,N)
        Band_Lum[0,*]= transpose(Lum_U)
        Band_Lum[1,*]= transpose(Lum_B)
        Band_Lum[2,*]= transpose(Lum_V)
        Band_Lum[3,*]= transpose(Lum_R)
        Band_Lum[4,*]= transpose(Lum_I)
        Band_Lum[5,*]= transpose(Lum_J)
        Band_Lum[6,*]= transpose(Lum_H)
        Band_Lum[7,*]= transpose(Lum_K)
        Band_Lum[8,*]= transpose(Lum_sdss_u)
        Band_Lum[9,*]= transpose(Lum_sdss_g)
        Band_Lum[10,*]= transpose(Lum_sdss_r)
        Band_Lum[11,*]= transpose(Lum_sdss_i)
        Band_Lum[12,*]= transpose(Lum_sdss_z)

	print, "Total Band_Lum= ", total(Band_Lum)

; ---------------------------------
;  Add BH Luminosity
; ---------------------------------
;include_BH_lum= 1
include_BH_lum= 0
if include_BH_lum eq 1 and fload_npart(5) gt 0 then begin

	; spectrum is just the current band
	; -----------------------------------
	N_spectrum= 10000   ; doesn't work well with low N
	c_light= 3.0e+8
	;log_nu_min= alog10(c_light/(Band_Lambda+0.5*Band_width)*1.0e+10)    ; 1e+10 converts to m from Angstroms
	;log_nu_max= alog10(c_light/(Band_Lambda-0.5*Band_width)*1.0e+10)
	log_nu_min= alog10(c_light/(max_Lambda)*1.0e+10)    ; 1e+10 converts to m from Angstroms
	log_nu_max= alog10(c_light/(min_Lambda)*1.0e+10)
	print, "max_lambda= ", max_lambda
	print, "min_lambda= ", min_lambda


	; get BH lum
	; -------------
        Nbh= fload_npart(5)
	bhids= fload_blackhole_id(1)
	bhbololum= fload_blackhole_lum(frun,Ttime,/bolometric)
	BH_Bolo_Lum[si]= alog10(total(bhbololum))

	OrigBand_Lum= Band_Lum
	NewBand_Lum= fltarr(NBands,N+Nbh)
	for bandi=0,NBands-1 do NewBand_Lum[bandi,Nbh:N+Nbh-1]= Band_Lum[bandi,*]
	Band_Lum= NewBand_Lum

	for i=1,fload_npart(5) do begin
		print, "Adding BH ",i, "  id= ",bhids[i-1]
		;BH_bolometric_lum= 12.0      ; bolometric luminosity in log
		BH_bolometric_lum= alog10(bhbololum[i-1])
		agn= fload_marconi_agn_spectrum(BH_bolometric_lum, 0.0, N_spectrum, log_nu_min, log_nu_max)
		nu= agn[0:N_spectrum-1]
		lambda_nu= c_light / (10.0^nu) * 1.0e+10
		l_nu= agn[N_spectrum:2*N_spectrum-1]

		; Add BH to list(s)
		; -------------------
		for bandi=0,NBands-1 do begin
			idx= where((lambda_nu gt Band_Lambda[bandi]-Band_width[bandi]) and $
					(lambda_nu lt Band_Lambda[bandi]+Band_width[bandi]))
			BH_lum_inband = mean(10.0^(l_nu(idx)) * 10.0^(nu(idx)))
			BH_lum_inband= BH_lum_inband * Bolometric_Solar_Luminosity / Band_Solar_Luminosity[bandi]
			;Band_Lum[bandi,*]= [BH_lum_inband, Band_Lum[bandi,*]]
			Band_Lum[bandi,i-1]= BH_lum_inband
			print, "Band= ", Band[bandi], "  L_bolo= ", 10.0^(BH_bolometric_lum), "  L_band= ", BH_lum_inband
		endfor

		;print, "BH Bolometric Luminosity = ",10.0^(BH_bolometric_lum)
		;print, "BH "+Band+"-band Luminosity (in bolometric solar units) = ",BH_lum_inband
		;print, "Bolometric Correction = ",10.0^(BH_bolometric_lum)/BH_lum_inband
		;BH_lum_inband= BH_lum_inband * Bolometric_Solar_Luminosity / Band_Solar_Luminosity
		;print, "BH "+Band+"-band Luminosity (in "+Band+"-band solar units) = ",BH_lum_inband

	endfor


endif


for bandi= 0, NBands-1 do Total_Band_Lums[si,bandi]= alog10(total(float(transpose(Band_Lum[bandi,*]))))





; ---------------------------------
;  Attenuate the Light
; ---------------------------------

attenuate_light= 0
;attenuate_light= 1
if keyword_set(attenuate_light) then begin


        Ngas= fload_npart(0)
        Nstars= fload_npart(2)+fload_npart(3)+fload_npart(4)
        Nbh= fload_npart(5)
if include_BH_lum eq 0 then Nbh= 0

        theta = (!PI / 180.0) * rotate_theta
        phi =(!PI / 180.0) * rotate_phi

        ; fields to pass
        Coord = fltarr(9,Ngas+Nstars+Nbh)

        ; gas
        Coord(0,0:Ngas-1) = fload_gas_xyz('x',center=[0,0,0])
        Coord(1,0:Ngas-1) = fload_gas_xyz('y',center=[0,0,0])
        Coord(2,0:Ngas-1) = fload_gas_xyz('z',center=[0,0,0])

        Coord(3,0:Ngas-1) = fload_gas_u(1)
        Coord(4,0:Ngas-1) = fload_gas_rho(1)
        Coord(5,0:Ngas-1) = fload_gas_hsml(1)
        Coord(6,0:Ngas-1) = fload_gas_numh(1)
        Coord(7,0:Ngas-1) = fload_gas_nume(1)
        Coord(8,0:Ngas-1) = fload_gas_metallicity(1)

        ; stars
        Coord(0,Ngas:Nstars+Ngas-1) = fload_allstars_xyz('x',center=[0,0,0])
        Coord(1,Ngas:Nstars+Ngas-1) = fload_allstars_xyz('y',center=[0,0,0])
        Coord(2,Ngas:Nstars+Ngas-1) = fload_allstars_xyz('z',center=[0,0,0])

        ; black holes
	if Nbh gt 0 then begin
        	Coord(0,Ngas+Nstars:Nstars+Ngas+Nbh-1) = fload_blackhole_xyz('x',center=[0,0,0])
        	Coord(1,Ngas+Nstars:Nstars+Ngas+Nbh-1) = fload_blackhole_xyz('y',center=[0,0,0])
        	Coord(2,Ngas+Nstars:Nstars+Ngas+Nbh-1) = fload_blackhole_xyz('z',center=[0,0,0])
	endif

        los_NH= fltarr(Nstars+Nbh)
        los_Z= fltarr(Nstars+Nbh)

        print, "PASSING: "
        print, "N_gas= ", Ngas
        print, "N_stars= ", Nstars
        print, "N_bh= ", Nbh
        print, "theta= ", theta
        print, "phi= ", phi
        help, Coord

	S = CALL_EXTERNAL('/n/home03/tcox/Tools/C-Routines_for_IDL/LOSColumn_RhoZ/getnh.so', $
                'getnh', $
                Ngas, $ 
                Nstars, $
                Nbh, $ 
                theta, $
                phi, $ 
                Coord, $
                los_NH, $
		los_Z)

        ; trap for really low NH values
        idx= where(los_NH lt 1.0e+10)
        if idx(0) ne -1 then los_NH(idx)= 1.0e+10
        idx= where(los_Z le 0.0)
        if idx(0) ne -1 then los_Z(idx)= 1.0e-5

	;  OK, now we've got the N_h values,
	; what is the attenuated luminosity

	print, "Now calling attenuation program."



	for bandi=0,NBands-1 do begin
		; pass in anstroms, but it gets converted
		; to microns in ism_absorption
		iBand_Lambda= float(Band_Lambda[bandi])
		print, "Band_Lambda= ",iBand_Lambda

	        AttLums= fltarr(Nstars+Nbh)      ;attenuated luminosities
		PreLums= float(transpose(Band_Lum[bandi,*]))

	        S = CALL_EXTERNAL('/n/home03/tcox/Tools/C-Routines_for_IDL/ISMAbsorption/ism_absorption', $
	                'calc_atten_lum', $
	                Nstars+Nbh, $
	                iBand_Lambda, $
	                ;ProjectThisLum, $
	                PreLums, $
	                los_NH, $
	                los_Z, $
	                AttLums)


		Total_Band_Lums[si,bandi]= alog10(total(PreLums))
		Atten_Band_Lums[si,bandi]= alog10(total(AttLums))


print, "total light= ", total(PreLums), total(PreLums)
print, "attenuated light= ", total(AttLums)
print, "    abs. f= ", total(AttLums)/total(PreLums)

print, " "
print, " Stars "
print, " ------"
print, "att. light = ", total(AttLums[Nbh:Nbh+Nstars-1])
print, "     abs. f= ", total(AttLums[Nbh:Nbh+Nstars-1])/total(PreLums[Nbh:Nbh+Nstars-1])

if Nbh gt 0 then begin
	print, " "
	print, " BH  "
	print, " ------"
	print, "att. light = ", total(AttLums[0:Nbh-1])
	print, "     abs. f= ", total(AttLums[0:Nbh-1])/total(PreLums[0:Nbh-1])
	print, "     los_NH= ", los_NH[0:Nbh-1]
	print, "     los_Z =  ", los_Z[0:Nbh-1]
	print, " "
endif

		if Band[bandi] eq 'B' then begin
			B_Lum[si]= alog10(total(PreLums))
			S_B_Lum[si]= alog10(total(PreLums[Nbh:Nbh+Nstars-1]))
			S_AttBLum[si]= alog10(total(AttLums[Nbh:Nbh+Nstars-1]))
			if Nbh gt 0 then begin
				BH_B_Lum[si]= alog10(total(PreLums[0:Nbh-1]))
				BH_AttBLum[si]= alog10(total(AttLums[0:Nbh-1]))
			endif else begin
				BH_B_Lum[si]= 0.0
				BH_AttBLum[si]= 0.0
			endelse
		endif
	endfor

	irlum= 10^(Total_Band_Lums[si,*]) - 10^(Atten_Band_Lums[si,*])
	IR_Lum[si]= alog10(total(irlum))

	print, "total Total_Band_Lums= ",total(10^Total_Band_Lums[si,*])
	print, "total Atten_Band_Lums= ",total(10^Atten_Band_Lums[si,*])
	print, "total IR_Lum= ",10^IR_Lum[si]

endif



endfor


; ------------------------------
;  Done - now write text file
; ------------------------------

openw, 1, frun+'/lum_att_B'+fadd+'.txt', ERROR=err

printf, 1, "#   lum_B_att.txt (+fadd, maybe)"
printf, 1, "#   B-band luminosity of stars and blackhole, from one direction, with and without attentuation"
printf, 1, "#        "
printf, 1, "# time      total    stellar    stellar         BH         BH "
printf, 1, "#(Gyrh)luminosity luminosity   att.lum. luminosity    att.lum."
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"  ",5(F9.4,"  "))', $
                time[i], B_Lum[i], S_B_Lum[i], S_AttBLum[i], BH_B_Lum[i], BH_AttBLum[i]
endfor
close, 1



openw, 1, frun+'/lum_bolo'+fadd+'.txt', ERROR=err

printf, 1, "#   lum_bolo.txt (+fadd, maybe)"
printf, 1, "#   Bolometric luminosity of stars and blackhole"
printf, 1, "#        "
printf, 1, "# time    stellar        BH  "
printf, 1, "#(Gyrh)luminosity luminosity "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"  ",2(F9.4,"  "))', $
                time[i], S_Bolo_Lum[i], BH_Bolo_Lum[i]
endfor
close, 1






;  total lums in all bands
; --------------------------
openw, 1, frun+'/lum_tot_allbnds'+fadd+'.txt', ERROR=err

printf, 1, "#   lum_tot_allbnds.txt (+fadd, maybe)"
printf, 1, "#   Total (unobscured) luminosity of stars and blackhole, from one direction"
printf, 1, "#        "
printf, 1, "# time      "
printf, 1, "#(Gyr/h)    U          B          V          R          I          J          H          K          u          g          r          i          z"
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"  ",13(F9.4,"  "))', $
                time[i], Total_Band_Lums[i,0], Total_Band_Lums[i,1], Total_Band_Lums[i,2], $
                	Total_Band_Lums[i,3], Total_Band_Lums[i,4], Total_Band_Lums[i,5], $
                	Total_Band_Lums[i,6], Total_Band_Lums[i,7], Total_Band_Lums[i,8], $
                	Total_Band_Lums[i,9], Total_Band_Lums[i,10], Total_Band_Lums[i,11], $
                	Total_Band_Lums[i,12]
endfor
close, 1



;  attenuated lums in all bands
; --------------------------
openw, 1, frun+'/lum_att_allbnds'+fadd+'.txt', ERROR=err

printf, 1, "#   lum_att_allbnds.txt (+fadd, maybe)"
printf, 1, "#   Attenuated luminosity of stars and blackhole, from one direction"
printf, 1, "#        "
printf, 1, "# time      "
printf, 1, "#(Gyr/h)    U          B          V          R          I          J          H          K          u          g          r          i          z"
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"  ",13(F9.4,"  "))', $
                time[i], Atten_Band_Lums[i,0], Atten_Band_Lums[i,1], Atten_Band_Lums[i,2], $
                	Atten_Band_Lums[i,3], Atten_Band_Lums[i,4], Atten_Band_Lums[i,5], $
                	Atten_Band_Lums[i,6], Atten_Band_Lums[i,7], Atten_Band_Lums[i,8], $
                	Atten_Band_Lums[i,9], Atten_Band_Lums[i,10], Atten_Band_Lums[i,11], $
                	Atten_Band_Lums[i,12]
endfor
close, 1


; Now, the total IR lum
; -----------------------
openw, 1, frun+'/lum_IR'+fadd+'.txt', ERROR=err

printf, 1, "#   lum_IR.txt (+fadd, maybe)"
printf, 1, "#   IR luminosity, assuming all absorbed radiation comes out in IR"
printf, 1, "#        "
printf, 1, "# time      total IR "
printf, 1, "#(Gyr/h)  luminosity "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"     ",1(F9.4,"  "))', $
                time[i], IR_Lum[i]
endfor
close, 1



end








; ================================================================================










; --------------------------------
;  Read Attenuated Lum File
; ----------------------------------
pro read_attlum_file, frun, time, bololum, s_bololum, s_attlum, bh_bololum, bh_attlum, $
				attfile=attfile

if keyword_set(attfile) then attfile= frun+'/'+attfile else attfile= frun+'/lum_att.txt'

spawn, "wc "+attfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then re_data= fltarr(6,lines)

openr, 1, attfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, re_data
close, 1


time= re_data[0,*]
bololum= re_data[1,*]
s_bololum= re_data[2,*]
s_attlum= re_data[3,*]
bh_bololum= re_data[4,*]
bh_attlum= re_data[5,*]


end






; --------------------------------
;  Read IR Lum File
; ----------------------------------
pro read_irlum_file, frun, time, irlum, $
		irfile=irfile

if not keyword_set(irfile) then irfile= frun+'/lum_IR.txt'

spawn, "wc "+irfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then re_data= fltarr(2,lines)

openr, 1, irfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, re_data
close, 1


time= re_data[0,*]
irlum= re_data[1,*]


end







; --------------------------------
;  Read Bolometric Lum File
; ----------------------------------
pro read_bololum_file, frun, time, sbololum, bhbololum, $
                bolofile=bolofile

if not keyword_set(bolofile) then bolofile= frun+'/lum_bolo.txt'

spawn, "wc "+bolofile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then re_data= fltarr(3,lines)

openr, 1, bolofile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, re_data
close, 1


time= re_data[0,*]
sbololum= re_data[1,*]
bhbololum= re_data[2,*]


end











; ================================================================================








;-------------------------------------------------
;
;-------------------------------------------------
pro plot_lumf, frun


if not keyword_set(frun) then begin
        print, " "
        print, " plot_lumf, frun"
        print, " "
        print, " "
        return
endif

filename='lumf.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;read_attlum_file, frun, time, bololum, s_bololum, s_attlum, bh_bololum, bh_attlum, attfile='lum_att_B.txt'
;S_bololum= s_attlum
;BH_bololum= bh_attlum
;bololum= alog10(10^(S_bololum) + 10^(BH_bololum))
;attlum= alog10(10^(s_attlum) + 10^(bh_attlum))

read_bololum_file, frun, time, sbololum, bhbololum, bolofile=bolofile
S_bololum= sbololum
BH_bololum= bhbololum
bololum= alog10(10^(sbololum) + 10^(bhbololum))

tmerg= 1.1


;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"
time= time/0.7
xmax = 2.0
xmin = 0.0

;yaxistitle = "!6Log L (L!DB,!9n!6!N)"
yaxistitle = "!6Log L (L!D!9n!6!N)"
;yaxistitle = "!6Log L!Dattenuated!N (L!D!9n!6!N)"
;yaxistitle = "!6Log L!Dbolometric!N (L!D!9n!6!N)"
ymax = 12.2
ymin = 10.6


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
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


        ;oplot, time, bololum, thick=4.0, color= 0, psym=-3
        ;oplot, time, S_bololum, thick=4.0, color= 150, psym=-3, linestyle=1
        oplot, time, S_bololum, thick=4.0, color= 150, psym=-3, linestyle=0
        ;oplot, time, BH_bololum, thick=4.0, color= 50, psym=-3, linestyle=2

        ;oplot, time, bololum, thick=4.0, color= 0, psym=-3
        ;oplot, time, S_attlum, thick=4.0, color= 150, psym=-3, linestyle=1
        ;oplot, time, BH_attlum, thick=4.0, color= 50, psym=-3, linestyle=2


	;xyouts, 0.76, 0.90, "Total", /normal, color= 0, charthick=3.0, charsize=1.3
	xyouts, 0.76, 0.86, "Stars", /normal, color= 150, charthick=3.0, charsize=1.3
	;xyouts, 0.76, 0.82, "BH", /normal, color= 50, charthick=3.0, charsize=1.3


        ;oplot, time, bololum, thick=4.0, color= 0, psym=-3
        ;oplot, time, attlum, thick=4.0, color= 150, psym=-3, linestyle=1
	;xyouts, 0.26, 0.23, "Bolometric", /normal, color= 0, charthick=3.0, charsize=1.3
	;xyouts, 0.26, 0.19, "Obscured", /normal, color= 150, charthick=3.0, charsize=1.3
	;xyouts, 0.45, 0.19, "(from one viewing angle)", /normal, color= 150, charthick=3.0, charsize=1.0
	
	xyouts, 0.22, 0.90, fload_fid(frun), /normal, color= 0, charthick=3.0, charsize=1.0

device, /close




end



;-------------------------------------------------
;
;-------------------------------------------------
pro att_lumf, frun


if not keyword_set(frun) then begin
        print, " "
        print, " att_lumf, frun"
        print, " "
        print, " "
        return
endif

filename='lumf.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;read_attlum_file, frun, time, bololum, s_bololum, s_attlum, bh_bololum, bh_attlum
read_bololum_file, frun, time, s_bololum, bh_bololum

tmerg= 1.1


xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xmax = 2.0
xmin = 0.0

yaxistitle = "Fraction of L!DBolometric!N from BH"
ymax = 1.0
ymin = 0.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
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

	Total_B_lum= 10^(S_bololum) + 10^(BH_bololum)
	BH_B_f= 10^(BH_bololum) / Total_B_lum

        oplot, time, BH_B_f, thick=4.0, color= 50, psym=-3


	; asterisks at times
	; --------------------
	write_asterisks= 0
	if write_asterisks eq 1 then begin
		thistime= 0.45
		idx= where(time eq thistime)
		thisf= BH_B_f(idx)
		oplot, [thistime], [thisf], thick=4.0, color=150, psym=2

		thistime= 0.85
		idx= where(time eq thistime)
		thisf= BH_B_f(idx)
		oplot, [thistime], [thisf], thick=4.0, color=150, psym=2

		thistime= 1.105
		idx= where(time eq thistime)
		thisf= BH_B_f(idx)
		oplot, [thistime], [thisf], thick=4.0, color=150, psym=2

		thistime= 1.215
		idx= where(time eq thistime)
		thisf= BH_B_f(idx)
		oplot, [thistime], [thisf], thick=4.0, color=150, psym=2

		thistime= 1.205
		idx= where(time eq thistime)
		thisf= BH_B_f(idx)
		oplot, [thistime], [thisf], thick=4.0, color=150, psym=2

		thistime= 1.310
		idx= where(time eq thistime)
		thisf= BH_B_f(idx)
		oplot, [thistime], [thisf], thick=4.0, color=150, psym=2
	endif


	; old stuff
	; -----------
        ;oplot, time, BH_B_f, thick=4.0, color= 150, psym=-3
	;xyouts, 0.22, 0.88, "Unobscured", /normal, color= 150, charthick=3.0

	;Obs_B_lum= 10^(S_attlum) + 10^(BH_attlum)
	;BH_ObsB_f= 10^(BH_attlum) / Obs_B_lum

        ;oplot, time, BH_ObsB_f, thick=4.0, color= 50, psym=-3
	;xyouts, 0.22, 0.82, "Observed", /normal, color= 50, charthick=3.0


device, /close




end



;-------------------------------------------------
;
;-------------------------------------------------
pro att_f, frun


if not keyword_set(frun) then begin
        print, " "
        print, " att_f, frun"
        print, " "
        print, " "
        return
endif

filename='lumf.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;read_attlum_file, frun, time, bololum, s_bololum, s_attlum, bh_bololum, bh_attlum
read_bololum_file, frun, time, s_bololum, bh_bololum

tmerg= 1.1


xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xmax = 2.0
xmin = 0.0

;yaxistitle = "Attenuated Luminosity from BH"
;ymax = 1.0
;ymin = 0.0
;yaxistitle = "Attenuated Luminosity (~L!DIR!N)"
yaxistitle = "!6Bolometric Luminosity (L!D!9n!6!N)"
ymax = 12.9
ymin = 9.4


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
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


	; attlum file
	; ---------------
	;Atten_Star_B_lum= 10^(S_bololum) - 10^(S_attlum)
	;Atten_BH_B_lum= 10^(BH_bololum) - 10^(BH_attlum)

	;Total_Atten_B_lum= Atten_Star_B_lum + Atten_BH_B_lum

	;Total_Atten_B_lum= alog10(Total_Atten_B_lum)
	;oplot, time, Total_Atten_B_lum, thick=4.0, color= 0, psym=-3

	;BH_Atten_B_f= Atten_BH_B_lum / Total_Atten_B_lum
        ;oplot, time, BH_Atten_B_f, thick=4.0, color= 50, psym=-3


	; bolo lum file
	; ----------------
	TotalLum= 10^(s_bololum) + 10^(bh_bololum)
	TotalLum= alog10(TotalLum)
	oplot, time, TotalLum, thick=4.0, color= 0, psym=-3

	oplot, time, s_bololum, thick=4.0, color= 50, psym=-3, linestyle=1

	oplot, time, bh_bololum, thick=4.0, color= 150, psym=-3, linestyle=2


	; key
	; -----
	oplot, [0.1,0.3], [12.5,12.5], thick=4.0, color=0, psym=-3
	xyouts, 0.32, 0.85, 'total', color=0, charthick=2.0, size=1.33, /normal
	oplot, [0.1,0.3], [12.3,12.3], thick=4.0, color=50, psym=-3, linestyle=1
	xyouts, 0.32, 0.805, 'stellar', color=50, charthick=2.0, size=1.33, /normal
	oplot, [0.1,0.3], [12.1,12.1], thick=4.0, color=150, psym=-3, linestyle=2
	xyouts, 0.32, 0.76, 'black hole', color=150, charthick=2.0, size=1.33, /normal





	; arrows
	; -------------
	;thistime= 1.105
	;idx= where(time eq thistime)
	;thislum= TotalLum(idx)
	;;oplot, [thistime], [thisf], thick=4.0, color=150, psym=2
	;thislum= 12.3
	;arrow, thistime, thislum+0.5, thistime, thislum+0.1, color= 200, thick=3.0, hthick=4.0, /data

	; vc3vc3e_2
	xyouts, 0.87, 0.87, '(e)', color=0, charthick=2.0, size=1.33, /normal
	xyouts, 0.82, 0.82, 'no BH', color=0, charthick=2.0, size=1.33, /normal
	draw_arrow, 1.070, time, TotalLum
	draw_arrow, 1.105, time, TotalLum
	draw_arrow, 1.125, time, TotalLum
	draw_arrow, 1.145, time, TotalLum
	draw_arrow, 1.165, time, TotalLum
	draw_arrow, 1.185, time, TotalLum
	draw_arrow, 1.205, time, TotalLum
	draw_arrow, 1.225, time, TotalLum
	draw_arrow, 1.245, time, TotalLum
	draw_arrow, 1.265, time, TotalLum
	draw_arrow, 1.285, time, TotalLum
	draw_arrow, 1.310, time, TotalLum


	; vc3vc3h_2
	;xyouts, 0.87, 0.87, '(h)', color=0, charthick=2.0, size=1.33, /normal
	;draw_arrow, 0.900, time, TotalLum
	;draw_arrow, 0.950, time, TotalLum
	;draw_arrow, 1.030, time, TotalLum
	;draw_arrow, 1.070, time, TotalLum
	;draw_arrow, 1.105, time, TotalLum
	;draw_arrow, 1.125, time, TotalLum
	;draw_arrow, 1.145, time, TotalLum
	;draw_arrow, 1.165, time, TotalLum
	;draw_arrow, 1.185, time, TotalLum
	;draw_arrow, 1.205, time, TotalLum
	;draw_arrow, 1.225, time, TotalLum
	;draw_arrow, 1.245, time, TotalLum



device, /close




end







pro draw_arrow, thistime, time, yquant

        ; arrows
        ; -------------
        ;thistime= 1.105
        idx= where(time eq thistime)
        thislum= yquant(idx)
        ;oplot, [thistime], [thisf], thick=4.0, color=150, psym=2
        thislum= 12.3
        arrow, thistime, thislum+0.5, thistime, thislum+0.1, color= 200, thick=3.0, hthick=4.0, /data

end




;============================================================================
;
;  We adapt this from the mags reading program.
;  Note that this reads in luminosities and then
;  converts to magnitudes.
;
;
;============================================================================

pro read_lums_file, frun, time, lums, mags=mags, $
                        totlum=totlum, $
			attlum=attlum, $
			fadd=fadd

cfile= frun+'/lum_tot_allbnds.txt'
if keyword_set(attlum) then cfile= frun+'/lum_att_allbnds.txt'
if keyword_set(fadd) then cfile= frun+'/lum_att_allbnds_'+fadd+'.txt'

spawn, "wc "+cfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then c_data= fltarr(14,lines)

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
lums= c_data[1:13,*]


solar_mags= [5.60, 5.51, 4.84, 4.48, 4.13, 3.70, $
		3.37, 3.33, 6.28, 4.95, 4.45, 4.35, 4.36]
if keyword_set(mags) then begin
    mags= lums
    for i=0,lines-1 do begin
	;mags[*,i]= -2.5*alog10(lums[*,i]) + solar_mags
	mags[*,i]= -2.5*lums[*,i] + solar_mags
    endfor
endif

end












;============================================================================
;============================================================================


;============================================================================
;============================================================================






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

;filename=frun+'/cc_diagram.eps'
filename='cc_diagram.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


xaxistitle = "B-V"
yaxistitle = "U-B"
xmax = 1.2
xmin = -0.1
ymax = -1.0
;ymax = 0.0
ymin = 0.8
;ymin = 1.2

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


do_one_color_color, "/raid4/tcox/vc3vc3e", linecolor=50

do_one_color_color, "/raid4/tcox/vc3vc3e", linecolor=150, /attlum


;  Done
; ------
device, /close




end







pro do_one_color_color, frun, linecolor=linecolor, $
			attlum=attlum

if not keyword_set(linecolor) then linecolor= 50

mags=1
read_lums_file, frun, time, lums, mags=mags, attlum=attlum

U= mags[0,*]
B= mags[1,*]
V= mags[2,*]
R= mags[3,*]
I= mags[4,*]
J= mags[5,*]
H= mags[6,*]
K= mags[7,*]
u_sdss= mags[8,*]
g_sdss= mags[9,*]
r_sdss= mags[10,*]
i_sdss= mags[11,*]
z_sdss= mags[12,*]


b_minus_v= B-V
u_minus_b= U-B

merger_time= 1.10

; pre-merger
;idx=where(time le (merger_time+0.1))
;b_m_v_PM= b_minus_v(idx)
;u_m_b_PM= u_minus_b(idx)
;oplot, b_m_v_PM, u_m_b_PM, thick=3.0, psym=-2, color= 150

; post-merger
;idx=where(time ge (merger_time-0.1))
;b_m_v_PM= b_minus_v(idx)
;u_m_b_PM= u_minus_b(idx)
;oplot, b_m_v_PM, u_m_b_PM, thick=3.0, psym=-7, color= 50


; whole thing
oplot, b_minus_v, u_minus_b, thick=3.0, psym=-2, color= linecolor


;xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0


end










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





;-------------------------------------------------
;  g-u color evolution during a merger
;   (motivated by figure 3. or springel
;     di Matteio, Hernquist red gals)
;-------------------------------------------------
pro plot_gmu_time, junk


if not keyword_set(junk) then begin
        print, " "
        print, " plot_gmu_time, junk"
        print, " "
        print, " "
        return
endif

;filename=frun+'/gmu.eps'
filename='gmu.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;---------------------------

yaxistitle = "u-r"
;xaxistitle = "Time (Gyr h!E-1!N)"
xaxistitle = "Time (Gyr)"
xmax = 4.3
xmin = 0.0

;yaxistitle = "u-r"
;ymax = -25.0
;ymin = -20.0

;yaxistitle = "U-V"
;ymax = 1.2
;ymin = -0.8

;yaxistitle = "J-K"
;ymax = 0.9
;ymin = 0.3

yaxistitle = "U-B"
yaxistitle = "U-B (unattenuated)"
ymax = 0.6
ymin =-0.9



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


x0= 0.65

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

;process_onesim_gmu, "/raid4/tcox/vc3vc3e_2", 4, x0, 0.48
;process_onesim_gmu, "/raid4/tcox/vc3vc3e_2", 6, x0, 0.44



; green vally
polyfill, [xmin,xmin,xmax,xmax,xmin], [0.0,0.25,0.25,0.0,0.0], color= 100, /line_fill, linestyle= 1
xyouts, 0.25, 0.68, 'Green Valley', /normal, charthick=4.0, size=1.8, color=100


process_onesim_gmu, "data1/lblecha/q1_fg0.4_tc_nokick", 4, x0, 0.48
xyouts, 0.6, 0.40, '!18v!Dk!N/v!Desc!N!3= 0', /normal, charthick=1.5, size=1.33, color=0

process_onesim_gmu, "data1/lblecha/q1_fg0.4_tc_v0.9", 6, x0, 0.43
xyouts, 0.6, 0.35, '!18v!Dk!N/v!Desc!N!3= 0.9', /normal, charthick=1.5, size=1.33, color=150

process_onesim_gmu, "data1/lblecha/q1_fg0.4_nobh", 5, x0, 0.43
xyouts, 0.6, 0.45, 'No BH', /normal, charthick=1.5, size=1.35, color=50


;  Done
; ------
device, /close




end




; do the dirty work
; --------------------
pro process_onesim_gmu, frun, pointselection, x0, y0, $
				funky=funky, attenlums=attenlums

if keyword_set(attenlums) then begin
	read_colors_file, frun, time, mags, funky=funky
	;bolo= mags[0,*]
	;U= mags[1,*]
	;Bmag= mags[2,*]
	;V= mags[3,*]
	;R= mags[4,*]
	;I= mags[5,*]
	J= mags[6,*]
	U= J
	;H= mags[7,*]
	Kmag= mags[8,*]
	V= Kmag
	;u_sdss= mags[9,*]
	;g_sdss= mags[10,*]
	;r_sdss= mags[11,*]
	;i_sdss= mags[12,*]
	;z_sdss= mags[13,*]
endif else begin
	mags=1
	read_lums_file, frun, time, lums, mags=mags, attlum=attlum  ; by default the total (non-atten) file is opened
	U= mags[0,*]
	B= mags[1,*]
	help, U, B
	print, "U_final= ", U[n_elements(U)-1]
	print, "B_final= ", B[n_elements(B)-1]
	;V= mags[2,*]
	;J= mags[5,*]
	;U= J
        ;Kmag= mags[7,*]
	;V= Kmag
endelse


time= time/0.7


;  open (or filled) square
; --------------------------
if pointselection eq 1 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 8
        symcolor= 150
endif   

;  x's
; -----
if pointselection eq 2 then begin
        symsel= 7
        symcolor= 100
endif   

;  open circle
; --------------
if pointselection eq 3 then begin
        symsize= 1.0
        ;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
        symsel= 8
        symcolor= 50
endif   

;  open triangle
; ----------------
if pointselection eq 4 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
        usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
        symsel= 3
        symcolor= 0
	thislinestyle= 2
endif

;  asterisk
; ----------
if pointselection eq 5 then begin
        symsel= 2
        symsel= 3
        symcolor= 50
	thislinestyle= 1
endif

;  filled square
; --------------------------
if pointselection eq 6 then begin
        symsize= 1.0
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 3
        symcolor= 150
endif


;  diamond
; ----------
if pointselection eq 7 then begin
        symsel= 3
        symcolor= 150
	thislinestyle= 3
endif

;  triangle (open)
; -----------------
if pointselection eq 8 then begin
        symsel= 5
        symcolor= 180
endif

;  plus sign
; ----------
if pointselection eq 9 then begin
        symsel= 1
        symcolor= 170
endif





;oplot, time, u_sdss-r_sdss, thick=3.0, psym=-symsel, color= symcolor
;oplot, time, U-V, thick=3.0, psym=-symsel, color= symcolor
oplot, time, U-B, thick=3.0, psym=-symsel, color= symcolor, linestyle= thislinestyle

;xyouts, x0, y0, fload_getid(frun), /normal, charthick=1, size=1.33, color=symcolor

;if not keyword_set(funky) then begin
;	xyouts, x0, y0, 'old=[0,5.0]', /normal, charthick=1, size=1.33, color=symcolor
;endif else begin
;	xyouts, x0, y0, 'old=[0,1.0]', /normal, charthick=1, size=1.33, color=symcolor
;endelse


if not keyword_set(attenlums) then begin
	;xyouts, x0, y0, 'Std', /normal, charthick=1, size=1.33, color=symcolor
endif else begin
	;xyouts, x0, y0, 'Obscured', /normal, charthick=1, size=1.33, color=symcolor
endelse


end







;===============================================================================








; ----------------------------------
;  Averge IR Lum from 3 directions
; ----------------------------------
pro avg_irlum_fromxyz, frun

if not keyword_set(frun) then irfile= frun+'/lum_IR.txt'


read_irlum_file, "junk", time, irlum_xy, irfile=frun+"/lum_IR_xy.txt"
read_irlum_file, "junk", time, irlum_yz, irfile=frun+"/lum_IR_yz.txt"
read_irlum_file, "junk", time, irlum_xz, irfile=frun+"/lum_IR_xz.txt"

nsnaps= n_elements(time)

irlum= (10^irlum_xy + 10^irlum_yz + 10^irlum_xz) / 3.0
irlum= alog10(irlum)


; Now, the total IR lum
; -----------------------
openw, 1, frun+'/lum_IR_avg.txt', ERROR=err

printf, 1, "#   lum_IR.txt "
printf, 1, "#   IR luminosity, averaged from xy, yz, and xz"
printf, 1, "#        "
printf, 1, "# time      total IR "
printf, 1, "#(Gyr)    luminosity "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"     ",1(F9.4,"  "))', $
                time[i], irlum[i]
endfor
close, 1



end







;===============================================================================

function get_stellar_lum, N, TTime, m, age, zmets


        ; get the luminosities
        ;  - in units of solar luminosities
        print, "load luminosities"
        load_all_stellar_luminosities, N, TTime, m, age, zmets, $
                                Lum_Bol, Lum_U, Lum_B, Lum_V, Lum_R, Lum_I, Lum_J, Lum_H, Lum_K, $
                                Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z   ;, $
                                ;/notsolar

        ; trap for any NaN's
        idx=where(finite(Lum_Bol) eq 0)
        if idx(0) ne -1 then begin
                Lum_Bol(idx)= 100.0
        	print, "Fixed ",n_elements(idx), " luminosities when trapping for NaN values."
        endif

	return, total(Lum_Bol)


end







; ----------------------------------
;  v2
; ----------------------------------

pro v2, frun



if not keyword_set(frun) then begin
   print, "  "
   print, "v2, frun"
   print, "  "
   print, "   needs:  .run determine_lums "
   print, "           .run sfr_multi"
   print, "           .run bh_multi "
   print, "  "
   return
endif



;  get sfr info
;------------------
open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass


Ni= n_elements(sfrtime)-1

Time_Lum= fltarr(Ni)
BH_Bolo_Lum= fltarr(Ni)
S_Bolo_Lum= fltarr(Ni)


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




;  get bh info
;------------------
open_blackhole_txt, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, $
                                bh_totalmass, bh_mdot_edd

Ni_bh= n_elements(bhtime)

; cgs units, cm, sec, g
L_solar= 3.826d+33                               ; in ergs per sec
cc= 2.9979d+10                                 ; speed of light in cm/sec
convert_sunyr_gsec= 6.30428d+25                  ; convert msun/yr -> g/sec
    
boloL_ergss= 0.1*bh_mdot_sunyr*cc*cc*convert_sunyr_gsec
boloL_sun= boloL_ergss / L_solar



;-------------------------------------------

for i= 1, Ni do begin

	Ttime= 0.5*(sfrtime[i] + sfrtime[i-1])   ; gyr/h
	Time_Lum[i-1]= Ttime

	;  bh
	; ====
	;uses black hole details - slow
        ;bhbololum= fload_blackhole_lum(frun,Ttime,/bolometric)
	;uses black hole txt file - fast
	if i ge Ni_bh then bhbololum= boloL_sun[Ni_bh-1] else bhbololum= 0.5*(boloL_sun[i]+boloL_sun[i-1])

        if bhbololum[0] gt 0.0 then BH_Bolo_Lum[i-1]= alog10(total(bhbololum))

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

	stellarlum= get_stellar_lum(iN, itime, iMass, iAge, iZ)
	S_Bolo_Lum[i-1]= alog10(stellarlum)
endfor



;-------------------------------------------



openw, 1, frun+'/luminosity.txt', ERROR=err

printf, 1, "#   luminosity.txt "
printf, 1, "#   Bolometric luminosity of stars and blackhole"
printf, 1, "#        "
printf, 1, "#  time    stellar        BH  "
printf, 1, "#(Gyr/h)  luminosity luminosity "
for i=0,Ni-1 do begin
        printf, 1, FORMAT= '(F8.4,"   ",2(F9.4,"  "))', $
                Time_Lum[i], S_Bolo_Lum[i], BH_Bolo_Lum[i]
endfor
close, 1



end






;=====================================================================================




