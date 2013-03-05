pro desika_angles, junk



Band= 'B'
;frun= '/home/tcox/yuexing'
;frun= '/home/tcox/yuexing2'
;frun= '/home/tcox/yuexing3'
frun= '/home/tcox/yuexing4'
;frun= '/home/tcox/yuexing5'
;frun= '/raid4/tcox/z6/a6'
angledir= '/los'
snapnum= 4
;snapnum= 8
use_calc_center= 1


;  Center
; --------------------------------
center=[0.0,0.0,0.0]




; --------------------------------
;  load desika angles
; --------------------------------

; full version
;desika_file= 'desika_thousandlos.snap8.check'
;palantiri /home/tcox> wc desika_thousandlos.snap8.check 
; 1000  2000 26000 desika_thousandlos.snap8.check
;dlines= 1000
; testing
;desika_file= 'desika_thousandlos.snap8.test' & tag= '.testheader'
;dlines= 5
; parse it up
;if junk eq  1 then desika_file= '/home/tcox/yuexing/los/loslist/desika_los.1.txt'  & tag= '.1'

junklbl= strcompress(string(junk),/remove_all)

desika_file= '/home/tcox/yuexing/los/loslist/desika_los.'+junklbl+'.txt'
print, "desika_file=",desika_file
tag= '.'+junklbl
print, "tag=",tag

spawn, "wc "+desika_file,result
dlines= long(result[0])
;dlines= 100


ddata= fltarr(2,dlines)

openr, 1, desika_file
readf, 1, ddata
close, 1

desika_theta= transpose(ddata[0,*])
desika_phi= transpose(ddata[1,*])




Stellar_Luminosity= fltarr(dlines)
BH_Luminosity= fltarr(dlines)





; --------------------------------
;  Load variables for smoothing
; --------------------------------
    if not keyword_set(loadedsnap) then begin
      ;if (fload_snapshot(frun, snapnum)) then begin
      if (fload_snapshot_bh(frun, snapnum)) then begin
        print, "PROBLEM: opening file"
        return
      endif
    endif

    if keyword_set(use_calc_center) then center=fload_center_alreadycomp(1)

    ndisk= fload_npart(2)				; disk particles
    nbulge= fload_npart(3)				; bulge particles
    nstars= fload_npart(4)
    npart= long(ndisk) + long(nbulge) + long(nstars)
    print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars

    x= fload_allstars_xyz('x',center=[0,0,0])
    y= fload_allstars_xyz('y',center=[0,0,0])
    z= fload_allstars_xyz('z',center=[0,0,0])
    ;m= fload_allstars_mass(1)

    ; if we ever want to seperate out just one galaxy
    ;if keyword_set(startid) and keyword_set(numpart) then begin
	;x= fload_1gal_allstars_xyz('x',startid,numpart,center=[0,0,0])
	;y= fload_1gal_allstars_xyz('y',startid,numpart,center=[0,0,0])
	;z= fload_1gal_allstars_xyz('z',startid,numpart,center=[0,0,0])
	;m= fload_1gal_allstars_xyz(startid,numpart)
    ;endif

    ;hsml= 0.0*x + 0.1    ; assume a standard smoothing size
    ;hsml= 0.5*fload_allstars_hsml(1)

    TTime= float(fload_time(1))

    N= npart
    m= 1.0e+10*fload_allstars_mass(1)
    age=fload_allstars_age(1)
    age=float(TTime-age)
    ;zmets=fload_allstars_z(1)
    zmets=fload_allstars_z(1) / 0.02    ; brant's new code requires it in solar


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


    Bolometric_Luminosity= transpose(Lum_Bol)

    print, "Bolometric Lum= ", total(Lum_Bol)
    print, "Total U-band Lum= ", total(Lum_U)
    print, "Total B-band Lum= ", total(Lum_B)
    print, "Total V-band Lum= ", total(Lum_V)
    print, "Total K-band Lum= ", total(Lum_K)
    print, "Total sdss_u-band Lum= ", total(Lum_sdss_u)






; ------------------------------------
;  Specify Band-related Information
; ------------------------------------


Bolometric_Solar_Luminosity= 3.826d+33    ; in ergs/sec

;
; lambda and width are in Angstroms
; where'd I get solar lums?
;
; -----------------------------------

; Number of Bands
iNBands= 13

; only used in this next 30 lines of code
;
iBand= strarr(iNBands)
iBand_Lambda= fltarr(iNBands)
iBand_width= fltarr(iNBands)
iBand_Solar_Luminosity= dblarr(iNBands)

iBand[0]= "U"  & iBand_Lambda[0] = 3650   &  iBand_width[0]= 700   & iBand_Solar_Luminosity[0]= 2.0d+32       ; U
iBand[1]= "B"  & iBand_Lambda[1] = 4400   &  iBand_width[1]= 980   & iBand_Solar_Luminosity[1]= 5.2d+32       ; B
iBand[2]= "V"  & iBand_Lambda[2] = 5500   &  iBand_width[2]= 980   & iBand_Solar_Luminosity[2]= 5.2d+32       ; V
iBand[3]= "R"  & iBand_Lambda[3] = 6500   &  iBand_width[3]= 1180  & iBand_Solar_Luminosity[3]= 7.7d+32       ; R
iBand[4]= "I"  & iBand_Lambda[4] = 8000   &  iBand_width[4]= 1400  & iBand_Solar_Luminosity[4]= 5.2d+32       ; I
iBand[5]= "J"  & iBand_Lambda[5] = 12150  &  iBand_width[5]= 2600  & iBand_Solar_Luminosity[5]= 2.8d+32       ; J
iBand[6]= "H"  & iBand_Lambda[6] = 16540  &  iBand_width[6]= 2900  & iBand_Solar_Luminosity[6]= 1.8d+32       ; H
iBand[7]= "K"  & iBand_Lambda[7] = 21790  &  iBand_width[7]= 4100  & iBand_Solar_Luminosity[7]= 0.8d+32       ; K
iBand[8]= "u"  & iBand_Lambda[8] = 3543   &  iBand_width[8]= 650   & iBand_Solar_Luminosity[8]= 2.0d+32       ; u
iBand[9]= "g"  & iBand_Lambda[9] = 4770   &  iBand_width[9]= 1480  & iBand_Solar_Luminosity[9]= 5.2d+32       ; g
iBand[10]= "r" & iBand_Lambda[10] = 6231  &  iBand_width[10]= 1385 & iBand_Solar_Luminosity[10]= 7.5d+32      ; r
iBand[11]= "i" & iBand_Lambda[11] = 7625  &  iBand_width[11]= 1565 & iBand_Solar_Luminosity[11]= 5.2d+32      ; i
iBand[12]= "z" & iBand_Lambda[12] = 9134  &  iBand_width[12]= 1130 & iBand_Solar_Luminosity[12]= 4.0d+32      ; z


idx= where(Band eq iBand)
if (idx(0) eq -1) or (n_elements(idx) gt 1) then begin
	print, 'what up with idx?'
	help, idx
	print, ' '
	print, ' Something is wrong with selecting the band,'
	print, ' now exiting.'
	print, ' '
	return
endif else begin
	Band_Lambda = iBand_Lambda[idx(0)]
	Band_width= iBand_width[idx(0)]
	Band_Solar_Luminosity= iBand_Solar_Luminosity[idx(0)]
endelse

case idx of
   0: ProjectThisLum= Lum_U
   1: ProjectThisLum= Lum_B
   2: ProjectThisLum= Lum_V
   3: ProjectThisLum= Lum_R
   4: ProjectThisLum= Lum_I
   5: ProjectThisLum= Lum_J
   6: ProjectThisLum= Lum_H
   7: ProjectThisLum= Lum_K
   8: ProjectThisLum= Lum_sdss_u
   9: ProjectThisLum= Lum_sdss_g
   10: ProjectThisLum= Lum_sdss_r
   12: ProjectThisLum= Lum_sdss_i
   12: ProjectThisLum= Lum_sdss_z
endcase


print, "----------------------------------------"
print, "using Band= ", Band
print, "using Band_Lambda= ", Band_Lambda
print, "using Band_width= ", Band_width
print, "using Solar_Lum= ", Bolometric_Solar_Luminosity
print, "using Band Luminosity= ", Band_Solar_Luminosity
print, "----------------------------------------"


if not keyword_set(msg) then msg=Band+'-band'



; ---------------------------------
;  Add BH Luminosity
; ---------------------------------
include_BH_lum= 1
;include_BH_lum= 0
if include_BH_lum eq 1 and fload_npart(5) gt 0 then begin

	; spectrum is just the current band
	; -----------------------------------
	min_Lambda= min(iBand_Lambda - iBand_width) - 100.0
	max_Lambda= max(iBand_Lambda - iBand_width) + 100.0
	N_spectrum= 20000L   ; doesn't work well with low N
	c_light= 3.0e+8
	;log_nu_min= alog10(c_light/(Band_Lambda+0.5*Band_width)*1.0e+10)    ; 1e+10 converts to m from Angstroms
	;log_nu_max= alog10(c_light/(Band_Lambda-0.5*Band_width)*1.0e+10)
	log_nu_min= alog10(c_light/(max_Lambda)*1.0e+10)    ; 1e+10 converts to m from Angstroms
        log_nu_max= alog10(c_light/(min_Lambda)*1.0e+10)
	print, "max_lambda= ", max_lambda
        print, "min_lambda= ", min_lambda

	; need this to (easily) concatenate the BH luminosities
	ProjectThisLum= transpose(ProjectThisLum)

	; get BH lum
	; -------------
	Nbh= fload_npart(5)
	bhids= fload_blackhole_id(1)
	bhbololum= fload_blackhole_lum(frun,Ttime,/bolometric)

	for i=1,fload_npart(5) do begin
		print, "Adding BH ",i, "  id= ",bhids[i-1]
		;BH_bolometric_lum= 12.0      ; bolometric luminosity in log
		BH_bolometric_lum= alog10(bhbololum[i-1])

		; updated by Phil (June 2006)
		; ---------------------------
		;agn= fload_marconi_agn_spectrum(BH_bolometric_lum, 0.0, N_spectrum, log_nu_min, log_nu_max)
		;nu= agn[0:N_spectrum-1]
		;lambda_nu= c_light / (10.0^nu) * 1.0e+10
		;l_nu= agn[N_spectrum:2*N_spectrum-1]
		;idx= where((lambda_nu gt (Band_Lambda-Band_width)) and $
                ;                        (lambda_nu lt (Band_Lambda+Band_width)))
		;BH_lum_inband = mean(10.0^(l_nu(idx)) * 10.0^(nu(idx)))
		; ---------------------------
                log_nu = alog10(c_light/Band_Lambda*1.0e+10)
                BH_lum_inband = attenuated_spectrum(log_nu,BH_bolometric_lum,18.)
                        BH_lum_inband = 10^(BH_lum_inband(0))


		print, "BH Bolometric Luminosity = ",10.0^(BH_bolometric_lum)
		print, "BH "+Band+"-band Luminosity (in bolometric solar units) = ",BH_lum_inband
		print, "Bolometric Correction = ",10.0^(BH_bolometric_lum)/BH_lum_inband
		BH_lum_inband= BH_lum_inband * Bolometric_Solar_Luminosity / Band_Solar_Luminosity
		print, "BH "+Band+"-band Luminosity (in "+Band+"-band solar units) = ",BH_lum_inband

		; Add BH to list
		; ---------------
		ProjectThisLum= [BH_lum_inband, ProjectThisLum]
		Bolometric_Luminosity= [10^BH_bolometric_lum, Bolometric_Luminosity]
	endfor

	; Add BH to other lists
	; -----------------------
	x= [fload_blackhole_xyz('x',center=[0,0,0]), x]
	y= [fload_blackhole_xyz('y',center=[0,0,0]), y]
	z= [fload_blackhole_xyz('z',center=[0,0,0]), z]

	;for i=1,fload_npart(5) do hsml= [0.1, hsml]

	set_bh_to_zero= 0

endif else begin
	set_bh_to_zero= 1
endelse



; ---------------------------------
;  Attenuate the Light
; ---------------------------------


        Ngas= long(fload_npart(0))
        Nstars= long(fload_npart(2)+fload_npart(3)+fload_npart(4))
        Nbh= long(fload_npart(5))
	if set_bh_to_zero eq 1 then Nbh= 0L

        ; fields to pass
        Coord = fltarr(10,Ngas+Nstars+Nbh)

        ; gas
        Coord(0,0:Ngas-1) = fload_gas_xyz('x', center=[0,0,0])
        Coord(1,0:Ngas-1) = fload_gas_xyz('y', center=[0,0,0])
        Coord(2,0:Ngas-1) = fload_gas_xyz('z', center=[0,0,0])

        Coord(3,0:Ngas-1) = fload_gas_u(1)
        Coord(4,0:Ngas-1) = fload_gas_rho(1)
        Coord(5,0:Ngas-1) = fload_gas_hsml(1)
        Coord(6,0:Ngas-1) = fload_gas_numh(1)
        Coord(7,0:Ngas-1) = fload_gas_nume(1)
        Coord(8,0:Ngas-1) = fload_gas_metallicity(1)
        Coord(9,0:Ngas-1) = fload_gas_mass(1)

        ; stars
        Coord(0,Ngas:Nstars+Ngas-1) = fload_allstars_xyz('x', center=[0,0,0])
        Coord(1,Ngas:Nstars+Ngas-1) = fload_allstars_xyz('y', center=[0,0,0])
        Coord(2,Ngas:Nstars+Ngas-1) = fload_allstars_xyz('z', center=[0,0,0])

        ; black holes
	if Nbh gt 0 then begin
        	Coord(0,Ngas+Nstars:Nstars+Ngas+Nbh-1) = fload_blackhole_xyz('x', center=[0,0,0])
        	Coord(1,Ngas+Nstars:Nstars+Ngas+Nbh-1) = fload_blackhole_xyz('y', center=[0,0,0])
        	Coord(2,Ngas+Nstars:Nstars+Ngas+Nbh-1) = fload_blackhole_xyz('z', center=[0,0,0])
	endif





; ------------------------------------------------
; ------------------------------------------------
;
;    Now actually calculate the luminosities
;
; ------------------------------------------------
; ------------------------------------------------
for id= 0, dlines-1 do begin

        los_NH= fltarr(Nstars+Nbh)
        los_Z= fltarr(Nstars+Nbh)

        print, "=================================================="
	print, id,id,id,id,id,id
        print, "=================================================="

	; need it in radians
	theta= (!PI / 180.0) * desika_theta[id]
	phi= (!PI / 180.0) * desika_phi[id]


        print, "PASSING: "
        print, "N_gas= ", Ngas
        print, "N_stars= ", Nstars
        print, "N_bh= ", Nbh
        print, "theta= ", theta
        print, "phi= ", phi
        help, Coord

	; old version
	;S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/LOSColumn_RhoZ/getnh.so', $
	; Phil's newest version
	S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/LOS_column_code/getnh.so', $
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
	;   and zero metallicity (make it small instead)
	idx= where(los_NH lt 1.0e+10)
	if idx(0) ne -1 then los_NH(idx)= 1.0e+10
	idx= where(los_Z le 0.0)
	if idx(0) ne -1 then los_Z(idx)= 1.0e-5


	;  OK, now we've got the N_h values,
	; what is the attenuated luminosity
	PreLums= ProjectThisLum
        band_nu = 3.0d8/(band_lambda * 1.0d-10)
        sigma   = cross_section(band_nu)
        AttLums = PreLums * EXP(-sigma*los_NH)

	print, "total light= ", total(ProjectThisLum), total(PreLums)
	print, "attenuated light= ", total(AttLums)
	print, "    abs. f= ", total(AttLums)/total(ProjectThisLum)

	print, " "
	print, " Stars "
	print, " ------"
	print, "att. light = ", total(AttLums[Nbh:Nbh+Nstars-1])
	print, "     abs. f= ", total(AttLums[Nbh:Nbh+Nstars-1])/total(ProjectThisLum[Nbh:Nbh+Nstars-1])

	Stellar_Luminosity[id]= total(AttLums[Nbh:Nbh+Nstars-1])

	if Nbh gt 0 then begin
		print, " "
		print, " BH  "
		print, " ------"
		print, "att. light = ", total(AttLums[0:Nbh-1])
		print, "     abs. f= ", total(AttLums[0:Nbh-1])/total(ProjectThisLum[0:Nbh-1])
		print, "     los_NH= ", los_NH[0:Nbh-1]
		print, "     los_Z =  ", los_Z[0:Nbh-1]
		print, " "
	endif else begin
		print, " NO Black Hole !!!"
	endelse

	BH_Luminosity[id]= total(AttLums[0:Nbh-1])


endfor


; ---------------------
blum_tot= strcompress(string(total(ProjectThisLum)),/remove_all)
blum_bh= strcompress(string(total(ProjectThisLum[0:Nbh-1])),/remove_all)
blum_sf= strcompress(string(total(ProjectThisLum[Nbh:Nbh+Nstars-1])),/remove_all)

alum_tot= strcompress(string(total(AttLums)),/remove_all)
alum_bh= strcompress(string(total(AttLums[0:Nbh-1])),/remove_all)
alum_sf= strcompress(string(total(AttLums[Nbh:Nbh+Nstars-1])),/remove_all)

bololum_tot= strcompress(string(total(Bolometric_Luminosity)),/remove_all)
bololum_bh= strcompress(string(total(Bolometric_Luminosity[0:Nbh-1])),/remove_all)
bololum_sf= strcompress(string(total(Bolometric_Luminosity[Nbh:Nbh+Nstars-1])),/remove_all)


; ---------------------
;  Done - print it up
; ---------------------
snapstr= strcompress(string(snapnum),/remove_all)
filename= frun+angledir+'/anglelums_'+snapstr+'.txt'+tag
openw, 1, filename, ERROR=err
        
printf, 1, "#  "+filename
printf, 1, "#        B-Band Luminosity=  "+blum_tot+"  (bh/sf= "+blum_bh+"/"+blum_sf+")"
printf, 1, "#  (att) B-Band Luminosity=  "+alum_tot+"  (bh/sf= "+alum_bh+"/"+alum_sf+")"
printf, 1, "#    Bolometric Luminosity=  "+bololum_tot+"  (bh/sf= "+bololum_bh+"/"+bololum_sf+")"
printf, 1, "#  theta      phi       Stellar        BH  "
printf, 1, "#  (deg)     (deg)     (L_solar)   (L_solar)"
for i=0,dlines-1 do begin
        printf, 1, FORMAT= '(F8.3,"  ",F8.3,"  ", 2(E10.3,"  "))', $
                desika_theta[i], desika_phi[i], Stellar_Luminosity[i], BH_Luminosity[i]
endfor
close, 1






end








;=============================================================================








pro da2, junk



Band= 'B'
;frun= '/home/tcox/yuexing'
;frun= '/home/tcox/yuexing2'
;frun= '/home/tcox/yuexing3'
;frun= '/home/tcox/yuexing4'
frun= '/home/tcox/yuexing5'
;frun= '/raid4/tcox/z6/a6'
angledir= '/los'
;snapnum= 4
snapnum= 5
;snapnum= 6
;snapnum= 10
use_calc_center= 1


;  Center
; --------------------------------
center=[0.0,0.0,0.0]




; --------------------------------
;  load desika angles
; --------------------------------

; full version
;desika_file= 'desika_thousandlos.snap8.check'
;palantiri /home/tcox> wc desika_thousandlos.snap8.check 
; 1000  2000 26000 desika_thousandlos.snap8.check
;dlines= 1000
; testing
;desika_file= 'desika_thousandlos.snap8.test' & tag= '.testheader'
;dlines= 5
; parse it up
;if junk eq  1 then desika_file= '/home/tcox/yuexing/los/loslist/desika_los.1.txt'  & tag= '.1'

junklbl= strcompress(string(junk),/remove_all)

desika_file= '/home/tcox/yuexing/los/loslist/desika_los.'+junklbl+'.txt'
print, "desika_file=",desika_file
tag= '.'+junklbl
print, "tag=",tag

spawn, "wc "+desika_file,result
dlines= long(result[0])
;dlines= 100


ddata= fltarr(2,dlines)

openr, 1, desika_file
readf, 1, ddata
close, 1

desika_theta= transpose(ddata[0,*])
desika_phi= transpose(ddata[1,*])




Stellar_Luminosity= fltarr(dlines)
BH_Luminosity= fltarr(dlines)





; --------------------------------
;  Load variables for smoothing
; --------------------------------
    if not keyword_set(loadedsnap) then begin
      ;if (fload_snapshot(frun, snapnum)) then begin
      if (fload_snapshot_bh(frun, snapnum)) then begin
        print, "PROBLEM: opening file"
        return
      endif
    endif

    if keyword_set(use_calc_center) then center=fload_center_alreadycomp(1)

    ndisk= fload_npart(2)				; disk particles
    nbulge= fload_npart(3)				; bulge particles
    nstars= fload_npart(4)
    npart= long(ndisk) + long(nbulge) + long(nstars)
    print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars

    x= fload_allstars_xyz('x',center=[0,0,0])
    y= fload_allstars_xyz('y',center=[0,0,0])
    z= fload_allstars_xyz('z',center=[0,0,0])
    ;m= fload_allstars_mass(1)

    ; if we ever want to seperate out just one galaxy
    ;if keyword_set(startid) and keyword_set(numpart) then begin
	;x= fload_1gal_allstars_xyz('x',startid,numpart,center=[0,0,0])
	;y= fload_1gal_allstars_xyz('y',startid,numpart,center=[0,0,0])
	;z= fload_1gal_allstars_xyz('z',startid,numpart,center=[0,0,0])
	;m= fload_1gal_allstars_xyz(startid,numpart)
    ;endif

    ;hsml= 0.0*x + 0.1    ; assume a standard smoothing size
    ;hsml= 0.5*fload_allstars_hsml(1)

    TTime= float(fload_time(1))

    N= npart
    m= 1.0e+10*fload_allstars_mass(1)
    age=fload_allstars_age(1)
    age=float(TTime-age)
    ;zmets=fload_allstars_z(1)
    zmets=fload_allstars_z(1) / 0.02    ; brant's new code requires it in solar


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


    Bolometric_Luminosity= transpose(Lum_Bol)

    print, "Bolometric Lum= ", total(Lum_Bol)
    print, "Total U-band Lum= ", total(Lum_U)
    print, "Total B-band Lum= ", total(Lum_B)
    print, "Total V-band Lum= ", total(Lum_V)
    print, "Total K-band Lum= ", total(Lum_K)
    print, "Total sdss_u-band Lum= ", total(Lum_sdss_u)






; ------------------------------------
;  Specify Band-related Information
; ------------------------------------


Bolometric_Solar_Luminosity= 3.826d+33    ; in ergs/sec

;
; lambda and width are in Angstroms
; where'd I get solar lums?
;
; -----------------------------------

; Number of Bands
iNBands= 13

; only used in this next 30 lines of code
;
iBand= strarr(iNBands)
iBand_Lambda= fltarr(iNBands)
iBand_width= fltarr(iNBands)
iBand_Solar_Luminosity= dblarr(iNBands)

iBand[0]= "U"  & iBand_Lambda[0] = 3650   &  iBand_width[0]= 700   & iBand_Solar_Luminosity[0]= 2.0d+32       ; U
iBand[1]= "B"  & iBand_Lambda[1] = 4400   &  iBand_width[1]= 980   & iBand_Solar_Luminosity[1]= 5.2d+32       ; B
iBand[2]= "V"  & iBand_Lambda[2] = 5500   &  iBand_width[2]= 980   & iBand_Solar_Luminosity[2]= 5.2d+32       ; V
iBand[3]= "R"  & iBand_Lambda[3] = 6500   &  iBand_width[3]= 1180  & iBand_Solar_Luminosity[3]= 7.7d+32       ; R
iBand[4]= "I"  & iBand_Lambda[4] = 8000   &  iBand_width[4]= 1400  & iBand_Solar_Luminosity[4]= 5.2d+32       ; I
iBand[5]= "J"  & iBand_Lambda[5] = 12150  &  iBand_width[5]= 2600  & iBand_Solar_Luminosity[5]= 2.8d+32       ; J
iBand[6]= "H"  & iBand_Lambda[6] = 16540  &  iBand_width[6]= 2900  & iBand_Solar_Luminosity[6]= 1.8d+32       ; H
iBand[7]= "K"  & iBand_Lambda[7] = 21790  &  iBand_width[7]= 4100  & iBand_Solar_Luminosity[7]= 0.8d+32       ; K
iBand[8]= "u"  & iBand_Lambda[8] = 3543   &  iBand_width[8]= 650   & iBand_Solar_Luminosity[8]= 2.0d+32       ; u
iBand[9]= "g"  & iBand_Lambda[9] = 4770   &  iBand_width[9]= 1480  & iBand_Solar_Luminosity[9]= 5.2d+32       ; g
iBand[10]= "r" & iBand_Lambda[10] = 6231  &  iBand_width[10]= 1385 & iBand_Solar_Luminosity[10]= 7.5d+32      ; r
iBand[11]= "i" & iBand_Lambda[11] = 7625  &  iBand_width[11]= 1565 & iBand_Solar_Luminosity[11]= 5.2d+32      ; i
iBand[12]= "z" & iBand_Lambda[12] = 9134  &  iBand_width[12]= 1130 & iBand_Solar_Luminosity[12]= 4.0d+32      ; z


idx= where(Band eq iBand)
if (idx(0) eq -1) or (n_elements(idx) gt 1) then begin
	print, 'what up with idx?'
	help, idx
	print, ' '
	print, ' Something is wrong with selecting the band,'
	print, ' now exiting.'
	print, ' '
	return
endif else begin
	Band_Lambda = iBand_Lambda[idx(0)]
	Band_width= iBand_width[idx(0)]
	Band_Solar_Luminosity= iBand_Solar_Luminosity[idx(0)]
endelse

case idx of
   0: ProjectThisLum= Lum_U
   1: ProjectThisLum= Lum_B
   2: ProjectThisLum= Lum_V
   3: ProjectThisLum= Lum_R
   4: ProjectThisLum= Lum_I
   5: ProjectThisLum= Lum_J
   6: ProjectThisLum= Lum_H
   7: ProjectThisLum= Lum_K
   8: ProjectThisLum= Lum_sdss_u
   9: ProjectThisLum= Lum_sdss_g
   10: ProjectThisLum= Lum_sdss_r
   12: ProjectThisLum= Lum_sdss_i
   12: ProjectThisLum= Lum_sdss_z
endcase


print, "----------------------------------------"
print, "using Band= ", Band
print, "using Band_Lambda= ", Band_Lambda
print, "using Band_width= ", Band_width
print, "using Solar_Lum= ", Bolometric_Solar_Luminosity
print, "using Band Luminosity= ", Band_Solar_Luminosity
print, "----------------------------------------"


if not keyword_set(msg) then msg=Band+'-band'



; ---------------------------------
;  Add BH Luminosity
; ---------------------------------
include_BH_lum= 1
;include_BH_lum= 0
if include_BH_lum eq 1 and fload_npart(5) gt 0 then begin

	; spectrum is just the current band
	; -----------------------------------
	min_Lambda= min(iBand_Lambda - iBand_width) - 100.0
	max_Lambda= max(iBand_Lambda - iBand_width) + 100.0
	N_spectrum= 20000L   ; doesn't work well with low N
	c_light= 3.0e+8
	;log_nu_min= alog10(c_light/(Band_Lambda+0.5*Band_width)*1.0e+10)    ; 1e+10 converts to m from Angstroms
	;log_nu_max= alog10(c_light/(Band_Lambda-0.5*Band_width)*1.0e+10)
	log_nu_min= alog10(c_light/(max_Lambda)*1.0e+10)    ; 1e+10 converts to m from Angstroms
        log_nu_max= alog10(c_light/(min_Lambda)*1.0e+10)
	print, "max_lambda= ", max_lambda
        print, "min_lambda= ", min_lambda

	; need this to (easily) concatenate the BH luminosities
	ProjectThisLum= transpose(ProjectThisLum)

	; get BH lum
	; -------------
	Nbh= fload_npart(5)
	bhids= fload_blackhole_id(1)
	bhbololum= fload_blackhole_lum(frun,Ttime,/bolometric)

	for i=1,fload_npart(5) do begin
		print, "Adding BH ",i, "  id= ",bhids[i-1]
		;BH_bolometric_lum= 12.0      ; bolometric luminosity in log
		BH_bolometric_lum= alog10(bhbololum[i-1])

		; updated by Phil (June 2006)
		; ---------------------------
		;agn= fload_marconi_agn_spectrum(BH_bolometric_lum, 0.0, N_spectrum, log_nu_min, log_nu_max)
		;nu= agn[0:N_spectrum-1]
		;lambda_nu= c_light / (10.0^nu) * 1.0e+10
		;l_nu= agn[N_spectrum:2*N_spectrum-1]
		;idx= where((lambda_nu gt (Band_Lambda-Band_width)) and $
                ;                        (lambda_nu lt (Band_Lambda+Band_width)))
		;BH_lum_inband = mean(10.0^(l_nu(idx)) * 10.0^(nu(idx)))
		; ---------------------------
                log_nu = alog10(c_light/Band_Lambda*1.0e+10)
                BH_lum_inband = attenuated_spectrum(log_nu,BH_bolometric_lum,18.)
                        BH_lum_inband = 10^(BH_lum_inband(0))


		print, "BH Bolometric Luminosity = ",10.0^(BH_bolometric_lum)
		print, "BH "+Band+"-band Luminosity (in bolometric solar units) = ",BH_lum_inband
		print, "Bolometric Correction = ",10.0^(BH_bolometric_lum)/BH_lum_inband
		BH_lum_inband= BH_lum_inband * Bolometric_Solar_Luminosity / Band_Solar_Luminosity
		print, "BH "+Band+"-band Luminosity (in "+Band+"-band solar units) = ",BH_lum_inband

		; Add BH to list
		; ---------------
		ProjectThisLum= [BH_lum_inband, ProjectThisLum]
		Bolometric_Luminosity= [10^BH_bolometric_lum, Bolometric_Luminosity]
	endfor

	; Add BH to other lists
	; -----------------------
	x= [fload_blackhole_xyz('x',center=[0,0,0]), x]
	y= [fload_blackhole_xyz('y',center=[0,0,0]), y]
	z= [fload_blackhole_xyz('z',center=[0,0,0]), z]

	;for i=1,fload_npart(5) do hsml= [0.1, hsml]

	set_bh_to_zero= 0

endif else begin
	set_bh_to_zero= 1
endelse



; ---------------------------------
;  Attenuate the Light
; ---------------------------------


        Ngas= long(fload_npart(0))
        Nstars= long(fload_npart(2)+fload_npart(3)+fload_npart(4))
        Nbh= long(fload_npart(5))
	if set_bh_to_zero eq 1 then Nbh= 0L

        ; fields to pass
        Coord = fltarr(10,Ngas+Nstars+Nbh)

        ; gas
        Coord(0,0:Ngas-1) = fload_gas_xyz('x', center=[0,0,0])
        Coord(1,0:Ngas-1) = fload_gas_xyz('y', center=[0,0,0])
        Coord(2,0:Ngas-1) = fload_gas_xyz('z', center=[0,0,0])

        Coord(3,0:Ngas-1) = fload_gas_u(1)
        Coord(4,0:Ngas-1) = fload_gas_rho(1)
        Coord(5,0:Ngas-1) = fload_gas_hsml(1)
        Coord(6,0:Ngas-1) = fload_gas_numh(1)
        Coord(7,0:Ngas-1) = fload_gas_nume(1)
        Coord(8,0:Ngas-1) = fload_gas_metallicity(1)
        Coord(9,0:Ngas-1) = fload_gas_mass(1)

        ; stars
        Coord(0,Ngas:Nstars+Ngas-1) = fload_allstars_xyz('x', center=[0,0,0])
        Coord(1,Ngas:Nstars+Ngas-1) = fload_allstars_xyz('y', center=[0,0,0])
        Coord(2,Ngas:Nstars+Ngas-1) = fload_allstars_xyz('z', center=[0,0,0])

        ; black holes
	if Nbh gt 0 then begin
        	Coord(0,Ngas+Nstars:Nstars+Ngas+Nbh-1) = fload_blackhole_xyz('x', center=[0,0,0])
        	Coord(1,Ngas+Nstars:Nstars+Ngas+Nbh-1) = fload_blackhole_xyz('y', center=[0,0,0])
        	Coord(2,Ngas+Nstars:Nstars+Ngas+Nbh-1) = fload_blackhole_xyz('z', center=[0,0,0])
	endif





; ------------------------------------------------
; ------------------------------------------------
;
;    Now actually calculate the luminosities
;
; ------------------------------------------------
; ------------------------------------------------


startangle= 0
endangle= dlines-1


; ---------------------
; ---------------------
snapstr= strcompress(string(snapnum),/remove_all)
filename= frun+angledir+'/anglelums_'+snapstr+'_data.txt'+tag


spawn, "wc "+filename,result
filename_lines= long(result[0])
; there is supposed to be 106 (6 header + 100 angles) lines

if filename_lines eq 106 then begin
	print, " "
	print, " filename= "+filename
	print, " already has 106 lines of data "
	print, " "
	return
endif

if filename_lines gt 0 then begin
;	openu, 1, filename, ERROR=err, /append
	startangle= filename_lines
endif ; else begin
;	openw, 1, filename, ERROR=err

        ; actually, we don't have these quantities yet
	;printf, 1, "#  "+filename
	;printf, 1, "#        B-Band Luminosity=  "+blum_tot+"  (bh/sf= "+blum_bh+"/"+blum_sf+")"
	;printf, 1, "#  (att) B-Band Luminosity=  "+alum_tot+"  (bh/sf= "+alum_bh+"/"+alum_sf+")"
	;printf, 1, "#    Bolometric Luminosity=  "+bololum_tot+"  (bh/sf= "+bololum_bh+"/"+bololum_sf+")"
	;printf, 1, "#  theta      phi       Stellar        BH  "
	;printf, 1, "#  (deg)     (deg)     (L_solar)   (L_solar)"
;endelse




; ---------------------
; ---------------------
for id= startangle, endangle do begin

        los_NH= fltarr(Nstars+Nbh)
        los_Z= fltarr(Nstars+Nbh)

        print, "=================================================="
	print, id,id,id,id,id,id
        print, "=================================================="

	; need it in radians
	theta= (!PI / 180.0) * desika_theta[id]
	phi= (!PI / 180.0) * desika_phi[id]


        print, "PASSING: "
        print, "N_gas= ", Ngas
        print, "N_stars= ", Nstars
        print, "N_bh= ", Nbh
        print, "theta= ", theta
        print, "phi= ", phi
        help, Coord

	; old version
	;S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/LOSColumn_RhoZ/getnh.so', $
	; Phil's newest version
	S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/LOS_column_code/getnh.so', $
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
	;   and zero metallicity (make it small instead)
	idx= where(los_NH lt 1.0e+10)
	if idx(0) ne -1 then los_NH(idx)= 1.0e+10
	idx= where(los_Z le 0.0)
	if idx(0) ne -1 then los_Z(idx)= 1.0e-5


	;  OK, now we've got the N_h values,
	; what is the attenuated luminosity
	PreLums= ProjectThisLum
        band_nu = 3.0d8/(band_lambda * 1.0d-10)
        sigma   = cross_section(band_nu)
        AttLums = PreLums * EXP(-sigma*los_NH)

	print, "sigma= ", sigma
	print, "los_NH= ", min(los_NH), max(los_NH), mean(los_NH)

	print, "total light= ", total(ProjectThisLum), total(PreLums)
	print, "attenuated light= ", total(AttLums)
	print, "    abs. f= ", total(AttLums)/total(ProjectThisLum)

	print, " "
	print, " Stars "
	print, " ------"
	print, "att. light = ", total(AttLums[Nbh:Nbh+Nstars-1])
	print, "     abs. f= ", total(AttLums[Nbh:Nbh+Nstars-1])/total(ProjectThisLum[Nbh:Nbh+Nstars-1])

	Stellar_Luminosity[id]= total(AttLums[Nbh:Nbh+Nstars-1])

	if Nbh gt 0 then begin
		print, " "
		print, " BH  "
		print, " ------"
		print, "att. light = ", total(AttLums[0:Nbh-1])
		print, "     abs. f= ", total(AttLums[0:Nbh-1])/total(ProjectThisLum[0:Nbh-1])
		print, "     los_NH= ", los_NH[0:Nbh-1]
		print, "     los_Z =  ", los_Z[0:Nbh-1]
		print, " "
	endif else begin
		print, " NO Black Hole !!!"
	endelse

	BH_Luminosity[id]= total(AttLums[0:Nbh-1])



	; ---------------------
	; ---------------------
	openu, 1, filename, ERROR=err, /append
	printf, 1, FORMAT= '(F8.3,"  ",F8.3,"  ", 2(E10.3,"  "))', $
                desika_theta[id], desika_phi[id], Stellar_Luminosity[id], BH_Luminosity[id]
	close, 1


endfor





; ---------------------
; ---------------------
blum_tot= strcompress(string(total(ProjectThisLum)),/remove_all)
blum_bh= strcompress(string(total(ProjectThisLum[0:Nbh-1])),/remove_all)
blum_sf= strcompress(string(total(ProjectThisLum[Nbh:Nbh+Nstars-1])),/remove_all)

alum_tot= strcompress(string(total(AttLums)),/remove_all)
alum_bh= strcompress(string(total(AttLums[0:Nbh-1])),/remove_all)
alum_sf= strcompress(string(total(AttLums[Nbh:Nbh+Nstars-1])),/remove_all)

bololum_tot= strcompress(string(total(Bolometric_Luminosity)),/remove_all)
bololum_bh= strcompress(string(total(Bolometric_Luminosity[0:Nbh-1])),/remove_all)
bololum_sf= strcompress(string(total(Bolometric_Luminosity[Nbh:Nbh+Nstars-1])),/remove_all)




; ---------------------
; ---------------------
tempfilename= frun+angledir+'/anglelums_'+snapstr+'_header.txt'+tag
finalfilename= frun+angledir+'/anglelums_'+snapstr+'.txt'+tag

openw, 1, tempfilename, ERROR=err
printf, 1, "#  "+finalfilename
printf, 1, "#        B-Band Luminosity=  "+blum_tot+"  (bh/sf= "+blum_bh+"/"+blum_sf+")"
printf, 1, "#  (att) B-Band Luminosity=  "+alum_tot+"  (bh/sf= "+alum_bh+"/"+alum_sf+")"
printf, 1, "#    Bolometric Luminosity=  "+bololum_tot+"  (bh/sf= "+bololum_bh+"/"+bololum_sf+")"
printf, 1, "#  theta      phi       Stellar        BH  "
printf, 1, "#  (deg)     (deg)     (L_solar)   (L_solar)"
close, 1




spawn, "cat "+tempfilename+" "+filename+" > "+finalfilename
spawn, "mv "+tempfilename+" "+frun+angledir+"/origs/"
spawn, "mv "+filename+" "+frun+angledir+"/origs/"




end




;=============================================================================








