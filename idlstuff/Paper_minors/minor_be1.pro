

;--------------------------------------------------------------------------
;
;  The standard, SPF01 formular
;
;   e= p1 * x^p0
;

function perform_sbefficiency_fit, mratios, burstes


mratios_tofit= mratios
burstes_tofit= burstes
weight= burstes_tofit*0.0 + 0.07
p1= [1.0,1.0]
;spec_sf_eff= MPFITFUN('specific_sf_law', mratios_tofit, burstes_tofit, weight, p1, BESTNORM=bestnorm, DOF=dof)
burst_eff_params= MPFITFUN('minor_func_samsbl', mratios_tofit, burstes_tofit, weight, p1, BESTNORM=bestnorm, DOF=dof)

redchi2= bestnorm/dof
print, "Reduced Chi^2 = ", redchi2
print, burst_eff_params

return, burst_eff_params


end



;--------------------------------------------------------------------------
;
;  This is our modified law that is
;
;
;   e =   p1 ( x - e_0)^p0    for x > e_0
;         0                   for x < e_0

function perform_sbefficiency_fit_v2, mratios, burstes


mratios_tofit= mratios
burstes_tofit= burstes
weight= burstes_tofit*0.0 + 0.07

	; initial guess
        guess= [1.0,1.0,0.1]
   
        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
        ; Normalization is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ; width is greater than 0
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.001
        ; mean is greater than 0
        pi[2].limited(0) = 1
        pi[2].limits(0) = 0.000000001

burst_eff_params= MPFITFUN('minor_func_samsbl2', mratios_tofit, burstes_tofit, weight, guess, BESTNORM=bestnorm, DOF=dof, PARINFO=pi)

redchi2= bestnorm/dof
print, "Reduced Chi^2 = ", redchi2
print, burst_eff_params

return, burst_eff_params


end



;--------------------------------------------------------------------------
;
;  Another trial law, which we assume is linear.
;
;     e= p1 * x + p0
;

function perform_sbefficiency_fit_v3, mratios, burstes


mratios_tofit= mratios
burstes_tofit= burstes
weight= burstes_tofit*0.0 + 0.07
p1= [1.0,1.0]
burst_eff_params= MPFITFUN('minor_func_samsbpl', mratios_tofit, burstes_tofit, weight, p1, BESTNORM=bestnorm, DOF=dof)

redchi2= bestnorm/dof
print, "Reduced Chi^2 = ", redchi2
print, burst_eff_params

return, burst_eff_params


end




;--------------------------------------------------------------------------
;
;  This is Avishai's kind of crazy formula that looks like
;
;    e = p1 ( x / (1 + x))
;
function perform_sbefficiency_fit_v4, mratios, burstes


mratios_tofit= mratios
burstes_tofit= burstes
weight= burstes_tofit*0.0 + 0.07
p1= [1.0]
burst_eff_params= MPFITFUN('minor_func_be', mratios_tofit, burstes_tofit, weight, p1, BESTNORM=bestnorm, DOF=dof)

redchi2= bestnorm/dof
print, "Reduced Chi^2 = ", redchi2
print, burst_eff_params

return, burst_eff_params


end



;--------------------------------------------------------------------------
;
;  This is our modified law, i.e., the same one as above, only we
;  manually fix e_0
;
;

function perform_sbefficiency_fit_v5, mratios, burstes


mratios_tofit= mratios
burstes_tofit= burstes
weight= burstes_tofit*0.0 + 0.07

; initial guess
guess= [1.0,1.0]

burst_eff_params= MPFITFUN('minor_func_samsbl3', mratios_tofit, burstes_tofit, weight, guess, BESTNORM=bestnorm, DOF=dof)

redchi2= bestnorm/dof
print, "Reduced Chi^2 = ", redchi2
print, burst_eff_params

return, burst_eff_params


end



;--------------------------------------------------------------------------
;
;  This is Patrik's linear suggestion.
;
;   e = p0 * alog10(x) + p1
;

function perform_sbefficiency_fit_v6, mratios, burstes


mratios_tofit= mratios
burstes_tofit= burstes
weight= burstes_tofit*0.0 + 0.07

; initial guess
guess= [1.0,1.0]

burst_eff_params= MPFITFUN('minor_func_samsblp', mratios_tofit, burstes_tofit, weight, guess, BESTNORM=bestnorm, DOF=dof)

redchi2= bestnorm/dof
print, "Reduced Chi^2 = ", redchi2
print, burst_eff_params

return, burst_eff_params


end




;--------------------------------------------------------------------------







pro calculate_burst_params, frun, MaxSFRTime= MaxSFRTime, AvgSFR= AvgSFR, $
            MassFormedStars= MassFormedStars, MaximumSFR= MaximumSFR, $
            OrigGasMass= OrigGasMass, FinalGasMass= FinalGasMass, $
            GasConsumption= GasConsumption, GasCFrac= GasCFrac


        ; get sfr data
            sfrfile= '/home2/tcox/upsand/sfr/'+fload_getid(frun)+'.sfr'

            openr, 1, sfrfile, ERROR=err

            if (err NE 0) then begin
                print, "  "
                print, "Problem: ",!ERR_STRING
                print, "  "
                close, 1
                ERR= 0
                sfrtime= [0]
                sfrsfr= [0]
                sfrmsg= "No Star Formation"
            endif else begin
                close, 1
                print, "opening: ",sfrfile
                sfrdata= read_ascii(sfrfile)
                sfrtime= sfrdata.field1[0,*]
                sfrsfr= sfrdata.field1[1,*]
                sfrmfs= sfrdata.field1[3,*]
                n_cols= n_elements(sfrdata.field1[*,0])
                n_rows= n_elements(sfrdata.field1[0,*])
                finalnewstarmass= sfrmfs[n_rows-1]
                if n_cols ge 6 then sfrgasmass= sfrdata.field1[6,*] else sfrgasmass=[20.0,20.0-finalnewstarmass]
                sfrmsg= ''
            endelse


            t1per= 0
            idx=where(sfrtime ge 1.0)
            if idx(0) ne -1 and n_cols gt 5 then begin
                gm= sfrgasmass[0]-sfrgasmass(idx(0))
                t1per= 100.0*gm/sfrgasmass[0]
            endif

            t4per= 0
            idx=where(sfrtime ge 4.0)
            if idx(0) ne -1 and n_cols gt 5 then begin
                gm= sfrgasmass[0]-sfrgasmass(idx(0))
                t4per= 100.0*gm/sfrgasmass[0]
            endif


            sfraftertmerg= 0.0
            smaftertmerg= 0.0
            taftertmerg= 0.0
            tmerger= 0.0
            ;if n_elements(tmerg) gt 1 then begin
            ;    idx= where(sfrtime ge tmerg[i]+0.2)
            ;    if idx(0) eq -1 then begin
            ;            sfraftertmerg= sfrsfr[n_rows-1]
            ;            smaftertmerg= sfrmfs[n_rows-1]
            ;            tmerger= sfrtime[n_rows-1]
            ;            taftertmerg= tmerger
            ;    endif else begin
            ;            sfraftertmerg= sfrsfr[idx(0)]
            ;            smaftertmerg= sfrmfs[idx(0)]
            ;            tmerger= tmerg[i]
            ;            taftertmerg= sfrtime[idx(0)]
            ;    endelse
            ;endif


            MaximumSFR= max(sfrsfr)
            idx= where(sfrsfr eq MaximumSFR)
            MaxSFRTime= sfrtime(idx)
            AvgSFR= total(sfrsfr)/n_rows
            n_mfs= n_elements(sfrgasmass)
            MassFormedStars= sfrmfs[n_mfs-1]
            print, "-------------------------------------"
            print, "maximum sfr= ", MaximumSFR, " occurs at ", MaxSFRTime
            print, "average sfr= ", AvgSFR
            print, "mass of new stars formed= ", MassFormedStars

            OrigGasMass= sfrgasmass[0]
            FinalGasMass= sfrgasmass[n_mfs-1]
            GasConsumption= OrigGasMass - FinalGasMass
            GasCFrac= GasConsumption / OrigGasMass
            print, "-------------------------------------"
            print, "original gas mass   =", OrigGasMass
            print, "remnant gas mass    =", FinalGasMass
            print, "gas consumed        =", GasConsumption
            print, "gas consump f       =", GasCFrac
            print, "-------------------------------------"




end



;--------------------------------------------------------------------------



;--------------------------------------------------------------------------


pro grab_sf_info, fruns, MaxSFRlist, MaxSFRTimelist, AvgSFRlist, $
                        MassFSlist, OrigGMlist, FinalGMlist, $
                        CGlist, GCFlist

   numfrun= n_elements(fruns)

   MaxSFRlist= fltarr(numfrun)
   MaxSFRTimelist= fltarr(numfrun)
   AvgSFRlist= fltarr(numfrun)
   MassFSlist= fltarr(numfrun)
   OrigGMlist= fltarr(numfrun)
   FinalGMlist= fltarr(numfrun)
   CGlist= fltarr(numfrun)
   GCFlist= fltarr(numfrun)

   for i=0, numfrun-1 do begin

        calculate_burst_params, fruns[i], MaxSFRTime= MaxSFRTime, AvgSFR= AvgSFR, $
            MassFormedStars= MassFormedStars, $
            OrigGasMass= OrigGasMass, FinalGasMass= FinalGasMass, $
            MaximumSFR= MaximumSFR, $
            GasConsumption= GasConsumption, GasCFrac= GasCFrac

        MaxSFRlist[i]= MaximumSFR
        MaxSFRTimelist[i]= MaxSFRTime
        AvgSFRlist[i]= AvgSFR
        MassFSlist[i]= MassFormedStars
        OrigGMlist[i]= OrigGasMass
        FinalGMlist[i]= FinalGasMass
        CGlist[i]= GasConsumption
        GCFlist[i]= GasCFrac
   endfor

end





;--------------------------------------------------------------------------




function calculate_be, mergerrun, gal1run, gal2run

        fruns= [mergerrun, gal1run, gal2run]
        grab_sf_info, fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist

        GasConsumption= GClist[1] + GClist[2]
        OrigGM= j5(1) + j5(2)
        isogasc= GasConsumption / OrigGM

        be= GCFlist[0] - isogasc

	OrigGM_diff= j5(0) - OrigGM
	if OrigGM_diff gt 0.001 then begin
		print, "Hey, why don't the gas masses add up?"
	endif

	return, be

end







;--------------------------------------------------------------------------








pro minor_be, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "minor_be, junk"
   return
endif


l0= 0
l1= 0
l2= 0
l3= 0
l4= 0





;   Extent of gaseous disk
;------------------------------------------
;do_gasdiskext= 1 
do_gasdiskext= 0 
if do_gasdiskext eq 1 then begin

        xvar_0= 1./[1.0, 2.3, 5.8, 22.7]
        xvar_1= 1./[1.0, 2.3, 5.8, 22.7]
        xvar_2= 1./[1.0, 2.3, 5.8, 22.7]
        xvar_3= 1./[1.0, 2.3, 5.8, 22.7]
        xvar_4= 1./[1.0, 2.3, 5.8, 22.7]


        ; ----------------------
        l0= 1
        msg0= '!9a!6= 1'
        ; get isolated gas consumption
        Isofruns= ["G3Rd4", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
        grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        GasConsumption= GClist + GClist[0]
        OrigGM= j5 + j5(0)
        isogasc= GasConsumption / OrigGM
        ; now get the merger gas consumption
        G3fruns= ["G3Rd4G3Rd4", "G3Rd4G2", "G3Rd4G1", "G3Rd4G0"]
        grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;
        be_0= GCFlist - isogasc

        ; ----------------------
        l1= 1
        msg1= '!9a!6= 3'
        ; get isolated gas consumption
        Isofruns= ["G3Rd4e", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
        grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        GasConsumption= GClist + GClist[0]
        OrigGM= j5 + j5(0)
        isogasc= GasConsumption / OrigGM
        ; now get the merger gas consumption
        G3fruns= ["G3Rd4eG3", "G3Rd4eG2", "G3Rd4eG1", "G3Rd4eG0"]
        grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;
        be_1= GCFlist - isogasc


	; isolated's
	;Isofruns= ["G3Rd4e", "G3Rd4f", "G3Rd4g", "G3Rd4h", "G2im-u1a"]
	;grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
	;GasConsumption= GClist + GClist[4]
        ;OrigGM= j5 + j5(4)
        ;isogasc= GasConsumption / OrigGM
	;isogasc= isogasc[0:3]
        ; now get the merger gas consumption
        ;G3fruns= ["G3Rd4eG2", "G3Rd4fG2", "G3Rd4gG2", "G3Rd4hG2"]
        ;grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;
        ;be_eg= GCFlist - isogasc


	; ----------------------
        l2= 1
        msg2= '!9a!6= 3 (n0)'
        ; get isolated gas consumption
        Isofruns= ["G3Rd4e_n0", "G2im-u2", "G1i-u2", "G0i-u2"]
        grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        GasConsumption= GClist + GClist[0]
        OrigGM= j5 + j5(0)
        isogasc= GasConsumption / OrigGM
        ; now get the merger gas consumption
        G3fruns= ["G3Rd4eG3_n0", "G3Rd4eG2_n0", "G3Rd4eG1_n0", "G3Rd4eG0_n0"]
        grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;
        be_2= GCFlist - isogasc



endif






;
;   Bulge-to-Disk Ratio
;------------------------------------------
;do_b_to_disk= 1
do_b_to_disk= 0
if do_b_to_disk eq 1 then begin

	xvar_0= 1./[1.0, 2.3, 5.8, 22.7]
	xvar_1= 1./[1.0, 2.3, 5.8, 22.7]
	xvar_2= 1./[1.0, 2.3, 5.8, 22.7]
	xvar_3= 1./[1.0, 2.3, 5.8, 22.7]
	xvar_4= 1./[1.0, 2.3, 5.8, 22.7]


	; ----------------------
	l0= 1
	msg0= 'B/D= 0'
	; get isolated gas consumption
	Isofruns= ["G3il-u1a", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
	grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
	GasConsumption= GClist + GClist[0]
        OrigGM= j5 + j5(0)
        isogasc= GasConsumption / OrigGM
	; now get the merger gas consumption
	G3fruns= ["G3blv5G3blv5-u1", "G3blv5G2-u3", "G3blv5G1-u3", "G3blv5G0-u3"]
	grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
	;
	be_0= GCFlist - isogasc


        ;l2= 1
        ;msg2= 'B/D= 0 (G4)'
        ; get isolated gas consumption
        ;Isofruns= ["Sb_n0", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
        ;grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;GasConsumption= GClist + GClist[0]
        ;OrigGM= j5 + j5(0)
        ;isogasc= GasConsumption / OrigGM
        ; now get the merger gas consumption
        ;G3fruns= ["SbSb_o2_n0", "SbG2_o2_n0", "SbG1_o2_n0", "SbG0_o2_n0"]
        ;grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;
        ;be_2= GCFlist - isogasc


	; ----------------------
	l1= 1
	msg1= '!6B/D= 0.2 (fiducial)'
	; get isolated gas consumption
	Isofruns= ["G3il-u1a", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
	grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
	GasConsumption= GClist + GClist[0]
        OrigGM= j5 + j5(0)
        isogasc= GasConsumption / OrigGM
	; now get the merger gas consumption
	G3fruns= ["G3G3b-u1", "G3G2-u3", "G3G1-u3", "G3G0e-u3"]
	grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
	;
	be_1= GCFlist - isogasc


	; ----------------------
	;l2= 1
	;msg2= 'B/D= 0.25'
        ; get isolated gas consumption
        ;Isofruns= ["G3BT1i-u1", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
        ;grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
	;GasConsumption= GClist + GClist[0]
        ;OrigGM= j5 + j5(0)
        ;isogasc= GasConsumption / OrigGM
	; now get the merger gas consumption
	;G3fruns= ["G3BT1G3BT1-u1", "G3BT1G2-u1", "G3BT1G1-u1", "G3BT1G0-u1"]
	;grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist 
        ;
        ;be_2= GCFlist - isogasc


	; ----------------------
	l2= 1
	msg2= 'B/D= 0.5'
        ; get isolated gas consumption
        Isofruns= ["G3BT2i-u1", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
        grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
	GasConsumption= GClist + GClist[0]
        OrigGM= j5 + j5(0)
        isogasc= GasConsumption / OrigGM
	; now get the merger gas consumption
	G3fruns= ["G3BT2G3BT2-u1", "G3BT2G2-u1", "G3BT2G1-u1", "G3BT2G0-u1"]
	grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist 
        ;
        be_2= GCFlist - isogasc


	; ----------------------
	;l4= 1
	;msg4= 'B/D= 1.0'
        ; get isolated gas consumption
        ;Isofruns= ["G3BT3i-u1", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
        ;grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
	;GasConsumption= GClist + GClist[0]
        ;OrigGM= j5 + j5(0)
        ;isogasc= GasConsumption / OrigGM
	; now get the merger gas consumption
	;G3fruns= ["G3BT3G3BT3-u1", "G3BT3G2-u1", "G3BT3G1-u1", "G3BT3G0-u1"]
	;grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist 
        ;
        ;be_4= GCFlist - isogasc

endif



;
;   Gas Fraction Ratio
;------------------------------------------
do_gf= 1
;do_gf= 0
if do_gf eq 1 then begin

	xvar= 1./[1.0, 2.3, 5.8, 22.7]

	; ----------------------
	l1= 1
	msg1= '0.2'
	; get isolated gas consumption
	Isofruns= ["G3il-u1a", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
	grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
	GasConsumption= GClist + GClist[0]
        OrigGM= j5 + j5(0)
        isogasc= GasConsumption / OrigGM
	; now get the merger gas consumption
	xvar_1= 1.0 / [1.0, 2.3, 5.8, 22.7]
	G3fruns= ["G3G3b-u1", "G3G2-u3", "G3G1-u3", "G3G0e-u3"]
	grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
	;
	be_1= GCFlist - isogasc


	; ----------------------
	l2= 1
	msg2= '0.5'
        ; get isolated gas consumption
        Isofruns= ["G3gf1i-u1", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
        grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
	GasConsumption= GClist + GClist[0]
        OrigGM= j5 + j5(0)
        isogasc= GasConsumption / OrigGM
	; now get the merger gas consumption
	xvar_2= 1.0 / [1.0, 2.3, 5.8, 22.7]
	G3fruns= ["G3gf1G3gf1b-u1", "G3gf1G2-u1", "G3gf1G1-u1", "G3gf1G0-u1"]
	grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist 
        ;
        be_2= GCFlist - isogasc


	; ----------------------
	l3= 1
	msg3= '0.78'
        ; get isolated gas consumption
        Isofruns= ["G3gf2i-u1", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
        grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
	GasConsumption= GClist + GClist[0]
        OrigGM= j5 + j5(0)
        isogasc= GasConsumption / OrigGM
	; now get the merger gas consumption
	xvar_3= 1.0 / [1.0, 2.3, 5.8, 22.7]
	G3fruns= ["G3gf2G3gf2b-u1", "G3gf2G2-u1", "G3gf2G1-u1", "G3gf2G0-u1"]
	grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist 
        ;
        be_3= GCFlist - isogasc


endif




;
;   Standard Mergers
;------------------------------------------
;do_std= 1
do_std= 0
if do_std eq 1 then begin

        ; ----------------------
        l1= 1
        msg1= 'G3Gx'
        ; get isolated gas consumption
        Isofruns= ["G3il-u1a", "G2im-u1a", "G1i-u1a", "G0i-u1a"]
        ;Isofruns= ["G3_n0", "G2im-u2", "G1i-u2", "G0i-u2"]
        grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;isogasc= GCFlist + GCFlist[0]
        GasConsumption= GClist + GClist[0]
	OrigGM= j5 + j5(0)
	isogasc= GasConsumption / OrigGM
        ; now get the merger gas consumption
        xvar_1= 1.0 / [1.0, 2.3, 5.8, 22.7]
        G3fruns= ["G3G3b-u1", "G3G2-u3", "G3G1-u3", "G3G0e-u3"]
        ;G3fruns= ["G3G3b-u2", "G3G2-u4", "G3G1-u4", "G3G0e-u4"]
        grab_sf_info, G3fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;
        be_1= GCFlist - isogasc
	;idx= where(be_1 lt 0.0)
	;if idx(0) ne -1 then be_1(idx)= 0

;IDL> print, GCFlist
;     0.696758     0.552270     0.313528     0.245250
;
;IDL> print, isogasc
;     0.465085     0.464670     0.373773     0.295343
;IDL> 



        ; ----------------------
        l2= 1
        msg2= 'G2Gx'
        ; get isolated gas consumption
        Isofruns= ["G2im-u1a", "G1i-u1a", "G0i-u1a"]
        grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;isogasc= GCFlist + GCFlist[0]
        GasConsumption= GClist + GClist[0]
	OrigGM= j5 + j5(0)
	isogasc= GasConsumption / OrigGM
        ; now get the merger gas consumption
        xvar_2= 1.0 / [1.0, 2.6, 10.0]
        G2fruns= ["G2G2-u1", "G2G1-u3", "G2G0-u3"]
        grab_sf_info, G2fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist 
        ;
        be_2= GCFlist - isogasc
	;idx= where(be_2 lt 0.0)
	;if idx(0) ne -1 then be_2(idx)= 0


        ; ----------------------
        l3= 1
        msg3= 'G1Gx'
        ; get isolated gas consumption
        Isofruns= ["G1i-u1a", "G0i-u1a"]
        grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;isogasc= GCFlist + GCFlist[0]
        GasConsumption= GClist + GClist[0]
	OrigGM= j5 + j5(0)
	isogasc= GasConsumption / OrigGM
        ; now get the merger gas consumption
        xvar_3= 1.0 / [1.0, 3.9]
        G1fruns= ["G1G1a-u1", "G1G0-u3"]
        grab_sf_info, G1fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;
        be_3= GCFlist - isogasc
	;idx= where(be_3 lt 0.0)
	;if idx(0) ne -1 then be_3(idx)= 0


        ; ----------------------
        l4= 1
        msg4= 'G0G0'
        ; get isolated gas consumption
        Isofruns= ["G0i-u1a"]
        grab_sf_info, Isofruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;isogasc= GCFlist + GCFlist[0]
        GasConsumption= GClist + GClist[0]
	OrigGM= j5 + j5(0)
	isogasc= GasConsumption / OrigGM
        ; now get the merger gas consumption
        xvar_4= 1.0 / [1.0]
        G0fruns= ["G0G0a-u1"]
        grab_sf_info, G0fruns, j1, j2, j3, j4, j5, j6, GClist, GCFlist
        ;
        be_4= GCFlist - isogasc
	;idx= where(be_4 lt 0.0)
	;if idx(0) ne -1 then be_4(idx)= 0


endif




;============================================
;============================================
;
;
;
;============================================
;============================================


xaxistitle= "!6Mass Ratio (M!Dsat!N/M!Dprimary!N)"
xmax= 2.0
;xmax= 1.05
xmin= 0.01
;xmin= 0.0

yaxistitle= "!6burst efficiency, !8e!6"
;ymax= 1.0
ymax= 0.8
ymin= 0.0
;ymin= -0.2

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename='minburste.eps', colortable=4

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, charthick=3.0, xthick=4.0, ythick=4.0, color= 0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /xlog


;----------------------------------------

        ;if l0 ne 0 then oplot, xvar_0, be_0, psym=-4, thick=3.0, color= 200, symsize=1.5
        ;if l1 ne 0 then oplot, xvar_1, be_1, psym=-2, thick=3.0, color= 150, symsize=1.5
        ;if l2 ne 0 then oplot, xvar_2, be_2, psym=-5, thick=3.0, color= 100, symsize=1.5
        ;if l3 ne 0 then oplot, xvar_3, be_3, psym=-7, thick=3.0, color= 50, symsize=1.5
        ;if l4 ne 0 then oplot, xvar_4, be_4, psym=-6, thick=3.0, color= 0, symsize=1.5
        if l0 ne 0 then oplot, xvar_0, be_0, psym=4, thick=3.0, color= 200, symsize=1.5
        if l1 ne 0 then oplot, xvar_1, be_1, psym=2, thick=3.0, color= 150, symsize=1.5
        if l2 ne 0 then oplot, xvar_2, be_2, psym=5, thick=3.0, color= 100, symsize=1.5
        if l3 ne 0 then oplot, xvar_3, be_3, psym=7, thick=3.0, color= 50, symsize=1.5
        if l4 ne 0 then oplot, xvar_4, be_4, psym=6, thick=3.0, color= 0, symsize=1.5


        x0= 0.30
        y0= 0.82
        if l0 ne 0 then xyouts, x0+0.05, y0+0.05, msg0, /normal, charthick=3.0, size=1.3, color= 200
        if l1 ne 0 then xyouts, x0+0.05, y0, msg1, /normal, charthick=3.0, size=1.3, color= 150
        if l2 ne 0 then xyouts, x0+0.05, y0-0.05, msg2, /normal, charthick=3.0, size=1.3, color= 100
        if l3 ne 0 then xyouts, x0+0.05, y0-0.10, msg3, /normal, charthick=3.0, size=1.3, color= 50
	if l4 ne 0 then xyouts, x0+0.05, y0-0.15, msg4, /normal, charthick=3.0, size=1.3, color= 0
	if l0 ne 0 then oplot, [0.018], [0.73], psym=4, thick=3.0, symsize=1.5, color=200
	if l1 ne 0 then oplot, [0.018], [0.68], psym=2, thick=3.0, symsize=1.5, color=150
	if l2 ne 0 then oplot, [0.018], [0.63], psym=5, thick=3.0, symsize=1.5, color=100
	if l3 ne 0 then oplot, [0.018], [0.58], psym=7, thick=3.0, symsize=1.5, color=50
	if l4 ne 0 then oplot, [0.018], [0.53], psym=6, thick=3.0, symsize=1.5, color= 0

	;if l0 ne 0 then oplot, [0.013,0.026], [0.73,0.73], psym=-3, thick=10.0, color=200, linestyle= 1
	;if l1 ne 0 then oplot, [0.013,0.026], [0.68,0.68], psym=-3, thick=10.0, color=150, linestyle= 0
	;if l2 ne 0 then oplot, [0.013,0.026], [0.63,0.63], psym=-3, thick=10.0, color=100, linestyle= 2

	;x0= 0.01 & xs= 0.21
	;y0= 0.92 & ys= 0.25
	;oplot, [x0,x0+xs,x0+xs,x0,x0],[y0,y0,y0-ys,y0-ys,y0], thick=3.0, color= 0, psym=-3

	oplot, [0.03], [0.4], psym=-3, color=0, thick= 3.0, symsize= 1.5
	oploterror, [0.03], [0.4], [0.07] , psym=-3, errcolor=0, color=0, thick= 3.0, errthick= 3.0



	; plot fits
	;------------
	;x= (alog10(xmax)-alog10(xmin)) * findgen(11)/10.0 + alog10(xmin)
	x= (alog10(1.0)-alog10(xmin)) * findgen(21)/20.0 + alog10(xmin)
	x= 10^x

	x=[x,0.09]
	x= x(sort(x))

	;do_allinone_fit= 1
	do_allinone_fit= 0
	if do_allinone_fit eq 1 then begin
		xvar= [xvar_1, xvar_2, xvar_3, xvar_4]
		be= [be_1, be_2, be_3, be_4]


		; standard SPF formulae: p0 x ^p2
		fitresult= perform_sbefficiency_fit(xvar, be)
		y= minor_func_samsbl(x,fitresult) 
                oplot, x, y, psym=-3, color= 0, thick=10.0, linestyle= 0

                ; standard SPF formulae w/ non-zero crossing: p0 ( x - 0.07 )^p2
                fitresult= perform_sbefficiency_fit_v5(xvar, be)
                y= minor_func_samsbl3(x,fitresult) 
                oplot, x, y, psym=-3, color= 0, thick=10.0, linestyle= 1

		; standard SPF formulae w/ non-zero crossing: p0 ( x - p1 )^p2
		;fitresult= perform_sbefficiency_fit_v2(xvar, be)
		;y= minor_func_samsbl2(x,fitresult) 
                ;oplot, x, y, psym=-3, color= 100, thick=10.0, linestyle= 1

		; broken power-law: p0 + x * p1
		;fitresult= perform_sbefficiency_fit_v3(xvar, be)
		;y= minor_func_samsbpl(x,fitresult) 
                ;oplot, x, y, psym=-3, color= 150, thick=10.0, linestyle= 2

		; avishai's suggestion: p0 x / (1+x)
		;fitresult= perform_sbefficiency_fit_v4(xvar, be)
		;y= minor_func_be(x,fitresult) 
                ;oplot, x, y, psym=-3, color= 200, thick=10.0, linestyle= 3

		; patrik's suggestion: p0 * alog10(x) + p1
		;fitresult= perform_sbefficiency_fit_v6(xvar, be)
		;y= minor_func_samsblp(x,fitresult) 
                ;oplot, x, y, psym=-3, color= 200, thick=10.0, linestyle= 3

		;oplot, [xmin, xmax], [0.0, 0.0], color= 0, thick=3.0, linestyle= 1
	endif


	;do_separate_fits= 0
	do_separate_fits= 1
	if do_separate_fits eq 1 then begin
	   if l0 ne 0 then begin
		fitresult= perform_sbefficiency_fit(xvar_0, be_0)
		y= minor_func_samsbl(x,fitresult)
		oplot, x, y, psym=-3, color= 200, thick=10.0, linestyle= 1
	   endif

	   if l1 ne 0 then begin
		fitresult= perform_sbefficiency_fit(xvar_1, be_1)
		y= minor_func_samsbl(x,fitresult)
		oplot, x, y, psym=-3, color= 150, thick=10.0, linestyle= 0
	   endif

	   if l2 ne 0 then begin
		fitresult= perform_sbefficiency_fit(xvar_2, be_2)
		y= minor_func_samsbl(x,fitresult)
		oplot, x, y, psym=-3, color= 100, thick=10.0, linestyle= 2
	   endif

	   if l3 ne 0 and n_elements(xvar_3) gt 1 then begin
		fitresult= perform_sbefficiency_fit(xvar_3, be_3)
		y= minor_func_samsbl(x,fitresult)
		oplot, x, y, psym=-3, color= 50, thick=10.0, linestyle= 0
	   endif

	   if l4 ne 0 and n_elements(xvar_4) gt 1 then begin
		fitresult= perform_sbefficiency_fit(xvar_4, be_4)
		y= minor_func_samsbl(x,fitresult)
		oplot, x, y, psym=-3, color= 0, thick=10.0, linestyle= 0
	   endif
	endif


;============================================

device, /close




end




