;------------------------------------------------------------------------
;
;    Random Procedures related to Velocity Stuff
;
;
;
;
;------------------------------------------------------------------------




;   average re
; --------------------------
function fload_avg_re, fruns, ddir, all=all, reerr=reerr
   ns= n_elements(fruns)
   allre= fltarr(ns)
   reerr= fltarr(ns)
   res= [-1]
   for i=0,ns-1 do begin
	read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir[i]
	res= [res, transpose(Re)]
	allre[i]= mean(Re)
	reerr[i]= sqrt(variance(Re))
   endfor
   res= res[1:n_elements(res)-1]
   if keyword_set(all) then return, res
   return, allre
end


;   average ellipticities
; --------------------------
function fload_avg_es, fruns, ddir, all=all, eserr=eserr
   ns= n_elements(fruns)
   alles= fltarr(ns)
   eserr= fltarr(ns)
   es= [-1]
   for i=0,ns-1 do begin
        read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir[i]
	es= [es,transpose(ellipticity)]
        alles[i]= mean(ellipticity)
	eserr[i]= sqrt(variance(ellipticity))
   endfor
   es= es[1:n_elements(es)-1]
   if keyword_set(all) then return, es
   return, alles
end


;   average vmax
; --------------------------
function fload_avg_vmax, fruns, ddir, all=all
   ns= n_elements(fruns)
   allvmax= fltarr(ns)
   vs= [-1]
   for i=0,ns-1 do begin
        read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir[i]
	vs= [vs, transpose(vmax)]
        allvmax[i]= mean(vmax)
   endfor
   vs= vs[1:n_elements(vs)-1]
   if keyword_set(all) then return, vs
   return, allvmax
end


;   average sigmas
; --------------------------
function fload_avg_sig, fruns, ddir, all=all
   ns= n_elements(fruns)
   allsig= fltarr(ns)
   sigs= [-1]
   for i=0,ns-1 do begin
        read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir[i]
	sigs= [sigs, transpose(sigma)]
        allsig[i]= mean(sigma)
   endfor
   sigs= sigs[1:n_elements(sigs)-1]
   if keyword_set(all) then return, sigs
   return, allsig
end


;   average vmax/sigmas
; --------------------------
function fload_avg_vs, fruns, ddir, all=all
   ns= n_elements(fruns)
   allvs= fltarr(ns)
   for i=0,ns-1 do begin
        read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir[i]
        allvs[i]= mean(vmax/sigma)
   endfor
   return, allvs
end


;   average (vmax/sigma)*
; --------------------------
function fload_avg_vsst, fruns, ddir, all=all, vserr=vserr
   ns= n_elements(fruns)
   allvsst= fltarr(ns)
   vserr= fltarr(ns)
   vsigst= [-1]
   for i=0,ns-1 do begin
        read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir[i]
	thisvsigst= (vmax/sigma)/(sqrt(ellipticity/(1-ellipticity)))
	;vsigst= [vsigst,(vmax/sigma)/(sqrt(ellipticity/(1-ellipticity)))]
	vsigst= [vsigst,transpose(thisvsigst)]
        allvsst[i]= mean(thisvsigst)
	vserr[i]= sqrt(variance(thisvsigst))
   endfor
   vsigst= vsigst[1:n_elements(vsigst)-1]
   if keyword_set(all) then return, vsigst
   return, allvsst
end



;   return Triaxiality
; --------------------------
function fload_Ts, fruns
   ns= n_elements(fruns)
   allTTs= fltarr(ns)
   for i=0,ns-1 do begin
	if strmid(fruns[i],0,1) eq 'A' then fruns[i]='As/'+fruns[i]
	if strmid(fruns[i],0,1) eq 'c' then fruns[i]='collisionless/'+fruns[i]
	read_moi_file, fruns[i], i_xx, i_yy, i_zz, i_xy, i_xz, i_yz, bb, cc, TT
        allTTs[i]=  TT
   endfor
   return, allTTs
end


;   return delta
; --------------------------
function fload_deltas, fruns
   ns= n_elements(fruns)
   alldeltas= fltarr(ns)
   for i=0,ns-1 do begin
	if strmid(fruns[i],0,1) eq 'c' then fruns[i]='collisionless/'+fruns[i]
        read_tvt_file, fruns[i], pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz
	;delta= 1.0 - pi_zz / sqrt(pi_xx*pi_xx + pi_yy*pi_yy)
	delta= 1.0 - 2.0*pi_zz / (pi_xx + pi_yy)
        alldeltas[i]=  abs(delta)
   endfor
   return, alldeltas
end


;   kinematic misalignments
; ---------------------------
function fload_kinematicmas, fruns, cless=cless, all=all
   ddir_major= 0
   ;ddir_major= 13
   ddir_minor= 1
   if keyword_set(cless) then begin
	ddir_major= 2
	ddir_minor= 3
   endif
   ns= n_elements(fruns)
   allkms= fltarr(ns)
   kms= [-1]
   for i=0,ns-1 do begin
	read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir_major
	vmaj= vmax
	read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir_minor
	vmin= vmax
	this_kms= abs(atan(vmin/vmaj))*180.0/!PI
	print, "mean kinematic misalignment, <Phi>= ",mean(this_kms)
	mu= vmin / sqrt(vmin*vmin + vmax*vmax)
	print, "                              <mu>= ",mean(mu)
	allkms[i]= mean(this_kms)
	kms= [kms,transpose(this_kms)]
   endfor
   kms= kms[1:n_elements(kms)-1]
   if keyword_set(all) then return, kms
   return, allkms
end












;------------------------------------------------------------------------
;
;------------------------------------------------------------------------






; ----------------------------
;  Read vsig.txt file
; ----------------------------
pro read_vsig_file, frun, phi, theta, Re, PA, ellipticity, vmax, sigma, $
			ddir=ddir

if not keyword_set(ddir) then ddir= 0

case ddir of
	0: datadir='vc3vc3_major'
	1: datadir='vc3vc3_minor'
	2: datadir='cvc3vc3_major'
	3: datadir='cvc3vc3_minor'
	4: datadir='masses'
	5: datadir='vc3vc3_new'
	6: datadir='vc3vc3_old'
	7: datadir='vc3vc3_20gas'
	8: datadir='vc3vc3_nobh'
	9: datadir='vc3vc3_orbits'
	10: datadir='vc3vc3_resolution'
	11: datadir='vc3vc3_emerger'
	12: datadir='vc3vc3_gasfrac'
	13: datadir='vc3vc3_mass'
	14: datadir='gfs'
	15: datadir='ds'
	else: break
endcase

junk= strsplit(frun,'/', /extract,count=count)
nfrun= junk(count-1)

;spawn, '/bin/ls /raid4/sdutta/projections/'+datadir+'/'+nfrun+'/*vc*.prof ',result
;spawn, '/bin/ls /raid4/sdutta/projections/'+datadir+'/'+nfrun+'/'+nfrun+'*.prof ',result
spawn, '/bin/ls /raid4/sdutta/projections/'+datadir+'/'+nfrun+'/'+strmid(nfrun,0,3)+'*.prof ',result

proffile=strcompress(result[0],/remove_all)

;filename= '/raid4/sdutta/projections/'+datadir+'/'+frun+'/'+frun+'.prof'
;filename= '/raid4/sdutta/projections/'+datadir+'/'+frun+'/'+proffile
filename= proffile

print, " --------------------------- "
print, 'opening: ',filename

spawn, "wc "+filename,result
lines=long(result)
lines=lines(0)-1
if lines GT 0 then vsig_data= fltarr(7,lines)

openr, 1, filename

junk=''
readf, 1, junk
;readf, 1, junk
;readf, 1, junk
;readf, 1, junk
;readf, 1, junk
readf, 1, vsig_data
close, 1


phi= vsig_data[0,*]
theta= vsig_data[1,*]
Re= vsig_data[2,*]
PA= vsig_data[3,*]
ellipticity= vsig_data[4,*]
vmax= vsig_data[5,*]
sigma= vsig_data[6,*]

print, "<r_e>=",mean(Re), "  +/-", sqrt(variance(Re))
print, "<e>=",mean(ellipticity), "  +/-", sqrt(variance(ellipticity))
print, "<v>=",mean(vmax), "  +/-", sqrt(variance(vmax))
print, "<sigma>=",mean(sigma), "  +/-", sqrt(variance(sigma))

end










;===================================================================
;===================================================================







;
;-----------------------------------------
pro process_vmin, frun, vmax_maj, vmax_min

	;ddir1= 0
	ddir1= 13
	ddir2= 1
	if strmid(frun,0,1) eq 'c' then begin
		ddir1= 2 & ddir2= 3
	endif

	read_vsig_file, frun, phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir1
	vmax_maj= vmax
	read_vsig_file, frun, phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir2
	vmax_min= vmax

vtot= sqrt(vmax_maj*vmax_maj + vmax_min*vmax_min)
vminfrac= vmax_min/vtot
print, " --------------------------- "
print, "<v_maj>=",mean(vmax_maj), "  +/-", sqrt(variance(vmax_maj))
print, "<v_min>=",mean(vmax_min), "  +/-", sqrt(variance(vmax_min))
print, "<v_tot>=",mean(vtot), "  +/-", sqrt(variance(vtot))
print, "<f_min>=",mean(vminfrac), "  +/-", sqrt(variance(vminfrac))
print, "<v_maj/sigma>=",mean(vmax_maj/sigma), "  +/-", sqrt(variance(vmax_maj/sigma))

	

end






;===================================================================================
;===================================================================================
;===================================================================================
;===================================================================================






;=================================================================================









pro rems_v_something, junk



if not keyword_set(junk) then begin
        print, " "
        print, " rems_v_something, junk"
        print, " "
        print, " "
        return
endif


; =====================
; =====================

secondgrp= 0
thirdgrp= 0

;------------------------------------------------------------------------------------


; gas fraction
; -------------
do_gf= 0
;do_gf= 1
if do_gf eq 1 then begin
	;xvar= [0,0.05,0.20,0.40,0.60,0.80,1.00]
	xvar= [0,0.05,0.20,0.40,0.60,0.80,1.00]
	xmax= 1.1 & xmin= -0.1
	xaxistitle='Progenitor Gas Fraction'


	; gas fractions
	cmt1= 'orbit e'
        fruns= ['cvc3vc3e','gfs/vc3vc3z_e','gfs/vc3vc3y_e','gfs/vc3vc3x2_e', $
		'gfs/vc3vc3w_e','gfs/vc3vc3v_e','gfs/vc3vc3u_e']
        ddir= [2,14,14,14,14,14,14]
	re_avgs=     fload_avg_re(fruns,ddir, reerr=reerr)
	re_avgs= re_avgs * 1.3 + 0.3
	es_avgs=     fload_avg_es(fruns,ddir, eserr=eserr)
	;vmax_avgs=   fload_avg_vmax(fruns,ddir)
	;sig_avgs=    fload_avg_sig(fruns,ddir)
	;vs_avgs=     fload_avg_vs(fruns,ddir)
	vs_avgs=     fload_avg_vsst(fruns,ddir, vserr=vserr)
	TTs_avgs=    fload_Ts(fruns)
	deltas_avgs= fload_deltas(fruns)


	; second group
	cmt2= 'orbit h'
	secondgrp= 1
	xvar2= xvar
        fruns= ['cvc3vc3h','gfs/vc3vc3z_h','gfs/vc3vc3y_h','gfs/vc3vc3x2_h', $
		'gfs/vc3vc3w_h','gfs/vc3vc3v_h','gfs/vc3vc3u_h']
        ddir= [2,14,14,14,14,14,14]
        re_avgs2=     fload_avg_re(fruns,ddir, reerr=reerr2)
	re_avgs2= re_avgs2 * 1.3 + 0.3
        es_avgs2=     fload_avg_es(fruns,ddir, eserr=eserr2)
        ;vmax_avgs2=   fload_avg_vmax(fruns,ddir)
        ;sig_avgs2=    fload_avg_sig(fruns,ddir)
        ;vs_avgs2=     fload_avg_vs(fruns,ddir)
        vs_avgs2=     fload_avg_vsst(fruns,ddir, vserr=vserr2)
        TTs_avgs2=    fload_Ts(fruns)
        deltas_avgs2= fload_deltas(fruns)


	; third group
	cmt3= 'orbit k'
	thirdgrp= 1
	xvar3= xvar
        fruns= ['cvc3vc3k','gfs/vc3vc3z_k','gfs/vc3vc3y_k','gfs/vc3vc3x2_k', $
		'gfs/vc3vc3w_k','gfs/vc3vc3v_k','gfs/vc3vc3u_k']
        ddir= [2,14,14,14,14,14,14]
        re_avgs3=     fload_avg_re(fruns,ddir, reerr=reerr3)
	re_avgs3= re_avgs3 * 1.3 + 0.3
        es_avgs3=     fload_avg_es(fruns,ddir, eserr=eserr3)
        ;vmax_avgs3=   fload_avg_vmax(fruns,ddir)
        ;sig_avgs3=    fload_avg_sig(fruns,ddir)
        ;vs_avgs3=     fload_avg_vs(fruns,ddir)
        vs_avgs3=     fload_avg_vsst(fruns,ddir, vserr=vserr3)
        TTs_avgs3=    fload_Ts(fruns)
        deltas_avgs3= fload_deltas(fruns)

endif



;------------------------------------------------------------------------------------

; sizes
; -----------
do_sizes= 0
;do_sizes= 1
if do_sizes eq 1 then begin
	; various sizes (sigma's)
	;xvar= [50,70,102.0,140.0,271.0,388.0,563.0] ; original disk v_200
	;xvar= [62.0,105.0,146.0,178.0,241.0,367.0,552.0] ; remnant sigma
	;xmax= 600.0 & xmin= 30.0
	xmax= 500.0 & xmin= 30.0
	xaxistitle='!7r!3 (km s!E-1!N)'


        ; vsig files
	cmt1= 'orbit e'
	;fruns= ['vc1vc1','vc2vc2','vc3vc3','vc4vc4a','vc5vc5a','vc6vc6a']
	;ddir= [4,4,7,4,4,4]
        ;fruns= ['ds/d0e','ds/d1e','ds/d2e','vc3vc3e','ds/d4e','ds/d5e','ds/d6e']
        ;ddir= [9,9,9,0,9,9,9]
        fruns= ['ds/d0e2_q','ds/d1e2_q','ds/d2e2_q','vc3vc3e','ds/d4e2_q','ds/d5e2_q','ds/d6e2_q']
        ddir= [15,15,15,0,15,15,15]

        re_avgs=     fload_avg_re(fruns,ddir, reerr=reerr)
	re_avgs= re_avgs * 1.3 + 0.3
        es_avgs=     fload_avg_es(fruns,ddir, eserr=eserr)
        ;vmax_avgs=   fload_avg_vmax(fruns,ddir)
        sig_avgs=    fload_avg_sig(fruns,ddir)
	xvar= sig_avgs
        ;vs_avgs=     fload_avg_vs(fruns,ddir, vserr=vserr)
        vs_avgs=     fload_avg_vsst(fruns,ddir, vserr=vserr)
        TTs_avgs=    fload_Ts(fruns)
        deltas_avgs= fload_deltas(fruns)
	fruns= '/raid4/tcox/' + fruns
	a4diva_avgs= 100.0*fload_a4diva(fruns, a4err=a4err)
	a4err= 100.0*a4err


	; second group
	secondgrp= 1
	cmt2= 'orbit h'
        ;fruns= ['ds/d0h','ds/d1h','ds/d2h','vc3vc3h','ds/d4h','ds/d5h','ds/d6h']
        ;ddir= [9,9,9,0,9,9,9]
        fruns= ['ds/d0h2_q','ds/d1h2_q','ds/d2h2_q','vc3vc3h','ds/d4h2_q','ds/d5h2_q','ds/d6h2_q']
        ddir= [15,15,15,0,15,15,15]
        re_avgs2=     fload_avg_re(fruns,ddir, reerr=reerr2)
	re_avgs2= re_avgs2 * 1.3 + 0.3
        es_avgs2=     fload_avg_es(fruns,ddir, eserr=eserr2)
        ;vmax_avgs2=   fload_avg_vmax(fruns,ddir)
        sig_avgs2=    fload_avg_sig(fruns,ddir)
	xvar2= sig_avgs2
        ;vs_avgs2=     fload_avg_vs(fruns,ddir, vserr=vserr2)
        vs_avgs2=     fload_avg_vsst(fruns,ddir, vserr=vserr2)
        TTs_avgs2=    fload_Ts(fruns)
        deltas_avgs2= fload_deltas(fruns)
	fruns= '/raid4/tcox/' + fruns
	a4diva_avgs2= 100.0*fload_a4diva(fruns, a4err=a4err2)
	a4err2= 100.0*a4err2


        ; third group
        thirdgrp= 1
	cmt3= 'orbit k'
        ;fruns= ['ds/d0k','ds/d1k','ds/d2k','vc3vc3k','ds/d4k','ds/d5k','ds/d6k']
        ;ddir= [9,9,9,0,9,9,9]
        fruns= ['ds/d0k2_q','ds/d1k2_q','ds/d2k2_q','vc3vc3k','ds/d4k2_q','ds/d5k2_q','ds/d6k2_q']
        ddir= [15,15,15,0,15,15,15]
        re_avgs3=     fload_avg_re(fruns,ddir, reerr=reerr3)
	re_avgs3= re_avgs3 * 1.3 + 0.3
        es_avgs3=     fload_avg_es(fruns,ddir, eserr=eserr3)
        ;vmax_avgs3=   fload_avg_vmax(fruns,ddir)
        sig_avgs3=    fload_avg_sig(fruns,ddir)
	xvar3= sig_avgs3
        ;vs_avgs3=     fload_avg_vs(fruns,ddir, vserr=vserr3)
        vs_avgs3=     fload_avg_vsst(fruns,ddir, vserr=vserr3)
        TTs_avgs3=    fload_Ts(fruns)
        deltas_avgs3= fload_deltas(fruns)
	fruns= '/raid4/tcox/' + fruns
	a4diva_avgs3= 100.0*fload_a4diva(fruns, a4err=a4err3)
	a4err3= 100.0*a4err3


        ; forth group
        fourthgrp= 1
        cmt4= 'orbit f'
        ;fruns= ['ds/d0f','ds/d1f','ds/d2f','vc3vc3f','ds/d4f','ds/d5f','ds/d6f']
        ;ddir= [9,9,9,0,9,9,9]
        fruns= ['ds/d0f2_q','ds/d1f2_q','ds/d2f2_q','vc3vc3f','ds/d4f2_q','ds/d5f2_q','ds/d6f2_q']
        ddir= [15,15,15,0,15,15,15]
        re_avgs4=     fload_avg_re(fruns,ddir, reerr=reerr4)
        re_avgs4= re_avgs4 * 1.3 + 0.3
        es_avgs4=     fload_avg_es(fruns,ddir, eserr=eserr4)
        ;vmax_avgs4=   fload_avg_vmax(fruns,ddir)
        sig_avgs4=    fload_avg_sig(fruns,ddir)
        xvar4= sig_avgs4
        ;vs_avgs4=     fload_avg_vs(fruns,ddir, vserr=vserr4)
        vs_avgs4=     fload_avg_vsst(fruns,ddir, vserr=vserr4)
        TTs_avgs4=    fload_Ts(fruns)
        deltas_avgs4= fload_deltas(fruns)
        fruns= '/raid4/tcox/' + fruns
        a4diva_avgs4= 100.0*fload_a4diva(fruns, a4err=a4err4)
	a4err4= 100.0*a4err4


	print, "sig_avgs"
	print, "--------"
	print, "e"
	print, sig_avgs
	print, "h"
	print, sig_avgs2
	print, "k"
	print, sig_avgs3
	print, "f"
	print, sig_avgs4

endif



;------------------------------------------------------------------------------------

; orbits
; -----------
do_orbits= 0
;do_orbits= 1
if do_orbits eq 1 then begin
        ; various sizes (sigma's)
        xvar= [2.5,5,10,15,20,30,40]
	xvar= xvar/0.7
        xmax= 60.0 & xmin= 0.0
        ;xaxistitle='R!Dperi!N (kpc h!E-1!N)'
        xaxistitle='!6R!Dperi!N (kpc)'


        ; vsig files
	cmt1= 'orbit e'
        fruns= ['ds/d3e1','vc3vc3e','ds/d3e2','ds/d3e3','ds/d3e4','ds/d3e6','ds/d3e5']
        ddir= [9,0,9,9,9,9,9]

        re_avgs=     fload_avg_re(fruns,ddir, reerr=reerr)
	re_avgs= re_avgs * 1.3 + 0.3
        es_avgs=     fload_avg_es(fruns,ddir, eserr=eserr)
        ;vmax_avgs=   fload_avg_vmax(fruns,ddir)
        ;sig_avgs=    fload_avg_sig(fruns,ddir)
        ;vs_avgs=     fload_avg_vs(fruns,ddir, vserr=vserr)
        vs_avgs=     fload_avg_vsst(fruns,ddir, vserr=vserr)
        TTs_avgs=    fload_Ts(fruns)
        deltas_avgs= fload_deltas(fruns)


	; second group
	secondgrp= 1
	xvar2= xvar
	cmt2= 'orbit h'
        fruns= ['ds/d3h1','vc3vc3h','ds/d3h2','ds/d3h3','ds/d3h4','ds/d3h6','ds/d3h5']
        ddir= [9,0,9,9,9,9,9]
        re_avgs2=     fload_avg_re(fruns,ddir, reerr=reerr2)
	re_avgs2= re_avgs2 * 1.3 + 0.3
        es_avgs2=     fload_avg_es(fruns,ddir, eserr=eserr2)
        ;vmax_avgs2=   fload_avg_vmax(fruns,ddir)
        ;sig_avgs2=    fload_avg_sig(fruns,ddir)
        ;vs_avgs2=     fload_avg_vs(fruns,ddir, vserr=vserr2)
        vs_avgs2=     fload_avg_vsst(fruns,ddir, vserr=vserr2)
        TTs_avgs2=    fload_Ts(fruns)
        deltas_avgs2= fload_deltas(fruns)


        ; third group
        thirdgrp= 1
	xvar3= xvar
	cmt3= 'orbit k'
        fruns= ['ds/d3k1','vc3vc3k','ds/d3k2','ds/d3k3','ds/d3k4','ds/d3k6','ds/d3k5']
        ;fruns= ['ds/d3k1','vc3vc3k','ds/d3k2','ds/d3k4','ds/d3k4','ds/d3k6','ds/d3k5']
        ddir= [9,0,9,9,9,9,9]
        re_avgs3=     fload_avg_re(fruns,ddir, reerr=reerr3)
	re_avgs3= re_avgs3 * 1.3 + 0.3
        es_avgs3=     fload_avg_es(fruns,ddir, eserr=eserr3)
        ;vmax_avgs3=   fload_avg_vmax(fruns,ddir)
        ;sig_avgs3=    fload_avg_sig(fruns,ddir)
        ;vs_avgs3=     fload_avg_vs(fruns,ddir, vserr=vserr3)
        vs_avgs3=     fload_avg_vsst(fruns,ddir, vserr=vserr3)
        TTs_avgs3=    fload_Ts(fruns)
        deltas_avgs3= fload_deltas(fruns)


endif




;------------------------------------------------------------------------------------


; Print thie mess up
; -------------------
filename='something.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable= 1, newxsize=12, newysize=22
;setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12, newysize=22
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12, newysize=30

x0= 0.15 & x1= 0.99
y0= 0.05 & ysize=0.156667
y1= y0+ysize
y2= y0+ysize+ysize
y3= y0+ysize+ysize+ysize
y4= y0+ysize+ysize+ysize+ysize
y5= y0+ysize+ysize+ysize+ysize+ysize
y6= y0+ysize+ysize+ysize+ysize+ysize+ysize


f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill



; zeroth: Re
; -------------
yaxistitle='!6<!17a!6>'
if do_sizes eq 1 then ymax = 15.8
if do_gf eq 1 then ymax = 6.9
if do_orbits eq 1 then ymax = 5.2
ymin= 0.8
if do_orbits eq 1 then ymin = 2.3
!p.position= [x0, y5, x1, y6]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, xvar, re_avgs, psym=-2, color=150, thick= 3.0, symsize= 1.5
if secondgrp eq 1 then oplot, xvar2, re_avgs2, psym=-8, color=50, thick= 3.0, symsize= 1.5
if thirdgrp eq 1 then oplot, xvar3, re_avgs3, psym=-5, color=100, thick= 3.0, symsize= 1.5
if fourthgrp eq 1 then oplot, xvar4, re_avgs4, psym=-7, color=200, thick= 3.0, symsize= 1.5

oploterror, xvar, re_avgs, reerr, psym=-3, errcolor=150, color=150, errthick= 3.0, thick= 3.0
if secondgrp eq 1 then oploterror, xvar2, re_avgs2, reerr2, psym=-3, errcolor=50, color=50, errthick= 3.0, thick=3.0
if thirdgrp eq 1 then oploterror, xvar3, re_avgs3, reerr3, psym=-3, errcolor=100, color=100, errthick= 3.0, thick=3.0
if fourthgrp eq 1 then oploterror, xvar4, re_avgs4, reerr4, psym=-3, errcolor=200, color=200, errthick= 3.0, thick=3.0


; gas fraction
if do_gf eq 1 then begin
    oplot, [0.65,0.75], [5.8,5.8], psym=-2, color=150, thick= 3.0, symsize= 1.5
    if secondgrp eq 1 then oplot, [0.65,0.75], [4.75,4.75], psym=-8, color=50, thick= 3.0, symsize= 1.5
    if thirdgrp eq 1 then oplot, [0.65,0.75], [3.73,3.73], psym=-5, color=100, thick= 3.0, symsize= 1.5
endif

; orbits
if do_orbits eq 1 then begin
    oplot, [6.5,10.7], [4.65,4.65], psym=-2, color=150, thick= 3.0, symsize= 1.5
    if secondgrp eq 1 then oplot, [6.5,10.7], [4.15,4.15], psym=-8, color=50, thick= 3.0, symsize= 1.5
    if thirdgrp eq 1 then oplot, [6.5,10.7], [3.7,3.7], psym=-5, color=100, thick= 3.0, symsize= 1.5
endif

; size
if do_sizes eq 1 then begin
    oplot, [55.0,90.0], [13.0,13.0], psym=-2, color=150, thick= 3.0, symsize= 1.5
    if secondgrp eq 1 then oplot, [55.0,90.0], [11.2,11.2], psym=-8, color=50, thick= 3.0, symsize= 1.5
    if thirdgrp eq 1 then oplot, [55.0,90.0], [9.4,9.4], psym=-5, color=100, thick= 3.0, symsize= 1.5
    if fourthgrp eq 1 then oplot, [55.0,90.0], [7.6,7.6], psym=-7, color=200, thick= 3.0, symsize= 1.5
endif




; zeroth: Ecc
; -------------
yaxistitle='!6<!7e!6>'
ymax = 0.9
ymin= 0.0
!p.position= [x0, y4, x1, y5]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, xvar, es_avgs, psym=-2, color=150, thick= 3.0, symsize= 1.5
if secondgrp eq 1 then oplot, xvar2, es_avgs2, psym=-8, color=50, thick= 3.0, symsize= 1.5
if thirdgrp eq 1 then oplot, xvar3, es_avgs3, psym=-5, color=100, thick= 3.0, symsize= 1.5
if fourthgrp eq 1 then oplot, xvar4, es_avgs4, psym=-7, color=200, thick= 3.0, symsize= 1.5

oploterror, xvar, es_avgs, eserr, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0
if secondgrp eq 1 then oploterror, xvar2, es_avgs2, eserr2, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0
if thirdgrp eq 1 then oploterror, xvar3, es_avgs3, eserr3, psym=-3, errcolor=100, color=100, thick= 3.0, errthick= 3.0
if fourthgrp eq 1 then oploterror, xvar4, es_avgs4, eserr4, psym=-3, errcolor=200, color=200, thick= 3.0, errthick= 3.0


; gas fraction
;oplot, [0.65,0.75], [0.75,0.75], psym=-2, color=150, thick= 3.0, symsize= 1.5
;if secondgrp eq 1 then oplot, [0.65,0.75], [0.65,0.65], psym=-8, color=50, thick= 3.0, symsize= 1.5
;if thirdgrp eq 1 then oplot, [0.65,0.75], [0.53,0.53], psym=-5, color=100, thick= 3.0, symsize= 1.5

; orbits
;oplot, [30.0,32.5], [0.75,0.75], psym=-2, color=150, thick= 3.0, symsize= 1.5
;if secondgrp eq 1 then oplot, [30.0,32.5], [0.65,0.65], psym=-8, color=50, thick= 3.0, symsize= 1.5
;if thirdgrp eq 1 then oplot, [30.0,32.5], [0.53,0.53], psym=-5, color=100, thick= 3.0, symsize= 1.5




; first: V_maj/sigma
; -------------------
;yaxistitle='<V!Dmaj!N/!7r!3>'
yaxistitle='!6<(V!Dmaj!N/!7r!6)!E*!N>'
ymax = 1.79
ymin= 0.0
!p.position= [x0, y3, x1, y4]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, xvar, vs_avgs, psym=-2, color=150, thick= 3.0, symsize= 1.5
if secondgrp eq 1 then oplot, xvar2, vs_avgs2, psym=-8, color=50, thick= 3.0, symsize= 1.5
if thirdgrp eq 1 then oplot, xvar3, vs_avgs3, psym=-5, color=100, thick= 3.0, symsize= 1.5
if fourthgrp eq 1 then oplot, xvar4, vs_avgs4, psym=-7, color=200, thick= 3.0, symsize= 1.5

oploterror, xvar, vs_avgs, vserr, psym=-3, errcolor=150, color=150, thick= 3.0, errthick=3.0
if secondgrp eq 1 then oploterror, xvar2, vs_avgs2, vserr2, psym=-3, errcolor=50, color=50, thick= 3.0, errthick=3.0
if thirdgrp eq 1 then oploterror, xvar3, vs_avgs3, vserr3, psym=-3, errcolor=100, color=100, thick= 3.0, errthick=3.0
if fourthgrp eq 1 then oploterror, xvar4, vs_avgs4, vserr4, psym=-3, errcolor=200, color=200, thick= 3.0, errthick=3.0


; second: Triaxiality
; ---------------------
yaxistitle='!6Triaxiality, T'
;ymax = 0.9
ymax = 1.1
ymin= 0.0
!p.position= [x0, y2, x1, y3]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, xvar, TTs_avgs, psym=-2, color=150, thick= 3.0, symsize= 1.5
if secondgrp eq 1 then oplot, xvar2, TTs_avgs2, psym=-8, color=50, thick= 3.0, symsize= 1.5
if thirdgrp eq 1 then oplot, xvar3, TTs_avgs3, psym=-5, color=100, thick= 3.0, symsize= 1.5
if fourthgrp eq 1 then oplot, xvar4, TTs_avgs4, psym=-7, color=200, thick= 3.0, symsize= 1.5



; second: delta
; ----------------
yaxistitle='!6Anisotropy, !7d!6'
ymax = 0.74
ymin= 0.0
!p.position= [x0, y1, x1, y2]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        ;xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, xvar, deltas_avgs, psym=-2, color=150, thick= 3.0, symsize= 1.5
if secondgrp eq 1 then oplot, xvar2, deltas_avgs2, psym=-8, color=50, thick= 3.0, symsize= 1.5
if thirdgrp eq 1 then oplot, xvar3, deltas_avgs3, psym=-5, color=100, thick= 3.0, symsize= 1.5
if fourthgrp eq 1 then oplot, xvar4, deltas_avgs4, psym=-7, color=200, thick= 3.0, symsize= 1.5





; last: a4/a
; ----------------
yaxistitle='!6100 a!D4!N/a'
ymax = 2.2
ymin= -2.0
!p.position= [x0, y0, x1, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase

oplot, xvar, a4diva_avgs, psym=-2, color=150, thick= 3.0, symsize= 1.5
if secondgrp eq 1 then oplot, xvar2, a4diva_avgs2, psym=-8, color=50, thick= 3.0, symsize= 1.5
if thirdgrp eq 1 then oplot, xvar3, a4diva_avgs3, psym=-5, color=100, thick= 3.0, symsize= 1.5
if fourthgrp eq 1 then oplot, xvar4, a4diva_avgs4, psym=-7, color=200, thick= 3.0, symsize= 1.5

oploterror, xvar, a4diva_avgs, a4err, psym=-3, errcolor=150, color=150, thick= 3.0, errthick=3.0
if secondgrp eq 1 then oploterror, xvar2, a4diva_avgs2, a4err2, psym=-3, errcolor=50, color=50, thick= 3.0, errthick=3.0
if thirdgrp eq 1 then oploterror, xvar3, a4diva_avgs3, a4err3, psym=-3, errcolor=100, color=100, thick= 3.0, errthick=3.0
if fourthgrp eq 1 then oploterror, xvar4, a4diva_avgs4, a4err4, psym=-3, errcolor=200, color=200, thick= 3.0, errthick=3.0



x=[xmin,xmax]
y=0.0 + 0*x
oplot, x, y, psym=-3, color=0, thick=2, linestyle= 1


; -----------------------------------------------------------------------------


if do_gf eq 1 then begin
    xyouts, 0.78, 0.95, cmt1, /normal, charthick=3.0, size=1.5, color=150
    if secondgrp eq 1 then xyouts, 0.78, 0.92, cmt2, /normal, charthick=3.0, size=1.5, color=50
    if thirdgrp eq 1 then xyouts, 0.78, 0.89, cmt3, /normal, charthick=3.0, size=1.5, color=100
endif

if do_orbits eq 1 then begin
    xyouts, 0.32, 0.95, cmt1, /normal, charthick=3.0, size=1.5, color=150
    if secondgrp eq 1 then xyouts, 0.32, 0.92, cmt2, /normal, charthick=3.0, size=1.5, color=50
    if thirdgrp eq 1 then xyouts, 0.32, 0.89, cmt3, /normal, charthick=3.0, size=1.5, color=100
endif

if do_sizes eq 1 then begin
    xyouts, 0.30, 0.9575, cmt1, /normal, charthick=3.0, size=1.5, color=150
    if secondgrp eq 1 then xyouts, 0.30, 0.9375, cmt2, /normal, charthick=3.0, size=1.5, color=50
    if thirdgrp eq 1 then xyouts, 0.30, 0.9175, cmt3, /normal, charthick=3.0, size=1.5, color=100
    if fourthgrp eq 1 then xyouts, 0.30, 0.8975, cmt3, /normal, charthick=3.0, size=1.5, color=200
endif



; -------------
;  Done
; -------------

device, /close



end







;================================================================================














