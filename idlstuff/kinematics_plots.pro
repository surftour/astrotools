; =================================================================
; =================================================================
;
; Read kinematics information (from files) and make plots.
;
; =================================================================
; =================================================================
;
;
;
;   1. first section of scripts read data (.txt) files
;   2. second set are floads
;   3. lastly come the plotting scripts
;       * a4_hist = histogram of a4's
;       * kb_map = (v/sig)* & vmin/vtot versus a4/a
;       * db_map = a4/a versus ellipticity
;       * vsig_map = (v/sig) versus ellipticity
;       * vsig_map_5x3 = 5x3 version of vsig_map
;











; =================================================================
;
;
;
;    read files in directories
;
;
; =================================================================



pro read_diskyboxy, frun, aa, bb, ellip, a4diva, $
			theta=theta, phi=phi

dbfile= frun+'/diskyboxy.txt'

spawn, "wc "+dbfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then db_data= fltarr(6,lines)

openr, 1, dbfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, db_data
close, 1


theta= db_data[0,*]
phi= db_data[1,*]
aa= db_data[2,*]
bb= db_data[3,*]
ellip= db_data[4,*]
a4diva= db_data[5,*]

end



; =================================================================


pro read_kinematics, frun, Re, Vrot_maj, Vrotavg_maj, $
			Vrot_min, Vrotavg_min, $
			Sig0, aa, bb, ellip, a4diva, pa, $
			theta=theta, phi=phi, $
			small_fov=small_fov, $
			med_fov=med_fov, $
			large_fov=large_fov

dbfile= frun+'/kinematics.txt'
if keyword_set(small_fov) then dbfile= frun+'/kinematics_small_fov.txt'
if keyword_set(med_fov) then dbfile= frun+'/kinematics_med_fov.txt'
if keyword_set(large_fov) then dbfile= frun+'/kinematics_large_fov.txt'

spawn, "wc "+dbfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then db_data= fltarr(13,lines)

print, "opening: ", dbfile
openr, 1, dbfile, ERROR= err

if (err NE 0) then begin
        err= 0 
        close, 1
        print, " "
        print, " error opening file: ", dbfile
        print, " "
        theta= [0,0]
        phi= [0,0]
        Re= [0,0]
        aa= [0,0]
        bb= [0,0]
        ellip= [0,0]
        pa= [0,0] 
        a4diva= [0,0]
        Sig0= [0,0]
	Vrot_maj= [0,0]
	Vrotavg_maj= [0,0]
	Vrot_min= [0,0]
	Vrotavg_min= [0,0]
        return
endif

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, db_data
close, 1


theta= transpose(db_data[0,*])
phi= transpose(db_data[1,*])
Re= transpose(db_data[2,*])
aa= transpose(db_data[3,*])
bb= transpose(db_data[4,*])
ellip= transpose(db_data[5,*])
pa= transpose(db_data[6,*])
a4diva= transpose(db_data[7,*])
Sig0= transpose(db_data[8,*])
Vrot_maj= transpose(db_data[9,*])
Vrotavg_maj= transpose(db_data[10,*])
Vrot_min= transpose(db_data[11,*])
Vrotavg_min= transpose(db_data[12,*])


; trap for any NaN's
found_bad= 0
idx=where(finite(ellip) eq 0)
if idx(0) ne -1 then begin & badidx= idx & found_bad= 1 & endif

if found_bad eq 1 then begin
	idx= where(finite(ellip) eq 1)
	theta= theta(idx)
	phi= phi(idx)
	Re= Re(idx)
	Vrot_maj= Vrot_maj(idx)
	Vrotavg_maj= Vrotavg_maj(idx)
	Vrot_min= Vrot_min(idx)
	Vrotavg_min= Vrotavg_min(idx)
	Sig0= Sig0(idx)
	aa= aa(idx)
	bb= bb(idx)
	ellip= ellip(idx)
	a4diva= a4diva(idx)
	pa= pa(idx)
        print, "Fixed ",n_elements(badidx), " kinematics measures when trapping for NaN values."
endif


print, "<Re>= ",mean(Re), "  +/-", sqrt(variance(Re))
print, "<Vrot_maj>= ",mean(Vrot_maj), "  +/-", sqrt(variance(Vrot_maj))
print, "<Vrotavg_maj>= ",mean(Vrotavg_maj), "  +/-", sqrt(variance(Vrotavg_maj))
print, "<Vrot_min>= ",mean(Vrot_min), "  +/-", sqrt(variance(Vrot_min))
print, "<Vrotavg_min>= ",mean(Vrotavg_min), "  +/-", sqrt(variance(Vrotavg_min))
print, "<Sig0>= ",mean(Sig0), "  +/-", sqrt(variance(Sig0))
print, "<aa>= ",mean(aa), "  +/-", sqrt(variance(aa))
print, "<bb>= ",mean(bb), "  +/-", sqrt(variance(bb))
print, "<ellip>= ",mean(ellip), "  +/-", sqrt(variance(ellip))
print, "<a4diva>= ",mean(a4diva), "  +/-", sqrt(variance(a4diva))

end





; =================================================================
;
;  Not currently used by any scripts - a relic of when suvendra
;  did the anlysis
;
;
;
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




; =================================================================



pro read_kinemetry, frun, angleidx, radius, pa, ellip, q, k1, k1_err, k5k1, k5k1_err

;dbfile= frun+'/kinematics.txt'
if angleidx ge 0 then begin
	anglelbl= strcompress(string(angleidx),/remove_all)
	dbfile= frun+'/'+anglelbl+'_lst_kinemetry.txt'
endif else begin
	dbfile= frun
endelse

spawn, "wc "+dbfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then db_data= fltarr(6,lines)
;if lines GT 0 then db_data= fltarr(8,lines)

print, "opening: ", dbfile
openr, 1, dbfile, ERROR= err

if (err NE 0) then begin
	err= 0
	close, 1
	print, " "
	print, " error opening file: ", dbfile
	print, " "
	radius= [0]
	pa= [0]
	ellip= [0]
	q= [0]
	k1= [0]
	k1_err= [0]
	k5k1= [0]
	k5k1_err= [0]
	return
endif

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, db_data
close, 1


radius= transpose(db_data[0,*])
pa= transpose(db_data[1,*])
ellip= transpose(db_data[2,*])
q= transpose(db_data[3,*])
k1= transpose(db_data[4,*])
k1_err= 0.1 * k1
;k1_err= transpose(db_data[5,*])
k5k1= transpose(db_data[5,*])
;k5k1= transpose(db_data[6,*])
k5k1_err= 0.1 * k5k1
;k5k1_err= transpose(db_data[7,*])


end





; =================================================================




pro read_isophotes, frun, angleidx, radius, sd, ellip, pa, a4diva

if angleidx ge 0 then begin
	anglelbl= strcompress(string(angleidx),/remove_all)
	dbfile= frun+'/'+anglelbl+'_lst_isophotes.txt'
endif else begin
	dbfile= frun
endelse

spawn, "wc "+dbfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then db_data= fltarr(5,lines)

print, "opening: ", dbfile
openr, 1, dbfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, db_data
close, 1


radius= transpose(db_data[0,*])
sd= transpose(db_data[1,*])
ellip= transpose(db_data[2,*])
pa= transpose(db_data[3,*])
a4diva= transpose(db_data[4,*])


end





; =================================================================




pro read_lambda, frun, angleidx, radius, lambda, sigma, v

if angleidx ge 0 then begin
	anglelbl= strcompress(string(angleidx),/remove_all)
	dbfile= frun+'/'+anglelbl+'_lst_lambda.txt'
endif else begin
	dbfile= frun
endelse

spawn, "wc "+dbfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then db_data= fltarr(4,lines)

print, "opening: ", dbfile
openr, 1, dbfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, db_data
close, 1


radius= transpose(db_data[0,*])
lambda= transpose(db_data[1,*])
sigma= transpose(db_data[2,*])
v= transpose(db_data[3,*])


end





; =================================================================




pro read_slit_info, frun, angleidx, radius, v, sigma, major=major, minor=minor

majormin= "major"
if keyword_set(major) then majormin="major"
if keyword_set(minor) then majormin="minor"

if angleidx ge 0 then begin
	anglelbl= strcompress(string(angleidx),/remove_all)
	dbfile= frun+'/'+anglelbl+'_lst_slit_'+majormin+'.txt'
endif else begin
	dbfile= frun
endelse

spawn, "wc "+dbfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then db_data= fltarr(3,lines)

print, "opening: ", dbfile
openr, 1, dbfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, db_data
close, 1


radius= transpose(db_data[0,*])
v= transpose(db_data[1,*])
sigma= transpose(db_data[2,*])


end




; =================================================================
;
;
;   Floads
;
;
; =================================================================



;   ellipticities
; --------------------------
function fload_ellips, fruns, all=all, eserr=eserr
   ns= n_elements(fruns)
   alles= fltarr(ns)
   eserr= fltarr(ns)
   es= [-1]
   for i=0,ns-1 do begin
        ;read_diskyboxy, fruns[i], aa, bb, ellip, a4diva
	;read_kinematics, fruns[i], Re, Vrot, Vrotavg, Sig0, aa, bb, ellip, a4diva
	;read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa
	read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa, /small_fov
        es= [es,transpose(ellip)]
        alles[i]= mean(ellip)
        eserr[i]= sqrt(variance(ellip))
   endfor
   es= es[1:n_elements(es)-1]
   if keyword_set(all) then return, es
   return, alles
end




;   semi-major axis
; --------------------------
function fload_aa, fruns, all=all, aaerr=aaerr
   ns= n_elements(fruns)
   allas= fltarr(ns)
   aaerr= fltarr(ns)
   as= [-1]
   for i=0,ns-1 do begin
        ;read_diskyboxy, fruns[i], aa, bb, ellip, a4diva
	;read_kinematics, fruns[i], Re, Vrot, Vrotavg, Sig0, aa, bb, ellip, a4diva
	;read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa
	read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa, /small_fov
	thisa= transpose(aa)
	thisb= transpose(bb)
	semimajor= thisa
	idx= where(thisb gt thisa)
	if idx(0) ne -1 then semimajor(idx)= thisb(idx)
        as= [as,semimajor]
        allas[i]= mean(semimajor)
        aaerr[i]= sqrt(variance(semimajor))
   endfor
   as= as[1:n_elements(as)-1]
   if keyword_set(all) then return, as
   return, allas
end




;   disky/boxy (a4/a)
; --------------------------
function fload_a4diva, fruns, all=all, a4err=a4err
   ns= n_elements(fruns)
   alla4= fltarr(ns)
   a4err= fltarr(ns)
   a4s= [-1]
   for i=0,ns-1 do begin
        ;read_diskyboxy, fruns[i], aa, bb, ellip, a4diva
	;read_kinematics, fruns[i], Re, Vrot, Vrotavg, Sig0, aa, bb, ellip, a4diva
	;read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa
	read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa, /small_fov
        a4s= [a4s,transpose(a4diva)]
        alla4[i]= mean(a4diva)
	print, fruns[i], "    <a4/a>= ", mean(100*a4diva), "   +/- ", sqrt(variance(100*a4diva))
        a4err[i]= sqrt(variance(a4diva))
   endfor
   a4s= a4s[1:n_elements(a4s)-1]
   if keyword_set(all) then return, a4s
   return, alla4
end



;   average (vmax/sigma)*
; --------------------------
function fload_avg_vsst, fruns, ddir, all=all, vserr=vserr
   ns= n_elements(fruns)
   allvsst= fltarr(ns)
   vserr= fltarr(ns)
   vsigst= [-1]
   for i=0,ns-1 do begin
        ;read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir[i]
	;read_kinematics, fruns[i], Re, Vrot, Vrotavg, Sig0, aa, bb, ellip, a4diva
	;read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa
	read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa, /small_fov
	vmax= Vrot
	sigma= Sig0
	ellipticity= ellip
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



;   average vmax
; --------------------------
function fload_avg_vmax, fruns, ddir, all=all
   ns= n_elements(fruns)
   allvmax= fltarr(ns)
   vs= [-1]
   for i=0,ns-1 do begin
        ;read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir[i]
	;read_kinematics, fruns[i], Re, Vrot, Vrotavg, Sig0, aa, bb, ellip, a4diva
	;read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa
	read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa, /small_fov
	vmax= Vrot
	sigma= Sig0
	ellipticity= ellip
        vs= [vs, transpose(vmax)]
        allvmax[i]= mean(vmax)
   endfor
   vs= vs[1:n_elements(vs)-1]
   if keyword_set(all) then return, vs
   return, allvmax
end



;   mu's
; ---------------------------
function fload_mus, fruns, cless=cless, all=all
   ddir_major= 0
   ;ddir_major= 13
   ddir_minor= 1
   if keyword_set(cless) then begin
        ddir_major= 2
        ddir_minor= 3
   endif
   ns= n_elements(fruns)
   allmus= fltarr(ns)
   mus= [-1]
   for i=0,ns-1 do begin
        ;read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir_major
	;read_kinematics, fruns[i], Re, Vrot, Vrotavg, Sig0, aa, bb, ellip, a4diva
	;read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa
	read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa, /small_fov
        ;this_mus= abs(atan(vmin/vmaj))*180.0/!PI
        ;print, "mean kinematic misalignment, <Phi>= ",mean(this_mus)
        this_mus= Vrot_min / sqrt(Vrot_min*Vrot_min + Vrot_maj*Vrot_maj)
        print, "                              <mu>= ",mean(this_mus)
        allmus[i]= mean(this_mus)
        mus= [mus,transpose(this_mus)]
   endfor
   mus= mus[1:n_elements(mus)-1]
   if keyword_set(all) then return, mus
   return, allmus
end



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













;===================================================================================
;===================================================================================
;
;
;
;   Plots
;
;   * a4_hist = histogram of a4's
;   * kb_map = (v/sig)* & vmin/vtot versus a4/a
;   * db_map = a4/a versus ellipticity
;   * vsig_map = (v/sig) versus ellipticity
;   * vsig_map_5x3 = 5x3 version of vsig_map
;
;
;
;
;===================================================================================
;===================================================================================





pro a4_hist, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "a4_hist, junk"
   print, "  "
   return
endif

filename='a4hist.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;--------------------------------------
;--------------------------------------

xaxistitle= '!6100 a4/a'
;xticklbls=[' ','-0.06',' ','-0.04',' ','-0.02',' ','0',' ','0.02',' ','0.04',' ','0.06',' ']
;xticknum= n_elements(xticklbls)-1
xmax= 7.00
xmin= -7.00


; --------------------


ymax = 1.2
ymin = 0.0


;----------------------------------------
; Generate plot
;----------------------------------------

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


; -----------------
;  std's
; -----------------

;fruns= ['vc3vc3h','vc3vc3b','vc3vc3c','vc3vc3d','vc3vc3e', $
;        'vc3vc3f','vc3vc3g','vc3vc3i','vc3vc3j','vc3vc3k', $
;        'vc3vc3l','vc3vc3m','vc3vc3n','vc3vc3o','vc3vc3p']
fruns= ['vc3vc3h','vc3vc3e']

for i=0, n_elements(fruns)-1 do fruns[i]= '/raid4/tcox/'+fruns[i]

alla4div=   100.0*fload_a4diva(fruns,/all)

print, "dissipational"
print, "-------------"
print, "a4/a   (max/min) = ", max(alla4div), min(alla4div)
print, "       (mean) = ", mean(alla4div)
print, "       (median) = ", median(alla4div)


; -----------------
;  collisionless
; -----------------

;fruns= ['cvc3vc3h','cvc3vc3b','cvc3vc3c','cvc3vc3d','cvc3vc3e', $
;        'cvc3vc3f','cvc3vc3g','cvc3vc3i','cvc3vc3j','cvc3vc3k', $
;        'cvc3vc3l','cvc3vc3m','cvc3vc3n','cvc3vc3o','cvc3vc3p']
fruns= ['d1h2_q','d1e2_q']

;for i=0, n_elements(fruns)-1 do fruns[i]= '/raid4/tcox/collisionless/'+fruns[i]
for i=0, n_elements(fruns)-1 do fruns[i]= '/raid4/tcox/ds/'+fruns[i]

allca4div=   100.0*fload_a4diva(fruns,/all)

print, "dissipationless"
print, "-------------"
print, "a4/a   (max/min) = ", max(allca4div), min(allca4div)
print, "       (mean) = ", mean(allca4div)
print, "       (median) = ", median(allca4div)


; --------------------



levels= 30.0


; two arrays to histogram

; a4/a
alltts= alla4div
allctts= allca4div


;xyouts, 0.62, 0.88, '!640% gas', /normal, charthick=3.0, size=1.5, color=150
;xyouts, 0.62, 0.82, '!6dissipationless', /normal, charthick=3.0, size=1.5, color=50
xyouts, 0.62, 0.88, '!6160', /normal, charthick=3.0, size=1.5, color=150
xyouts, 0.62, 0.82, '!680', /normal, charthick=3.0, size=1.5, color=50


idx=where(alltts gt 0.0)
print, '40: percent above 0.0', 100.0*n_elements(idx)/n_elements(alltts)
idx=where(allctts gt 0.0)
print, 'c: percent above 0.0', 100.0*n_elements(idx)/n_elements(allctts)

; --------------------

step= (xmax-xmin)/(levels)
bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
;bins=[bins,xmax]
bins=bins+(step*0.5)


; std histogram
hist_alltts= histogram(alltts, binsize=step, max=xmax, min=xmin)


; collisionless histogram
hist_allctts= histogram(allctts, binsize=step, max=xmax, min=xmin)

; make sure it extends to 0
;bins= [xmin,bins, bins(levels-1)]
bins= [xmin,bins,xmax]
hist_alltts= [hist_alltts(0),hist_alltts,hist_alltts(levels-1)]
hist_allctts= [hist_allctts(0),hist_allctts,hist_allctts(levels-1)]


normalization= float(max([hist_alltts,hist_allctts]))
print, "histogram normalization= ", normalization
hist_alltts= hist_alltts/normalization
hist_alltts_area = total(hist_alltts*step)
hist_allctts= hist_allctts/normalization

oplot, bins, hist_alltts, psym=10, color=150, thick=4.0
;oplot, bins, hist_allctts, psym=10, color=50, thick=4.0
oplot, bins, hist_allctts, psym=10, color=50, thick=8.0


; fill in histogram
; ------------------
nbins= bins+(step*0.5)          ; make x coord
;nbins[0]= 0.0
nbins[0]= xmin
nbins=[nbins,nbins]
nbins= nbins(sort(nbins))
;nbins= nbins[nbins, xmax]

ntts= fltarr(2.*levels + 2)
nctts= fltarr(2.*levels + 2)
for i=1,levels do begin
   ntts[2.*i-1]= hist_alltts[i]
   ntts[2.*i]= hist_alltts[i]
   nctts[2.*i-1]= hist_allctts[i]
   nctts[2.*i]= hist_allctts[i]
endfor

; collisionless
;
;polyfill, nbins, nctts, /data, color= 50, /line_fill, linestyle=0, $
;                               thick=3.0
;polyfill, nbins, nctts, /data, color= 50, /fill, linestyle=0, $
;                               thick=3.0
;polyfill, nbins, nctts, /data, color= 220, /fill, linestyle=0, $
;                               thick=3.0

; 40% gas
;polyfill, nbins, ntts, /data, color= 150, /line_fill, linestyle=0, $
;                               thick=3.0, orientation=90.0
polyfill, nbins, ntts, /data, color= 150, /line_fill, linestyle=0, $
                                thick=3.0, orientation=45.0



; ------------------------------------------------------------------


        ; --------------
        ; read bbf '92
        bbffile= '/home/tcox/RotvAnisSupport/bbf92_spheroids.txt'
        read_bbf_data, bbffile, mb, vsigst, a4diva

        idx= where(a4diva lt 3.5)
        if idx(0) ne -1 then a4diva= a4diva(idx)

	print, "BBH '92   N= ", n_elements(a4diva)
        print, "a4/a   (max/min) = ", max(a4diva), min(a4diva)
	print, "       (mean) = ", mean(a4diva)
	print, "       (median) = ", median(a4diva)
        idx= where(a4diva gt 0.0) 
        print, "a4/a percent above 0.0= ", 100.0*n_elements(idx)/n_elements(a4diva)

        ; now, make histogram
        levels= 10.0
        step= (xmax-xmin)/(levels)
        bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
        bins=bins+(step*0.5)

        ; std histogram
        hist_alla4diva= histogram(a4diva, binsize=step, max=xmax, min=xmin)

        ; make sure it extends to 0
        bins= [xmin,bins,xmax]
        hist_alla4diva= [hist_alla4diva(0),hist_alla4diva,hist_alla4diva(levels-1)]

        normalization= hist_alltts_area / total(hist_alla4diva*step)
        hist_alla4diva= hist_alla4diva*normalization
        ;hist_alla4diva= hist_alla4diva/normalization

        oplot, bins, hist_alla4diva, psym=10, color=0, thick=12.0, linestyle= 0

        oplot, [1.6,2.5], [0.90,0.90], psym=-3, color= 0, thick=12.0, linestyle= 0
        xyouts, 0.70, 0.76, 'BBF (1992)', /normal, charthick=4.0, size=1.5, color=0





        ; --------------------------
        ; read bender current data
	read_bender_current, name, ellip, a4a, vmaj, vmin, sig, pK, dPA, M_t, SB_e, L_R, L_X

        idx= where(a4a lt 50.0)
        if idx(0) ne -1 then a4a= a4a(idx)

        print, "Bender Current Data N= ", n_elements(a4a)
        print, "a4/a   (max/min) = ", max(a4a), min(a4a)
        print, "       (mean) = ", mean(a4a)
        print, "       (median) = ", median(a4a)
        idx= where(a4a gt 0.0)
        print, "a4/a percent above 0.0= ", 100.0*n_elements(idx)/n_elements(a4a)

        ; now, make histogram
        levels= 10.0
        step= (xmax-xmin)/(levels)
        bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
        bins=bins+(step*0.5)

        ; std histogram
        hist_alla4diva= histogram(a4a, binsize=step, max=xmax, min=xmin)

        ; make sure it extends to 0
        bins= [xmin,bins,xmax]
        hist_alla4diva= [hist_alla4diva(0),hist_alla4diva,hist_alla4diva(levels-1)]

        normalization= hist_alltts_area / total(hist_alla4diva*step)
        hist_alla4diva= hist_alla4diva*normalization
        ;hist_alla4diva= hist_alla4diva/normalization

        oplot, bins, hist_alla4diva, psym=10, color=0, thick=12.0, linestyle= 1

        oplot, [1.6,2.5], [0.85,0.85], psym=-3, color= 0, thick=12.0, linestyle= 1
        xyouts, 0.70, 0.71, 'Bender (prv.)', /normal, charthick=4.0, size=1.5, color=0



; ------------------------------------------------------------------



device, /close


end





;====================================================================





;------------------------------------------------------------------------
;
;     Density Map of (v/sig)* and vmin/vtot versus a4/a
;     --------------------------------------------------
;
;
;------------------------------------------------------------------------

pro kb_map, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "kb_map, junk"
   print, "  "
   return 
endif

filename='kbmap.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3



;--------------------------------------
;--------------------------------------


   fruns_40gas= ['data/ds/vc3vc3b', $
           'data/ds/vc3vc3c', $
           'data/ds/vc3vc3d', $
           'data/ds/vc3vc3e', $
           'data/ds/vc3vc3f', $
           'data/ds/vc3vc3g', $
           'data/ds/vc3vc3h', $
           'data/ds/vc3vc3i', $
           'data/ds/vc3vc3j', $
           'data/ds/vc3vc3k', $
           'data/ds/vc3vc3l', $
           'data/ds/vc3vc3m', $
           'data/ds/vc3vc3n', $
           'data/ds/vc3vc3o', $
           'data/ds/vc3vc3p']


   fruns_dless= ['data/collisionless/cvc3vc3b', $
           'data/collisionless/cvc3vc3c', $
           'data/collisionless/cvc3vc3d', $
           'data/collisionless/cvc3vc3e', $
           'data/collisionless/cvc3vc3f', $
           'data/collisionless/cvc3vc3g', $
           'data/collisionless/cvc3vc3h', $
           'data/collisionless/cvc3vc3i', $
           'data/collisionless/cvc3vc3j', $
           'data/collisionless/cvc3vc3k', $
           'data/collisionless/cvc3vc3l', $
           'data/collisionless/cvc3vc3m', $
           'data/collisionless/cvc3vc3n', $
           'data/collisionless/cvc3vc3o', $
           'data/collisionless/cvc3vc3p']


  fruns_remerg_b= ['data/remergers/b2eb2e', $
           'data/remergers/b2eb2h', $
           'data/remergers/b2eb2k', $
           'data/remergers/b3eb3e', $
           ;'data/remergers/b3eb3f', $
           'data/remergers/b3eb3h', $
           'data/remergers/b4eb3e', $
           'data/remergers/b4eb4e', $
           ;'data/remergers/b4eb4k', $
           'data/remergers/b4fb4h']


  fruns_remerg_v= ['data/remergers/vc3fvc3f', $
           'data/remergers/vc3fvc3f_polar', $
           'data/remergers/vc3hvc3h_1', $
           'data/remergers/vc3hvc3h_2', $
           'data/remergers/vc3hvc3h_3', $
           ;'data/remergers/vc3hvc3h_4', $
           'data/remergers/vc3mvc3m', $
           'data/remergers/vc3mvc3m_polar', $
           'data/remergers/vc3vc3i_rem3', $
           'data/remergers/vc3vc3i_rem4']


  fruns_remerg_e= ['data/remergers/e2ee2e', $
           'data/remergers/e2ke2k', $
           'data/remergers/e3ee3e', $
           'data/remergers/e3ke3e', $
           'data/remergers/e3ke3k', $
           'data/remergers/e4ee4e', $
           'data/remergers/e4ke4k']


  fruns_remerg_z= ['data/remergers/z3fz3f', $
           'data/remergers/z3fz3f_polar', $
           'data/remergers/z3hz3h', $
           'data/remergers/z3hz3h_polar', $
           'data/remergers/z3mz3m', $
           'data/remergers/z3mz3m_polar']



;  now plot things up
; ---------------------------


x0= 0.15
x1= 0.98

y0= 0.12
y1= 0.55
y2= 0.98


; ---------------------------

do_40gas= 0
;do_40gas= 1
if do_40gas eq 1 then begin
   process_one_map, 'kb1', fruns_40gas, x0, y1, x1, y2, /showyaxis, /show_x_eq_0, /show_y_eq_0
   oplot_kb1_data, 1
   process_one_map, 'kb2', fruns_40gas, x0, y0, x1, y1, /showyaxis, /showxaxis, /show_x_eq_0
   oplot_kb2_data, 1

   xyouts, 0.72, 0.63, /normal, '40% gas', size=1.5, charthick= 3.0, color=0
   xyouts, 0.70, 0.58, /normal, '(all orient)', size=1.5, charthick= 3.0, color=0

   ; --------------------
   ;process_one_map, 'kb1', fruns_40gas, x0, y0, x1, y1, /showyaxis, /showxaxis, /docontour
   ;xyouts, 0.25, 0.90, /normal, 'contours - 40% gas (all orient)', size=1.2, charthick= 3.0, color=0
endif


;do_dless= 0
do_dless= 1
if do_dless eq 1 then begin
   process_one_map, 'kb1', fruns_dless, x0, y1, x1, y2, /showyaxis, /show_x_eq_0, /show_y_eq_0
   oplot_kb1_data, 1
   process_one_map, 'kb2', fruns_dless, x0, y0, x1, y1, /showyaxis, /showxaxis, /show_x_eq_0
   oplot_kb2_data, 1
   xyouts, 0.62, 0.63, /normal, 'dissipationless', size=1.5, charthick= 3.0, color=0
   xyouts, 0.72, 0.58, /normal, '(all orient)', size=1.5, charthick= 3.0, color=0
endif


do_remerg= 0
;do_remerg= 1
if do_remerg eq 1 then begin
   ;fruns_remerg= [fruns_remerg_e(*),fruns_remerg_z(*)]
   ;fruns_remerg= fruns_remerg_b & msg='remergers (b)'
   ;fruns_remerg= fruns_remerg_v & msg='remergers (v)'
   fruns_remerg= fruns_remerg_e & msg='remergers (e)'
   ;fruns_remerg= fruns_remerg_z & msg='remergers (z)'
   process_one_map, 'kb1', fruns_remerg, x0, y1, x1, y2, /showyaxis, /show_x_eq_0, /show_y_eq_0
   oplot_kb1_data, 1
   process_one_map, 'kb2', fruns_remerg, x0, y0, x1, y1, /showyaxis, /showxaxis, /show_x_eq_0
   oplot_kb2_data, 1
   xyouts, 0.65, 0.58, /normal, msg, size=1.2, charthick= 3.0, color=0
endif



; --------------------

device, /close


end






;====================================================================












;------------------------------------------------------------------------
;
;     Density Map of a4/a vs. Ellipticity
;     ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro db_map, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "db_map, junk"
   print, "  "
   return
endif

filename='db.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3



;--------------------------------------
;--------------------------------------


x0= 0.18
x1= 0.98

y0= 0.15
y1= 0.98


;--------------------------------------
;--------------------------------------


   ;fruns_40gas= ['data/ds/vc3vc3b', $
   fruns_40gas= ['data/ds/vc3vc3c', $
           ;'data/ds/vc3vc3c', $
           'data/ds/vc3vc3d', $
           ;'data/ds/vc3vc3e', $
           ;'data/ds/vc3vc3f', $
           'data/ds/vc3vc3g', $
           'data/ds/vc3vc3h', $
           'data/ds/vc3vc3i', $
           'data/ds/vc3vc3j', $
           'data/ds/vc3vc3k', $
           'data/ds/vc3vc3l', $
           'data/ds/vc3vc3m', $
           'data/ds/vc3vc3n', $
           'data/ds/vc3vc3o', $
           'data/ds/vc3vc3p'] 


   fruns_dless= ['data/collisionless/cvc3vc3b', $
	   ;'data/collisionless/cvc3vc3c', $
	   'data/collisionless/cvc3vc3d', $
	   'data/collisionless/cvc3vc3e', $
	   'data/collisionless/cvc3vc3f', $
	   'data/collisionless/cvc3vc3g', $
	   'data/collisionless/cvc3vc3h', $
	   'data/collisionless/cvc3vc3i', $
	   'data/collisionless/cvc3vc3j', $
	   'data/collisionless/cvc3vc3k', $
	   'data/collisionless/cvc3vc3l', $
	   'data/collisionless/cvc3vc3m', $
	   'data/collisionless/cvc3vc3n', $
	   ;'data/collisionless/cvc3vc3o', $
	   'data/collisionless/cvc3vc3p']



  ;fruns_remerg_b= ['data/remergers/b2eb2e', $
           ;'data/remergers/b2eb2h', $
           ;'data/remergers/b2eb2k', $
           ;'data/remergers/b3eb3e', $
           ;'data/remergers/b3eb3f', $
           ;'data/remergers/b3eb3h', $
           ;'data/remergers/b4eb3e', $
           ;'data/remergers/b4eb4e', $
           ;'data/remergers/b4eb4k', $
           ;'data/remergers/b4fb4h']
  fruns_remerg_b= ['data/remergers/b3eb3f', $
           ;'data/remergers/b3eb3h', $
           'data/remergers/b4eb3e', $
           'data/remergers/b4eb4e', $
           'data/remergers/b4eb4k', $
           'data/remergers/b4fb4h']



  fruns_remerg_v= ['data/remergers/vc3fvc3f', $
           'data/remergers/vc3fvc3f_polar', $
           'data/remergers/vc3hvc3h_1', $
           'data/remergers/vc3hvc3h_2', $
           'data/remergers/vc3hvc3h_3', $
           'data/remergers/vc3hvc3h_4', $
           'data/remergers/vc3mvc3m', $
           'data/remergers/vc3mvc3m_polar', $
           'data/remergers/vc3vc3i_rem3', $
           'data/remergers/vc3vc3i_rem4']


  fruns_remerg_e= ['data/remergers/e2ee2e', $
           'data/remergers/e2ke2k', $
           'data/remergers/e3ee3e', $
           'data/remergers/e3ke3e', $
           'data/remergers/e3ke3k', $
           'data/remergers/e4ee4e', $
           'data/remergers/e4ke4k']


  fruns_remerg_z= ['data/remergers/z3fz3f', $
           'data/remergers/z3fz3f_polar', $
           'data/remergers/z3hz3h', $
           'data/remergers/z3hz3h_polar', $
           'data/remergers/z3mz3m', $
           'data/remergers/z3mz3m_polar']




;  now plot things up
; ---------------------------

do_one= 0
;do_one= 1
if do_one eq 1 then begin
   ;fruns= "data/ds/vc3vc3e_2" & color= 150
   ;fruns= "data/ds/vc3vc3d" & color= 150
   ;fruns= "data/ds/vc3vc3f" & color= 150
   fruns= "data/ds/vc3vc3i" & color= 150
   process_one_map, 'db', fruns, x0, y0, x1, y1, /showyaxis, /showxaxis, /dopoints, color= color
   xyouts, 0.26, 0.88, /normal, fruns, size=1.5, charthick= 3.0, color=0
endif


do_40gas= 0
;do_40gas= 1
if do_40gas eq 1 then begin
   process_one_map, 'db', fruns_40gas, x0, y0, x1, y1, /showyaxis, /showxaxis
   xyouts, 0.26, 0.88, /normal, '40% gas', size=1.5, charthick= 3.0, color=0
   xyouts, 0.24, 0.83, /normal, '(all orient)', size=1.5, charthick= 3.0, color=0
   ;process_one_map, 'db', fruns_40gas, x0, y0, x1, y1, /showyaxis, /showxaxis, /docontour
   ;xyouts, 0.25, 0.90, /normal, 'contours - 40% gas (all orient)', size=1.2, charthick= 3.0, color=0
endif


do_dless= 0
;do_dless= 1
if do_dless eq 1 then begin
   process_one_map, 'db', fruns_dless, x0, y0, x1, y1, /showyaxis, /showxaxis
   xyouts, 0.24, 0.88, /normal, 'dissipationless', size=1.5, charthick= 3.0, color=0
   xyouts, 0.26, 0.83, /normal, '(all orient)', size=1.5, charthick= 3.0, color=0
endif


;do_remerg= 0
do_remerg= 1
if do_remerg eq 1 then begin
   ;fruns_remerg= [fruns_remerg_e(*),fruns_remerg_z(*)]
   fruns_remerg= fruns_remerg_b & msg='remergers (b)'
   ;fruns_remerg= fruns_remerg_v & msg='remergers (v)'
   ;fruns_remerg= fruns_remerg_e & msg='remergers (e)'
   ;fruns_remerg= fruns_remerg_z & msg='remergers (z)'
   process_one_map, 'db', fruns_remerg, x0, y0, x1, y1, /showyaxis, /showxaxis
   ;xyouts, 0.25, 0.85, /normal, 'map - remergers', size=1.2, charthick= 3.0, color=0
   xyouts, 0.25, 0.85, /normal, msg, size=1.2, charthick= 3.0, color=0
endif


;do_40gas= 0
do_40gas= 1
if do_40gas eq 1 then begin
   process_one_map, 'db', fruns_40gas, x0, y0, x1, y1, /docontour
   xyouts, 0.25, 0.90, /normal, 'contours - 40% gas (all orient)', size=1.2, charthick= 3.0, color=0
endif



; --------------------

; line of isotropic rotator
;x=findgen(50)*(0.6/50.)
;y=sqrt(x/(1-x))
;oplot, x, y, psym=-3, color=0, thick=3.0


; --------------------

; add some real data!
;oplot_bender_vsig_data, 1
;oplot_davies_vsig_data, 1
;oplot_dezeeuw_vsig_data, 1

;oplot_bender_vsig_data, 1, /addkey
;oplot_davies_vsig_data, 1, /addkey
;oplot_deZeeuw_vsig_data, 1, /addkey




; --------------------

device, /close



end



;====================================================================












;------------------------------------------------------------------------
;
;     Density Map of v/sigma vs. Ellipticity
;     ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro vsig_map, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "vsig_map, junk"
   print, "  "
   return
endif

filename='vsig.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3



;--------------------------------------
;--------------------------------------


x0= 0.18
x1= 0.98

y0= 0.15
y1= 0.98


;--------------------------------------
;--------------------------------------


   ;fruns_40gas= ['data/ds/vc3vc3b', $
   fruns_40gas= ['data/ds/vc3vc3c', $
           ;'data/ds/vc3vc3c', $
           'data/ds/vc3vc3d', $
           ;'data/ds/vc3vc3e', $
           ;'data/ds/vc3vc3f', $
           'data/ds/vc3vc3g', $
           'data/ds/vc3vc3h', $
           'data/ds/vc3vc3i', $
           'data/ds/vc3vc3j', $
           'data/ds/vc3vc3k', $
           'data/ds/vc3vc3l', $
           'data/ds/vc3vc3m', $
           'data/ds/vc3vc3n', $
           'data/ds/vc3vc3o', $
           'data/ds/vc3vc3p'] 


   fruns_dless= ['data/collisionless/cvc3vc3b', $
	   ;'data/collisionless/cvc3vc3c', $
	   'data/collisionless/cvc3vc3d', $
	   'data/collisionless/cvc3vc3e', $
	   'data/collisionless/cvc3vc3f', $
	   'data/collisionless/cvc3vc3g', $
	   'data/collisionless/cvc3vc3h', $
	   'data/collisionless/cvc3vc3i', $
	   'data/collisionless/cvc3vc3j', $
	   'data/collisionless/cvc3vc3k', $
	   'data/collisionless/cvc3vc3l', $
	   'data/collisionless/cvc3vc3m', $
	   'data/collisionless/cvc3vc3n', $
	   ;'data/collisionless/cvc3vc3o', $
	   'data/collisionless/cvc3vc3p']


  ;fruns_remerg_b= ['data/remergers/b2eb2e', $
           ;'data/remergers/b2eb2h', $
           ;'data/remergers/b2eb2k', $
           ;'data/remergers/b3eb3e', $
           ;'data/remergers/b3eb3f', $
  fruns_remerg_b= ['data/remergers/b3eb3f', $
           ;'data/remergers/b3eb3h', $
           'data/remergers/b4eb3e', $
           'data/remergers/b4eb4e', $
           'data/remergers/b4eb4k', $
           'data/remergers/b4fb4h']


  fruns_remerg_v= ['data/remergers/vc3fvc3f', $
           'data/remergers/vc3fvc3f_polar', $
           'data/remergers/vc3hvc3h_1', $
           'data/remergers/vc3hvc3h_2', $
           'data/remergers/vc3hvc3h_3', $
           'data/remergers/vc3hvc3h_4', $
           'data/remergers/vc3mvc3m', $
           'data/remergers/vc3mvc3m_polar', $
           'data/remergers/vc3vc3i_rem3', $
           'data/remergers/vc3vc3i_rem4']


  fruns_remerg_e= ['data/remergers/e2ee2e', $
           'data/remergers/e2ke2k', $
           'data/remergers/e3ee3e', $
           'data/remergers/e3ke3e', $
           'data/remergers/e3ke3k', $
           'data/remergers/e4ee4e', $
           'data/remergers/e4ke4k']


  fruns_remerg_z= ['data/remergers/z3fz3f', $
           'data/remergers/z3fz3f_polar', $
           'data/remergers/z3hz3h', $
           'data/remergers/z3hz3h_polar', $
           'data/remergers/z3mz3m', $
           'data/remergers/z3mz3m_polar']




;  now plot things up
; ---------------------------

do_40gas= 0
;do_40gas= 1
if do_40gas eq 1 then begin
   process_one_map, 'vsig', fruns_40gas, x0, y0, x1, y1, /showyaxis, /showxaxis
   xyouts, 0.26, 0.88, /normal, '40% gas', size=1.5, charthick= 3.0, color=0
   xyouts, 0.24, 0.83, /normal, '(all orient)', size=1.5, charthick= 3.0, color=0
   ;process_one_map, 'vsig', fruns_40gas, x0, y0, x1, y1, /showyaxis, /showxaxis, /docontour
   ;xyouts, 0.25, 0.90, /normal, 'contours - 40% gas (all orient)', size=1.2, charthick= 3.0, color=0
endif


do_dless= 0
;do_dless= 1
if do_dless eq 1 then begin
   process_one_map, 'vsig', fruns_dless, x0, y0, x1, y1, /showyaxis, /showxaxis
   xyouts, 0.24, 0.88, /normal, 'dissipationless', size=1.5, charthick= 3.0, color=0
   xyouts, 0.26, 0.83, /normal, '(all orient)', size=1.5, charthick= 3.0, color=0
endif


;do_remerg= 0
do_remerg= 1
if do_remerg eq 1 then begin
   ;fruns_remerg= [fruns_remerg_e(*),fruns_remerg_z(*)]
   fruns_remerg= fruns_remerg_b & msg='map-remergers (b)'
   ;fruns_remerg= fruns_remerg_v & msg='map-remergers (v)'
   ;fruns_remerg= fruns_remerg_e & msg='map-remergers (e)'
   ;fruns_remerg= fruns_remerg_z & msg='map-remergers (z)'
   process_one_map, 'vsig', fruns_remerg, x0, y0, x1, y1, /showyaxis, /showxaxis
   ;xyouts, 0.25, 0.85, /normal, 'map - remergers', size=1.2, charthick= 3.0, color=0
   xyouts, 0.25, 0.85, /normal, msg, size=1.2, charthick= 3.0, color=0
endif


;do_40gas= 0
do_40gas= 1
if do_40gas eq 1 then begin
   process_one_map, 'vsig', fruns_40gas, x0, y0, x1, y1, /docontour
   xyouts, 0.25, 0.90, /normal, 'contours - 40% gas (all orient)', size=1.2, charthick= 3.0, color=0
endif



; --------------------

; line of isotropic rotator
x=findgen(50)*(0.6/50.)
y=sqrt(x/(1-x))
oplot, x, y, psym=-3, color=0, thick=3.0


; --------------------

; add some real data!
oplot_bender_vsig_data, 1
oplot_davies_vsig_data, 1
oplot_dezeeuw_vsig_data, 1

;oplot_bender_vsig_data, 1, /addkey
;oplot_davies_vsig_data, 1, /addkey
;oplot_deZeeuw_vsig_data, 1, /addkey




; --------------------

device, /close



end






;===================================================================








;------------------------------------------------------------------------
;
;     Separate Map of ellip. versus disky/boxy
;       for each orbit
;   ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro vsig_map_5x3, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "vsig_map_5x3, junk"
   print, "  "
   return
endif

filename='vsig_map_5x3.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3, newxsize=22, newysize=14


;--------------------------------------
;--------------------------------------

;
;  WARNING: xmax, xmin, ymax, ymin, bins,
;   xaxistitle, yaxistitle, are all set
;   in process_one_vsig_file
;


;  -------------------------------
;  |     |     |     |     |     |
;  |     |     |     |     |     |
;  |     |     |     |     |     |
;  -------------------------------
;  |     |     |     |     |     |
;  |     |     |     |     |     |
;  |     |     |     |     |     |
;  -------------------------------
;  |     |     |     |     |     |
;  |     |     |     |     |     |
;  |     |     |     |     |     |
;  -------------------------------


x0= 0.07
xsize= 0.184
x1= x0+xsize
x2= x0+xsize+xsize
x3= x0+xsize+xsize+xsize
x4= x0+xsize+xsize+xsize+xsize
x5= x0+xsize+xsize+xsize+xsize+xsize

y0= 0.11
ysize= 0.293333
y1= y0+ysize
y2= y0+ysize+ysize
y3= y0+ysize+ysize+ysize



; orientation vc3's (40% gas)
; ---------------------------
;do_orient= 0
do_orient= 1
ddir= 13
if do_orient eq 1 then begin
	;process_one_map, 'vsig', 'data/ds/vc3vc3h', x0, y2, x1, y3, ddir=ddir, /showyaxis, msg='h', secmsg='(pro-pro)', /showdelta
	;process_one_map, 'vsig', 'data/ds/vc3vc3b', x1, y2, x2, y3, ddir=ddir, msg='b', secmsg='(pro-ret)'
	;process_one_map, 'vsig', 'data/ds/vc3vc3c', x2, y2, x3, y3, ddir=ddir, msg='c', secmsg='(ret-ret)'
	;process_one_map, 'vsig', 'data/ds/vc3vc3d', x3, y2, x4, y3, ddir=ddir, msg='d', secmsg='(polar1)'
	;process_one_map, 'vsig', 'data/ds/vc3vc3e', x4, y2, x5, y3, ddir=ddir, msg='e', secmsg='(tilted1)'
	;process_one_map, 'vsig', 'data/ds/vc3vc3f', x0, y1, x1, y2, ddir=ddir, /showyaxis, msg='f', secmsg='(polar2)'
	;process_one_map, 'vsig', 'data/ds/vc3vc3g', x1, y1, x2, y2, ddir=ddir, msg='g', secmsg='(tilted2)'
	;process_one_map, 'vsig', 'data/ds/vc3vc3i', x2, y1, x3, y2, ddir=ddir, msg='i'
	;process_one_map, 'vsig', 'data/ds/vc3vc3j', x3, y1, x4, y2, ddir=ddir, msg='j'
	;process_one_map, 'vsig', 'data/ds/vc3vc3k', x4, y1, x5, y2, ddir=ddir, msg='k'
	;process_one_map, 'vsig', 'data/ds/vc3vc3l', x0, y0, x1, y1, ddir=ddir, /showyaxis, /showxaxis, msg='l'
	;process_one_map, 'vsig', 'data/ds/vc3vc3m', x1, y0, x2, y1, ddir=ddir, /showxaxis, msg='m'
	;process_one_map, 'vsig', 'data/ds/vc3vc3n', x2, y0, x3, y1, ddir=ddir, /showxaxis, msg='n'
	;process_one_map, 'vsig', 'data/ds/vc3vc3o', x3, y0, x4, y1, ddir=ddir, /showxaxis, msg='o'
	;process_one_map, 'vsig', 'data/ds/vc3vc3p', x4, y0, x5, y1, ddir=ddir, /showxaxis, msg='p'
	process_one_map, 'vsig', 'data/ds/vc3vc3h', x0, y2, x1, y3, /dopoints, ddir=ddir, /showyaxis, msg='h', secmsg='(pro-pro)', /showdelta
	process_one_map, 'vsig', 'data/ds/vc3vc3b', x1, y2, x2, y3, /dopoints, ddir=ddir, msg='b', secmsg='(pro-ret)'
	process_one_map, 'vsig', 'data/ds/vc3vc3c', x2, y2, x3, y3, /dopoints, ddir=ddir, msg='c', secmsg='(ret-ret)'
	process_one_map, 'vsig', 'data/ds/vc3vc3d', x3, y2, x4, y3, /dopoints, ddir=ddir, msg='d', secmsg='(polar1)'
	process_one_map, 'vsig', 'data/ds/vc3vc3e', x4, y2, x5, y3, /dopoints, ddir=ddir, msg='e', secmsg='(tilted1)'
	process_one_map, 'vsig', 'data/ds/vc3vc3f', x0, y1, x1, y2, /dopoints, ddir=ddir, /showyaxis, msg='f', secmsg='(polar2)'
	process_one_map, 'vsig', 'data/ds/vc3vc3g', x1, y1, x2, y2, /dopoints, ddir=ddir, msg='g', secmsg='(tilted2)'
	process_one_map, 'vsig', 'data/ds/vc3vc3i', x2, y1, x3, y2, /dopoints, ddir=ddir, msg='i'
	process_one_map, 'vsig', 'data/ds/vc3vc3j', x3, y1, x4, y2, /dopoints, ddir=ddir, msg='j'
	process_one_map, 'vsig', 'data/ds/vc3vc3k', x4, y1, x5, y2, /dopoints, ddir=ddir, msg='k'
	process_one_map, 'vsig', 'data/ds/vc3vc3l', x0, y0, x1, y1, /dopoints, ddir=ddir, /showyaxis, /showxaxis, msg='l'
	process_one_map, 'vsig', 'data/ds/vc3vc3m', x1, y0, x2, y1, /dopoints, ddir=ddir, /showxaxis, msg='m'
	process_one_map, 'vsig', 'data/ds/vc3vc3n', x2, y0, x3, y1, /dopoints, ddir=ddir, /showxaxis, msg='n'
	process_one_map, 'vsig', 'data/ds/vc3vc3o', x3, y0, x4, y1, /dopoints, ddir=ddir, /showxaxis, msg='o'
	process_one_map, 'vsig', 'data/ds/vc3vc3p', x4, y0, x5, y1, /dopoints, ddir=ddir, /showxaxis, msg='p'
endif



;   *collisionless*
;  orientation vc3's
; ---------------------------
do_orient= 0
;do_orient= 1
if do_orient eq 1 then begin
	;process_one_map, 'vsig', 'cvc3vc3h', x0, y2, x1, y3, ddir=2, /showyaxis, msg='h', secmsg='(pro-pro)', /showdelta
	;process_one_map, 'vsig', 'cvc3vc3b', x1, y2, x2, y3, ddir=2, msg='b', secmsg='(pro-ret)'
	;process_one_map, 'vsig', 'cvc3vc3c', x2, y2, x3, y3, ddir=2, msg='c', secmsg='(ret-ret)'
	;process_one_map, 'vsig', 'cvc3vc3d', x3, y2, x4, y3, ddir=2, msg='d', secmsg='(polar1)'
	;process_one_map, 'vsig', 'cvc3vc3e', x4, y2, x5, y3, ddir=2, msg='e', secmsg='(tilted1)'
	;process_one_map, 'vsig', 'cvc3vc3f', x0, y1, x1, y2, ddir=2, /showyaxis, msg='f', secmsg='(polar2)'
	;process_one_map, 'vsig', 'cvc3vc3g', x1, y1, x2, y2, ddir=2, msg='g', secmsg='(tilted2)'
	;process_one_map, 'vsig', 'cvc3vc3i', x2, y1, x3, y2, ddir=2, msg='i'
	;process_one_map, 'vsig', 'cvc3vc3j', x3, y1, x4, y2, ddir=2, msg='j'
	;process_one_map, 'vsig', 'cvc3vc3k', x4, y1, x5, y2, ddir=2, msg='k'
	;process_one_map, 'vsig', 'cvc3vc3l', x0, y0, x1, y1, ddir=2, /showyaxis, /showxaxis, msg='l'
	;process_one_map, 'vsig', 'cvc3vc3m', x1, y0, x2, y1, ddir=2, /showxaxis, msg='m'
	;process_one_map, 'vsig', 'cvc3vc3n', x2, y0, x3, y1, ddir=2, /showxaxis, msg='n'
	;process_one_map, 'vsig', 'cvc3vc3o', x3, y0, x4, y1, ddir=2, /showxaxis, msg='o'
	;process_one_map, 'vsig', 'cvc3vc3p', x4, y0, x5, y1, ddir=2, /showxaxis, msg='p'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3h', x0, y2, x1, y3, /dopoints, ddir=2, /showyaxis, msg='h', secmsg='(pro-pro)', /showdelta
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3b', x1, y2, x2, y3, /dopoints, ddir=2, msg='b', secmsg='(pro-ret)'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3c', x2, y2, x3, y3, /dopoints, ddir=2, msg='c', secmsg='(ret-ret)'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3d', x3, y2, x4, y3, /dopoints, ddir=2, msg='d', secmsg='(polar1)'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3e', x4, y2, x5, y3, /dopoints, ddir=2, msg='e', secmsg='(tilted1)'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3f', x0, y1, x1, y2, /dopoints, ddir=2, /showyaxis, msg='f', secmsg='(polar2)'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3g', x1, y1, x2, y2, /dopoints, ddir=2, msg='g', secmsg='(tilted2)'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3i', x2, y1, x3, y2, /dopoints, ddir=2, msg='i'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3j', x3, y1, x4, y2, /dopoints, ddir=2, msg='j'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3k', x4, y1, x5, y2, /dopoints, ddir=2, msg='k'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3l', x0, y0, x1, y1, /dopoints, ddir=2, /showyaxis, /showxaxis, msg='l'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3m', x1, y0, x2, y1, /dopoints, ddir=2, /showxaxis, msg='m'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3n', x2, y0, x3, y1, /dopoints, ddir=2, /showxaxis, msg='n'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3o', x3, y0, x4, y1, /dopoints, ddir=2, /showxaxis, msg='o'
	process_one_map, 'vsig', 'data/collisionless/cvc3vc3p', x4, y0, x5, y1, /dopoints, ddir=2, /showxaxis, msg='p'
endif



print, "WARNING: msg's have been manually moved around - check their status"


device, /close



end





; ========================================================================
; ========================================================================
;
;
;
;
;     process_ stuff
;
;
;
; ========================================================================
; ========================================================================




;
;
;
; ---------------------------------------------------
pro process_one_map, maptype, fruns, x0, y0, x1, y1, $
				showyaxis=showyaxis, showxaxis=showxaxis, $
				msg=msg, secmsg=secmsg, $
				showdelta=showdelta, $
				docontour=docontour, dopoints=dopoints, $
				show_x_eq_0=show_x_eq_0, show_y_eq_0=show_y_eq_0, $
				ddir=ddir, $   ; not currently used
				color= color


if not keyword_set(color) then thiscolor= 150 else thiscolor= color

bins= 40

for i=0,n_elements(fruns)-1 do begin

	if maptype eq 'vsig' then begin
		; ellipticity
		;xaxistitle= 'ellipticity'
		xaxistitle= '!7e!6'
		xmax= 0.7
		xmin= 0.0

		; v/sigma
		;yaxistitle= 'v/sigma'
		;yaxistitle= '!7t!3/!7r!3'
		yaxistitle= '!6V!Dmaj!N/!7r!6'
		;yaxistitle= 'V!Dmin!N/!7r!3'
		ymax = 1.7
		ymin = 0.0

		;read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa
		read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa, /small_fov

		this_xvar= ellip
		;this_yvar= Vrot_maj/Sig0
		this_yvar= Vrotavg_maj/Sig0
	endif

	if maptype eq 'db' then begin
		; ellipticity
		xaxistitle= 'ellipticity, !7e!3'
		xmax= 0.7
		xmin= 0.0

		; disky/boxy
		yaxistitle= '!6100 a!D4!N/a'
		ymax = 4.2
		ymin =-4.2

		;read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa
		read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa, /small_fov

		this_yvar= 100.0 * a4diva
		this_xvar= ellip
	endif

	if maptype eq 'kb1' then begin
		; a4/a
		xaxistitle= '!6100 a!D4!N/a'
		xmax= 3.5
		xmin= -2.0

		; (v/sig)*
		yaxistitle= '!6Log (V/!7r!6)!E*!N'
		ymax = 0.5
		ymin = -1.6

		;read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa
		read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa, /small_fov

		this_xvar= 100.0 * a4diva
		this_yvar= alog10((Vrotavg_maj/Sig0)/sqrt(ellip/(1-ellip)))
	endif

	if maptype eq 'kb2' then begin
		; a4/a
		xaxistitle= '!6100 a!D4!N/a'
		xmax= 3.5
		xmin= -2.0

		; Vmin/(|V|)
		yaxistitle= '!6V!Dmin!N/(V!Dmaj!N!E2!N+V!Dmin!N!E2!N)!E1/2!N'
		ymax = 1.05
		ymin = -0.05

		;read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa
		read_kinematics, fruns[i], Re, Vrot_maj, Vrotavg_maj, Vrot_min, Vrotavg_min, Sig0, aa, bb, ellip, a4diva, pa, /small_fov

		this_xvar= 100.0 * a4diva
		;this_yvar=  Vrotavg_min / sqrt(Vrotavg_min*Vrotavg_min + Vrotavg_maj*Vrotavg_maj)
		this_yvar=  Vrotavg_min / sqrt(Vrot_min*Vrot_min + Vrot_maj*Vrot_maj)
	endif


	; load variables
	;
	if i eq 0 then begin
		xvar= this_xvar
		yvar= this_yvar
	endif else begin
		xvar= [xvar,this_xvar]
		yvar= [yvar,this_yvar]
	endelse

endfor


print, xaxistitle+"  <avg>=",mean(xvar), "  +/-", sqrt(variance(xvar))
print, yaxistitle+"  <avg>=",mean(yvar), "  +/-", sqrt(variance(yvar))

; ------------------
; compute histogram
; ------------------

contour_makegeneralpic, xvar, yvar, xmax, xmin, ymax, ymin, $
				pixels= bins, $
				NxNImage=NxNImage


;print, min(NxNImage), max(NxNImage)


x_size= (x1-x0)
y_size= (y1-y0)


if not keyword_set(docontour) and not keyword_set(dopoints) then begin
	tv, NxNImage, x0, y0, xsize=x_size, ysize=y_size, /normal
endif




;
;   generate axes
; ------------------

!p.position=[x0, y0, x0+x_size,y0+y_size]

tnames=['0',' ','0.2',' ','0.4',' ','0.6',' ']

if keyword_set(showyaxis) and keyword_set(showxaxis) then begin
   plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
	xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
	xticks= 7, xtickname=tnames, $
	xtitle=xaxistitle, ytitle=yaxistitle
endif else begin
   if keyword_set(showyaxis) then begin
        plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
	xticks= 7, $
        ytitle=yaxistitle, xtickformat='(a1)'
   endif

   if keyword_set(showxaxis) then begin
        plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
	xticks= 7, xtickname=tnames, $
        xtitle=xaxistitle, ytickformat='(a1)'
   endif

   if ((not keyword_set(showyaxis)) and (not keyword_set(showxaxis))) then begin
	plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
	xticks= 7, $
        xtickformat='(a1)', ytickformat='(a1)'
   endif

endelse



;
;  now do the actual plotting
; -----------------------------

if keyword_set(docontour) then begin
	oplot_contours, NxNImage, x0, y0, x1, y1, pixels= bins, clr=0   ; black (default)
	;oplot_contours, NxNImage, x0, y0, x1, y1, pixels= bins, clr=1   ; white
endif

if keyword_set(dopoints) then begin
	oplot, xvar, yvar, psym=7, color= thiscolor, symsize=1.2
endif



if keyword_set(show_x_eq_0) then begin
        y= [ymin,ymax]
        x= [0.0,0.0]
        oplot, x, y, psym=-3,color=0,thick=4.0,linestyle=2
endif
if keyword_set(show_y_eq_0) then begin
        y= [0.0,0.0]
        x= [xmin,xmax]
        oplot, x, y, psym=-3,color=0,thick=4.0,linestyle=2
endif


if keyword_set(msg) then begin
   ;xyouts, 0.65, 0.90, /normal, fload_fid(1), size=1.33, color=0
   xyouts, 0.08, 1.4, /data, msg, size=1.2, charthick= 3.0, color=0    ; default one
   ;xyouts, 0.05, 1.4, /data, msg, size=1.2, charthick= 3.0, color=0    ; vsig masses
endif

if keyword_set(secmsg) then begin
   ;xyouts, 0.65, 0.90, /normal, fload_fid(1), size=1.33, color=0
   xyouts, 0.03, 1.25, /data, secmsg, size=1.1, color=0           ; default
   ;xyouts, 0.09, 1.25, /data, secmsg, size=1.2, color=0           ; vsig masses
endif





; --------------------
; line for isotropic rotator
;show_isorot_line= 0
show_isorot_line= 1
if keyword_set(show_isorot_line) then begin
	x=findgen(50)*(0.6/50.)
	y=sqrt(x/(1-x))
	oplot, x, y, psym=-3, color=0, thick=3.0
endif



; --------------------
; line for disky/boxy
if keyword_set(show_db_line) then begin
	x=[0,0]
	y=[-0.5,1.5]
	oplot, x, y, psym=-3, color=0, thick=3.0, linestyle= 2
endif



end


; ------------------------------------------------------------------



; ------------------------------------------------------------------


