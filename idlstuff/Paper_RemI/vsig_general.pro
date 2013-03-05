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


bins= 40


;filename='/home/sdutta/Ellipse/DirectionProj/vc3vc3e.prof'


; orientation vc3's (40% gas)
; ---------------------------
do_orient= 0
;do_orient= 1
ddir= 0    ; vc3vc3_major
;ddir= 1    ; vc3vc3_minor
;ddir= 5    ; vc3vc3_new
;ddir= 6    ; vc3vc3_old
ddir= 13    ; vc3vc3_mass
if do_orient eq 1 then begin
	read_vsig_file, 'vc3vc3b', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= vmax/sigma
	allvmax= vmax
	allsig= sigma
	ell= ellipticity

	read_vsig_file, 'vc3vc3c', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'vc3vc3d', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'vc3vc3e', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]
	;vsig= vmax/sigma
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	;ell= ellipticity

	read_vsig_file, 'vc3vc3f', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'vc3vc3g', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'vc3vc3h', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]
	;vsig= vmax/sigma
	;ell= ellipticity

	read_vsig_file, 'vc3vc3i', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'vc3vc3j', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'vc3vc3k', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'vc3vc3l', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'vc3vc3m', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'vc3vc3n', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'vc3vc3o', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'vc3vc3p', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]
endif



;   *collisionless*
; orientation vc3's (40% gas)
; ---------------------------
;do_orient= 0
do_orient= 1
ddir= 2    ; cvc3vc3_major
;ddir= 3    ; cvc3vc3_minor
if do_orient eq 1 then begin
	read_vsig_file, 'cvc3vc3b', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= vmax/sigma
	allvmax= vmax
	allsig= sigma
	ell= ellipticity

	read_vsig_file, 'cvc3vc3c', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'cvc3vc3d', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'cvc3vc3e', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]
	;;vsig= vmax/sigma
	;;ell= ellipticity

	read_vsig_file, 'cvc3vc3f', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'cvc3vc3g', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'cvc3vc3h', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]
	;vsig= vmax/sigma
	;ell= ellipticity

	read_vsig_file, 'cvc3vc3i', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'cvc3vc3j', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'cvc3vc3k', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'cvc3vc3l', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'cvc3vc3m', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'cvc3vc3n', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'cvc3vc3o', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]

	read_vsig_file, 'cvc3vc3p', phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
	vsig= [vsig,vmax/sigma]
	allvmax= [allvmax,vmax]
	allsig= [allsig,sigma]
	ell= [ell,ellipticity]
endif



print, "<v>=",mean(allvmax), "  +/-", sqrt(variance(allvmax))
print, "<sigma>=",mean(allsig), "  +/-", sqrt(variance(allsig))


; ------------------
; compute histogram
; ------------------

contour_makegeneralpic, ell, vsig, xmax, xmin, ymax, ymin, $
				pixels= bins, $
				NxNImage=NxNImage



;---------------------------------------
;  Print Tv image
;-----------------------------------------

x0= 0.15
y0= 0.15
x_size= 0.80
y_size= 0.80
           
                
tv, NxNImage, x0, y0, xsize=x_size, ysize=y_size, /normal



;----------------------------------------
; Generate plot
;----------------------------------------

!p.position=[x0, y0, x0+x_size,y0+y_size]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, $
        xrange=[xmin,xmax],$
        yrange=[ymin,ymax],$
        color=0, $
        xcharsize=1.5, ycharsize=1.5, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	;xticks= 14, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /noerase, $
        /nodata



;oplot, ell, vsig, psym=2, color=0


;xyouts, 0.22, 0.90, /normal, fload_timelbl(1), size=1.0, color=0
;xyouts, 0.65, 0.90, /normal, fload_fid(1), size=1.33, color=0

;xyouts, 0.30, 0.80, 'old stars', /normal, size=1.5, charthick=3.0, color= 0
;xyouts, 0.70, 0.30, 'new stars', /normal, size=1.5, charthick=3.0, color= 0

;xyouts, 0.72, 0.83, '!640% gas', /normal, size=1.5, charthick=3.0, color= 0
xyouts, 0.56, 0.83, '!6dissipationless', /normal, size=1.5, charthick=3.0, color= 0


; --------------------

; line of isotropic rotator
x=findgen(50)*(0.6/50.)
y=sqrt(x/(1-x))
oplot, x, y, psym=-3, color=0, thick=3.0


; --------------------

; add some real data!
readandplot_bender_data, 1
readandplot_davies_data, 1
readandplot_deZeeuw_data, 1

;readandplot_bender_data, 1, /addkey
;readandplot_davies_data, 1, /addkey
;readandplot_deZeeuw_data, 1, /addkey




; --------------------

device, /close



end








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








;------------------------------------------------------------------------
;
;     Separate Map of v/sigma vs. Ellipticity
;       for each orbit
;   ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro vsig_map_orb, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "vsig_map_orb, junk"
   print, "  "
   return
endif

filename='vsig_orb.eps'


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
do_orient= 0
;do_orient= 1
ddir= 13
if do_orient eq 1 then begin
	process_one_vsig_file, 'vc3vc3h', x0, y2, x1, y3, ddir=ddir, /showyaxis, msg='h', secmsg='(pro-pro)', /showdelta
	process_one_vsig_file, 'vc3vc3b', x1, y2, x2, y3, ddir=ddir, msg='b', secmsg='(pro-ret)'
	process_one_vsig_file, 'vc3vc3c', x2, y2, x3, y3, ddir=ddir, msg='c', secmsg='(ret-ret)'
	process_one_vsig_file, 'vc3vc3d', x3, y2, x4, y3, ddir=ddir, msg='d', secmsg='(polar1)'
	process_one_vsig_file, 'vc3vc3e', x4, y2, x5, y3, ddir=ddir, msg='e', secmsg='(tilted1)'
	process_one_vsig_file, 'vc3vc3f', x0, y1, x1, y2, ddir=ddir, /showyaxis, msg='f', secmsg='(polar2)'
	process_one_vsig_file, 'vc3vc3g', x1, y1, x2, y2, ddir=ddir, msg='g', secmsg='(tilted2)'
	process_one_vsig_file, 'vc3vc3i', x2, y1, x3, y2, ddir=ddir, msg='i'
	process_one_vsig_file, 'vc3vc3j', x3, y1, x4, y2, ddir=ddir, msg='j'
	process_one_vsig_file, 'vc3vc3k', x4, y1, x5, y2, ddir=ddir, msg='k'
	process_one_vsig_file, 'vc3vc3l', x0, y0, x1, y1, ddir=ddir, /showyaxis, /showxaxis, msg='l'
	process_one_vsig_file, 'vc3vc3m', x1, y0, x2, y1, ddir=ddir, /showxaxis, msg='m'
	process_one_vsig_file, 'vc3vc3n', x2, y0, x3, y1, ddir=ddir, /showxaxis, msg='n'
	process_one_vsig_file, 'vc3vc3o', x3, y0, x4, y1, ddir=ddir, /showxaxis, msg='o'
	process_one_vsig_file, 'vc3vc3p', x4, y0, x5, y1, ddir=ddir, /showxaxis, msg='p'
endif



;   *collisionless*
;  orientation vc3's
; ---------------------------
;do_orient= 0
do_orient= 1
if do_orient eq 1 then begin
	process_one_vsig_file, 'cvc3vc3h', x0, y2, x1, y3, ddir=2, /showyaxis, msg='h', secmsg='(pro-pro)', /showdelta
	process_one_vsig_file, 'cvc3vc3b', x1, y2, x2, y3, ddir=2, msg='b', secmsg='(pro-ret)'
	process_one_vsig_file, 'cvc3vc3c', x2, y2, x3, y3, ddir=2, msg='c', secmsg='(ret-ret)'
	process_one_vsig_file, 'cvc3vc3d', x3, y2, x4, y3, ddir=2, msg='d', secmsg='(polar1)'
	process_one_vsig_file, 'cvc3vc3e', x4, y2, x5, y3, ddir=2, msg='e', secmsg='(tilted1)'
	process_one_vsig_file, 'cvc3vc3f', x0, y1, x1, y2, ddir=2, /showyaxis, msg='f', secmsg='(polar2)'
	process_one_vsig_file, 'cvc3vc3g', x1, y1, x2, y2, ddir=2, msg='g', secmsg='(tilted2)'
	process_one_vsig_file, 'cvc3vc3i', x2, y1, x3, y2, ddir=2, msg='i'
	process_one_vsig_file, 'cvc3vc3j', x3, y1, x4, y2, ddir=2, msg='j'
	process_one_vsig_file, 'cvc3vc3k', x4, y1, x5, y2, ddir=2, msg='k'
	process_one_vsig_file, 'cvc3vc3l', x0, y0, x1, y1, ddir=2, /showyaxis, /showxaxis, msg='l'
	process_one_vsig_file, 'cvc3vc3m', x1, y0, x2, y1, ddir=2, /showxaxis, msg='m'
	process_one_vsig_file, 'cvc3vc3n', x2, y0, x3, y1, ddir=2, /showxaxis, msg='n'
	process_one_vsig_file, 'cvc3vc3o', x3, y0, x4, y1, ddir=2, /showxaxis, msg='o'
	process_one_vsig_file, 'cvc3vc3p', x4, y0, x5, y1, ddir=2, /showxaxis, msg='p'
endif



print, "WARNING: msg's have been manually moved around - check their status"


device, /close



end





; ========================================================================








;------------------------------------------------------------------------
;
;     Separate Map of v/sigma vs. Ellipticity
;       for masses vc0-vc6
;   ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro vsig_map_mass, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "vsig_map_mass, junk"
   print, "  "
   return
endif

filename='vsig_mass.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3, newxsize=21, newysize=14


;--------------------------------------
;--------------------------------------

;
;  WARNING: xmax, xmin, ymax, ymin, bins,
;   xaxistitle, yaxistitle, are all set
;   in process_one_vsig_file
;


;  -------------------
;  |     |     |     |
;  |     |     |     |
;  |     |     |     |
;  -------------------
;  |     |     |     |
;  |     |     |     |
;  |     |     |     |
;  -------------------
;

x0= 0.08
xsize= 0.30
x1= x0+xsize
x2= x0+xsize+xsize
x3= x0+xsize+xsize+xsize
x4= x0+xsize+xsize+xsize+xsize

y0= 0.08
ysize= 0.45
y1= y0+ysize
y2= y0+ysize+ysize



; vc masses  (most 20% gas)
; ---------------------------
;do_masses= 0
do_masses= 1
if do_masses eq 1 then begin
	process_one_vsig_file, 'vc1vc1', x0, y1, x1, y2, ddir= 4, /showyaxis, $
	;	msg='M!Dtot!N=2.4 x 10!E11!N M!D!9n!3!N', secmsg='!7r!3= 56 km s!E-1!N', /showdelta
		msg='!7r!3= 70 km s!E-1!N', /showdelta
	process_one_vsig_file, 'vc2vc2', x1, y1, x2, y2, ddir= 4, $
	;	msg='M!Dtot!N=7.1 x 10!E11!N M!D!9n!3!N', secmsg='!7r!3= 81 km s!E-1!N'
		msg='!7r!3=102 km s!E-1!N'
	;process_one_vsig_file, 'vc3vc3h', x2, y1, x3, y2, ddir= 0, $
	process_one_vsig_file, 'vc3vc3', x2, y1, x3, y2, ddir= 7, $
	;	msg='M!Dtot!N=1.9 x 10!E12!N M!D!9n!3!N', secmsg='!7r!3=156 km s!E-1!N'
		msg='!7r!3=140 km s!E-1!N'
	process_one_vsig_file, 'vc4vc4a', x0, y0, x1, y1, ddir= 4, /showyaxis, /showxaxis, $
	;	msg='M!Dtot!N=5.4 x 10!E12!N M!D!9n!3!N', secmsg='!7r!3=226 km s!E-1!N'
		msg='!7r!3=271 km s!E-1!N'
	process_one_vsig_file, 'vc5vc5a', x1, y0, x2, y1, ddir= 4, /showxaxis, $
	;	msg='M!Dtot!N=1.5 x 10!E13!N M!D!9n!3!N', secmsg='!7r!3=378 km s!E-1!N'
		msg='!7r!3=388 km s!E-1!N'
	process_one_vsig_file, 'vc6vc6a', x2, y0, x3, y1, ddir= 4, /showxaxis, $
	;	msg='M!Dtot!N=5.8 x 10!E13!N M!D!9n!3!N', secmsg='!7r!3=531 km s!E-1!N'
		msg='!7r!3=563 km s!E-1!N'
endif



;   *collisionless*
;      vc masses
; ---------------------------
do_masses= 0
;do_masses= 1
if do_masses eq 1 then begin
endif



device, /close



end




; ========================================================================


;------------------------------------------------------------------------
;
;     Separate Map of v/sigma vs. Ellipticity
;       for various gas fractions
;   ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro vsig_map_gf, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "vsig_map_gf, junk"
   print, "  "
   return
endif

filename='vsig_gf.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3, newxsize=21, newysize=14


;--------------------------------------
;--------------------------------------

;
;  WARNING: xmax, xmin, ymax, ymin, bins,
;   xaxistitle, yaxistitle, are all set
;   in process_one_vsig_file
;


;  -------------------
;  |     |     |     |
;  |     |     |     |
;  |     |     |     |
;  -------------------
;  |     |     |     |
;  |     |     |     |
;  |     |     |     |
;  -------------------
;

x0= 0.08
xsize= 0.30
x1= x0+xsize
x2= x0+xsize+xsize
x3= x0+xsize+xsize+xsize
x4= x0+xsize+xsize+xsize+xsize

y0= 0.08
ysize= 0.45
y1= y0+ysize
y2= y0+ysize+ysize



; vc masses  (most 20% gas)
; ---------------------------
;do_masses= 0
do_masses= 1
if do_masses eq 1 then begin
	process_one_vsig_file, 'cvc3vc3e', x0, y1, x1, y2, ddir= 2, /showyaxis, $
		msg='dissipationless', /showdelta
	process_one_vsig_file, 'vc3vc3_e', x1, y1, x2, y2, ddir= 11, $
		msg='20%'
	process_one_vsig_file, 'vc3vc3e', x2, y1, x3, y2, ddir= 0, $
		msg='40%'
	process_one_vsig_file, 'vc3vc3e_gf6', x0, y0, x1, y1, ddir= 12, /showyaxis, /showxaxis, $
		msg='60%'
	process_one_vsig_file, 'vc3vc3_gf8', x1, y0, x2, y1, ddir= 12, /showxaxis, $
		msg='80%'
	process_one_vsig_file, 'A3e', x2, y0, x3, y1, ddir= 12, /showxaxis, $
		msg='100 %'
endif



;   *collisionless*
;      vc masses
; ---------------------------
do_masses= 0
;do_masses= 1
if do_masses eq 1 then begin
endif



device, /close



end






; ========================================================================


;------------------------------------------------------------------------
;
;     Separate Map of v/sigma vs. Ellipticity
;       for various merger orbits
;   ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro vsig_map_orbs, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "vsig_map_orbs, junk"
   print, "  "
   return
endif

filename='vsig_orbs.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3, newxsize=21, newysize=14


;--------------------------------------
;--------------------------------------

;
;  WARNING: xmax, xmin, ymax, ymin, bins,
;   xaxistitle, yaxistitle, are all set
;   in process_one_vsig_file
;


;  -------------------
;  |     |     |     |
;  |     |     |     |
;  |     |     |     |
;  -------------------
;  |     |     |     |
;  |     |     |     |
;  |     |     |     |
;  -------------------
;

x0= 0.08
xsize= 0.30
x1= x0+xsize
x2= x0+xsize+xsize
x3= x0+xsize+xsize+xsize
x4= x0+xsize+xsize+xsize+xsize

y0= 0.08
ysize= 0.45
y1= y0+ysize
y2= y0+ysize+ysize



; vc masses  (most 20% gas)
; ---------------------------
;do_masses= 0
do_masses= 1
if do_masses eq 1 then begin
	;process_one_vsig_file, 'vc3vc3h', x0, y1, x1, y2, ddir= 0, /showyaxis, $
	process_one_vsig_file, 'vc3vc3e', x0, y1, x1, y2, ddir= 0, /showyaxis, $
		msg='R!Dperi!N= 5 kpc/h', /showdelta
	;process_one_vsig_file, 'ds/d3h2', x1, y1, x2, y2, ddir= 9, $
	process_one_vsig_file, 'ds/d3e2', x1, y1, x2, y2, ddir= 9, $
		msg='R!Dperi!N= 10 kpc/h'
	;process_one_vsig_file, 'ds/d3h3', x2, y1, x3, y2, ddir= 9, $
	process_one_vsig_file, 'ds/d3e3', x2, y1, x3, y2, ddir= 9, $
		msg='R!Dperi!N= 15 kpc/h'
	;process_one_vsig_file, 'ds/d3h4', x0, y0, x1, y1, ddir= 9, /showyaxis, /showxaxis, $
	process_one_vsig_file, 'ds/d3e4', x0, y0, x1, y1, ddir= 9, /showyaxis, /showxaxis, $
		msg='R!Dperi!N= 20 kpc/h'
	;process_one_vsig_file, 'ds/d3h5', x1, y0, x2, y1, ddir= 9, /showxaxis, $
	process_one_vsig_file, 'ds/d3e5', x1, y0, x2, y1, ddir= 9, /showxaxis, $
		msg='R!Dperi!N= 40 kpc/h'
	;process_one_vsig_file, 'A3e', x2, y0, x3, y1, ddir= 12, /showxaxis, $
	;	msg='100 %'
endif


device, /close


end





; ========================================================================
; ========================================================================


;------------------------------------------------------------------------
;
;     Separate Map of v/sigma vs. Ellipticity
;       for and two quantities
;   ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro vsig_map_2, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "vsig_map_2, junk"
   print, "  "
   return
endif

filename='vsig_2.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3, newxsize=16, newysize=8


;--------------------------------------
;--------------------------------------

;
;  WARNING: xmax, xmin, ymax, ymin, bins,
;   xaxistitle, yaxistitle, are all set
;   in process_one_vsig_file
;


;  -------------
;  |     |     |
;  |     |     |
;  |     |     |
;  -------------
;

x0= 0.13
xsize= 0.43
x1= x0+xsize
x2= x0+xsize+xsize

y0= 0.12
ysize= 0.86
y1= y0+ysize



do_bulge= 0
;do_bulge= 1
if do_bulge eq 1 then begin
	frun1='vc3vc3' & msg1= ' ' & ddir1= 7
	frun2='vc3bvc3b' & msg2= 'bulge' & ddir2= 7
endif


do_gf= 0
;do_gf= 1
if do_gf eq 1 then begin
	frun1='vc3vc3e_gf8' & msg1= '8' & ddir1= 12
	frun2='vc3vc3e_gf8a' & msg2= '8a' & ddir2= 12
endif


do_res= 0
;do_res= 1
if do_res eq 1 then begin
	frun1='vc3vc3' & msg1= '1x' & ddir1= 7
	frun2='vc3vc3_8x' & msg2= '8x' & ddir2= 10
	;frun1='vc3bvc3b' & msg1= '1x' & ddir1= 7
	;frun2='vc3bvc3b_8x' & msg2= '8x' & ddir2= 10
endif


;do_bh= 0
do_bh= 1
if do_bh eq 1 then begin
	;frun1='vc3vc3' & msg1= 'black hole' & ddir1= 7
	;frun2='vc3vc3_no' & msg2= 'no black hole' & ddir2= 7
	;frun1='vc3bvc3b' & msg1= 'black hole' & ddir1= 7
	;frun2='vc3bvc3b_no' & msg2= 'no black hole' & ddir2= 7
	frun1='vc3vc3h' & msg1= 'black hole' & ddir1= 0
	frun2='vc3vc3h_no' & msg2= 'no black hole' & ddir2= 8
	;frun1='vc3vc3e' & msg1= 'black hole' & ddir1= 0
	;frun2='vc3vc3e_no' & msg2= 'no black hole' & ddir2= 8
	;frun1='vc3vc3f' & msg1= 'black hole' & ddir1= 0
	;frun2='vc3vc3f_no' & msg2= 'no black hole' & ddir2= 8
	;frun1='vc3vc3j' & msg1= 'black hole' & ddir1= 0
	;frun2='vc3vc3j_no' & msg2= 'no black hole' & ddir2= 8
	;frun1='vc4vc4a' & msg1= 'black hole' & ddir1= 4
	;frun2='vc4vc4_no' & msg2= 'no black hole' & ddir2= 4
	;frun1='vc5vc5a' & msg1= 'black hole' & ddir1= 4
	;frun2='vc5vc5_no' & msg2= 'no black hole' & ddir2= 4
	;frun1='vc6vc6a' & msg1= 'black hole' & ddir1= 4
	;frun2='vc6vc6_no' & msg2= 'no black hole' & ddir2= 4
endif



process_one_vsig_file, frun1, x0, y0, x1, y1, ddir=ddir1, /showyaxis, /showdelta, msg=msg1
process_one_vsig_file, frun2, x1, y0, x2, y1, ddir=ddir2, msg=msg2


device, /close



end




; ========================================================================












;
;
;
; ---------------------------------------------------
pro process_one_vsig_file, frun, x0, y0, x1, y1, $
				showyaxis=showyaxis, showxaxis=showxaxis, $
				msg=msg, $
				secmsg=secmsg, $
				ddir=ddir, $
				showdelta=showdelta



; ellipticity
;xaxistitle= 'ellipticity'
xaxistitle= '!7e!3'
xmax= 0.7
;xmax= 0.69
xmin= 0.0

; v/sigma
;yaxistitle= 'v/sigma'
;yaxistitle= '!7t!3/!7r!3'
yaxistitle= 'V!Dmaj!N/!7r!3'
;yaxistitle= 'V!Dmin!N/!7r!3'
ymax = 1.7
ymin = 0.0


bins= 40




read_vsig_file, frun, phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
vsig= vmax/sigma
ell= ellipticity

; ------------------
; compute histogram
; ------------------

contour_makegeneralpic, ell, vsig, xmax, xmin, ymax, ymin, $
				pixels= bins, $
				NxNImage=NxNImage


;print, min(NxNImage), max(NxNImage)


x_size= (x1-x0)
y_size= (y1-y0)


tv, NxNImage, x0, y0, xsize=x_size, ysize=y_size, /normal



;----------------------------------------
; Generate plot
;----------------------------------------

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


;oplot, ell, vsig, psym=2, color=0


;xyouts, 0.22, 0.90, /normal, fload_timelbl(1), size=1.0, color=0
;xyouts, 0.65, 0.90, /normal, fload_fid(1), size=1.33, color=0

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
; dispersion and ke 
; tensor components
addtvt= 1
;addtvt= 0
if keyword_set(addtvt) then begin
	read_tvt_file, frun, pi_xx, pi_yy, pi_zz, pi_xy, pi_xz, pi_yz, $
                        t_xx, t_yy, t_zz, t_xy, t_xz, t_yz

	Trace_PI= pi_xx + pi_yy + pi_zz
	xzratio= pi_xx / pi_zz
	yzratio= pi_yy / pi_zz
	Trace_T= t_xx + t_yy + t_zz
	Trace_ratio= Trace_T/Trace_PI
	Trace_Total= Trace_T+Trace_PI

	delta_x= 1.0 - (pi_zz/pi_xx)
	delta_y= 1.0 - (pi_zz/pi_yy)

	delta= 1.0 - 2.0 * pi_zz / (pi_xx + pi_yy)

	;xzr = strmid(strcompress(string(xzratio),/remove_all),0,4)   ; 0.xx 
	;xyouts, 0.55, 0.78, xzr, /data, size=0.9, charthick=3.0, color= 0

	;yzr = strmid(strcompress(string(yzratio),/remove_all),0,4)   ; 0.xx 
	;xyouts, 0.55, 0.66, yzr, /data, size=0.9, charthick=3.0, color= 0

	;xzr = strmid(strcompress(string(Trace_ratio),/remove_all),0,4)   ; 0.xx 
	;xyouts, 0.55, 0.54, xzr, /data, size=0.9, charthick=3.0, color= 0

	;xzr = strmid(strcompress(string(Trace_Total),/remove_all),0,6)   ; 0.xx 
	;xyouts, 0.55, 0.42, xzr, /data, size=0.9, charthick=3.0, color= 0

	; delta x
	; ---------
	;if keyword_set(showdelta) then xyouts, 0.45, 0.66, '!7d!3!Dx!N=', /data, size=0.9, charthick=3.0, color= 0
	;if delta_x lt 0.0 then begin
	;	delta_x= -1.0*delta_x
	;	xyouts, 0.51, 0.66, '-', /data, size=0.9, charthick=3.0, color= 0
	;endif
	;xzr = strmid(strcompress(string(delta_x),/remove_all),0,4)   ; 0.xx
	;xyouts, 0.55, 0.66, xzr, /data, size=0.9, charthick=3.0, color= 0

	; delta y
	; ---------
	;if keyword_set(showdelta) then xyouts, 0.45, 0.54, '!7d!3!Dy!N=', /data, size=0.9, charthick=3.0, color= 0
	;if delta_y lt 0.0 then begin
	;	delta_y= -1.0*delta_y
	;	xyouts, 0.51, 0.54, '-', /data, size=0.9, charthick=3.0, color= 0
	;endif
	;xzr = strmid(strcompress(string(delta_y),/remove_all),0,4)   ; 0.xx
	;xyouts, 0.55, 0.54, xzr, /data, size=0.9, charthick=3.0, color= 0

	; delta
	; ------
	if keyword_set(showdelta) then xyouts, 0.45, 0.66, '!7d!3=', /data, size=0.9, charthick=3.0, color= 0
	if delta lt -0.01 then begin
		delta= -1.0*delta
		;xyouts, 0.51, 0.66, '-', /data, size=0.9, charthick=3.0, color= 0
	endif
	if delta lt 0.0 then delta= -1.0*delta
	xzr = strmid(strcompress(string(delta),/remove_all),0,4)   ; 0.xx
	xyouts, 0.55, 0.66, xzr, /data, size=0.9, charthick=3.0, color= 0

endif


; --------------------
; line for isotropic rotator
x=findgen(50)*(0.6/50.)
y=sqrt(x/(1-x))
oplot, x, y, psym=-3, color=0, thick=3.0



end






; =======================================================================================






;------------------------------------------------------------------------
;
;     Density Map of 
;         some measure of minor axis
;         velocity, maybe v_major versus
;         v_tot, (v_maj^2 + v_min^2)^1/2
;     ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro vmin_map, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "vmin_map, junk"
   print, "  "
   return
endif

filename='vmin.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3


;--------------------------------------
;--------------------------------------


; v_maj
;yaxistitle= 'V!Dmaj!N/!7r!3'
;yaxistitle= 'V!Dmaj!N'
yaxistitle= '!6V!Dtot!N=(V!Dmaj!N!E2!N+V!Dmin!N!E2!N)!E1/2!N'
;ymax= 225.0
ymax= 160.0
ymin= 0.0

; v_tot
;xaxistitle= '(V!Dmaj!N!E2!N+V!Dmin!N!E2!N)!E1/2!N'
;xaxistitle= '!6f!Dmin!N=V!Dmin!N/(V!Dmaj!N!E2!N+V!Dmin!N!E2!N)!E1/2!N'
xaxistitle= '!7l!6!N=V!Dmin!N/(V!Dmaj!N!E2!N+V!Dmin!N!E2!N)!E1/2!N'
;ymax = 200.0
xmax = 1.0
xmin = 0.0


bins= 40




; orientation vc3's (40% gas)
; ---------------------------
;do_orient= 0
do_orient= 1
if do_orient eq 1 then begin
	process_vmin, 'vc3vc3b', vmax_maj, vmax_min
	allvmaj= vmax_maj(*) & allvmin= vmax_min(*)
	process_vmin, 'vc3vc3c', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3d', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3e', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3f', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3g', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3h', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3i', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3j', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3k', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3l', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3m', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3n', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3o', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'vc3vc3p', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
endif



;   *collisionless*
; orientation vc3's
; ---------------------------
do_orient= 0
;do_orient= 1
if do_orient eq 1 then begin
	process_vmin, 'cvc3vc3b', vmax_maj, vmax_min
	allvmaj= vmax_maj(*) & allvmin= vmax_min(*)
	process_vmin, 'cvc3vc3c', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3d', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3e', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3f', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3g', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3h', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3i', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3j', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3k', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3l', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3m', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3n', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3o', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
	process_vmin, 'cvc3vc3p', vmax_maj, vmax_min
	allvmaj= [allvmaj, vmax_maj(*)] & allvmin= [allvmin, vmax_min(*)]
endif



print, " --------------------------- "
print, "<v_maj>=",mean(allvmaj), "  +/-", sqrt(variance(allvmaj))
print, "<v_min>=",mean(allvmin), "  +/-", sqrt(variance(allvmin))


; ------------------
; compute histogram
; ------------------

allvtot= sqrt(allvmaj*allvmaj + allvmin*allvmin)
allvminfrac= allvmin/allvtot

print, "<v_tot>=",mean(allvtot), "  +/-", sqrt(variance(allvtot))
print, "<f_min>=",mean(allvminfrac), "  +/-", sqrt(variance(allvminfrac))

contour_makegeneralpic, allvminfrac, allvtot, xmax, xmin, ymax, ymin, $
				pixels= bins, $
				NxNImage=NxNImage



;---------------------------------------
;  Print Tv image
;-----------------------------------------

x0= 0.18
y0= 0.15
x_size= 0.78
y_size= 0.83
           
                
tv, NxNImage, x0, y0, xsize=x_size, ysize=y_size, /normal



;----------------------------------------
; Generate plot
;----------------------------------------

!p.position=[x0, y0, x0+x_size,y0+y_size]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, color= 0, /noerase, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle=xaxistitle, ytitle=yaxistitle, $
        xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, charthick=3.0


; --------------------

y= [0.0, 105.0]
x= y
x(*)= 1./sqrt(2.0)
oplot, x, y, psym=-3,color=0,thick=4.0,linestyle=2

; --------------------


;xyouts, 0.30, 0.80, 'old stars', /normal, size=1.5, charthick=3.0, color= 0
;xyouts, 0.70, 0.30, 'new stars', /normal, size=1.5, charthick=3.0, color= 0

xyouts, 0.50, 0.80, '40% gas', /normal, size=1.5, charthick=3.0, color= 0
;xyouts, 0.30, 0.80, 'dissipationless', /normal, size=1.5, charthick=3.0, color= 0


; --------------------

device, /close


end








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
;print, " --------------------------- "
;print, "<v_maj>=",mean(vmax_maj), "  +/-", sqrt(variance(vmax_maj))
;print, "<v_min>=",mean(vmax_min), "  +/-", sqrt(variance(vmax_min))
;print, "<v_tot>=",mean(vtot), "  +/-", sqrt(variance(vtot))
;print, "<f_min>=",mean(vminfrac), "  +/-", sqrt(variance(vminfrac))
;print, "<v_maj/sigma>=",mean(vmax_maj/sigma), "  +/-", sqrt(variance(vmax_maj/sigma))

	

end






;===================================================================================
;===================================================================================
;===================================================================================
;===================================================================================







pro random_hist, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "random_hist, junk"
   print, "  "
   return
endif

filename='rhist.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;--------------------------------------
;--------------------------------------


xaxistitle= '!6Ellipticity, !7e!6'
;;xaxistitle= 'f!Dmin!N'
;xaxistitle= '!6Anisotropy, !7d!6'
xticklbls=['0',' ','0.2',' ','0.4',' ','0.6',' ','0.8',' ','1']
xticknum= 10
xmax= 1.0
xmin= 0.0

;xaxistitle= '!6R!D!8a!6!N (h!E-1!N kpc)'
;xaxistitle= '!8a!6 (h!E-1!N kpc)'
;xaxistitle= '!8a!6 (kpc)'
;;xticklbls=['1.0','2.0','3.0','4.0','5.0','6.0']
;xticklbls=['2','3','4','5','6','7','8']
;xticklbls=['3','4','5','6','7','8','9','10','11']
;xticknum= n_elements(xticklbls)-1
;xmax= 11.0
;xmin= 3.0

;xaxistitle= '!7r!6 (km s!E-1!N)'
;xticklbls=[' ','120',' ','140',' ','160',' ','180',' ','200',' ',' ']
;xticknum= n_elements(xticklbls)-1
;xmax= 220.0
;xmin= 110.0

;xaxistitle= '!6V!Dmaj!N (km s!E-1!N)'
;;xticklbls=['0','50','100','150','200']
;xticklbls=['0','20','40','60','80','100','120','140','160']
;xticknum= n_elements(xticklbls)-1
;;xticknum= 10
;;xmax= 200.0
;;xticknum= 8
;xmax= 160.0
;xmin= 0.0

;xaxistitle= 'V!Dmaj!N/!7r!3'
;xticklbls=['0',' ','0.2',' ','0.4',' ','0.6',' ','0.8',' ','1.0',' ','1.2',' ','1.4',' ',' ']
;xticknum= n_elements(xticklbls)-1
;xmax= 1.6
;xmin= 0.0

;xaxistitle= '!6(V!Dmaj!N/!7r!3)!E*!N'
;;xticklbls=['0',' ','0.2',' ','0.4',' ','0.6',' ','0.8',' ','1.0',' ','1.2',' ',' ']
;xticklbls=['0',' ','0.2',' ','0.4',' ','0.6',' ','0.8',' ','1.0',' ','1.2',' ','1.4',' ','1.6',' ',' ']
;xticknum= n_elements(xticklbls)-1
;xmax= 1.8
;;xmax= 1.4
;;xmax = 3.0
;xmin= 0.0

;xaxistitle= 'Kinematic Misalignment, !7W!6'
;xticklbls=['0',' ','20',' ','40',' ','60',' ','80',' ']
;xticknum= n_elements(xticklbls)-1
;xmax= 90.0
;xmin= 0.0


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


; --------------------

;  std's
; -----------------

fruns= ['vc3vc3h','vc3vc3b','vc3vc3c','vc3vc3d','vc3vc3e', $
        'vc3vc3f','vc3vc3g','vc3vc3i','vc3vc3j','vc3vc3k', $
        'vc3vc3l','vc3vc3m','vc3vc3n','vc3vc3o','vc3vc3p']
ddir= intarr(15)
;ddir(*)= 0
ddir(*)= 13

;fruns= ['vc3vc3i']
;ddir= [13]


;allre=     fload_avg_re(fruns,ddir,/all)
;allre= allre * 1.3 + 0.3
;allre= allre/0.7
alle=      fload_avg_es(fruns,ddir,/all)
;allvmax=   fload_avg_vmax(fruns,ddir,/all)
;allsig=    fload_avg_sig(fruns,ddir,/all)
;allvs=     fload_avg_vs(fruns,ddir,/all)
;allvsst=   fload_avg_vsst(fruns,ddir,/all)
;allts=     fload_Ts(fruns)
;allds=     fload_deltas(fruns)
;allkms=    fload_kinematicmas(fruns,/all)
;allfmins=  [0.59, 0.64, 0.52, 0.93, 0.78,$ 
;		0.93, 0.70, 0.90, 0.77, 0.87, $
;		0.82, 0.94, 0.91, 0.75, 0.94]


; --------------------


;  collisionless
; -----------------

fruns= ['cvc3vc3h','cvc3vc3b','cvc3vc3c','cvc3vc3d','cvc3vc3e', $
        'cvc3vc3f','cvc3vc3g','cvc3vc3i','cvc3vc3j','cvc3vc3k', $
        'cvc3vc3l','cvc3vc3m','cvc3vc3n','cvc3vc3o','cvc3vc3p']
ddir= intarr(15)
ddir(*)= 2

;fruns= ['cvc3vc3i']
;ddir= [2]


;allcre=     fload_avg_re(fruns,ddir,/all)
;allcre=     allcre * 1.425 + 0.1    ; fix the screwy major axis numbers
;allcre=     allcre / 0.7
allce=      fload_avg_es(fruns,ddir,/all)
;allcvmax=   fload_avg_vmax(fruns,ddir,/all)
;allcsig=    fload_avg_sig(fruns,ddir,/all)
;allcvs=     fload_avg_vs(fruns,ddir,/all)
;allcvsst=   fload_avg_vsst(fruns,ddir,/all)
;allcts=     fload_Ts(fruns)
;allcds=     fload_deltas(fruns)
;allckms=     fload_kinematicmas(fruns,/all,/cless)
;allcfmins=  [0.59, 0.64, 0.52, 0.93, 0.78,$ 
;		0.93, 0.70, 0.90, 0.77, 0.87, $
;		0.82, 0.94, 0.91, 0.75, 0.94]



; --------------------



levels= 50.0
;levels= 20.0
;levels= 10.0


; two arrays to histogram

; re
;alltts= allre
;allctts= allcre
; e
alltts= alle
allctts= allce
; vmax
;alltts= allvmax
;allctts= allcvmax
; sig
;alltts= allsig
;allctts= allcsig
; v/sig
;alltts= allvs
;allctts= allcvs
; v/sig^star
;alltts= allvsst
;allctts= allcvsst
; anisotropy
;alltts= allds
;allctts= allcds
; kinematic misalignment
;alltts= allkms
;allctts= allckms


xyouts, 0.57, 0.85, '!640% gas', /normal, charthick=3.0, size=1.5, color=150
xyouts, 0.57, 0.79, '!6dissipationless', /normal, charthick=3.0, size=1.5, color=50


idx=where(alltts gt 1.0)
print, '40: percent above 1.0', 1.0*n_elements(idx)/n_elements(alltts)
idx=where(allctts gt 1.0)
print, 'c: percent above 1.0', 1.0*n_elements(idx)/n_elements(allctts)

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
;bins= [0,bins, bins(levels-1)]
bins= [0,bins,xmax]
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
;				thick=3.0
;polyfill, nbins, nctts, /data, color= 50, /fill, linestyle=0, $
;				thick=3.0
;polyfill, nbins, nctts, /data, color= 220, /fill, linestyle=0, $
;				thick=3.0

; 40% gas
;polyfill, nbins, ntts, /data, color= 150, /line_fill, linestyle=0, $
;				thick=3.0, orientation=90.0
polyfill, nbins, ntts, /data, color= 150, /line_fill, linestyle=0, $
				thick=3.0, orientation=45.0




; add the data on kinematic misalignment angles
;   from Franx, Illingworth, and Heckman (1989)
;   or newer information from Franx, Ill, and de Z (1991)
; ----------------------------------------------
;add_kinma= 1
add_kinma= 0
if add_kinma eq 1 then begin
	;hi= [14, 14, 4, 0, 2, 0, 0, 0, 1, 1, 1]/16.0     ; this is FIH '89
	hi= [26, 26, 8, 0, 3, 1, 1, 0, 1, 4, 4]/29.0     ; this is FIZ '91
	bins= [0, 5, 15, 25, 35, 45, 55, 65, 75, 85, 90]
	oplot, bins, hi, psym=10, color= 0, thick=12.0

	;xyouts, 0.67, 0.73, 'FIH (1989)', /normal, charthick=4.0, size=1.5, color=0
	xyouts, 0.57, 0.73, 'FIZ (1991)', /normal, charthick=4.0, size=1.5, color=0
endif


; add an arrow for the initial R_e
; ----------------------------------
;add_rearrow= 1
add_rearrow= 0
if add_rearrow eq 1 then begin
	Initial_Re= 3.10 / 0.7
	arrow, Initial_Re, 1.08, Initial_Re, 0.90, COLOR=0, THICK=3.0, hthick=3.0, /data
	x0= 0.29
	y0= 0.85
	xyouts, x0, y0, '!6initial', /normal, charthick=4.0, size=1.2, color=0
	xyouts, x0-0.04, y0-0.04, '!6half-mass', /normal, charthick=4.0, size=1.2, color=0
	xyouts, x0, y0-0.08, '!6radius', /normal, charthick=4.0, size=1.2, color=0
endif


; histogram of ellipticities from Franx, Illingworth, & de Zeeuw
; ------------------------------------------------------------------
add_ehist= 1
;add_ehist= 0
if add_ehist eq 1 then begin
	ellips_ccd= [16,16,22,33,34,27,23,27,14,11,6,4,0,0,0]/40.0
	ellips_7sam= [38,38,30,41,39,52,58,30,26,16,15,5,2,2,0]/65.0
	bins= [0,0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475, $
		0.525,0.575,0.625,0.675]  ;,0.725,0.775,0.825,0.875,0.925,0.975]

	oplot, bins, ellips_ccd, psym=10, color= 0, thick=12.0, linestyle= 1
	oplot, [0.58,0.65], [0.86,0.86], psym=-3, color= 0, thick=12.0, linestyle= 1
	xyouts, 0.67, 0.73, '!6FIZ (1991)', /normal, charthick=4.0, size=1.5, color=0

	oplot, bins, ellips_7sam, psym=10, color= 0, thick=12.0
	oplot, [0.58,0.65], [0.77,0.77], psym=-3, color= 0, thick=12.0
	xyouts, 0.67, 0.67, '!6Faber (1989)', /normal, charthick=4.0, size=1.5, color=0
endif




; histogram of (v/sig)* from all our compiled datasets
; ------------------------------------------------------
add_vsigstar_hist= 0
;add_vsigstar_hist= 1
if add_vsigstar_hist eq 1 then begin
	; davies-ellipticals
	benderfile= '/home/tcox/RotvAnisSupport/davies83_ellipticals.txt'
	read_bender_data, benderfile, ellip, mb, vmax, sigma
	vsigst= (vmax/sigma)/(sqrt(ellip/(1-ellip)))
	; davies-bulges
	benderfile= '/home/tcox/RotvAnisSupport/davies83_bulges.txt'
	read_bender_data, benderfile, ellip, mb, vmax, sigma
	vsigst= [vsigst,(vmax/sigma)/(sqrt(ellip/(1-ellip)))]
	; --------
	; bender '88
	benderfile= '/home/tcox/RotvAnisSupport/bender88_ellipticals.txt'
	read_bender_data, benderfile, ellip, mb, vmax, sigma
	vsigst= [vsigst,(vmax/sigma)/(sqrt(ellip/(1-ellip)))]
	benderfile= '/home/tcox/RotvAnisSupport/bender90_dwarfe_ellipticals.txt'
	read_bender_data, benderfile, ellip, mb, vmax, sigma
	vsigst= [vsigst,(vmax/sigma)/(sqrt(ellip/(1-ellip)))]
	benderfile= '/home/tcox/RotvAnisSupport/bender90_lowlum_ellipticals.txt'
	read_bender_data, benderfile, ellip, mb, vmax, sigma
	vsigst= [vsigst,(vmax/sigma)/(sqrt(ellip/(1-ellip)))]
	; --------
	; de Zeeuw - SAURON
	benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_cluster_es.txt'
	read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
	vsigst= [vsigst,(vmax/sigma)/(sqrt(ellip/(1-ellip)))]
	;benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_cluster_len.txt'
	;read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
	;vsigst= [vsigst,(vmax/sigma)/(sqrt(ellip/(1-ellip)))]
	;benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_cluster_sp.txt'
	;read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
	;vsigst= [vsigst,(vmax/sigma)/(sqrt(ellip/(1-ellip)))]
	benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_field_es.txt'
	read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
	vsigst= [vsigst,(vmax/sigma)/(sqrt(ellip/(1-ellip)))]
	;benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_field_len.txt'
	;read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
	;vsigst= [vsigst,(vmax/sigma)/(sqrt(ellip/(1-ellip)))]
	;benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_field_sp.txt'
	;read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
	;vsigst= [vsigst,(vmax/sigma)/(sqrt(ellip/(1-ellip)))]
	; --------
	; Halliday - low lum e's
	benderfile= '/home/tcox/RotvAnisSupport/halliday01_lowlum_ellipticals.txt'
	read_bender_data, benderfile, ellip, mb, vmax, sigma
	vsigst= [vsigst,(vmax/sigma)/(sqrt(ellip/(1-ellip)))]


	idx= where(vsigst gt 0.0)
	if idx(0) ne -1 then vsigst= vsigst(idx)

	print, "max (v/sig)* = ", max(vsigst)
	idx= where(vsigst gt 1.0)
	print, "percent greater than 1.0= ", 1.0*n_elements(idx)/n_elements(vsigst)

	; now, make histogram
	levels= 10.0
	step= (xmax-xmin)/(levels)
	bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
	bins=bins+(step*0.5)

	; std histogram
	hist_allvsigsts= histogram(vsigst, binsize=step, max=xmax, min=xmin)

	; make sure it extends to 0
	bins= [0,bins,xmax]
	hist_allvsigsts= [hist_allvsigsts(0),hist_allvsigsts,hist_allvsigsts(levels-1)]

	normalization= hist_alltts_area / total(hist_allvsigsts*step)
	hist_allvsigsts= hist_allvsigsts*normalization
	;hist_allvsigsts= hist_allvsigsts/normalization

	oplot, bins, hist_allvsigsts, psym=10, color=0, thick=12.0

	oplot, [1.05,1.18], [0.77,0.77], psym=-3, color= 0, thick=12.0
	xyouts, 0.67, 0.67, 'Fig. 4', /normal, charthick=4.0, size=1.5, color=0



	; --------------
	; read bbf '92
	bbffile= '/home/tcox/RotvAnisSupport/bbf92_spheroids.txt'
        read_bbf_data, bbffile, mb, vsigst, a4diva


        idx= where(vsigst gt 0.0)
        if idx(0) ne -1 then vsigst= vsigst(idx)
   
        print, "max (v/sig)* = ", max(vsigst)
        idx= where(vsigst gt 1.0) 
        print, "percent greater than 1.0= ", 1.0*n_elements(idx)/n_elements(vsigst)

        ; now, make histogram
        levels= 10.0
        step= (xmax-xmin)/(levels)
        bins= (IndGen(levels)/levels*(xmax-xmin)) + xmin
        bins=bins+(step*0.5)

        ; std histogram
        hist_allvsigsts= histogram(vsigst, binsize=step, max=xmax, min=xmin)

        ; make sure it extends to 0
        bins= [0,bins,xmax]
        hist_allvsigsts= [hist_allvsigsts(0),hist_allvsigsts,hist_allvsigsts(levels-1)]

        normalization= hist_alltts_area / total(hist_allvsigsts*step)
        hist_allvsigsts= hist_allvsigsts*normalization
        ;hist_allvsigsts= hist_allvsigsts/normalization

        oplot, bins, hist_allvsigsts, psym=10, color=0, thick=12.0, linestyle= 1

	oplot, [1.05,1.18], [0.86,0.86], psym=-3, color= 0, thick=12.0, linestyle= 1
	xyouts, 0.67, 0.73, 'BBF (1992)', /normal, charthick=4.0, size=1.5, color=0

endif


device, /close


end









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
;do_sizes= 0
do_sizes= 1
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


	print, "sig_avgs"
	print, "--------"
	print, "e"
	print, sig_avgs
	print, "h"
	print, sig_avgs2
	print, "k"
	print, sig_avgs3

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
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12, newysize=27

x0= 0.15 & x1= 0.99
y0= 0.07 & ysize=0.184
y1= y0+ysize
y2= y0+ysize+ysize
y3= y0+ysize+ysize+ysize
y4= y0+ysize+ysize+ysize+ysize
y5= y0+ysize+ysize+ysize+ysize+ysize


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
!p.position= [x0, y4, x1, y5]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, xvar, re_avgs, psym=-2, color=150, thick= 3.0, symsize= 1.5
if secondgrp eq 1 then oplot, xvar2, re_avgs2, psym=-8, color=50, thick= 3.0, symsize= 1.5
if thirdgrp eq 1 then oplot, xvar3, re_avgs3, psym=-5, color=100, thick= 3.0, symsize= 1.5
oploterror, xvar, re_avgs, reerr, psym=-3, errcolor=150, color=150, errthick= 3.0, thick= 3.0
if secondgrp eq 1 then oploterror, xvar2, re_avgs2, reerr2, psym=-3, errcolor=50, color=50, errthick= 3.0, thick=3.0
if thirdgrp eq 1 then oploterror, xvar3, re_avgs3, reerr3, psym=-3, errcolor=100, color=100, errthick= 3.0, thick=3.0


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
    if secondgrp eq 1 then oplot, [55.0,90.0], [10.5,10.5], psym=-8, color=50, thick= 3.0, symsize= 1.5
    if thirdgrp eq 1 then oplot, [55.0,90.0], [8.0,8.0], psym=-5, color=100, thick= 3.0, symsize= 1.5
endif




; zeroth: Ecc
; -------------
yaxistitle='!6<!7e!6>'
ymax = 0.9
ymin= 0.0
!p.position= [x0, y3, x1, y4]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, xvar, es_avgs, psym=-2, color=150, thick= 3.0, symsize= 1.5
if secondgrp eq 1 then oplot, xvar2, es_avgs2, psym=-8, color=50, thick= 3.0, symsize= 1.5
if thirdgrp eq 1 then oplot, xvar3, es_avgs3, psym=-5, color=100, thick= 3.0, symsize= 1.5

oploterror, xvar, es_avgs, eserr, psym=-3, errcolor=150, color=150, thick= 3.0, errthick= 3.0
if secondgrp eq 1 then oploterror, xvar2, es_avgs2, eserr2, psym=-3, errcolor=50, color=50, thick= 3.0, errthick= 3.0
if thirdgrp eq 1 then oploterror, xvar3, es_avgs3, eserr3, psym=-3, errcolor=100, color=100, thick= 3.0, errthick= 3.0


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
ymax = 1.09
ymin= 0.0
!p.position= [x0, y2, x1, y3]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, xvar, vs_avgs, psym=-2, color=150, thick= 3.0, symsize= 1.5
if secondgrp eq 1 then oplot, xvar2, vs_avgs2, psym=-8, color=50, thick= 3.0, symsize= 1.5
if thirdgrp eq 1 then oplot, xvar3, vs_avgs3, psym=-5, color=100, thick= 3.0, symsize= 1.5

oploterror, xvar, vs_avgs, vserr, psym=-3, errcolor=150, color=150, thick= 3.0, errthick=3.0
if secondgrp eq 1 then oploterror, xvar2, vs_avgs2, vserr, psym=-3, errcolor=50, color=50, thick= 3.0, errthick=3.0
if thirdgrp eq 1 then oploterror, xvar3, vs_avgs3, vserr, psym=-3, errcolor=100, color=100, thick= 3.0, errthick=3.0


; second: Triaxiality
; ---------------------
yaxistitle='!6Triaxiality, T'
;ymax = 0.9
ymax = 1.1
ymin= 0.0
!p.position= [x0, y1, x1, y2]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, xvar, TTs_avgs, psym=-2, color=150, thick= 3.0, symsize= 1.5
if secondgrp eq 1 then oplot, xvar2, TTs_avgs2, psym=-8, color=50, thick= 3.0, symsize= 1.5
if thirdgrp eq 1 then oplot, xvar3, TTs_avgs3, psym=-5, color=100, thick= 3.0, symsize= 1.5



; second: delta
; ----------------
yaxistitle='!6Anisotropy, !7d!6'
ymax = 0.74
ymin= 0.0
!p.position= [x0, y0, x1, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase

oplot, xvar, deltas_avgs, psym=-2, color=150, thick= 3.0, symsize= 1.5
if secondgrp eq 1 then oplot, xvar2, deltas_avgs2, psym=-8, color=50, thick= 3.0, symsize= 1.5
if thirdgrp eq 1 then oplot, xvar3, deltas_avgs3, psym=-5, color=100, thick= 3.0, symsize= 1.5



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
    xyouts, 0.28, 0.95, cmt1, /normal, charthick=3.0, size=1.5, color=150
    if secondgrp eq 1 then xyouts, 0.28, 0.92, cmt2, /normal, charthick=3.0, size=1.5, color=50
    if thirdgrp eq 1 then xyouts, 0.28, 0.89, cmt3, /normal, charthick=3.0, size=1.5, color=100
endif



; -------------
;  Done
; -------------

device, /close



end







;================================================================================














;=================================================================================
;
;   Global Relationships
;
;=================================================================================








pro global_relations, junk



if not keyword_set(junk) then begin
        print, " "
        print, " global_relations, junk"
        print, " "
        print, " "
        return
endif


; =====================
; =====================


;  std's
; -------------

fruns= ['vc3vc3h','vc3vc3b','vc3vc3c','vc3vc3d','vc3vc3e', $
	'vc3vc3f','vc3vc3g','vc3vc3i','vc3vc3j','vc3vc3k', $
	'vc3vc3l','vc3vc3m','vc3vc3n','vc3vc3o','vc3vc3p']
ddir= intarr(15)
;ddir(*)= 0
ddir(*)= 13


re_avgs=     fload_avg_re(fruns,ddir)
es_avgs=     fload_avg_es(fruns,ddir)
vmax_avgs=   fload_avg_vmax(fruns,ddir)
sig_avgs=    fload_avg_sig(fruns,ddir)
vs_avgs=     fload_avg_vs(fruns,ddir)
vsst_avgs=   fload_avg_vsst(fruns,ddir)
TTs_avgs=    fload_Ts(fruns)
deltas_avgs= fload_deltas(fruns)
fmin_avgs= [0.30, 0.37, 0.10, 0.29, 0.70, $
		0.23, 0.31, 0.31, 0.50, 0.32, $
		0.51, 0.57, 0.33, 0.36, 0.34]




;  collisionless
; -----------------

fruns= ['cvc3vc3h','cvc3vc3b','cvc3vc3c','cvc3vc3d','cvc3vc3e', $
        'cvc3vc3f','cvc3vc3g','cvc3vc3i','cvc3vc3j','cvc3vc3k', $
        'cvc3vc3l','cvc3vc3m','cvc3vc3n','cvc3vc3o','cvc3vc3p']
ddir= intarr(15)
ddir(*)= 2



Cre_avgs=     fload_avg_re(fruns,ddir)
Ces_avgs=     fload_avg_es(fruns,ddir)
Cvmax_avgs=   fload_avg_vmax(fruns,ddir)
Csig_avgs=    fload_avg_sig(fruns,ddir)
Cvs_avgs=     fload_avg_vs(fruns,ddir)
Cvsst_avgs=   fload_avg_vsst(fruns,ddir)
CTTs_avgs=    fload_Ts(fruns)
Cdeltas_avgs= fload_deltas(fruns)
Cfmin_avgs= [0.59, 0.64, 0.52, 0.93, 0.78,$ 
		0.93, 0.70, 0.90, 0.77, 0.87, $
		0.82, 0.94, 0.91, 0.75, 0.94]



namelabels= ['!6h','b','c','d','e','f','g','i','j','k','l','m','n','o','p']

;------------------------------------------------------------------------------------
;xaxistitle='V!Dmaj!N/!7r!6'
;xmax = 1.2
;xmin= 0.0
;xvar= vs_avgs
;Cxvar= Cvs_avgs
;------------------------------------------------------------------------------------
xaxistitle='(V!Dmaj!N/!7r!6)!E*!N'
xmax = 1.8
xmin= 0.0
xvar= vsst_avgs
Cxvar= Cvsst_avgs
;------------------------------------------------------------------------------------
;xaxistitle='Triaxiality, T'
;xmax = 1.2
;xmin= 0.0
;xvar= TTs_avgs
;Cxvar= CTTs_avgs
;------------------------------------------------------------------------------------
yaxistitle='Triaxiality, T'
ymax = 0.9
ymin= 0.0
yvar= TTs_avgs
Cyvar= CTTs_avgs
;------------------------------------------------------------------------------------
;yaxistitle='!7r!3'
;ymax = 190.0
;ymin= 120.0
;yvar= sig_avgs
;Cyvar= Csig_avgs
;------------------------------------------------------------------------------------
;yaxistitle='Anisotropy, !7d!3'
;ymax = 0.85
;ymin= -0.05
;yvar= deltas_avgs
;Cyvar= Cdeltas_avgs
;------------------------------------------------------------------------------------
;xaxistitle='Anisotropy, !7d!3'
;xmax = 1.0
;xmin= 0.0
;xvar= deltas_avgs
;Cxvar= Cdeltas_avgs
;------------------------------------------------------------------------------------
;yaxistitle='<!7e!3>'
;ymax = 0.9
;ymin= 0.0
;yvar= es_avgs
;Cyvar= Ces_avgs
;------------------------------------------------------------------------------------
;yaxistitle='f!Dmin!N'
;ymax = 1.0
;ymin= 0.0
;yvar= fmin_avgs
;Cyvar= Cfmin_avgs
;------------------------------------------------------------------------------------



; Print thie mess up
; -------------------
filename='something.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12, newysize=12

x0= 0.15 & x1= 0.98
y0= 0.15 & y1= 0.98



!p.position= [x0, y0, x1, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase



;oplot, xvar, yvar, psym=2, color=150, thick= 3.0
symsize= 2.0
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
oplot, xvar, yvar, psym=8, color=150, thick= 3.0
for i=0,n_elements(namelabels)-1 do begin
	;xyouts, xvar[i], yvar[i], namelabels[i], /data, color=150, charthick=3.0
	;xyouts, xvar[i]-0.01, yvar[i]-1.0, namelabels[i], /data, color=1, charthick=3.0, size=1.2   ; sigma-T
	xyouts, xvar[i]-0.01, yvar[i]-0.01, namelabels[i], /data, color=1, charthick=3.0, size=1.1   ; e-T
endfor

;oplot, Cxvar, Cyvar, psym=5, color=50, thick= 3.0
for i=0,n_elements(namelabels)-1 do begin
	xyouts, Cxvar[i], Cyvar[i], namelabels[i], /data, color=50, charthick=3.0, size=1.2
endfor




;------------------------------------------------------------------------------------


; Print thie mess up
; -------------------
;filename='something.eps'

;initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12, newysize=22

;x0= 0.15 & x1= 0.99
;y0= 0.07 & ysize=0.23
;y1= y0+ysize
;y2= y0+ysize+ysize
;y3= y0+ysize+ysize+ysize
;y4= y0+ysize+ysize+ysize+ysize


;  Ecc
; -------------
;yaxistitle='<!7e!N>'
;ymax = 0.9
;ymin= 0.0
;!p.position= [x0, y3, x1, y4]
;plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
;        xrange=[xmin,xmax], yrange=[ymin,ymax], $
;        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
;        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

;oplot, xvar, es_avgs, psym=-2, color=150, thick= 3.0
;if secondgrp eq 1 then oplot, xvar, es_avgs2, psym=-5, color=50, thick= 3.0



; first: V_maj/sigma
; -------------------
;yaxistitle='<V!Dmaj!N/!7r!3>'
;ymax = 0.9
;ymin= 0.0
;!p.position= [x0, y2, x1, y3]
;plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
;        xrange=[xmin,xmax], yrange=[ymin,ymax], $
;        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
;        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

;oplot, xvar, vs_avgs, psym=-2, color=150, thick= 3.0
;if secondgrp eq 1 then oplot, xvar, vs_avgs2, psym=-5, color=50, thick= 3.0


; second: Triaxiality
; ---------------------
;yaxistitle='Triaxiality, T'
;ymax = 1.1
;ymin= 0.0
;!p.position= [x0, y1, x1, y2]
;plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
;        xrange=[xmin,xmax], yrange=[ymin,ymax], $
;        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
;        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

;oplot, xvar, TTs_avgs, psym=-2, color=150, thick= 3.0
;if secondgrp eq 1 then oplot, xvar, TTs_avgs2, psym=-5, color=50, thick= 3.0



; second: delta
; ----------------
;yaxistitle='Anisotropy, !7d!3'
;ymax = 0.9
;ymin= 0.2
;!p.position= [x0, y0, x1, y1]
;plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
;        xrange=[xmin,xmax], yrange=[ymin,ymax], $
;        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
;        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase

;oplot, xvar, deltas_avgs, psym=-2, color=150, thick= 3.0
;if secondgrp eq 1 then oplot, xvar, deltas_avgs2, psym=-5, color=50, thick= 3.0



; -----------------------------------------------------------------------------



xyouts, 0.67, 0.85, '40% gas', /normal, charthick=3.0, size=1.5, color=150
xyouts, 0.67, 0.79, 'dissipationless', /normal, charthick=3.0, size=1.5, color=50
;xyouts, 0.67, 0.35, '40% gas', /normal, charthick=3.0, size=1.5, color=150
;xyouts, 0.67, 0.29, 'dissipationless', /normal, charthick=3.0, size=1.5, color=50
;xyouts, 0.24, 0.26, '40% gas', /normal, charthick=3.0, size=1.5, color=150
;xyouts, 0.24, 0.20, 'dissipationless', /normal, charthick=3.0, size=1.5, color=50


;xyouts, 0.65, 0.95, fload_fid(frun), /normal, charthick=3.0, size=1.33, color=0
;xyouts, 0.65, 0.91, 'xy projection', /normal, charthick=3.0, size=1.33, color=0



; -------------
;  Done
; -------------

device, /close



end























;================================================================================
;
;
;     V/sigma  Data
;  -------------------
;
;
;
;
;
;
;================================================================================





;========================================
;
;   Read the Bender data
;
;========================================

; This is the Bender (1988) data:
; 
; 
;
pro readandplot_bender_data, junk, addkey=addkey

; this first file has spaces in the name which 
; screws things up, use  the second one.
;osullivanfile= '/home/tcox/osullivan_01_all.txt'


; --------

benderfile= '/home/tcox/RotvAnisSupport/bender88_ellipticals.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=2, symsize=1.0, color= 0
; --------


benderfile= '/home/tcox/RotvAnisSupport/bender90_dwarfe_ellipticals.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=1, symsize=1.0, color= 0
; --------


benderfile= '/home/tcox/RotvAnisSupport/bender90_lowlum_ellipticals.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=1, symsize=1.0, color= 0
; --------


;---------------------------------------------------------
; normal coordinates
x0=0.25
y0=0.93

; data coordinates
x0d= 0.05
y0d= 1.58

if keyword_set(addkey) then begin
	xyouts, x0+.01, y0-0.07, 'Bender (1988)', color= 0, size=0.8, charthick=2.0, /normal
	oplot, [x0d], [y0d-0.06], psym=2, color= 0

	xyouts, x0+.01, y0-0.11, 'Bender & Nieto (1990)', color= 0, size=0.8, charthick=2.0, /normal
	oplot, [x0d], [y0d-0.14], psym=1, color= 0
endif


;---------------------------------------------------------





end






;========================================
;
;   Read the Davies data
;
;========================================

; This is the Davies et al. (1983) data:
; 
; 
;
pro readandplot_davies_data, junk, addkey=addkey

; this has the exact same file format
; as the bender data, so we'll use the 
; identical procedure
;
;


benderfile= '/home/tcox/RotvAnisSupport/davies83_ellipticals.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma

; low luminosity ellipticals
symsize= 1.0
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
idx= where(mb gt -20.5)
oplot, ellip(idx), vmax(idx)/sigma(idx), psym=8, color= 0
; --------

; high luminosity ellipticals
symsize= 0.8
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
idx= where(mb le -20.5)
oplot, ellip(idx), vmax(idx)/sigma(idx), psym=8, color= 0
; --------


; bulges
benderfile= '/home/tcox/RotvAnisSupport/davies83_bulges.txt'
read_bender_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=7, symsize=0.8, color= 0, thick=2.0
; --------


; ----------------------------------------------------------
; normal coordinates
x0=0.25
y0=0.93

; data coordinates
x0d= 0.05
y0d= 1.61

if keyword_set(addkey) then begin
	xyouts, x0+.01, y0-0.03, 'Davies et al. (1983)', color= 0, size=0.8, charthick=2.0, /normal
	symsize=1.0
	usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	oplot, [x0d-0.025], [y0d], psym=8, color= 0
	xyouts, x0-.0615, y0-0.03, ',', color= 0, size=0.8, charthick=2.0, /normal
	symsize= 0.8
	usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
	oplot, [x0d], [y0d], psym=8, color= 0
	xyouts, x0-.0325, y0-0.03, ',', color= 0, size=0.8, charthick=2.0, /normal
	oplot, [x0d+0.025], [y0d], psym=7, color= 0,thick=2.0
endif
; ----------------------------------------------------------



end







;================================================================================







pro read_bender_data, filename, ellip, mb, vmax, sigma, $
				no_neg_trap=no_neg_trap


;
;  Format of this file is:
; -------------------------
;#
;#
;#
;# 		Ellip	M_b		V_max	Sigma
;#					(km/s)	(km/s)
;NGC205		0.40	-15.7		20	50
;M32		0.20	-15.5		15	60
;Fornax		0.33	-12.6		1.7	6
; etc...
;

spawn, "wc "+filename,result
lines=long(result)
datalines=lines(0)-5
ellip= fltarr(datalines)
mb= fltarr(datalines)
vmax= fltarr(datalines)
sigma= fltarr(datalines)

openr, 1, filename
print, "opening: ", filename
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
        readf, 1, junk
        ;print, junk
        tempjunk= strsplit(junk,/extract,count=count)
        name= tempjunk(0)
        ellip(i)= float(tempjunk(1))
        mb(i)= float(tempjunk(2))
        vmax(i)= float(tempjunk(3))
        sigma(i)= float(tempjunk(4))
endfor

close, 1


; trap for -1 values
if not keyword_set(no_neg_trap) then begin
   idx= where((vmax lt 0.0) or (sigma lt 0.0))
   if idx(0) ne -1 then begin
	vmax(idx)= -1.0
	sigma(idx)= 1.0
	print,"Number of data points: ", n_elements(ellip)
	print,"Fixing ", n_elements(idx), " for negative values."
   endif
endif


end







;================================================================================






pro read_bender_5, filename, ellip, mb, vmax_min, vmax_maj, sigma, $
				no_neg_trap=no_neg_trap


;
;  Format of this file is:
; -------------------------
;#   5x of these
;NGC205		0.40	-15.7		20	50    0.0
; etc...
;

spawn, "wc "+filename,result
lines=long(result)
datalines=lines(0)-5
ellip= fltarr(datalines)
mb= fltarr(datalines)
vmax_maj= fltarr(datalines)
vmax_min= fltarr(datalines)
sigma= fltarr(datalines)

openr, 1, filename
print, "opening: ", filename
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
        readf, 1, junk
        ;print, junk
        tempjunk= strsplit(junk,/extract,count=count)
        name= tempjunk(0)
        ellip(i)= float(tempjunk(1))
        mb(i)= float(tempjunk(2))
        vmax_min(i)= float(tempjunk(3))
        vmax_maj(i)= float(tempjunk(4))
        sigma(i)= float(tempjunk(5))
endfor

close, 1


; trap for -1 values
if not keyword_set(no_neg_trap) then begin
   idx= where((vmax lt 0.0) or (sigma lt 0.0))
   if idx(0) ne -1 then begin
	vmax(idx)= -1.0
	sigma(idx)= 1.0
	print,"Number of data points: ", n_elements(ellip)
	print,"Fixing ", n_elements(idx), " for negative values."
   endif
endif


end







;================================================================================








;========================================
;
;   Read the deZeeuw data
;
;========================================

; This is the de Zeeuw et al. (2002) SAURON data:
; 
; 
;
pro readandplot_deZeeuw_data, junk, addkey=addkey

; this first file has spaces in the name which 
; screws things up, use  the second one.
;osullivanfile= '/home/tcox/osullivan_01_all.txt'

benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_cluster_es.txt'
read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

;benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_cluster_len.txt'
;read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
;oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

;benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_cluster_sp.txt'
;read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
;oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0


benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_field_es.txt'
read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

;benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_field_len.txt'
;read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
;oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

;benderfile= '/home/tcox/RotvAnisSupport/deZeeuw02_field_sp.txt'
;read_deZeeuw_data, benderfile, ellip, mb, vmax, sigma
;oplot, ellip, vmax/sigma, psym=5, symsize=1.0, color= 0

; normal coordinates
x0=0.25
y0=0.93

; data coordinates
x0d= 0.05
y0d= 1.58

if keyword_set(addkey) then begin
	xyouts, x0+.01, y0-0.15, 'de Zeeuw (2002)', color= 0, size=0.8, charthick=2.0, /normal
	oplot, [x0d], [y0d-0.23], psym=5, color= 0
endif


;---------------------------------------------------------

end





;
; ----------------------------------------------------------
pro read_deZeeuw_data, filename, ellip, mb, vmax, sigma


;
;  Format of this file is:
; -------------------------
;#  "Field"
;#  Elliptical, Lenticular, and Spiral data from de Zeeuw et al. (2002) - SAURON II
;#
;#Galaxy Type            T-type  Vsys    d_m     M_b     (B-V)   Mg2     R_e     e       mu_e    V_max   Sigma   Gamma   M_BH
;#                               (km/s)  (mag)   (mag)   (mag)   mag     (arcsec)        (mag)   (km/s)  (km/s)          a-b(10^c) Msol
;Ellipticals
;NGC821  E6?             4.2     1742    31.86   20.44   1.020   0.316   50      0.32    22.02   91      208     -1      3.0-7.0(7)
;NGC2699 E:              5.0     1825    31.83   18.85   0.980   0.282   -1      0.06    -1      -1      -1      -1      -1 
;NGC2768 E6:             3.1     1324    31.66   21.15   0.960   0.276   64      0.42    21.94   148     188     -1      -1 
;NGC2974 E4              3.6     1983    31.93   20.32   1.005   0.305   24      0.39    20.74   207     229     -1      -1 
; etc...
;#
;#
;

spawn, "wc "+filename,result
lines=long(result)
datalines=lines(0)-5
ellip= fltarr(datalines)
mb= fltarr(datalines)
vmax= fltarr(datalines)
sigma= fltarr(datalines)

openr, 1, filename
print, "opening: ", filename
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
        readf, 1, junk
        ;print, junk
        tempjunk= strsplit(junk,/extract,count=count)
        name= tempjunk(0)
        ellip(i)= float(tempjunk(9))
        mb(i)= float(tempjunk(5))
        vmax(i)= float(tempjunk(11))
        sigma(i)= float(tempjunk(12))
endfor

close, 1


end








;
; ----------------------------------------------------------
pro read_bbf_data, filename, mb, vsigst, a4diva


;
; # E's, cE's, dE's and bulges with known internal kinematics and/or Mg_2
; # ----------------------------------------------------------------------
; #
; # obj       l     b    Type     D   S  lgsig  S  lg re  SB_e     M_T   S  A_B   lg vds*   Mg2  S  B-Vo  S   a4   lg rc  R.C.
; # (1)      (2)   (3)    (4)    (5) (6)  (7)  (8)  (9)   (10)    (11) (12) (13)   (14)    (15)(16) (17)(18) (19)   (20)  (21)
; N0315   124.6 -32.5   LA    107.2 1  2.546  1  1.486  22.36  -23.61  1  0.26  -1.046   0.283 1  1.00  1 -0.31   9.999  3 - giant
; N0584   149.8 -67.6   E4     39.6 1  2.337  1  0.724  20.44  -21.72  1  0.13   0.190   0.283 1  0.93  1  1.50  -0.609  3   E's
;
;

spawn, "wc "+filename,result
lines=long(result)
;datalines=lines(0)-5
datalines=114               ; set this manually
ellip= fltarr(datalines)
mb= fltarr(datalines)
vmax= fltarr(datalines)
sigma= fltarr(datalines)
vsigst= fltarr(datalines)
a4diva= fltarr(datalines)

openr, 1, filename
print, "opening: ", filename
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
        readf, 1, junk
        ;print, junk
        tempjunk= strsplit(junk,/extract,count=count)
        name= tempjunk(0)

	if name eq 'stop' then break

        ;ellip(i)= float(tempjunk(x))
        mb(i)= float(tempjunk(10))
        ;vmax(i)= float(tempjunk(x))
        ;sigma(i)= float(tempjunk(x))
	vsigst(i)= 10^(float(tempjunk(13)))
	a4diva(i)= float(tempjunk(18))
endfor

close, 1


end







;
; ----------------------------------------------------------
pro read_bender_current, name, ellip, a4a, vmaj, vmin, sig, pK, dPA, M_t, SB_e, L_R, L_X



filename='/home/tcox/RotvAnisSupport/bender_current_data.txt'
headerlen= 6
;#  TABLE : e_kinxr.dat
;#  upper limits for velocities are given by negative values,
;#  no data available or a4 not classifiable: 999; if dPA = 999: galaxy
;#  re-classified as `face-on' SB0, if a4>4 then a4:=4;  pK = pec.core kinem., strong mi.ax.rot (no SB0)
;#   NGC 1-b/a   a4/a  Vmaj Vmin <sig> pK dPA   M_t   SB_e   L_R    L_X
;#                     km/s km/s km/s     deg   mag  mag/""  W/Hz   erg/s
;  N067A 0.53   3.70  999  999   999  0   1  -99.99 -99.99  99.99  99.99
;  N227  0.36   1.30  999  999   999  0   8  -21.90  21.33  99.99  99.99
;  N315  0.25  -0.31   30  -30   290  0   2  -23.61  22.36  24.31  41.93


spawn, "wc "+filename,result
lines=long(result)
datalines=lines(0)-headerlen
;datalines=114               ; set this manually

name= strarr(datalines)

ellip= fltarr(datalines)     ; 1
a4a= fltarr(datalines)       ; 2
vmaj= fltarr(datalines)      ; 3
vmin= fltarr(datalines)      ; 4
sig= fltarr(datalines)       ; 5
pK= fltarr(datalines)        ; 6
dPA= fltarr(datalines)       ; 7
M_t= fltarr(datalines)       ; 8
SB_e= fltarr(datalines)      ; 9
L_R= fltarr(datalines)       ; 10
L_X= fltarr(datalines)       ; 11

openr, 1, filename
print, "opening: ", filename
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,(lines(0)-headerlen-1) do begin
        readf, 1, junk
        ;print, junk
        tempjunk= strsplit(junk,/extract,count=count)

        name(i)= tempjunk(0)
	if name(i) eq 'stop' then break

        ellip(i)= float(tempjunk(1))
        a4a(i)= float(tempjunk(2))
        vmaj(i)= float(tempjunk(3))
        vmin(i)= float(tempjunk(4))
        sig(i)= float(tempjunk(5))
	pK(i)= float(tempjunk(6))
	dPA(i)= float(tempjunk(7))
	M_t(i)= float(tempjunk(8))
	SB_e(i)= float(tempjunk(9))
	L_R(i)= float(tempjunk(10))
	L_X(i)= float(tempjunk(11))
endfor

close, 1


end







