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
spawn, '/bin/ls /raid4/sdutta/projections/'+datadir+'/'+nfrun+'/'+nfrun+'*.prof ',result

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
;     Separate Map of ellip. versus disky/boxy
;       for each orbit
;   ---------------------------------------------
;
;
;------------------------------------------------------------------------

pro dichmap, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "dichmap, junk"
   print, "  "
   return
endif

filename='dichmap.eps'


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
	process_one_dichmap, 'vc3vc3h', x0, y2, x1, y3, ddir=ddir, /showyaxis, msg='h', secmsg='(pro-pro)', /showdelta
	process_one_dichmap, 'vc3vc3b', x1, y2, x2, y3, ddir=ddir, msg='b', secmsg='(pro-ret)'
	process_one_dichmap, 'vc3vc3c', x2, y2, x3, y3, ddir=ddir, msg='c', secmsg='(ret-ret)'
	process_one_dichmap, 'vc3vc3d', x3, y2, x4, y3, ddir=ddir, msg='d', secmsg='(polar1)'
	process_one_dichmap, 'vc3vc3e', x4, y2, x5, y3, ddir=ddir, msg='e', secmsg='(tilted1)'
	process_one_dichmap, 'vc3vc3f', x0, y1, x1, y2, ddir=ddir, /showyaxis, msg='f', secmsg='(polar2)'
	process_one_dichmap, 'vc3vc3g', x1, y1, x2, y2, ddir=ddir, msg='g', secmsg='(tilted2)'
	process_one_dichmap, 'vc3vc3i', x2, y1, x3, y2, ddir=ddir, msg='i'
	process_one_dichmap, 'vc3vc3j', x3, y1, x4, y2, ddir=ddir, msg='j'
	process_one_dichmap, 'vc3vc3k', x4, y1, x5, y2, ddir=ddir, msg='k'
	process_one_dichmap, 'vc3vc3l', x0, y0, x1, y1, ddir=ddir, /showyaxis, /showxaxis, msg='l'
	process_one_dichmap, 'vc3vc3m', x1, y0, x2, y1, ddir=ddir, /showxaxis, msg='m'
	process_one_dichmap, 'vc3vc3n', x2, y0, x3, y1, ddir=ddir, /showxaxis, msg='n'
	process_one_dichmap, 'vc3vc3o', x3, y0, x4, y1, ddir=ddir, /showxaxis, msg='o'
	process_one_dichmap, 'vc3vc3p', x4, y0, x5, y1, ddir=ddir, /showxaxis, msg='p'
endif



;   *collisionless*
;  orientation vc3's
; ---------------------------
do_orient= 0
;do_orient= 1
if do_orient eq 1 then begin
	process_one_dichmap, 'cvc3vc3h', x0, y2, x1, y3, ddir=2, /showyaxis, msg='h', secmsg='(pro-pro)', /showdelta
	process_one_dichmap, 'cvc3vc3b', x1, y2, x2, y3, ddir=2, msg='b', secmsg='(pro-ret)'
	process_one_dichmap, 'cvc3vc3c', x2, y2, x3, y3, ddir=2, msg='c', secmsg='(ret-ret)'
	process_one_dichmap, 'cvc3vc3d', x3, y2, x4, y3, ddir=2, msg='d', secmsg='(polar1)'
	process_one_dichmap, 'cvc3vc3e', x4, y2, x5, y3, ddir=2, msg='e', secmsg='(tilted1)'
	process_one_dichmap, 'cvc3vc3f', x0, y1, x1, y2, ddir=2, /showyaxis, msg='f', secmsg='(polar2)'
	process_one_dichmap, 'cvc3vc3g', x1, y1, x2, y2, ddir=2, msg='g', secmsg='(tilted2)'
	process_one_dichmap, 'cvc3vc3i', x2, y1, x3, y2, ddir=2, msg='i'
	process_one_dichmap, 'cvc3vc3j', x3, y1, x4, y2, ddir=2, msg='j'
	process_one_dichmap, 'cvc3vc3k', x4, y1, x5, y2, ddir=2, msg='k'
	process_one_dichmap, 'cvc3vc3l', x0, y0, x1, y1, ddir=2, /showyaxis, /showxaxis, msg='l'
	process_one_dichmap, 'cvc3vc3m', x1, y0, x2, y1, ddir=2, /showxaxis, msg='m'
	process_one_dichmap, 'cvc3vc3n', x2, y0, x3, y1, ddir=2, /showxaxis, msg='n'
	process_one_dichmap, 'cvc3vc3o', x3, y0, x4, y1, ddir=2, /showxaxis, msg='o'
	process_one_dichmap, 'cvc3vc3p', x4, y0, x5, y1, ddir=2, /showxaxis, msg='p'
endif



print, "WARNING: msg's have been manually moved around - check their status"


device, /close



end





; ========================================================================















;
;
;
; ---------------------------------------------------
pro process_one_dichmap, frun, x0, y0, x1, y1, $
				showyaxis=showyaxis, showxaxis=showxaxis, $
				msg=msg, $
				secmsg=secmsg, $
				ddir=ddir, $
				showdelta=showdelta



; ellipticity
yaxistitle= 'ellipticity, !7e!3'
ymax= 0.7
ymin= 0.0

; disky/boxy
xaxistitle= '!6100 a!D4!N/a'
xmax = 4.2
xmin =-2.0


bins= 40


;read_vsig_file, frun, phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir
;ell=ellipticity

if strmid(frun,0,1) eq 'c' then frun='/raid4/tcox/collisionless/'+frun else frun='/raid4/tcox/'+frun
read_diskyboxy, frun, aa, bb, ellip, a4diva
db= 100.0 * a4diva
ell= ellip

; ------------------
; compute histogram
; ------------------

contour_makegeneralpic, db, ell, xmax, xmin, ymax, ymin, $
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

;tnames=['0',' ','0.2',' ','0.4',' ','0.6',' ']

if keyword_set(showyaxis) and keyword_set(showxaxis) then begin
   plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
	xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
	;xticks= 7, xtickname=tnames, $
	xtitle=xaxistitle, ytitle=yaxistitle
endif else begin
   if keyword_set(showyaxis) then begin
        plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
	;xticks= 7, $
        ytitle=yaxistitle, xtickformat='(a1)'
   endif

   if keyword_set(showxaxis) then begin
        plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
	;xticks= 7, xtickname=tnames, $
        xtitle=xaxistitle, ytickformat='(a1)'
   endif

   if ((not keyword_set(showyaxis)) and (not keyword_set(showxaxis))) then begin
	plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
	;xticks= 7, $
        xtickformat='(a1)', ytickformat='(a1)'
   endif

endelse


;xyouts, 0.22, 0.90, /normal, fload_timelbl(1), size=1.0, color=0
;xyouts, 0.65, 0.90, /normal, fload_fid(1), size=1.33, color=0

if keyword_set(msg) then begin
   xyouts, 0.9, 0.62, /data, msg, size=1.2, charthick= 3.0, color=0    ; default one
endif

if keyword_set(secmsg) then begin
   xyouts, 0.5, 0.53, /data, secmsg, size=1.1, color=0           ; default
endif


; --------------------
; line for disky/boxy
x=[0,0]
y=[-0.5,1.5]
oplot, x, y, psym=-3, color=0, thick=3.0, linestyle= 2



end






; =======================================================================================








