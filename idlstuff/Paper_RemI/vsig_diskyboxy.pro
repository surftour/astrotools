pro diskyboxy, frun, snapnum, xlen


;frun= '/raid4/tcox/vc3vc3b'
;snapnum= 30
;xlen= 10.0



; open the file up
;   and load variables
; ---------------------
ok= fload_snapshot_bh(frun, snapnum)
center=fload_center_alreadycomp(1)

ndisk= fload_npart(2)                               ; disk particles
nbulge= fload_npart(3)                              ; bulge particles
nstars= fload_npart(4)
npart= long(ndisk) + long(nbulge) + long(nstars)
print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars

x= fload_allstars_xyz('x',center=[0,0,0])
y= fload_allstars_xyz('y',center=[0,0,0])
z= fload_allstars_xyz('z',center=[0,0,0])
m= fload_allstars_mass(1)
hsml= fload_allstars_hsml(1)

hsml= hsml*2.0
print, "WARNING: multiplying the stellar hsml by 2.0"


; --------------------------------------------------------

;  read the list of theta's and phi's

anglefile= '/home/tcox/unitsphereangles.txt'

spawn, "wc "+anglefile,result
lines=long(result)
Nangles=lines(0)-5
if Nangles GT 0 then angle_data= fltarr(2,Nangles)

openr, 1, anglefile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, angle_data
close, 1


theta= angle_data[0,*]
phi= angle_data[1,*]


; --------------------------------------------------------
;
;   Loop through the angles, project and calculate
;  the disky/boxy - ness
;

fit_a= fltarr(Nangles)
fit_b= fltarr(Nangles)
ellip= fltarr(Nangles)
a4diva= fltarr(Nangles)



for i= 0, Nangles-1 do begin

	rotate_phi= phi[i]
	rotate_theta= theta[i]


	; ------------------------------------------------------------------

	set_maxden=10.0
	set_dynrng= 1.0e+5

	istring= strcompress(string(i), /remove_all)
	fitsfilename= frun+'/fits/massdensity'+istring+'.eps'    ; kind of funky becuase
							; contour_makepic assumes
							; the ending is eps

        contour_makepic, x, y, z, m, xlen, $
                        hsml=hsml, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        center=center, filename=fitsfilename, /fitstoo, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage, HalfMassSB=HalfMassSB

        Pic= NxNImage
        PIXELS= n_elements(Pic[0,*])



		; do some slight processing
		; --------------------------
		idx=where(Pic eq 1)
                Pic= 256-Pic
                if idx(0) ne -1 then Pic(idx)= 1
                print, "before: min/max= ",min(Pic), max(Pic)
                print, "n > 240= ", n_elements(where(Pic gt 240))
                print, "n < 20= ", n_elements(where(Pic lt 20))
                ; smooth?
                width= 10      ; try other things?
                newMap= smooth(Pic, width, /edge_truncate)
                Pic= newMap
                print, "after: min/max= ",min(Pic), max(Pic)



	; ------------------------------------------------------------------



        ; --------------------------------
        ;
        ;  Fit isophote & calc deviations
        ;
        ; --------------------------------
        calc_isophote_boxydisky= 1
        if calc_isophote_boxydisky eq 1 then begin

            Thresh= 256-HalfMassSB
            SurfaceBrightness= Pic
            Scale= 2.0*xlen

            Xfit=fltarr(360)
            YFit=fltarr(360)

            Xcon=fltarr(360)
            Ycon=fltarr(360)

            AFit=0.0
            BFit=0.0
            X0Fit=0.0
            Y0Fit=0.0
            PhiFit=0.0
            a4Fit=0.0


            S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/ComputeIsoPhotFitIDL/isophotfit.so', $
                      'fit_isophot', $
                      PIXELS, $
                      SurfaceBrightness,$
                      Thresh,$
                      Scale, $
                      Xfit, Yfit,$
                      PhiFit, AFit, BFit, X0Fit, Y0Fit, a4Fit, $
                      Xcon, Ycon)

           ; disky/boxy
           print, "a4/a= ", a4Fit/AFit
	   fit_a[i]= AFit
	   fit_b[i]= BFit
	   if AFit lt BFit then ellip[i]= 1-(AFit/BFit) else ellip[i]= 1-(BFit/AFit)
	   a4diva[i]= a4Fit/AFit


	endif


	; ------------------------------------------------------------------

endfor


; ------------------------------------------------------------------

openw, 1, frun+'/diskyboxy.txt', ERROR=err
        
printf, 1, "#   diskyboxy.txt"
printf, 1, "#   "
printf, 1, "#   "
printf, 1, "#  theta      phi        a       b   "
printf, 1, "#  (deg)    (deg)   (kpc/h) (kpc/h)     ellip     a4/a  "

for i=0,Nangles-1 do begin
	printf, 1, FORMAT= '(2(F8.3," "),4(F8.5,"  "))', $
                theta[i], phi[i], fit_a[i], fit_b[i], ellip[i], a4diva[i]
endfor
close, 1




; ------------------------------------------------------------------



end









; =================================================================


pro read_diskyboxy, frun, aa, bb, ellip, a4diva

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


aa= db_data[2,*]
bb= db_data[3,*]
ellip= db_data[4,*]
a4diva= db_data[5,*]

end



; =================================================================



;   ellipticities
; --------------------------
function fload_ellips, fruns, all=all, eserr=eserr
   ns= n_elements(fruns)
   alles= fltarr(ns)
   eserr= fltarr(ns)
   es= [-1]
   for i=0,ns-1 do begin
        read_diskyboxy, fruns[i], aa, bb, ellip, a4diva
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
        read_diskyboxy, fruns[i], aa, bb, ellip, a4diva
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
        read_diskyboxy, fruns[i], aa, bb, ellip, a4diva
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
        read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir_major
        vmaj= vmax
        read_vsig_file, fruns[i], phi, theta, Re, PA, ellipticity, vmax, sigma, ddir=ddir_minor
        vmin= vmax
        ;this_mus= abs(atan(vmin/vmaj))*180.0/!PI
        ;print, "mean kinematic misalignment, <Phi>= ",mean(this_mus)
        this_mus= vmin / sqrt(vmin*vmin + vmaj*vmaj)
        print, "                              <mu>= ",mean(this_mus)
        allmus[i]= mean(this_mus)
        mus= [mus,transpose(this_mus)]
   endfor
   mus= mus[1:n_elements(mus)-1]
   if keyword_set(all) then return, mus
   return, allmus
end






; =================================================================
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

pro db_map1, junk


if not keyword_set(junk) then begin
   print, "  "
   print, "db_map1, junk"
   print, "  "
   return 
endif

filename='dbmap1.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=3


;--------------------------------------
;--------------------------------------

; -----------------
;  std's
; -----------------
;do_regs= 1
do_regs= 0
if do_regs eq 1 then begin
	fruns= ['vc3vc3h','vc3vc3b','vc3vc3c','vc3vc3d','vc3vc3e', $
	        'vc3vc3f','vc3vc3g','vc3vc3i','vc3vc3j','vc3vc3k', $
	        'vc3vc3l','vc3vc3m','vc3vc3n','vc3vc3o','vc3vc3p']

	ddir= intarr(15)
	;ddir(*)= 0
	ddir(*)= 13

	allvsst=   fload_avg_vsst(fruns,ddir,/all)
	allvsigst_log= alog10(allvsst)

	ddir(*)= 0
	allvmaj=   fload_avg_vmax(fruns,ddir,/all)

	ddir(*)= 1
	allvmin=   fload_avg_vmax(fruns,ddir,/all)
	allvmin_div_allvtot= allvmin / sqrt(allvmin*allvmin + allvmaj*allvmaj)

	;allvmin_div_allvtot= fload_mus(fruns, /all)


	for i=0, n_elements(fruns)-1 do fruns[i]= '/raid4/tcox/'+fruns[i]

	alla4diva=   100.0 * fload_a4diva(fruns,/all)
endif



; -----------------
;  collisionless
; -----------------
do_cless= 1
;do_cless= 0
if do_cless eq 1 then begin
	fruns= ['cvc3vc3h','cvc3vc3b','cvc3vc3c','cvc3vc3d','cvc3vc3e', $
	        'cvc3vc3f','cvc3vc3g','cvc3vc3i','cvc3vc3j','cvc3vc3k', $
	        'cvc3vc3l','cvc3vc3m','cvc3vc3n','cvc3vc3o','cvc3vc3p']

	ddir= intarr(15)
	ddir(*)= 2

	allvsst=   fload_avg_vsst(fruns,ddir,/all)
	allvsigst_log= alog10(allvsst)

	allvmaj=   fload_avg_vmax(fruns,ddir,/all)

	ddir(*)= 3
	allvmin=   fload_avg_vmax(fruns,ddir,/all)

	allvtot= sqrt(allvmin*allvmin + allvmaj*allvmaj)
	allvmin_div_allvtot= allvmin / allvtot

	print, " --------------------------- "
	print, "<v_maj>=",mean(allvmaj), "  +/-", sqrt(variance(allvmaj))
	print, "<v_min>=",mean(allvmin), "  +/-", sqrt(variance(allvmin))
	print, "<v_tot>=",mean(allvtot), "  +/-", sqrt(variance(allvtot))
	print, "<f_min>=",mean(allvmin_div_allvtot), "  +/-", sqrt(variance(allvmin_div_allvtot))



	;allvmin_div_allvtot= fload_mus(fruns, /all, /cless)
	;testvar= allvmin_div_allvtot
	;print, min(testvar), max(testvar), mean(testvar), median(testvar)

	for i=0, n_elements(fruns)-1 do fruns[i]= '/raid4/tcox/collisionless/'+fruns[i]

	alla4diva= 100.0 * fload_a4diva(fruns,/all)
endif


; --------------------



;----------------------------------------
; Plot #1   (v/sig)* versus a4/a
;----------------------------------------

xaxistitle= '!6100 a!D4!N/a'
xmax= 3.5
xmin= -2.0

yaxistitle= '!6Log (V/!7r!6)!E*!N'
ymax = 0.5
ymin = -1.6

bins= 30

contour_makegeneralpic, alla4diva, allvsigst_log, xmax, xmin, ymax, ymin, $
                                pixels= bins, $
                                NxNImage=NxNImage

; --------------------

x0= 0.15 & y0= 0.55
x1= 0.98 & y1= 0.98

tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

!p.position=[x0, y0, x1, y1]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, color= 0, /noerase, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], xtickformat='(a1)', ytitle=yaxistitle, $
        ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0


; --------------------

y= [ymin,ymax]
x= [0.0,0.0]
oplot, x, y, psym=-3,color=0,thick=4.0,linestyle=2

y= [0.0,0.0]
x= [xmin,xmax]
oplot, x, y, psym=-3,color=0,thick=4.0,linestyle=2




; ------------------------------------------------------------------


        ; --------------
        ; read bbf '92
        ;bbffile= '/home/tcox/RotvAnisSupport/bbf92_spheroids.txt'
        ;read_bbf_data, bbffile, mb, vsigst, a4diva

        ;idx= where(a4diva lt 3.5)
        ;if idx(0) ne -1 then begin
	;	a4diva= a4diva(idx)
	;	vsigst= vsigst(idx)
	;endif

	;vsigst= alog10(vsigst)

        ;print, "a4/a   (max/min) = ", max(a4diva), min(a4diva)
        ;idx= where(a4diva ge 0.0) 
        ;print, "a4/a percent >= 0.0: ", 100.0*n_elements(idx)/n_elements(a4diva)
        ;idx= where(a4diva gt 0.0) 
        ;print, "a4/a percent >  0.0: ", 100.0*n_elements(idx)/n_elements(a4diva)

	;f = (2.0*!pi/16.0)*findgen(17)
	;usersym,0.7*cos(f),0.7*sin(f),/fill
        ;oplot, a4diva, vsigst, psym=8, color=0, thick=3.0

        ;oplot, [1.05,1.18], [0.86,0.86], psym=-3, color= 0, thick=12.0, linestyle= 1
        ;xyouts, 0.67, 0.73, 'BBF (1992)', /normal, charthick=4.0, size=1.5, color=0



        ; -------------------------
        ; read bender current data
	read_bender_current, name, ellip, a4a, vmaj, vmin, sig, pK, dPA, M_t, SB_e, L_R, L_X

        ; make sure there is a major axis velocity
        idx= where(vmaj lt 800.0)
        if idx(0) ne -1 then begin
                vmaj= vmaj(idx)
                vmin= vmin(idx)
                a4a= a4a(idx)
                ellip= ellip(idx)
		sig= sig(idx)
        endif

        ; make lower limits for vmaj
        idx= where(vmaj lt 0.0)
        if idx(0) ne -1 then vmaj(idx)= -1.0 * vmaj(idx)


	vsigst= vmaj/sig/sqrt(ellip/(1-ellip))
        vsigst= alog10(vsigst)

        print, "a4/a   (max/min) = ", max(a4a), min(a4a)
        idx= where(a4a ge 0.0)
        print, "a4/a percent >= 0.0: ", 100.0*n_elements(idx)/n_elements(a4a)
        idx= where(a4a gt 0.0) 
        print, "a4/a percent >  0.0: ", 100.0*n_elements(idx)/n_elements(a4a)

        f = (2.0*!pi/16.0)*findgen(17)
        usersym,0.7*cos(f),0.7*sin(f),/fill
        oplot, a4a, vsigst, psym=8, color=0, thick=3.0

        ;oplot, [1.05,1.18], [0.86,0.86], psym=-3, color= 0, thick=12.0, linestyle= 1
        ;xyouts, 0.67, 0.73, 'BBF (1992)', /normal, charthick=4.0, size=1.5, color=0









;----------------------------------------
; Plot #2   v_min/v_tot versus a4/a
;----------------------------------------

xaxistitle= '!6100 a!D4!N/a'
xmax= 3.5
xmin= -2.0

yaxistitle= '!6V!Dmin!N/(V!Dmaj!N!E2!N+V!Dmin!N!E2!N)!E1/2!N'
ymax = 1.05
ymin = -0.05

bins= 30

contour_makegeneralpic, alla4diva, allvmin_div_allvtot, xmax, xmin, ymax, ymin, $
                                pixels= bins, $
                                NxNImage=NxNImage

; --------------------

x0= 0.15 & y0= 0.12
x1= 0.98 & y1= 0.55

tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

!p.position=[x0, y0, x1, y1]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, color= 0, /noerase, /nodata, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle=xaxistitle, ytitle=yaxistitle, $
        xcharsize=1.2, ycharsize=1.2, xthick=4.0, ythick=4.0, charthick=3.0


; --------------------

y= [ymin,ymax]
x= [0.0,0.0]
oplot, x, y, psym=-3,color=0,thick=4.0,linestyle=2



; ------------------------------------------------------------------


        ; --------------
        ; read b,s & g '94
        ;bbffile= '/home/tcox/RotvAnisSupport/bender94_minoraxis.txt'
        ;read_bender_5, bbffile, ellip, mb, vmin, vmaj, a4diva, /no_neg_trap

	;mu= vmin / sqrt(vmin*vmin + vmaj*vmaj)

        ;idx= where(a4diva lt 3.5)
        ;if idx(0) ne -1 then begin
	;	a4diva= a4diva(idx)
	;	mu= mu(idx)
	;endif

	;f = (2.0*!pi/16.0)*findgen(17)
	;usersym,0.7*cos(f),0.7*sin(f),/fill
        ;oplot, a4diva, mu, psym=8, color=0, thick=3.0



        ; -------------------------
        ; read bender current data
	read_bender_current, name, ellip, a4a, vmaj, vmin, sig, pK, dPA, M_t, SB_e, L_R, L_X

	; make sure there is a major axis velocity
	idx= where(vmaj lt 800.0)
	if idx(0) ne -1 then begin
		vmaj= vmaj(idx)
		vmin= vmin(idx)
		a4a= a4a(idx)
	endif

	; make sure there is a minor axis velocity
	idx= where(vmin lt 800.0)
	if idx(0) ne -1 then begin
		vmaj= vmaj(idx)
		vmin= vmin(idx)
		a4a= a4a(idx)
	endif

	; ditch if both velocities are lower limits
	;idx= where((vmin gt 0.0) and (vmaj gt 0.0))
	;if idx(0) ne -1 then begin
	;	vmaj= vmaj(idx)
	;	vmin= vmin(idx)
	;	a4a= a4a(idx)
	;endif

	; make lower limits for vmaj
	idx= where(vmaj lt 0.0)
	if idx(0) ne -1 then begin
	   vmaj(idx)= -1.0 * vmaj(idx)
	   ;for i=0,n_elements(idx)-1 do begin
	;	vmaj(idx(i))= randomu(seed)*10.0
	   ;endfor
	endif

	; make lower limits for vmin
	idx= where(vmin lt 0.0)
	if idx(0) ne -1 then begin
	    vmin(idx)= -1.0 * vmin(idx)
	   ;for i=0,n_elements(idx)-1 do begin
	;	vmin(idx(i))= -1.0*randomu(seed)*vmin(idx(i))
	   ;endfor
	endif

        mu= vmin / sqrt(vmin*vmin + vmaj*vmaj)

        idx= where(a4a lt 50.0)
        if idx(0) ne -1 then begin
                a4a= a4a(idx)
                mu= mu(idx)
        endif

        f = (2.0*!pi/16.0)*findgen(17)
        usersym,0.7*cos(f),0.7*sin(f),/fill
        oplot, a4a, mu, psym=8, color=0, thick=3.0




; ---------------------------------------------------------------------

xyouts, 0.60, 0.60, '40% gas', /normal, size=1.5, charthick=3.0, color= 0
;xyouts, 0.60, 0.60, 'dissipationless', /normal, size=1.5, charthick=3.0, color= 0

; --------------------

device, /close


end









