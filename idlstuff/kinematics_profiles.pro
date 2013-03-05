


pro testit, junk

;do_one_angle, "data/ds/vc3vc3i", 0.0, 0.0, snapnum=19, filename="test.eps", xlen=7.0, $
;					R_e= R_e, $
;                                        Vrot_maj=Vrot_maj, avgVrot_maj=avgVrot_maj, $
;                                        Vrot_min=Vrot_min, avgVrot_min=avgVrot_min, $
;                                        pa=pa, Sig=sig, $
;                                        fit_a=fit_a, fit_b=fit_b, $
;                                        fit_ellip=fit_ellip, a4diva=a4diva
;
;
;print, "R_e= ", R_e
;print, "Vrot_maj= ", Vrot_maj
;print, "avgVrot_maj= ", avgVrot_maj
;print, "Vrot_min= ", Vrot_min
;print, "avgVrot_min= ", avgVrot_min
;print, "pa= ", pa
;print, "Sig= ", Sig
;print, "fit_a= ", fit_a
;print, "fit_b= ", fit_b
;print, "a4diva= ", fit_ellip
;print, "a4diva= ", a4diva


end



; -----------------------------------------------------------------------------



pro fit_isophote, map, level, xlen, $
		x0, x1, y0, y1, $
		fit_a=fit_a, fit_b=fit_b, $
		fit_ellip=fit_ellip, $
		a4diva=a4diva, pa=pa
		

	Thresh= float(level)+0.1
	SurfaceBrightness= map
	help, map
	print, "max(SurfaceBrightness)= ", max(SurfaceBrightness)
	print, "min(SurfaceBrightness)= ", min(SurfaceBrightness)
	Scale= 2.0*xlen          ; allows the conversion to physical units

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

	PIXELS= long((size(map))[1])

	spawn, 'echo $HOME', result
	homedir= strcompress(result,/remove_all)
	libfile= homedir+'/Tools/C-Routines_for_IDL/ComputeIsoPhotFitIDL/isophotfit.so'
	S = CALL_EXTERNAL(libfile, $
		'fit_isophot', $
		;'main', $
		PIXELS, $
		SurfaceBrightness,$
		Thresh,$
		Scale, $
		Xfit, Yfit,$
		PhiFit, AFit, BFit, X0Fit, Y0Fit, a4Fit, $
		Xcon, Ycon)

	; disky/boxy
	print, " "
	fit_a= max([AFit,BFit])
	print, "fit_a= ", fit_a
	fit_b= min([AFit,BFit])
	print, "fit_b= ", fit_b
	fit_ellip= 1-(fit_b/fit_a)
	print, "fit_ellip= ", fit_ellip
	a4diva= a4Fit/AFit
	print, "a4/a (a4diva)= ", a4diva
	pa= PhiFit * 180.0 / !PI   + 90.0
	if AFit gt BFit then pa= pa - 90.0     ; the fit doesn't care which axis it is, just the closest one
	print, "pa= ", pa
	print, " "


	; for some reason we need axes to be reset
	!p.position=[x0,y0,x1,y1]
	!p.ticklen=0.03
	plot,[0],[0], psym=3,  xrange=[-xlen,xlen], yrange=[-xlen,xlen], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xthick=4.0, ythick=4.0, charthick=3.0, $
                      xtickformat='(a1)', ytickformat='(a1)', xstyle=1, ystyle=1, /nodata

	oplot, Xfit, Yfit, color= 0, thick= 1.0


end



; -----------------------------------------------------------------------------


;##############################################################################################
;##############################################################################################
;##############################################################################################


; -----------------------------------------------------------------------------
;
;
;   profiles:
;
;   surface brightness, and ecc, a4/a, PA versus r
;
;
;
;


pro do_isophote_profile, frun, kinnum, snapnum=snapnum, loadedsnap=loadedsnap, $
					filename=filename, $
					fitsdir=fitsdir, $
					sdfitsf=sdfitsf, vfitsf=vfitsf, dfitsf=dfitsf, $
					xlen=xlen, $
					R_e= R_e, $
					ellip_Re=ellip_Re, $
					pa_Re= pa_Re, $
					a4diva_Re= a4diva_Re, $
					ellip_5Re=ellip_5Re, $
					pa_5Re= pa_5Re, $
					a4diva_5Re= a4diva_5Re


if not keyword_set(frun) then begin
        print, " "
        print, " do_isophote_profile, frun,"
        print, " "
        ;return
endif


if not keyword_set(snapnum) then begin
        print, " "
        print, " do_isophote_profile, frun,"
        print, " "
        ;return
endif






; --------------------------
;  Setup stuff
; --------------------------

if not keyword_set(filename) then filename='profile_isophotes.eps'
if not keyword_set(kinnum) then kinnum= 0



; panel dimensions
; ------------------


x0= 0.17
x1= 0.97

y0= 0.07 & ys= 0.12
y1= y0+ys
y2= y0+ys+ys
y3= y0+ys+ys+ys    ; 0.48

ys= 0.28
y4= y3+ys
y5= y3+ys+ys
    



; ------------------
;  Set up file info
; ------------------

if not keyword_set(loadedsnap) then ok=fload_snapshot_bh(frun,snapnum, /nopot_in_snap)


; setup directory information
if not keyword_set(fitsdir) then fitsdir= frun+'/snap_'+fload_finfo_exts(1)+'_info/'

istring= strcompress(string(kinnum), /remove_all)

;if keyword_set(sdfitsf) then surfacedensity_fits_file= sdfitsf else surfacedensity_fits_file= fitsdir+'/kinmap_'+istring+'.fits'
;if keyword_set(vfitsf) then velocity_fits_file= vfitsf else velocity_fits_file= fitsdir+'/kinmap_'+istring+'_vel.fits'
;if keyword_set(dfitsf) then dispersion_fits_file= dfitsf else dispersion_fits_file= fitsdir+'/kinmap_'+istring+'_dis.fits'
if keyword_set(sdfitsf) then surfacedensity_fits_file= sdfitsf else surfacedensity_fits_file= fitsdir+'/'+istring+'_slitmap_sd.fits'
if keyword_set(vfitsf) then velocity_fits_file= vfitsf else velocity_fits_file= fitsdir+'/'+istring+'_slitmap_vel.fits'
if keyword_set(dfitsf) then dispersion_fits_file= dfitsf else dispersion_fits_file= fitsdir+'/'+istring+'_slitmap_dis.fits'


; testing
;fitsdir= "./"
;surfacedensity_fits_file= 'temp_kin.fits'
;velocity_fits_file= 'temp_kin_vel.fits'
;dispersion_fits_file= 'temp_kin_dis.fits'



; --------------
;   Plot it Up
; --------------


; get print stuff ready
; ----------------------
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 0, newxsize=12, newysize=30
print, "writing: ", filename





; read in and display surface density
; -------------------------------------
surfacedensity=readfits(surfacedensity_fits_file,header)
xlen= sxpar(header,'XLEN')   &   if xlen le 0 then xlen= 10.0
print, "xlen= ",xlen
print, header
nx= sxpar(header,'NAXIS1')
ny= sxpar(header,'NAXIS2')
HalfMassContourLevel= sxpar(header,'HALFMASS')
print, "HalfMassContourLevel= ", HalfMassContourLevel
n= nx * ny
xbin= fltarr(n)
ybin= fltarr(n)
surfacedensity_list= fltarr(n)
dx= 2. * xlen / nx
dy= 2. * xlen / ny
idx= 0L
for i=0, nx-1 do begin
  for j=0, ny-1 do begin
        xbin[idx]= dx*(i+0.5) - xlen
        ybin[idx]= dy*(j+0.5) - xlen
        surfacedensity_list[idx]= surfacedensity[i,j]
        idx= idx + 1
  endfor
endfor





; scale it to 256 and plot
; -----------------------

MaxDensXY= max(surfacedensity)
DynRange=1.0e4

ma=MaxDensXY
mi=MaxDensXY/DynRange

print, "surfacedensity     max: ", max(surfacedensity), "   min: ", min(surfacedensity)
print, "Clipping at   ma= ", ma, " mi= ", mi

; do some clipping
surfacedensity_clipped= surfacedensity
ind=where(surfacedensity lt mi)
if ind(0) ne -1 then surfacedensity_clipped(ind)=mi
ind=where(surfacedensity gt ma)
if ind(0) ne -1 then surfacedensity_clipped(ind)=ma

; number of colors
Cols=255
; map it to 0->Cols
Pic=(alog10(surfacedensity_clipped)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2
HalfMCnt=(alog10(HalfMassContourLevel)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2


; these should never occur
print, "max(Pic)= ", max(Pic)
print, "min(Pic)= ", min(Pic)
ind=where(Pic ge 256)
if ind(0) ne -1 then Pic(ind)=255
ind=where(Pic le 0)
if ind(0) ne -1 then Pic(ind)=1


; I like my colors inverted
invertcolors= 1
;invertcolors= 0
if invertcolors eq 1 then begin
   Pic_toshow=256-Pic                          ; invert color table
   HalfMCnt= 256-HalfMCnt
   idx= where(Pic EQ 254)               ; set background to white
endif else begin
   Pic_toshow= Pic
   idx= where(Pic EQ 2)
endelse

if idx(0) ne -1 then Pic_toshow(idx)= 1       ; this is setting it to white
;if idx(0) ne -1 then Pic_toshow(idx)= 0       ; this sets background to black

fload_newcolortable, 0

tv, Pic_toshow, x0, y4, xsize=(x1-x0), ysize=(y5-y4), /normal

print, "max(Pic_toshow)= ", max(Pic_toshow)
print, "min(Pic_toshow)= ", min(Pic_toshow)



; creates axes and plot style
; ---------------------------
xmax= +xlen
xmin= -xlen
ymax= +xlen
ymin= -xlen
!p.position=[x0,y4,x1,y5]
!p.ticklen=0.03
plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xthick=4.0, ythick=4.0, charthick=3.0, $
                      xtickformat='(a1)', ytickformat='(a1)', xstyle=1, ystyle=1, /nodata



show_xlen= 1
;show_xlen= 0
if keyword_set(show_xlen) then begin
        xlenlbl= strcompress(string(2.0*xlen),/remove_all)
        if (2.0*xlen) ge 1.0 then digs= 1
        if (2.0*xlen) ge 10.0 then digs= 2
        if (2.0*xlen) ge 100.0 then digs= 3
        xlenlbl = strmid(xlenlbl,0,digs)        ; T=0.x (4+digits after decimal)
        xlenlbl= '!94!6'+xlenlbl+' kpc/h !96!6'
        xyouts, x1-0.26, y5-0.02, xlenlbl, /normal, size= 1.2, charthick=3.0, color= 0
endif








;  now, calculate and display isophotal information
; ----------------------------------------------------

; smooth map slightly
;tempsd= surfacedensity
tempsd= Pic
;width= 2
;tempsd= smooth(Pic, width, /edge_truncate)
;thislevel= max(tempsd) * 0.5

levels= 100   ; every 0.2 dex

;sdfactor= 0.8     ; where to do the fits?
dlevel= (max(tempsd) - min(tempsd))/25.0

;step= (256)/levels
;userLevels= IndGen(levels) * step + step/2 
;print, "userLevels= ",userLevels
;level_to_fit= 0

; test
;userLevels= [75, 100, 120]
;levels= 7


prof_surfacedensity= fltarr(levels)
prof_ellip= fltarr(levels)
prof_pa= fltarr(levels)
prof_a4diva= fltarr(levels)
prof_radius= fltarr(levels)

for i= 0, levels-1 do begin

	;thislevel= userLevels[i]
	;thislevel= userLevels[levels-1-i]
	;thislevel= thislevel*sdfactor

        ;print, "i = ", i

	; good for pure sd
        ;if i eq 0 then high_sd= max(tempsd)*0.5 else high_sd= low_sd
        ;low_sd= high_sd * 0.9
	;prof_surfacedensity[i]= low_sd
	; good for sd scaled to color scale (linear)
        if i eq 0 then high_sd= max(tempsd) else high_sd= low_sd
        low_sd= high_sd - dlevel
	prof_surfacedensity[i]= (low_sd-2.0) / (Cols-3.0) * alog10(DynRange) + alog10(mi)

        ; if we're below the map, then get out of here
        if low_sd lt min(tempsd) then break

	thislevel= low_sd


	go_for_fit= 1
	;if thislevel gt max(tempsd) then go_for_fit= -1
	;if thislevel lt min(tempsd) then go_for_fit= -1

	;
	ind= where(tempsd gt thislevel)
	print, "i= ", i
	print, "thislevel= ", thislevel
	print, "n_elements(ind)= ", n_elements(ind)
	if n_elements(ind) lt 1500 then begin
		print, " "
		print, " *** skipping this level ***"
		print, " "
		go_for_fit= 0
	endif

	if go_for_fit eq 1 then begin
	     fit_isophote, tempsd, thislevel, xlen, $
		x0, x1, y4, y5, $
		fit_a=fit_a, fit_b=fit_b, $
		fit_ellip=fit_ellip, $
		a4diva=a4diva, pa=pa
	     prof_ellip[i]= fit_ellip
	     prof_pa[i]= pa
	     prof_a4diva[i]= a4diva
	     prof_radius[i]= 0.5 * (fit_a + fit_b)
	endif else begin
	     prof_radius[i]= -1
	endelse
		
	if (HalfMassContourLevel gt low_sd) and (HalfMassContourLevel lt high_sd) then begin
		R_e= prof_radius[i]
	endif

endfor


idx= where(prof_radius gt 0.0)
if idx(0) ne -1 then begin
	prof_radius= prof_radius(idx)
	prof_surfacedensity= prof_surfacedensity(idx)
	prof_ellip= prof_ellip(idx)
	prof_pa= prof_pa(idx)
	prof_a4diva= prof_a4diva(idx)
endif


idx= where(prof_radius lt sqrt(2)*xlen)
if idx(0) ne -1 then begin
	prof_radius= prof_radius(idx)
	prof_surfacedensity= prof_surfacedensity(idx)
	prof_ellip= prof_ellip(idx)
	prof_pa= prof_pa(idx)
	prof_a4diva= prof_a4diva(idx)
endif


idx= where(prof_surfacedensity lt alog10(HalfMassContourLevel))
if n_elements(idx) le 1 then begin
	R_e= max(prof_radius)
endif else begin
	psd= 10^(prof_surfacedensity)
	ii= idx(0)
	if ii eq 0 then ii= 1
	slp= (psd(ii-1)-psd(ii))/(prof_radius(ii-1)-prof_radius(ii))
	R_e= (HalfMassContourLevel - psd(ii))/slp + prof_radius(ii)
endelse

xaxistitle= "!6Radius (kpc/h)"
xmin= 0.0
xmax = max(prof_radius) * 1.1


;
;   surface density (it should be in log)
;
;-------------------------------------
ymax= alog10(3.0 * ma * 0.7 * 1.0e4)
ymin= alog10(0.2 * mi * 0.7 * 1.0e4)
prof_surfacedensity= prof_surfacedensity + alog10(0.7 * 1.0e4)

yaxistitle= "Log !7R!6 (M!9n!6 pc!E-2!N)"   ; & ymax = alog10(ma*3)  & ymin= ymax-5.0
!p.position= [x0, y3, x1, y4]
!p.ticklen=0.03
plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
                xstyle=1, ystyle=1, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, $
                charthick=3.0, xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, prof_radius, prof_surfacedensity, psym= -2, color= 0, thick= 1.0, symsize= 2.0

x= [R_e, R_e]
y= [ymin, ymax]
oplot, x, y, psym=-3, color= 0, thick= 1.0, linestyle= 1


;
;   ellipticity
;
;-------------------------------------
yaxistitle= "ellip."
ymax = max(prof_ellip) + 0.05
ymin= -0.05
!p.position= [x0, y2, x1, y3]
!p.ticklen=0.03
plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
                xstyle=1, ystyle=1, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, $
                charthick=3.0, xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, prof_radius, prof_ellip, psym= -2, color= 0, thick= 1.0, symsize= 2.0

x= [R_e, R_e]
y= [ymin, ymax]
oplot, x, y, psym=-3, color= 0, thick= 1.0, linestyle= 1



;
;   pa
;
;-------------------------------------
yaxistitle= "PA" & ymax = 190  & ymin= -10.0
!p.position= [x0, y1, x1, y2]
!p.ticklen=0.03
plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
                xstyle=1, ystyle=1, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, $
                charthick=3.0, xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

idx= where(prof_pa lt 0.0)
if idx(0) ne -1 then prof_pa(idx) = 180.0 + prof_pa(idx)
oplot, prof_radius, prof_pa, psym= -2, color= 0, thick= 1.0, symsize= 2.0

x= [R_e, R_e]
y= [ymin, ymax]
oplot, x, y, psym=-3, color= 0, thick= 1.0, linestyle= 1




;
;   a4/4
;
;-------------------------------------
yaxistitle= '100 a!D4!N/a'
;ymax = 3.3 
;ymin= -3.3
maxa4= max([prof_a4diva,abs(prof_a4diva)])
ymax= maxa4*105.0
ymin= -ymax
!p.position= [x0, y0, x1, y1]
!p.ticklen=0.03
plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
                xstyle=1, ystyle=1, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, $
                charthick=3.0, xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase

oplot, prof_radius, 100.0*prof_a4diva, psym= -2, color= 0, thick= 1.0, symsize= 2.0

x= [xmin, xmax]
y= [0.0, 0.0]
oplot, x, y, psym=-3, color= 0, thick= 1.0, linestyle= 2

x= [R_e, R_e]
y= [ymin, ymax]
oplot, x, y, psym=-3, color= 0, thick= 1.0, linestyle= 1



; ----------------------------------------------------
;
;
;  Calculate quantities and pass back
;    ( not currently in use -
;          now we just write a file)
;
;
;
;idx= where(prof_radius gt R_e)
;ellip_Re= prof_ellip(idx(0))
;pa_Re= prof_pa(idx(0))
;a4diva_Re= prof_a4diva(idx(0))
;
;idx= where(prof_radius gt 5*R_e)
;if idx(0) ne -1 then begin
;	ellip_5Re= prof_ellip(idx(0))
;	pa_5Re= prof_pa(idx(0))
;	a4diva_5Re= prof_a4diva(idx(0))
;endif
;
;
; ----------------------------------------------------

;fname= fitsdir+'/isophotes_'+istring+'.txt'
fname= fitsdir+'/'+istring+'_lst_isophotes.txt'

openw, 1, fname, ERROR=err

printf, 1, "#    isophote information, projection "+istring
printf, 1, "# "
printf, 1, "#          SurfDen  "
printf, 1, "# radius    Sigma    "
printf, 1, "#  (kpc)  (Mo/pc2)    Ellip.      PA   100.*a4/a"
for i=0, n_elements(prof_radius)-1 do begin
	printf, 1, FORMAT='(5(F8.4,"  "))', $
		prof_radius[i], $
		prof_surfacedensity[i], $
		prof_ellip[i], $
		prof_pa[i], $
		prof_a4diva[i]
endfor

close, 1





; ----------------------------------------------------



device, /close

; -------------
;  Done
; -------------


end






; -----------------------------------------------------------------------------
;
;
;   profiles:
;
;   surface brightness, velocity and lambda_r
;
;
;
;   two versions of this plot:
;        1. three panel (sb, vel, lambda)
;        2. four panel (sb, vel, disp, lambda)
;
;
;


;
;    3 panel version
;

pro do_lambda_profile_3, frun, kinnum, snapnum=snapnum, loadedsnap=loadedsnap, $
					filename=filename, $
					fitsdir=fitsdir, $
					sdfitsf=sdfitsf, vfitsf=vfitsf, dfitsf=dfitsf, $
					xlen=xlen, $
					lam_Re=lam_Re, $
					lam_5Re=lam_5Re


if not keyword_set(frun) then begin
        print, " "
        print, " do_isophote_profile_3, frun,"
        print, " "
        ;return
endif


if not keyword_set(snapnum) then begin
        print, " "
        print, " do_isophote_profile_3, frun,"
        print, " "
        ;return
endif






; --------------------------
;  Setup stuff
; --------------------------

if not keyword_set(filename) then filename='profile_kinematics.eps'
if not keyword_set(kinnum) then kinnum= 0



; panel dimensions
; ------------------


x0= 0.15
x1= 0.98

y0= 0.08 & ys= 0.20
y1= y0+ys

ys= 0.35
y2= y1+ys
y3= y1+ys+ys
    



; ------------------
;  Set up file info
; ------------------

if not keyword_set(loadedsnap) then ok=fload_snapshot_bh(frun,snapnum, /nopot_in_snap)


; setup directory information
if not keyword_set(fitsdir) then fitsdir= frun+'/snap_'+fload_finfo_exts(1)+'_info/'

istring= strcompress(string(kinnum), /remove_all)

;surfacedensity_fits_file= fitsdir+'/kinmap_'+istring+'.fits'
;velocity_fits_file= fitsdir+'/kinmap_'+istring+'_vel.fits'
;dispersion_fits_file= fitsdir+'/kinmap_'+istring+'_dis.fits'
;if keyword_set(sdfitsf) then surfacedensity_fits_file= sdfitsf else surfacedensity_fits_file= fitsdir+'/kinmap_'+istring+'.fits'
;if keyword_set(vfitsf) then velocity_fits_file= vfitsf else velocity_fits_file= fitsdir+'/kinmap_'+istring+'_vel.fits'
;if keyword_set(dfitsf) then dispersion_fits_file= dfitsf else dispersion_fits_file= fitsdir+'/kinmap_'+istring+'_dis.fits'
if keyword_set(sdfitsf) then surfacedensity_fits_file= sdfitsf else surfacedensity_fits_file= fitsdir+'/'+istring+'_slitmap_sd.fits'
if keyword_set(vfitsf) then velocity_fits_file= vfitsf else velocity_fits_file= fitsdir+'/'+istring+'_slitmap_vel.fits'
if keyword_set(dfitsf) then dispersion_fits_file= dfitsf else dispersion_fits_file= fitsdir+'/'+istring+'_slitmap_dis.fits'


; testing
;fitsdir= "./"
;surfacedensity_fits_file= 'temp_kin.fits'
;velocity_fits_file= 'temp_kin_vel.fits'
;dispersion_fits_file= 'temp_kin_dis.fits'



; --------------
;   Plot it Up
; --------------


; get print stuff ready
; ----------------------
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 0, newxsize=15, newysize=32
print, "writing: ", filename




;-----------------------------------------------------



; read in and display surface density
; -------------------------------------
surfacedensity=readfits(surfacedensity_fits_file,header)
xlen= sxpar(header,'XLEN')  &  if xlen le 0 then xlen= 10.0
print, "xlen= ", xlen
print, header
nx= sxpar(header,'NAXIS1')
ny= sxpar(header,'NAXIS2')
HalfMassContourLevel= sxpar(header,'HALFMASS')
print, "HalfMassContourLevel= ", HalfMassContourLevel
n= nx * ny
xbin= fltarr(n)
ybin= fltarr(n)
surfacedensity_list= fltarr(n)
dx= 2. * xlen / nx
dy= 2. * xlen / ny
idx= 0L
for i=0, nx-1 do begin
  for j=0, ny-1 do begin
        xbin[idx]= dx*(i+0.5) - xlen
        ybin[idx]= dy*(j+0.5) - xlen
        surfacedensity_list[idx]= surfacedensity[i,j]
        idx= idx + 1
  endfor
endfor





; scale it to 256 and plot
; -----------------------

MaxDensXY= max(surfacedensity)
DynRange=1.0e4

ma=MaxDensXY
mi=MaxDensXY/DynRange

print, "surfacedensity     max: ", max(surfacedensity), "   min: ", min(surfacedensity)
print, "Clipping at   ma= ", ma, " mi= ", mi

; do some clipping
surfacedensity_clipped= surfacedensity
ind=where(surfacedensity lt mi)
if ind(0) ne -1 then surfacedensity_clipped(ind)=mi
ind=where(surfacedensity gt ma)
if ind(0) ne -1 then surfacedensity_clipped(ind)=ma

; number of colors
Cols=255
; map it to 0->Cols
Pic=(alog10(surfacedensity_clipped)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2
HalfMCnt=(alog10(HalfMassContourLevel)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2


; these should never occur
print, "max(Pic)= ", max(Pic)
print, "min(Pic)= ", min(Pic)
ind=where(Pic ge 256)
if ind(0) ne -1 then Pic(ind)=255
ind=where(Pic le 0)
if ind(0) ne -1 then Pic(ind)=1


; I like my colors inverted
invertcolors= 1
;invertcolors= 0
if invertcolors eq 1 then begin
   Pic=256-Pic                          ; invert color table
   HalfMCnt= 256-HalfMCnt
   idx= where(Pic EQ 254)               ; set background to white
endif else begin
   idx= where(Pic EQ 2)
endelse

if idx(0) ne -1 then Pic(idx)= 1       ; this is setting it to white
;if idx(0) ne -1 then Pic(idx)= 0       ; this sets background to black

fload_newcolortable, 0

tv, Pic, x0, y2, xsize=(x1-x0), ysize=(y3-y2), /normal

print, "max(Pic)= ", max(Pic)
print, "min(Pic)= ", min(Pic)


; creates axes and plot style
; ---------------------------
xmax= +xlen
xmin= -xlen
ymax= +xlen
ymin= -xlen
!p.position=[x0,y2,x1,y3]
!p.ticklen=0.03
plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xthick=4.0, ythick=4.0, charthick=3.0, $
                      xtickformat='(a1)', ytickformat='(a1)', xstyle=1, ystyle=1, /nodata


; overplot contours, but not ellipse
oplot_contours, Pic, x0, y2, x1, y3, xlen=xlen, hmsb= HalfMCnt
SBPic= Pic




show_xlen= 1
;show_xlen= 0
if keyword_set(show_xlen) then begin
        xlenlbl= strcompress(string(2.0*xlen),/remove_all)
        if (2.0*xlen) ge 1.0 then digs= 1
        if (2.0*xlen) ge 10.0 then digs= 2
        if (2.0*xlen) ge 100.0 then digs= 3
        xlenlbl = strmid(xlenlbl,0,digs)        ; T=0.x (4+digits after decimal)
        xlenlbl= '!94!6'+xlenlbl+'!96!6'
        xyouts, x1-0.10, y3-0.05, xlenlbl, /normal, size= 1.2, charthick=3.0, color= 0
endif






;-----------------------------------------------------





; read in and display velocity map
; -------------------------------------
velocitymap=readfits(velocity_fits_file,header)
xlen= sxpar(header,'XLEN')  & if xlen le 0 then xlen= 10.0
print, "xlen= ", xlen
print, header
nx= sxpar(header,'NAXIS1')
ny= sxpar(header,'NAXIS2')
n= nx * ny
xbin= fltarr(n)
ybin= fltarr(n)
xgrid= velocitymap
ygrid= velocitymap
velocity_list= fltarr(n)
dx= 2. * xlen / nx
dy= 2. * xlen / ny
idx= 0L
for i=0, nx-1 do begin
  for j=0, ny-1 do begin
        xbin[idx]= dx*(i+0.5) - xlen
	xgrid[i,j]= dx*(i+0.5) - xlen
        ybin[idx]= dy*(j+0.5) - xlen
	ygrid[i,j]= dy*(j+0.5) - xlen
        velocity_list[idx]= velocitymap[i,j]
        idx= idx + 1
  endfor
endfor





; scale it to 256 and plot
; -----------------------


;ma= 250.0
ma= 150.0
mi= -ma

print, "velocitymap     max: ", max(velocitymap), "   min: ", min(velocitymap)
print, "Clipping at   ma= ", ma, " mi= ", mi

; do some clipping
velocitymap_clipped= velocitymap
ind=where(velocitymap lt mi)
if ind(0) ne -1 then velocitymap_clipped(ind)=mi
ind=where(velocitymap gt ma)
if ind(0) ne -1 then velocitymap_clipped(ind)=ma

; number of colors
Cols=255
; map it to 0->Cols
Pic= (velocitymap_clipped-mi)/(ma-mi) * (cols-3) + 2


; these should never occur
print, "max(Pic)= ", max(Pic)
print, "min(Pic)= ", min(Pic)
ind=where(Pic ge 256)
if ind(0) ne -1 then Pic(ind)=255
ind=where(Pic le 0)
if ind(0) ne -1 then Pic(ind)=1


; I like my colors inverted
invertcolors= 1
;invertcolors= 0
if invertcolors eq 1 then begin
   Pic=256-Pic                          ; invert color table
   idx= where(Pic EQ 254)               ; set background to white
endif else begin
   idx= where(Pic EQ 2)
endelse

if idx(0) ne -1 then Pic(idx)= 1       ; this is setting it to white
;if idx(0) ne -1 then Pic(idx)= 0       ; this sets background to black

;fload_newcolortable, 4
;fload_newcolortable, 12
fload_newcolortable, 25

tv, Pic, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal

print, "max(Pic)= ", max(Pic)
print, "min(Pic)= ", min(Pic)



; creates axes and plot style
; ---------------------------
xmax= +xlen
xmin= -xlen
ymax= +xlen
ymin= -xlen
!p.position=[x0,y1,x1,y2]
!p.ticklen=0.03
plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xthick=4.0, ythick=4.0, charthick=3.0, $
                      xtickformat='(a1)', ytickformat='(a1)', xstyle=1, ystyle=1, /nodata


; overplot contours, but not ellipse
oplot_contours, SBPic, x0, y1, x1, y2, xlen=xlen, hmsb=HalfMCnt






;
; read in dispersion map
; -------------------------------------
;dispersionmap= velocitymap
;dispersionmap[*,*]= 155.0
dispersionmap=readfits(dispersion_fits_file,header)
print, "dispersionmap     max: ", max(dispersionmap), "   min: ", min(dispersionmap)
xlen= sxpar(header,'XLEN')  &  if xlen le 0 then xlen= 10.0
print, "xlen= ", xlen
print, header







;  now calculate lambda
; ------------------------

; smooth map slightly
tempsd= surfacedensity
;width= 2
;tempMap= smooth(Pic, width, /edge_truncate)

tempv= abs(velocitymap)

idx= where(surfacedensity eq max(surfacedensity))
xg= xgrid - xgrid(idx(0))
yg= ygrid - ygrid(idx(0))
radius= sqrt(xg * xg + yg * yg)

; for the moment, we don't have the dispersion map
totalvel= sqrt(velocitymap*velocitymap + dispersionmap*dispersionmap)

levels= 100

prof_radius= fltarr(levels) & prof_radius(*)= -1
prof_lambda= fltarr(levels) & prof_lambda(*)= -1
prof_sigma= fltarr(levels) & prof_sigma(*)= -1
prof_meanabsvel= fltarr(levels) & prof_meanabsvel(*)= -1


for i= 0, levels-1 do begin

	;print, "i = ", i

	if i eq 0 then high_sd= max(tempsd) else high_sd= low_sd
	low_sd= high_sd * 0.9
	idx= where((tempsd gt low_sd) and (tempsd le high_sd))

	; if we're below the map, then get out of here
	if high_sd lt min(tempsd) then break

	if idx(0) ne -1 then begin
	     num= total(tempsd(idx) * radius(idx) * tempv(idx))
	     denom= total(tempsd(idx) * radius(idx) * totalvel(idx))
	     prof_lambda[i]= num / denom
	     prof_radius[i]= mean(radius(idx))
	     prof_meanabsvel[i]= mean(abs(velocitymap(idx)))
	     prof_sigma[i]= mean(dispersionmap(idx))
	endif

	if (HalfMassContourLevel gt low_sd) and (HalfMassContourLevel lt high_sd) then begin
		R_e= prof_radius[i]
	endif
endfor



idx= where(prof_radius lt 0.0)
if idx(0) ne -1 then begin
	ind= where(prof_radius gt 0.0)
	prof_radius= prof_radius(ind)
	prof_lambda= prof_lambda(ind)
	prof_sigma= prof_sigma(ind)
	prof_meanabsvel= prof_meanabsvel(ind)
endif



xaxistitle= "!6Radius (kpc)"
xmin= 0.0
xmax = max(prof_radius) * 1.1



;
;   lambda 
;
;-------------------------------------
yaxistitle= "!7k!6!DR!E" & ymax = 1.1  & ymin= 0.0
!p.position= [x0, y0, x1, y1]
!p.ticklen=0.03
plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
                xstyle=1, ystyle=1, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, $
                charthick=3.0, xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase

oplot, prof_radius, prof_lambda, psym= -2, color= 0, thick= 1.0, symsize= 2.0


x= [R_e, R_e]
y= [ymin, ymax]
oplot, x, y, psym=-3, color= 0, thick= 1.0, linestyle= 1


; ----------------------------------------------------
;
;   Old Method - grab values and return
;
;
;
;idx= where(prof_radius gt R_e)
;lam_Re= prof_lambda(idx(0))
;
;idx= where(prof_radius gt 5*R_e)
;lam_5Re= prof_lambda(idx(0))
;
;
; ----------------------------------------------------

;fname= fitsdir+'/lambda_'+istring+'.txt'
fname= fitsdir+'/'+istring+'_lst_lambda.txt'

openw, 1, fname, ERROR=err

printf, 1, "#    lambda_"+istring+".txt"
printf, 1, "# "
printf, 1, "# "
printf, 1, "# radius            Sigma    <|V|>   "
printf, 1, "#  (kpc)  Lambda    (km/s)   (km/s) "
for i=0, n_elements(plot_radius)-1 do begin
        printf, 1, FORMAT='(4(F8.4,"  "))', $
		prof_radius[i], $
		prof_lambda[i], $
		prof_sigma[i], $
		prof_meanabsvel
endfor

close, 1





; ----------------------------------------------------



device, /close

; -------------
;  Done
; -------------


end



; -----------------------------------------------------------------------------





;
;    4 panel version
;

pro do_lambda_profile_4, frun, kinnum, snapnum=snapnum, loadedsnap=loadedsnap, $
					filename=filename, $
					fitsdir=fitsdir, $
					sdfitsf=sdfitsf, vfitsf=vfitsf, dfitsf=dfitsf, $
					xlen=xlen, $
					lam_Re=lam_Re, $
					lam_5Re=lam_5Re


if not keyword_set(frun) then begin
        print, " "
        print, " do_lambda_profile_4, frun,"
        print, " "
        ;return
endif


if not keyword_set(snapnum) then begin
        print, " "
        print, " do_lambda_profile_4, frun,"
        print, " "
        ;return
endif






; --------------------------
;  Setup stuff
; --------------------------

if not keyword_set(filename) then filename='profile_kinematics.eps'
if not keyword_set(kinnum) then kinnum= 0



; panel dimensions
; ------------------


x0= 0.08
x1= 0.49
x2= 0.90

y0= 0.08
y1= 0.53
y2= 0.98
    



; ------------------
;  Set up file info
; ------------------

if not keyword_set(loadedsnap) then ok=fload_snapshot_bh(frun,snapnum, /nopot_in_snap)


; setup directory information
if not keyword_set(fitsdir) then fitsdir= frun+'/snap_'+fload_finfo_exts(1)+'_info/'

istring= strcompress(string(kinnum), /remove_all)

;surfacedensity_fits_file= fitsdir+'/kinmap_'+istring+'.fits'
;velocity_fits_file= fitsdir+'/kinmap_'+istring+'_vel.fits'
;dispersion_fits_file= fitsdir+'/kinmap_'+istring+'_dis.fits'
;if keyword_set(sdfitsf) then surfacedensity_fits_file= sdfitsf else surfacedensity_fits_file= fitsdir+'/kinmap_'+istring+'.fits'
;if keyword_set(vfitsf) then velocity_fits_file= vfitsf else velocity_fits_file= fitsdir+'/kinmap_'+istring+'_vel.fits'
;if keyword_set(dfitsf) then dispersion_fits_file= dfitsf else dispersion_fits_file= fitsdir+'/kinmap_'+istring+'_dis.fits'
if keyword_set(sdfitsf) then surfacedensity_fits_file= sdfitsf else surfacedensity_fits_file= fitsdir+'/'+istring+'_slitmap_sd.fits'
if keyword_set(vfitsf) then velocity_fits_file= vfitsf else velocity_fits_file= fitsdir+'/'+istring+'_slitmap_vel.fits'
if keyword_set(dfitsf) then dispersion_fits_file= dfitsf else dispersion_fits_file= fitsdir+'/'+istring+'_slitmap_dis.fits'


; testing
;fitsdir= "./"
;surfacedensity_fits_file= 'temp_kin.fits'
;velocity_fits_file= 'temp_kin_vel.fits'
;dispersion_fits_file= 'temp_kin_dis.fits'



; --------------
;   Plot it Up
; --------------


; get print stuff ready
; ----------------------
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 0, newxsize=32, newysize=30
print, "writing: ", filename




;-----------------------------------------------------



; read in and display surface density
; -------------------------------------
surfacedensity=readfits(surfacedensity_fits_file,header)
xlen= sxpar(header,'XLEN')  &  if xlen le 0 then xlen= 10.0
print, "xlen= ", xlen
print, header
nx= sxpar(header,'NAXIS1')
ny= sxpar(header,'NAXIS2')
HalfMassContourLevel= sxpar(header,'HALFMASS')
print, "HalfMassContourLevel= ", HalfMassContourLevel
n= nx * ny
xbin= fltarr(n)
ybin= fltarr(n)
surfacedensity_list= fltarr(n)
dx= 2. * xlen / nx
dy= 2. * xlen / ny
idx= 0L
for i=0, nx-1 do begin
  for j=0, ny-1 do begin
        xbin[idx]= dx*(i+0.5) - xlen
        ybin[idx]= dy*(j+0.5) - xlen
        surfacedensity_list[idx]= surfacedensity[i,j]
        idx= idx + 1
  endfor
endfor





; scale it to 256 and plot
; -----------------------

MaxDensXY= max(surfacedensity)
DynRange=1.0e4

ma=MaxDensXY
mi=MaxDensXY/DynRange

print, "surfacedensity     max: ", max(surfacedensity), "   min: ", min(surfacedensity)
print, "Clipping at   ma= ", ma, " mi= ", mi

; do some clipping
surfacedensity_clipped= surfacedensity
ind=where(surfacedensity lt mi)
if ind(0) ne -1 then surfacedensity_clipped(ind)=mi
ind=where(surfacedensity gt ma)
if ind(0) ne -1 then surfacedensity_clipped(ind)=ma

; number of colors
Cols=255
; map it to 0->Cols
Pic=(alog10(surfacedensity_clipped)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2
HalfMCnt=(alog10(HalfMassContourLevel)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2


; these should never occur
print, "max(Pic)= ", max(Pic)
print, "min(Pic)= ", min(Pic)
ind=where(Pic ge 256)
if ind(0) ne -1 then Pic(ind)=255
ind=where(Pic le 0)
if ind(0) ne -1 then Pic(ind)=1


; I like my colors inverted
invertcolors= 1
;invertcolors= 0
if invertcolors eq 1 then begin
   Pic=256-Pic                          ; invert color table
   HalfMCnt= 256-HalfMCnt
   idx= where(Pic EQ 254)               ; set background to white
endif else begin
   idx= where(Pic EQ 2)
endelse

if idx(0) ne -1 then Pic(idx)= 1       ; this is setting it to white
;if idx(0) ne -1 then Pic(idx)= 0       ; this sets background to black

fload_newcolortable, 0

tv, Pic, x0, y1, xsize=(x1-x0), ysize=(y2-y1), /normal

print, "max(Pic)= ", max(Pic)
print, "min(Pic)= ", min(Pic)

;
; creates axes and plot style
; ---------------------------
xmax= +xlen
xmin= -xlen
ymax= +xlen
ymin= -xlen
!p.position=[x0,y1,x1,y2]
!p.ticklen=0.03
plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xthick=4.0, ythick=4.0, charthick=3.0, $
                      xtickformat='(a1)', ytickformat='(a1)', xstyle=1, ystyle=1, /nodata


;
; overplot contours, but not ellipse
oplot_contours, Pic, x0, y1, x1, y2, xlen=xlen, hmsb=HalfMCnt
SBPic= Pic


;
;
show_xlen= 1
;show_xlen= 0
if keyword_set(show_xlen) then begin
        xlenlbl= strcompress(string(2.0*xlen),/remove_all)
        if (2.0*xlen) ge 1.0 then digs= 1
        if (2.0*xlen) ge 10.0 then digs= 2
        if (2.0*xlen) ge 100.0 then digs= 3
        xlenlbl = strmid(xlenlbl,0,digs)        ; T=0.x (4+digits after decimal)
        xlenlbl= '!94!6'+xlenlbl+' kpc/h !96!6'
        xyouts, x1-0.12, y2-0.03, xlenlbl, /normal, size= 1.4, charthick=3.0, color= 0
endif


;
;  put scale on right
bar= REPLICATE(1B, 10)#BINDGEN(256)
if invertcolors eq 1 then bar= 255-bar
idx= where(bar eq 0 or bar eq 1)
if idx(0) ne -1 then bar(idx)= 2   ; many of my color tables have 0 and 1 as black and white
barwidth= 0.025

tv, bar, x0-0.01-barwidth, y1+0.04, xsize=barwidth, ysize=(y2-y1-0.08), /normal

ma= alog10(ma * 0.7 * 1.0e4)
mi= alog10(mi * 0.7 * 1.0e4)

!p.position= [x0-0.01-barwidth,y1+0.04,x0-0.01,y2-0.04]
!p.ticklen= 0.20
plot, [0],[0], psym=3, xrange=[0,1], yrange=[mi,ma], /normal, /noerase, color= 0, /nodata, $
                        xstyle=1, ystyle=1, xthick=4.0, ythick=4.0, $
                        xticks=1, xtickformat='(a1)', ytickformat='(a1)'

ytit= 'Log !7R!6 (M!D!9n!6!N pc!E-2!N)' 
axis, yaxis=0, yrange=[mi,ma], ystyle=1, /normal, ycharsize=1.1, charthick=3.0, ytitle= ytit





;-----------------------------------------------------





; read in and display velocity map
; -------------------------------------
velocitymap=readfits(velocity_fits_file,header)
xlen= sxpar(header,'XLEN')  & if xlen le 0 then xlen= 10.0
print, "xlen= ", xlen
print, header
nx= sxpar(header,'NAXIS1')
ny= sxpar(header,'NAXIS2')
n= nx * ny
xbin= fltarr(n)
ybin= fltarr(n)
xgrid= velocitymap
ygrid= velocitymap
velocity_list= fltarr(n)
dx= 2. * xlen / nx
dy= 2. * xlen / ny
idx= 0L
for i=0, nx-1 do begin
  for j=0, ny-1 do begin
        xbin[idx]= dx*(i+0.5) - xlen
	xgrid[i,j]= dx*(i+0.5) - xlen
        ybin[idx]= dy*(j+0.5) - xlen
	ygrid[i,j]= dy*(j+0.5) - xlen
        velocity_list[idx]= velocitymap[i,j]
        idx= idx + 1
  endfor
endfor





; scale it to 256 and plot
; -----------------------


;ma= 250.0
ma= 150.0
mi= -ma

print, "velocitymap     max: ", max(velocitymap), "   min: ", min(velocitymap)
print, "Clipping at   ma= ", ma, " mi= ", mi

; do some clipping
velocitymap_clipped= velocitymap
ind=where(velocitymap lt mi)
if ind(0) ne -1 then velocitymap_clipped(ind)=mi
ind=where(velocitymap gt ma)
if ind(0) ne -1 then velocitymap_clipped(ind)=ma

; number of colors
;Cols=255
Cols=254
; map it to 0->Cols
;Pic= (velocitymap_clipped-mi)/(ma-mi) * (cols-3) + 2
Pic= (velocitymap_clipped-mi)/(ma-mi) * Cols + 2


; these should never occur
print, "max(Pic)= ", max(Pic)
print, "min(Pic)= ", min(Pic)
ind=where(Pic ge 254)
if ind(0) ne -1 then Pic(ind)=254
ind=where(Pic le 1)
if ind(0) ne -1 then Pic(ind)=2


; I like my colors inverted
invertcolors= 1
;invertcolors= 0
if invertcolors eq 1 then begin
   Pic=256-Pic                          ; invert color table
   idx= where(Pic gt 255)               ; set background to white
   if idx(0) ne -1 then Pic(idx)= 255
endif else begin
   ;idx= where(Pic EQ 2)
endelse

;if idx(0) ne -1 then Pic(idx)= 1       ; this is setting it to white
;if idx(0) ne -1 then Pic(idx)= 0       ; this sets background to black

;fload_newcolortable, 4
;fload_newcolortable, 12
fload_newcolortable, 25

tv, Pic, x1, y1, xsize=(x2-x1), ysize=(y2-y1), /normal

print, "max(Pic)= ", max(Pic)
print, "min(Pic)= ", min(Pic)



; creates axes and plot style
; ---------------------------
xmax= +xlen
xmin= -xlen
ymax= +xlen
ymin= -xlen
!p.position=[x1,y1,x2,y2]
!p.ticklen=0.03
plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xthick=4.0, ythick=4.0, charthick=3.0, $
                      xtickformat='(a1)', ytickformat='(a1)', xstyle=1, ystyle=1, /nodata


;
; overplot contours, but not ellipse
oplot_contours, SBPic, x1, y1, x2, y2, xlen=xlen, hmsb=HalfMCnt


;
;  put scale on right
bar= REPLICATE(1B, 10)#BINDGEN(256)
if invertcolors eq 1 then bar= 255-bar 
idx= where(bar eq 0 or bar eq 1)
if idx(0) ne -1 then bar(idx)= 2   ; many of my color tables have 0 and 1 as black and white
barwidth= 0.025

tv, bar, x2+0.01, y1+0.04, xsize=barwidth, ysize=(y2-y1-0.08), /normal

!p.position= [x2+0.01,y1+0.04,x2+0.01+barwidth,y2-0.04]
!p.ticklen= 0.20
plot, [0],[0], psym=3, xrange=[0,1], yrange=[mi,ma], /normal, /noerase, color= 0, /nodata, $
                        xstyle=1, ystyle=1, xthick=4.0, ythick=4.0, $
                        xticks=1, xtickformat='(a1)', ytickformat='(a1)'

ytit= '!6Velocity (km s!E-1!N)'        
axis, yaxis=1, yrange=[mi,ma], ystyle=1, /normal, ycharsize=1.1, charthick=3.0, ytitle= ytit






;
;
;
;
;
; read in dispersion map
; -------------------------------------

;ma= 250.0
ma= 200.0
mi= 0.0

dispersionmap=readfits(dispersion_fits_file,header)
print, "dispersionmap     max: ", max(dispersionmap), "   min: ", min(dispersionmap)
xlen= sxpar(header,'XLEN')  &  if xlen le 0 then xlen= 10.0
print, "xlen= ", xlen
print, header

print, "Clipping at   ma= ", ma, " mi= ", mi

; do some clipping
dispersionmap_clipped= dispersionmap
ind=where(dispersionmap lt mi)
if ind(0) ne -1 then dispersionmap_clipped(ind)=mi
ind=where(dispersionmap gt ma)
if ind(0) ne -1 then dispersionmap_clipped(ind)=ma

; number of colors
Cols=255
; map it to 0->Cols
Pic= (dispersionmap_clipped-mi)/(ma-mi) * (cols-3) + 2


; these should never occur
print, "max(Pic)= ", max(Pic)
print, "min(Pic)= ", min(Pic)
ind=where(Pic ge 256)
if ind(0) ne -1 then Pic(ind)=255
ind=where(Pic le 0)
if ind(0) ne -1 then Pic(ind)=1


; I like my colors inverted
invertcolors= 1
;invertcolors= 0
if invertcolors eq 1 then begin 
   Pic=256-Pic                          ; invert color table
   idx= where(Pic EQ 254)               ; set background to white
endif else begin 
   idx= where(Pic EQ 2)
endelse

if idx(0) ne -1 then Pic(idx)= 1       ; this is setting it to white
;if idx(0) ne -1 then Pic(idx)= 0       ; this sets background to black

;fload_newcolortable, 4
fload_newcolortable, 5
;fload_newcolortable, 7
;fload_newcolortable, 12
;fload_newcolortable, 25
;fload_newcolortable, 33

tv, Pic, x1, y0, xsize=(x2-x1), ysize=(y1-y0), /normal

print, "max(Pic)= ", max(Pic)
print, "min(Pic)= ", min(Pic)



; creates axes and plot style
; ---------------------------
xmax= +xlen
xmin= -xlen
ymax= +xlen
ymin= -xlen
!p.position=[x1,y0,x2,y1]
!p.ticklen=0.03
plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, color= 0, $
                      xcharsize=0.01, ycharsize=0.01, xthick=4.0, ythick=4.0, charthick=3.0, $
                      xtickformat='(a1)', ytickformat='(a1)', xstyle=1, ystyle=1, /nodata


;
; overplot contours, but not ellipse
oplot_contours, SBPic, x1, y0, x2, y1, xlen=xlen, hmsb=HalfMCnt


;
;  put scale on right
bar= REPLICATE(1B, 10)#BINDGEN(256)
if invertcolors eq 1 then bar= 255-bar
idx= where(bar eq 0 or bar eq 1)
if idx(0) ne -1 then bar(idx)= 2   ; many of my color tables have 0 and 1 as black and white
barwidth= 0.025

tv, bar, x2+0.01, y0+0.04, xsize=barwidth, ysize=(y1-y0-0.08), /normal

!p.position= [x2+0.01,y0+0.04,x2+0.01+barwidth,y1-0.04]
!p.ticklen= 0.20
plot, [0],[0], psym=3, xrange=[0,1], yrange=[mi,ma], /normal, /noerase, color= 0, /nodata, $
                        xstyle=1, ystyle=1, xthick=4.0, ythick=4.0, $
                        xticks=1, xtickformat='(a1)', ytickformat='(a1)'

ytit= '!7r!6 (km s!E-1!N)'
axis, yaxis=1, yrange=[mi,ma], ystyle=1, /normal, ycharsize=1.1, charthick=3.0, ytitle= ytit








;
;
;
;
;
;  now calculate lambda
; ------------------------

; smooth map slightly
tempsd= surfacedensity
;width= 2
;tempMap= smooth(Pic, width, /edge_truncate)

tempv= abs(velocitymap)

idx= where(surfacedensity eq max(surfacedensity))
xg= xgrid - xgrid(idx(0))
yg= ygrid - ygrid(idx(0))
radius= sqrt(xg * xg + yg * yg)

; for the moment, we don't have the dispersion map
totalvel= sqrt(velocitymap*velocitymap + dispersionmap*dispersionmap)

levels= 100

prof_radius= fltarr(levels) & prof_radius(*)= -1
prof_lambda= fltarr(levels) & prof_lambda(*)= -1
prof_sigma= fltarr(levels) & prof_sigma(*)= -1
prof_meanabsvel= fltarr(levels) & prof_meanabsvel(*)= -1


for i= 0, levels-1 do begin

	;print, "i = ", i

	if i eq 0 then high_sd= max(tempsd) else high_sd= low_sd
	low_sd= high_sd * 0.9
	idx= where((tempsd gt low_sd) and (tempsd le high_sd))

	; if we're below the map, then get out of here
	if high_sd lt min(tempsd) then break

	if idx(0) ne -1 then begin
	     num= total(tempsd(idx) * radius(idx) * tempv(idx))
	     denom= total(tempsd(idx) * radius(idx) * totalvel(idx))
	     prof_lambda[i]= num / denom
	     prof_radius[i]= mean(radius(idx))
	     prof_meanabsvel[i]= mean(abs(velocitymap(idx)))
	     prof_sigma[i]= mean(dispersionmap(idx))
	endif

	if (HalfMassContourLevel gt low_sd) and (HalfMassContourLevel lt high_sd) then begin
		R_e= prof_radius[i]
	endif
endfor



idx= where(prof_radius lt 0.0)
if idx(0) ne -1 then begin
	ind= where(prof_radius gt 0.0)
	prof_radius= prof_radius(ind)
	prof_lambda= prof_lambda(ind)
	prof_sigma= prof_sigma(ind)
	prof_meanabsvel= prof_meanabsvel(ind)
endif



xaxistitle= "!6Radius (kpc/h)"
xmin= 0.0
xmax = max(prof_radius) * 1.1



;
;   lambda 
;
;-------------------------------------
yaxistitle= "!7k!6!DR!E" & ymax = 1.1  & ymin= 0.0
!p.position= [x0, y0, x1, y1]
!p.ticklen=0.03
plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
                xstyle=1, ystyle=1, xcharsize=1.5, ycharsize=1.5, xthick=4.0, ythick=4.0, $
                charthick=3.0, xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase

oplot, prof_radius, prof_lambda, psym= -2, color= 0, thick= 1.0, symsize= 2.0


x= [R_e, R_e]
y= [ymin, ymax]
oplot, x, y, psym=-3, color= 0, thick= 1.0, linestyle= 1


; ----------------------------------------------------
;
;   Old Method - grab values and return
;
;
;
;idx= where(prof_radius gt R_e)
;lam_Re= prof_lambda(idx(0))
;
;idx= where(prof_radius gt 5*R_e)
;lam_5Re= prof_lambda(idx(0))
;
;
; ----------------------------------------------------

;fname= fitsdir+'/lambda_'+istring+'.txt'
fname= fitsdir+'/'+istring+'_lst_lambda.txt'

openw, 1, fname, ERROR=err

printf, 1, "#    lambda_"+istring+".txt"
printf, 1, "# "
printf, 1, "# "
printf, 1, "# radius              Sigma     <|V|>   "
printf, 1, "#  (kpc)    Lambda    (km/s)    (km/s) "
for i=0, n_elements(prof_radius)-1 do begin
        printf, 1, FORMAT='(4(F8.4,"  "))', $
		prof_radius[i], $
		prof_lambda[i], $
		prof_sigma[i], $
		prof_meanabsvel[i]
endfor

close, 1





; ----------------------------------------------------



device, /close

; -------------
;  Done
; -------------


end



; -----------------------------------------------------------------------------
