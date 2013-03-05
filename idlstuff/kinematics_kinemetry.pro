;
;
;  This goes through a directory and generates the kinemetry
;  from the kinematics info (which was run while skipping the kinemetry).
;
;
pro process_kinemetry, frun

if not keyword_set(frun) then begin
        print, "  "
        print, " process_kinemetry, frun"
        print, "  "
        return
endif


small_fov= 1
if keyword_set(small_fov) then begin

	spawn, "head -n 1 "+frun+"/kinematics_small_fov.txt ", result
	tempjunk= strsplit(result[0],/extract,count=count)
	tempjunk='000'+tempjunk(4)
	snapwithkin= strmid(tempjunk,strlen(tempjunk)-3,3)

	fitsdir= frun+'/snap_'+snapwithkin+'_small_fov/'
	print, fitsdir

	spawn, "wc "+frun+"/kinematics_small_fov.txt ", result
	Nangles= long(result[0])-5
	print, Nangles
endif


;
; now, actually produce the kinemetry information
;
for i= 0, Nangles-1 do begin
	
	istring= strcompress(string(i), /remove_all)
	;print, istring
	kinematics_kinemetry, frun, nbins= 30, istring=istring, fitsdir=fitsdir
endfor


end




;################################################################
;
;
;
; This is based after the kinemetry_example used in the 
; kinemetry subdirectory.  
;
; 
;###############################################################

pro kinematics_kinemetry, frun, nbins= nbins, $
				fitsdir=fitsdir, $
				istring= istring


spawn, "mkdir "+fitsdir+"/"+istring

;
vfitsf= fitsdir+'/'+istring+'_slitmap_48_vel.fits'

;test=readfits("/n/home/tcox/temp_kin_vel.fits",header)
;test=readfits("/Users/tcox/temp_kin_vel.fits",header)
;test=readfits("/Users/tcox/test_48_vel.fits",header)
;test=readfits("/n/home/tcox/test_48_vel.fits",header)
test=readfits(vfitsf,header)
;stop
xlen= sxpar(header,'XLEN')   &   if xlen le 0 then xlen= 10.0
nx= sxpar(header,'NAXIS1')
ny= sxpar(header,'NAXIS2')
n= nx * ny
xbin= fltarr(n)
ybin= fltarr(n)
velbin= fltarr(n)
dx= 2. * xlen / nx
dy= 2. * xlen / ny
idx= 0L
for i=0, nx-1 do begin
  for j=0, ny-1 do begin
	xbin[idx]= dx*(i+0.5) - xlen
	ybin[idx]= dy*(j+0.5) - xlen
	velbin[idx]= test[i,j]
	idx= idx + 1
  endfor
endfor

;
; kinemetry on velocity map
;
t=systime(1)

; manually set the radial bins
if not keyword_set(nbins) then nbins= 30
radius= xlen * findgen(nbins)/float(nbins-1)
radius= radius[1:nbins-1]

KINEMETRY, xbin, ybin, velbin, rad, pa, q, cf, ntrm=6, scale=1.0, $
	ERROR=er_velbin, name='Kinemetry',er_cf=er_cf, er_pa=er_pa, $
	er_q=er_q, /plot, /verbose, RADIUS=radius, $   ;, npa= 5, nq= 5
	fitsdir=fitsdir, istring=istring
print, systime(1) -t, 'seconds'

;
; kinemetry parameters as defined in Krajnovic et al. (2006)
; 
k0 = cf[*,0]
k1 = SQRT(cf[*,1]^2 + cf[*,2]^2)
k5 = SQRT(cf[*,5]^2 + cf[*,6]^2)
k51 = k5/k1
erk1 = (SQRT( (cf[*,1]*er_cf[*,1])^2 + (cf[*,2]*er_cf[*,2])^2 ))/k1
erk5 = (SQRT( (cf[*,5]*er_cf[*,5])^2 + (cf[*,6]*er_cf[*,6])^2 ))/k5
erk51 = ( SQRT( ((k5/k1) * erk1)^2 + erk5^2  ) )/k1 

;
; plot coeffs.
;

epsfilename= fitsdir+'/'+istring+'_kinemetry.eps'

;r = GET_SCREEN_SIZE()
;window, 1, xsize=r[0]*0.3, ysize=r[1]*0.8
initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename="tj_test_again.eps", colortable= 0, newxsize= 18.0, newysize= 30.0
setup_plot_stuff, 'ps', filename=epsfilename, colortable= 0, newxsize= 18.0, newysize= 30.0
;!y.style=1
;!p.multi=[0,1,4]
;!Y.MARGIN=[0,0] ; show plots with shared X axis
;!Y.OMARGIN=[6,5] ; allow for space for the axis labels
print, "pa= ", min(pa), max(pa)
!p.position= [0.14, 0.72, 0.96, 0.94]
!p.ticklen= 0.03
;ploterror, rad, pa, er_pa, PSYM=-5, TITLE='TJs Test', xtickformat = '(A1)', YTITLE='!7C!X!N [degrees]', YRANGE=[30,70], color= 0, charsize= 2.5
ploterror, rad, pa, er_pa, PSYM=-5, TITLE=frun+", "+istring, xtickformat = '(A1)', YTITLE='!7C!X!N [degrees]', YRANGE=[min(pa)-10,max(pa)+10], color= 0, charsize= 1.3
;ploterror, rad, q, er_q, PSYM=-5, YRANGE=[0,1.1], xtickformat = '(A1)', YTITLE='q', color= 0, charsize= 2.5
!p.position= [0.14, 0.50, 0.96, 0.72]
ellip= 1.0 - q
ploterror, rad, ellip, er_q, PSYM=-5, YRANGE=[0,1.1], xtickformat = '(A1)', YTITLE='ellipticity', color= 0, charsize= 1.3, /noerase
!p.position= [0.14, 0.28, 0.96, 0.50]
ploterror, rad, k1, erk1, PSYM=-5, xtickformat = '(A1)', YTITLE='k1 [km/s]',YRANGE=[0,245], color= 0, charsize= 1.3, /noerase
!p.position= [0.14, 0.06, 0.96, 0.28]
ploterror, rad, k51, erk51, PSYM=-5, XTITLE='R [kpc/h]', YTITLE='k5/k1', YRANGE=[0,0.25], color= 0, charsize= 1.3, /noerase
;!P.MULTI=0
;!Y.MARGIN=[4,2] ; back to default values
;!Y.OMARGIN=[0,0]
device,/close



;
;    write information to file
;
filename= fitsdir+'/'+istring+'_lst_kinemetry.txt'
openw, 1, filename, ERROR=err

printf, 1, "#   ###_lst_kinemetry.txt, "+frun
printf, 1, "#   "
printf, 1, "#                                                  "
printf, 1, "#  radius      pa                            k1            "
printf, 1, "# (kpc/h)     (deg) ellipticity     q      (km/s)     k5/k1" 


nbins= n_elements(rad)

for i=0,nbins-1 do begin
        printf, 1, FORMAT= '(" ",6(F8.3,"  "))', $
                rad[i], pa[i], ellip[i], q[i], k1[i], k51[i]
endfor
close, 1







END



;
;
;======================================================================
;
;









;
;
;======================================================================
;
;





pro plot_compiled_radial_data, fitsdir, istring, filename=filename




if not keyword_set(fitsdir) then begin
   print, "  "
   print, "plot_compiled_radial_data, fitsdir, istring, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename=fitsdir+'/x_radial.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize= 18.0, newysize= 30.0



;--------------------------------------
;--------------------------------------


x0= 0.13
x1= 0.97

y0= 0.08
yend= 0.98

npanels= 4

ys= (yend-y0)/npanels
y1= y0 + ys
y2= y0 + ys + ys
y3= y0 + ys + ys + ys
y4= y0 + ys + ys + ys + ys
y5= y0 + ys + ys + ys + ys + ys



;--------------------------------------
;--------------------------------------


xaxistitle= "R [kpc/h]"

;xmax= max(radius)
xmax= 20.0
xmin=  0.0



;--------------------------------------
;--------------------------------------



frun= fitsdir+"/"+istring+"_lst_kinemetry.txt"
read_kinemetry, frun, -1, k_radius, k_pa, k_ellip, k_q, k_k1, k_k1_err, k_k5k1, k_k5k1_err



frun= fitsdir+"/"+istring+"_lst_isophotes.txt"
read_isophotes, frun, -1, i_radius, i_sd, i_ellip, i_pa, i_a4diva



frun= fitsdir+"/"+istring+"_lst_lambda.txt"
read_lambda, frun, -1, l_radius, l_lambda, l_sigma, l_v



frun= fitsdir+"/"+istring+"_lst_slit_major.txt"
read_slit_info, frun, -1, radius_major, vrot_major, sigma_major



frun= fitsdir+"/"+istring+"_lst_slit_minor.txt"
read_slit_info, frun, -1, radius_minor, vrot_minor, sigma_minor



xmax= max([k_radius, i_radius, l_radius, radius_major, radius_minor])
print, "xmax= ", xmax


;--------------------------------------
;--------------------------------------


yaxistitle= "PA [degrees]"
ymax=  180.0
ymin= -180.0

!p.position=[x0, y3, x1,y4]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.4, ycharsize=1.4, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
        ytitle=yaxistitle, xtickformat='(a1)'


oplot, i_radius, i_pa, psym= -5, color= 50, symsize= 1.2, thick= 2.0
;ploterror, i_radius, i_pa, i_pa_err, PSYM=-5, color= 50, /noerase
xyouts, 0.75, y3+0.06, /normal, "isophotes", size=1.4, charthick= 3.0, color=50

oplot, k_radius, k_pa, psym= -2, color= 150, symsize= 1.2, thick= 2.0
;ploterror, k_radius, k_pa, k_pa_err, PSYM=-2, color= 150, /noerase
xyouts, 0.75, y3+0.03, /normal, "kinematic", size=1.4, charthick= 3.0, color=150





;--------------------------------------


yaxistitle= "Ellipticity"
ymax= 1.1
ymin= 0.0

!p.position=[x0, y2, x1, y3]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.4, ycharsize=1.4, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
        ytitle=yaxistitle, xtickformat='(a1)'


oplot, i_radius, i_ellip, psym= -5, color= 50, symsize= 1.2, thick= 2.0
;ploterror, i_radius, i_ellip, i_ellip_err, PSYM=-5, color= 50, /noerase
;xyouts, 0.7, y4-0.08, /normal, "isophotes", size=1.2, charthick= 3.0, color=50

oplot, k_radius, k_ellip, psym= -2, color= 150, symsize= 1.2, thick= 2.0
;ploterror, k_radius, k_ellip, k_ellip_err, PSYM=-2, color= 150, /noerase
;xyouts, 0.7, y4-0.15, /normal, "kinematic", size=1.2, charthick= 3.0, color=150



   

; --------------------


yaxistitle= "Velocity [km/s]"
ymax= 185.0
ymin= 0.0

!p.position=[x0, y1, x1, y2]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.4, ycharsize=1.4, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
        ytitle=yaxistitle, xtickformat='(a1)'



nr= long(n_elements(radius_major))
nrhalf= long(nr/2)
rm= fltarr(nrhalf)
vm= fltarr(nrhalf)
for i= 0, nrhalf-1 do begin
   rm[i]= 0.5 * (abs(radius_major[i]) + abs(radius_major[nr-1-i]))
   vm[i]= 0.5 * (abs(vrot_major[i]) + abs(vrot_major[nr-1-i]))
endfor
oplot, rm, vm, psym= -3, color= 0, thick= 4.0
;ploterror, radius_major, vrot_major, vrot_major_err, psym=-3, color= 0, /noerase, thick= 4.0
xyouts, 0.2, y2-0.03, /normal, "major", size=1.4, charthick= 5.0, color=0

nr= long(n_elements(radius_minor))
nrhalf= long(nr/2)
rm= fltarr(nrhalf)
vm= fltarr(nrhalf)
for i= 0, nrhalf-1 do begin
   rm[i]= 0.5 * (abs(radius_minor[i]) + abs(radius_minor[nr-1-i]))
   vm[i]= 0.5 * (abs(vrot_minor[i]) + abs(vrot_minor[nr-1-i]))
endfor
oplot, rm, vm, psym= -3, color= 0, thick= 1.0, linestyle= 1
;ploterror, radius_minor, vrot_minor, vrot_minor_err, psym=-3, color= 0, /noerase
xyouts, 0.2, y2-0.055, /normal, "minor", size=1.4, charthick= 1.0, color=0

oplot, l_radius, 200.0*l_lambda, psym= -7, color= 100, thick= 1.0
xyouts, 0.75, y1+0.06, /normal, "lambda (x200)", size=1.4, charthick= 1.0, color=100

oplot, k_radius, k_k1, psym= -2, color= 150, symsize= 1.2, thick= 2.0
;ploterror, k_radius, k_k1, k_k1_err, PSYM=-2, color= 150, /noerase
xyouts, 0.75, y1+0.03, /normal, "k1", size=1.4, charthick= 3.0, color=150




; --------------------

yaxistitle= " "
ymax=  0.44
ymin= -0.44

!p.position=[x0, y0, x1, y1]

plot, [0], [0], psym=3,xstyle=1,ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax], color=0, $
        xcharsize=1.4, ycharsize=1.4, xthick=4.0, ythick=4.0, charthick=3.0, /noerase, /nodata, $
        ytitle=yaxistitle, xtitle=xaxistitle

x= [xmin, xmax]
y= [0.0, 0.0]
oplot, x, y, psym=-3, color= 0, linestyle= 1

oplot, i_radius, 10.0*i_a4diva, psym= -5, color= 50, symsize= 1.2, thick= 2.0
;ploterror, k_radius, i_a4diva, i_a4diva_err, PSYM=-2, color= 50, /noerase
xyouts, 0.75, y0+0.06, /normal, '10 x a4/a', size=1.2, charthick= 3.0, color=50

oplot, k_radius, k_k5k1, psym= -2, color= 150, symsize= 1.2, thick= 2.0
;ploterror, k_radius, k_k5k1, k_k5k1_err, PSYM=-2, color= 150, /noerase
xyouts, 0.75, y0+0.03, /normal, "k5/k1", size=1.2, charthick= 3.0, color=150



; --------------------

device, /close






end



