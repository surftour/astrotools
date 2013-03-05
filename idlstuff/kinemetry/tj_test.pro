;################################################################
; This is an example routine which calls KINEMETRY 
; to analyse SAURON velocity map of NGC2974 as presented in
; (Emsellem et al. 2004 MNRAS, 352, 271). It makes a plot
; of kinemetric coefficients as in the first column of Fig.7 
; in Krajnovic et al. (2006).
;
; This routine uses RDFLOAT.PRO and PLOTERROR.PRO avaiable from
; IDL Astronomer Unser's Library.
; 
; Davor Krajnovic, Oxford, 07.12.2005.
;###############################################################

PRO kinemetry_example

;
; read in all data
;
;file = 'NGC2974_SAURON_kinematics.dat'
;rdfloat, file, num, xbin, ybin, velbin, er_velbin,  SKIPLIN=1

;test=readfits("/n/home/tcox/temp_kin_vel.fits",header)
;test=readfits("/Users/tcox/temp_kin_vel.fits",header)
;test=readfits("/Users/tcox/test_48_vel.fits",header)
test=readfits("/n/home/tcox/test_48_vel.fits",header)
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
nbins= 30
radius= xlen * findgen(nbins)/float(nbins-1)
radius= radius[1:nbins-1]

TJ_KINEMETRY, xbin, ybin, velbin, rad, pa, q, cf, ntrm=6, scale=0.8, $
	ERROR=er_velbin, name='TJs Test',er_cf=er_cf, er_pa=er_pa, $
	er_q=er_q, /plot, /verbose, RADIUS=radius   ;, npa= 5, nq= 5
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
;r = GET_SCREEN_SIZE()
;window, 1, xsize=r[0]*0.3, ysize=r[1]*0.8
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename="tj_test_again.eps", colortable= 0, newxsize= 18.0, newysize= 30.0
!y.style=1
!p.multi=[0,1,4]
!Y.MARGIN=[0,0] ; show plots with shared X axis
!Y.OMARGIN=[6,5] ; allow for space for the axis labels
print, "pa= ", min(pa), max(pa)
;ploterror, rad, pa, er_pa, PSYM=-5, TITLE='TJs Test', xtickformat = '(A1)', YTITLE='!7C!X!N [degrees]', YRANGE=[30,70], color= 0, charsize= 2.5
ploterror, rad, pa, er_pa, PSYM=-5, TITLE='TJs Test', xtickformat = '(A1)', YTITLE='!7C!X!N [degrees]', YRANGE=[min(pa)-10,max(pa)+10], color= 0, charsize= 2.5
;ploterror, rad, q, er_q, PSYM=-5, YRANGE=[0,1.1], xtickformat = '(A1)', YTITLE='q', color= 0, charsize= 2.5
ellip= 1.0 - q
ploterror, rad, ellip, er_q, PSYM=-5, YRANGE=[0,1.1], xtickformat = '(A1)', YTITLE='ellipticity', color= 0, charsize= 2.5
ploterror, rad, k1, erk1, PSYM=-5, xtickformat = '(A1)', YTITLE='k1 [km/s]',YRANGE=[0,245], color= 0, charsize= 2.5
ploterror, rad, k51, erk51, PSYM=-5, XTITLE='R [kpc/h]', YTITLE='k5/k1', YRANGE=[0,0.13], color= 0, charsize= 2.5
!P.MULTI=0
!Y.MARGIN=[4,2] ; back to default values
!Y.OMARGIN=[0,0]
device,/close

END
