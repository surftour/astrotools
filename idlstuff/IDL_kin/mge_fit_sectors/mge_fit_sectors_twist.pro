;######################################################################
;
; Copyright (C) 1999-2002, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; For details on the method see:
;   Cappellari M., 2002, MNRAS, 333, 400
;
; Updated versions of the software are available from my web page
; http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
; If you have found this software useful for your
; research, I would appreciate an acknowledgment to use of
; `the MGE fitting software developed by Cappellari (2002)'.
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;######################################################################
;+
; NAME:
;     MGE_FIT_SECTORS_TWIST
;
; AUTHOR:
;       Michele Cappellari, Astrophysics Sub-department, University of Oxford, UK
;
; PURPOSE:
;       Fit a Multi Gaussian Expansion (MGE) model to a set of galaxy surface
;       brightness measurements. The measurements are usually taken along
;       sectors with a previous call to the routine SECTORS_PHOTOMETRY_TWIST.
;       The MGE model is intended to be used as a parametrization for
;       the galaxy surface brightness. All measurements within this program
;       are in the instrumental units of PIXELS and COUNTS.
;       This routine fits MGE models with varying position angle and common center.
;
; EXPLANATION:
;       Further information on MGE_FIT_SECTORS_TWIST is available in
;       Cappellari M., 2002, MNRAS, 333, 400
;
; CALLING SEQUENCE:
;       MGE_FIT_SECTORS_TWIST, Radius, Angle, Counts, Eps,
;           NGAUSS=ngauss, SIGMAPSF=[sigma1,sigma2,...], NORMPSF=[norm1,norm2,...],
;           SCALE=scale, RBOUNDS=[rmin,rmax], QBOUNDS=[qmin,qmax],
;           PABOUNDS=[PAmin,PAmax], /PRINT, /FASTNORM, /NEGATIVE, SOL=sol,
;           OUTER_SLOPE=outer_slope
;
; INPUTS:
;       Radius = Vector containing the radius of the surface brightness
;               measurements, taken from the galaxy center. This is given
;               in units of PIXELS (!) of the high resolution image.
;       Angle = Vector containing the polar angle of the surface brightness
;               measurements, taken from the galaxy major axis.
;       Counts = Vector containing the actual surface brightness measurements
;               in COUNTS (!) at the polar coordinates specified by the vectors
;               Radius and Angle. These three vectors need to have the same
;               number of elements.
;       Eps = Estimate for the galaxy `average' ellipticity Eps = 1-b/a = 1-q'
;
; OUTPUTS:
;       No output parameters. The results are printed on the screen, plotted
;       in a PS file with the /PRINT keyword and passed with the optional
;       keyword SOL
;
; OPTIONAL INPUT KEYWORDS:
;       /NEGATIVE - Set this keyword to allow for negative Gaussians in the fit.
;               Use this only if everything else didn't work or if there is clear
;               evidence that negative Gaussians are actually needed.
;       NGAUSS - Number of Gaussians desired in the MGE fit.
;               Typical values are in the range 10-20 (default: 12).
;       NORMPSF - This is optional if only a scalar is given for SIGMAPSF,
;               otherwise it must contain the normalization of each MGE component
;               of the PSF, whose sigma is given by SIGMAPSF. The vector needs to
;               have the same number of elements of SIGMAPSF and the condition
;               TOTAL(normpsf) = 1 must be verified. In other words the MGE PSF
;               needs to be normalized. (default: 1).
;       OUTER_SLOPE - This scalar forces the surface brightness profile of
;               the MGE model to decrease at least as fast as R^(-OUTER_SLOPE)
;               at the largest measured radius (Default: OUTER_SLOPE=2).
;       PABOUNDS - Two elements vector giving the minimum and maximum
;               position angle PA allowed in the MGE fit.
;       /PRINT - Set this keyword to print the best-fitting MGE profiles in IDL.PS
;       QBOUNDS - Two elements vector giving the minimum and maximum
;               axial ratio Q allowed in the MGE fit.
;       RBOUNDS - Two elements vector giving the minimum and maximum sigma
;               allowed in the MGE fit. This is in PIXELS units.
;       SCALE - the pixel scale in arcsec/pixels. This is only needed for the plots.
;               (default: 1)
;       SIGMAPSF - Scalar giving the sigma of the PSF of the high resolution
;               image (see Cappellari [2002, pg. 406] for details), or vector
;               with the sigma of an MGE model for the circular PSF.
;               This has to be in pixels, as the vector RADIUS above.
;               (Default: no convolution)
;
; OPTIONAL OUTPUT KEYWORDS:
;       SOL - Output keyword containing a 4xNgauss array with the
;               best-fitting solution:
;               1) sol[0,*] = TotalCounts, of the Gaussians components.
;                   The relation TotalCounts = Height*(2*!PI*Sigma^2*qObs)
;                   can be used compute the Gaussian central surface
;                   brightness (Height)
;               2) sol[1,*] = Sigma, is the dispersion of the best-fitting
;                   Gaussians in pixels.
;               3) sol[2,*] = qObs, is the observed axial ratio of the
;                   best-fitting Gaussian components.
;               4) sol[3,*] = PA, is the observed position angle of the
;                   best-fitting Gaussian components. It is measured
;                   with respect to the same (arbitrary) system of angular
;                   coordinates one has adopted for the input vector Angle.
;
; EXAMPLE:
;       The sequence of commands below was used to generate the complete
;       MGE model of the Figures 11-12 in Cappellari (2002).
;       1) The FITS file is read, the bad regions are masked and
;           the image is sky subtracted;
;       2) the photometry along sectors is performed with SECTORS_PHOTOMETRY_TWIST;
;       3) the resulting measurements are fitted with MGE_FIT_SECTORS_TWIST;
;       4) the contour are printed on a PS file with MGE_PRINT_CONTOURS_TWIST.
;       The geometric parameters of the galaxy (eps,ang,xc,yc) were
;       previously determined using FIND_GALAXY.
;
;           fits_read, 'ngc5831_f702w_mosaic.fits', img, h
;
;           badpixels = where(img le 0) ; mask the unexposed part of the mosaic
;           ;
;           ; Mask a galaxy and two bright stars
;           ;
;           s = SIZE(img)
;           dist_circle, r, s[1:2], 1408.09, 357.749
;           badpixels = [badpixels,WHERE(r LT 200)]
;           dist_circle, r, s[1:2], 539.243, 161.264
;           badpixels = [badpixels,WHERE(r LT 50)]
;           dist_circle, r, s[1:2], 393.028, 713.681
;           badpixels = [badpixels,WHERE(r LT 50)]
;
;           skylevel = 13.0 ; In counts
;           img = img - skyLevel
;           scale = 0.100   ; WFPC2/WF CCD pixel scale arcsec/pixels
;
;           sigmaPSF = 0.4  ; WF pixels: here I use one single Gaussian PSF
;           eps = 0.28      ; These values were measured with FIND_GALAXY
;           ang = 141.0     ; major axis PA in the inner regions
;           xc = 969
;           yc = 974
;
;           sectors_photometry_twist, img, ang, xc, yc, radius, counts, angle, $
;               BADPIXELS=badpixelsMosaic, MINLEVEL=2.0
;
;           MGE_fit_sectors_twist, radius, angle, counts, eps, $
;               NGAUSS=11, SIGMAPSF=sigmaPSF, SOL=sol, /PRINT, $
;               SCALE=scale, QBOUNDS=[0.2,1]
;
;           MGE_print_contours_twist, img, ang, xc, yc, sol, $
;               FILE='NGC5831_Mosaic.ps', BINNING=5, MAGRANGE=9.5, $
;               SIGMAPSF=sigmaPSF, SCALE=scale
;
; PROCEDURES USED:
;       The following procedures are contained in the main MGE_FIT_SECTORS program.
;           MGE_FIT_SECTORS_TWIST_PRINT  -- called during the fit to show the fit profiles
;           FITFUNC_MGE_SECTORS_TWIST    -- returns the residuals between data and model
;
;       Other IDL routines needed:
;           BVLS  -- Michele Cappellari porting of Lawson & Hanson generalized NNLS
;                    http://www.strw.leidenuniv.nl/~mcappell/idl/
;           MGE_MPFIT -- Craig Markwardt porting of Levenberg-Marquardt MINPACK-1
;
; MODIFICATION HISTORY:
;       V1.0 First implementation, Padova, February 1999, Michele Cappellari
;       V2.0 Major revisions, Leiden, January 2000, MC
;       V3.0 Significant changes, Padova, July 2000, MC
;       V3.1 More robust definition of err in FITFUNC_MGE_SECTORS,
;           Leiden, 27 April 2001, MC
;       V3.2 Graphical changes: always show about 7 sectors on the screen,
;           and print plots with shared axes. Leiden, 8 July 2001, MC
;       V3.3 Added MGE PSF convolution, central pixel integration and changed
;           program input parameters to make it independent from SECTORS_PHOTOMETRY
;           Leiden, 26 July 2001, MC
;       V3.4 Added varying position angles fit, August 2001, MC
;       V3.5 Added /FASTNORM keyword, Leiden, 20 September 2001, MC
;       V3.6 Updated documentation, Leiden, 8 October 2001, MC
;       V3.7 Modified implementation with /NEGATIVE keyword
;           Leiden, 23 October 2001, MC
;       V3.8 Added explicit stepsize (STEP) of numerical derivative in
;           parinfo structure, after suggestion by Craig B. Markwardt.
;           Leiden, 23 February 2002, MC
;       V3.81 Added compilation options, Leiden 20 May 2002, MC
;       V3.82 Corrected checking of input parameters when SIGMAPSF is a vector.
;           Added ERRMSG keyword to MPFIT call. Leiden, 13 October 2002, MC
;       V3.83 Force the input parameters to the given bounds if they
;           fall outside the required range before starting the fit.
;           Leiden, 7 March 2003, MC
;       V3.84: Use N_ELEMENTS instead of KEYWORD_SET to test
;           non-logical keywords. Leiden, 9 May 2003, MC
;       V3.85: Use updated calling sequence for BVLS. Leiden, 20 March 2004, MC
;       V3.9: Force the surface brightness of the MGE model to decrease at
;           least as R^-2 at the largest measured radius. Leiden, 8 May 2004, MC
;       V3.91: Make sure this routine uses the Nov/2000 version of Craig Markwardt
;           MPFIT which was renamed MGE_MPFIT to prevent potential conflicts with
;           more recent versions of the same routine. Vicenza, 23 August 2004, MC.
;       V3.92: Allow forcing the outer slope of the surface brightness profile of
;           the MGE model to decrease at least as R^-n at the largest measured
;           radius (cfr. version 3.9). Leiden, 23 October 2004, MC
;       V3.93 Replaced LOGRANGE keyword in example with the new MAGRANGE.
;           MC, Leiden, 1 May 2005
;       V3.94: Changed axes labels in plots. Leiden, 18 October 2005, MC
;       V3.95: Use more robust la_least_squares (IDL 5.6) instead of SVDC with
;           /NEGATIVE keyword. MC, Oxford, 16 May 2008
;       V3.96: Force Gaussians smaller than the PSF, which have a degenerate
;           axial ratio, to have the same axial ratio as the mean of the first 
;           two well determined Gaussians. MC, Oxford, 24 September 2008
;-
;----------------------------------------------------------------------------
PRO mge_fit_sectors_twist_print, MYFUNCT, pars, iter, fnorm, $
    FUNCTARGS=fcnargs, PARINFO=parinfo, QUIET=quiet, DOF=dof
;
; This is a plotting routine that is called every NPRINT iterations of MPFIT
;
COMPILE_OPT IDL2, HIDDEN

COMMON mge_fit_sectors_twist, radius, counts, angle, gauss, soluz, $
    yfit, err, sigmaPSF, normPSF, sectors, negative, scale

PRINT, 'Iteration: ', STRTRIM(iter,2), ';  chi2: ', STRTRIM(fnorm,2)
;
; Select an x and y plot range that is the same for all plots
;
minrad = MIN(radius*scale,MAX=maxrad)
mincnt = MIN(counts,MAX=maxcnt)
xrange = minrad * (maxrad/minrad)^[-0.02,+1.02]
yrange = mincnt * (maxcnt/mincnt)^[-0.05,+1.05]
xtxt = minrad * (maxrad/minrad)^0.96
ytxt = mincnt * (maxcnt/mincnt)^0.87

xtickformat = '(A1)'
xtitle = ' '
ngauss = N_ELEMENTS(pars)/3
n = N_ELEMENTS(sectors)
dn = ROUND(n/6.0)

!P.MULTI=[0,2,(n-1)/dn+1]
!Y.MARGIN=[0,0] ; show plots with shared X axis
!Y.OMARGIN=[5,5] ; allow for space for the axis labels

FOR j=0,n-1,dn DO BEGIN
    w = WHERE(angle EQ sectors[j])
    w = w[SORT(radius[w])]
    if (j eq (n-1)/dn*dn) then begin ; Integer division here!
        xtickformat = !x.tickformat
        xtitle = 'arcsec'
    endif
    PLOT, radius[w]*scale, counts[w], /YLOG, /XLOG, /XSTYLE, /YSTYLE, $
        XTITLE=xtitle, YTITLE='counts', XTICKFORMAT=xtickformat, $
        XRANGE=xrange, YRANGE=yrange, CHARSIZE=2, TICKLEN=0.04, PSYM=6, SYMSIZE=0.5
    OPLOT, radius[w]*scale, yfit[w]
    str = STRING(sectors[j],FORMAT='(G0.3)')
    XYOUTS, xtxt, ytxt, str, ALIGN=1
    FOR k=0,ngauss-1 DO OPLOT, radius[w]*scale, gauss[w,k]*soluz[k]
    PLOT, radius[w]*scale, err[w]*100, /XLOG, /YSTYLE, /XSTYLE, $
        XTITLE=xtitle, YTITLE='error (%)', XTICKFORMAT=xtickformat, $
        XRANGE=xrange, YRANGE=[-20,19.9], CHARSIZE=2, TICKLEN=0.04, YMINOR=5, $
        PSYM=6, SYMSIZE=0.5
    PLOTS, xrange, [0,0], LINE=1
ENDFOR

!P.MULTI=0
!Y.MARGIN=[4,2] ; back to default values
!Y.OMARGIN=[0,0]

END
;----------------------------------------------------------------------------
FUNCTION fitfunc_mge_sectors_twist, pars
;
COMPILE_OPT IDL2, HIDDEN

COMMON mge_fit_sectors_twist

ngauss = N_ELEMENTS(pars)/3
npoints = N_ELEMENTS(radius)
logsigma = pars[0:ngauss-1]
q = pars[ngauss:2*ngauss-1]
pa = pars[2*ngauss:*]
gauss = DBLARR(npoints,ngauss)
A = gauss
sigma2 = 10d^(2d*logsigma)

w = WHERE(radius LT 0.5, m) ; Will perform flux integration inside the central pixel
FOR j=0,ngauss-1 DO BEGIN ; loop over the galaxy MGE Gaussians
    ang = (angle-pa[j])/!RADEG
    FOR k=0,N_ELEMENTS(sigmaPSF)-1 DO BEGIN ; loop over the PSF Gaussians
        ;
        ; Analytic convolution with a single Gaussian circular PSF
        ; Equations (4,5) in Cappellari (2002)
        ;
        sigmaX = SQRT(sigma2[j] + sigmaPSF[k]^2)
        sigmaY = SQRT(sigma2[j]*q[j]^2 + sigmaPSF[k]^2)
        ;
        ; Evaluate the normalized (volume=1) 2d Gaussian in polar coordinates
        ;
        g = EXP( -0.5d*radius^2 * ((COS(ang)/sigmaX)^2 + (SIN(ang)/sigmaY)^2) ) $
             / ( 2d*!DPI*sigmaX*sigmaY )
        ;
        ; Analytic integral of the Gaussian on the central pixel.
        ; Below we assume the central pixel is aligned with the galaxy axes.
        ; This is generally not the case, but the error due to this
        ; approximation is negligible in realistic situations.
        ;
        IF m NE 0 THEN $
            g[w] = ERRORF(2d^(-1.5d)/sigmaX) * ERRORF(2d^(-1.5d)/sigmaY)
        ;
        ; Each convolved MGE Gaussian is weighted with the
        ; normalization of the corresponding PSF component
        ;
        gauss[*,j] = gauss[*,j] + normPSF[k] * g
    ENDFOR
    A[*,j] = gauss[*,j]/counts ; gauss*SQRT(weights) = gauss/y
ENDFOR
B = REPLICATE(1d,npoints)  ; y*SQRT(weights) = 1 <== weights = 1/sigma^2 = 1/y^2

IF negative THEN $  ; Solution by LAPACK linear least-squares
    soluz = la_least_squares(transpose(A),B,METHOD=3,STATUS=ierr) $
ELSE BEGIN        ; Solution by Bounded-Variables Least-Squares
    BND = DBLARR(2,ngauss)
    BND[0,*] = 0d                       ; The usual positivity constraint
    BND[1,*] = (MACHAR()).XMAX
    ;NNLS_EXTERNAL, A, b, soluz, rnorm & ierr=0 ; If available call the compiled Fortran 77 routine
    BVLS, A, B, BND, soluz, IERR=ierr, ITMAX=ngauss*30, FASTNORM=fastNorm
ENDELSE
IF ierr NE 0 THEN MESSAGE, 'BVLS/la_least_squares Error n. ' + STRTRIM(ierr,2) + ': try to decrease NGAUSS'

yfit = gauss # soluz   ; Evaluate predictions by matrix multiplications
err = 1d - yfit/counts ; relative error: yfit, counts are positive quantities
chi2 = TOTAL(err^2) ; rnorm^2 = TOTAL(err^2) (this value is only used with the /LINEAR keyword)

RETURN, FLOAT(err) ; err is a vector. Solve BVLS in DOUBLE but MPFIT in FLOAT
END
;----------------------------------------------------------------------------
PRO MGE_fit_sectors_twist, radius1, angle1, counts1, eps, $
    NGAUSS=ngauss, NEGATIVE=negative1, SIGMAPSF=sigmaPSF1, NORMPSF=normPSF1, $
    SCALE=scale1, RBOUNDS=rbounds, QBOUNDS=qbounds, PABOUNDS=PAbounds, $
    SOL=sol, PRINT=print, OUTER_SLOPE=outer_slope
;
; This is the main routine. It is the only one that has to be called from outside.
;
COMPILE_OPT IDL2
ON_ERROR, 2

COMMON mge_fit_sectors_twist

t = SYSTIME(1)
;
; Starting from here the biggest part of this routine is devoted to input
; checking and loading of shared parameters into a common block
;
IF (WHERE(counts1 LE 0))[0] NE -1 THEN MESSAGE, 'Error: negative input counts'
n = N_ELEMENTS(radius1)
IF (N_ELEMENTS(angle1) NE n OR N_ELEMENTS(counts1) NE n) THEN $
    MESSAGE, 'Error: Input vectors must have the same length'
IF N_ELEMENTS(outer_slope) EQ 0 THEN outer_slope = 2.0 $
    ELSE outer_slope = outer_slope > 1.0 < 4.0

; load data vectors into the COMMON block
;
radius = radius1
counts = counts1
angle = angle1

IF N_ELEMENTS(ngauss) EQ 0 THEN ngauss = 12
IF KEYWORD_SET(negative1) THEN negative = 1 ELSE negative = 0
IF N_ELEMENTS(sigmaPSF1) EQ 0 THEN sigmaPSF = 0.0 ELSE sigmaPSF = sigmaPSF1
nPsf = N_ELEMENTS(sigmaPSF)
IF nPsf EQ 1 THEN normPSF = 1.0 ELSE BEGIN
    normPSF = normPSF1
    IF N_ELEMENTS(normPSF) NE nPsf THEN $
        MESSAGE, 'sigmaPSF and normPSF must have the same length'
    IF ROUND(TOTAL(normPSF)*100.0) NE 100 THEN MESSAGE, 'Error: PSF not normalized'
ENDELSE
IF N_ELEMENTS(scale1) EQ 0 THEN scale = 1.0 ELSE scale = scale1

sectors = angle[UNIQ(angle, SORT(angle))] ; Finds the different position angles
;
; Open grid in the range [rminLog,rmaxLog]
; The logarithmic slope of a Gaussian at R = sqrt(n)*sigma is -n.
; Below we constrain the largest Gaussian to have sigma < rmax/sqrt(n)
; to force the surface brightness of the MGE model to decrease
; at least as fast as R^-n at the largest measured radius.
;
IF N_ELEMENTS(rbounds) EQ 0 THEN BEGIN
    rmin = MIN(radius, MAX=rmax)
    rminLog = ALOG10(rmin)
    rmaxLog = ALOG10(rmax/sqrt(outer_slope))
ENDIF ELSE BEGIN
    rminlog = ALOG10(rbounds[0])
    rmaxlog = ALOG10(rbounds[1])
ENDELSE
IF N_ELEMENTS(qBounds) EQ 0 THEN qbounds=[0.05,1.0] ; no gaussians flatter than q=0.05
IF N_ELEMENTS(paBounds) EQ 0 THEN paBounds=[-90.0,90.0]
;
; If the smallest intrinsic Gaussian has sigma=0.75*sigmaPSF it will produce an
; observed Gaussian with sigmaObs=SQRT(sigma^2+sigmaPSF^2)=1.25*sigmaPSF.
; We require the sigma of the Gaussians to be larger than 0.75*sigmaPSF,
; or the smallest measured radius, whichever is larger.
;
IF TOTAL(sigmaPSF) GT 0 THEN rminLog = rminLog > ALOG10(0.75*MIN(sigmaPSF))

IF N_ELEMENTS(sol) EQ 0 THEN BEGIN
    logsigma = ((rmaxLog-rminLog)*(0.5+FINDGEN(ngauss))/ngauss + rminLog) ; open grid
    pars = [logsigma, REPLICATE(1.0-eps > qBounds[0] < qBounds[1], ngauss), $
            REPLICATE(0.0 > paBounds[0] < paBounds[1], ngauss)]
ENDIF ELSE BEGIN
    ngauss = N_ELEMENTS(sol)/4
    pars = FLTARR(ngauss*3)
    pars[0:ngauss-1] = ALOG10(sol[1,*]) > rminLog < rmaxLog      ; Log(sigma)
    pars[ngauss:2*ngauss-1] = sol[2,*] > qBounds[0] < qBounds[1] ; qObs
    pars[2*ngauss:*] = sol[3,*] > paBounds[0] < paBounds[1]      ; PA
ENDELSE
parinfo = REPLICATE({limited:[1,1], limits:[rminLog,rmaxLog], step:0.01}, 3*ngauss)
parinfo[ngauss:2*ngauss-1].limits = qBounds
parinfo[2*ngauss:*].limits = paBounds
parinfo[2*ngauss:*].step = 1.0 ; degrees: step for numerical derivative

sol = mge_mpfit('fitfunc_mge_sectors_twist', pars, BESTNORM=bestnorm, NITER=niter, $
    ITERPROC='mge_fit_sectors_twist_print', NPRINT=5, PARINFO=parinfo, ERRMSG=errmsg)
IF errmsg NE '' THEN MESSAGE, errmsg
;
; Print a PS file with the best fit to the profiles
;
IF KEYWORD_SET(print) THEN BEGIN
    LOADCT, 2
    SET_PLOT, 'PS'
    DEVICE, XSIZE=19, YSIZE=28, YOFFSET=1, XOFFSET=1, /COLOR, BITS=8
    mge_fit_sectors_twist_print, MYFUNCT, sol, niter, bestnorm
    DEVICE,/CLOSE
    CASE !VERSION.OS_FAMILY OF
        'Windows': SET_PLOT, 'WIN'
        'MacOS': SET_PLOT, 'MAC'
    ELSE: SET_PLOT, 'X'
    ENDCASE
ENDIF
;
; Print the results for the nonzero Gaussians
;
logSigma = sol[0:ngauss-1]
q = sol[ngauss:2*ngauss-1]
pa = sol[2*ngauss:*]
w = WHERE(soluz NE 0, m)
logSigma = logSigma[w]
q = q[w]
pa = pa[w]
soluz = soluz[w]
j = SORT(logSigma) ; Sort by increasing sigma

sol = FLTARR(4,m)
sol[0,*] = soluz[j]
sol[1,*] = 10.0^logSigma[j]
sol[2,*] = q[j]
sol[3,*] = pa[j]

; Force Gaussians smaller than the PSF, which have a degenerate
; axial ratio, to have the same axial ratio as the mean of the 
; first two well determined Gaussians. Flux is conserved by PSF
; convolution so no other changes are required
;
w1 = where(sol[1,*] le min(sigmaPSF), COMPLEMENT=w2,p)
if p gt 0 then sol[2,w1] = mean(sol[2,w2[0:1]])

PRINT, '####################################################'
PRINT, '  Computation time: ', STRTRIM(SYSTIME(1)-t,2), ' seconds'
PRINT, '  Total Iterations: ', STRTRIM(niter,2)
PRINT, ' Nonzero Gaussians: ', STRTRIM(m,2)
PRINT, '  Unused Gaussians: ', STRTRIM(ngauss-m,2)
PRINT, ' Sectors used in the fit: ', STRTRIM(N_ELEMENTS(sectors),2)
PRINT, ' Total number of points fitted: ', STRTRIM(n,2)
PRINT, ' Chi2: ', STRTRIM(bestnorm,2)
PRINT, ' STDEV: ', STRTRIM(STDDEV(err),2)
PRINT, ' MEANABSDEV: ', STRTRIM(MEANABSDEV(err),2)
PRINT, ' Total_Counts   Sigma_Pixels    qObs          PA'
PRINT, '####################################################'
PRINT, sol
PRINT, '++++++++++++++++++++++++++++++++++++++++++++++++++++'

END
;----------------------------------------------------------------------------
