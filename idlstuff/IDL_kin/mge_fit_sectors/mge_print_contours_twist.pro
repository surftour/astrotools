;######################################################################
;
; Copyright (C) 1999-2005, Michele Cappellari
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
;     MGE_PRINT_CONTOURS_TWIST
;
; AUTHOR:
;       Michele Cappellari, Astrophysics Sub-department, University of Oxford, UK
;
; PURPOSE:
;       Produces a contour plot comparing a convolved
;       MGE model to the original fitted image.
;
; CALLING SEQUENCE:
;       MGE_PRINT_CONTOURS_TWIST, Img, Ang, Xc, Yc, Sol,
;           FILE=file, MAGRANGE=magRange, SCALE=scale, BINNING=binning,
;           SIGMAPSF=[sigma1,sigma2,...], NORMPSF=[norm1,norm2,...], MODEL=model
;
; INPUTS:
;       Img = array containing the image that was fitted by MGE_FIT_SECTORS
;       Ang = Scalar defining the direction of the zero Position Angle of the
;               Gaussians. This is measured counterclockwise from the image Y
;               axis to the Gaussians major axis, as measured by FIND_GALAXY.
;       Xc = Scalar giving the common X coordinate in pixels of the
;               center of the Gaussians.
;       Yc = Scalar giving the common Y coordinate in pixels of the
;               center of the Gaussians.
;       SOL - Array containing a 4xNgauss array with the MGE best-fitting
;               solution as produced by MGE_FIT_SECTORS_TWIST:
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
;                   with respect to the PA given by the input scalar Ang.
;
; OUTPUTS:
;       No output parameters. The results are plotted in a PS file.
;       But see also the optional MODEL keyword.
;
; OPTIONAL KEYWORDS:
;       FILE - String: name of the output PS file.
;       MAGRANGE - Scalar giving the range in magnitudes of the equally
;               spaced contours to plot, in steps of 0.5 mag/arcsec^2,
;               starting from the model maximum value.
;               (default: from maximum to minimum of the model)
;       SCALE - The pixel scale in arcsec/pixels used for the plot axes.
;               (default: 1)
;       BINNING - Pixels to bin together before plotting.
;               Helps producing MUCH smaller PS files (default: no binning).
;       /CONVOL - Set this keyword to use Gaussian convolution instead of the 
;               default SMOOTH routine to smooth the galaxy image before rebining 
;               (when the BINNING keyword is used). Noise suppression is in principle 
;               more effective when this keyword is set but the computation 
;               time and memory requirement can be significantly larger. 
;       SIGMAPSF - Scalar giving the sigma of the PSF, or vector with the
;               sigma of an MGE model for the circular PSF. (Default: no convolution)
;       NORMPSF - This is optional if only a scalar is given for SIGMAPSF,
;               otherwise it must contain the normalization of each MGE component
;               of the PSF, whose sigma is given by SIGMAPSF. The vector needs to
;               have the same number of elements of SIGMAPSF and the condition
;               TOTAL(normpsf) = 1 must be verified. In other words the MGE PSF
;               needs to be normalized. (default: 1).
;       MODEL -- Named variable that will contain in output an image with
;               the MGE model, convolved with the given PSF.
;       _EXTRA - Additional parameters can be passed to the IDL CONTOUR
;               procedure via the _EXTRA mechanism
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
;       The following procedures are contained in the main MGE_PRINT_CONTOURS program.
;           GAUSS2D_MGE  -- Returns a 2D Gaussian image
;           MULTI_GAUSS_TWIST  -- Returns a 2D MGE expansion by calling GAUSS2D_MGE
;       Astronomy User's Library (http://idlastro.gsfc.nasa.gov/) routines needed:
;           CONVOLVE  -- Fast convolution using FFT
;
; MODIFICATION HISTORY:
;       V1.0 First implementation, Padova, February 1999, Michele Cappellari
;       V2.0 Major revisions, Leiden, January 2000, MC
;       V2.1 Updated documentation, Leiden, 8 October 2001, MC
;       V2.2 Implemented MGE PSF, Leiden, 29 October 2001, MC
;       V2.3 Added MODEL keyword, Leiden, 30 October 2001, MC
;       V2.31 Added compilation options, Leiden 20 May 2002, MC
;       V2.32: Use N_ELEMENTS instead of KEYWORD_SET to test
;           non-logical keywords. Leiden, 9 May 2003, MC
;       V2.4: Convolve image with a Gaussian kernel instead of using
;           the SMOOTH function before binning. Always shows contours
;           in steps of 0.5 mag/arcsec^2. Replaced LOGRANGE and NLEVELS
;           keywords with MAGRANGE. Leiden, 30 April 2005, MC
;       V2.41: Added /CONVOL keyword. MC, Oxford, 23 September 2008
;-
;----------------------------------------------------------------------------
FUNCTION gauss2d_mge, n, xc, yc, sx, sy, pos_ang
;
; Returns a 2D Gaussian image with size N[0]xN[1], center (XC,YC),
; sigma (SX,SY) along the principal axes and position angle POS_ANG, measured
; from the positive Y axis to the Gaussian major axis (positive counter-clockwise).
;
COMPILE_OPT IDL2, HIDDEN

ang = (pos_ang-90.0)/!RADEG

x = FINDGEN(n[0]) - xc
y = FINDGEN(n[1]) - yc

xcosang = COS(ang)/(SQRT(2.0)*sx) * x
xsinang = SIN(ang)/(SQRT(2.0)*sy) * x
ycosang = COS(ang)/(SQRT(2.0)*sy) * y
ysinang = SIN(ang)/(SQRT(2.0)*sx) * y

im = FLTARR(n[0], n[1], /NOZERO)

FOR i=0,n[1]-1 DO BEGIN
   xtemp =  xcosang + ysinang[i]
   ytemp = -xsinang + ycosang[i]
   im[0,i] = xtemp^2 + ytemp^2
ENDFOR

RETURN, EXP(-im)
END
;----------------------------------------------------------------------------
FUNCTION multi_gauss_twist, pars, img, sigmaPSF, normPSF, xpeak, ypeak, theta
;
COMPILE_OPT IDL2, HIDDEN

s = SIZE(img)
u = 0.0

ngauss = N_ELEMENTS(pars)/4
lum = pars[0,*]
sigma = pars[1,*]
q = pars[2,*]
pa = pars[3,*]
;
; Analytic convolution with an MGE circular Gaussian
; Eq.(4,5) in Cappellari (2002)
;
FOR j=0,ngauss-1 DO BEGIN
    FOR k=0,N_ELEMENTS(sigmaPSF)-1 DO BEGIN
        qx = SQRT(sigma[j]^2 + sigmaPSF[k]^2)
        qy = SQRT((sigma[j]*q[j])^2 + sigmaPSF[k]^2)
        g = gauss2d_mge(s[1:2], xpeak, ypeak, qx, qy, theta-pa[j])
        u = TEMPORARY(u) + lum[j]*normPSF[k]/(2d*!DPI*qx*qy) * g
    ENDFOR
ENDFOR

RETURN, u
END
;----------------------------------------------------------------------------
PRO MGE_print_contours_twist, img, ang, xc, yc, sol, $
        FILE=file, MAGRANGE=magRange, SCALE=scale, BINNING=binning, $
        SIGMAPSF=sigmaPSF, NORMPSF=normPSF, MODEL=model, $
        CONVOL=convol, _EXTRA=extra

COMPILE_OPT IDL2
ON_ERROR, 2

IF N_ELEMENTS(sigmaPSF) EQ 0 THEN sigmaPSF = 0.0 ; Default no convolution
n = N_ELEMENTS(sigmaPSF)
IF n EQ 1 THEN normPSF = 1.0 ELSE BEGIN
    IF N_ELEMENTS(normPSF) NE n THEN $
        MESSAGE, 'sigmaPSF and normPSF must have the same length'
    IF ROUND(TOTAL(normPSF)*100.0) NE 100 THEN MESSAGE, 'Error: PSF not normalized'
ENDELSE
IF N_ELEMENTS(file) EQ 0 THEN file = 'idl.ps'

SET_PLOT,'ps'
DEVICE, XSIZE=20, YSIZE=28, XOFFSET=0.5, YOFFSET=1, FONT_SIZE=12, FILE=file

s = SIZE(img)
model = multi_gauss_twist(sol, img, sigmaPSF, normPSF, xc, yc, ang)
IF N_ELEMENTS(magRange) EQ 0 THEN magRange = 10.0
levels = REVERSE(2.5*ALOG10(model[xc,yc]) - 0.5*FINDGEN(magRange*2)) ; 0.5 mag/arcsec^2 steps

IF N_ELEMENTS(scale) EQ 0 THEN scale = 1.0
IF N_ELEMENTS(binning) GT 0 THEN BEGIN
    dx = 3*binning ; Size of the Gaussian kernel for convolution. Adopt FWHM=binning
    kernel = gauss2d_mge([dx,dx], dx/2.0, dx/2.0, binning/2.35, binning/2.35, 0.0)
    kernel = kernel/total(kernel) ; Normalize the kernel
    if keyword_set(convol) then begin
        model1 = frebin(convolve(model,kernel),s[1]/binning,s[2]/binning)
        img1 = frebin(convolve(img,kernel),s[1]/binning,s[2]/binning)
    endif else begin
        model1 = frebin(smooth(model,binning),s[1]/binning,s[2]/binning)
        img1 = frebin(smooth(img,binning),s[1]/binning,s[2]/binning)
    endelse
ENDIF ELSE BEGIN
    binning = 1.0
    model1 = model
    img1 = img
ENDELSE

CONTOUR, 2.5*ALOG10(model1), XTITLE='arcsec', YTITLE='arcsec', $
    FINDGEN(s[1]/binning)*scale*binning, FINDGEN(s[2]/binning)*scale*binning, $
    /ISOTROPIC, /XSTYLE, /YSTYLE, LEVELS=levels, _EXTRA=extra
CONTOUR, 2.5*ALOG10(img1), LEVELS=levels, /OVERPLOT, $
    FINDGEN(s[1]/binning)*scale*binning, FINDGEN(s[2]/binning)*scale*binning
IF s[1] EQ 1600 AND s[2] EQ 1600 THEN $ ; assumes WFPC2 and overplot the mosaic
    PLOTS,[3,155,155,114,114,80,80,3,3], [4,4,80,80,113.5,113.5,155,155,4]

DEVICE,/CLOSE

CASE !VERSION.OS_FAMILY OF
    'Windows': SET_PLOT, 'WIN'
    'MacOS': SET_PLOT, 'MAC'
ELSE: SET_PLOT, 'X'
ENDCASE

END
;----------------------------------------------------------------------------
