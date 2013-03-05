;######################################################################
;
; Copyright (C) 1999-2001, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; For details on the method see:
;   Cappellari M., 2002, MNRAS, 333, 400
;
; Updated versions of the software are available from my web page
; http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
; If you have found this software useful for your
; research, I would appreciate an acknowledgment to
; `use of the MGE fitting software developed by Cappellari (2002)'.
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
;       WFPC2_MGE_FIT
;
; AUTHOR:
;       Michele Cappellari, Astrophysics Sub-department, University of Oxford, UK
;
; PURPOSE:
;       Perform an MGE fit to a galaxy observed with the WFPC2 camera
;       of the Hubble Space Telescope. It is assumed that the galaxy nucleus
;       was centered on the PC1 CCD and the MGE model is determined by fitting
;       both the PC1 CCD and the whole WFPC2 mosaic at the same time.
;       This program is provided as an usage example of the MGE_FIT_SECTORS
;       package. It essentially passes parameters to the proper function
;       of the package in the correct sequence.
;
; EXPLANATION:
;       Further information is available in
;       Cappellari M., 2002, MNRAS, 333, 400
;
; CALLING SEQUENCE:
;       WFPC2_MGE_FIT, Eps, Ang, Xc, Yc,
;           NGAUSS=Ngauss, N_SECTORS=n_sectors,
;           SIGMAPSF=[sigma1,sigma2,...], NORMPSF=[norm1,norm2,...],
;           SKYLEVEL=Skylevel, MINLEVEL=Minlevel,
;           BADPIXELSPC=badpixelsPC, BADPIXELSMOSAIC=badpixelsMosaic,
;           IMGPC=imgPC, IMGMOSAIC=imgMosaic, MAGRANGE=magRange,
;           FILENAME=filename, HEADER=header, SOL=sol
;
; INPUTS:
;       Eps = The galaxy "average" ellipticity Eps = 1-b/a = 1-q'.
;           Photometry will be measured along elliptical annuli with
;           constant axial ellipticity Eps. The four quantities
;           (Eps,Ang,Xc,Yc) can be measured with the routine FIND_GALAXY.
;       Ang = Position angle measured counterclockwise from
;           the image Y axis, to the galaxy major axis.
;       Xc = X coordinate of the galaxy center in pixels, as measured on
;           the PC1 trimmed image (read with WFPC2_READ, filename, img ,/TRIM)
;       Yc = Y coordinate of the galaxy center in pixels on the PC1 CCD.
;
; OUTPUTS:
;       No output parameters. The results are printed on the screen,
;       and plotted in a PS file.
;
; OPTIONAL INPUT KEYWORDS:
;       NGAUSS - Number of Gaussians desired in the MGE fit.
;               Typical values are in the range 10-20 when the /LINEAR
;               keyword is NOT set (default: 12) and 20^2-40^2 when the
;               /LINEAR keyword is set (default: 30^2).
;       N_SECTORS - Number of the sectors, equally spaced in angle, from
;               the galaxy major axis to the minor axis (one quadrant).
;               (default: 19 to cover the whole image with 5 degrees width).
;       SIGMAPSF - Scalar giving the sigma of the PSF of the high resolution
;               image (see Cappellari [2002] for details), or vector with the
;               sigma of an MGE model for the circular PSF. (Default: no convolution)
;       NORMPSF - This is optional if only a scalar is given for SIGMAPSF,
;               otherwise it must contain the normalization of each MGE component
;               of the PSF, whose sigma is given by SIGMAPSF. The vector needs to
;               have the same number of elements of SIGMAPSF and the condition
;               TOTAL(normpsf) = 1 must be verified. In other words the MGE PSF
;               needs to be normalized. (default: 1).
;       SKYLEVEL - Sky level in units of counts of the PC1 (!) image.
;       MINLEVEL - Scalar giving the minimum COUNTS level to include
;               in the photometry, in units of the PC1 CCD. The measurement
;               along one profile will stop when the counts go below this level.
;       BADPIXELSPC - Vector containing the indices of the pixels in the
;               PC1 image that have to be excluded from the photometry
;               calculation.
;       BADPIXELSMOSAIC - Vector containing the indices of the pixels in the
;               mosaic image that have to be excluded from the photometry
;               calculation.
;       IMGPC - The PC1 image as read with 'WFPC2_READ, filename, img ,/TRIM'
;       IMGMOSAIC - The WFPC2 mosaic image as read with
;               'WFPC2_READ, filename, img ,/BATWING'
;       MAGRANGE - Scalar giving the range in magnitudes of the equally
;               spaced contours to plot, in steps of 0.5 mag/arcsec^2,
;               starting from the model maximum value.
;               (default: from maximum to minimum of the model)
;       FILENAME - Name of the FITS file containing the WFPC2 image in
;               the standard STScI 4x800x800 format. This keyword has to be
;               supplied if IMGPC and IMGMOSAIC are not given.
;       HEADER - If FILENAME is not given this parameter has to be supplied,
;               containing the header of the WFPC2 FITS image.
;       SOL - Array containing a 3xNgauss array with the MGE best-fitting
;               solution as produced by MGE_FIT_SECTORS:
;               1) sol[0,*] = TotalCounts, of the Gaussians components.
;                   the relation TotalCounts = Height*(2*!PI*Sigma^2*qObs)
;                   can be used compute the Gaussian central surface
;                   brightness (Height)
;               2) sol[1,*] = sigma, is the dispersion of the best-fitting
;                   Gaussians in pixels.
;               3) sol[2,*] = qObs, is the observed axial ratio of the
;                   best-fitting Gaussian components.
;       _EXTRA - Any additional accepted parameter can be passed to the
;               procedures SECTORS_PHOTOMETRY, MGE_FIT_SECTORS and
;               MGE_PRINT_CONTOURS via the IDL _EXTRA mechanism.
;
; EXAMPLE:
;       The following commands generate a fit to a WFPC2 image of NGC4473.
;       Both the PC and the WFPC2 mosaic are fitted at the same time.
;
;       These parameters are for the PC1 trimmed image as given by FIND_GALAXY:
;
;           skylevel = 3.6 ; PC1 counts
;           sigmaPSF = 0.8 ; pixels
;           eps = 0.42
;           ang = -79.4
;           xc = 378
;           yc = 376
;           ngauss = 15
;           minlevel = 0.0
;
;           WFPC2_MGE_fit, eps, ang, xc, yc, NGAUSS=ngauss, SIGMAPSF=sigmaPSF, $
;               SKYLEVEL=skylevel, MINLEVEL=minlevel, FILENAME='ngc4473_f814w.fits'
;
; PROCEDURES USED:
;       No additional routines are included in the main WFPC2_MGE_FIT program.
;
;       Other routines needed from the MGE_FIT_SECTORS package by
;       Michele Cappellari (http://www.strw.leidenuniv.nl/~mcappell/idl)
;           SECTORS_PHOTOMETRY       -- perform photometry along sectors
;           MGE_FIT_SECTORS          -- do the actual MGE fit
;           MGE_PRINT_CONTOURS       -- plot the results
;
;       Astronomy User's Library routines needed (http://idlastro.gsfc.nasa.gov):
;           WFPC2_READ
;           SXPAR
;
; MODIFICATION HISTORY:
;       V1.0 Michele Cappellari, Vicenza, 9 November 2000
;       V2.0 Updated MC, Leiden July 2001
;       V2.01 Added compilation options, MC, Leiden 20 May 2002
;       V2.02 Allow passing of *any* allowed keyword to SECTORS_PHOTOMETRY,
;           MGE_FIT_SECTORS and MGE_PRINT_CONTOURS via the _EXTRA keyword.
;           Use N_ELEMENTS instead of KEYWORD_SET to test non-logical keywords.
;           MC, Leiden, 9 May 2003
;       V2.03 Included documentation of new the MAGRANGE keyword.
;           MC, Leiden, 30 April 2005
;-
;----------------------------------------------------------------------------
PRO WFPC2_MGE_FIT, Eps, Ang, Xc1, Yc1, $
    BADPIXELSPC=badpixelsPC, BADPIXELSMOSAIC=badpixelsMosaic, $
    IMGPC=imgPC, IMGMOSAIC=imgMosaic, FILENAME=filename, SKYLEVEL=skylevel, $
    SIGMAPSF=sigmaPSF, HEADER=header, SOL=sol, _EXTRA=extra
;
; This procedure is specialized for the case of an MGE model of a galaxy
; that has been observed with the WFPC2 by placing the galaxy nucleus on the
; PC1 CCD. This is the most common situation for galaxy photometry with WFPC2.
;
; NB: (xc,yc) have to be measured on the trimmed PC1 image. With the /TRIM
; keyword the PC1 image has the same orientation as the whole WFPC2 mosaic.
; SKYLEVEL and MINLEVEL are also in the units for the PC pixels.
; The final photometry will be in PC pixels and PC counts/pixels.
;
COMPILE_OPT IDL2
ON_ERROR, 2

scale1 = 0.0455     ; (arcsec/pixel) PC1
scale2 = 0.100      ; (arcsec/pixel) WF
ratio = scale2/scale1
;
; Perform photometry on the PC1 image
;
IF N_ELEMENTS(imgPC) EQ 0 THEN $
    wfpc2_read, filename, imgPC, header, /trim
IF N_ELEMENTS(skyLevel) EQ 0 THEN skyLevel = 0.0
imgPC1 = imgPC - skylevel ; the input image is not modified
sectors_photometry, imgPC1, eps, ang, xc1, yc1, radius1, angle1, cnt1, $
    BADPIXELS=badpixelsPC, _EXTRA=extra
;
; Perform photometry on the whole mosaic WFPC2 image.
; Convert the coordinate of the nucleus from the PC1 to the mosaic coordinates.
;
xc2 = ROUND(799.5 + 0.4572 * xc1)
yc2 = ROUND(799.5 + 0.4572 * yc1)
IF N_ELEMENTS(imgMosaic) EQ 0 THEN $
    wfpc2_read, filename, imgMosaic, /batwing
IF N_ELEMENTS(badpixelsMosaic) EQ 0 THEN $ ; Do not change the actual input variable
    badpixelsMosaic1 = WHERE(imgMosaic LE 0.0) $
ELSE $
    badpixelsMosaic1 = [badpixelsMosaic, WHERE(imgMosaic LE 0.0)]
imgMosaic1 = imgMosaic/ratio^2 - skylevel ; the input image is not modified
sectors_photometry, imgMosaic1, eps, ang, xc2, yc2, radius2, angle2, cnt2, $
    BADPIXELS=badpixelsMosaic1, _EXTRA=extra
radius2 = radius2 * ratio ; Put all radii on the same scale
;
; Uses PC1 photometry only inside a 15" radius from the nucleus
; and Mosaic photometry outside this radius.
;
w1 = WHERE(radius1 LE 15.0/scale1)
w2 = WHERE(radius2 GT 15.0/scale1)
radius = [radius1[w1],radius2[w2]]
angle = [angle1[w1],angle2[w2]]
counts = [cnt1[w1],cnt2[w2]]

MGE_fit_sectors, radius, angle, counts, eps, $
    SIGMAPSF=sigmaPSF, SCALE=scale1, SOL=sol, _EXTRA=extra

exptime = sxpar(header, 'EXPTIME')
photmode = STRTRIM(sxpar(header, 'PHOTMODE'),2)
targname = STRTRIM(sxpar(header, 'TARGNAME'),2)
PRINT, 'EXPTIME:     ', exptime
PRINT, 'PHOTMODE:    ', photmode
PRINT, 'TARGNAME:    ', targname

sol1 = sol ; Do not change SOL since this has to be returned from this procedure
mge_print_contours, imgPC1, ang, xc1, yc1, sol1, $
    FILE=targname+'_PC.ps', SCALE=scale1, BINNING=2, $
    SIGMAPSF=sigmaPSF, TITLE=targname, _EXTRA=extra
ngauss = N_ELEMENTS(sol)/3
j = INDGEN(ngauss)*3
sol1[j] = sol1[j]/ratio^2    ; totalCounts scales as sigma^2
sol1[j+1] = sol1[j+1]/ratio  ; sigma
mge_print_contours, imgMosaic1, ang, xc2, yc2, sol1, $
    FILE=targname+'_Mosaic.ps', SCALE=scale2, BINNING=6, $
    SIGMAPSF=sigmaPSF/ratio, TITLE=targname, _EXTRA=extra

END
;----------------------------------------------------------------------------
