;#############################################################################
;
; Copyright (C) 2004-2008, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; Updated versions of the software are available from my web page
; http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
; If you have found this software useful for your research,
; I would appreciate an acknowledgment to the use of the
; `JAM modelling method of Cappellari (2008)'.
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;#############################################################################
;+
; NAME:
;   JAM_SPHERICAL_RMS
;
; PURPOSE:
;    This procedure calculates a prediction for the projected second velocity
;    moments V_RMS = sqrt(V^2 + sigma^2) for an anisotropic spherical galaxy
;    model. It implements the solution of the anisotropic Jeans equations
;    presented in equation (50) of Cappellari (2008, MNRAS, 390, 71). PSF
;    convolution is done as described in the Appendix of that paper.
;
; CALLING SEQUENCE:
;    JAM_SPHERICAL_RMS, $
;        surf_lum, sigma_lum, surf_pot, sigma_pot, mbh, distance, rad, sigModel, $
;        BETA=beta, NORMPSF=normPsf, NRAD=nrad, PIXSIZE=pixSize, $
;        SIGMAPSF=sigmaPsf, STEP=step
;
; INPUT PARAMETERS:
;   SURF_LUM: vector of length N containing the peak surface brightness of the
;       MGE Gaussians describing the galaxy surface brightness in units of
;       Lsun/pc^2 (solar luminosities per parsec^2).
;   SIGMA_LUM: vector of length N containing the dispersion in arcseconds of
;       the MGE Gaussians describing the galaxy surface brightness.
;   SURF_POT: vector of length M containing the peak value of the MGE Gaussians
;       describing the galaxy surface density in units of Msun/pc^2 (solar
;       masses per parsec^2). This is the MGE model from which the model
;       potential is computed.
;     - In a common usage scenario, with a self-consistent model, one will have
;       the same Gaussians for both the surface brightness and the potential.
;       This implies SURF_POT = SURF_LUM*(M/L) and SIGMA_POT = SIGMA_LUM, where
;       M/L is the desired global mass-to-light ratio of the model.
;       Alternatively one can adopt SURF_POT = SURF_POT as input and scale the
;       output SIGMODEL to fit the data with equation (51) after the
;       computation.
;     - IMPORTANT: when one scales SIGMODEL by a given M/L to fit the data, the
;       MBH that was adopted for the calculation is also effectively scaled by
;       the same M/L!
;   SIGMA_POT: vector of length M containing the dispersion in arcseconds of
;       the MGE Gaussians describing the galaxy surface density.
;   MBH: Mass of a nuclear supermassive black hole in solar masses.
;   DISTANCE: distance of the galaxy in Mpc.
;   RAD: Vector of length P with the radius from the galaxy center in
;       arcseconds of the bins (or pixels) at which one wants to compute the
;       model predictions.
;
; KEYWORDS:
;   BETA: Vector of length N with the anisotropy
;       beta = 1 - (sigma_theta/sigma_R)^2 of the individual MGE Gaussians. A
;       scalar can be used if the model has constant anisotropy.
;   NORMPSF: Vector of length Q with the fraction of the total PSF flux
;       contained in the various circular Gaussians describing the PSF of the
;       observations. It has to be total(NORMPSF) = 1. The PSF will be used for
;       seeing convolution of the model kinematics.
;   NRAD: Number of logarithmically spaced radial positions for which the
;       models is evaluated before interpolation and PSF convolution. One may
;       want to increase this value if the model has to be evaluated over many
;       orders of magnitutes in radius (default: NRAD=50).
;   PIXSIZE: Size in arcseconds of the (square) spatial elements at which the
;       kinematics is obtained. This may correspond to the size of the spaxel
;       or lenslets of an integral-field spectrograph. This size is used to
;       compute the kernel for the seeing and aperture convolution.
;     - If this is not set, or PIXSIZE = 0, then convolution is not performed.
;   SIGMAPSF: Vector of length Q with the dispersion in arcseconds of the
;       circular Gaussians describing the PSF of the observations.
;     - If this is not set, or SIGMAPSF = 0, then convolution is not performed.
;     - IMPORTANT: PSF convolution is done by creating a 2D image, with pixels
;       size given by STEP = MAX(SIGMAPSF,PIXSIZE/2)/4, and convolving it with
;       the PSF + aperture. If the input radii RAD are very large with respect
;       to STEP, the 2D image may require a too large amount of memory. If this
;       is the case one may compute the model predictions at small radii
;       separately from those at large radii, where PSF convolution is not
;       needed.
;   STEP: Spatial step for the model calculation and PSF convolution in arcsec.
;       This value is automatically computed by default as
;       STEP = MAX(SIGMAPSF,PIXSIZE/2)/4. It is assumed that when PIXSIZE or
;       SIGMAPSF are big, high resolution calculations are not needed. In some
;       cases however, e.g. to accurately estimate the central Vrms in a very
;       cuspy galaxy inside a large aperture, one may want to override the
;       default value to force smaller spatial pixels using this keyword.
;
; OUTPUT PARAMETER:
;   SIGMODEL: Vector of length P with the model predictions for the velocity
;       second moments (sigma in the spherical case) of each bin.
;
; USAGE EXAMPLE:
;    A simple usage example is given in the procedure TEST_JAM_SPHERICAL_RMS at
;    the end of this file.
;
; REQUIRED ROUTINES:
;       By M. Cappellari (included in the JAM distribution):
;       - ANY
;       - DIFF
;       - HYPERGEOMETRIC2F1
;       - IBETAM
;       - MESHGRID
;       - QUADVA
;       - RANGE
;
; MODIFICATION HISTORY:
; V1.0: Written and tested isotropic case.
;    Michele Cappellari, Vicenza, 10 August 2004
; V2.0: Include anisotropic case with 1D integral. MC, Oxford, 4 April 2008
; V3.1: First released version. MC, Oxford, 12 August 2008
; V3.2: Updated documentation. MC, Oxford, 14 August 2008
; V4.0: Implemented PSF convolution using interpolation on polar grid. Dramatic
;     speed-up of calculation. Further documentation.
;     MC, Oxford, 11 September 2008
; V4.01: Included keyword STEP. MC, Windhoek, 29 September 2008
; V4.02: Added keywords NRAD. Thanks to Michael Williams for reporting possible
;     problems with too coarse interpolation. MC, Oxford, 21 November 2008
; -
;#############################################################################
FUNCTION sphani_integrand_spherical_jeans, r, $
    SIG_L=sig_l, SIG_M=sig_m, LUM=lum, MASS=mass, MBH=Mbh, RMIN=rmin, BETA=beta
compile_opt idl2, hidden
;
; This routine implements the integrand of equation (50) of Cappellari (2008).
; The routine tries to speed up the calculation by treating differently
; the isotropic or constant-anisotropy cases.

mass_r = Mbh
FOR j=0L,n_elements(mass)-1 DO $
    mass_r += mass[j]*( ERF(r/(SQRT(2d)*sig_m[j])) $
            - r*SQRT(2d/!dpi)*EXP(-0.5d*(r/sig_m[j])^2)/sig_m[j] ) ; equation (49)

if array_equal(beta,beta[0]) then begin   ; Faster constant-anisotropy model
    if beta[0] eq 0 then $                ; Isotropic case
        er = sqrt(r^2-rmin^2) $           ; equation (44)
    else begin                            ; Anisotropic case
        rat = (rmin/r)^2
        er = 0.5d*rmin/rat^beta[0]* $     ; equation (43)
            ( beta[0]*ibetam(0.5d + beta[0],0.5d,rat) - ibetam(beta[0] - 0.5d,0.5d,rat) $
            + sqrt(!dpi)*(1.5d - beta[0])*gamma(beta[0] - 0.5d)/gamma(beta[0]) )
    endelse
    lum_dens = 0d
    FOR j=0L,n_elements(lum)-1 DO $
        lum_dens += lum[j]*EXP(-0.5d*(r/sig_l[j])^2)/(SQRT(2d*!dpi)*sig_l[j])^3
    fun = er * lum_dens
endif else begin                          ; Slower variable-anisotropy model
    fun = 0d
    rat = (rmin/r)^2
    for j=0,n_elements(lum)-1 DO begin
        if beta[j] eq 0 then $            ; Isotropic case
            er = sqrt(r^2-rmin^2) $       ; equation (44)
        else $                            ; Anisotropic case
            er = 0.5d*rmin/rat^beta[j]* $ ; equation (43)
                ( beta[j]*ibetam(0.5d + beta[j],0.5d,rat) - ibetam(beta[j] - 0.5d,0.5d,rat) $
                + sqrt(!dpi)*(1.5d - beta[j])*gamma(beta[j] - 0.5d)/gamma(beta[j]) )
        fun += er * lum[j]*EXP(-0.5d*(r/sig_l[j])^2)/(SQRT(2d*!dpi)*sig_l[j])^3
    endfor
endelse

; This routine returns a vector of values computed at different values of r
;
G = 0.00430237d    ; (km/s)^2 pc/Msun [6.674e-11 SI units]

return, 2d*G*fun*mass_r/r^2
END
;------------------------------------------------------------------
FUNCTION sphani_weighted_sigma2, R, s
compile_opt idl2, hidden
;
; Integration of equation (50)

s.rmin = R
rmax = 3d*max(s.sig_l)
if R ge rmax then message, 'R > rmax'
quadva, 'sphani_integrand_spherical_jeans', [R, rmax], int, $
	RELTOL=1d-5, ABSTOL=0d, FUNCTARGS=s

return, int
END
;------------------------------------------------------------------
function shpani_second_moment, R, s
compile_opt idl2, hidden
;
; This routine gives the second V moment after convolution with a PSF.
; The convolution is done using interpolation of the model on a
; polar grid, as described in Appendix A of Cappellari (2008).

if max(s.sigmaPsf) gt 0 and s.pixSize gt 0 then begin ; PSF convolution

    ; Kernel step is 1/4 of largest value between sigma(min) and 1/2 pixel side.
    ; Kernel half size is the sum of 3*sigma(max) and 1/2 pixel diagonal.
    ;
    if s.step gt 0 then step = s.step $
        else step = (s.pixSize/2d > min(s.sigmaPsf))/4d
    mx = 3d*max(s.sigmaPsf) + s.pixSize/sqrt(2d)

    ; Make grid linear in log of radius RR
    ;
    rmax = max(R) + mx ; Radius of circle containing all data + convolution
    nrad = s.NRAD ; radial elements
    logRad = range(alog(step),alog(rmax),nrad) ; Linear grid in log(RR)
    rr = exp(logRad)

    ; The model Vrms computation is only performed on the radial grid
    ; which is then used to interpolate the values at any other location
    ;
    wm2Pol = rr*0
    mgePol = wm2Pol
    for j=0,n_elements(rr)-1 do begin
        wm2Pol[j] = sphani_weighted_sigma2(rr[j],s)
        mgePol[j] = total( s.surf_l * exp(-0.5d*(rr[j]/s.sig_l)^2) )
    endfor

    nx = ceil(rmax/step/64d)*64  ; Make 2*nx dimensions divisible by 2^7
    x1 = range(-nx,nx,2*nx)*step ; for much faster FFT computation.
    meshgrid, x1, x1, xCar, yCar ; Cartesian grid for convolution

    ; Interpolate MGE model and Vrms over cartesian grid
    ;
    r1 = 0.5d*alog(xCar^2 + yCar^2) ; Log radius of cartesian grid
    indx = (r1 - logRad[0])/(logRad[1] - logRad[0]) ; Convert to indices for INTERPOLATE
    wm2Car = interpolate(wm2Pol,indx,CUBIC=-0.5)
    mgeCar = interpolate(mgePol,indx,CUBIC=-0.5)

    nk = ceil(mx/step)
    kgrid = range(-nk,nk,2*nk)*step
    meshgrid, kgrid, kgrid, xgrid, ygrid ; Kernel is square

    ; Compute kernel with equation (A6) of Cappellari (2008).
    ; Normaliztion is irrelevant here as it cancels out.
    ;
    kernel = 0d
    dx = s.pixSize/2d
    sp = sqrt(2d)*s.sigmaPsf
    for j=0,n_elements(s.sigmapsf)-1 do $
        kernel += s.normPsf[j] $
            * (erf((dx-xgrid)/sp[j]) + erf((dx+xgrid)/sp[j])) $
            * (erf((dx-ygrid)/sp[j]) + erf((dx+ygrid)/sp[j]))

    ; Seeing and aperture convolution with equation (A3)
    ;
    muCar = sqrt(convolve(wm2Car,kernel)/convolve(mgeCar,kernel))

    ; Interpolate convolved image at observed apertures.
    ; Aperture integration was already included in the kernel.
    ;
    indx = (R/sqrt(2) - x1[0])/(x1[1] - x1[0]) ; Convert to indices for INTERPOLATE
    mu = interpolate(muCar,indx,indx)

endif else begin ; No PSF convolution: just compute values

    mu = R*0
    for j=0,n_elements(R)-1 do begin
        wm2Pol = sphani_weighted_sigma2(R[j],s)
        mgePol = total( s.surf_l * exp(-0.5d*(R[j]/s.sig_l)^2) )
        mu[j] = sqrt(wm2Pol/mgePol)
    endfor

endelse

return, mu
end
;----------------------------------------------------------------------
pro jam_spherical_rms, $
    surf_lum, sigma_lum, surf_pot, sigma_pot, mbh, distance, rad, sigp, $
    BETA=beta, NORMPSF=normPsf, PIXSIZE=pixSize1, SIGMAPSF=sigmaPsf1, $
    STEP=step, NRAD=nrad
compile_opt idl2

if n_elements(beta) eq 0 then beta = 0d
if n_elements(sigmaPsf1) eq 0 then sigmaPsf1 = 0d
if n_elements(normPsf) eq 0 then normPsf = 1d
if n_elements(pixSize1) eq 0 then pixSize1 = 0d
if n_elements(mbh) eq 0 then mbh = 0d
if n_elements(step) eq 0 then step = 0
if n_elements(nrad) eq 0 then nrad = 50

pc = distance*!dpi/0.648d ; Constant factor to convert arcsec --> pc

sigmaPsf = sigmaPsf1*pc
pixSize = pixSize1*pc
step_pc = step*pc

sigma_lum_pc = sigma_lum*pc     ; Convert from arcsec to pc
lum = surf_lum*(sigma_lum_pc*sqrt(2d*!dpi))^2

sigma_pot_pc = sigma_pot*pc     ; Convert from arcsec to pc
mass = surf_pot*(sigma_pot_pc*sqrt(2d*!dpi))^2

s = {SIG_L:sigma_lum_pc, SIG_M:sigma_pot_pc, LUM:lum, MASS:mass, MBH:mbh, $
    RMIN:0d, SURF_L:surf_lum, SIGMAPSF:sigmapsf, NORMPSF:normpsf, $
    PIXSIZE:pixSize, BETA:beta, STEP:step_pc, NRAD:nrad}

sigp = shpani_second_moment(rad*pc, s)

END
;--------------------------------------------------------------------
pro test_jam_spherical_rms
;
; This example takes 1.1s on a 2GHz computer

; Realistic MGE galaxy surface brightness
;
surf_pc = [20818.4, 6762.15, 26875.8, 8136.50, 4321.21, $
           3485.05, 1164.02, 305.010, 61.6719] ; Lsun/pc^2
sigma_arcsec = [0.0172900, 0.0636382, 0.299341, 0.919666, 2.08990, $
                4.30922, 12.4481, 32.2990, 124.410]

distance = 24d    ; Mpc
mbh = 2d8 ; Msun
rad = 10d^range(-1,1,20)

; Assume self-consistency: same MGE for luminosity and potential
;
surf_lum = surf_pc
surf_pot = surf_pc
sigma_lum = sigma_arcsec
sigma_pot = sigma_arcsec

t = systime(1)

; Plot V_RMS profiles assuming isotropy (beta=0) for the Gaussians
; with sigma > 4 and three different anisotropies for the smaller
; Gaussians: (i) isotropy (2) radial (3) tangential

loadct, 12
jam_spherical_rms, surf_lum, sigma_lum, surf_pot, sigma_pot, $
    mbh, distance, rad, sigp, BETA=0, SIGMAPSF=0.5, PIXSIZE=0.3
plot, rad, sigp, PSYM=-4, /YNOZERO, YRANGE=[115,135], $
    XTITLE='R (arcsec)', YTITLE='!4r!3 (km/s)'

beta = sigma_arcsec*0 ; Initializes beta to zero
beta[where(sigma_arcsec lt 4)] = 0.3
jam_spherical_rms, surf_lum, sigma_lum, surf_pot, sigma_pot, $
    mbh, distance, rad, sigp, BETA=beta, SIGMAPSF=0.5, PIXSIZE=0.3
oplot, rad, sigp, PSYM=-4, COLOR=100

beta[where(sigma_arcsec lt 4)] = -0.3
jam_spherical_rms, surf_lum, sigma_lum, surf_pot, sigma_pot, $
    mbh, distance, rad, sigp, BETA=beta, SIGMAPSF=0.5, PIXSIZE=0.3
oplot, rad, sigp, PSYM=-4, COLOR=200

print, 'Time (s):', systime(1) - t

END
;--------------------------------------------------------------------
