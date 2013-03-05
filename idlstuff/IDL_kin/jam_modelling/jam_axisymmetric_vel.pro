;###############################################################################
;
; Copyright (C) 2008, Michele Cappellari
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
;###############################################################################
;+
; NAME:
;   JAM_AXISYMMETRIC_VEL
;
; PURPOSE:
;    This procedure calculates a prediction for the projected mean velocity V
;    for an anisotropic axisymmetric galaxy model. It implements the solution
;    of the anisotropic Jeans equations presented in equation (38) of
;    Cappellari (2008, MNRAS, 390, 71). PSF convolution is done as described
;    in the Appendix of that paper.
;
; CALLING SEQUENCE:
;    JAM_AXISYMMETRIC_VEL, $
;        surf_lum, sigma_lum, qObs_lum, surf_pot, sigma_pot, qObs_pot, $
;        inc_deg, mbh, distance, xbin, ybin, velModel, $
;        BETA=beta, CHI2=chi2, EVEL=evel, FLUX=flux, GAMMA=gamma, $
;        GOODBINS=goodBins, KAPPA=kappa, NORMPSF=normPsf, NANG=nang, NRAD=nrad, $
;        PIXANG=pixAng, /PLOT, /QUIET, SIGMAPSF=sigmaPsf, STEP=step, VEL=vel
;
; INPUT PARAMETERS:
;   SURF_LUM: vector of length N containing the peak surface brightness of the
;       MGE Gaussians describing the galaxy surface brightness in units of
;       Lsun/pc^2 (solar luminosities per parsec^2).
;   SIGMA_LUM: vector of length N containing the dispersion in arcseconds of
;       the MGE Gaussians describing the galaxy surface brightness.
;   QOBS_LUM: vector of length N containing the observed axial ratio of the MGE
;       Gaussians describing the galaxy surface brightness.
;   SURF_POT: vector of length M containing the peak value of the MGE Gaussians
;       describing the galaxy surface density in units of Msun/pc^2 (solar
;       masses per parsec^2). This is the MGE model from which the model
;       potential is computed.
;     - In a common usage scenario, with a self-consistent model, one will use
;       the same Gaussians for both the surface brightness and the potential.
;       One will use JAM_AXISYMMETRIC_RMS to fit the global M/L and then call
;       JAM_AXISYMMETRIC_VEL with SURF_POT = SURF_LUM*(M/L),
;       SIGMA_POT = SIGMA_LUM and QOBS_POT = QOBS_LUM.
;   SIGMA_POT: vector of length M containing the dispersion in arcseconds of
;       the MGE Gaussians describing the galaxy surface density.
;   QOBS_POT: vector of length M containing the observed axial ratio of the MGE
;       Gaussians describing the galaxy surface density.
;   INC_DEG: inclination in degrees (90 being edge-on).
;   MBH: Mass of a nuclear supermassive black hole in solar masses.
;   DISTANCE: distance in Mpc.
;   XBIN: Vector of length P with the X coordinates in arcseconds of the bins
;       (or pixels) at which one wants to compute the model predictions. The
;       X-axis is assumed to coincide with the galaxy projected major axis. The
;       galaxy center is at (0,0).
;   YBIN: Vector of length P with the Y coordinates in arcseconds of the bins
;       (or pixels) at which one wants to compute the model predictions. The
;       Y-axis is assumed to concide with the projected galaxy symmetry axis.
;
; KEYWORDS:
;   BETA: Vector of length N with the vertical anisotropy
;       beta_z = 1 - (sigma_z/sigma_R)^2 of the individual MGE Gaussians.
;       A scalar can be used if the model has constant anisotropy.
;       (Default: BETA=0)
;   CHI2: Reduced chi^2 describing the quality of the fit
;       chi^2 = total( ((vel[goodBins]-velModel[goodBins])/evel[goodBins])^2 )
;             / n_elements(goodBins)
;     - IMPORTANT: by default the rotation parameter KAPPA is determined by
;       requiring the projected angular momentum to be the same in the data and
;       in the model. The chi^2 is not necessarily minimized (Sec.4.2 of
;       Cappellari 2008 for details).
;   EVEL: Vector of length P with the 1sigma error associated to the VEL
;       measurements in each bin. (Default EVEL=1)
;   FLUX: In output this contains a vector of length P with the unconvolved MGE
;       surface brightness of each bin, used to plot the isophotes on the model
;       results.
;   GAMMA: Vector of length N with the tangential anisotropy
;       gamma = 1 - (sigma_phi/sigma_R)^2 of the individual MGE Gaussians.
;       A scalar can be used if the model has constant anisotropy.
;       (Default: GAMMA=0)
;     - IMPORTANT: In a common usage this keyword is *not* needed and the
;       tangential anisotropy is only parameterized with the much more
;       efficient KAPPA parameter below. However this form, which uses the
;       definition of tangential anisotropy of equation (33) of Cappellari
;       (2008) can be sometimes useful as an alternative (Sec.3.1.5 of
;       Cappellari 2008 for details).
;     - When this keyword is set and nonzero, one should generally set KAPPA=1
;       to have a model with a well defined shape of the velocity ellipsoid.
;       However KAPPA can still be fitted to allow for further generality.
;   GOODBINS: Vector of length <= P with the indices of the bins which have to
;       be included in the fit (if requested) and chi^2 calculation.
;       (Default: fit all bins).
;   KAPPA: Rotation parameter (Default: KAPPA=1). Vector of length N defining
;       by how much the model rotation of each individual MGE Gaussian is
;       scaled, with respect to the velocity of a model with an oblate (GAMMA=0)
;       velocity ellipsoid (Sec.3.1.5 of Cappellari 2008 for details).
;     - If the observed velocities are passed via the keyword VEL then a single
;		KAPPA is fitted to the data and returned in output.
;     - If GAMMA is set and nonzero then KAPPA defines by how much the model
;       velocity field is scaled with respect to a model with tangential
;       anisotropy GAMMA.
;   NORMPSF: Vector of length Q with the fraction of the total PSF flux
;       contained in the circular Gaussians describing the PSF of the
;       observations. It has to be total(NORMPSF) = 1. The PSF will be used for
;       seeing convolution of the model kinematics.
;   NRAD: Number of logarithmically spaced radial positions for which the
;       models is evaluated before interpolation and PSF convolution. One may
;       want to increase this value if the model has to be evaluated over many
;       orders of magnitudes in radius (default: NRAD=50). The computation time
;       scales as NRAD*NANG.
;   NANG: Same as for NRAD, but for the number of angular intervals
;       (default: NANG=10).
;   PIXANG: angle between the observed spaxels and the galaxy major axis X.
;   PIXSIZE: Size in arcseconds of the (square) spatial elements at which the
;       kinematics is obtained. This may correspond to the side of the spaxel
;       or lenslets of an integral-field spectrograph. This size is used to
;       compute the kernel for the seeing and aperture convolution.
;     - If this is not set, or PIXSIZE=0, then convolution is not performed.
;   /PLOT: Set this keyword to produce a plot at the end of the calculation.
;   /QUIET: Set this keyword not to print values on the screen.
;   SIGMAPSF: Vector of length Q with the dispersion in arcseconds of the
;       circular Gaussians describing the PSF of the observations.
;     - If this is not set, or SIGMAPSF=0, then convolution is not performed.
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
;   VEL: Vector of length P with the input observed mean stellar velocity V at
;        the coordinates positions given by the vectors XBIN and YBIN.
;
; OUTPUT PARAMETER:
;   VELMODEL: Vector of length P with the model predictions for the mean
;       velocity of each bin.
;
; USAGE EXAMPLE:
;    A simple usage example is given in the procedure TEST_JAM_AXISYMMETRIC_VEL
;    at the end of this file.
;
; REQUIRED ROUTINES:
;       By M. Cappellari (included in the JAM distribution):
;       - ANY
;       - DIFF
;       - MESHGRID
;       - QUADVA
;       - RANGE
;       - SIGN
;       - ROTATE_POINTS
;       - SYMMETRIZE_VELFIELD
;
;       By M. Cappellari available here
;       http://www-astro.physics.ox.ac.uk/~mxc/idl/#jam
;       - PLOT_VELFIELD
;       - SAURON_COLORMAP
;
; MODIFICATION HISTORY:
; V1.0: Written and tested by Michele Cappellari, Oxford, 24 April 2008
; V1.1: First released version. MC, Oxford, 12 August 2008
; V1.2: Updated documentation. MC, Oxford, 14 August 2008
; V2.0: Implemented PSF convolution using interpolation on polar grid. Dramatic
;     speed-up of calculation. Further documentation.
;     MC, Oxford, 11 September 2008;
; V2.01: Bug fix: when EVEL was not given, the default was not properly set.
;     Included keyword STEP. The keyword FLUX is now only used for output: the
;     surface brightness for plotting is computed from the MGE model.
;     MC, Windhoek, 29 September 2008
; V2.02: Bug fix: Velocity was not multiplied by sign(x) without convolution.
;     MC, Oxford, 22 October 2008
; V2.03: Added keywords NRAD and NANG. Thanks to Michael Williams for reporting
;     possible problems with too coarse interpolation.
;     MC, Oxford, 21 November 2008
; V2.1: Allow for alternative definition of tangential anisotropy via the
;     keyword GAMMA. MC, Oxford, 6 December 2008
; -
;#############################################################################
function janis1_jeans_mge_integrand, u, $
    DENS_lum=dens_lum, SIGMA_lum=sigma_lum, Q_lum=q_lum, $
    DENS_pot=dens_pot, SIGMA_pot=sigma_pot, Q_pot=q_pot, $
    R=R, Z=z, INC=inc, BETA=beta, KAPPA=kappa, GAMMA=gamma
compile_opt idl2, hidden
;
; This function returns the inner integrand in equation (38) of Cappellari (2008)

bani = 1d/(1d - beta) ; Anisotropy ratio b = (sig_R/sig_z)^2
R2 = R^2
z2 = z^2
u2 = u^2

s2_lum = sigma_lum^2
q2_lum = q_lum^2
s2q2_lum = s2_lum*q2_lum

s2_pot = sigma_pot^2
e2_pot = 1d - q_pot^2

; The next line gives an extra term for equation (38) using the alternative
; definition of tangential anisotropy of equation (33) instead of the default
; form of equation (35). This term is zero by default (gamma=0).
;
dani = gamma*bani*s2q2_lum/R2

; Double summation over (j,k) of eq.(38) vectorized over integration variable u.
; The j-index refers to the Gaussians describing the total mass,
; from which the potential is derived, while the k-index is used
; for the MGE components describing the galaxy stellar luminosity.
;
sum = 0d
nu = kappa^2 * dens_lum * exp(-0.5d/s2_lum * (R2 + z2/q2_lum))
for j=0,n_elements(dens_pot)-1 do begin ; loop over mass Gaussians
    p2 = 1d - e2_pot[j]*u2
    hj = exp(-0.5d/s2_pot[j]*u2*(r2+z2/p2))/sqrt(p2) ; equation (17)
    e = q_pot[j]*dens_pot[j]*hj*u2
    for k=0,n_elements(dens_lum)-1 do begin ; loop over luminous Gaussians
        c = e2_pot[j] - s2q2_lum[k]/s2_pot[j] ; equation (22)
        d = 1d - bani[k]*q2_lum[k] - ((1d - bani[k])*c + e2_pot[j]*bani[k])*u2 ; equation (23)
        sum += nu[k]*e*(d + dani[k]) / (1d - c*u2) ; equation (38)
    endfor
endfor

return, sum ; sum is a vector with the same size as the input u
end
;----------------------------------------------------------------------
function janis1_jeans_mge_los_integrand, z1, $
    DENS_lum=dens_lum, SIGMA_lum=sigma_lum, Q_lum=q_lum, $
    DENS_pot=dens_pot, SIGMA_pot=sigma_pot, Q_pot=q_pot, $
    X1=x1, Y1=y1, INC=inc, BETA=beta, KAPPA=kappa, GAMMA=gamma
compile_opt idl2, hidden
;
; This function returns the outer integrand in equation (38) of Cappellari (2008)

n = n_elements(z1)
int = dblarr(n,/NOZERO)

for j=0,n-1 do begin ; loop over all z1 values along the LOS

    R = sqrt( (z1[j]*sin(inc) - y1*cos(inc))^2 + x1^2 ) ; Equation (25)
    z = sqrt( (z1[j]*cos(inc) + y1*sin(inc))^2 )

    functargs = {DENS_lum:dens_lum, SIGMA_lum:sigma_lum, Q_lum:q_lum, $
                 DENS_pot:dens_pot, SIGMA_pot:sigma_pot, Q_pot:q_pot, $
                 R:r, Z:z, INC:inc, BETA:beta, KAPPA:kappa, GAMMA:gamma}
    quadva, 'janis1_jeans_mge_integrand', [0d,1d], tmp, $
        RELTOL=1d-5, ABSTOL=0, FUNCTARGS=functargs

    nu = total(dens_lum * exp(-0.5d/sigma_lum^2 * (R^2 + (z/q_lum)^2)))
    int[j] = sqrt(nu*tmp>0) ; Prevents this from becoming nearly zero but negative due to rounding errors

endfor

; This routine returns a vector of values computed at different values of z1
;
G = 0.00430237d    ; (km/s)^2 pc/Msun [6.674e-11 SI units]

return, 2d*sqrt(!dpi*G)*x1*sin(inc) * int
end
;----------------------------------------------------------------------------
function janis1_weighted_first_moment, x, y, inc_deg, s
compile_opt idl2, hidden
;
; This routine gives the projected weighted first moment \Sigma*<V_los>

; Axisymmetric deprojection of both luminous and total mass.
; See equation (12)-(14) of Cappellari (2008)

inc = inc_deg/!radeg

ms = 3d*max(s.sigma_lum)

qintr_lum = s.qobs_lum^2 - cos(inc)^2
if any(qintr_lum le 0d) then message, 'Inclination too low q < 0'
qintr_lum = sqrt(qintr_lum)/sin(inc)
if any(qintr_lum lt 0.05d) then message, 'q < 0.05 components'
dens_lum = s.surf_lum*s.qobs_lum / (s.sigma_lum*qintr_lum*sqrt(2d*!dpi))

qintr_pot = s.qobs_pot^2 - cos(inc)^2
if any(qintr_pot le 0d) then message, 'Inclination too low q < 0'
qintr_pot = sqrt(qintr_pot)/sin(inc)
if any(qintr_pot lt 0.05d) then message, 'q < 0.05 components'
dens_pot = s.surf_pot*s.qobs_pot / (s.sigma_pot*qintr_pot*sqrt(2d*!dpi))

functargs = {DENS_lum:dens_lum, SIGMA_lum:s.sigma_lum, Q_lum:qintr_lum, $
             DENS_pot:dens_pot, SIGMA_pot:s.sigma_pot, Q_pot:qintr_pot, $
             X1:x, Y1:y, INC:inc, BETA:s.beta, KAPPA:s.kappa, GAMMA:s.gamma}
quadva, 'janis1_jeans_mge_los_integrand', [-ms,ms], sb_mu1, $
    RELTOL=1d-5, ABSTOL=0d, FUNCTARGS=functargs

return, sb_mu1
end
;----------------------------------------------------------------------
function janis1_first_moment, x, y, inc_deg, s
compile_opt idl2, hidden
;
; This routine gives the first V moment after convolution with a PSF.
; The convolution is done using interpolation of the model on a
; polar grid, as described in Appendix A of Cappellari (2008).

; Define parameters of polar grid for interpolation
;
w = where(s.sigma_lum lt max(x),m) ; Characteristic MGE axial ratio in observed range
if m lt 3 then qmed = median(s.qobs_lum) else qmed = median(s.qobs_lum[w])
rell = sqrt(x^2 + (y/qmed)^2) ; Elliptical radius of input (x,y)

; Kernel step is 1/4 of largest value between sigma(min) and 1/2 pixel size.
; Kernel half size is the sum of 3*sigma(max) and 1/2 pixel diagonal.
;
if max(s.sigmaPsf) gt 0 and s.pixSize gt 0 then begin ; PSF convolution
    if s.step gt 0 then step = s.step $
        else step = (s.pixSize/2d > min(s.sigmaPsf))/4d
    mx = 3d*max(s.sigmaPsf) + s.pixSize/sqrt(2d)
endif else begin ; No convolution
    step = min(rell)
    mx = 0d
endelse

; Make linear grid in log of elliptical radius RAD and eccentric anomaly ANG
; See Appendix A
;
rmax = max(rell) + mx ; Major axis of ellipse containing all data + convolution
nRad = s.NRAD ; radial elements
nAng = s.NANG ; angular sectors
logRad = range(alog(step),alog(rmax),nRad) ; Linear grid in log(rell)
ang = range(0d,!dpi/2d,nAng) ; Linear grid in eccentric anomaly
meshgrid, exp(logRad), ang, radGrid, angGrid
xPol = radGrid*cos(angGrid)
yPol = radGrid*sin(angGrid) * qmed

; The model Vrms computation is only performed on the polar grid
; which is then used to interpolate the values at any other location
;
wm2Pol = xPol*0
mgePol = wm2Pol
for j=0,n_elements(xPol)-1 do begin
    wm2Pol[j] = janis1_weighted_first_moment(xPol[j], yPol[j], inc_deg, s)
    mgePol[j] = total(s.surf_lum * exp(-0.5d/s.sigma_lum^2 * (xPol[j]^2 + (yPol[j]/s.qobs_lum)^2)))
endfor

if max(s.sigmaPsf) gt 0 and s.pixSize gt 0 then begin ; PSF convolution

    nx = ceil(rmax/step/64d)*64       ; Make 2*nx dimensions divisible by 2^7
    ny = ceil(rmax*qmed/step/64d)*64  ; for much faster FFT computation
    x1 = range(-nx,nx,2*nx)*step
    y1 = range(-ny,ny,2*ny)*step
    meshgrid, x1, y1, xCar, yCar ; Cartesian grid for convolution

    ; Interpolate MGE model and V over cartesian grid
    ;
    r1 = 0.5d*alog(xCar^2 + (yCar/qmed)^2) ; Log elliptical radius of cartesian grid
    e1 = atan(abs(yCar/qmed),abs(xCar))    ; Eccentric anomaly of cartesian grid
    indx = (r1 - logRad[0])/(logRad[1] - logRad[0]) ; Convert to indices for INTERPOLATE
    indy = (e1 - ang[0])/(ang[1] - ang[0])
    wm2Car = interpolate(wm2Pol,indx,indy,CUBIC=-0.5)*sign(xCar) ; V is anti-symmetric
    mgeCar = interpolate(mgePol,indx,indy,CUBIC=-0.5)

    nk = ceil(mx/step)
    kgrid = range(-nk,nk,2*nk)*step
    meshgrid, kgrid, kgrid, xx, yy ; Kernel is square
    rotate_points, xx, yy, s.pixAng, xgrid, ygrid

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

    ; Seeing and aperture convolution with equation (A2)
    ;
    muCar = convolve(wm2Car,kernel)/convolve(mgeCar,kernel)

    ; Interpolate convolved image at observed apertures.
    ; Aperture integration was already included in the kernel.
    ;
    indx = (x - x1[0])/(x1[1] - x1[0]) ; Convert to indices for INTERPOLATE
    indy = (y - y1[0])/(y1[1] - y1[0])
    mu = interpolate(muCar,indx,indy)

endif else begin ; No PSF convolution: just interpolate values

    muPol = wm2Pol/mgePol
    r1 = 0.5d*alog(x^2 + (y/qmed)^2) ; Log elliptical radius of input (x,y)
    e1 = atan(abs(y/qmed),abs(x))    ; Eccentric anomaly of input (x,y)
    indx = (r1 - logRad[0])/(logRad[1] - logRad[0]) ; Convert to indices for INTERPOLATE
    indy = (e1 - ang[0])/(ang[1] - ang[0])
    mu = interpolate(muPol,indx,indy,CUBIC=-0.5)*sign(x)

endelse

return, mu
end
;----------------------------------------------------------------------
pro jam_axisymmetric_vel, $
    surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot1, $
    inc, mbh, distance, xbin, ybin, velModel, $
    NORMPSF=normPsf, PIXANG=pixAng, PIXSIZE=pixSize, PLOT=plot, $
    VEL=vel, EVEL=evel, CHI2=chi2, SIGMAPSF=sigmaPsf, GOODBINS=goodBins, $
    QUIET=quiet, FLUX=flux, BETA=beta1, KAPPA=kappa1, GAMMA=gamma1, $
    STEP=step, NRAD=nrad, NANG=nang
compile_opt idl2
on_error, 2

ngau = n_elements(surf_lum)
if n_elements(beta1) eq 0 then beta1 = 0d   ; Anisotropy parameter beta = 1 - (sig_z/sig_R)^2 (beta=0-->circle)
if n_elements(beta1) eq 1 then beta = sigma_lum*0+beta1 ; All components have the same anisotropy
if n_elements(beta1) eq ngau then beta = beta1

if n_elements(kappa1) eq 0 then kappa1 = 1d ; Anisotropy parameter V_obs = kappa * V_iso (kappa=1-->circle)
if n_elements(kappa1) eq 1 then kappa = sigma_lum*0+kappa1 ; All components have the same anisotropy
if n_elements(kappa1) eq ngau then kappa = kappa1

if n_elements(gamma1) eq 0 then gamma = 0d  ; Anisotropy parameter gamma = 1 - (sig_phi/sig_R)^2 (gamma=0-->circle)
if n_elements(gamma1) eq 1 then gamma = sigma_lum*0+gamma1 ; All components have the same anisotropy
if n_elements(gamma1) eq ngau then gamma = gamma1

if n_elements(evel) eq 0 and n_elements(vel) gt 0 then evel = vel*0+1 ; Constant errors
if n_elements(sigmaPsf) eq 0 then sigmaPsf = 0d
if n_elements(normPsf) eq 0 then normPsf = 1d
if n_elements(pixAng) eq 0 then pixAng = 0d
if n_elements(pixSize) eq 0 then pixSize = 0d
if n_elements(step) eq 0 then step = 0
if n_elements(nrad) eq 0 then nrad = 50
if n_elements(nang) eq 0 then nang = 10

t = systime(1)

pc = distance*!dpi/0.648d ; Constant factor to convert arcsec --> pc

surf_lum_pc = surf_lum
surf_pot_pc = surf_pot
sigma_lum_pc = sigma_lum*pc         ; Convert from arcsec to pc
sigma_pot_pc = sigma_pot*pc         ; Convert from arcsec to pc
xbin_pc = xbin*pc                   ; Convert all distances to pc
ybin_pc = ybin*pc
pixSize_pc = pixSize*pc
sigmaPsf_pc = sigmaPsf*pc
step_pc = step*pc

; Add a Gaussian with small sigma and the same total mass as the BH.
; The Gaussian provides an excellent representation of the second moments
; of a point-like mass, to 1% accuracy out to a radius 2*sigmaBH.
; The error increses to 14% at 1*sigmaBH; Independently of the BH mass.
;
if mbh gt 0 then begin
    sigmaBH_pc = 0.01d*pc ; Adopt for the BH just a very small size sigma=0.01"
    surfBH_pc = mbh/(2d*!dpi*sigmaBH_pc^2)
    surf_pot_pc = [surfBH_pc,surf_pot_pc[*]] ; Add Gaussian to potential only!
    sigma_pot_pc = [sigmaBH_pc,sigma_pot_pc[*]]
    qobs_pot = [1d,qobs_pot1[*]]  ; Make sure vectors do not have extra dimensions
endif else qobs_pot = qobs_pot1 ; Important: do not change input vector qobs_pot!

s = {SURF_LUM:surf_lum_pc, SIGMA_LUM:sigma_lum_pc, QOBS_LUM:qobs_lum, $
     SURF_POT:surf_pot_pc, SIGMA_POT:sigma_pot_pc, QOBS_POT:qobs_pot, $
     SIGMAPSF:sigmaPsf_pc, NORMPSF:normPsf, BETA:beta, KAPPA:kappa, GAMMA:gamma, $
     PIXSIZE:pixSize_pc, PIXANG:pixAng, STEP:step_pc, NRAD:nrad, NANG:nang}

velModel = janis1_first_moment(xbin_pc, ybin_pc, inc, s)

if ~keyword_set(quiet) then print, 'Elapsed time sec:', systime(1) - t

;###### Output and optional KAPPA anisotropy fit
; If VEL keyword is not given all this section is skipped

if n_elements(vel) gt 0 then begin

    ; Only consider the good bins for the chi^2 estimation
    ;
    if n_elements(goodBins) eq 0 then goodBins = indgen(n_elements(xbin))

    ; Scale by having the same angular momentum
    ; in the model and in the galaxy (equation 52)
    ;
    kappa1 = total(abs(vel[goodBins]*xbin[goodBins])) $
           / total(abs(velModel[goodBins]*xbin[goodBins]))

    ; Measure the scaling one would have from a standard chi^2 fit of the V field.
    ; This value is only used to get proper sense of rotation for the model.
    ; y1 = rms; dy1 = erms (y1 are the data, y2 the model)
    ; scale = total(y1*y2/dy1^2)/total(y2^2/dy1^2)  (equation 51)
    ;
    kappa2 = total(vel[goodBins]*velModel[goodBins]/evel[goodBins]^2) $
           / total(velModel[goodBins]^2/evel[goodBins]^2)

    velModel *= kappa1*sign(kappa2)

    chi2 = total( ((vel[goodBins]-velModel[goodBins])/evel[goodBins])^2 ) / n_elements(goodBins)
    if ~keyword_set(quiet) then print, inc, chi2, kappa1, beta[0], mbh, $
        FORMAT='("inc: ", f0.1, "; chi^2: ", g0.3, "; kappa: ", f0.2, "; beta: ", f0.2, "; M_BH: ", e0.1)'

    if keyword_set(plot) then begin
        str = 'i=' + STRING(inc,FORMAT='(f0.1)') + $
              ' !7b!X!Dz!N=' + STRING(beta[0],FORMAT='(f0.2)') + $
              ' !7j!X=' + STRING(abs(kappa1),FORMAT='(f0.2)') + $
              ' BH=' + STRING(mbh,FORMAT='(e0.1)')

        ; The X-axis was aligned with the major axis as input
        ;
        dx = randomn(seed,n_elements(xbin))*0.001 ; Avoids collinear-points IDL bug in TRIGRID
        symmetrize_velfield, xbin[goodBins]+dx, ybin[goodBins]+dx, vel[goodBins], tmp, SYM=1, PA=90d
        vel1 = vel
        vel1[goodBins] = tmp ; Only symmetrize good bins

        flux = 0d ; Total MGE surface brightness for plotting
        for j=0,n_elements(surf_lum)-1 do $
            flux += surf_lum[j]*exp(-0.5d/sigma_lum[j]^2*(xbin^2 + (ybin/qObs_lum[j])^2))

        !p.multi = [0,1,2]
        sauron_colormap
        mx = max(sigrange(vel1[goodBins]),MIN=mn)
        mx = abs(mn) < mx
        range = [-mx,mx]
        plot_velfield, xbin, ybin, vel1, RANGE=range, FLUX=flux

        ; Overplot bad bins on the data
        ;
        tmp = xbin*0+1
        tmp[goodBins] = 0
        badBins = where(tmp,m)
        if m gt 0 then begin
            oplot, xbin[badBins], ybin[badBins], PSYM=4, SYMSIZE=0.5, COLOR=0, THICK=1
            oplot, xbin[badBins], ybin[badBins], PSYM=4, SYMSIZE=0.25, COLOR=255, THICK=1
        endif
        plot_velfield, xbin, ybin, velModel, RANGE=range, TITLE=str, FLUX=flux
        !p.multi = 0
    endif

endif

end
;----------------------------------------------------------------------
pro test_jam_axisymmetric_vel
;
; This example takes 300s on a 2GHz computer

; Realistic MGE galaxy surface brightness
;
surf = [3146.61, 1110.0924, 315.25342, 257.74897, 179.70656, 51.292491,$
    21.554115, 39.289108, 8.5841701, 5.4990117, 1.7729374, 0.10307876] ; Lsun/pc^2
sigma = [0.0273000, 0.108019, 0.278020, 0.623963, 1.32313, 3.63162,$
    4.88739, 5.96427, 14.8843, 26.8627, 53.0694, 82.2861] ; arcsec
qObs = [0.751920, 0.728844, 0.774995, 0.800000, 0.729212, 0.522413,$
    0.792778, 0.255372, 0.547606, 0.250000, 0.300203, 0.735538]

; Assume self-consistency: same MGE for both light and mass
;
surf_lum = surf
sigma_lum = sigma
qobs_lum = qObs
surf_pot = surf
sigma_pot = sigma
qobs_pot = qObs

inc = 80d ; Inclination in degrees
mbh = 1d6 ; BH mass in solar masses
distance = 10d ; Mpc

nbins = 1000
xbin = 10*(randomu(seed,nbins)-0.5) ; Random between [-5,5] arcsec
ybin = 10*(randomu(seed,nbins)-0.5) ; Random between [-5,5] arcsec
vel = 4*xbin ; Adopt nonsense input rotation field

t = systime(1)

; Note that the model in this example does not look in any way like the data!

jam_axisymmetric_vel, $
    surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, $
    inc, mbh, distance, xbin, ybin, velModel, $
    BETA=0.1, VEL=vel, /PLOT, SIGMAPSF=0.5, PIXSIZE=0.3, NRAD=20

print, 'Time (s):', systime(1) - t

end
;----------------------------------------------------------------------
