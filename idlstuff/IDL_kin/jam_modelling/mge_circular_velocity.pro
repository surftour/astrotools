;#############################################################################
;
; Copyright (C) 2003-2008, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; Updated versions of the software are available from my web page
; http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
; If you have found this software useful for your research,
; I would appreciate an acknowledgment to the use of the
; `JAM modelling package of Cappellari (2008)'.
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
;   MGE_CIRCULAR_VELOCITY
;
; PURPOSE:
;    This procedure calculates the circular velocity in the equatorial plane of
;    an axisymmetric galaxy model described by a Multi-Gaussian Expansion
;    parametrization. This implementation follows the approach described in
;    Appendix A of Cappellari et al. (2002, ApJ, 578, 787), which allows for
;    quick and accurate calculations also at very small
;    and very large radii.
;
; CALLING SEQUENCE:
;    MGE_CIRCULAR_VELOCITY, surf_pot, sigma_pot, qObs_pot, $
;        inc_deg, mbh, distance, rad, vcirc, SOFT=soft
;
; INPUT PARAMETERS:
;   SURF_POT: vector of length M containing the peak value of the MGE Gaussians
;       describing the galaxy surface density in units of Msun/pc^2 (solar
;       masses per parsec^2). This is the MGE model from which the model
;       potential is computed.
;   SIGMA_POT: vector of length M containing the dispersion in arcseconds of
;       the MGE Gaussians describing the galaxy surface density.
;   QOBS_POT: vector of length M containing the observed axial ratio of the MGE
;       Gaussians describing the galaxy surface density.
;   INC_DEG: inclination in degrees (90 being edge-on).
;   MBH: Mass of a nuclear supermassive black hole in solar masses.
;   DISTANCE: distance of the galaxy in Mpc.
;   RAD: Vector of length P with the radius in arcseconds, measured from the
;       galaxy centre, at which one wants to compute the model predictions.
;
; KEYWORDS:
;   SOFT: Softening length in arcsec for the Keplerian potential of the black
;       hole. When this keyword is nonzero the black hole potential will be
;       replaced by a Plummer potential with the given scale length.
;
; OUTPUT PARAMETER:
;   VCIRC: Vector of length P with the model predictions for the circular
;       velocity at the given input radii RAD.
;
; USAGE EXAMPLE:
;    A simple usage example is given in the procedure
;    TEST_MGE_CIRCULAR_VELOCITY at the end of this file.
;
; REQUIRED ROUTINES:
;       By M. Cappellari (included in the JAM distribution):
;       - ANY
;       - DIFF
;       - QUADVA
;       - RANGE
;
; MODIFICATION HISTORY:
; V1.0: Written and tested as part of the implementation of the orbit-based
;     numerical superposition Schwarzschild's method described in
;     Cappellari et al. (2006). Michele Cappellari, Leiden, 3 February 2003
; V3.0: This version retains only the few routines required for the computation
;     of the circular velocity. All other unnecessary modelling routines have
;     been removed. MC, Leiden, 22 November 2005
; V3.01: Minor code polishing. MC, Oxford, 9 November 2006
; V3.02: First released version. Included documentation.
;     MC, Windhoek, 1 October 2008
;-
;#############################################################################
;
; The following set of routines computes the R acceleration
; for a density parametrized via the Multi-Gaussian Expansion method.
; The routines are designed to GUARANTEE a maximum relative error of
; 1e-4 in the case of positive Gaussians. This maximum error is reached
; only at the extremes of the usable radial range and only for a very
; flattened Gaussian (q=0.1). Inside the radial range normally adopted
; during orbit integration the error is instead <1e-6.
;
function mgevc_accelerationR_dRRcapitalh, u, R2=r2, Z2=z2, E2=e2, S2=s2
compile_opt idl2, hidden
;
; Computes: -D[H[R,z,u],R]/R

    u2 = u^2
    p2 = 1d - e2*u2
    us2 = u2/s2
    return, exp(-0.5d*us2*(r2+z2/p2))*us2/sqrt(p2) ; Cfr. equation (A3)
end
;----------------------------------------------------------------------
function mgevc_accR, r, z, dens, sigma, qintr, bhMass, soft
compile_opt idl2, hidden

nrad = n_elements(R)
ngauss = n_elements(dens)
mgepot = dblarr(nrad,/NOZERO)
pot = dblarr(ngauss,/NOZERO)
e2 = 1d - qintr^2
s2 = sigma^2
r2 = R^2
z2 = z^2
d2 = r2 + z2

for k=0,nrad-1 do begin
    for j=0,ngauss-1 do begin
        if (d2[k] lt s2[j]/240d^2) then begin
            e = sqrt(e2[j]) ; pot is Integral in {u,0,1} of -D[H[R,z,u],R]/R at (R,z)=0
            pot[j] = (asin(e)/e - qintr[j])/(2d*e2[j]*s2[j]) ; Cfr. equation (A5)
        endif else if (d2[k] lt s2[j]*245d^2) then begin
            functargs = {R2:r2[k],Z2:z2[k],E2:e2[j],S2:s2[j]}
            quadva, 'mgevc_accelerationR_dRRcapitalh', [0d,1d], int, $
                RELTOL=1d-6, ABSTOL=0d, FUNCTARGS=functargs
            pot[j] = int
        endif else $ ; R acceleration in Keplerian limit (Cappellari et al. 2002)
           pot[j] = sqrt(!dpi/2d)*sigma[j]/d2[k]^1.5d ; Cfr. equation (A4)
    endfor
    mgepot[k] = total(s2*qintr*dens*pot)
endfor

G = 0.00430237d    ; (km/s)^2 pc/Msun [6.674e-11 SI units]

return, -r*(4d*!dpi*G*mgepot + G*bhMass/(d2 + soft^2)^1.5d)
end
;----------------------------------------------------------------------
pro mge_circular_velocity, surf_pc, sigma_arcsec, qobs, $
    inc_deg, mbh, distance, rad, vcirc, SOFT=softl_arcsec
compile_opt idl2
on_error, 2

if n_elements(softl_arcsec) eq 0 then softl_arcsec = 0d

pc = distance*!dpi/0.648d ; Constant factor to convert arcsec --> pc

soft = softl_arcsec*pc      ; Convert from arcsec to pc
rcirc = rad*pc              ; Convert from arcsec to pc
sigma = sigma_arcsec*pc     ; Convert from arcsec to pc

; Axisymmetric deprojection of total mass.
; See equation (12)-(14) of Cappellari (2008)
;
inc = inc_deg / !radeg      ; Convert inclination to radians
qintr = qobs^2 - cos(inc)^2
if any(qintr le 0.0) then message, 'Inclination too low for deprojection'
qintr = sqrt(qintr)/sin(inc)
if any(qintr le 0.05) then message, 'q < 0.05 components'
dens = surf_pc*qobs/(qintr*sigma*sqrt(2d*!dpi)) ; MGE deprojection

; Equality of gravitational and centrifugal acceleration accR at z=0
; R Vphi^2 == accR --> R (vcirc/R)^2 == accR
;
accR = mgevc_accR(Rcirc, Rcirc*0d, dens, sigma, qintr<0.999d, mbh, soft)
vcirc = sqrt(rcirc*abs(accR))  ; circular velocity at rcirc

end
;----------------------------------------------------------------------------
pro test_mge_circular_velocity
;
; This example takes a fraction of a second on a 2GHz computer

; Realistic MGE galaxy surface brightness
;
surf = [3146.61, 1110.0924, 315.25342, 257.74897, 179.70656, 51.292491,$
    21.554115, 39.289108, 8.5841701, 5.4990117, 1.7729374, 0.10307876] ; Lsun/pc^2
sigma = [0.0273000, 0.108019, 0.278020, 0.623963, 1.32313, 3.63162,$
    4.88739, 5.96427, 14.8843, 26.8627, 53.0694, 82.2861] ; arcsec
qObs = [0.751920, 0.728844, 0.774995, 0.800000, 0.729212, 0.522413,$
    0.792778, 0.255372, 0.547606, 0.250000, 0.300203, 0.735538]

inc = 80d ; Inclination in degrees
mbh = 1d6 ; BH mass in solar masses
distance = 10d ; Mpc
rad = 10d^range(-1,1.1,20) ; Radii in arscec where Vcirc has to be computed
ml = 5.0 ; Adopted M/L ratio

mge_circular_velocity, surf*ml, sigma, qobs, inc, mbh, distance, rad, vcirc

plot, rad, vcirc, PSYM=-4, XTITLE='R (arcsec)', YTITLE='V!Dcirc!N (km/s)'

end
;----------------------------------------------------------------------
