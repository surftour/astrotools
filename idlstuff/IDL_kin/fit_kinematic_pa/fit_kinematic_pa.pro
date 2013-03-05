;#############################################################################
;
; Copyright (C) 2005-2008, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; Updated versions of the software are available from my web page
; http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
; If you have found this software useful for your research,
; we would appreciate an acknowledgment to use of the
; `Method described in Appendix C of Krajnovic et al. (2006)'.
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
;   FIT_KINEMATIC_PA
;
; PURPOSE:
;   Determine the global kinematic position angle of a
;   galaxy with the method described in Appendix C of
;   Krajnovic, Cappellari, de Zeeuw, & Copin 2006, MNRAS, 366, 787
;
; CALLING SEQUENCE:
;   FIT_KINEMATIC_PA, Xbin, Ybin, Vel, AngleBest, AngleError, VelSyst, $
;        DEBUG=debug, NSTEPS=nsteps, QUIET=quiet
;
; INPUT PARAMETERS:
;   XBIN, YBIN: vectors with the coordinates of the bins (or pixels)
;        measured from the centre of rotation (typically the galaxy centre).
;   VEL: measured velocity at the position (XBIN,YBIN).
;
; INPUT KEYWORDS:
;   NSTEPS: number of steps along which the angle is sampled.
;       Default is 361 steps which gives a 0.5 degr accuracy.
;       Decrease this number to limit computation time during testing.
;
; OUTPUT PARAMETER:
;   ANGLEBEST: kinematical PA. Note that this is the angle along which
;       |Vel| is maximum (note modulus!). If one reverses the sense of
;       rotation in a galaxy ANGLEBEST does not change. The sense of
;       rotation can be trivially determined by looking at the map of Vel.
;   ANGLEERROR: corresponding error to assign to ANGLEBEST.
;   VELSYST: best-fitting systemic velocity of the galaxy.
;
; REQUIRED ROUTINES:
;   The following five additional routines are needed:
;   - 1. SYMMETRIZE_VELFIELD and 2. RANGE: by Michele Cappellari
;     (included in this FIT_KINEMATIC_PA distribution)
;   - 3. SAURON_COLORMAP and 4. PLOT_VELFIELD: can be obtained from:
;     http://www.strw.leidenuniv.nl/~mcappell/idl/#binning
;   - 5. SIGRANGE: from IDL astro library http://idlastro.gsfc.nasa.gov/
;
; MODIFICATION HISTORY:
;   V1.0 -- Created by Michele Cappellari, Leiden, 30 May 2005
;   V1.1 -- Written documentation. MC, Oxford, 9 October 2007
;   V1.11 -- Corrected handling of NSTEPS keyword. Thanks to Roland Jesseit.
;       MC, Oxford, 17 October 2007
;   V1.12 -- Force error to be smaller than 1/2 of the angular step size.
;       MC, Oxford, 19 October 2007
;   V1.13 -- Determine plotting ranges from velSym instead of vel.
;       Thanks to Anne-Marie Weijmans. Leiden, 25 May 2008
;-
;----------------------------------------------------------------------------
pro fit_kinematic_pa, x, y, vel, angBest, angErr, vSyst, DEBUG=debug, NSTEPS=nsteps, QUIET=quiet
compile_opt idl2
on_error, 2

sauron_colormap

dvel = vel*0+10.0 ; Adopt here constant 10 km/s errors!

nbins = n_elements(x)
if n_elements(nsteps) eq 0 then n = 361 else n = nsteps
ang = range(0,180,n) ; 0.5 degrees steps
chi2 = ang*0
TRIANGULATE, x, y, tr
for j=0,n-1 do begin
    symmetrize_velfield, x, y, vel, velSym, SYM=1, PA=ang[j], TRIANG=tr
    chi2[j] = total(((vel-velSym)/dvel)^2)
    if keyword_set(debug) then begin
        print, 'Ang, chi2/DOF:', ang[j], chi2[j]/nbins
        plot_velfield, x, y, velSym, TITLE=chi2[j]
    endif
endfor
tmp = min(chi2,k)
angBest = ang[k]

; Compute fit at the best position
;
symmetrize_velfield, x, y, vel, velSym, SYM=1, PA=angBest
if angBest lt 0 then angBest += 180

f = where(chi2 - chi2[k] le 15.1) ; 99.9% confidence limit, to allow for systematics
AngErr = (max(ang[f]) - min(ang[f]))/2.0
if AngErr ge 45 then begin
    good = atan(tan(ang[f]/!radeg))*!radeg
    AngErr = (max(good) - min(good))/2.0
endif

angErr = angErr > (ang[1]-ang[0])/2.0 > 0.5 ; Force errors to be larger than 0.5 deg
vSyst = median(vel - velSym)

if not keyword_set(quiet) then begin
    print, '  Kin PA:', angBest, ' +/- ', angErr, FORMAT='(a,f6.1,a,f6.1)'
    print, 'Systemic Velocity:', vSyst
endif

; Plot results
;
mx = max(sigrange(velSym,FRAC=0.95),MIN=mn)
mx = mx < (-mn)
!p.multi=[0,2,1]
plot_velfield, x, y, velSym, range=[-mx,mx], TITLE='Symmetrized'
plot_velfield, x, y, vel-vSyst, range=[-mx,mx], TITLE='Data'
!p.multi=0

end
;----------------------------------------------------------------------
