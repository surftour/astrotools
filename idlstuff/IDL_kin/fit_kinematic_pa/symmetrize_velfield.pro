;######################################################################
;
; Copyright (C) 2004-2007, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;######################################################################
pro symm_rotate_points, x, y, ang, xNew, yNew
compile_opt idl2
;
; Michele cappellari, Leiden, 15 August 2003

theta = (ang - 90.0)/!RADEG  ; Measure from positive Y axis
xNew = x*COS(theta) - y*SIN(theta)
yNew = x*SIN(theta) + y*COS(theta)

end
;----------------------------------------------------------------------
function symm_sign, a
compile_opt idl2
;
; Michele Cappellari, Leiden, 29 May 2003

return, (a gt 0)*2 - 1
end
;----------------------------------------------------------------------
pro symmetrize_velfield, xbin, ybin, velBin, velSym, SYM=sym, PA=pa, TRIANG=tr
compile_opt idl2
;
; This routine generates a bi-symmetric ('axisymmetric') or point-symmetric
; ('triaxial') version of a given set of kinematical measurements.
; PA: is the angle in degrees, measured counter-clockwise,
;   from the vertical axis (Y axis) to the galaxy major axis.
; SYM: by-simmetry: is 1 for (V,h3,h5) and 2 for (sigma,h4,h6)
;   point-simmetry: is 3 for (V,h3,h5) and 4 for (sigma,h4,h6)
;
; HISTORY:
;
; V1.0: Michele Cappellari, Vicenza, 21 May 2004
; V1.01: Added MISSING keyword to TRIGRID call. Flipped velocity sign.
;   Written basic documentation. MC, Leiden, 25 May 2004
; V1.1: Included point-symmetric case. Remco van den Bosch, Leiden, 18 January 2005
; V1.11: Minor code revisions. MC, Leiden, 23 May 2005
; V1.12: Important: changed definition of PA to be measured counterclockwise
;   with respect to the positive Y axis, as in astronomical convention and
;   consistently with my FIND_GALAXY routine. MC, Leiden, 1 June 2005
; V1.13: Added optional keyword TRIANG. Corrected rare situation with w=-1.
;   MC, Leiden, 2 June 2005
; V1.14: Added prefix SYMM_ to internal functions to prevent conflicts
;   with external functions with the same name. MC, Oxford, 11 May 2007
;
on_error, 2

n = n_elements(xbin)
if n_elements(pa) eq 0 then pa = 0.0
symm_rotate_points, xbin, ybin, -pa, x, y  ; Negative PA for counter-clockwise
if n_elements(tr) eq 0 then TRIANGULATE, x, y, tr
velSym = velBin

for j=0,n-1 do begin
    xtmp = abs(x[j]) > 1e-4 ; Avoid degenerate case x=0
    ytmp = abs(y[j]) > 1e-4 ; Avoid degenerate case y=0
    vel = TRIGRID(x, y, velBin, tr, XOUT=[-xtmp,xtmp], YOUT=[-ytmp,ytmp], MISSING=1234567)
    w = where(vel ne 1234567,p)
    if p gt 1 then begin
        if sym ge 3 then begin
            vec = (x[j]*y[j] gt 0) ? [-1,0,0,1] : [0,1,-1,0]
            tmp = where(vec[w], m)
        endif
        case sym of
            1: velSym[j] = mean((vel*[-1,1,-1,1])[w])*symm_sign(x[j])
            2: velSym[j] = mean(vel[w])
            3: velSym[j] = total((vel*vec)[w])*symm_sign(x[j])/m
            4: velSym[j] = total((vel*abs(vec))[w])/m
        endcase
    endif
endfor

end
;----------------------------------------------------------------------
