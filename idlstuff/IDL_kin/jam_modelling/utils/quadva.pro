;#############################################################################
;
; Copyright (C) 2007-2008, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; Updated versions of the software are available from my web page
; http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
; If you have found this software useful for your research,
; I would appreciate an acknowledgment and a link to the website.
; 
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;#############################################################################
;
; function [Ifx,errbnd] = quadva(fun,interval,reltol,abstol)
; When INTERVAL = [A,B], computes the integral of a continuous function
; f(x) from A to B for A < B.  A can be -Inf and/or B can be Inf. FUN
; accepts a row vector X and returns a vector Y with Y(m) = f(X(m))
; for m = 1,...,length(X). FUN can be a string that QUADVA will convert
; to a vectorized, inline function.
;
; QUADVA returns an approximation Ifx to the integral and optionally an
; approximate bound, ERRBND, on the error |integral - Ifx|. It attempts
; to compute Ifx such that |Ifx - integral| <= max(ABSTOL,RELTOL*|Ifx|).
; If QUADVA is unsuccessful, Ifx and ERRBND are still meaningful, so a
; warning is issued that includes ERRBND.
;
; If the interval is infinite, say [a,Inf), then for the integral to
; exist, f(x) must decay as x -> Inf and QUADVA requires it to decay
; rapidly. Special methods should be used for oscillatory functions on
; infinite intervals, but QUADVA can be used if f(x) decays fast enough.
;
; QUADVA will integrate f(x) that are singular at finite end points if
; the singularities are not too strong. For example, it will integrate
; f(x) that behave like log|x-c| or |x-c|^p for p >= -1/2 with c = a
; and/or c = b. If f(x) is singular at points inside (A,B), write the
; integral as a sum of integrals over subintervals with the singular
; points as end points, compute them with QUADVA, and add the results.
;
; QUADVA starts with samples of f(x) at 150 points in (A,B). It must be
; able to recognize the behavior of f(x) from these samples, so if f(x)
; oscillates very rapidly or has sharp peaks, it may be necessary to
; subdivide the interval. To do this, make INTERVAL an array with entries
; (breakpoints) that increase from A = INTERVAL(1) to B = INTERVAL(end).
; If f(x) is only piecewise smooth, the points of discontinuity should
; be breakpoints.
; 
; Based on the algorithm "Vectorized Adaptive Quadrature in Matlab"
; L.F. Shampine, Journal of Computational and Applied Mathematics 211 (2008) 131-140
;
; MODIFICATION HISTORY:
; V1.0: Written and tested against the corresponding MATLAB version,
;       Michele Cappellari, Oxford, 22 October, 2007
; V1.1: Allow function parameters to be passed via the FUNCTARGS keyword.
;       MC, Oxford, 23 October 2007
; V1.11: Added STATUS keyword. Provide more informative error messages on failure.
;       MC, Windoek, 5 October 2008
;
;----------------------------------------------------------------------
function qva_split, interval
compile_opt idl2, hidden

; If breakpoints are specified, split subintervals in
; half as needed to get a minimum of 10 subintervals.
v = interval
while 1 do begin
    npts = n_elements(interval)
    if npts ge 11 then break
    v = dblarr(npts*2-1,/NOZERO)
    v[0:2*npts-2:2] = interval
    v[1:2*npts-2:2] = interval[0:npts-2] + 0.5d*diff(interval)
    interval = v
endwhile

return, v
end ; qva_split
;----------------------------------------------------------------------
function qva_check_spacing, x
compile_opt idl2, hidden

ax = abs(x)
n = n_elements(x)
too_close = any(diff(x) le 100d*(machar(/DOUBLE)).eps*(ax[0:n-2]>ax[1:n-1]))

return, too_close
end ; too_close
;----------------------------------------------------------------------
pro qva_f1,fun,t,a,b,y,too_close,FUNCTARGS=fargs
compile_opt idl2, hidden

; Transform to weaken singularities at both ends: [a,b] -> [-1,1]
Tt = 0.25d*(b-a)*t*(3d - t^2) + 0.5d*(b+a)
too_close = qva_check_spacing(Tt)
if ~too_close then begin
    if n_elements(fargs) GT 0 then $
        y = call_function(fun,Tt,_EXTRA=fargs) $
    else $
        y = call_function(fun,Tt)
    y = 0.75d*(b-a)*y*(1d - t^2)
endif

end ; qva_f1
;----------------------------------------------------------------------
pro qva_f2,fun,t,a,b,y,too_close,FUNCTARGS=fargs
compile_opt idl2, hidden

; Transform to weaken singularity at left end: [a,Inf) -> [0,Inf).
; Then transform to finite interval: [0,Inf) -> [0,1].
Tt = t / (1d - t)
T2t = a + Tt^2
too_close = qva_check_spacing(T2t)
if ~too_close then begin
    if n_elements(fargs) GT 0 then $
        y = call_function(fun,T2t,_EXTRA=fargs) $
    else $
        y = call_function(fun,T2t)
    y =  2d*Tt * y / (1d - t)^2
endif

end ; qva_f2
;----------------------------------------------------------------------
pro qva_f3,fun,t,a,b,y,too_close,FUNCTARGS=fargs
compile_opt idl2, hidden

; Transform to weaken singularity at right end: (-Inf,b] -> (-Inf,b].
; Then transform to finite interval: (-Inf,b] -> (-1,0].
Tt = t / (1d + t)
T2t = b - Tt^2
too_close = qva_check_spacing(T2t)
if ~too_close then begin
    if n_elements(fargs) GT 0 then $
        y = call_function(fun,T2t,_EXTRA=fargs) $
    else $
        y = call_function(fun,T2t)
    y = -2d*Tt * y / (1d + t)^2
endif

end ; qva_f3
;----------------------------------------------------------------------
pro qva_f4,fun,t,a,b,y,too_close,FUNCTARGS=fargs
compile_opt idl2, hidden

; Transform to finite interval: (-Inf,Inf) -> (-1,1).
Tt = t / (1d - t^2)
too_close = qva_check_spacing(Tt)
if ~too_close then begin
    if n_elements(fargs) GT 0 then $
        y = call_function(fun,Tt,_EXTRA=fargs) $
    else $
        y = call_function(fun,Tt)
    y = y * (1d + t^2) / (1d - t^2)^2
endif

end ; qva_f4
;----------------------------------------------------------------------
pro qva_Vadapt,f,tinterval,rtol,atol, samples,nodes,wt,ewt,fun,a,b, Ifx,errbnd,OK,FUNCTARGS=fargs
compile_opt idl2, hidden

nint = n_elements(tinterval)
tbma = abs(tinterval[nint-1] - tinterval[0]) ; length of transformed interval

; Initialize array of subintervals of [a,b].
subs = [[tinterval[0:nint-2]],[tinterval[1:nint-1]]] ; Two columns array[n,2]

; Initialize partial sums.
IfxOK = 0d
errOK = 0d

; Initialize main loop
OK = 1 ; true
first = 1 ; true
Ifx = 0d
errbnd = 0d
iter = 0
while 1 do begin
    ; SUBS contains subintervals of [a,b] where the integral is not
    ; sufficiently accurate.  The first row of SUBS holds the left end
    ; points and the second row, the corresponding right end points.
    midpt = total(subs,2)/2d   ; midpoints of the subintervals
    halfh = diff(subs)/2d      ; half the lengths of the subintervals
    x = nodes # halfh + replicate(1d,samples) # midpt
    call_procedure,f,fun,x[*],a,b,fx,too_close,FUNCTARGS=fargs
    ;plot, x, fx, psym=4 & wait, 0.5

    ; Quit if mesh points are too close or too close to a
    ; singular point or got into trouble on first evaluation.
    not_finite = any(~finite(fx))
    if too_close or not_finite then break
    fx = reform(fx,samples,n_elements(fx)/samples)

    ; Quantities for subintervals.
    Ifxsubs = wt # fx * halfh
    errsubs = ewt # fx * halfh

    ; Quantities for all of [a,b].
    Ifx = total(Ifxsubs) + IfxOK
    errbnd = abs(total(errsubs) + errOK)
    ;print, ++iter, ifx, errbnd, format='(i5,": ",g0.7," +/- ",g0.2)'

    ; Test for convergence:
    tol = atol > rtol*abs(Ifx)
    if errbnd le tol then return

    ; Locate subintervals where the approximate integrals are
    ; sufficiently accurate and use them to update partial sums.
    ndx = where(abs(errsubs) le (2d/tbma)*halfh*tol, COMPLEMENT=bad)
    if ndx[0] ne -1 then begin
        errOK += total(errsubs[ndx])
        IfxOK += total(Ifxsubs[ndx])
        if bad[0] eq -1 then return ; all intervals are accurate
    endif
    ; Only keep subintervals still without accurate approximations.
    subs = subs[bad,*]

    ; Split the remaining subintervals in half. Quit if splitting
    ; results in too many subintervals.
    many_subint = 2d*n_elements(bad) gt 650 ; multiplied limit by 10x MC 26/FEB/2008 
    if many_subint then break 
    midpt = total(subs,2)/2d
    tmp = [[subs[*,0]], [midpt], [midpt], [subs[*,1]]]
    subs = transpose(reform(transpose(tmp),2,n_elements(tmp)/2)) ; ---> subs[n,2]
    first = 0
endwhile
if first then begin
    if too_close then print, '***Sub intervals too close.' $
    else if not_finite then print, '***Infinite values in integrand.' $
    else if many_subint then print, '***Too many sub intervals.'
    OK = 0
endif

end ;qva_Vadapt
;----------------------------------------------------------------------
pro quadva, fun, interval, Ifx, errbnd, $
    RELTOL=reltol, ABSTOL=abstol, FUNCTARGS=fargs, STATUS=status
compile_opt idl2
on_error, 2

nint = n_elements(interval)
if nint lt 2 then $
    message, 'INTERVAL must be a real vector of at least two entries.'
if any(diff(interval) le 0) then $
    message, 'Entries of INTERVAL must strictly increase.'

a = interval[0]
b = interval[nint-1]

; Generally the error test is a mixed one, but pure absolute error
; and pure relative error are allowed.  If a pure relative error
; test is specified, the tolerance must be at least 100*EPS. Defaults
; are 1e-5 for relative error, 1e-10 for absolute error.
;
if n_elements(reltol) le 0 then $
    rtol = 1d-5 $
else if reltol lt 0 then $
    rtol = 0d $
else $
    rtol = reltol > 100d*(machar(/DOUBLE)).eps
if n_elements(abstol) le 0 then $
    atol = 1d-10 $
else $
    atol = abstol > 0d
if atol+rtol eq 0 then begin
    rtol = 1d-5
    atol = 1d-10
endif

; Gauss-Kronrod (7,15) pair. Use symmetry in defining nodes and weights.
;
samples = 15
pnodes = [0.2077849550078985d, 0.4058451513773972d, 0.5860872354676911d, $
          0.7415311855993944d, 0.8648644233597691d, 0.9491079123427585d, $
          0.9914553711208126d]
nodes = [-reverse(pnodes), 0d, pnodes]
pwt = [0.2044329400752989d, 0.1903505780647854d, 0.1690047266392679d, $
       0.1406532597155259d, 0.1047900103222502d, 0.06309209262997855d, $
       0.02293532201052922d]
wt = [reverse(pwt), 0.2094821410847278d, pwt]
pwt7 = [0d, 0.3818300505051189d, 0d, 0.2797053914892767d, $
        0d, 0.1294849661688697d, 0d]
ewt = wt - [reverse(pwt7), 0.4179591836734694d, pwt7]

; Identify the task. If breakpoints are specified, work out
; how they map into the standard interval.
;
if finite(a) and finite(b) then begin
    if nint gt 2 then begin
        ; Analytical transformation suggested by K.L. Metlov:
        alpha = 2d*sin( asin((a + b - 2d*interval[1:nint-2])/(a - b))/3d )
        tinterval = [-1d, alpha, 1d]
        tinterval = qva_split(tinterval)
    endif else $
        tinterval = range(-1d,1d,11)
    qva_Vadapt,'qva_f1',tinterval,rtol,atol,samples,nodes,wt,ewt,fun,a,b,Ifx,errbnd,OK,FUNCTARGS=fargs
endif else if finite(a) and finite(b,/INFINITY) then begin
    if nint gt 2 then begin
        alpha = sqrt(interval[1:nint-2] - a)
        tinterval = [0d, alpha/(1d + alpha), 1d]
        tinterval = qva_split(tinterval)
    endif else $
        tinterval = range(0d,1d,11)
    qva_Vadapt,'qva_f2',tinterval,rtol,atol,samples,nodes,wt,ewt,fun,a,b,Ifx,errbnd,OK,FUNCTARGS=fargs
endif else if finite(a,/INFINITY) and finite(b) then begin
    if nint gt 2 then begin
        alpha = sqrt(b - interval[1:nint-2])
        tinterval = [-1d, -alpha/(1d + alpha), 0d]
        tinterval = qva_split(tinterval)
    endif else $
        tinterval = range(-1d,0d,11)
    qva_Vadapt,'qva_f3',tinterval,rtol,atol,samples,nodes,wt,ewt,fun,a,b,Ifx,errbnd,OK,FUNCTARGS=fargs
endif else if finite(a,/INFINITY) and finite(b,/INFINITY) then begin
    if nint gt 2 then begin
        ; Analytical transformation suggested by K.L. Metlov:
        alpha = tanh( asinh(2d*interval[1:nint-2])/2d )
        tinterval = [-1d, alpha, +1d]
        tinterval = qva_split(tinterval)
    endif else $
        tinterval = range(-1d,1d,11)
    qva_Vadapt,'qva_f4',tinterval,rtol,atol,samples,nodes,wt,ewt,fun,a,b,Ifx,errbnd,OK,FUNCTARGS=fargs
endif

if ~OK then begin
   print, '***Integral does not satisfy error test.'
   print, '***Approximate bound on error is', errbnd
   status = 1
endif else status = 0 ; success   

end ; quadva
;######################################################################
; All the routines below this line are just examples
function quadva_test1, x
; Gladwell's problem no1. limits x=[0,8]
; Precise result: 0.33333333332074955152

return, exp(-3d*x)-cos(5d*!dpi*x)
end
;----------------------------------------------------------------------
function quadva_test2, x
; Gladwell's problem no2. limits x=[-1,2]
; Precise result: 5.9630898453302550932

return, abs(x-1d/sqrt(3d)) + abs(x+1d/sqrt(2))
end
;----------------------------------------------------------------------
function quadva_test3, x, a=a
; Gladwell's problem no3. limits x=[0,1]
; Precise result: 3 with x^(-2d/3d)

return, x^a
end
;----------------------------------------------------------------------
function quadva_test4, x
; Gladwell's problem no3. limits x=[0,1]
; Precise result: 2.5066282746310005024

return, exp(-x^2/2d)
end
;----------------------------------------------------------------------
pro quadva_examples

print, 'test 1 ###############################'
quadva, 'quadva_test1', [0d,8d], RELTOL=0d, ABSTOL=1d-12
print, 'test 2 ###############################'
quadva, 'quadva_test2', [-1d,2d], RELTOL=0d, ABSTOL=1d-12
print, 'test 2 with breakpoints ##############'
quadva, 'quadva_test2', [-1,-1/sqrt(2),1/sqrt(3),2], RELTOL=0d, ABSTOL=1d-12
print, 'test 3 ###############################'
quadva, 'quadva_test3', [0d,1d], RELTOL=0d, ABSTOL=1d-12, FUNCTARGS={a:-2d/3d}

quadva, 'quadva_test4', [-!VALUES.F_INFINITY,!VALUES.F_INFINITY], int, err, RELTOL=1d-5, ABSTOL=0
quadva, 'quadva_test4', [-6d,6d], int, err, RELTOL=1d-5, ABSTOL=0
end
;----------------------------------------------------------------------
