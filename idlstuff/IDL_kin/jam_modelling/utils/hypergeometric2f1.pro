;#############################################################################
;
; Copyright (C) 2008, Michele Cappellari
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
; hf = Hypergeometric2F1(a, b, c, x)
;
; From the book "Computation of Special Functions"
;    by Shanjie Zhang and Jianming Jin
;    Copyright 1996 by John Wiley & Sons, Inc.
; The authors state:
;   "We give permission to the reader who purchases this book 
;    to incorporate any of these programs into his or her 
;    programs provided that the copyright is acknowledged."
;
;       ====================================================
;       Purpose: Compute hypergeometric function F(a,b,c,x)
;       Input :  a --- Parameter
;                b --- Parameter
;                c --- Parameter, c <> 0,-1,-2,...
;                x --- Argument   ( x < 1 )
;       Output:  HF --- F(a,b,c,x)
;       Routines called:
;            (1) GAMMA for computing gamma function
;            (2) PSI for computing psi function
;       ====================================================
;
;------------------------------------------------------------------------------
function hygfx_psi, x
compile_opt idl2, hidden
;
;       ======================================
;       Purpose: Compute PolyGamma Psi function
;       Input :  x  --- Argument of psi(x)
;       Output:  PS --- psi(x)
;       ======================================

el = 0.57721566490153286061d ; N[EulerGamma,20] in Mathematica

xa = ABS(x)
s = 0.0D0
IF (x eq fix(x) and x LE 0.0) THEN begin
  ps = 1.0D+300
  RETURN, ps
endif ELSE IF (xa eq fix(xa)) THEN begin
  n = xa
  for  k = 1, n - 1 do s = s + 1.0D0 / k
  ps = -el + s
endif ELSE IF (xa+0.5 eq fix(xa+0.5)) THEN begin
  n = xa - 0.5d 
  for  k = 1, n do s = s + 1.0d / (2*k-1)
  ps = -el + 2.0D0 * s - 1.386294361119891D0
endif else begin
  IF (xa LT 10.0) THEN begin
    n = 10 - fix(xa)
    for  k = 0, n - 1 do s = s + 1.0D0 / (xa + k)
    xa = xa + n
  endif
  x2 = 1.0D0 / (xa*xa)
  a1 = -.83333333333333333D-01
  a2 =  .83333333333333333D-02
  a3 = -.39682539682539683D-02
  a4 =  .41666666666666667D-02
  a5 = -.75757575757575758D-02
  a6 =  .21092796092796093D-01
  a7 = -.83333333333333333D-01
  a8 =  .4432598039215686D0
  ps = aLOG(xa) - .5D0 / xa + x2 * (((((((a8*x2 + a7)*x2 + a6)*x2 + a5)*  $
       x2 + a4)*x2 + a3)*x2 + a2)*x2 + a1)
  ps = ps - s
endelse
IF (x LT 0.0) then ps = ps - pi * COS(pi*x) / SIN(pi*x) - 1.0D0 / x
RETURN, ps

END
;------------------------------------------------------------------------------
function hygfx_Hypergeometric2F1, a, b, c, x
compile_opt idl2
on_error, 2

el = 0.57721566490153286061d ; N[EulerGamma,20] in Mathematica

l0 = c eq fix(c) and c LT 0.0
l1 = 1.0D0 - x LT 1.0D-15 and c - a - b LE 0.0
l2 = a eq fix(a) and a LT 0.0
l3 = b eq fix(b) and b LT 0.0
l4 = c - a eq fix(c-a) and c - a LE 0.0
l5 = c - b eq fix(c-b) and c - b LE 0.0
IF (l0 or l1) THEN message, 'The hypergeometric series is divergent'
eps = 1.0D-15
IF (x GT 0.95) then eps = 1.0D-8
IF (x eq 0.0 or a eq 0.0 OR b eq 0.0) THEN begin
  hf = 1.0D0
  return, hf
endif ELSE IF (1.0D0-x eq eps and c-a-b GT 0.0) THEN begin
  gc = gamma(c)
  gcab = gamma(c-a-b)
  gca = gamma(c-a)
  gcb = gamma(c-b)
  hf = gc * gcab / (gca*gcb)
  return, hf
endif ELSE IF (1.0D0+x LE eps and ABS(c-a+b-1.0) LE eps) THEN begin
  g0 = SQRT(!dpi) * 2.0D0^(-a)
  g1 = gamma(c)
  g2 = gamma(1.0D0+a/2.0-b)
  g3 = gamma(0.5D0+0.5*a)
  hf = g0 * g1 / (g2*g3)
  return, hf
endif ELSE IF (l2 or l3) THEN begin
  IF (l2) then nm = fix(ABS(a))
  IF (l3) then nm = fix(ABS(b))
  hf = 1.0D0
  r = 1.0D0
  for  k = 1, nm do begin
    r = r * (a+k-1.0D0) * (b+k-1.0D0) / (k*(c+k-1.0D0)) * x
    hf = hf + r
  endfor
  return, hf
endif ELSE IF (l4 or l5) THEN begin
  IF (l4) then nm = fix(ABS(c-a))
  IF (l5) then nm = fix(ABS(c-b))
  hf = 1.0D0
  r = 1.0D0
  for  k = 1, nm do begin
    r = r * (c-a+k-1.0D0) * (c-b+k-1.0D0) / (k*(c+k-1.0D0)) * x
    hf = hf + r
  endfor
  hf = (1.0D0-x)^(c-a-b) * hf
  return, hf
endif
aa = a
bb = b
x1 = x
IF (x LT 0.0D0) THEN begin
  x = x / (x-1.0D0)
  IF (c GT a and b LT a and b GT 0.0) THEN begin
    a = bb
    b = aa
  endif
  b = c - b
endif
IF (x GE 0.75D0) THEN begin
  gm = 0.0D0
  IF (ABS(c-a-b-fix(c-a-b)) LT 1.0D-15) THEN begin
    m = fix(c-a-b)
    ga = gamma(a)
    gb = gamma(b)
    gc = gamma(c)
    gam = gamma(a+m)
    gbm = gamma(b+m)
    pa = hygfx_psi(a)
    pb = hygfx_psi(b)
    IF (m ne 0) then gm = 1.0D0
    for  j = 1, ABS(m) - 1 do gm = gm * j
    rm = 1.0D0
    for  j = 1, ABS(m) do rm = rm * j
    f0 = 1.0D0
    r0 = 1.0D0
    r1 = 1.0D0
    sp0 = 0.0d0
    sp = 0.0D0
    IF (m ge 0) THEN begin
      c0 = gm * gc / (gam*gbm)
      c1 = -gc * (x-1.0D0)^m / (ga*gb*rm)
      for  k = 1, m - 1 do begin
        r0 = r0 * (a+k-1.0D0) * (b+k-1.0) / (k*(k-m)) * (1.0-x)
        f0 = f0 + r0
      endfor
      for  k = 1, m do $
        sp0 = sp0 + 1.0D0 / (a+k-1.0) + 1.0 / (b+k-1.0) - 1.0 / k
      f1 = pa + pb + sp0 + 2.0D0 * el + aLOG(1.0D0-x)
      hw = f1
      for  k = 1, 250 do begin
        sp = sp + (1.0D0-a) / (k*(a+k-1.0)) + (1.0-b) / (k*(b+k-1.0))
        sm = 0.0D0
        for  j = 1, m do $
          sm = sm + (1.0D0-a) / ((j+k)*(a+j+k-1.0)) + 1.0 / (b+j+k-1.0)
        rp = pa + pb + 2.0D0 * el + sp + sm + aLOG(1.0D0-x)
        r1 = r1 * (a+m+k-1.0D0) * (b+m+k-1.0) / (k*(m+k)) * (1.0-x)
        f1 = f1 + r1 * rp
        IF (ABS(f1-hw) LT ABS(f1)*eps) then break 
        hw = f1
      endfor
      hf = f0 * c0 + f1 * c1
    endif ELSE IF (m LT 0) THEN begin
      m = -m
      c0 = gm * gc / (ga*gb*(1.0D0-x)^m)
      c1 = -(-1)^m * gc / (gam*gbm*rm)
      for  k = 1, m - 1 do begin
        r0 = r0 * (a-m+k-1.0D0) * (b-m+k-1.0) / (k*(k-m)) * (1.0-x)
        f0 = f0 + r0
      endfor
      for  k = 1, m do sp0 = sp0 + 1.0D0 / k
      f1 = pa + pb - sp0 + 2.0D0 * el + aLOG(1.0D0-x)
      hw = f1
      for  k = 1, 250 do begin
        sp = sp + (1.0D0-a) / (k*(a+k-1.0)) + (1.0-b) / (k*(b+k-1.0))
        sm = 0.0D0
        for  j = 1, m do sm = sm + 1.0D0 / (j+k)
        rp = pa + pb + 2.0D0 * el + sp - sm + aLOG(1.0D0-x)
        r1 = r1 * (a+k-1.0D0) * (b+k-1.0) / (k*(m+k)) * (1.0-x)
        f1 = f1 + r1 * rp
        IF (ABS(f1-hw) LT ABS(f1)*eps) then break 
        hw = f1
      endfor
      hf = f0 * c0 + f1 * c1
    endif
  endif ELSE begin
    ga = gamma(a)
    gb = gamma(b)
    gc = gamma(c)
    gca = gamma(c-a)
    gcb = gamma(c-b)
    gcab = gamma(c-a-b)
    gabc = gamma(a+b-c)
    c0 = gc * gcab / (gca*gcb)
    c1 = gc * gabc / (ga*gb) * (1.0D0-x)^(c-a-b)
    hf = 0.0D0
    hw = hf
    r0 = c0
    r1 = c1
    for  k = 1, 250 do begin
      r0 = r0 * (a+k-1.0D0) * (b+k-1.0) / (k*(a+b-c+k)) * (1.0-x)
      r1 = r1 * (c-a+k-1.0D0) * (c-b+k-1.0) / (k*(c-a-b+k)) * (1.0-x)
      hf = hf + r0 + r1
      IF (ABS(hf-hw) LT ABS(hf)*eps) then break
      hw = hf
    endfor
    hf = hf + c0 + c1
  endelse
endif ELSE begin
  a0 = 1.0D0
  IF (c GT a and c LT 2.0D0*a and c GT b and c LT 2.0D0*b) THEN begin
    a0 = (1.0D0-x)^(c-a-b)
    a = c - a
    b = c - b
  endif
  hf = 1.0D0
  hw = hf
  r = 1.0D0
  for  k = 1, 250 do begin
    r = r * (a+k-1.0D0) * (b+k-1.0D0) / (k*(c+k-1.0D0)) * x
    hf = hf + r
    IF (ABS(hf-hw) LE ABS(hf)*eps) then break 
    hw = hf
  endfor
  hf = a0 * hf
endelse
IF (x1 LT 0.0D0) THEN begin
  x = x1
  c0 = 1.0D0 / (1.0D0-x)^aa
  hf = c0 * hf
endif
a = aa
b = bb
IF (k GT 120) then print, ' Warning; You should check the accuracy'
return, hf

END
;------------------------------------------------------------------------------
function Hypergeometric2F1, a, b, c, x
compile_opt idl2
on_error, 2
;
; This is just a wrapper for the actual function, 
; as the computation routine is not (yet) vectorized.

n = n_elements(x)
hyp = dblarr(n,/NOZERO)
for j=0,n-1 do hyp[j] = hygfx_Hypergeometric2F1(a,b,c,x[j])

return, hyp
END
;------------------------------------------------------------------------------
pro hypergeometric_test
;
;    ============================================================
;    Purpose: This program tests the hypergeometric function
;             F(a,b,c,x) using subroutine HYGFX
;    Input :  a --- Parameter
;             b --- Parameter
;             c --- Parameter, c <> 0,-1,-2,...
;             x --- Argument ( x � 1 )
;    Output:  HF --- F(a,b,c,x)
;    Example:
;           b = 3.30,  c = 6.70
;           a     F(a,b,c,.25)     F(a,b,c,.55)    F(a,b,c,.85)
;         ------------------------------------------------------
;         -2.5   .72356129D+00    .46961432D+00   .29106096D+00
;         -0.5   .93610145D+00    .85187390D+00   .75543187D+00
;          0.5   .10689695D+01    .11795358D+01   .13510497D+01
;          2.5   .14051563D+01    .23999063D+01   .57381566D+01
;
;           a = 3.30,  b = 6.70
;           c     F(a,b,c,.25)     F(a,b,c,.55)    F(a,b,c,.85)
;         ------------------------------------------------------
;         -5.5   .15090670D+05    .10170778D+11   .58682088D+19
;         -0.5  -.21631479D+04   -.30854772D+07  -.10217370D+13
;          0.5   .26451677D+03    .11967860D+06   .92370648D+10
;          4.5   .41946916D+01    .58092729D+02   .20396914D+05
;    ============================================================

while 1 do begin
  READ, a, b, c, x, PROMPT='Please enter a,b,c and x: '
  Hypergeometric2F1, a, b, c, x, hf
  print, 'a,b,c,x=', a, b, c, x
  print, 'F(a,b,c,x)=', hf
endwhile

END
;------------------------------------------------------------------------------
 