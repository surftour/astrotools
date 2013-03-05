;----------------------------------------------------------------------
function ibetam, a, b, x
;
; Incomplete beta function defined in the same 
; way as the Mathematica function Beta[x,a,b].
; This routine works for x<1 and was tested against Mathematica.
;
; V1.0: Michele Cappellari, Oxford, 01/APR/2008
; V2.0: Use Hypergeometric function for negative parameters.
;    From equation (6.6.8) of Abramoviz & Stegun
;    http://www.nrbook.com/abramowitz_and_stegun/page_263.htm 
;    MC, Oxford, 04/APR/2008
 
if a gt 0 && b gt 0 then $
    ib = ibeta(a,b,x)*beta(a,b) $
else $ 
    ib = x^a/a * hypergeometric2F1(a, 1d - b, a + 1d, x)

return, ib 
end
;----------------------------------------------------------------------