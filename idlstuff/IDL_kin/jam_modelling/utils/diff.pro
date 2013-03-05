;----------------------------------------------------------------------
function diff, a
;
; Emulates the matlab function diff(a).
; This function works for both vectors and arrays.
; Michele Cappellari, Oxford, 23 October 2007

s = size(a)
case s[0] of
1 : d = (a-shift(a,1))[1:*]
2 : d = (a-shift(a,0,1))[*,1:*]
endcase

return, d
end
;----------------------------------------------------------------------
