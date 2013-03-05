pro doit, junk


snapnum= 52
center=[0.0, 0.0, 0.0]
;xlen= 10.0
xlen= 20.0


; manual tasks to do
; -------------------
;
; the eta, and v_w need to
; be manually set in
; contour_makeplot



; ----------------------

; eta= 2.0, v_w= 800
frun= '/raid4/tcox/sbw/sb8'
snapnum= 76
contour_gas, frun, snapnum, xlen, 'ps', filename='f2b.eps', center=center, /nolabels

; eta= 0.5, v_w= 800
frun= '/raid4/tcox/sbw/sb10'
contour_gas, frun, snapnum, xlen, 'ps', filename='f2d.eps', center=center, /nolabels

; eta= 0.05, v_w= 800
frun= '/raid4/tcox/sbw/sb9'
snapnum= 42
contour_gas, frun, snapnum, xlen, 'ps', filename='f2f.eps', center=center, /nolabels

; ----------------------

; eta= 2.0, v_w= 105
frun= '/raid4/tcox/sbw/sb14'
contour_gas, frun, snapnum, xlen, 'ps', filename='f2a.eps', center=center, /nolabels

; eta= 0.5, v_w= 105
frun= '/raid4/tcox/sbw/sb17'
contour_gas, frun, snapnum, xlen, 'ps', filename='f2c.eps', center=center, /nolabels

; eta= 0.05, v_w= 105
frun= '/raid4/tcox/sbw/sb7'
snapnum= 39   ;42
contour_gas, frun, snapnum, xlen, 'ps', filename='f2e.eps', center=center, /nolabels


end








pro doit2, junk


frun= '/raid4/tcox/sbw/sb10'
snapnum= 52    ; T= 1.2  (near the merger)
center=[0.0, 0.0, 0.0]
xlen= 20.0


;========================================
;   5 panel version
;========================================

; t= 0.200 Gyr/h
;contour_gas, frun, 2, xlen, 'ps', filename='sb10gas2.eps', center=center, /nolabels
; t= 0.400 Gyr/h
;contour_gas, frun, 9, xlen, 'ps', filename='sb10gas9.eps', center=center, /nolabels
; t= 1.100 Gyr/h
;contour_gas, frun, 31, xlen, 'ps', filename='sb10gas31.eps', center=center, /nolabels
; t= 1.200 Gyr/h
;contour_gas, frun, 51, xlen, 'ps', filename='sb10gas51.eps', center=center, /nolabels
; t= 1.700 Gyr/h
;contour_gas, frun, 93, xlen, 'ps', filename='sb10gas93.eps', center=center, /nolabels

; ----------------------

; t= 0.200 Gyr/h
;contour_allstars, frun, 2, xlen, 'ps', filename='sb10stars2.eps', center=center, /nolabels
; t= 0.400 Gyr/h
;contour_allstars, frun, 9, xlen, 'ps', filename='sb10stars9.eps', center=center, /nolabels
; t= 1.100 Gyr/h
;contour_allstars, frun, 31, xlen, 'ps', filename='sb10stars31.eps', center=center, /nolabels
; t= 1.200 Gyr/h
;contour_allstars, frun, 51, xlen, 'ps', filename='sb10stars51.eps', center=center, /nolabels
; t= 1.700 Gyr/h
;contour_allstars, frun, 93, xlen, 'ps', filename='sb10stars93.eps', center=center, /nolabels

; ----------------------




;========================================
;   6 panel version
;========================================

; approximate times (in Gyr/h)
;
; 0.180833     0.542500     0.904169      1.26583      1.62750      1.98917
;
;

;contour_gas, frun, 2, xlen, 'ps', filename='sb10gas2.eps', center=center, /nolabels
;contour_gas, frun, 12, xlen, 'ps', filename='sb10gas12.eps', center=center, /nolabels
contour_gas, frun, 19, xlen, 'ps', filename='sb10gas19.eps', center=center, /nolabels
contour_gas, frun, 64, xlen, 'ps', filename='sb10gas64.eps', center=center, /nolabels
;contour_gas, frun, 92, xlen, 'ps', filename='sb10gas92.eps', center=center, /nolabels
;contour_gas, frun, 96, xlen, 'ps', filename='sb10gas96.eps', center=center, /nolabels

; ----------------------

;contour_allstars, frun, 2, xlen, 'ps', filename='sb10stars2.eps', center=center, /nolabels
;contour_allstars, frun, 12, xlen, 'ps', filename='sb10stars12.eps', center=center, /nolabels
;contour_allstars, frun, 19, xlen, 'ps', filename='sb10stars19.eps', center=center, /nolabels
;contour_allstars, frun, 64, xlen, 'ps', filename='sb10stars64.eps', center=center, /nolabels
;contour_allstars, frun, 92, xlen, 'ps', filename='sb10stars92.eps', center=center, /nolabels
;contour_allstars, frun, 96, xlen, 'ps', filename='sb10stars96.eps', center=center, /nolabels

; ----------------------



end








