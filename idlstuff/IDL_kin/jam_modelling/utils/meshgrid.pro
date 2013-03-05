;----------------------------------------------------------------------------
pro meshgrid, x, y, xx, yy
;
; Emulate MATLAB function MESHGRID.
; By Michele Cappellari, Leiden, 14 February 2004
;
COMPILE_OPT IDL2

nx = n_elements(x)
ny = n_elements(y)
xx = x # replicate(1,ny)
yy = replicate(1,nx) # y

END
;----------------------------------------------------------------------------
