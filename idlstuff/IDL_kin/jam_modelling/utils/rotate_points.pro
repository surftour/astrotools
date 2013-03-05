;----------------------------------------------------------------------
pro rotate_points, x, y, ang, xNew, yNew
;
; Rotates points conter-clockwise by an angle ANG in degrees.
; Michele cappellari, Leiden, 15 August 2003

theta = ang/!RADEG
xNew = x*COS(theta) - y*SIN(theta)
yNew = x*SIN(theta) + y*COS(theta)

end
;----------------------------------------------------------------------
