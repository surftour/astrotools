pro doit, junk


frun= '/raid4/tcox/vc3vc3e_2'

; manually set scale to
;  set_maxden= 10.0
;  set_dynrng= 1.0e+4
;
; and set timelbl to have 3 sig. digs

snapnum= 10
center=[17.4,24.2,0.4]
contour_gas, frun, snapnum, 2.0, 'ps', filename='snap_10_innergas.eps', /showbhs, center=center
contour_gas, frun, snapnum, 0.2, 'ps', filename='snap_10_veryinnergas.eps', /showbhs, center=center

snapnum= 18
center=[14.3,18.3,0.3]
contour_gas, frun, snapnum, 2.0, 'ps', filename='snap_18_innergas.eps', /showbhs, center=center
contour_gas, frun, snapnum, 0.2, 'ps', filename='snap_18_veryinnergas.eps', /showbhs, center=center

snapnum= 32
center=[-2.4,-2.3,1.9]
contour_gas, frun, snapnum, 2.0, 'ps', filename='snap_32_innergas.eps', /showbhs, center=center
contour_gas, frun, snapnum, 0.2, 'ps', filename='snap_32_veryinnergas.eps', /showbhs, center=center

snapnum= 52
contour_gas, frun, snapnum, 2.0, 'ps', filename='snap_52_innergas.eps', /showbhs
contour_gas, frun, snapnum, 0.2, 'ps', filename='snap_52_veryinnergas.eps', /showbhs, /use_calc_center

snapnum= 54
contour_gas, frun, snapnum, 2.0, 'ps', filename='snap_54_innergas.eps', /showbhs
contour_gas, frun, snapnum, 0.2, 'ps', filename='snap_54_veryinnergas.eps', /showbhs, /use_calc_center

snapnum= 72
contour_gas, frun, snapnum, 2.0, 'ps', filename='snap_72_innergas.eps', /showbhs
contour_gas, frun, snapnum, 0.2, 'ps', filename='snap_72_veryinnergas.eps', /showbhs, /use_calc_center


end









pro doit2, junk


frun= '/raid4/tcox/vc3vc3e_2'

; manually set scale to
;  set_maxden= 10.0
;  set_dynrng= 1.0e+6
;
; and set timelbl to have 3 sig. digs

snapnum= 10
center=[17.4,24.2,0.4]
;contour_allstars, frun, snapnum, 2.0, 'ps', filename='snap_10_allstars.eps', /showbhs, center=center
contour_newstars, frun, snapnum, 2.0, 'ps', filename='snap_10_newstars.eps', /showbhs, center=center

snapnum= 18
center=[14.3,18.3,0.3]
;contour_allstars, frun, snapnum, 2.0, 'ps', filename='snap_18_allstars.eps', /showbhs, center=center
contour_newstars, frun, snapnum, 2.0, 'ps', filename='snap_18_newstars.eps', /showbhs, center=center

snapnum= 32
center=[-2.4,-2.3,1.9]
;contour_allstars, frun, snapnum, 2.0, 'ps', filename='snap_32_allstars.eps', /showbhs, center=center
contour_newstars, frun, snapnum, 2.0, 'ps', filename='snap_32_newstars.eps', /showbhs, center=center

snapnum= 52
;contour_allstars, frun, snapnum, 2.0, 'ps', filename='snap_52_allstars.eps', /showbhs
contour_newstars, frun, snapnum, 2.0, 'ps', filename='snap_52_newstars.eps', /showbhs, /use_calc_center

snapnum= 54
;contour_allstars, frun, snapnum, 2.0, 'ps', filename='snap_54_allstars.eps', /showbhs
contour_newstars, frun, snapnum, 2.0, 'ps', filename='snap_54_newstars.eps', /showbhs, /use_calc_center

snapnum= 72
;contour_allstars, frun, snapnum, 2.0, 'ps', filename='snap_72_allstars.eps', /showbhs
contour_newstars, frun, snapnum, 2.0, 'ps', filename='snap_72_newstars.eps', /showbhs, /use_calc_center


end


