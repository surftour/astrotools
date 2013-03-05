

pro doit, junk


;contour_allstars, "/raid4/tcox/vc3vc3i", 30, 50.0, 'ps', filename='/raid4/tcox/vc3vc3i/shells.eps', /particlesonly
;contour_allstars, "/raid4/tcox/vc3vc3i_rem1", 440, 50.0, 'ps', filename='/raid4/tcox/vc3vc3i_rem1/shells.eps', /particlesonly
;contour_allstars, "/raid4/tcox/vc3vc3i_rem2", 600, 50.0, 'ps', filename='/raid4/tcox/vc3vc3i_rem2/shells.eps', /particlesonly

;contour_allstars, "/raid4/tcox/vc3vc3h", 30, 50.0, 'ps', filename='/raid4/tcox/vc3vc3h/shells.eps', /particlesonly, /use_calc_center, /pubstyle
;contour_allstars, "/raid4/tcox/vcs_remergers/vc3rem_vc3rem", 40, 50.0, 'ps', filename='/raid4/tcox/vcs_remergers/vc3rem_vc3rem/shells.eps', /particlesonly, /use_calc_center, /pubstyle

;contour_allstars, "/raid4/tcox/vc3vc3h", 30, 50.0, 'ps', filename='/raid4/tcox/vc3vc3h/shells2.eps', /use_calc_center, /pubstyle
;contour_allstars, "/raid4/tcox/vcs_remergers/vc3rem_vc3rem", 40, 50.0, 'ps', filename='/raid4/tcox/vcs_remergers/vc3rem_vc3rem/shells2.eps', /use_calc_center, /pubstyle

contour_allstars, "/raid4/tcox/vc3vc3h", 30, 75.0, 'ps', filename='/raid4/tcox/vc3vc3h/shells3.eps', /use_calc_center, /pubstyle, /nolabels
contour_allstars, "/raid4/tcox/vcs_remergers/vc3rem_vc3rem", 40, 75.0, 'ps', filename='/raid4/tcox/vcs_remergers/vc3rem_vc3rem/shells3.eps', /use_calc_center, /pubstyle, /nolabels



; .run shells_find
doshells, "/raid4/tcox/vc3vc3h", 30, filename='/raid4/tcox/vc3vc3h/shells_p.eps'
doshells, "/raid4/tcox/vcs_remergers/vc3rem_vc3rem", 40, filename='/raid4/tcox/vcs_remergers/vc3rem_vc3rem/shells_p.eps'


end



