pro doit, junk


;frun="/raid4/tcox/vc3vc3"
;frun="/raid4/tcox/vc3vc3_no"
;contour_xrays, frun, 2, 50.0, 'ps', filename=frun+"/cntst_x2.eps", /nolabels, /pubstyle, msg=' '
;contour_xrays, frun, 4, 50.0, 'ps', filename=frun+"/cntst_x4.eps", /nolabels, /pubstyle, msg=' '
;contour_xrays, frun, 6, 50.0, 'ps', filename=frun+"/cntst_x6.eps", /nolabels, /pubstyle, msg=' '
;contour_xrays, frun, 8, 50.0, 'ps', filename=frun+"/cntst_x8.eps", /nolabels, /pubstyle, msg=' '
;contour_xrays, frun, 10, 50.0, 'ps', filename=frun+"/cntst_x10.eps", /nolabels, /pubstyle, msg=' '
;contour_xrays, frun, 11, 50.0, 'ps', filename=frun+"/cntst_x11.eps", /nolabels, /pubstyle, msg=' '
;contour_xrays, frun, 12, 50.0, 'ps', filename=frun+"/cntst_x12.eps", /nolabels, /pubstyle, msg=' '
;contour_xrays, frun, 15, 50.0, 'ps', filename=frun+"/cntst_x15.eps", /nolabels, /pubstyle, msg=' '
;contour_xrays, frun, 20, 50.0, 'ps', filename=frun+"/cntst_x20.eps", /nolabels, /pubstyle, msg=' '


;frun="pool/vc3vc3_wBH"
;contour_xrays, frun, 2, 50.0, 'ps', filename=frun+"/cnt_x2.eps"
;contour_xrays, frun, 4, 50.0, 'ps', filename=frun+"/cnt_x4.eps"
;contour_xrays, frun, 6, 50.0, 'ps', filename=frun+"/cnt_x6.eps"
;contour_xrays, frun, 8, 50.0, 'ps', filename=frun+"/cnt_x8.eps"
;contour_xrays, frun, 10, 50.0, 'ps', filename=frun+"/cnt_x10.eps"
;contour_xrays, frun, 11, 50.0, 'ps', filename=frun+"/cnt_x11.eps"
;contour_xrays, frun, 12, 50.0, 'ps', filename=frun+"/cnt_x12.eps"
;contour_xrays, frun, 15, 50.0, 'ps', filename=frun+"/cnt_x15.eps"
;contour_xrays, frun, 20, 50.0, 'ps', filename=frun+"/cnt_x20.eps"


;frun="pool/vc3bvc3b"
;contour_xrays, frun, 2, 50.0, 'ps', filename=frun+"/cnt_x2.eps"
;contour_xrays, frun, 4, 50.0, 'ps', filename=frun+"/cnt_x4.eps"
;contour_xrays, frun, 6, 50.0, 'ps', filename=frun+"/cnt_x6.eps"
;contour_xrays, frun, 8, 50.0, 'ps', filename=frun+"/cnt_x8.eps"
;contour_xrays, frun, 10, 50.0, 'ps', filename=frun+"/cnt_x10.eps"
;contour_xrays, frun, 11, 50.0, 'ps', filename=frun+"/cnt_x11.eps"
;contour_xrays, frun, 12, 50.0, 'ps', filename=frun+"/cnt_x12.eps"
;contour_xrays, frun, 15, 50.0, 'ps', filename=frun+"/cnt_x15.eps"
;contour_xrays, frun, 20, 50.0, 'ps', filename=frun+"/cnt_x20.eps"


;frun="pool/vc3bvc3b_wBH"
;contour_xrays, frun, 2, 50.0, 'ps', filename=frun+"/cnt_x2.eps"
;contour_xrays, frun, 4, 50.0, 'ps', filename=frun+"/cnt_x4.eps"
;contour_xrays, frun, 6, 50.0, 'ps', filename=frun+"/cnt_x6.eps"
;contour_xrays, frun, 8, 50.0, 'ps', filename=frun+"/cnt_x8.eps"
;contour_xrays, frun, 10, 50.0, 'ps', filename=frun+"/cnt_x10.eps"
;contour_xrays, frun, 11, 50.0, 'ps', filename=frun+"/cnt_x11.eps"
;contour_xrays, frun, 12, 50.0, 'ps', filename=frun+"/cnt_x12.eps"
;contour_xrays, frun, 15, 50.0, 'ps', filename=frun+"/cnt_x15.eps"
;contour_xrays, frun, 20, 50.0, 'ps', filename=frun+"/cnt_x20.eps"



end








pro doit2, junk


frun="pool/vc3vc3"
contour_phasediagram, frun, 2, 'ps', filename=frun+"/phase2.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 4, 'ps', filename=frun+"/phase4.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 6, 'ps', filename=frun+"/phase6.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 8, 'ps', filename=frun+"/phase8.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 10, 'ps', filename=frun+"/phase10.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 11, 'ps', filename=frun+"/phase11.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 12, 'ps', filename=frun+"/phase12.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 15, 'ps', filename=frun+"/phase15.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 20, 'ps', filename=frun+"/phase20.eps", sfrhothresh=0.00854924


frun="pool/vc3vc3_wBH"
contour_phasediagram, frun, 2, 'ps', filename=frun+"/phase2.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 4, 'ps', filename=frun+"/phase4.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 6, 'ps', filename=frun+"/phase6.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 8, 'ps', filename=frun+"/phase8.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 10, 'ps', filename=frun+"/phase10.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 11, 'ps', filename=frun+"/phase11.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 12, 'ps', filename=frun+"/phase12.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 15, 'ps', filename=frun+"/phase15.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 20, 'ps', filename=frun+"/phase20.eps", sfrhothresh=0.00854924


frun="pool/vc3bvc3b"
contour_phasediagram, frun, 2, 'ps', filename=frun+"/phase2.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 4, 'ps', filename=frun+"/phase4.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 6, 'ps', filename=frun+"/phase6.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 8, 'ps', filename=frun+"/phase8.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 10, 'ps', filename=frun+"/phase10.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 11, 'ps', filename=frun+"/phase11.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 12, 'ps', filename=frun+"/phase12.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 15, 'ps', filename=frun+"/phase15.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 20, 'ps', filename=frun+"/phase20.eps", sfrhothresh=0.00854924


frun="pool/vc3bvc3b_wBH"
contour_phasediagram, frun, 2, 'ps', filename=frun+"/phase2.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 4, 'ps', filename=frun+"/phase4.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 6, 'ps', filename=frun+"/phase6.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 8, 'ps', filename=frun+"/phase8.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 10, 'ps', filename=frun+"/phase10.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 11, 'ps', filename=frun+"/phase11.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 12, 'ps', filename=frun+"/phase12.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 15, 'ps', filename=frun+"/phase15.eps", sfrhothresh=0.00854924
contour_phasediagram, frun, 20, 'ps', filename=frun+"/phase20.eps", sfrhothresh=0.00854924



end











pro doit3, junk

;frun="/raid2/tcox/vc3bvc3b"
frun="/raid2/tcox/vc3bvc3b_no"

;contour_xrays, frun, 2, 50.0, 'ps', filename=frun+"/cnt_x2.eps"
contour_xrays, frun, 2, 40.0, 'ps', filename="cnt_x2.eps", /nolabels, /pubstyle, msg=' '
;contour_xrays, frun, 4, 50.0, 'ps', filename=frun+"/cnt_x4.eps"
;contour_xrays, frun, 6, 50.0, 'ps', filename=frun+"/cnt_x6.eps"
;contour_xrays, frun, 8, 50.0, 'ps', filename=frun+"/cnt_x8.eps"
;contour_xrays, frun, 10, 50.0, 'ps', filename=frun+"/cnt_x10.eps"
;contour_xrays, frun, 11, 50.0, 'ps', filename=frun+"/cnt_x11.eps"
contour_xrays, frun, 11, 40.0, 'ps', filename="cnt_x11.eps", /nolabels, /pubstyle, msg=' '
;contour_xrays, frun, 12, 50.0, 'ps', filename=frun+"/cnt_x12.eps"
;contour_xrays, frun, 15, 50.0, 'ps', filename=frun+"/cnt_x15.eps"
;contour_xrays, frun, 20, 50.0, 'ps', filename=frun+"/cnt_x20.eps"
contour_xrays, frun, 20, 40.0, 'ps', filename="cnt_x20.eps", /nolabels, /pubstyle, msg=' ', center=[4.0,0.0,4.0]

end



pro bhcomp_junk, junk

;plot_devacouleur_profile_avg_comparison, 'ps', smoothlen=0.1, filename='dev15.eps', snapnum=15
;plot_devacouleur_profile_avg_comparison, 'ps', smoothlen=0.1, filename='dev20.eps', snapnum=20
;plot_devacouleur_profile_avg_comparison, 'ps', smoothlen=0.1, filename='dev25.eps', snapnum=25
;plot_devacouleur_profile_avg_comparison, 'ps', smoothlen=0.1, filename='dev30.eps', snapnum=30

;plot_metallicity_profile_avg_comparison, 'ps', smoothlen=0.1, filename='zcomp15.eps', snapnum=15
;plot_metallicity_profile_avg_comparison, 'ps', smoothlen=0.1, filename='zcomp20.eps', snapnum=20
;plot_metallicity_profile_avg_comparison, 'ps', smoothlen=0.1, filename='zcomp25.eps', snapnum=25
;plot_metallicity_profile_avg_comparison, 'ps', smoothlen=0.1, filename='zcomp30.eps', snapnum=30

;plot_surface_brightness_comparison, 'ps', smoothlen=0.1, filename='bcomp15.eps', snapnum=15
;plot_surface_brightness_comparison, 'ps', smoothlen=0.1, filename='bcomp20.eps', snapnum=20
;plot_surface_brightness_comparison, 'ps', smoothlen=0.1, filename='bcomp25.eps', snapnum=25
;plot_surface_brightness_comparison, 'ps', smoothlen=0.1, filename='bcomp30.eps', snapnum=30

;plot_surface_colors_comparison, 'ps', smoothlen=0.1, filename='bmv_comp15.eps', snapnum=15
;plot_surface_colors_comparison, 'ps', smoothlen=0.1, filename='bmv_comp20.eps', snapnum=20
;plot_surface_colors_comparison, 'ps', smoothlen=0.1, filename='bmv_comp25.eps', snapnum=25
;plot_surface_colors_comparison, 'ps', smoothlen=0.1, filename='bmv_comp30.eps', snapnum=30

end



pro remstuff, junk

initialize_plotinfo, 1

;plot_devacouleur_profile_avg_comparison, 'ps', smoothlen=0.1, filename='diskrem.eps', snapnum=25
plot_devacouleur_profile_avg_comparison, 'ps', smoothlen=0.1, filename='remrem.eps', snapnum=25

end



pro remimages, junk


frun="/raid4/tcox/vc3rem_vc3"
contour_allstars, frun, 2, 50.0, 'ps', filename=frun+"/stars2.eps"
contour_allstars, frun, 4, 50.0, 'ps', filename=frun+"/stars4.eps"
contour_allstars, frun, 6, 50.0, 'ps', filename=frun+"/stars6.eps"
contour_allstars, frun, 8, 50.0, 'ps', filename=frun+"/stars8.eps"
contour_allstars, frun, 10, 50.0, 'ps', filename=frun+"/stars10.eps"
contour_allstars, frun, 12, 50.0, 'ps', filename=frun+"/stars11.eps"
contour_allstars, frun, 14, 50.0, 'ps', filename=frun+"/stars12.eps"
contour_allstars, frun, 16, 50.0, 'ps', filename=frun+"/stars15.eps"
contour_allstars, frun, 20, 50.0, 'ps', filename=frun+"/stars20.eps"
contour_gas, frun, 2, 50.0, 'ps', filename=frun+"/gas2.eps"
contour_gas, frun, 4, 50.0, 'ps', filename=frun+"/gas4.eps"
contour_gas, frun, 6, 50.0, 'ps', filename=frun+"/gas6.eps"
contour_gas, frun, 8, 50.0, 'ps', filename=frun+"/gas8.eps"
contour_gas, frun, 10, 50.0, 'ps', filename=frun+"/gas10.eps"
contour_gas, frun, 12, 50.0, 'ps', filename=frun+"/gas11.eps"
contour_gas, frun, 14, 50.0, 'ps', filename=frun+"/gas12.eps"
contour_gas, frun, 16, 50.0, 'ps', filename=frun+"/gas15.eps"
contour_gas, frun, 20, 50.0, 'ps', filename=frun+"/gas20.eps"


frun="/raid4/tcox/vc3rem_vc3rem"
contour_allstars, frun, 2, 50.0, 'ps', filename=frun+"/stars2.eps"
contour_allstars, frun, 4, 50.0, 'ps', filename=frun+"/stars4.eps"
contour_allstars, frun, 6, 50.0, 'ps', filename=frun+"/stars6.eps"
contour_allstars, frun, 8, 50.0, 'ps', filename=frun+"/stars8.eps"
contour_allstars, frun, 10, 50.0, 'ps', filename=frun+"/stars10.eps"
contour_allstars, frun, 12, 50.0, 'ps', filename=frun+"/stars11.eps"
contour_allstars, frun, 14, 50.0, 'ps', filename=frun+"/stars12.eps"
contour_allstars, frun, 16, 50.0, 'ps', filename=frun+"/stars15.eps"
contour_allstars, frun, 20, 50.0, 'ps', filename=frun+"/stars20.eps"
contour_gas, frun, 2, 50.0, 'ps', filename=frun+"/gas2.eps"
contour_gas, frun, 4, 50.0, 'ps', filename=frun+"/gas4.eps"
contour_gas, frun, 6, 50.0, 'ps', filename=frun+"/gas6.eps"
contour_gas, frun, 8, 50.0, 'ps', filename=frun+"/gas8.eps"
contour_gas, frun, 10, 50.0, 'ps', filename=frun+"/gas10.eps"
contour_gas, frun, 12, 50.0, 'ps', filename=frun+"/gas11.eps"
contour_gas, frun, 14, 50.0, 'ps', filename=frun+"/gas12.eps"
contour_gas, frun, 16, 50.0, 'ps', filename=frun+"/gas15.eps"
contour_gas, frun, 20, 50.0, 'ps', filename=frun+"/gas20.eps"

end





; this is the Marijn stuff
; --------------------------
pro bands, junk

contour_lum_b, "/raid4/tcox/vc3vc3e", 2, 40.0, 'ps', filename='bimg1.eps', /fitstoo
contour_lum_b, "/raid4/tcox/vc3vc3e", 4, 40.0, 'ps', filename='bimg2.eps', /fitstoo
contour_lum_b, "/raid4/tcox/vc3vc3e", 8, 40.0, 'ps', filename='bimg3.eps', /fitstoo
contour_lum_b, "/raid4/tcox/vc3vc3e", 9, 40.0, 'ps', filename='bimg4.eps', /fitstoo
contour_lum_b, "/raid4/tcox/vc3vc3e", 10, 40.0, 'ps', filename='bimg5.eps', /fitstoo
contour_lum_b, "/raid4/tcox/vc3vc3e", 11, 40.0, 'ps', filename='bimg6.eps', /fitstoo
contour_lum_b, "/raid4/tcox/vc3vc3e", 12, 40.0, 'ps', filename='bimg7.eps', /fitstoo
contour_lum_b, "/raid4/tcox/vc3vc3e", 13, 40.0, 'ps', filename='bimg8.eps', /fitstoo
contour_lum_b, "/raid4/tcox/vc3vc3e", 14, 40.0, 'ps', filename='bimg9.eps', /fitstoo
contour_lum_b, "/raid4/tcox/vc3vc3e", 20, 40.0, 'ps', filename='bimg10.eps', /fitstoo
contour_lum_b, "/raid4/tcox/vc3vc3e", 30, 40.0, 'ps', filename='bimg11.eps', /fitstoo

contour_lum_u, "/raid4/tcox/vc3vc3e", 2, 40.0, 'ps', filename='uimg1.eps', /fitstoo
contour_lum_u, "/raid4/tcox/vc3vc3e", 4, 40.0, 'ps', filename='uimg2.eps', /fitstoo
contour_lum_u, "/raid4/tcox/vc3vc3e", 8, 40.0, 'ps', filename='uimg3.eps', /fitstoo
contour_lum_u, "/raid4/tcox/vc3vc3e", 9, 40.0, 'ps', filename='uimg4.eps', /fitstoo
contour_lum_u, "/raid4/tcox/vc3vc3e", 10, 40.0, 'ps', filename='uimg5.eps', /fitstoo
contour_lum_u, "/raid4/tcox/vc3vc3e", 11, 40.0, 'ps', filename='uimg6.eps', /fitstoo
contour_lum_u, "/raid4/tcox/vc3vc3e", 12, 40.0, 'ps', filename='uimg7.eps', /fitstoo
contour_lum_u, "/raid4/tcox/vc3vc3e", 13, 40.0, 'ps', filename='uimg8.eps', /fitstoo
contour_lum_u, "/raid4/tcox/vc3vc3e", 14, 40.0, 'ps', filename='uimg9.eps', /fitstoo
contour_lum_u, "/raid4/tcox/vc3vc3e", 20, 40.0, 'ps', filename='uimg10.eps', /fitstoo
contour_lum_u, "/raid4/tcox/vc3vc3e", 30, 40.0, 'ps', filename='uimg11.eps', /fitstoo

contour_lum_v, "/raid4/tcox/vc3vc3e", 2, 40.0, 'ps', filename='vimg1.eps', /fitstoo
contour_lum_v, "/raid4/tcox/vc3vc3e", 4, 40.0, 'ps', filename='vimg2.eps', /fitstoo
contour_lum_v, "/raid4/tcox/vc3vc3e", 8, 40.0, 'ps', filename='vimg3.eps', /fitstoo
contour_lum_v, "/raid4/tcox/vc3vc3e", 9, 40.0, 'ps', filename='vimg4.eps', /fitstoo
contour_lum_v, "/raid4/tcox/vc3vc3e", 10, 40.0, 'ps', filename='vimg5.eps', /fitstoo
contour_lum_v, "/raid4/tcox/vc3vc3e", 11, 40.0, 'ps', filename='vimg6.eps', /fitstoo
contour_lum_v, "/raid4/tcox/vc3vc3e", 12, 40.0, 'ps', filename='vimg7.eps', /fitstoo
contour_lum_v, "/raid4/tcox/vc3vc3e", 13, 40.0, 'ps', filename='vimg8.eps', /fitstoo
contour_lum_v, "/raid4/tcox/vc3vc3e", 14, 40.0, 'ps', filename='vimg9.eps', /fitstoo
contour_lum_v, "/raid4/tcox/vc3vc3e", 20, 40.0, 'ps', filename='vimg10.eps', /fitstoo
contour_lum_v, "/raid4/tcox/vc3vc3e", 30, 40.0, 'ps', filename='vimg11.eps', /fitstoo


end






pro threepro, junk


contour_three_proj, "/raid4/tcox/vc3vc3b", 30, 50.0, 'ps', filename='vc3vc3b.eps'
contour_three_proj, "/raid4/tcox/vc3vc3c", 30, 50.0, 'ps', filename='vc3vc3c.eps'
contour_three_proj, "/raid4/tcox/vc3vc3d", 30, 50.0, 'ps', filename='vc3vc3d.eps'
contour_three_proj, "/raid4/tcox/vc3vc3e", 30, 50.0, 'ps', filename='vc3vc3e.eps'
contour_three_proj, "/raid4/tcox/vc3vc3f", 30, 50.0, 'ps', filename='vc3vc3f.eps'
contour_three_proj, "/raid4/tcox/vc3vc3g", 30, 50.0, 'ps', filename='vc3vc3g.eps'
contour_three_proj, "/raid4/tcox/vc3vc3h", 30, 50.0, 'ps', filename='vc3vc3h.eps'
contour_three_proj, "/raid4/tcox/vc3vc3i", 30, 50.0, 'ps', filename='vc3vc3i.eps'
contour_three_proj, "/raid4/tcox/vc3vc3j", 30, 50.0, 'ps', filename='vc3vc3j.eps'
contour_three_proj, "/raid4/tcox/vc3vc3k", 30, 50.0, 'ps', filename='vc3vc3k.eps'
contour_three_proj, "/raid4/tcox/vc3vc3l", 30, 50.0, 'ps', filename='vc3vc3l.eps'
contour_three_proj, "/raid4/tcox/vc3vc3m", 30, 50.0, 'ps', filename='vc3vc3m.eps'
contour_three_proj, "/raid4/tcox/vc3vc3n", 30, 50.0, 'ps', filename='vc3vc3n.eps'
contour_three_proj, "/raid4/tcox/vc3vc3o", 30, 50.0, 'ps', filename='vc3vc3o.eps'
contour_three_proj, "/raid4/tcox/vc3vc3p", 30, 50.0, 'ps', filename='vc3vc3p.eps'

end




pro velfields, junk


velocity_field, "/raid4/tcox/vc3vc3e", 3, 50.0, 'ps', filename='ve_xz3a.eps', /xz
velocity_field, "/raid4/tcox/vc3vc3e", 6, 50.0, 'ps', filename='ve_xz6a.eps', /xz
velocity_field, "/raid4/tcox/vc3vc3e", 8, 50.0, 'ps', filename='ve_xz8a.eps', /xz
;velocity_field, "/raid4/tcox/vc3vc3e", 10, 50.0, 'ps', filename='ve_xz10.eps', /xz
;velocity_field, "/raid4/tcox/vc3vc3e", 11, 50.0, 'ps', filename='ve_xz11.eps', /xz
;velocity_field, "/raid4/tcox/vc3vc3e", 12, 50.0, 'ps', filename='ve_xz12.eps', /xz
;velocity_field, "/raid4/tcox/vc3vc3e", 13, 50.0, 'ps', filename='ve_xz13.eps', /xz
;velocity_field, "/raid4/tcox/vc3vc3e", 14, 50.0, 'ps', filename='ve_xz14.eps', /xz
;velocity_field, "/raid4/tcox/vc3vc3e", 15, 50.0, 'ps', filename='ve_xz15.eps', /xz
;velocity_field, "/raid4/tcox/vc3vc3e", 16, 50.0, 'ps', filename='ve_xz16.eps', /xz
;velocity_field, "/raid4/tcox/vc3vc3e", 20, 50.0, 'ps', filename='ve_xz20.eps', /xz
;velocity_field, "/raid4/tcox/vc3vc3e", 30, 50.0, 'ps', filename='ve_xz30.eps', /xz

end






pro doit_multi_slits, junk


;multi_slit_plot_2, "/raid4/tcox/vc3vc3b",  0.0,  0.0, 75.0, filename="kdc/vc3b_xy.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3b", 90.0,  0.0, 91.0, filename="kdc/vc3b_yz.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3b", 90.0, 90.0, 85.0, filename="kdc/vc3b_xz.eps"

multi_slit_plot_2, "/raid4/tcox/vc3vc3c",  0.0,  0.0, 78.0, filename="kdc/vc3c_xy.eps"
multi_slit_plot_2, "/raid4/tcox/vc3vc3c", 90.0,  0.0, 91.0, filename="kdc/vc3c_yz.eps"
multi_slit_plot_2, "/raid4/tcox/vc3vc3c", 90.0, 90.0, 88.0, filename="kdc/vc3c_xz.eps"

multi_slit_plot_2, "/raid4/tcox/vc3vc3d",  0.0,  0.0, 73.0, filename="kdc/vc3d_xy.eps"
multi_slit_plot_2, "/raid4/tcox/vc3vc3d", 90.0,  0.0, 101.0, filename="kdc/vc3d_yz.eps"
multi_slit_plot_2, "/raid4/tcox/vc3vc3d", 90.0, 90.0, -4.0, filename="kdc/vc3d_xz.eps"

multi_slit_plot_2, "/raid4/tcox/vc3vc3e",  0.0,  0.0, 65.0, filename="kdc/vc3e_xy.eps"
multi_slit_plot_2, "/raid4/tcox/vc3vc3e", 90.0,  0.0, 75.0, filename="kdc/vc3e_yz.eps"
multi_slit_plot_2, "/raid4/tcox/vc3vc3e", 90.0, 90.0, 37.0, filename="kdc/vc3e_xz.eps"

;multi_slit_plot_2, "/raid4/tcox/vc3vc3f",  0.0,  0.0, 137.0, filename="kdc/vc3f_xy.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3f", 90.0,  0.0, 48.0, filename="kdc/vc3f_yz.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3f", 90.0, 90.0, 124.0, filename="kdc/vc3f_xz.eps"

;multi_slit_plot_2, "/raid4/tcox/vc3vc3g",  0.0,  0.0, 40.0, filename="kdc/vc3g_xy.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3g", 90.0,  0.0, -28.0, filename="kdc/vc3g_yz.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3g", 90.0, 90.0, 98.0, filename="kdc/vc3g_xz.eps"

;multi_slit_plot_2, "/raid4/tcox/vc3vc3h",  0.0,  0.0, 72.0, filename="kdc/vc3h_xy.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3h", 90.0,  0.0, 90.0, filename="kdc/vc3h_yz.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3h", 90.0, 90.0, 89.0, filename="kdc/vc3h_xz.eps"

;multi_slit_plot_2, "/raid4/tcox/vc3vc3i",  0.0,  0.0, 87.0, filename="kdc/vc3i_xy.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3i", 90.0,  0.0, 127.0, filename="kdc/vc3i_yz.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3i", 90.0, 90.0, 9.0, filename="kdc/vc3i_xz.eps"

;multi_slit_plot_2, "/raid4/tcox/vc3vc3j",  0.0,  0.0, 103.0, filename="kdc/vc3j_xy.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3j", 90.0,  0.0, 110.0, filename="kdc/vc3j_yz.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3j", 90.0, 90.0, 52.0, filename="kdc/vc3j_xz.eps"

;multi_slit_plot_2, "/raid4/tcox/vc3vc3k",  0.0,  0.0, 71.0, filename="kdc/vc3k_xy.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3k", 90.0,  0.0, 80.0, filename="kdc/vc3k_yz.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3k", 90.0, 90.0, 88.0, filename="kdc/vc3k_xz.eps"

;multi_slit_plot_2, "/raid4/tcox/vc3vc3l",  0.0,  0.0, 17.0, filename="kdc/vc3l_xy.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3l", 90.0,  0.0, 50.0, filename="kdc/vc3l_yz.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3l", 90.0, 90.0, 83.0, filename="kdc/vc3l_xz.eps"

;multi_slit_plot_2, "/raid4/tcox/vc3vc3m",  0.0,  0.0, 49.0, filename="kdc/vc3m_xy.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3m", 90.0,  0.0, 48.0, filename="kdc/vc3m_yz.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3m", 90.0, 90.0, 52.0, filename="kdc/vc3m_xz.eps"

;multi_slit_plot_2, "/raid4/tcox/vc3vc3n",  0.0,  0.0, 27.0, filename="kdc/vc3n_xy.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3n", 90.0,  0.0, 10.0, filename="kdc/vc3n_yz.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3n", 90.0, 90.0, 78.0, filename="kdc/vc3n_xz.eps"

;multi_slit_plot_2, "/raid4/tcox/vc3vc3o",  0.0,  0.0, 107.0, filename="kdc/vc3o_xy.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3o", 90.0,  0.0, 91.0, filename="kdc/vc3o_yz.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3o", 90.0, 90.0, 145.0, filename="kdc/vc3o_xz.eps"

;multi_slit_plot_2, "/raid4/tcox/vc3vc3p",  0.0,  0.0, 26.0, filename="kdc/vc3p_xy.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3p", 90.0,  0.0, 65.0, filename="kdc/vc3p_yz.eps"
;multi_slit_plot_2, "/raid4/tcox/vc3vc3p", 90.0, 90.0, 77.0, filename="kdc/vc3p_xz.eps"




end



