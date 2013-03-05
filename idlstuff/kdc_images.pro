

pro doit, junk

;  three panel figs
; -------------------
;.run time_m_sigma_re
;.run kdc_slitvel

;multi_slit_plot_3, "/raid4/tcox/vc3vc3i", 0.0, 0.0, snapnum= 30, filename='kdc_i_0.eps'

;multi_slit_plot_3, "/raid4/tcox/collisionless/cvc3vc3i", 0.0, 0.0, snapnum= 30, filename='kdc_ci_0.eps'


;multi_slit_plot_3, "/raid4/tcox/vc3vc3i_rem1", 0.0, 0.0, snapnum= 440, filename='/raid4/tcox/vc3vc3i_rem1/kdc_0.eps'
;multi_slit_plot_3, "/raid4/tcox/vc3vc3i_rem1", 90.0, 0.0, snapnum= 440, filename='/raid4/tcox/vc3vc3i_rem1/kdc_1.eps'
;multi_slit_plot_3, "/raid4/tcox/vc3vc3i_rem1", 90.0, 90.0, snapnum= 440, filename='/raid4/tcox/vc3vc3i_rem1/kdc_2.eps'

;multi_slit_plot_3, "/raid4/tcox/vc3vc3i_rem2", 0.0, 0.0, snapnum= 600, filename='/raid4/tcox/vc3vc3i_rem2/kdc_0.eps'
;multi_slit_plot_3, "/raid4/tcox/vc3vc3i_rem2", 90.0, 0.0, snapnum= 600, filename='/raid4/tcox/vc3vc3i_rem2/kdc_1.eps'
;multi_slit_plot_3, "/raid4/tcox/vc3vc3i_rem2", 90.0, 90.0, snapnum= 600, filename='/raid4/tcox/vc3vc3i_rem2/kdc_2.eps'

;multi_slit_plot_3, "/raid4/tcox/vc3vc3i_s1", 0.0, 0.0, snapnum= 30, filename='/raid4/tcox/vc3vc3i_s1/kdc_0.eps'
;multi_slit_plot_3, "/raid4/tcox/vc3vc3i_s1", 90.0, 0.0, snapnum= 30, filename='/raid4/tcox/vc3vc3i_s1/kdc_1.eps'
;multi_slit_plot_3, "/raid4/tcox/vc3vc3i_s1", 90.0, 90.0, snapnum= 30, filename='/raid4/tcox/vc3vc3i_s1/kdc_2.eps'



frun="/raid4/tcox/bs/b3e_i3"
snapnum= 13
multi_slit_plot_3, frun, 0.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_13_0.eps'
multi_slit_plot_3, frun, 90.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_13_1.eps'
multi_slit_plot_3, frun, 90.0, 90.0, snapnum= snapnum, filename=frun+'/kdc_13_2.eps'
snapnum= 14
multi_slit_plot_3, frun, 0.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_14_0.eps'
multi_slit_plot_3, frun, 90.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_14_1.eps'
multi_slit_plot_3, frun, 90.0, 90.0, snapnum= snapnum, filename=frun+'/kdc_14_2.eps'
snapnum= 25
multi_slit_plot_3, frun, 0.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_25_0.eps'
multi_slit_plot_3, frun, 90.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_25_1.eps'
multi_slit_plot_3, frun, 90.0, 90.0, snapnum= snapnum, filename=frun+'/kdc_25_2.eps'


frun="/raid4/tcox/bs/b3e"
snapnum= 13
multi_slit_plot_3, frun, 0.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_13_0.eps'
multi_slit_plot_3, frun, 90.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_13_1.eps'
multi_slit_plot_3, frun, 90.0, 90.0, snapnum= snapnum, filename=frun+'/kdc_13_2.eps'
snapnum= 14
multi_slit_plot_3, frun, 0.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_14_0.eps'
multi_slit_plot_3, frun, 90.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_14_1.eps'
multi_slit_plot_3, frun, 90.0, 90.0, snapnum= snapnum, filename=frun+'/kdc_14_2.eps'
snapnum= 25
multi_slit_plot_3, frun, 0.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_25_0.eps'
multi_slit_plot_3, frun, 90.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_25_1.eps'
multi_slit_plot_3, frun, 90.0, 90.0, snapnum= snapnum, filename=frun+'/kdc_25_2.eps'




end




pro do_all_three_slits, frun, snapnum


snaplbl= '0000' + strcompress(string(snapnum),/remove_all)
snaplbl= strmid(snaplbl,strlen(snaplbl)-3,3)


multi_slit_plot_3, frun, 0.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_'+snaplbl+'_0.eps'
multi_slit_plot_3, frun, 90.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_'+snaplbl+'_1.eps'
multi_slit_plot_3, frun, 90.0, 90.0, snapnum= snapnum, filename=frun+'/kdc_'+snaplbl+'_2.eps'


end



pro do_galbreakdown, frun, snapnum


snaplbl= '0000' + strcompress(string(snapnum),/remove_all)
snaplbl= strmid(snaplbl,strlen(snaplbl)-3,3)

multi_slit_plot_4, frun, 0.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_gals_'+snaplbl+'_0.eps'
multi_slit_plot_4, frun, 90.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_gals_'+snaplbl+'_1.eps'
multi_slit_plot_4, frun, 90.0, 90.0, snapnum= snapnum, filename=frun+'/kdc_gals_'+snaplbl+'_2.eps'


end



pro do_starbreakdown, frun, snapnum

snaplbl= '0000' + strcompress(string(snapnum),/remove_all)
snaplbl= strmid(snaplbl,strlen(snaplbl)-3,3)


multi_slit_plot_5, frun, 0.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_newold_'+snaplbl+'_0.eps'
multi_slit_plot_5, frun, 90.0, 0.0, snapnum= snapnum, filename=frun+'/kdc_newold_'+snaplbl+'_1.eps'
multi_slit_plot_5, frun, 90.0, 90.0, snapnum= snapnum, filename=frun+'/kdc_newold_'+snaplbl+'_2.eps'


end

