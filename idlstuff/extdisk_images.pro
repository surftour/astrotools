

;
;
;
;

;gals= {frun:"data2/Sbc_q1", snapnum: 0, xlen: 120.0}      ; q= 1.0
;gals= replicate(gals,8)
;gals[*].snapnum= [0, 2, 3, 5, 10, 15, 20, 28]
;plot_contour, gals, "gas", panels='4x2', filename="gas_Sbc_q1.eps", msg=' '



;
;  youngstar maps with gas contours overlaid.
;
; needed to manually change "do_contour" in image/img_do_one_panel
;
;
;
;

;gals= {frun:"data2/Sbc201_q1", snapnum: 70, xlen: 120.0}      ; q= 1.0
;plot_contour, gals, "youngstars", panels='1x1', filename="map_Sbc_q1_70.eps", msg=' ', colortable=3
;plot_contour, gals, "gas", panels='1x1', filename="map_Sbc_q1_gas70.eps", msg=' ', colortable=1

;gals= {frun:"data2/Sbcm1201", snapnum: 70, xlen: 120.0}      ; q= 1.0
;plot_contour, gals, "youngstars", panels='1x1', filename="map_Sbcm1_70.eps", msg=' ', colortable=3
;plot_contour, gals, "gas", panels='1x1', filename="map_Sbcm1_gas70.eps", msg=' ', colortable=1

;gals= {frun:"data2/Sbcm2201", snapnum: 70, xlen: 120.0}      ; q= 1.0
;plot_contour, gals, "youngstars", panels='1x1', filename="map_Sbcm2_70.eps", msg=' ', colortable=3
;plot_contour, gals, "gas", panels='1x1', filename="map_Sbcm2_gas70.eps", msg=' ', colortable=1

;gals= {frun:"data2/Sbcm3201", snapnum: 70, xlen: 120.0}      ; q= 1.0
;plot_contour, gals, "youngstars", panels='1x1', filename="map_Sbcm3_70.eps", msg=' ', colortable=3
;plot_contour, gals, "gas", panels='1x1', filename="map_Sbcm3_gas70.eps", msg=' ', colortable=1

;gals= {frun:"data2/Sbccut201", snapnum: 70, xlen: 120.0}      ; q= 1.0
;plot_contour, gals, "youngstars", panels='1x1', filename="map_Sbccut_70.eps", msg=' ', colortable=3
;plot_contour, gals, "gas", panels='1x1', filename="map_Sbccut_gas70.eps", msg=' ', colortable=1




;
;
;  finding the tidal dwarfs
;
;
;

;gals= {frun:"data2/Sbc201_q1", snapnum: 300, xlen: 120.0}
;plot_contour, gals, "newstars", filename="map_Sbc_ns300.eps", msg=' ', /threeprojection

;gals= {frun:"data2/Sbccut_q1", snapnum: 300, xlen: 120.0}
;plot_contour, gals, "newstars", filename="map_Sbccut_ns300.eps", msg=' ', /threeprojection

;gals= {frun:"data2/Sbcm2201", snapnum: 300, xlen: 120.0}
;plot_contour, gals, "newstars", filename="map_Sbcm2_ns300.eps", msg=' ', /threeprojection

;gals= {frun:"data2/Sbcm1201", snapnum: 300, xlen: 120.0}
;plot_contour, gals, "newstars", filename="map_Sbcm1_ns300.eps", msg=' ', /threeprojection

;gals= {frun:"data2/Sbcm3201", snapnum: 300, xlen: 120.0}
;plot_contour, gals, "newstars", filename="map_Sbcm3_ns300.eps", msg=' ', /threeprojection





end



