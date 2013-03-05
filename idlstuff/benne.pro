

;
;
;
;

;gals= {frun:"data1/Sbc/Sbc_10x_wBH", snapnum: 0, xlen: 42.0}      ; q= 0.25
;gals= {frun:"data1/Sbc/Sbc_10x_wBHq", snapnum: 0, xlen: 42.0}      ; q= 0.5
gals= {frun:"data1/Sbc/Sbc_q1", snapnum: 0, xlen: 42.0}      ; q= 1.0
gals= replicate(gals,4)
gals[*].snapnum= [5, 10, 20, 30]
plot_contour, gals, "HI", panels='4x1', filename="HI_Sbc10xwBHq.eps", msg='exp/HI'
plot_contour, gals, "gas", panels='4x1', filename="gas_Sbc10xwBHq.eps", msg='exp/total'

gals= {frun:"data1/Sbc/Sbccut", snapnum: 0, xlen: 42.0}
gals= replicate(gals,4)
gals[*].snapnum= [5, 10, 20, 30]
plot_contour, gals, "HI", panels='4x1', filename="HI_Sbccut.eps", msg='expcut/HI'
plot_contour, gals, "gas", panels='4x1', filename="gas_Sbccut.eps", msg='expcut/total'

gals= {frun:"data1/Sbc/Sbcm1_q1", snapnum: 0, xlen: 42.0}
gals= replicate(gals,4)
gals[*].snapnum= [5, 10, 20, 30]
plot_contour, gals, "HI", panels='4x1', filename="HI_Sbcm1q1.eps", msg='m1/HI'
plot_contour, gals, "gas", panels='4x1', filename="gas_Sbcm1q1.eps", msg='m1/total'

gals= {frun:"data1/Sbc/Sbcm2_q1", snapnum: 0, xlen: 42.0}
gals= replicate(gals,4)
gals[*].snapnum= [5, 10, 20, 30]
plot_contour, gals, "HI", panels='4x1', filename="HI_Sbcm2q1.eps", msg='m2/HI'
plot_contour, gals, "gas", panels='4x1', filename="gas_Sbcm2q1.eps", msg='m2/total'

gals= {frun:"data1/Sbc/Sbcm3_q1", snapnum: 0, xlen: 42.0}
gals= replicate(gals,4)
gals[*].snapnum= [5, 10, 20, 30]
plot_contour, gals, "HI", panels='4x1', filename="HI_Sbcm3q1.eps", msg='m3/HI'
plot_contour, gals, "gas", panels='4x1', filename="gas_Sbcm3q1.eps", msg='m3/total'












end



