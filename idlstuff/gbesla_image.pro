
pro doit, junk




frun= "/n/scratch/hernquist_lab/gbesla/runs/LMC_SMC/LS_Orbit1_L1_8e11Mar5/O1_r90_Azy20"
xlen= 70.0
thiscenter= [0.0, 0.0, 0.0]

;contour_allstars, frun, 420, xlen, center=thiscenter, filename="lmcsmc_420_stars.eps"
;contour_gas, frun, 420, xlen, center=thiscenter, filename="lmcsmc_420_gas.eps"
;contour_dm, frun, 420, xlen, center=thiscenter, filename="lmcsmc_420_dm.eps"



;LMC :  startid = 0 , numpart = 200000 :  DM 100000,  Stars 100000
;SMC :  startid=200000, numpart= 500000 :  Gas 300000, DM 100000, Stars 100000


contour_allstars, frun,   0, xlen, center=thiscenter, filename="lmcsmc_0_starsdm.eps"
contour_allstars, frun, 120, xlen, center=thiscenter, filename="lmcsmc_120_starsdm.eps"
contour_allstars, frun, 320, xlen, center=thiscenter, filename="lmcsmc_320_starsdm.eps"
contour_allstars, frun, 420, xlen, center=thiscenter, filename="lmcsmc_420_starsdm.eps"



end





;==============================================================================





pro do_nature_fig1b, junk

frun= "/n/scratch/hernquist_lab/gbesla/runs/LMC_SMC/LS_Orbit1_L1_8e11Mar5/O1Ssjb_Azy2"
;snapnum= 320
snapnum= 363
thiscenter= [0.0, 0.0, 0.0]
;xlen=130.0
xlen=90.0

contour_gas, frun, snapnum, xlen, center=thiscenter, filename="lmcsmc2_gas.eps"

; now turn on the gas contours in contour_makeplot (above turned on dm contours)

contour_allstars, frun, snapnum, xlen, center=thiscenter, filename="lmcsmc2_starsHIcnt.eps"



end






;==============================================================================






;==============================================================================

