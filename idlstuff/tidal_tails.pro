
pro doit, junk
;
; with velocity field turned on
;
;contour_allstars, "data1/tides/disk_pro",  0, 30.0, 'ps', filename="test_pro.eps"
;contour_allstars, "data1/tides/disk_pro1", 0, 30.0, 'ps', filename="test_pro1.eps"
;contour_allstars, "data1/tides/disk_pro2", 0, 30.0, 'ps', filename="test_pro2.eps"
;contour_allstars, "data1/tides/disk_pro3", 0, 30.0, 'ps', filename="test_pro3.eps"
;contour_allstars, "data1/tides/disk_pro4", 0, 30.0, 'ps', filename="test_pro4.eps"
;contour_allstars, "data1/tides/disk_pro5", 0, 30.0, 'ps', filename="test_pro5.eps"
;contour_allstars, "data1/tides/disk_pro6", 0, 30.0, 'ps', filename="test_pro6.eps"


;contour_multi_3x4, "data1/tides/disk_pro",  filename="seq_pro.eps",  snaparr=[0,1,2,3,4,6,8,10,12,14,16,20]
;contour_multi_3x4, "data1/tides/disk_pro1", filename="seq_pro1.eps", snaparr=[0,1,2,3,4,6,8,10,12,14,16,20]
;contour_multi_3x4, "data1/tides/disk_pro2", filename="seq_pro2.eps", snaparr=[0,1,2,3,4,6,8,10,12,14,16,20]
;contour_multi_3x4, "data1/tides/disk_pro3", filename="seq_pro3.eps", snaparr=[0,1,2,3,4,6,8,10,12,14,16,20]
;contour_multi_3x4, "data1/tides/disk_pro4", filename="seq_pro4.eps", snaparr=[0,1,2,3,4,6,8,10,12,14,16,20], xxlen=70.0
;contour_multi_3x4, "data1/tides/disk_pro5", filename="seq_pro5.eps", snaparr=[0,1,2,3,4,6,8,10,12,14,16,20], xxlen=70.0
;contour_multi_3x4, "data1/tides/disk_pro6", filename="seq_pro6.eps", snaparr=[0,1,2,3,4,6,8,10,12,14,16,20], xxlen=70.0


;contour_multi_3x4, "data1/tides/4x_pro", filename="data1/tides/4x_pro/panel_stars_perturb.eps", snaparr=[49,50,51,52,53,55,57,59,61,63,65,69], xxlen=60.0



end





pro check_tail_material, junk

; .run id_tidalmaterial
; do_the_following, 1
;generate_idlist, "data1/tides/4x_pro", 55
;generate_idlist, "data1/tides/4x_pro", 61  ; a bit messy to pick out particles



;
; with velocity field turned off and idoverplotting turned on
;
;contour_allstars, "data1/tides/4x_pro", 22, 70.0, 'ps', filename="data1/tides/4x_pro/stars_022_arm1.eps"
;contour_allstars, "data1/tides/4x_pro", 49, 70.0, 'ps', filename="data1/tides/4x_pro/stars_049_arm1.eps"
;contour_allstars, "data1/tides/4x_pro", 55, 70.0, 'ps', filename="data1/tides/4x_pro/stars_055_arm1.eps"
;contour_allstars, "data1/tides/4x_pro", 58, 70.0, 'ps', filename="data1/tides/4x_pro/stars_058_arm1.eps"
;contour_allstars, "data1/tides/4x_pro", 61, 70.0, 'ps', filename="data1/tides/4x_pro/stars_061_arm1.eps"
;contour_allstars, "data1/tides/4x_pro", 71, 70.0, 'ps', filename="data1/tides/4x_pro/stars_071_arm1.eps"
contour_allstars, "data1/tides/4x_pro", 22, 70.0, 'ps', filename="data1/tides/4x_pro/stars_022_arm2.eps"
contour_allstars, "data1/tides/4x_pro", 49, 70.0, 'ps', filename="data1/tides/4x_pro/stars_049_arm2.eps"
contour_allstars, "data1/tides/4x_pro", 55, 70.0, 'ps', filename="data1/tides/4x_pro/stars_055_arm2.eps"
contour_allstars, "data1/tides/4x_pro", 58, 70.0, 'ps', filename="data1/tides/4x_pro/stars_058_arm2.eps"
contour_allstars, "data1/tides/4x_pro", 61, 70.0, 'ps', filename="data1/tides/4x_pro/stars_061_arm2.eps"
contour_allstars, "data1/tides/4x_pro", 71, 70.0, 'ps', filename="data1/tides/4x_pro/stars_071_arm2.eps"


;
; with idoverplotting turned on
;
;contour_multi_3x4, "data1/tides/4x_pro", filename="data1/tides/4x_pro/stars_seq_tides.eps"


end


