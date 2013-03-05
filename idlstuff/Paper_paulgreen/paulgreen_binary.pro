
pro doit, junk


;
;  .run sfr_multi
;
;sfr_one, "ds/d3h1", filename= "sfr_d3h1.eps", show_specific_times=[0.9/0.7], lcolor= 50
;sfr_one, "ds/d3h2", filename= "sfr_d3h2.eps", show_specific_times=[1.0/0.7, 1.1/0.7], lcolor= 50
;sfr_one, "ds/d3h2", filename= "sfr_d3h2a.eps", show_specific_times=[1.0/0.7], lcolor= 50
;sfr_one, "ds/d3h2", filename= "sfr_d3h2b.eps", show_specific_times=[1.1/0.7], lcolor= 50
;sfr_one, "ds/d3h5", filename= "sfr_d3h5.eps", show_specific_times=[3.6/0.7, 3.7/0.7], lcolor= 50, xmax= 5.0/0.7
;sfr_one, "ds/d3h5", filename= "sfr_d3h5a.eps", show_specific_times=[3.6/0.7], lcolor= 50, xmax= 5.0/0.7
;sfr_one, "ds/d3h5", filename= "sfr_d3h5b.eps", show_specific_times=[3.7/0.7], lcolor= 50, xmax= 5.0/0.7
;sfr_one, "ds/d3h6", filename= "sfr_d3h6.eps", show_specific_times=[2.7/0.7], lcolor= 50, xmax= 5.0/0.7
;sfr_one, "ds/d3h7", filename= "sfr_d3h7.eps", show_specific_times=[1.0/0.7], lcolor= 50


;contour_allstars, "data/ds/d3h1",  9, 70.0, filename="cnt_d3h1.eps", /pubstyle, /nolabels, /show_rsep
;contour_allstars, "data/ds/d3h2", 10, 70.0, filename="cnt_d3h2a.eps", /pubstyle, /nolabels, /show_rsep
;contour_allstars, "data/ds/d3h2", 11, 70.0, filename="cnt_d3h2b.eps", /pubstyle, /nolabels, /show_rsep
;contour_allstars, "data/ds/d3h5", 36, 70.0, filename="cnt_d3h5a.eps", /pubstyle, /nolabels, /show_rsep
;contour_allstars, "data/ds/d3h5", 37, 70.0, filename="cnt_d3h5b.eps", /pubstyle, /nolabels, /show_rsep
;contour_allstars, "data/ds/d3h6", 27, 70.0, filename="cnt_d3h6.eps", /pubstyle, /nolabels, /show_rsep
;contour_allstars, "data/ds/d3h7", 10, 70.0, filename="cnt_d3h7.eps", /pubstyle, /nolabels, /show_rsep


;
; .run movie_image_plus_plot
; .run bh_multi
; .run sfr_multi
; .run bh_details
; .run time_centers
;
img4_3, "data/ds/d3h1", 8, filename="d3h1_8.eps"
img4_3, "data/ds/d3h1", 9, filename="d3h1_9.eps"

img4_3, "data/ds/d3h2", 10, filename="d3h2_10.eps"
img4_3, "data/ds/d3h2", 11, filename="d3h2_11.eps"

img4_3, "data/ds/d3h5", 35, filename="d3h5_35.eps"
img4_3, "data/ds/d3h5", 36, filename="d3h5_36.eps"
img4_3, "data/ds/d3h5", 37, filename="d3h5_37.eps"
img4_3, "data/ds/d3h5", 38, filename="d3h5_38.eps"

img4_3, "data/ds/d3h6", 26, filename="d3h6_26.eps"
img4_3, "data/ds/d3h6", 27, filename="d3h6_27.eps"
img4_3, "data/ds/d3h6", 28, filename="d3h6_28.eps"

img4_3, "data/ds/d3h7",  8, filename="d3h7_8.eps"
img4_3, "data/ds/d3h7",  9, filename="d3h7_9.eps"
img4_3, "data/ds/d3h7", 10, filename="d3h7_10.eps"


end


