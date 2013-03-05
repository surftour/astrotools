pro do_image, junk


; in contour_makeplot.pro 
;
; make file read in the id lists;
;     hotgas_idlist.txt
;     wind_idlist.txt
;  
; i also turned on the contours.
; maybe we should turn off the pixel
; image? 
; -yes, turned it off
; -commented out timelbl, and added
;   the black hole, and no black hole
;   labels
;


contour_gas, "/raid4/tcox/vc3vc3e_2", 13, 50.0, 'ps', filename='gas.eps', /nolabels, /pubstyle


;
; need to 
;  - change comment to say "no" black hole
;  - change idlist to correct read
;

contour_gas, "/raid4/tcox/vc3vc3e_no", 13, 50.0, 'ps', filename='gas.eps', /nolabels, /pubstyle




end



