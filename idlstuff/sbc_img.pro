
pro doit, junk

doone, "Sbc201_q1"
doone, "Sbccut201"
doone, "Sbcm1201"
doone, "Sbcm2201"
doone, "Sbcm3201"


end


pro doone, runid

frun="data1/Sbc/"+runid

contour_youngstars, frun, 70, 70.0, filename="img_"+runid+"_young.eps", /nolabels
contour_gas, frun, 70, 70.0, filename="img_"+runid+"_gas.eps", /nolabels
contour_newstars, frun, 70, 70.0, 'ps', filename="img_"+runid+"_new.eps", /nolabels
contour_allstars, frun, 70, 70.0, filename="img_"+runid+"_all.eps", /nolabels


end


