
pro dosfr, junk

; .run sfr_multi
;sfr_test, "data1/ds/vc3vc3c_1", filename="sfr_vc3vc3c_1.eps"
;sfr_test, "data1/ds/vc3vc3c_2", filename="sfr_vc3vc3c_2.eps"
;sfr_test, "data1/ds/vc3vc3c_3", filename="sfr_vc3vc3c_3.eps"
;sfr_test, "data1/ds/vc3vc3c_4", filename="sfr_vc3vc3c_4.eps"
;sfr_test, "data1/ds/vc3vc3c_5", filename="sfr_vc3vc3c_5.eps"
sfr_test, "data1/ds/vc3vc3c_6", filename="sfr_vc3vc3c_6.eps"
sfr_test, "data1/ds/vc3vc3c_7", filename="sfr_vc3vc3c_7.eps"

;sfr_test, "data/ds/vc3vc3c", filename="sfr_vc3vc3c.eps"
;sfr_test, "data/z3/b5e", filename="sfr_b5e.eps"

end


;
;
; ----------------------------------------------------------
;
;


pro doit, junk


frun= "data/ds/vc3vc3c"

grab_young_stars, frun, 8
grab_young_stars, frun, 9



frun= "data/z3/b5e"

grab_young_stars, frun, 38
grab_young_stars, frun, 39


end

;
;
; ----------------------------------------------------------
;
;


pro grab_young_stars, frun, snapnum

ok=fload_snapshot_bh(frun,snapnum)

id=fload_newstars_id(1)
time=fload_time(1)
age=fload_newstars_age(1)
m=fload_newstars_mass(1)
help, m, age, id, time
print, min(age), max(age)


age= (time-age)/0.7
print, min(age), max(age)


mappingsparts_idx= where(age lt 0.01)
mappingsparts_age= age(mappingsparts_idx)
help, mappingsparts_age

end



