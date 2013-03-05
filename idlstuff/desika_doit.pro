
pro doit1, junk


dirlist=["/raid4/tcox/brantsruns/full_model/z6/v4g1q0z6_v4g1q0z6_p", $
	 "/raid4/tcox/brantsruns/full_model/z6/v4g1q1z6_v4g1q1z6_p", $
	 "/raid4/tcox/brantsruns/full_model/z6/v4g2q0z6_v4g2q0z6_p", $
	 "/raid4/tcox/brantsruns/full_model/z6/v4g2q1z6_v4g2q1z6_p", $
	 "/raid4/tcox/brantsruns/full_model/z6/v5g1q0z6_v5g1q0z6_p", $
	 "/raid4/tcox/brantsruns/full_model/z6/v5g1q1z6_v5g1q1z6_p", $
	 "/raid4/tcox/brantsruns/full_model/z6/v5g2q0z6_v5g2q0z6_p", $
	 "/raid4/tcox/brantsruns/full_model/z6/v5g2q1z6_v5g2q1z6_p", $
	 "/raid4/tcox/brantsruns/full_model/z6/v6g1q0z6_v6g1q0z6_p", $
	 "/raid4/tcox/brantsruns/full_model/z6/v6g1q1z6_v6g1q1z6_p", $
	 "/raid4/tcox/brantsruns/full_model/z6/v6g2q0z6_v6g2q0z6_p", $
	 "/raid4/tcox/brantsruns/full_model/z6/v6g2q1z6_v6g2q1z6_p"]

nd= n_elements(dirlist)-1


for i=0, nd do begin

	print, "==============================================================="
	print, "==============================================================="
	print, dirlist[i]

	print, fload_bh_mergertime(dirlist[i])
	find_rdiffs, dirlist[i]                   ; this needs .run time_centers

endfor

end




;====================================================================


pro doit2, junk



; need .run grid_make



desika_grid, "/raid4/tcox/brantsruns/full_model/z6/v4g1q0z6_v4g1q0z6_p", 1, /use_calc_center, /stellar_ages

desika_grid, "/raid4/tcox/brantsruns/full_model/z6/v4g1q1z6_v4g1q1z6_p", 1, /use_calc_center, /stellar_ages

desika_grid, "/raid4/tcox/brantsruns/full_model/z6/v4g2q0z6_v4g2q0z6_p", 1, /use_calc_center, /stellar_ages

desika_grid, "/raid4/tcox/brantsruns/full_model/z6/v4g2q1z6_v4g2q1z6_p", 1, /use_calc_center, /stellar_ages

desika_grid, "/raid4/tcox/brantsruns/full_model/z6/v5g1q0z6_v5g1q0z6_p", 1, /use_calc_center, /stellar_ages

desika_grid, "/raid4/tcox/brantsruns/full_model/z6/v5g1q1z6_v5g1q1z6_p", 1, /use_calc_center, /stellar_ages

desika_grid, "/raid4/tcox/brantsruns/full_model/z6/v5g2q0z6_v5g2q0z6_p", 1, /use_calc_center, /stellar_ages

desika_grid, "/raid4/tcox/brantsruns/full_model/z6/v5g2q1z6_v5g2q1z6_p", 1, /use_calc_center, /stellar_ages

desika_grid, "/raid4/tcox/brantsruns/full_model/z6/v6g1q0z6_v6g1q0z6_p", 1, /use_calc_center, /stellar_ages

desika_grid, "/raid4/tcox/brantsruns/full_model/z6/v6g1q1z6_v6g1q1z6_p", 1, /use_calc_center, /stellar_ages

desika_grid, "/raid4/tcox/brantsruns/full_model/z6/v6g2q0z6_v6g2q0z6_p", 1, /use_calc_center, /stellar_ages

desika_grid, "/raid4/tcox/brantsruns/full_model/z6/v6g2q1z6_v6g2q1z6_p", 1, /use_calc_center, /stellar_ages


end




