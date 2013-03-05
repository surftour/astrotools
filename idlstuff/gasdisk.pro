
pro doit, junk

do_gasdisk, "/raid4/tcox/As/A3e"
do_gasdisk, "/raid4/tcox/As/A3"

do_gasdisk, "/raid4/tcox/bs/b3e"
do_gasdisk, "/raid4/tcox/bs/b3h"
do_gasdisk, "/raid4/tcox/bs/b3f"
do_gasdisk, "/raid4/tcox/bs/b3k"

do_gasdisk, "/raid4/tcox/cs/c3e"
do_gasdisk, "/raid4/tcox/cs/c3h"
do_gasdisk, "/raid4/tcox/cs/c3f"
do_gasdisk, "/raid4/tcox/cs/c3k"

do_gasdisk, "/raid4/tcox/ds/d3e7"
do_gasdisk, "/raid4/tcox/ds/d3h7"
do_gasdisk, "/raid4/tcox/ds/d3f7"
do_gasdisk, "/raid4/tcox/ds/d3k7"

do_gasdisk, "/raid4/tcox/ds/vc3vc3b"
do_gasdisk, "/raid4/tcox/ds/vc3vc3c"
do_gasdisk, "/raid4/tcox/ds/vc3vc3d"
do_gasdisk, "/raid4/tcox/ds/vc3vc3e"
do_gasdisk, "/raid4/tcox/ds/vc3vc3f"
do_gasdisk, "/raid4/tcox/ds/vc3vc3g"
do_gasdisk, "/raid4/tcox/ds/vc3vc3h"
do_gasdisk, "/raid4/tcox/ds/vc3vc3i"
do_gasdisk, "/raid4/tcox/ds/vc3vc3j"
do_gasdisk, "/raid4/tcox/ds/vc3vc3k"
do_gasdisk, "/raid4/tcox/ds/vc3vc3l"
do_gasdisk, "/raid4/tcox/ds/vc3vc3m"
do_gasdisk, "/raid4/tcox/ds/vc3vc3n"
do_gasdisk, "/raid4/tcox/ds/vc3vc3o"
do_gasdisk, "/raid4/tcox/ds/vc3vc3p"

do_gasdisk, "/raid4/tcox/es/e3e"
do_gasdisk, "/raid4/tcox/es/e3h"
do_gasdisk, "/raid4/tcox/es/e3f"
do_gasdisk, "/raid4/tcox/es/e3k"

do_gasdisk, "/raid4/tcox/zs/z3e"
do_gasdisk, "/raid4/tcox/zs/z3h"
do_gasdisk, "/raid4/tcox/zs/z3f"
do_gasdisk, "/raid4/tcox/zs/z3k"

end




;
;
;   find image when remnant is relaxed
; --------------------------------------
pro do_gasdisk, frun

	premergtime= 0.3
	;relaxationtime= 0.4
	relaxationtime= 0.8

        spawn, "/bin/ls "+frun+"/snap* | wc ",result
        totnumsnap=long(result[0])-1
        snapnum=long(result[0])-4

        bhmtime= fload_bh_mergertime(frun)
        if bhmtime le 0 then begin
                print, " "
                print, " WARNING: can't determine bhmtime, setting to be 2.0"
                print, " "
                bhmtime= 2.0
        endif

        for i=0,totnumsnap do begin
            ok= fload_snapshot_bh(frun,i,/header_only)

            ; grab snapnum when system is "relaxed"
            if fload_time(1) lt (bhmtime-premergtime) then premergersnapnum= i


            ; grab snapnum when system is "relaxed"
            if (fload_time(1) gt (bhmtime+relaxationtime)) or (i eq totnumsnap) then begin

		postmergersnapnum= i

		contour_gas, frun, postmergersnapnum, 10.0, 'ps', filename=frun+'/gas.eps', /use_calc_center
		contour_gas, frun, postmergersnapnum, 10.0, 'ps', filename=frun+'/gasxz.eps', /use_calc_center, /xz

                break
            endif
        endfor


end


