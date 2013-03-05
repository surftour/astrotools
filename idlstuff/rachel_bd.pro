;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Determine D/B As a Function of gas f and M_star
;   -----------------------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro add_one_sim, frun, returnvars=returnvars, $
		bd_ratio=bd_ratio, $
		bt_ratio=bt_ratio, $
		gasf_pre=gasf_pre, $
		gasf_post=gasf_post, $
		fsb=fsb, $
		m_star=m_star


premergtime= 0.3
relaxationtime= 0.4

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

	    ; grab snapnum when system is "pre-merger"
            if fload_time(1) lt (bhmtime-premergtime) then premergersnapnum= i


	    ; grab snapnum when system is "relaxed"
            if (fload_time(1) gt (bhmtime+relaxationtime)) or (i eq totnumsnap) then begin

		; PRE= get gas fraction.
		ok= fload_snapshot_bh(frun,premergersnapnum,/skip_center)
		gasm= total(fload_gas_mass(1))
		stellarmass_pre= total(fload_allstars_mass(1))
		gasf_pre = gasm / (gasm+stellarmass_pre)

		; POST= get gas fraction.
		ok= fload_snapshot_bh(frun,i,/skip_center)
		gasm= total(fload_gas_mass(1))
		stellarmass_post= total(fload_allstars_mass(1))
		gasf_post = gasm / (gasm+stellarmass_post)

		exts='0000'
		exts=exts+strcompress(string(i),/remove_all)
		exts=strmid(exts,strlen(exts)-3,3)
                if not file_test(frun+'/bulgedisk/'+exts+'.txt') then compute_bd, frun, i

                break
            endif
        endfor



;read_bulge_file, frun+'/bulge.txt', m, m_thin, m_other, m_nonrot, $
read_bulge_file, frun, i, m, m_thin, m_other, m_nonrot, $
			g1osm, g1_os_m_thin, g1_os_m_other, g1_os_m_nonrot, $
			g1nsm, g1_ns_m_thin, g1_ns_m_other, g1_ns_m_nonrot, $
			g2osm, g2_os_m_thin, g2_os_m_other, g2_os_m_nonrot, $
			g2nsm, g2_ns_m_thin, g2_ns_m_other, g2_ns_m_nonrot


; now we can get the B/D ratio and M_*
bd_ratio = m_nonrot / (m_thin + m_other)
bt_ratio = m_nonrot / m
m_star = alog10(m*1.0e10/0.7)


print, "bd_ratio= ", bd_ratio
print, "bt_ratio= ", bt_ratio
print, "m_star= ", m_star
print, "gasf_pre= ", gasf_pre
print, "gasf_post= ", gasf_post

print, "stellarmass_pre= ", stellarmass_pre
print, "stellarmass_post= ", stellarmass_post
fsb= (stellarmass_post - stellarmass_pre)/stellarmass_post
print, "fsb= ", fsb


if keyword_set(returnvars) then return

; now plot this up
;=================

;if bd_ratio gt 0.0  then begin & symsize= 4.0 & thiscolor=  50 & endif
;if bd_ratio gt 0.5  then begin & symsize= 3.0 & thiscolor= 100 & endif
;if bd_ratio gt 0.75 then begin & symsize= 2.0 & thiscolor= 150 & endif
;if bd_ratio gt 0.9  then begin & symsize= 1.0 & thiscolor= 200 & endif

if bt_ratio gt 0.0  then begin & symsize= 4.0 & thiscolor=  50 & endif
if bt_ratio gt 0.2  then begin & symsize= 4.0 & thiscolor=  50 & endif
if bt_ratio gt 0.4  then begin & symsize= 3.0 & thiscolor=  70 & endif
if bt_ratio gt 0.5  then begin & symsize= 2.5 & thiscolor=  90 & endif
if bt_ratio gt 0.65 then begin & symsize= 2.2 & thiscolor= 110 & endif
if bt_ratio gt 0.7  then begin & symsize= 2.0 & thiscolor= 130 & endif
if bt_ratio gt 0.8  then begin & symsize= 1.5 & thiscolor= 150 & endif
if bt_ratio gt 0.9  then begin & symsize= 1.2 & thiscolor= 170 & endif
if bt_ratio gt 0.95 then begin & symsize= 1.0 & thiscolor= 190 & endif

symsize= symsize/2.0

;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=1.0, /fill
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=8.0
;oplot, [m_star], [gasf], psym=8, color= thiscolor, thick= 4.0

oplot, [m_star], [gasf], psym=2, color= thiscolor, thick= 8.0


end




;====================================================================================




pro bd_one_sim_time, frun, returnvars=returnvars, $
		time=time, $
		bd_ratio=bd_ratio, $
		bt_ratio=bt_ratio, $
		gasf=gasf, $
		m_star=m_star, $
		pointt= pointt


        spawn, "/bin/ls "+frun+"/snap* | wc ",result
        totnumsnap=long(result[0])-1

        spawn, "/bin/ls "+frun+"/bulgedisk/*.txt | wc ",result
        totbulgedisk=long(result[0])
	time= fltarr(totbulgedisk)
	;gasf= fltarr(totbulgedisk)
	bd_ratio= fltarr(totbulgedisk)
	bt_ratio= fltarr(totbulgedisk)
	m_star= fltarr(totbulgedisk)

        bhmtime= fload_bh_mergertime(frun)
	if bhmtime le 0 then begin
		print, " "
		print, " WARNING: can't determine bhmtime, setting to be 2.0"
		print, " "
		bhmtime= 2.0
	endif

	ii= 0
        for i=0,totnumsnap do begin

	    exts='0000'
	    exts=exts+strcompress(string(i),/remove_all)
	    exts=strmid(exts,strlen(exts)-3,3)
            if file_test(frun+'/bulgedisk/'+exts+'.txt') then begin

		ok= fload_snapshot_bh(frun,i,/header_only)
		time[ii]= fload_time(1)

		; gas fraction.
		;ok= fload_snapshot_bh(frun,i)
		;gasm= total(fload_gas_mass(1))
		;stellarmass= total(fload_allstars_mass(1))
		;gasf[ii] = gasm / (gasm+stellarmass)

		; get bulge/disk info
		read_bulge_file, frun, i, m, m_thin, m_other, m_nonrot, $
			g1osm, g1_os_m_thin, g1_os_m_other, g1_os_m_nonrot, $
			g1nsm, g1_ns_m_thin, g1_ns_m_other, g1_ns_m_nonrot, $
			g2osm, g2_os_m_thin, g2_os_m_other, g2_os_m_nonrot, $
			g2nsm, g2_ns_m_thin, g2_ns_m_other, g2_ns_m_nonrot


		; now we can get the B/D ratio and M_*
		bd_ratio[ii] = m_nonrot / (m_thin + m_other)
		bt_ratio[ii] = m_nonrot / m
		m_star[ii] = alog10(m*1.0e10/0.7)


		ii= ii+1
            endif
        endfor


	;select_thispoint, pointt, thispsym, thiscolor
	;oplot, time, bt_ratio, psym=-thispsym, color= thiscolor, thick= 8.0

	oplot, time, bt_ratio, psym=-3, color= pointt, thick= 8.0


end




;====================================================================================






pro bd_for_one_sim, frun

        spawn, "/bin/ls "+frun+"/bulgedisk/*.txt | wc ",result
        if long(result[0]) gt 1 then begin
		print, " "
		print, " =================================== "
		print, " "
		print, " bulgedisk already calculated"
		print, " "
		;print, " directory: ", frun
		;print, " num of bulgedisk files: ", long(result[0])
		;
		print, " but recalculating anyway"
		print, " "
		print, " =================================== "
		print, " "
		;return
	endif

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
            if (fload_time(1) gt bhmtime) then begin

		exts='0000'
		exts=exts+strcompress(string(i),/remove_all)
		exts=strmid(exts,strlen(exts)-3,3)
                ;if not file_test(frun+'/bulgedisk/'+exts+'.txt') then compute_bd, frun, i
                compute_bd, frun, i

            endif
        endfor

end









;====================================================================================



pro write_text_file, junk

;======================================================================
do_one, "/raid4/tcox/As/A0" & do_one, "/raid4/tcox/As/A0e"
do_one, "/raid4/tcox/As/A1" & do_one, "/raid4/tcox/As/A1e"
do_one, "/raid4/tcox/As/A2" & do_one, "/raid4/tcox/As/A2e"
do_one, "/raid4/tcox/As/A3" & do_one, "/raid4/tcox/As/A3e"
do_one, "/raid4/tcox/As/A4" & do_one, "/raid4/tcox/As/A4e"
do_one, "/raid4/tcox/As/A5" & do_one, "/raid4/tcox/As/A5e"
do_one, "/raid4/tcox/As/A6" & do_one, "/raid4/tcox/As/A6e"
;======================================================================
do_one, "/raid4/tcox/bs/b0e" & do_one, "/raid4/tcox/bs/b0h"
do_one, "/raid4/tcox/bs/b1e" & do_one, "/raid4/tcox/bs/b1h"
do_one, "/raid4/tcox/bs/b2e" & do_one, "/raid4/tcox/bs/b2h"
do_one, "/raid4/tcox/bs/b3e" & do_one, "/raid4/tcox/bs/b3h"
do_one, "/raid4/tcox/bs/b4e" & do_one, "/raid4/tcox/bs/b4h"
do_one, "/raid4/tcox/bs/b5e" & do_one, "/raid4/tcox/bs/b5h"
do_one, "/raid4/tcox/bs/b6e" & do_one, "/raid4/tcox/bs/b6h"
;----
do_one, "/raid4/tcox/bs/b0k" & do_one, "/raid4/tcox/bs/b0f"
do_one, "/raid4/tcox/bs/b1k" & do_one, "/raid4/tcox/bs/b1f"
do_one, "/raid4/tcox/bs/b2k" & do_one, "/raid4/tcox/bs/b2f"
do_one, "/raid4/tcox/bs/b3k" & do_one, "/raid4/tcox/bs/b3f"
do_one, "/raid4/tcox/bs/b4k" & do_one, "/raid4/tcox/bs/b4f"
do_one, "/raid4/tcox/bs/b5k" & do_one, "/raid4/tcox/bs/b5f"
do_one, "/raid4/tcox/bs/b6k" & do_one, "/raid4/tcox/bs/b6f"
;======================================================================
do_one, "/raid4/tcox/ds/d0e"    & do_one, "/raid4/tcox/ds/d0h"
do_one, "/raid4/tcox/ds/d0e2"   & do_one, "/raid4/tcox/ds/d0h2" 
do_one, "/raid4/tcox/ds/d0e2_q" & do_one, "/raid4/tcox/ds/d0h2_q"
do_one, "/raid4/tcox/ds/d0f"    & do_one, "/raid4/tcox/ds/d0k"
do_one, "/raid4/tcox/ds/d0f2_q" & do_one, "/raid4/tcox/ds/d0k2_q"
;---
do_one, "/raid4/tcox/ds/d1e"    & do_one, "/raid4/tcox/ds/d1h"
do_one, "/raid4/tcox/ds/d1e2"   & do_one, "/raid4/tcox/ds/d1h2" 
do_one, "/raid4/tcox/ds/d1e2_q" & do_one, "/raid4/tcox/ds/d1h2_q"
do_one, "/raid4/tcox/ds/d1f"    & do_one, "/raid4/tcox/ds/d1k"
do_one, "/raid4/tcox/ds/d1f2_q" & do_one, "/raid4/tcox/ds/d1k2_q"
do_one, "/raid4/tcox/ds/d1k2_q_i1"
;---
do_one, "/raid4/tcox/ds/d2e"    & do_one, "/raid4/tcox/ds/d2h"
do_one, "/raid4/tcox/ds/d2e2"   & do_one, "/raid4/tcox/ds/d2h2" 
do_one, "/raid4/tcox/ds/d2e2_q" & do_one, "/raid4/tcox/ds/d2h2_q"
do_one, "/raid4/tcox/ds/d2f"    & do_one, "/raid4/tcox/ds/d2k"
do_one, "/raid4/tcox/ds/d2f2_q" & do_one, "/raid4/tcox/ds/d2k2_q"
;---
do_one, "/raid4/tcox/ds/d3e1" & do_one, "/raid4/tcox/ds/d3h1" & do_one, "/raid4/tcox/ds/d3f1" & do_one, "/raid4/tcox/ds/d3k1" 
do_one, "/raid4/tcox/ds/d3e2" & do_one, "/raid4/tcox/ds/d3h2" & do_one, "/raid4/tcox/ds/d3f2" & do_one, "/raid4/tcox/ds/d3k2" 
do_one, "/raid4/tcox/ds/d3e3" & do_one, "/raid4/tcox/ds/d3h3" & do_one, "/raid4/tcox/ds/d3f3" & do_one, "/raid4/tcox/ds/d3k3" 
do_one, "/raid4/tcox/ds/d3e4" & do_one, "/raid4/tcox/ds/d3h4" & do_one, "/raid4/tcox/ds/d3f4" & do_one, "/raid4/tcox/ds/d3k4" 
do_one, "/raid4/tcox/ds/d3e5" & do_one, "/raid4/tcox/ds/d3h5" & do_one, "/raid4/tcox/ds/d3f5" & do_one, "/raid4/tcox/ds/d3k5" 
do_one, "/raid4/tcox/ds/d3e6" & do_one, "/raid4/tcox/ds/d3h6" & do_one, "/raid4/tcox/ds/d3f6" & do_one, "/raid4/tcox/ds/d3k6" 
do_one, "/raid4/tcox/ds/d3e7" & do_one, "/raid4/tcox/ds/d3h7" & do_one, "/raid4/tcox/ds/d3f7" & do_one, "/raid4/tcox/ds/d3k7" 
;---
do_one, "/raid4/tcox/ds/d4e"    & do_one, "/raid4/tcox/ds/d4h"
do_one, "/raid4/tcox/ds/d4e2"   & do_one, "/raid4/tcox/ds/d4h2" 
do_one, "/raid4/tcox/ds/d4e2_q" & do_one, "/raid4/tcox/ds/d4h2_q"
do_one, "/raid4/tcox/ds/d4f"    & do_one, "/raid4/tcox/ds/d4k"
do_one, "/raid4/tcox/ds/d4f2_q" & do_one, "/raid4/tcox/ds/d4k2_q"
;---
do_one, "/raid4/tcox/ds/d5e"    & do_one, "/raid4/tcox/ds/d5h"
do_one, "/raid4/tcox/ds/d5e2"   & do_one, "/raid4/tcox/ds/d5h2" 
do_one, "/raid4/tcox/ds/d5e2_q" & do_one, "/raid4/tcox/ds/d5h2_q"
do_one, "/raid4/tcox/ds/d5f"    & do_one, "/raid4/tcox/ds/d5k"
do_one, "/raid4/tcox/ds/d5f2_q" & do_one, "/raid4/tcox/ds/d5k2_q"
;---
do_one, "/raid4/tcox/ds/d6e"    & do_one, "/raid4/tcox/ds/d6h"
do_one, "/raid4/tcox/ds/d6e2"   & do_one, "/raid4/tcox/ds/d6h2" 
do_one, "/raid4/tcox/ds/d6e2_q" & do_one, "/raid4/tcox/ds/d6h2_q"
do_one, "/raid4/tcox/ds/d6f"    & do_one, "/raid4/tcox/ds/d6k"
do_one, "/raid4/tcox/ds/d6f2_q" & do_one, "/raid4/tcox/ds/d6k2_q"
;---
do_one, "/raid4/tcox/ds/vc3vc3b"
do_one, "/raid4/tcox/ds/vc3vc3c"
do_one, "/raid4/tcox/ds/vc3vc3d"
do_one, "/raid4/tcox/ds/vc3vc3e"
do_one, "/raid4/tcox/ds/vc3vc3f"
do_one, "/raid4/tcox/ds/vc3vc3g"
do_one, "/raid4/tcox/ds/vc3vc3h"
do_one, "/raid4/tcox/ds/vc3vc3i"
do_one, "/raid4/tcox/ds/vc3vc3j"
do_one, "/raid4/tcox/ds/vc3vc3k"
do_one, "/raid4/tcox/ds/vc3vc3l"
do_one, "/raid4/tcox/ds/vc3vc3m"
do_one, "/raid4/tcox/ds/vc3vc3n"
do_one, "/raid4/tcox/ds/vc3vc3o"
do_one, "/raid4/tcox/ds/vc3vc3p"
;======================================================================
do_one, "/raid4/tcox/es/e0e" & do_one, "/raid4/tcox/es/e0h"
do_one, "/raid4/tcox/es/e1e" & do_one, "/raid4/tcox/es/e1h"
do_one, "/raid4/tcox/es/e2e" & do_one, "/raid4/tcox/es/e2h"
do_one, "/raid4/tcox/es/e3e" & do_one, "/raid4/tcox/es/e3h"
do_one, "/raid4/tcox/es/e4e" & do_one, "/raid4/tcox/es/e4h"
do_one, "/raid4/tcox/es/e5e" & do_one, "/raid4/tcox/es/e5h"
do_one, "/raid4/tcox/es/e6e" & do_one, "/raid4/tcox/es/e6h"
;----
do_one, "/raid4/tcox/es/e0k" & do_one, "/raid4/tcox/es/e0f"
do_one, "/raid4/tcox/es/e1k" & do_one, "/raid4/tcox/es/e1f"
do_one, "/raid4/tcox/es/e2k" & do_one, "/raid4/tcox/es/e2f"
do_one, "/raid4/tcox/es/e3k" & do_one, "/raid4/tcox/es/e3f"
do_one, "/raid4/tcox/es/e4k" & do_one, "/raid4/tcox/es/e4f"
do_one, "/raid4/tcox/es/e5k" & do_one, "/raid4/tcox/es/e5f"
do_one, "/raid4/tcox/es/e6k" & do_one, "/raid4/tcox/es/e6f"
;======================================================================
do_one, "/raid4/tcox/zs/z0e" & do_one, "/raid4/tcox/zs/z0h"
do_one, "/raid4/tcox/zs/z1e" & do_one, "/raid4/tcox/zs/z1h"
do_one, "/raid4/tcox/zs/z2e" & do_one, "/raid4/tcox/zs/z2h"
do_one, "/raid4/tcox/zs/z3e" & do_one, "/raid4/tcox/zs/z3h"
do_one, "/raid4/tcox/zs/z4e" & do_one, "/raid4/tcox/zs/z4h"
do_one, "/raid4/tcox/zs/z5e" & do_one, "/raid4/tcox/zs/z5h"
do_one, "/raid4/tcox/zs/z6e" & do_one, "/raid4/tcox/zs/z6h"
;----
do_one, "/raid4/tcox/zs/z0k" & do_one, "/raid4/tcox/zs/z0f"
do_one, "/raid4/tcox/zs/z1k" & do_one, "/raid4/tcox/zs/z1f"
do_one, "/raid4/tcox/zs/z2k" & do_one, "/raid4/tcox/zs/z2f"
do_one, "/raid4/tcox/zs/z3k" & do_one, "/raid4/tcox/zs/z3f"
do_one, "/raid4/tcox/zs/z4k" & do_one, "/raid4/tcox/zs/z4f"
do_one, "/raid4/tcox/zs/z5k" & do_one, "/raid4/tcox/zs/z5f"
do_one, "/raid4/tcox/zs/z6k" & do_one, "/raid4/tcox/zs/z6f"
;======================================================================
do_one, "/raid4/tcox/collisionless/cvc0vc0e" & do_one, "/raid4/tcox/collisionless/cvc0vc0h"
do_one, "/raid4/tcox/collisionless/cvc1vc1e" & do_one, "/raid4/tcox/collisionless/cvc1vc1h"
do_one, "/raid4/tcox/collisionless/cvc2vc2e" & do_one, "/raid4/tcox/collisionless/cvc2vc2h"
do_one, "/raid4/tcox/collisionless/cvc3vc3b"
do_one, "/raid4/tcox/collisionless/cvc3vc3c"
do_one, "/raid4/tcox/collisionless/cvc3vc3d"
do_one, "/raid4/tcox/collisionless/cvc3vc3e"
do_one, "/raid4/tcox/collisionless/cvc3vc3f"
do_one, "/raid4/tcox/collisionless/cvc3vc3g"
do_one, "/raid4/tcox/collisionless/cvc3vc3h"
do_one, "/raid4/tcox/collisionless/cvc3vc3i"
do_one, "/raid4/tcox/collisionless/cvc3vc3j"
do_one, "/raid4/tcox/collisionless/cvc3vc3k"
do_one, "/raid4/tcox/collisionless/cvc3vc3l"
do_one, "/raid4/tcox/collisionless/cvc3vc3m"
do_one, "/raid4/tcox/collisionless/cvc3vc3n"
do_one, "/raid4/tcox/collisionless/cvc3vc3o"
do_one, "/raid4/tcox/collisionless/cvc3vc3p"
;do_one, "/raid4/tcox/collisionless/cvc3h1"
;do_one, "/raid4/tcox/collisionless/cvc3h2"
;do_one, "/raid4/tcox/collisionless/cvc3h3"
; do_one, "/raid4/tcox/collisionless/cvc3h4"
; do_one, "/raid4/tcox/collisionless/cvc3h5"
; do_one, "/raid4/tcox/collisionless/cvc3h6"
;do_one, "/raid4/tcox/collisionless/cvc4vc4e" & do_one, "/raid4/tcox/collisionless/cvc4vc4h"
;do_one, "/raid4/tcox/collisionless/cvc5vc5e" & do_one, "/raid4/tcox/collisionless/cvc5vc5h"
;do_one, "/raid4/tcox/collisionless/cvc6vc6e" & do_one, "/raid4/tcox/collisionless/cvc6vc6h"
;do_one, "/raid4/tcox/collisionless/cvc4nbmme"
;do_one, "/raid4/tcox/collisionless/cvc4nbmmh"
;do_one, "/raid4/tcox/collisionless/cvc5nbmme"
;do_one, "/raid4/tcox/collisionless/cvc5nbmmh"
;do_one, "/raid4/tcox/collisionless/cvc6nbmme"
;do_one, "/raid4/tcox/collisionless/cvc6nbmmh"
;======================================================================



end


pro testit, junk

;do_one, "/raid4/tcox/As/A3"
;do_one, "/raid4/tcox/bs/b3e" & do_one, "/raid4/tcox/bs/b3h"
;do_one, "/raid4/tcox/ds/d6e2_q" & do_one, "/raid4/tcox/ds/d6h2_q"

;do_one, "/raid4/tcox/es/e4e"
;do_one, "/raid4/tcox/es/e4k"

;======================================================================
do_one, "/raid4/tcox/As/A0" & do_one, "/raid4/tcox/As/A0e"
do_one, "/raid4/tcox/As/A1" & do_one, "/raid4/tcox/As/A1e"
do_one, "/raid4/tcox/As/A2" & do_one, "/raid4/tcox/As/A2e"
do_one, "/raid4/tcox/As/A3" & do_one, "/raid4/tcox/As/A3e"
do_one, "/raid4/tcox/As/A4" & do_one, "/raid4/tcox/As/A4e"
do_one, "/raid4/tcox/As/A5" & do_one, "/raid4/tcox/As/A5e"
do_one, "/raid4/tcox/As/A6" & do_one, "/raid4/tcox/As/A6e"


end



;-------------------------
pro do_one, frun

	;add_frun_tolist, frun
	add_frun_and_data_tolist, frun
	;bd_one_sim_time, frun
	;bd_for_one_sim, frun

end



;-------------------------
pro add_frun_tolist, frun

   rfilename = '/home2/tcox/rachellist.txt'

   if not file_test(rfilename) then begin
	openw, 1, rfilename
   endif else begin
	openu, 1, rfilename, /append
   endelse

   printf, 1, frun

   close, 1
end



;-------------------------
pro add_frun_and_data_tolist, frun

    rfilename = '/home2/tcox/rachellist.txt'

    ; check to see if we already have the file information
    spawn, "grep "+frun+" "+rfilename+" | head -n 1",result
    resultarr= strsplit(result, ' ', /extract)
    totret=n_elements(resultarr)
    if totret eq 7 then return

    add_one_sim, frun, /returnvars, $
                bd_ratio=bd_ratio, $
                bt_ratio=bt_ratio, $
                gasf_pre=gasf_pre, $
                gasf_post=gasf_post, $
                fsb=fsb, $
                m_star=m_star


   if not file_test(rfilename) then begin
        openw, 1, rfilename
   endif else begin
        openu, 1, rfilename, /append
   endelse

   tempfrun= strmid(frun+"                             ",0,40)

   ;printf, 1, FORMAT= '(A, "   ", 4(F10.3,"  "))',$
   ;		tempfrun, m_star, gasf, bt_ratio, bd_ratio
   printf, 1, FORMAT= '(A, "   ", 6(F10.3,"  "))',$
		tempfrun, m_star, fsb, gasf_pre, gasf_post, bt_ratio, bd_ratio

   close, 1
end




;----------------------------------

pro read_in_data, mstar_all=mstar_all, fsb_all=fsb_all, $
		  fgas_pre_all=fgas_pre_all, fgas_post_all=fgas_post_all, $
		  bt_all=bt_all, frun_all=frun_all

    rfilename = '/home2/tcox/rachellist.txt'

    ; read to open
    print, "opening: ",rfilename
    rdata= read_ascii(rfilename)

    frun_all= rdata.field1[0,*]
    mstar_all= rdata.field1[1,*]
    ;fgas_all= rdata.field1[2,*]
    ;bt_all= rdata.field1[3,*]
    ;bd_all= rdata.field1[4,*]
    fsb_all= rdata.field1[2,*]
    fgas_pre_all= rdata.field1[3,*]
    fgas_post_all= rdata.field1[4,*]
    bt_all= rdata.field1[5,*]
    bd_all= rdata.field1[6,*]

end


;----------------------------------

pro read_onesim_data, frun, mstar, fsb, fgas_pre, fgas_post, bt, bd

    rfilename = '/home2/tcox/rachellist.txt'

    ; check to see if we already have the file information
    ;spawn, "grep "+frun+" "+rfilename+" | head -n 1",result
    spawn, "grep "+frun+" "+rfilename,result
    if n_elements(result) gt 1 then begin
	print, " "
	print, " "
	print, " WARNING:  grep found two entries in "+rfilename+" but we're only using the first one"
	print, " "
	print, " "
    endif
    resultarr= strsplit(result[0], ' ', /extract)
    totret=n_elements(resultarr)
    if totret ne 7 then return

    frun_full= resultarr[0]
    mstar= float(resultarr[1])
    fsb= float(resultarr[2])
    fgas_pre= float(resultarr[3])
    fgas_post= float(resultarr[4])
    bt= float(resultarr[5])
    bd= float(resultarr[6])

end







;====================================================================================


pro determine_relation, junk

;read_in_data, mstar_all=mstar_all, fgas_all=fgas_all, bt_all=bt_all
read_in_data, mstar_all=mstar_all, fsb_all=fsb_all, $
		  fgas_pre_all=fgas_pre_all, fgas_post_all=fgas_post_all, $
		  bt_all=bt_all

y_tofit= bt_all
weight_tofit= y_tofit*0.0 + 0.07

	idx= where(fgas_pre_all le 0.0)
	if idx(0) ne -1 then fgas_pre_all(idx)= 0.01
	idx= where(fgas_post_all le 0.0)
	if idx(0) ne -1 then fgas_post_all(idx)= 0.01

	;x_tofit= 10^(mstar_all - 10.0)  ; v1
	;x_tofit = 10^(mstar_all - 10.0) * (fgas_all^(-0.3)) ; v2
	x_tofit = fgas_pre_all ; v3
        ; initial guess
        guess= [0.1,1.0]
   
        ; markwardt mpfit procedure
        ;fit_result = MPFITFUN('func_v1', x_tofit, y_tofit, weight_tofit, guess, $
        ;fit_result = MPFITFUN('func_v2', x_tofit, y_tofit, weight_tofit, guess, $
        fit_result = MPFITFUN('func_v3', x_tofit, y_tofit, weight_tofit, guess, $
                                BESTNORM=bestnorm, DOF=dof)


	redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2

end



; trial 1, simplest one
function func_v1, x, p

        ; p[0] = Normalization, i.e., the B/T(0)
        ; p[1] = m_star scaling, beta
        return, p[0]*(x^p[1])

end

function func_v2, x, p

        ; p[0] = Normalization, i.e., the B/T(0)
        ; p[1] = x is now a combination of fgas and mstar
        return, p[0]*(x^p[1])

end

function func_v3, x, p

        ; p[0] = Intercept, i.e., B/T at f_gas = 0
        ; p[1] = Slope
        return, p[0] + x*p[1]

end




;====================================================================================
;====================================================================================




pro bt_mstar, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "bt_mstar, junk"
   print, "  "
   return
endif


;----------------------------------------------------

read_in_data, mstar_all=mstar_all, fsb_all=fsb_all, $
		  fgas_pre_all=fgas_pre_all, fgas_post_all=fgas_post_all, $
		  bt_all=bt_all

xvar= mstar_all
print, min(xvar), max(xvar)

yvar= bt_all

;-----------------------------------------------------



filename='bt_mstar.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


; -----------------------------
; set up constants
; -----------------------------

xaxistitle= "!6Log Stellar Mass (M!D!9n!6!N)"
xmax = 12.8
xmin =  9.0

yaxistitle= '!6B/T'
ymax = 1.0
ymin = 0.0


; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
	color= 0, xcharsize=1.50, ycharsize=1.50, $
	xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle



oplot, xvar, yvar, psym=2, symsize= 1.0, color= 150, thick= 3.0



; show fit
; -------------
; v1
; ----
; bad b/t ratio
;
;;pfit= [0.525814, 0.0788186]   ; first round, smaller subset of data
;pfit= [0.513685, 0.0819290]   
;x= [xmin,xmax]
;xnotlog= 10^(x-10.0)
;y= func_v1(xnotlog, pfit)
;oplot, x, y, psym=-3, color= 0, linestyle= 2, thick= 3.0
;xyouts, 0.35, 0.20, '0.514 (M!D*!N/10!E10!N M!D!9n!6!N)!E0.082!N', size=2.0, color=0, /normal, charthick=4.0

pfit= [0.55765, 0.00485162]   
x= [xmin,xmax]
xnotlog= 10^(x-10.0)
y= func_v1(xnotlog, pfit)
oplot, x, y, psym=-3, color= 0, linestyle= 1, thick= 5.0
xyouts, 0.35, 0.20, '0.56 (M!D*!N/10!E10!N M!D!9n!6!N)!E0.004!N', size=2.0, color=0, /normal, charthick=4.0



; done
; -------------
device, /close


end








;====================================================================================




pro bt_fgas, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "bt_fgas, junk"
   print, "  "
   return
endif



filename='bt_fgas.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;----------------------------------------------------

read_in_data, mstar_all=mstar_all, fsb_all=fsb_all, $
		  fgas_pre_all=fgas_pre_all, fgas_post_all=fgas_post_all, $
		  bt_all=bt_all, frun_all=frun_all


yvar= bt_all
print, min(yvar), max(yvar)

yaxistitle= '!6B/T'
ymax = 1.0
ymin = 0.0



xvar= fgas_pre_all & xaxistitle= "!18f!6!Dgas!N (pre starburst)"
;xvar= fgas_post_all & xaxistitle= "!18f!6!Dgas!N (post starburst)"
;xvar= fsb_all & xaxistitle= "!18f!6!Dsb!N "
print, min(xvar), max(xvar)

xmax =  0.9
xmin =  0.0



; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
	color= 0, xcharsize=1.50, ycharsize=1.50, $
	xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle



oplot, xvar, yvar, psym=2, symsize= 1.0, color= 150, thick= 3.0



; show fit
; -------------

; v3
; ----
;pfit= [0.675566, -0.148056]  ; first, smaller amount of data, bad b/t
;pfit= [0.659938, -0.152539]  ; bad b/t
pfit= [0.594414, -0.160330]
x= [xmin,xmax]
y= pfit[0] + pfit[1]*x
oplot, x, y, psym=-3, color= 0, linestyle= 2, thick= 3.0
xyouts, 0.37, 0.20, '0.59 - 0.16!18f!6!Dgas!N', size=2.0, color=0, /normal, charthick=4.0


; done
; -------------
device, /close


end








;====================================================================================





pro bt_combo, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "bt_combo, junk"
   print, "  "
   return
endif


;------------------------------------------------------------------

;read_in_data, mstar_all=mstar_all, fgas_all=fgas_all, bt_all=bt_all
read_in_data, mstar_all=mstar_all, fsb_all=fsb_all, $
		  fgas_pre_all=fgas_pre_all, fgas_post_all=fgas_post_all, $
		  bt_all=bt_all

;idx= where(fgas_all le 0.0)
;if idx(0) ne -1 then fgas_all(idx)= 0.01

;xvar= (mstar_all-10.0) * fgas_all
xvar= (mstar_all-10.0) * (fgas_all^(0.2))
;xvar= alog10((10^(mstar_all-10.0)^(0.076)) * (fgas_all^(-0.023)))
print, min(xvar), max(xvar)

yvar= bt_all

;------------------------------------------------------------------



filename='bt_combo.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


; -----------------------------
; set up constants
; -----------------------------

xaxistitle= "!6Log (M!D*!N!E0.0786!N  !18f!6!Dgas!N!E-0.024!N) "
xmax =  0.3
xmin = -0.1

yaxistitle= '!6B/T'
ymax = 1.0
ymin = 0.0


; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
	color= 0, xcharsize=1.50, ycharsize=1.50, $
	xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle



oplot, xvar, yvar, psym=2, symsize= 1.0, color= 150, thick= 3.0



; show fit
; -------------

; v2
; ----
;pfit= [0.505061, 0.0756883]   ; first round, small subset of data
;pfit= [0.491373, 0.0756883]
;x= [xmin,xmax]
;y= pfit[0] + x   ; slope, pfit[1], is taken care of in x-axis scaling
;oplot, x, y, psym=-3, color= 0, linestyle= 2, thick= 3.0
;xyouts, 0.24, 0.20, '0.491 (M!D*!N/10!E10!N M!D!9n!6!N)!E0.0786!N !18f!6!Dgas!N!E-0.024!N', size=1.7, color=0, /normal, charthick=4.0


; done
; -------------
device, /close


end







;====================================================================================





pro bt_time, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "bt_time, junk"
   print, "  "
   return
endif


filename='bt_time.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


; -----------------------------
; set up constants
; -----------------------------

yaxistitle= '!6B/T'
ymax =  1.0
ymin = 0.0

xaxistitle= '!6Time (Gyr/h)'
xmax = 3.1
xmin = 1.0


; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
	color= 0, xcharsize=1.50, ycharsize=1.50, $
	xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle



;------------------------------------------------

; no real trend with time

;bd_one_sim_time, "/raid4/tcox/es/e1e", pointt= 50
;bd_one_sim_time, "/raid4/tcox/es/e1h", pointt= 50
;bd_one_sim_time, "/raid4/tcox/es/e1k", pointt= 50
;bd_one_sim_time, "/raid4/tcox/es/e1f", pointt= 50

;bd_one_sim_time, "/raid4/tcox/es/e3e", pointt= 100
;bd_one_sim_time, "/raid4/tcox/es/e3h", pointt= 100
;bd_one_sim_time, "/raid4/tcox/es/e3k", pointt= 100
;bd_one_sim_time, "/raid4/tcox/es/e3f", pointt= 100

;bd_one_sim_time, "/raid4/tcox/es/e5e", pointt= 150
;bd_one_sim_time, "/raid4/tcox/es/e5h", pointt= 150
;bd_one_sim_time, "/raid4/tcox/es/e5k", pointt= 150
;bd_one_sim_time, "/raid4/tcox/es/e5f", pointt= 150

;------------------------------------------------

; no real difference with integration

;bd_one_sim_time, "/raid4/tcox/es/e2k", pointt= 150
;bd_one_sim_time, "/raid4/tcox/es/e2k_i1", pointt= 150
;bd_one_sim_time, "/raid4/tcox/es/e3k", pointt= 50
;bd_one_sim_time, "/raid4/tcox/es/e3k_i1", pointt= 50
;bd_one_sim_time, "/raid4/tcox/es/e4k", pointt= 200
;bd_one_sim_time, "/raid4/tcox/es/e4k_i1", pointt= 200

;bd_one_sim_time, "/raid4/tcox/es/e2e", pointt= 150
;bd_one_sim_time, "/raid4/tcox/es/e2e_i1", pointt= 150
;bd_one_sim_time, "/raid4/tcox/es/e3e", pointt= 50
;bd_one_sim_time, "/raid4/tcox/es/e3e_i1", pointt= 50
;bd_one_sim_time, "/raid4/tcox/es/e4e", pointt= 200
;bd_one_sim_time, "/raid4/tcox/es/e4e_i1", pointt= 200

;------------------------------------------------

bd_one_sim_time, "/raid4/tcox/es/e5k", pointt= 150
bd_one_sim_time, "/raid4/tcox/es/e5k_badbd", pointt= 100

;------------------------------------------------


; done
; -------------
device, /close


end




;====================================================================================




pro bt_mstar_fgas, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "bt_mstar_fgas, junk"
   print, "  "
   return
endif


filename='bt_mstar_fgas.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;----------------------------------------------------

read_in_data, mstar_all=mstar_all, fsb_all=fsb_all, $
		  fgas_pre_all=fgas_pre_all, fgas_post_all=fgas_post_all, $
		  bt_all=bt_all, frun_all=frun_all


xvar= mstar_all
print, min(xvar), max(xvar)

xaxistitle= "!6Log Stellar Mass (M!D!9n!6!N)"
xmax = 12.8
xmin =  9.0

yvar= fgas_pre_all & yaxistitle= "!18f!6!Dgas!N (pre starburst)"
;yvar= fgas_post_all & yaxistitle= "!18f!6!Dgas!N (post starburst)"
;yvar= fsb_all & yaxistitle= "!18f!6!Dsb!N "
print, min(yvar), max(yvar)

ymax =  0.9
ymin =  0.0


;idx= where(yvar gt 0.6)
;print, yvar(idx)
;print, xvar(idx)
;print, frun_all(idx)


; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
	color= 0, xcharsize=1.50, ycharsize=1.50, $
	xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle



oplot, xvar, yvar, psym=2, symsize= 1.0, color= 150, thick= 3.0





; done
; -------------
device, /close


end




;====================================================================================





pro bt_mstar_select, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "bt_mstar_select, junk"
   print, "  "
   return
endif


filename='bt_mstar_select.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


; -----------------------------
; set up constants
; -----------------------------

xaxistitle= "!6Log Stellar Mass (M!D!9n!6!N)"
xmax = 12.8
xmin =  9.0

yaxistitle= '!6B/T'
ymax = 1.0
ymin = 0.0


; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
	color= 0, xcharsize=1.50, ycharsize=1.50, $
	xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
	xtitle=xaxistitle, ytitle=yaxistitle




process_one_simgroup, ["b0e", "b1e", "b2e", "b3e", "b4e", "b5e", "b6e"], 100, 3
process_one_simgroup, ["e0e", "e1e", "e2e", "e3e", "e4e", "e5e", "e6e"], 110, 2
process_one_simgroup, ["z0e", "z1e", "z2e", "z3e", "z4e", "z5e", "z6e"], 120, 5

process_one_simgroup, ["b0f", "b1f", "b2f", "b3f", "b4f", "b5f", "b6f"], 150, 3
process_one_simgroup, ["e0f", "e1f", "e2f", "e3f", "e4f", "e5f", "e6f"], 160, 2
process_one_simgroup, ["z0f", "z1f", "z2f", "z3f", "z4f", "z5f", "z6f"], 170, 5

;process_one_simgroup, ["b0h", "b1h", "b2h", "b3h", "b4h", "b5h", "b6h"], 40, 3
;process_one_simgroup, ["e0h", "e1h", "e2h", "e3h", "e4h", "e5h", "e6h"], 50, 2
;process_one_simgroup, ["z0h", "z1h", "z2h", "z3h", "z4h", "z5h", "z6h"], 60, 5

;process_one_simgroup, ["b0k", "b1k", "b2k", "b3k", "b4k", "b5k", "b6k"], 200, 3
;process_one_simgroup, ["e0k", "e1k", "e2k", "e3k", "e4k", "e5k", "e6k"], 210, 2
;process_one_simgroup, ["z0k", "z1k", "z2k", "z3k", "z4k", "z5k", "z6k"], 220, 5



; a little key
xx1= 9.2
xx2= 9.5
oplot, [xx1,xx2], [0.9, 0.9], psym= -3, color= 0, thick=3.0
oplot, [xx1 + 0.5*(xx2-xx1)], [0.9], psym=3, color= 0, thick= 3.0
xyouts, [xx2+0.2], [0.885], '80% gas fraction', size=1.3, color= 0, charthick=3.0, /data

oplot, [xx1,xx2], [0.85, 0.85], psym= -3, color= 0, thick=3.0
oplot, [xx1 + 0.5*(xx2-xx1)], [0.85], psym=2, color= 0, thick= 3.0
xyouts, [xx2+0.2], [0.835], '20% ', size=1.3, color= 0, charthick=3.0, /data

oplot, [xx1,xx2], [0.8, 0.8], psym= -3, color= 0, thick=3.0
oplot, [xx1 + 0.5*(xx2-xx1)], [0.8], psym=5, color= 0, thick= 3.0
xyouts, [xx2+0.2], [0.785], '5% ', size=1.3, color= 0, charthick=3.0, /data


xyouts, 0.35, 0.62, 'e-orbit', size=1.3, color= 100, charthick=3.0, /normal
xyouts, 0.35, 0.32, 'f-orbit', size=1.3, color= 150, charthick=3.0, /normal


; done
; -------------
device, /close


end



;--------------------------------------


pro process_one_simgroup, fruns, pptc, ppts

n_i= n_elements(fruns)
massarr= fltarr(n_i)
btarr= fltarr(n_i)

for i=0,n_i-1 do begin
	read_onesim_data, fruns[i], mstar, fsb, fgas_pre, fgas_post, bt, bd
	massarr[i]= mstar
	btarr[i]= bt
endfor

oplot, massarr, btarr, psym=-ppts, symsize= 1.0, color= pptc, thick= 3.0

end





;====================================================================================


