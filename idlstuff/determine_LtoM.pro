;------------------------------------------------------------------------
;
;    Random Procedures related to ULIRG L/M (Sukanya' quantity)
;
;
;
;------------------------------------------------------------------------




function grab_m, frun, snapnum, radius=radius

; load snapshot
ok=fload_snapshot_bh(frun,snapnum)



; -----------------
center= [0,0,0]
center_bh= fload_center_alreadycomp(1)

; two black holes
if fload_npart(5) gt 1 then begin
        bhid= fload_blackhole_id(1)
        bhid1= bhid[0]
        bhid2= bhid[1]
        print, "Blackhole ID's: ", bhid1, bhid2
        center1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
        center2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
        center_bh= center1
        ;center_bh= center2
endif

; one black hole
if fload_npart(5) eq 1 then begin
        bhid= fload_blackhole_id(1)
        center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)
endif



; ----------------
print, "using center= ", center_bh
x= fload_gas_xyz('x',center=center_bh)
y= fload_gas_xyz('y',center=center_bh)
z= fload_gas_xyz('z',center=center_bh)
rr= sqrt(x*x + y*y + z*z)

m= fload_gas_mass(1)
ids= fload_gas_id(1)

print, "Total Gas Mass= ", total(m)

coldf= fload_gas_coldfraction(1)
print, "coldf max/min: ", max(coldf), min(coldf)

ColdGasMass= m*coldf
totalcoldmass= total(ColdGasMass)
print, "Cold Gas Mass= ", totalcoldmass


sortr= sort(rr)
rr= rr(sortr)
m= m(sortr)
ColdGasMass= ColdGasMass(sortr)


;
; find radius which encloses 
; this mass fraction
enc_mf= 0.8
print, "Desired Enclosed Mass= ", enc_mf * totalcoldmass

n= n_elements(ColdGasMass)
n_idx= long(0.7 * enc_mf * n)
print, "n/n_idx= ",n, n_idx

;i= n_idx
i= 2000L
repeat begin
	enclosed_mass= total(ColdGasMass[0:i])
	i= i+1
endrep until ((i ge n) or (enclosed_mass gt enc_mf*totalcoldmass))
print, "i= ", i, enclosed_mass, enc_mf*totalcoldmass

radius= rr(i-2)
print, " Radius = ", radius
print, " Total gas mass within sphere = ", total(m[0:i-2])
print, " Cold gas mass within Sphere = ", total(ColdGasMass[0:i-2])


return, total(ColdGasMass[0:i-2])


end




pro do_just_one, frun, snapnum= snapnum

m_time= fload_bh_mergertime(frun)
print, "merger time= ", m_time, " Gyr/h"
if not keyword_set(snapnum) then snapnum= 1 + fix(10.0*m_time)
;if not keyword_set(snapnum) then snapnum= 2 + fix(10.0*m_time)    ; one snap later
;if not keyword_set(snapnum) then snapnum= 3 + fix(10.0*m_time)    ; two snaps later

; gas mass
M_gas= grab_m(frun,snapnum,radius=radius)
M_gas= M_gas * 1.0e+10
print, "M_gas (within ",radius," kpc/h)  = ", M_gas

; luminosity
read_bololum_file, frun, time, sbololum, bhbololum
L_bh_bolo= 10^(bhbololum[snapnum])
L_s_bolo= 10^(sbololum[snapnum])
L_bolo= L_bh_bolo + L_s_bolo
print, "Bolometric Luminosity= ", L_bolo
print, "                (BH) = ", L_bh_bolo
print, "             (Stars) = ", L_s_bolo

; L/M ratio
print, "xxxx ", frun, " ", fload_timelbl(1,3,/noteq), "  (", L_bolo, " / ", M_gas, " )  = ", L_bolo / M_gas


end




;==========================================================================




pro ltom, junk


;do_just_one, "/raid4/tcox/vc3vc3e_2"
; -------------------------------------------------
do_just_one, "/raid4/tcox/bs/b3e"
do_just_one, "/raid4/tcox/bs/b3h"
do_just_one, "/raid4/tcox/bs/b3f"
do_just_one, "/raid4/tcox/bs/b3k"
; -------------------------------------------------
do_just_one, "/raid4/tcox/bs/b4e"
do_just_one, "/raid4/tcox/bs/b4h"
do_just_one, "/raid4/tcox/bs/b4f"
do_just_one, "/raid4/tcox/bs/b4k"
; -------------------------------------------------
do_just_one, "/raid4/tcox/bs/b5e"
do_just_one, "/raid4/tcox/bs/b5h"
do_just_one, "/raid4/tcox/bs/b5f"
do_just_one, "/raid4/tcox/bs/b5k"
; -------------------------------------------------
;do_just_one, "/raid4/tcox/ds/d4e2_q"
;do_just_one, "/raid4/tcox/ds/d4h2_q"
;do_just_one, "/raid4/tcox/ds/d5e2_q"
;

end





;==========================================================================




pro ltom_time, junk


;frun= "/raid4/tcox/vc3vc3e_2" & start_snap= 42
;frun= "/raid4/tcox/vc3vc3h_2" & start_snap= 27
;frun= "/raid4/tcox/vc3vc3h_bulge" & start_snap= 29
;frun= "/raid4/tcox/vc3vc3e_bulge" & start_snap= 29
;frun= "/raid4/tcox/vc3vc3f_bulge" & start_snap= 52
;frun= "/raid4/tcox/vc3vc3k_bulge" & start_snap= 47
;start_snap= 10 + fix(fload_bh_mergertime(frun))
end_snap= 94

for i= start_snap, end_snap do begin
	do_just_one, frun, snapnum=i
endfor




end







