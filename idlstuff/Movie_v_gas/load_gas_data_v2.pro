pro load_gas_data_v2, Num, Path, Snapbase, N, CurrentTime, Pos, Vel, Hsml, Mass, Temp, Id


ok=fload_snapshot_bh(Path, Num, /nopot_in_snap)

CurrentTime = fload_time(1)

cen= [0,0,0]

;use_bh= 0
use_bh= 1
if use_bh eq 1 then begin
        bhid= fload_blackhole_id(1)
        ;bhid= bhid[0]
        ;bhid= bhid[1]
        ;bhid= 200001L
        ;bhid= 280002L   ; used for z3/b4e
        ;bhid= 400002L   ; used for ds/vc3vc3e_2
        bhid= long(min(bhid))   ; uses minimum bhid - then you don't need to know the BH ID
        cen= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)

        print, "Blackhole ID: ", bhid
        print, "Blackhole center: ", cen
endif

use_calc_cen= 0
;use_calc_cen= 1
if use_calc_cen eq 1 then begin
	cen= fload_center_alreadycomp(1)
endif

N= fload_npart(0)

Pos= fltarr(3,N)
Pos[0,0:N-1]= fload_gas_xyz('x', center= cen)
Pos[1,0:N-1]= fload_gas_xyz('y', center= cen)
Pos[2,0:N-1]= fload_gas_xyz('z', center= cen)

Vel= fltarr(3,N)
Vel[0,0:N-1]= fload_gas_v('x', center= cen)
Vel[1,0:N-1]= fload_gas_v('y', center= cen)
Vel[2,0:N-1]= fload_gas_v('z', center= cen)

Hsml= fload_gas_hsml(1)

Id = fload_gas_id(1)

Mass= fload_gas_mass(1)

Temp= fload_gas_temperature(1, /keV)

end
