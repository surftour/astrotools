function fload_gas_cold_j, dummy


    COMMON GalaxyHeader
    COMMON GasData


    ;sfr_rhocutoff= 0.00170994      ; gadget units
    sfr_rhocutoff= 0.00854924  ; gadget units
    cold_cutoff= 150.0      ; gadget units = 1.2 e4 T(K)


    if dummy NE -101 then begin

        sf_idx= where(rho GE sfr_rhocutoff, count_sf)

        rest_idx= where(rho LT sfr_rhocutoff, count_r)
        rest_u= u(rest_idx)


        if((count_r+count_sf) NE n_elements(mgas)) then begin
                print, "  "
                print, "PROBLEM: not all gas particles are above or below rho cutoff"
                print, "  "
                return, 0
        endif

        cold_idx= where((rest_u) LT 150.0)

        j_x= fload_gas_j(11)
        j_y= fload_gas_j(12)
        j_z= fload_gas_j(13)

        j_cold= [total(j_x(cold_idx)), total(j_y(cold_idx)), total(j_z(cold_idx))]

        return, j_cold
    endif

end

