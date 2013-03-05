function fload_cosmology, dummy

    hubble = 0.7
    omega_lambda = 0.7
    omega_matter = 0.3

    if dummy EQ 'h' then return, hubble
    if dummy EQ 'h0' then return, hubble
    if dummy EQ 'omega_lambda' then return, omega_lambda
    if dummy EQ 'omega_matter' then return, omega_matter


end




