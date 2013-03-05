function fload_rho_crit, dummy


  ;---------------------------
  ; H0 and G in gadget units
  ;---------------------------
  h= fload_cosmology('h')
  H0= 0.1*h
  G= 43007.1


  if dummy EQ 1 then begin

    ;rhocrit= 3*H0*H0/(8.0*(3.14159)*G)
    rhocrit= 3*H0*H0/(8.0*!PI*G)

    print, "rho_crit = ", rhocrit, "    (for h=",h,")"

    return, rhocrit

  endif

end


