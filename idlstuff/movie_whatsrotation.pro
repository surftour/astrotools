pro wr, frun, i



	ok=fload_snapshot_bh(frun,i,/nopot_in_snap,/skip_center)


	   com=fload_center(1)
	   r= fload_allstars_xyz('r', center=com)
	   comvel= fload_all_comvel(1, center=com, justcenter=10.0)
	   jx= fload_allstars_j(21, center=com, vcom=comvel)   ; specific j_x
	   jy= fload_allstars_j(22, center=com, vcom=comvel)   ; specific j_y
	   jz= fload_allstars_j(23, center=com, vcom=comvel)   ; specific j_z

	   ; angular momentum within some radius 
	   radius= 5.0
	   idx= where(r lt radius)
	   if idx(0) eq -1 then stop
	   jtot_x= total(jx(idx))
	   jtot_y= total(jy(idx))
	   jtot_z= total(jz(idx))
	   jtot= sqrt(jtot_x*jtot_x + jtot_y*jtot_y + jtot_z*jtot_z)
	   nx= jtot_x/jtot
	   ny= jtot_y/jtot
	   nz= jtot_z/jtot

	   n_hat= [nx, ny, nz]
	   n_hat_tot= sqrt(nx*nx + ny*ny + nz*nz)

	   theta= acos(nz)*180./!PI
	   phi= atan(ny,nx)*180./!PI

	   print, "angular momentum vector= ", nx, ny, nz
	   print, "              n_hat_tot= ", n_hat_tot
	   print, "                  theta= ", theta
	   print, "                    phi= ", phi






	; gas
           r= fload_gas_xyz('r', center=com)
           jx= fload_gas_j(21, center=com, vcom=comvel)   ; specific j_x
           jy= fload_gas_j(22, center=com, vcom=comvel)   ; specific j_y
           jz= fload_gas_j(23, center=com, vcom=comvel)   ; specific j_z

           ; angular momentum within some radius 
           radius= 5.0
           idx= where(r lt radius)
           if idx(0) eq -1 then stop
           jtot_x= total(jx(idx))
           jtot_y= total(jy(idx))
           jtot_z= total(jz(idx))
           jtot= sqrt(jtot_x*jtot_x + jtot_y*jtot_y + jtot_z*jtot_z)
           nx= jtot_x/jtot
           ny= jtot_y/jtot
           nz= jtot_z/jtot

           n_hat= [nx, ny, nz]
           n_hat_tot= sqrt(nx*nx + ny*ny + nz*nz)

           theta= acos(nz)*180./!PI
           phi= atan(ny,nx)*180./!PI

           print, "angular momentum vector= ", nx, ny, nz
           print, "              n_hat_tot= ", n_hat_tot
           print, "                  theta= ", theta
           print, "                    phi= ", phi




end



