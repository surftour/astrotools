pro process_dir, frun, $
		g1_sid=g1_sid, g1_npart=g1_npart, g1_gsid=g1_gsid, g1_gnpart=g1_gnpart, $
		g2_sid=g2_sid, g2_npart=g2_npart, g2_gsid=g2_gsid, g2_gnpart=g2_gnpart

if not keyword_set(frun) then begin
   print, "  "
   print, "process_dir, frun, "
   print, "             g1_sid=g1_sid, g1_npart=g1_npart, g1_gsid=g1_gsid, g1_gnpart=g1_gnpart,"
   print, "             g2_sid=g2_sid, g2_npart=g2_npart, g2_gsid=g2_gsid, g2_gnpart=g2_gnpart"
   print, "  "
   print, "  "
   return
endif




	if galaxysimtype eq 'merger' then begin
		radius= 50.0
		if not keyword_set(xlen) then xlen = 100.0

		if not keyword_set(g1_sid) then begin
		   ;g1_sid=      1
		   ;g1_npart=    81920
		   ;g1_gsid=     1
		   ;g1_gnpart=   16384

		   ;g2_sid=      81921
		   ;g2_npart=    81920
		   ;g2_gsid=     81921
		   ;g2_gnpart=   16384

		   g1_sid=      1
		   g1_npart=    70000
		   g1_gsid=     1
		   g1_gnpart=   20000

		   g2_sid=      70001
		   g2_npart=    70000
		   g2_gsid=     70001
		   g2_gnpart=   20000

                   ;g1_sid=      1
                   ;g1_npart=    180000
                   ;g1_gsid=     1
                   ;g1_gnpart=   40000

                   ;g2_sid=      180001
                   ;g2_npart=    180000
                   ;g2_gsid=     180001
                   ;g2_gnpart=   40000

                   ;g1_sid=      1
                   ;g1_npart=    700000L
                   ;g1_gsid=     1
                   ;g1_gnpart=   100000L

                   ;g2_sid=      700001L
                   ;g2_npart=    199000L
                   ;g2_gsid=     0
                   ;g2_gnpart=   0
		endif
	endif

	if galaxysimtype eq 'isolated' then begin
		radius= 20.0
		if not keyword_set(xlen) then xlen= 40.0

		if not keyword_set(g1_sid) then begin
		; mh gal 1
;               g1_sid=      1
;               g1_npart=    81920
;               g1_gsid=     1
;               g1_gnpart=   16384
		; mh gal 2
;               g1_sid=      81921
;               g1_npart=    81920
;               g1_gsid=     81921
;               g1_gnpart=   16384
		; std gal 1
;                g1_sid=      1
;                g1_npart=    70000
;                g1_gsid=     1
;                g1_gnpart=   20000
		; std gal 2
;               g2_sid=      70001
;               g2_npart=    70000
;               g2_gsid=     70001
;               g2_gnpart=   20000
                ; std gal w/ hg
                ;g1_sid=      1L
                ;g1_npart=    190000L
                ;g1_gsid=     1L
                ;g1_gnpart=   60000L
		; fly-by isolated
		g1_sid=	   1L
		g1_npart=  700000L
		g1_gsid=   1L
		g1_gnpart= 100000L

		g2_sid=    1L
		g2_npart=  1L
		g2_gsid=   1L
		g2_gnpart= 1L
		endif
	endif


g1_sid= 	long(g1_sid)
g1_npart=	long(g1_npart)
g1_gsid=	long(g1_gsid)
g1_gnpart=	long(g1_gnpart)

g2_sid=		long(g2_sid)
g2_npart=	long(g2_npart)
g2_gsid=	long(g2_gsid)
g2_gnpart=	long(g2_gnpart)

print, 'galaxy1: ',g1_sid,g1_npart,g1_gsid,g1_gnpart
print, 'galaxy2: ',g2_sid,g2_npart,g2_gsid,g2_gnpart





; ------------------------------------------------------------------
; xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
; ------------------------------------------------------------------
;
;   Loop through snaps, parameterize and save, print
;
; ------------------------------------------------------------------
; xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
; ------------------------------------------------------------------


;   general
; -----------------
time=           fltarr(endsnap-startsnap+1)



;   angular momentum
; --------------------
j=		fltarr(endsnap-startsnap+1)
jx=             fltarr(endsnap-startsnap+1)
jy=             fltarr(endsnap-startsnap+1)
jz=		fltarr(endsnap-startsnap+1)
j_orbit=	fltarr(endsnap-startsnap+1)
jz_orbit=       fltarr(endsnap-startsnap+1)
j_orbit_com=    fltarr(endsnap-startsnap+1)
jz_orbit_com=   fltarr(endsnap-startsnap+1)
j_spin=		fltarr(endsnap-startsnap+1)
jz_spin=        fltarr(endsnap-startsnap+1)
j_spin_com=     fltarr(endsnap-startsnap+1)
jz_spin_com=    fltarr(endsnap-startsnap+1)
j_tot=		fltarr(endsnap-startsnap+1)
jz_tot=		fltarr(endsnap-startsnap+1)
j_tot_com=      fltarr(endsnap-startsnap+1)
jz_tot_com=     fltarr(endsnap-startsnap+1)
lambda=		fltarr(endsnap-startsnap+1)



;  galaxy 1 quantities
; ---------------------
center_1=	fltarr(endsnap-startsnap+1,3)
com_1=		fltarr(endsnap-startsnap+1,3)
comvel_1=	fltarr(endsnap-startsnap+1,3)
centervel_1=	fltarr(endsnap-startsnap+1,3)

mass_1=		fltarr(endsnap-startsnap+1)

jtot1=		fltarr(endsnap-startsnap+1)
jztot1=		fltarr(endsnap-startsnap+1)
j_tot_1=	fltarr(endsnap-startsnap+1)
jz_tot_1=	fltarr(endsnap-startsnap+1)
j_tot_com_1=    fltarr(endsnap-startsnap+1)
jz_tot_com_1=   fltarr(endsnap-startsnap+1)
j_1=		fltarr(endsnap-startsnap+1)
j_gas_1=	fltarr(endsnap-startsnap+1)
j_halo_1=      fltarr(endsnap-startsnap+1)
j_disk_1=      fltarr(endsnap-startsnap+1)
j_bulge_1=     fltarr(endsnap-startsnap+1)
j_stars_1=	fltarr(endsnap-startsnap+1)
j_orbit_1=	fltarr(endsnap-startsnap+1)
jz_orbit_1=	fltarr(endsnap-startsnap+1)
j_orbit_com_1=  fltarr(endsnap-startsnap+1)
jz_orbit_com_1= fltarr(endsnap-startsnap+1)
js_1=		fltarr(endsnap-startsnap+1)
js_gas_1=	fltarr(endsnap-startsnap+1)
js_halo_1=	fltarr(endsnap-startsnap+1)
js_disk_1=	fltarr(endsnap-startsnap+1)
js_bulge_1=	fltarr(endsnap-startsnap+1)
js_stars_1=	fltarr(endsnap-startsnap+1)
jx_1=           fltarr(endsnap-startsnap+1)
jy_1=           fltarr(endsnap-startsnap+1)
jz_1=		fltarr(endsnap-startsnap+1)
j_com_1=	fltarr(endsnap-startsnap+1)
jx_com_1=       fltarr(endsnap-startsnap+1)
jy_com_1=       fltarr(endsnap-startsnap+1)
jz_com_1=       fltarr(endsnap-startsnap+1)
jz_gas_1=	fltarr(endsnap-startsnap+1)
jz_halo_1=	fltarr(endsnap-startsnap+1)
jz_disk_1=      fltarr(endsnap-startsnap+1)
jz_bulge_1=     fltarr(endsnap-startsnap+1)
jz_stars_1=	fltarr(endsnap-startsnap+1)
jsz_1=          fltarr(endsnap-startsnap+1)
jsz_gas_1=      fltarr(endsnap-startsnap+1)
jsz_halo_1=     fltarr(endsnap-startsnap+1)
jsz_disk_1=     fltarr(endsnap-startsnap+1)
jsz_bulge_1=    fltarr(endsnap-startsnap+1)
jsz_stars_1=    fltarr(endsnap-startsnap+1)



;  galaxy 2 quantities
; ---------------------
center_2=	fltarr(endsnap-startsnap+1,3)
com_2=		fltarr(endsnap-startsnap+1,3)
comvel_2=	fltarr(endsnap-startsnap+1,3)
centervel_2=	fltarr(endsnap-startsnap+1,3)

mass_=		fltarr(endsnap-startsnap+1)

jtot2=		fltarr(endsnap-startsnap+1)
jztot2=		fltarr(endsnap-startsnap+1)
j_tot_2=	fltarr(endsnap-startsnap+1)
jz_tot_2=	fltarr(endsnap-startsnap+1)
j_tot_com_2=    fltarr(endsnap-startsnap+1)
jz_tot_com_2=   fltarr(endsnap-startsnap+1)
j_2=		fltarr(endsnap-startsnap+1)
j_gas_2=	fltarr(endsnap-startsnap+1)
j_halo_2=      fltarr(endsnap-startsnap+1)
j_disk_2=      fltarr(endsnap-startsnap+1)
j_bulge_2=     fltarr(endsnap-startsnap+1)
j_stars_2=	fltarr(endsnap-startsnap+1)
j_gasstars_2=   fltarr(endsnap-startsnap+1)
j_orbit_2=	fltarr(endsnap-startsnap+1)
jz_orbit_2=	fltarr(endsnap-startsnap+1)
j_orbit_com_2=  fltarr(endsnap-startsnap+1)
jz_orbit_com_2= fltarr(endsnap-startsnap+1)
js_2=           fltarr(endsnap-startsnap+1)
js_gas_2=       fltarr(endsnap-startsnap+1)
js_halo_2=      fltarr(endsnap-startsnap+1)
js_disk_2=      fltarr(endsnap-startsnap+1)
js_bulge_2=     fltarr(endsnap-startsnap+1)
js_stars_2=     fltarr(endsnap-startsnap+1)
js_gasstars_2=  fltarr(endsnap-startsnap+1)
jx_2=		fltarr(endsnap-startsnap+1)
jy_2=		fltarr(endsnap-startsnap+1)
jz_2=           fltarr(endsnap-startsnap+1)
j_com_2=	fltarr(endsnap-startsnap+1)
jx_com_2=       fltarr(endsnap-startsnap+1)
jy_com_2=       fltarr(endsnap-startsnap+1)
jz_com_2=       fltarr(endsnap-startsnap+1)
jz_gas_2=       fltarr(endsnap-startsnap+1)
jz_halo_2=      fltarr(endsnap-startsnap+1)
jz_disk_2=      fltarr(endsnap-startsnap+1)
jz_bulge_2=     fltarr(endsnap-startsnap+1)
jz_stars_2=     fltarr(endsnap-startsnap+1)
jz_gasstars_2=  fltarr(endsnap-startsnap+1)
jsz_2=          fltarr(endsnap-startsnap+1)
jsz_gas_2=      fltarr(endsnap-startsnap+1)
jsz_halo_2=     fltarr(endsnap-startsnap+1)
jsz_disk_2=     fltarr(endsnap-startsnap+1)
jsz_bulge_2=    fltarr(endsnap-startsnap+1)
jsz_stars_2=    fltarr(endsnap-startsnap+1)
jsz_gasstars_2= fltarr(endsnap-startsnap+1)



; --------------------------------
;  Now cycle through directory
; --------------------------------

i=0 
for snapnum= startsnap,endsnap do begin
	print, "===================================================="
	snaptoopen= snapnum*snapeveryX
	if snaptoopen gt lastsnapnum then snaptoopen= lastsnapnum
	if keyword_set(snaps_invh) then begin
		ok= fload_snapshot(frun, snaptoopen, /h) 
	endif else begin
		ok= fload_snapshot(frun, snaptoopen)
	endelse
	print, "===================================================="

	if ok eq 0 then begin


		ngas= fload_npart(0)
		nhalo= fload_npart(1)
		ndisk= fload_npart(2)
		nbulge= fload_npart(3)
		nstars= fload_npart(4)



		; centers
		; ----------
		center_1[i,*]= fload_1gal_center(g1_sid,g1_npart)
		center_b_1[i,*]= fload_1gal_baryonic_center(g1_sid,g1_npart)
		center_d_1[i,*]= fload_1gal_dm_center(g1_sid,g1_npart)
		com_1[i,*]= fload_1gal_com(g1_sid,g1_npart)
		centervel_1[i,*]= fload_1gal_dm_comvel(g1_sid,g1_npart,center=center_1[i,*],rfact=0.05)
		comvel_1[i,*]= fload_1gal_comvel(g1_sid,g1_npart,center=center_1[i,*],rfact=1.0)

		if galaxysimtype eq 'merger' then begin
			center_2[i,*]= fload_1gal_center(g2_sid,g2_npart)
			center_b_2[i,*]= fload_1gal_baryonic_center(g2_sid,g2_npart)
			center_d_2[i,*]= fload_1gal_dm_center(g2_sid,g2_npart)
			com_2[i,*]= fload_1gal_com(g2_sid,g2_npart)
			centervel_2[i,*]= fload_1gal_dm_comvel(g2_sid,g2_npart,center=center_2[i,*],rfact=0.05)
			comvel_2[i,*]= fload_1gal_comvel(g2_sid,g2_npart,center=center_2[i,*],rfact=1.0)
		endif





		; angular momentum
		; ---------------------

		; total
		j[i]= fload_all_j(0, center=[0,0,0], vcom=[0,0,0])
		jz[i]= fload_all_j(3, center=[0,0,0], vcom=[0,0,0])


		; total, galaxy 1
		jtot1[i] = fload_1gal_all_j(0, g1_sid, g1_npart, center=[0,0,0], vcom=[0,0,0])

		; spin
		;   (by component)
		j_1[i] = fload_1gal_all_j(0, g1_sid, g1_npart, center=center_1[i,*], vcom=centervel_1[i,*])
		j_gas_1[i] = fload_1gal_gas_j(0, g1_gsid, g1_gnpart, center= center_1[i,*], vcom=centervel_1[i,*])
		j_halo_1[i] = fload_1gal_halo_j(0, g1_sid, g1_npart, center= center_1[i,*], vcom=centervel_1[i,*])
		j_disk_1[i] = fload_1gal_disk_j(0, g1_sid, g1_npart, center= center_1[i,*], vcom=centervel_1[i,*])
		j_bulge_1[i] = fload_1gal_bulge_j(0, g1_sid, g1_npart, center= center_1[i,*], vcom=centervel_1[i,*])
		j_stars_1[i]= fload_1gal_newstars_j(0, g1_gsid, g1_gnpart, center= center_1[i,*], vcom=centervel_1[i,*])
		j_gasstars_1= j_gas_1 + j_stars_1

                ; spin (specific, by component)
                ;js_1[i] = fload_1gal_all_j(44, g1_sid, g1_npart, center=center_1[i,*], vcom=centervel_1[i,*])
                ;js_gas_1[i] = fload_1gal_gas_j(44, g1_gsid, g1_gnpart, center= center_1[i,*], vcom=centervel_1[i,*])
		;js_halo_1[i] = fload_1gal_halo_j(44, g1_sid, g1_npart, center= center_1[i,*], vcom=centervel_1[i,*])
		;js_disk_1[i] = fload_1gal_disk_j(44, g1_sid, g1_npart, center= center_1[i,*], vcom=centervel_1[i,*])
		;js_bulge_1[i] = fload_1gal_bulge_j(44, g1_sid, g1_npart, center= center_1[i,*], vcom=centervel_1[i,*])
                ;js_stars_1[i]= fload_1gal_newstars_j(44, g1_gsid, g1_gnpart, center= center_1[i,*], vcom=centervel_1[i,*])
                ;js_gasstars_1= js_gas_1 + js_stars_1

		; spin, azimuthal
		;   (by component)
		jztot1[i] = fload_1gal_all_j(3, g1_sid, g1_npart, center=[0,0,0], vcom=[0,0,0])
		jz_1[i] = fload_1gal_all_j(3, g1_sid, g1_npart, center=center_1[i,*], vcom=centervel_1[i,*])
                jz_gas_1[i] = fload_1gal_gas_j(3, g1_gsid, g1_gnpart, center= center_1[i,*], vcom=centervel_1[i,*])
		jz_halo_1[i] = fload_1gal_halo_j(3, g1_sid, g1_npart, center= center_1[i,*], vcom=centervel_1[i,*])
                jz_disk_1[i] = fload_1gal_disk_j(3, g1_sid, g1_npart, center= center_1[i,*], vcom=centervel_1[i,*])
                jz_bulge_1[i] = fload_1gal_bulge_j(3, g1_sid, g1_npart, center= center_1[i,*], vcom=centervel_1[i,*])
                jz_stars_1[i]= fload_1gal_newstars_j(3, g1_gsid, g1_gnpart, center= center_1[i,*], vcom=centervel_1[i,*])
                jz_gasstars_1= jz_gas_1 + jz_stars_1

		jx_1[i] = fload_1gal_all_j(1, g1_sid, g1_npart, center=center_1[i,*], vcom=centervel_1[i,*])
		jy_1[i] = fload_1gal_all_j(2, g1_sid, g1_npart, center=center_1[i,*], vcom=centervel_1[i,*])


		; spin, in true c.o.m.
		j_com_1[i] = fload_1gal_all_j(0, g1_sid, g1_npart, center=com_1[i,*], vcom=comvel_1[i,*])
		jx_com_1[i] = fload_1gal_all_j(1, g1_sid, g1_npart, center=com_1[i,*], vcom=comvel_1[i,*])
                jy_com_1[i] = fload_1gal_all_j(2, g1_sid, g1_npart, center=com_1[i,*], vcom=comvel_1[i,*])
		jz_com_1[i] = fload_1gal_all_j(3, g1_sid, g1_npart, center=com_1[i,*], vcom=comvel_1[i,*])


                ; spin, azimuthal (specific, by component)
                ;jsz_1[i] = fload_1gal_all_j(45, g1_sid, g1_npart, center=center_1[i,*], vcom=centervel_1[i,*])
                ;jsz_gas_1[i] = fload_1gal_gas_j(45, g1_gsid, g1_gnpart, center= center_1[i,*], vcom=centervel_1[i,*])
                ;jsz_halo_1[i] = fload_1gal_halo_j(45, g1_sid, g1_npart, center= center_1[i,*], vcom=centervel_1[i,*])
                ;jsz_disk_1[i] = fload_1gal_disk_j(45, g1_sid, g1_npart, center= center_1[i,*], vcom=centervel_1[i,*])
                ;jsz_bulge_1[i] = fload_1gal_bulge_j(45, g1_sid, g1_npart, center= center_1[i,*], vcom=centervel_1[i,*])
                ;jsz_stars_1[i]= fload_1gal_newstars_j(45, g1_gsid, g1_gnpart, center= center_1[i,*], vcom=centervel_1[i,*])
                ;jsz_gasstars_1= jsz_gas_1 + jsz_stars_1

		; can't do lambda until we calculate potential of only 1 galaxy
		;lambda_1[i] = j_1[i] * sqrt(-total_e_1[i]) / ( 43007.1 * mass_1[i] * mass_1[i] * sqrt(mass_1[i]) )

		if galaxysimtype eq 'merger' then begin
			jtot2[i] = fload_1gal_all_j(0, g2_sid, g2_npart, center=[0,0,0], vcom=[0,0,0])

			j_2[i] = fload_1gal_all_j(0, g2_sid, g2_npart, center=center_2[i,*], vcom=centervel_2[i,*])
			j_gas_2[i] = fload_1gal_gas_j(0, g2_gsid, g2_gnpart, center= center_2[i,*], vcom=centervel_2[i,*])
			j_halo_2[i] = fload_1gal_halo_j(0, g2_sid, g2_npart, center= center_2[i,*], vcom=centervel_2[i,*])
			j_disk_2[i] = fload_1gal_disk_j(0, g2_sid, g2_npart, center= center_2[i,*], vcom=centervel_2[i,*])
			j_bulge_2[i] = fload_1gal_bulge_j(0, g2_sid, g2_npart, center= center_2[i,*], vcom=centervel_2[i,*])
			j_stars_2[i]= fload_1gal_newstars_j(0, g2_gsid, g2_gnpart, center= center_2[i,*], vcom=centervel_2[i,*])
			j_gasstars_2= j_gas_2 + j_stars_2

	                ;js_2[i] = fload_1gal_all_j(44, g2_sid, g2_npart, center=center_2[i,*], vcom=centervel_2[i,*])
	                ;js_gas_2[i] = fload_1gal_gas_j(44, g2_gsid, g2_gnpart, center= center_2[i,*], vcom=centervel_2[i,*])
	                ;js_halo_2[i] = fload_1gal_halo_j(44, g2_sid, g2_npart, center= center_2[i,*], vcom=centervel_2[i,*])
	                ;js_disk_2[i] = fload_1gal_disk_j(44, g2_sid, g2_npart, center= center_2[i,*], vcom=centervel_2[i,*])
	                ;js_bulge_2[i] = fload_1gal_bulge_j(44, g2_sid, g2_npart, center= center_2[i,*], vcom=centervel_2[i,*])
	                ;js_stars_2[i]= fload_1gal_newstars_j(44, g2_gsid, g2_gnpart, center= center_2[i,*], vcom=centervel_2[i,*])
	                ;js_gasstars_2= js_gas_2 + js_stars_2

			jztot2[i]= fload_1gal_all_j(3, g2_sid, g2_npart, center=[0,0,0], vcom=[0,0,0])
	                jz_2[i] = fload_1gal_all_j(3, g2_sid, g2_npart, center=center_2[i,*], vcom=centervel_2[i,*])
	                jz_gas_2[i] = fload_1gal_gas_j(3, g2_gsid, g2_gnpart, center= center_2[i,*], vcom=centervel_2[i,*])
			jz_halo_2[i] = fload_1gal_halo_j(3, g2_sid, g2_npart, center= center_2[i,*], vcom=centervel_2[i,*])
                        jz_disk_2[i] = fload_1gal_disk_j(3, g2_sid, g2_npart, center= center_2[i,*], vcom=centervel_2[i,*])
                        jz_bulge_2[i] = fload_1gal_bulge_j(3, g2_sid, g2_npart, center= center_2[i,*], vcom=centervel_2[i,*])
	                jz_stars_2[i]= fload_1gal_newstars_j(3, g2_gsid, g2_gnpart, center= center_2[i,*], vcom=centervel_2[i,*])
	                jz_gasstars_2= jz_gas_2 + jz_stars_2

			jx_2[i] = fload_1gal_all_j(1, g2_sid, g2_npart, center=center_2[i,*], vcom=centervel_2[i,*])
			jy_2[i] = fload_1gal_all_j(2, g2_sid, g2_npart, center=center_2[i,*], vcom=centervel_2[i,*])

			j_com_2[i] = fload_1gal_all_j(0, g2_sid, g2_npart, center=com_2[i,*], vcom=comvel_2[i,*])
			jx_com_2[i] = fload_1gal_all_j(1, g2_sid, g2_npart, center=com_2[i,*], vcom=comvel_2[i,*])
                        jy_com_2[i] = fload_1gal_all_j(2, g2_sid, g2_npart, center=com_2[i,*], vcom=comvel_2[i,*])
			jz_com_2[i] = fload_1gal_all_j(3, g2_sid, g2_npart, center=com_2[i,*], vcom=comvel_2[i,*])

	                ;jsz_2[i] = fload_1gal_all_j(45, g2_sid, g2_npart, center=center_2[i,*], vcom=centervel_2[i,*])
	                ;jsz_gas_2[i] = fload_1gal_gas_j(45, g2_gsid, g2_gnpart, center= center_2[i,*], vcom=centervel_2[i,*])
	                ;jsz_halo_2[i] = fload_1gal_halo_j(45, g2_sid, g2_npart, center= center_2[i,*], vcom=centervel_2[i,*])
	                ;jsz_disk_2[i] = fload_1gal_disk_j(45, g2_sid, g2_npart, center= center_2[i,*], vcom=centervel_2[i,*])
	                ;jsz_bulge_2[i] = fload_1gal_bulge_j(45, g2_sid, g2_npart, center= center_2[i,*], vcom=centervel_2[i,*])
	                ;jsz_stars_2[i]= fload_1gal_newstars_j(45, g2_gsid, g2_gnpart, center= center_2[i,*], vcom=centervel_2[i,*])
	                ;jsz_gasstars_2= jsz_gas_2 + jsz_stars_2

			;lambda_2[i] = j_2[i] * sqrt(-total_e_2[i]) / ( 43007.1 * mass_2[i] * mass_2[i] * sqrt(mass_2[i]) )



			; calculate total (gal 1 + 2) spins
			; -----------------------------------
			jx_spin= jx_1[i] + jx_2[i]
			jy_spin= jy_1[i] + jy_2[i]
			jz_spin[i]= jz_1[i] + jz_2[i]

			j_spin[i]= sqrt(jx_spin*jx_spin + jy_spin*jy_spin + jz_spin[i]*jz_spin[i])


			; com spins
                        jx_spin_com= jx_com_1[i] + jx_com_2[i]
                        jy_spin_com= jy_com_1[i] + jy_com_2[i]
                        jz_spin_com[i]= jz_com_1[i] + jz_com_2[i]

                        j_spin_com[i]= sqrt(jx_spin_com*jx_spin_com + jy_spin_com*jy_spin_com + jz_spin_com[i]*jz_spin_com[i])



			; orbital angular momentum
			; --------------------------

			;center
			m1= fload_1gal_mass(g1_sid,g1_npart)
			jx1= m1*(center_1[i,1]*centervel_1[i,2] - center_1[i,2]*centervel_1[i,1])
			jy1=-m1*(center_1[i,0]*centervel_1[i,2] - center_1[i,2]*centervel_1[i,0])
			jz_orbit_1[i]= m1*(center_1[i,0]*centervel_1[i,1] - center_1[i,1]*centervel_1[i,0])

			m2= fload_1gal_mass(g2_sid,g2_npart)
                        jx2= m2*(center_2[i,1]*centervel_2[i,2] - center_2[i,2]*centervel_2[i,1])
                        jy2=-m2*(center_2[i,0]*centervel_2[i,2] - center_2[i,2]*centervel_2[i,0])
                        jz_orbit_2[i]= m2*(center_2[i,0]*centervel_2[i,1] - center_2[i,1]*centervel_2[i,0])

			jx_orbit= jx1 + jx2
			jy_orbit= jy1 + jy2
			jz_orbit[i]= jz_orbit_1[i] + jz_orbit_2[i]
			j_orbit[i]= sqrt(jx_orbit*jx_orbit + jy_orbit*jy_orbit + jz_orbit[i]*jz_orbit[i])
			j_orbit_1[i]= sqrt(jx1*jx1 + jy1*jy1 + jz_orbit_1[i]*jz_orbit_1[i])
                        j_orbit_2[i]= sqrt(jx2*jx2 + jy2*jy2 + jz_orbit_2[i]*jz_orbit_2[i])


			;c.o.m.
                        jx1_com= m1*(com_1[i,1]*comvel_1[i,2] - com_1[i,2]*comvel_1[i,1])
                        jy1_com=-m1*(com_1[i,0]*comvel_1[i,2] - com_1[i,2]*comvel_1[i,0])
                        jz_orbit_com_1[i]= m1*(com_1[i,0]*comvel_1[i,1] - com_1[i,1]*comvel_1[i,0])

                        jx2_com= m2*(com_2[i,1]*comvel_2[i,2] - com_2[i,2]*comvel_2[i,1])
                        jy2_com=-m2*(com_2[i,0]*comvel_2[i,2] - com_2[i,2]*comvel_2[i,0])
                        jz_orbit_com_2[i]= m2*(com_2[i,0]*comvel_2[i,1] - com_2[i,1]*comvel_2[i,0])

                        jx_orbit_com= jx1_com + jx2_com
                        jy_orbit_com= jy1_com + jy2_com
                        jz_orbit_com[i]= jz_orbit_com_1[i] + jz_orbit_com_2[i]
                        j_orbit_com[i]= sqrt(jx_orbit_com*jx_orbit_com + jy_orbit_com*jy_orbit_com + jz_orbit_com[i]*jz_orbit_com[i])
                        j_orbit_com_1[i]= sqrt(jx1_com*jx1_com + jy1_com*jy1_com + jz_orbit_com_1[i]*jz_orbit_com_1[i])
                        j_orbit_com_2[i]= sqrt(jx2_com*jx2_com + jy2_com*jy2_com + jz_orbit_com_2[i]*jz_orbit_com_2[i])


			; calculate totals
			; -------------------
			jx_tot= jx_spin+jx_orbit
			jy_tot= jy_spin+jy_orbit
			jz_tot[i]= jz_spin[i]+jz_orbit[i]
			j_tot[i]= sqrt(jx_tot*jx_tot + jy_tot*jy_tot + jz_tot[i]*jz_tot[i])


                	; total, by hand - gal 1
                	jx_tot_1= jx_1[i] + jx1
                	jy_tot_1= jy_1[i] + jy1
                	jz_tot_1[i]= jz_1[i] + jz_orbit_1[i]
			j_tot_1[i]= sqrt(jx_tot_1*jx_tot_1 + jy_tot_1*jy_tot_1 + jz_tot_1[i]*jz_tot_1[i])

                        ; total, by hand - gal 2
                        jx_tot_2= jx_2[i] + jx2
                        jy_tot_2= jy_2[i] + jy2
                        jz_tot_2[i]= jz_2[i] + jz_orbit_2[i]
			j_tot_2[i]= sqrt(jx_tot_2*jx_tot_2 + jy_tot_2*jy_tot_2 + jz_tot_2[i]*jz_tot_2[i])


                        ; c.o.m. total
                        jx_tot_com= jx_spin_com+jx_orbit_com
                        jy_tot_com= jy_spin_com+jy_orbit_com
                        jz_tot_com[i]= jz_spin_com[i]+jz_orbit_com[i]
                        j_tot_com[i]= sqrt(jx_tot_com*jx_tot_com + jy_tot_com*jy_tot_com + jz_tot_com[i]*jz_tot_com[i])
			; com total, by hand - gal 1
                        jx_tot_com_1= jx_com_1[i] + jx1_com
                        jy_tot_com_1= jy_com_1[i] + jy1_com
                        jz_tot_com_1[i]= jz_com_1[i] + jz_orbit_com_1[i]
                        j_tot_com_1[i]= sqrt(jx_tot_com_1*jx_tot_com_1 + jy_tot_com_1*jy_tot_com_1 + jz_tot_com_1[i]*jz_tot_com_1[i])
                        ; com total, by hand - gal 2
                        jx_tot_com_2= jx_com_2[i] + jx2_com
                        jy_tot_com_2= jy_com_2[i] + jy2_com
                        jz_tot_com_2[i]= jz_com_2[i] + jz_orbit_com_2[i]
                        j_tot_com_2[i]= sqrt(jx_tot_com_2*jx_tot_com_2 + jy_tot_com_2*jy_tot_com_2 + jz_tot_com_2[i]*jz_tot_com_2[i])




		endif



                ; orbital energies
                ; ------------------
		if galaxysimtype eq 'merger' then begin
                   rsep2= (center_1[i,*]-center_2[i,*])*(center_1[i,*]-center_2[i,*])
                   rsep= sqrt(total(rsep2))
		   if i eq 0 then rstart= rsep
                   orbit_ke[i]= 0.5*m1*total(comvel_1[i,*]*comvel_1[i,*])  + 0.5*m2*total(comvel_2[i,*]*comvel_2[i,*])
                   if rsep lt rstart*0.05 then rsep= rstart*0.05       ; things blow up if we plot this all the way
		   orbit_pe[i]= -43007.1 * m1 * m2 / rsep
		endif





		; time
		; -------
                time[i]= fload_time(1)




        endif
        i=i+1
endfor



; -------------------
;  Angular Momentum
; -------------------

if galaxysimtype eq 'merger' then begin

	; using our 'center'
	; --------------------
	plot_angular_momentum, time, j, j_spin, j3=j_orbit, j4=j_tot, $
			lbl1="Total", lbl2="Spin", lbl3="Orbital", $
			ytit='|J|, |S|, |L|', $
                        frun, sendto, /percentage, filename=frun+"/jtot.eps"

	; azimuthal component of total's
	plot_angular_momentum, time, jz, jz_spin, j3=jz_orbit, j4=jz_tot, $ 
                        lbl1="Total", lbl2="Spin", lbl3="Orbital", $
                        ytit='J!Dz!N, S!Dz!N, L!Dz!N', $
                        frun, sendto, /percentage, filename=frun+"/jztot.eps"



	; now COM
	; --------
        plot_angular_momentum, time, j, j_spin_com, j3=j_orbit_com, j4=j_tot_com, $
                        lbl1="Total", lbl2="Spin", lbl3="Orbital", $
                        ytit='|J|, |S|, |L|', $
			msg='COM', $
                        frun, sendto, /percentage, filename=frun+"/jcomtot.eps"

        ; azimuthal component of total's
        plot_angular_momentum, time, jz, jz_spin_com, j3=jz_orbit_com, j4=jz_tot_com, $
                        lbl1="Total", lbl2="Spin", lbl3="Orbital", $
                        ytit='J!Dz!N, S!Dz!N, L!Dz!N', $
			msg='COM', $
                        frun, sendto, /percentage, filename=frun+"/jzcomtot.eps"

endif


;same as above for each galaxy separately
plot_angular_momentum, time, jtot1, j_1, j3=j_orbit_1, j4=j_tot_1, $
	lbl1="Total", lbl2="Spin", lbl3="Orbital", $
	ytit='|J|, |S|, |L|', $
	msg='Galaxy 1', $
	frun, sendto, /percentage, filename=frun+"/j1.eps"

plot_angular_momentum, time, jztot1, jz_1, j3=jz_orbit_1, j4=jz_tot_1, $
        lbl1="Total", lbl2="Spin", lbl3="Orbital", $
        ytit='J!Dz!N, S!Dz!N, L!Dz!N', $
        msg='Galaxy 1', $
        frun, sendto, /percentage, filename=frun+"/jz1.eps"



if galaxysimtype eq 'merger' then begin
        plot_angular_momentum, time, jtot2, j_2, j3=j_orbit_2, j4=j_tot_2, $
                        lbl1="Total", lbl2="Spin", lbl3="Orbital", $
                        ytit='|J|, |S|, |L|', $
			msg='Galaxy 2', $
                        frun, sendto, /percentage, filename=frun+"/j2.eps"

        plot_angular_momentum, time, jztot2, jz_2, j3=jz_orbit_2, j4=jz_tot_2, $
                        lbl1="Total", lbl2="Spin", lbl3="Orbital", $
			ytit='J!Dz!N, S!Dz!N, L!Dz!N', $
                        msg='Galaxy 2', $
                        frun, sendto, /percentage, filename=frun+"/jz2.eps"
endif




; spin of each galaxy (by component)
; -----------------------------------
plot_angularmomentum_comp, time, j_1, j_gas_1, j_halo_1, j_disk_1, j_bulge_1, j_stars_1, $
                       frun, sendto, /percentage, msg='Galaxy 1', filename=frun+"/j1comp.eps"
plot_angularmomentum_comp_delta, time, j_1, j_gas_1, j_halo_1, j_disk_1, j_bulge_1, j_stars_1, $
                        frun, sendto, /percentage, msg='Galaxy 1', filename=frun+"/j1comp_delta.eps"
if galaxysimtype eq 'merger' then begin
       plot_angularmomentum_comp, time, j_2, j_gas_2, j_halo_2, j_disk_2, j_bulge_2, j_stars_2, $
                       frun, sendto, /percentage, msg='Galaxy 2', filename=frun+"/j2comp.eps"
       plot_angularmomentum_comp_delta, time, j_2, j_gas_2, j_halo_2, j_disk_2, j_bulge_2, j_stars_2, $
                        frun, sendto, /percentage, msg='Galaxy 2', filename=frun+"/j2comp_delta.eps"


	; gas angular momentum for each galaxy
	if ngas gt 0 then begin
		plot_angular_momentum, time, j_gas_1, j_gas_2, msg="Gas", $
			lbl1='Galaxy 1', lbl2='Galaxy 2', $
			ytit= '|S|', $
                        frun, sendto, /percentage, filename=frun+"/jgas.eps"
		plot_angular_momentum, time, jz_gas_1, jz_gas_2, msg="Gas", $
			lbl1='Galaxy 1', lbl2='Galaxy 2', $
			ytit= 'S!Dz!N', $
                        frun, sendto, /percentage, filename=frun+"/jzgas.eps"

		; gas+new stars angular momentum
		plot_angular_momentum, time, j_gasstars_1, j_gasstars_2, msg='Gas+Stars', $
			lbl1='Galaxy 1', lbl2='Galaxy 2', $
			ytit= '|S|', $
                        frun, sendto, /percentage, filename=frun+"/jgasstars.eps"
		plot_angular_momentum, time, jz_gasstars_1, jz_gasstars_2, msg='Gas+Stars', $
			lbl1='Galaxy 1', lbl2='Galaxy 2', $
			ytit= 'S!Dz!N', $
                        frun, sendto, /percentage, filename=frun+"/jzgasstars.eps"
	endif


        ; halo angular momentum for each galaxy
        if nhalo gt 0 then begin
                plot_angular_momentum, time, j_halo_1, j_halo_2, msg="Halo", $
			lbl1='Galaxy 1', lbl2='Galaxy 2', $
                        ytit= '|S|', $
                        frun, sendto, /percentage, filename=frun+"/jhalo.eps"
                plot_angular_momentum, time, jz_halo_1, jz_halo_2, msg="Halo", $
			lbl1='Galaxy 1', lbl2='Galaxy 2', $
                        ytit= 'S!Dz!N', $
                        frun, sendto, /percentage, filename=frun+"/jzhalo.eps"
        endif

        ; disk angular momentum for each galaxy
        if ndisk gt 0 then begin
                plot_angular_momentum, time, j_disk_1, j_disk_2, msg="Disk", $
			lbl1='Galaxy 1', lbl2='Galaxy 2', $
                        ytit= '|S|', $
                        frun, sendto, /percentage, filename=frun+"/jdisk.eps"
                plot_angular_momentum, time, jz_disk_1, jz_disk_2, msg="Disk", $
			lbl1='Galaxy 1', lbl2='Galaxy 2', $
                        ytit= 'S!Dz!N', $
                        frun, sendto, /percentage, filename=frun+"/jzdisk.eps"
        endif

        ; bulge angular momentum for each galaxy
        if nbulge gt 0 then begin
                plot_angular_momentum, time, j_bulge_1, j_bulge_2, msg="Bulge", $
			lbl1='Galaxy 1', lbl2='Galaxy 2', $
                        ytit= '|S|', $
                        frun, sendto, /percentage, filename=frun+"/jbulge.eps"
                plot_angular_momentum, time, jz_bulge_1, jz_bulge_2, msg="Bulge", $
			lbl1='Galaxy 1', lbl2='Galaxy 2', $
                        ytit= 'S!Dz!N', $
                        frun, sendto, /percentage, filename=frun+"/jzbulge.eps"
        endif


endif






; ---------------------------------------------------------------------------------



print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"



end

























;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Angular Momentum
;     ---------------------------------------------
;
;   plots total angular momentum (magnitude of J) as a function of time
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro plot_angular_momentum, time, j1, j2, j3=j3, j4=j4, $
			percentage=percentage, msg=msg, $
			lbl1=lbl1, lbl2=lbl2, lbl3=lbl3, lbl4=lbl4, $
			ytit=ytit, $
                        frun, sendto, filename=filename


if not keyword_set(sendto) then sendto= 'ps'
if n_elements(time) lt 2 then begin
   print, "  "
   print, "PROBLEM:  plot_angular_momentum"
   print, "  "
   return
endif


setup_plot_stuff, sendto, filename=filename
   


if not keyword_set(lbl1) then lbl1= 'Galaxy 1'
if not keyword_set(lbl2) then lbl2= 'Galaxy 2'

if not keyword_set(lbl3) then lbl3= ''
if not keyword_set(lbl4) then lbl4= ''

if not keyword_set(ytit) then ytit= '|J|'




;--------------------------------------
;  Print the Shit
;--------------------------------------


xmax = max(time)
xmin = 0


; resize j's
; -------------
j1temp=j1/1e5
j2temp=j2/1e5
if keyword_set(j3) then j3temp=j3/1e5
if keyword_set(j4) then j4temp=j4/1e5



; linear y
; ---------
js= [j1temp,j2temp]
if keyword_set(j3) then js= [js, j3temp]
if keyword_set(j4) then js= [js, j4temp]

ymax= max(js)
ymax= fix(1.2*ymax)
if ymax eq 0 then ymax= 1.2*max(js)


ymin= min(js)
ymin=fix(ymin-1.0)
if ymin gt 0 then ymin= 0
if (min(js) gt -0.2 and ymin lt 0) then ymin=min(js)



; log y
; ------
;js= [abs(j1),abs(j2)]
;if keyword_set(j3) then js= [js, abs(j3)]

;ymax= 10.0^(long(alog10(max(js))+1))
;ymin= 10.0^(long(alog10(min(js)+0.1)))

;if ymax eq 0 then ymax= ymin*10


;idx= where(j1 le 0)
;if idx(0) ne -1 then begin
;	idx= where(j1 eq 0)
;	if idx(0) ne -1 then j1(idx)= 0.1
;	idx= where(j1 lt 0)
;	if idx(0) ne -1 then j1neg= abs(j1(idx))
;	if idx(0) ne -1 then j1(idx)= abs(j1(idx))
;endif

;idx= where(j2 le 0)
;if idx(0) ne -1 then begin
;        idx= where(j2 eq 0)
;        if idx(0) ne -1 then j2(idx)= 0.1
;        idx= where(j2 lt 0)
;        if idx(0) ne -1 then j2neg= abs(j2(idx))
;        if idx(0) ne -1 then j2(idx)= abs(j2(idx))
;endif


;ymax = 10
;ymin = 0 
;js= [j1,j2]
;ymax= 10.0^(long(alog10(max(js))+1))
;ymin= 10.0^(long(alog10(min(js))))



;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        ;/ylog, $
        color= 0, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle="Time (Gyr)", $
        ytitle=ytit, $
        /nodata



if sendto EQ 'x' then begin
        ;oplot, time, j1, psym=-2, color= getcolor('red')
        ;oplot, time, j2, psym=-5, color= getcolor('green')
	;if n_elements(j1neg) gt 0 then oplot, time, j1neg, psym=-2, color= getcolor('blue')
	;if n_elements(j2neg) gt 0 then oplot, time, j2neg, psym=-5, color= getcolor('blue')
endif else begin
        oplot, time, j1temp, psym=-3, linestyle= 0, color= 50, thick=3.0
        oplot, time, j2temp, psym=-3, color= 150, linestyle=1, thick=3.0
	if keyword_set(j3) then oplot, time, j3temp, psym=-3, linestyle= 3, color= 100, thick=3.0
        ;if n_elements(j1neg) gt 0 then oplot, time, j1neg, psym=-2, color= 100
        ;if n_elements(j2neg) gt 0 then oplot, time, j2neg, psym=-5, color= 100
	if keyword_set(j4) then oplot, time, j4temp, psym=-3, linestyle= 2, color= 50, thick=3.0
endelse



if ymin lt 0 then begin
	oplot, [xmin,xmax], [0,0], psym=-3, linestyle=2, color= 0
endif

;xyouts, 0.7, 0.90, fload_fid(1), /normal, charthick=1, size=1.33, color=0
xyouts, 0.7, 0.88, lbl1, /normal, charthick=2.0, size=1.0, color=50
xyouts, 0.7, 0.83, lbl2, /normal, charthick=2.0, size=1.0, color=150
if keyword_set(j3) then xyouts, 0.7, 0.78, lbl3, /normal, charthick=2.0, size=1.0, color=100
if keyword_set(j4) then xyouts, 0.7, 0.73, lbl4, /normal, charthick=2.0, size=1.0, color=50


if keyword_set(msg) then begin
	xyouts, 0.25, 0.90, msg, /normal, color= 0, size= 1.0
endif

if keyword_set(percentage) then begin
	lastidx= n_elements(j1)-1
	if j1[0] gt 0 then totper= 100.0*j1[lastidx]/j1[0] else totper= 0
	totlbl= strcompress(string(totper),/remove_all)
	totlbl= strmid(totlbl,0,5)+'%'
	totlbl= '('+totlbl+')'
	xyouts, 0.55, 0.88, totlbl, /normal, size= 1.0, color=0

        lastidx= n_elements(j2)-1
        if j2[0] gt 0 then totper= 100.0*j2[lastidx]/j2[0] else totper= 0
        totlbl= strcompress(string(totper),/remove_all)
        totlbl= strmid(totlbl,0,5)+'%'
        totlbl= '('+totlbl+')'
        xyouts, 0.55, 0.83, totlbl, /normal, size= 1.0, color=0
endif

;--------------------------------------
;--------------------------------------

if (sendto EQ 'ps') then device, /close


end















;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Angular Momentum
;     -------------------------------------------
;  plot angular momentum broken down by galactic component
;  
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro plot_angularmomentum_comp, time, totj, gasj, haloj, diskj, bulgej, newstarsj, $
			percentage=percentage, msg=msg, $
                        frun, sendto, filename=filename


if not keyword_set(sendto) then sendto='ps'
if n_elements(time) lt 2 then begin
   print, "  "
   print, "PROBLEM: plot_angularmomentum_comp"
   print, "  "
   return
endif


setup_plot_stuff, sendto, filename=filename



; ----------------------------------
;  Set Parameters
; ----------------------------------

xmax= max(time)
xmin= min(time)

ymax= -1000
if max(totj) gt ymax then ymax= max(totj)
if max(gasj) gt ymax then ymax= max(gasj)
if max(haloj) gt ymax then ymax= max(haloj)
if max(diskj) gt ymax then ymax= max(diskj)
if max(bulgej) gt ymax then ymax= max(bulgej)
if max(gasj+newstarsj) gt ymax then ymax= max(gasj+newstarsj)

ymax= ymax*5.0


ymin= 1.0e+10
if min(totj) lt ymin and min(totj) gt 1 then ymin= min(totj)
if min(gasj) lt ymin and min(gasj) gt 1 then ymin= min(gasj)
if min(haloj) lt ymin and min(haloj) gt 1 then ymin= min(haloj)
if min(diskj) lt ymin and min(diskj) gt 1 then ymin= min(diskj)
if min(bulgej) lt ymin and min(bulgej) gt 1 then ymin= min(bulgej)
if min(gasj+newstarsj) lt ymin and min(gasj+newstarsj) gt 1 then ymin= min(gasj+newstarsj)

ymin= ymin/5.0



; take out zero's and store a negative part
idx= where(totj le 0)
if idx(0) ne -1 then begin
        idx= where(totj eq 0)
        if idx(0) ne -1 then totj(idx)= 0.1
        idx= where(totj lt 0)
        if idx(0) ne -1 then totjneg= abs(totj(idx))
        if idx(0) ne -1 then totj(idx)= abs(totj(idx))
endif





;---------------------------
;  Print it
;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]


plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        color= 0, $
	/ylog, $
        xstyle=1, ystyle=1, $ 
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle="Time (Gyr)", $
        ytitle="j (Gadget Units)", $
        /nodata



if sendto EQ 'x' then begin
        oplot, time, totj, psym=-3, color= 0
        oplot, time, gasj, psym=-3, color= getcolor('red')
        oplot, time, haloj, psym=-3, color= getcolor('cyan')
        oplot, time, diskj, psym=-3, color= getcolor('green')
        oplot, time, bulgej, psym=-3, color= getcolor('blue')
	oplot, time, newstarsj, psym=-3, color= getcolor('navy')

        xyouts, 0.7, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0
        xyouts, 0.7, 0.86, 'Total', /normal, charthick=1, size=1.33, color= 0
        xyouts, 0.7, 0.82, 'Gas', /normal, charthick=1, size=1.33, color= getcolor('blue')
        xyouts, 0.7, 0.78, 'Halo', /normal, charthick=1, size=1.33, color= getcolor('green')
        xyouts, 0.7, 0.74, 'Disk', /normal, charthick=1, size=1.33, color= getcolor('red')
        xyouts, 0.7, 0.70, 'Bulge', /normal, charthick=1, size=1.33, color= getcolor('cyan')
	xyouts, 0.7, 0.66, 'G+NS', /normal, charthick=1, size=1.33, color= getcolor('cyan')

endif else begin
        if totj[0] ne 0 then oplot, time, totj, psym=-3, color= 0, thick=3.0
        if gasj[0] ne 0 then oplot, time, gasj, psym=-3, color= 100, thick=3.0
        if haloj[0] ne 0 then oplot, time, haloj, psym=-3, color= 150, thick=3.0
        if diskj[0] ne 0 then oplot, time, diskj, psym=-3, color= 50, thick=3.0
        if bulgej[0] ne 0 then oplot, time, bulgej, psym=-3, color= 80, thick=3.0
	if (newstarsj[0]+gasj[0]) ne 0 then oplot, time, newstarsj+gasj, psym=-3, color= 180, thick=3.0

	x0= 0.80
	y0= 0.88
        if totj[0] ne 0 then xyouts, x0, y0, 'Total', /normal, charthick=3.0, size=1.33, color= 0
        if gasj[0] ne 0 then xyouts, x0, y0-0.04, 'Gas', /normal, charthick=3.0, size=1.33, color= 100
        if haloj[0] ne 0 then xyouts, x0, y0-0.08, 'Halo', /normal, charthick=3.0, size=1.33, color= 150
        if diskj[0] ne 0 then xyouts, x0, y0-0.12, 'Disk', /normal, charthick=3.0, size=1.33, color= 50
        if bulgej[0] ne 0 then xyouts, x0, y0-0.16, 'Bulge', /normal, charthick=3.0, size=1.33, color= 80
	if (newstarsj[0]+gasj[0]) ne 0 then xyouts, x0, y0-0.20, 'G+NS', /normal, charthick=3.0, size=1.33, color= 180
endelse


; print message if set
if keyword_set(msg) then xyouts, 0.22, 0.85, msg, /normal, size=1, color= 0




;--------------------------------------
;--------------------------------------

if (sendto EQ 'ps') then device, /close


end




;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     change in gas energy vs. time
;     -------------------------------------------
;    we'll break this down by component: ke, pe, thermal and total
;  
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro plot_angularmomentum_comp_delta, time, totj, gasj, haloj, diskj, bulgej, newstarsj, $
			percentage=percentage, msg=msg, $
                        frun, sendto, filename=filename



if not keyword_set(sendto) then sendto='ps'
if n_elements(time) lt 2 then begin
   print, "  "
   print, "PROBLEM: plot_angularmomentum_comp_delta"
   print, "  "
   return
endif


setup_plot_stuff, sendto, filename=filename



; ----------------------------------
;  Set Parameters
; ----------------------------------

xmax= max(time)
xmin= min(time)

ymax= 3
ymin= -3



; delta specific j
deltatot= totj-totj[0]
deltagas= gasj-gasj[0]
deltahalo= haloj-haloj[0]
deltadisk= diskj-diskj[0]
deltabulge= bulgej-bulgej[0]

gnsj=gasj+newstarsj
deltagns= gnsj-gnsj[0]



maxs= [deltatot, deltagas, deltahalo, deltadisk, deltabulge, deltagns]

ymax= max(maxs)
ymin= min(maxs)


; percentage changes 
lastidx= n_elements(totj)-1
if totj[0] ne 0 then totper=100.0*deltatot[lastidx]/totj[0]
if gasj[0] ne 0 then gasper=100.0*deltagas[lastidx]/gasj[0]
if haloj[0] ne 0 then haloper=100.0*deltahalo[lastidx]/haloj[0]
if diskj[0] ne 0 then diskper=100.0*deltadisk[lastidx]/diskj[0]
if bulgej[0] ne 0 then bulgeper=100.0*deltabulge[lastidx]/bulgej[0]
if gnsj[0] ne 0 then gnsper=100.0*deltagns[lastidx]/gnsj[0]

if totj[0] ne 0 then dtot=deltatot[lastidx]
if gasj[0] ne 0 then dgas=deltagas[lastidx]
if haloj[0] ne 0 then dhalo=deltahalo[lastidx]
if diskj[0] ne 0 then ddisk=deltadisk[lastidx]
if bulgej[0] ne 0 then dbulge=deltabulge[lastidx]
if gnsj[0] ne 0 then dgns=deltagns[lastidx]




; ----------------------------------
;   Try this plot thingy
; ----------------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle="Time (Gyr)", $
        ytitle="!4D!3j (Gadget Units)", $
        /nodata



if sendto EQ 'x' then begin

endif else begin
        if totj[0] ne 0 then oplot, time, deltatot, thick=2.0, psym=-3, color= 0
        if gasj[0] ne 0 then oplot, time, deltagas, thick=1.0, psym=-3, color= 100
        if haloj[0] ne 0 then oplot, time, deltahalo, thick=2.0, psym=-3, color= 150
        if diskj[0] ne 0 then oplot, time, deltadisk, thick=2.0, psym=-3, color= 50
        if bulgej[0] ne 0 then oplot, time, deltabulge, thick=2.0, psym=-3, color= 80
	if gnsj[0] ne 0 then oplot, time, deltagns, thick=2.0, psym=-3, color= 180

        if totj[0] ne 0 then xyouts, 0.82, 0.90, 'Total ', /normal, charthick=1, size=1.33, color= 0
        if gasj[0] ne 0 then xyouts, 0.82, 0.86, 'Gas', /normal, charthick=1, size=1.33, color= 100
        if haloj[0] ne 0 then xyouts, 0.82, 0.82, 'Halo', /normal, charthick=1, size=1.33, color= 150
        if diskj[0] ne 0 then xyouts, 0.82, 0.78, 'Disk', /normal, charthick=1, size=1.33, color= 50
        if bulgej[0] ne 0 then xyouts, 0.82, 0.74, 'Bulge', /normal, charthick=1, size=1.33, color= 80
	if gnsj[0] ne 0 then xyouts, 0.82, 0.74, 'G+NS', /normal, charthick=1, size=1.33, color= 180

endelse


; plot zero
x=findgen(10)
y=0.0*x
oplot, x, y, psym=-3, color= 0


; print message if set
if keyword_set(msg) then xyouts, 0.22, 0.85, msg, /normal, size=1, color= 0


; print percentage of energy change
if totj[0] ne 0 then begin
  totlbl1= strcompress(string(dtot),/remove_all)
  totlbl2= strcompress(string(totper),/remove_all)
  totlbl= strmid(totlbl1,0,4)+' ('+strmid(totlbl2,0,5)+'%)'
  xyouts, 0.52, 0.90, totlbl, /normal,size= 1.2, color=0
endif

if gasj[0] ne 0 then begin
  glbl1= strcompress(string(dgas),/remove_all)
  glbl2= strcompress(string(gasper),/remove_all)
  glbl= strmid(glbl1,0,4)+' ('+strmid(glbl2,0,5)+'%)'
  xyouts, 0.52, 0.86, glbl, /normal,size= 1.2, color=0
endif

if haloj[0] ne 0 then begin
  halolbl1= strcompress(string(dhalo),/remove_all)
  halolbl2= strcompress(string(haloper),/remove_all)
  halolbl= strmid(halolbl1,0,4)+' ('+strmid(halolbl2,0,5)+'%)'
  xyouts, 0.52, 0.82, halolbl, /normal,size= 1.2, color=0
endif

if diskj[0] ne 0 then begin
  disklbl1= strcompress(string(ddisk),/remove_all)
  disklbl2= strcompress(string(diskper),/remove_all)
  disklbl= strmid(disklbl1,0,4)+' ('+strmid(disklbl2,0,5)+'%)'
  xyouts, 0.52, 0.78, disklbl, /normal,size= 1.2, color=0
endif

if bulgej[0] ne 0 then begin
  bulgelbl1= strcompress(string(dbulge),/remove_all)
  bulgelbl2= strcompress(string(bulgeper),/remove_all)
  bulgelbl= strmid(bulgelbl1,0,4)+' ('+strmid(bulgelbl2,0,5)+'%)'
  xyouts, 0.52, 0.74, bulgelbl, /normal,size= 1.2, color=0
endif

if gnsj[0] ne 0 then begin
  gnslbl1= strcompress(string(dgns),/remove_all)
  gnslbl2= strcompress(string(gnsper),/remove_all)
  gnslbl= strmid(gnslbl1,0,4)+' ('+strmid(gnslbl2,0,5)+'%)'
  xyouts, 0.52, 0.74, gnslbl, /normal,size= 1.2, color=0
endif




; -------------------------------

if (sendto EQ 'ps') then device, /close




end





















;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Angular Momentum
;     ---------------------------------------------
;
;   plots histogram of angular momentum  (i.e. the j distribution)
;   we will just do whatever is passed, so could either pass j per
;   particle, specific j or j_z
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------




; actually will do this with a routine we have already processed,
; angmo_hist, so look there.














;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Angular Momentum
;     ---------------------------------------------
;
;   plots radial distribution of angular momentum
;   again we'll just plot whatever is passed, so we could plot total
;   j per particle, specific j or j_z
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro plot_angmo_radial_distribution, rad, j, $
			msg=msg, $
			frun, sendto, filename=filename



if not keyword_set(sendto) then sendto='ps'
if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_angmo_radial_distribution"
   print, "  "
   return
endif


setup_plot_stuff, sendto, filename=filename



; -----------------------------
; set up constants
; -----------------------------

bins = 100

xmax = 3
xmin = -2
ymax = 1e+12
ymin = 1


; ---------------------------
;  Load it Up
; ---------------------------

r_log= fltarr(bins)
jdist= fltarr(bins)


;r= fload_gas_xyz('r')    ; or fload_all_xyz('r')    ; or fload_halo_xyz('r')
;j= fload_gas_j(33)       ; or fload_all_j(33)       ; or fload_halo_j(33)

radius= alog10(rad)

binsize = float((xmax-xmin))/bins

for i=1,bins do begin
  lg_r = i*binsize + xmin
  sm_r = (i-1)*binsize + xmin
  r_log(i-1) = 0.5*(lg_r + sm_r)

  idx= where((radius GE sm_r) AND (radius LT lg_r))
  if idx(0) lt 0 then begin
	m= [0.0]
	jdist(i-1)= 0
  endif else begin
	m= j(idx)
	if n_elements(idx) eq 1 then m=[m,m]
	minfo= moment(m)
	jdist(i-1)= minfo(0)
  endelse
endfor


; take out any zero's
; --------------------
idx= where(jdist GT 0)
r_log= r_log(idx)
jdist= jdist(idx)


ymin=min(jdist)
ymax=max(jdist)


;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        charsize=1.7, $
        xtitle="Log R [(h!E-1!Nkpc)]", $
        ytitle="J", $
        /nodata


        oplot, r_log, jdist, psym=-3, linestyle=0, color= 50

        ;x0= 0.65
        ;y0= 0.9
        ;xyouts, x0, y0, 'Dark Matter', /normal, charthick=1.7, size= 1.7, color= 50

	if keyword_set(msg) then begin
	        xyouts, 0.25, 0.80, msg, /normal, color= 0, size= 1.5
	endif



;--------------------------------------
;--------------------------------------

if (sendto EQ 'ps') then device, /close


end










