; fn='/raid5/yxli/Runs/merger_tree/third_runs/Gas/h13t939_0.5rd_1r200_mbh_z/eos0.5/z_8.49_0/z_8.49_0_104'

pro load_center_and_BHs, fname, bhlist, time, cm, center1, nbh, idbh, mbh, posbh

npart=lonarr(6)	
massarr=dblarr(6)
time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedbacktp=0L
npartTotal=lonarr(6)	
flag_cooling=0L
numfiles=0L
cosmocrap=dblarr(4)
flag_stellarage=0L
flag_metals=0L
bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4 - 2*4 - 4*8 - 2*4
la=intarr(bytesleft/2)

fname=strcompress(fname,/remove_all)
openr,1,fname,/f77_unformatted

readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedbacktp,npartTotal, $
  flag_cooling,numfiles,cosmocrap,flag_stellarage,flag_metals,la


; volker isn't using a flag, so do it manually
if keyword_set(nopot_in_snap) then flag_snaphaspot= 0L else flag_snaphaspot= 1L
;flag_snaphaspot= 0L
flag_stargens= 1L
flag_energydetail= 0L
flag_parentid= 0L
flag_starorig= 0L


N=total(npart)
pos=fltarr(3,N)
vel=fltarr(3,N)
id=lonarr(N)

ind=where((npart gt 0) and (massarr eq 0)) 
if ind(0) ne -1 then begin
    Nwithmass= total(npart(ind))
    mass=fltarr(Nwithmass)
endif else begin	
    Nwithmass= 0
endelse

readu,1,pos
readu,1,vel
readu,1,id
if Nwithmass gt 0 then begin
    readu,1,mass
    print, "Nwithmasses= ",Nwithmass
endif

Ngas=npart(0)
Nhalo=npart(1)
Ndisk=npart(2)
Nbulge=npart(3)
Nstars=npart(4)
Nbh=npart(5)

print, 'N, ngas, nhalo, nstars, nbh', N, ngas, nhalo, nstars, nbh
if Ngas gt 0 then begin
    u=fltarr(Ngas)
    readu,1,u
endif

if Ngas gt 0 then begin
    rho=fltarr(Ngas)
    readu,1,rho
endif

if flag_cooling gt 0 and Ngas gt 0 then begin
    nume=fltarr(Ngas)
    numh=fltarr(Ngas)
    readu,1,nume
    readu,1,numh
endif

if Ngas gt 0 then begin
    hsml=fltarr(Ngas)
    readu,1,hsml
endif

if flag_sfr gt 0 and Ngas gt 0 then begin
    sfr=fltarr(Ngas)
    readu,1,sfr
    print, "sfr= ", total(sfr)
endif

if flag_sfr gt 0 then begin
    if flag_stellarage gt 0 and Nstars gt 0 then begin
        stellage= fltarr(Nstars)
        readu,1,stellage
;        print, 'mean stellage', mean(stellage)
    endif
endif

if flag_metals gt 0 then begin
    if (Ngas+Nstars) gt 0 then begin
        mets= fltarr(Ngas+Nstars)
        readu,1,mets
        gmetals= fltarr(Ngas)
        gmetals(*)= mets(0:Ngas-1)
        if Nstars gt 0 then begin
            smetals= fltarr(Nstars)
            smetals(*)= mets(Ngas:Ngas+Nstars-1)
;            print, 'mean gmetals, smetals', mean(gmetals), mean(smetals)
        endif
    endif
endif

if flag_snaphaspot gt 0 then begin
    pot= fltarr(N)
    readu,1,pot
endif

close, 1

cosmos=0
hubble=0.7
;    if keyword_set(h) then begin
;	hubble= fload_cosmology('h')
if (cosmos) then begin
    print, "snapshot values corrected for h= ", hubble
    
    time= time/hubble
    pos= pos/hubble
    mass= mass/hubble
    massarr= massarr/hubble
    rho=rho*hubble*hubble
    mfs=mfs/hubble
    if flag_multiphase gt 0 then mclouds=mclouds/hubble
    hsml=hsml/hubble
    if flag_stellarage gt 0 and Nstars gt 0 then stellage= stellage/hubble
    if flag_metals gt 0 then gmetals=gmetals/hubble
    if flag_metals gt 0 and Nstars gt 0 then smetals=smetals/hubble
    if flag_starorig gt 0 and Nstars gt 0 then origmass=origmass/hubble
    if flag_starorig gt 0 and Nstars gt 0 then orighsml=orighsml/hubble
endif

cm=fltarr(3)
cm[0]=total(pos[0,*]*mass[*])/total(mass)
cm[1]=total(pos[1,*]*mass[*])/total(mass)
cm[2]=total(pos[2,*]*mass[*])/total(mass)

if Nstars gt 0 then begin               
    n1=Nhalo+Ngas+Ndisk+Nbulge
    n2=n1+Nstars
    idstars=intarr(Nstars)
    idstars=id(n1:n2-1)
    s=0L
    while (s lt Nstars) do begin
        if idstars[s] gt (2UL)^31 then idstars[s]=idstars[s]-(2UL)^31
        s=s+1
    endwhile
    
    id(n1:n2-1)=idstars(*)
endif

center1=fltarr(3)
w=where(id ge 1L and id le bhlist[0], c)
if c gt 0 then begin
    center1[0]=total(pos[0, w]*mass[w])/total(mass[w])
    center1[1]=total(pos[1, w]*mass[w])/total(mass[w])
    center1[2]=total(pos[2, w]*mass[w])/total(mass[w])
endif
   
print, 'gal center= ', cm
    
;;;

    ; ----------------------------------
    ;   now start parsing file
    ; ----------------------------------

if Nbh gt 0 then begin
    xbh=fltarr(Nbh) &  ybh=fltarr(Nbh)  & zbh=fltarr(Nbh)  & mbh=fltarr(Nbh) & idbh=intarr(Nbh)
    xbh(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
    ybh(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
    zbh(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
    
    vxbh=fltarr(Nbh) &  vybh=fltarr(Nbh) &  vzbh=fltarr(Nbh)
    vxbh(*)=vel(0,Nhalo+Ndisk+Ngas+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
    vybh(*)=vel(1,Nhalo+Ndisk+Ngas+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
    vzbh(*)=vel(2,Nhalo+Ndisk+Ngas+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
    
    if massarr(5) eq 0 then begin
        skip=0L
        for t=0,4 do begin
            if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                skip=skip + npart(t)
            endif
        endfor
        print, 'read from mass'
        mbh(*)=mass(0+skip:Nbh-1+skip)
    endif else begin
        print, 'read from massarr'
        mbh(*)= massarr(5)
    endelse
    
    if flag_snaphaspot gt 0 then begin
        pbh=fltarr(Nbh)
        pbh(*)=pot(Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
    endif

    idbh=id(Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)

    print, 'nbh, idbh', nbh, idbh
;    print, 'mbh', mbh
;    print, 'xbh, ybh, zbh', xbh, ybh, zbh

    w=where(mbh eq max(mbh), c)
    cbx=xbh[w[0]]
    cby=ybh[w[0]]
    cbz=zbh[w[0]]

;    print, 'BH center', cbx, cby, cbz

    posbh=fltarr(3, nbh)
    posbh[0,*]=xbh
    posbh[1,*]=ybh
    posbh[2,*]=zbh
    
endif

if Ngas gt 0 then begin
    xgas=fltarr(Ngas) &  ygas=fltarr(Ngas)  & zgas=fltarr(Ngas) & mgas=fltarr(Ngas) & idgas=intarr(Ngas)
    xgas(*)=pos(0,0:Ngas-1)
    ygas(*)=pos(1,0:Ngas-1)
    zgas(*)=pos(2,0:Ngas-1)
    
    vxgas=fltarr(Ngas) &  vygas=fltarr(Ngas)  & vzgas=fltarr(Ngas)
    vxgas(*)=vel(0,0:Ngas-1)
    vygas(*)=vel(1,0:Ngas-1)
    vzgas(*)=vel(2,0:Ngas-1)
    
    if massarr(0) eq 0 then begin
        mgas(*)=mass(0:Ngas-1)	
    endif else begin
        mgas(*)= massarr(0)
    endelse
    
    idgas=id(0:Ngas-1)
    
    if flag_snaphaspot gt 0 then begin
        pgas=fltarr(Ngas)
        pgas(*)=pot(0:Ngas-1)
    endif

endif

if Nhalo gt 0 then begin
    xhalo=fltarr(Nhalo) &  yhalo=fltarr(Nhalo) & zhalo=fltarr(Nhalo) & mhalo=fltarr(Nhalo) & idhalo=intarr(Nhalo)
    xhalo(*)=pos(0,0+Ngas:Nhalo+Ngas-1)
    yhalo(*)=pos(1,0+Ngas:Nhalo+Ngas-1)
    zhalo(*)=pos(2,0+Ngas:Nhalo+Ngas-1)
    
    vxhalo=fltarr(Nhalo) &  vyhalo=fltarr(Nhalo) & vzhalo=fltarr(Nhalo)
    vxhalo(*)=vel(0,0+Ngas:Nhalo+Ngas-1)
    vyhalo(*)=vel(1,0+Ngas:Nhalo+Ngas-1)
    vzhalo(*)=vel(2,0+Ngas:Nhalo+Ngas-1)
    
    if massarr(1) eq 0 then begin
        skip=0L
        for t=0,0 do begin	
            if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                skip=skip + npart(t)
            endif
        endfor
        mhalo(*)=mass(0+skip:Nhalo-1+skip)	
    endif else begin
        mhalo(*)= massarr(1)
    endelse
    
    if flag_snaphaspot gt 0 then begin
        phalo=fltarr(Nhalo)
        phalo(*)=pot(0+Ngas:Nhalo+Ngas-1)
    endif
    
    idhalo=id(0+Ngas:Nhalo+Ngas-1)
    
endif

 ; Disk
 ; ------
mdisk= 0
if Ndisk gt 0 then begin
    xdisk=fltarr(Ndisk) &  ydisk=fltarr(Ndisk) &  zdisk=fltarr(Ndisk) & mdisk=fltarr(Ndisk)
    xdisk(*)=pos(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
    ydisk(*)=pos(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
    zdisk(*)=pos(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
    
    vxdisk=fltarr(Ndisk) &  vydisk=fltarr(Ndisk) &  vzdisk=fltarr(Ndisk)
    vxdisk(*)=vel(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
    vydisk(*)=vel(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
    vzdisk(*)=vel(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
    
    if massarr(2) eq 0 then begin
        skip=0L
        for t=0,1 do begin	
            if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                skip=skip + npart(t)
            endif
        endfor
        mdisk(*)=mass(0+skip:Ndisk-1+skip)	
    endif else begin
        mdisk(*)= massarr(2)
    endelse
    
    if flag_snaphaspot gt 0 then begin
        pdisk=fltarr(Ndisk)
        pdisk(*)=pot(Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
    endif
    
    diskage= -5.0*randomu(seed,Ndisk)
 ;diskage= -1.0*randomu(seed,Ndisk)    ; changed this for Marijn
    diskmetals= 10^(1.5*randomu(seed,Ndisk) - 3.0)
    
endif

 ; Bulge
 ; -------
mbulge= 0
if Nbulge gt 0 then begin
    xbulge=fltarr(Nbulge) &  ybulge=fltarr(Nbulge) &  zbulge=fltarr(Nbulge) & mbulge=fltarr(Nbulge)
    xbulge(*)=pos(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
    ybulge(*)=pos(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
    zbulge(*)=pos(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
    
    vxbulge=fltarr(Nbulge) &  vybulge=fltarr(Nbulge) &  vzbulge=fltarr(Nbulge)
    vxbulge(*)=vel(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
    vybulge(*)=vel(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
    vzbulge(*)=vel(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
    
    if massarr(3) eq 0 then begin
        skip=0L
        for t=0,2 do begin	
            if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                skip=skip + npart(t)
            endif
        endfor
        mbulge(*)=mass(0+skip:Nbulge-1+skip)	
    endif else begin
        mbulge(*)= massarr(3)
    endelse
    
    if flag_snaphaspot gt 0 then begin
        pbulge=fltarr(Nbulge)
        pbulge(*)=pot(Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
    endif
    
    bulgeage= -7.0 + 0.0*mbulge
    bulgemetals= 0.001 + 0.0*mbulge
endif


if Nstars gt 0 then begin
    xstars=fltarr(Nstars) &  ystars=fltarr(Nstars)  & zstars=fltarr(Nstars)  & $
      mstars=fltarr(Nstars)  & idstars=intarr(Nstars)

    xstars(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
    ystars(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
    zstars(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
    
    vxstars=fltarr(Nstars) &  vystars=fltarr(Nstars) &  vzstars=fltarr(Nstars)
    vxstars(*)=vel(0,Nhalo+Ndisk+Ngas+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
    vystars(*)=vel(1,Nhalo+Ndisk+Ngas+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
    vzstars(*)=vel(2,Nhalo+Ndisk+Ngas+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
    
    if massarr(4) eq 0 then begin
        skip=0L
        for t=0,3 do begin	
            if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                skip=skip + npart(t)
            endif
        endfor
        mstars(*)=mass(0+skip:Nstars-1+skip)	
    endif else begin
        mstars(*)= massarr(4)
    endelse
    
    if flag_snaphaspot gt 0 then begin
        pstars=fltarr(Nstars)
        pstars(*)=pot(Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
    endif

    idstars=id(Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)

    s=0L
    while (s lt Nstars) do begin
        if idstars[s] gt (2UL)^31 then idstars[s]=idstars[s]-(2UL)^31
        s=s+1
    endwhile

endif


;    com=fload_center(1)
;    alternate_com= fload_center_felix(1)


end


