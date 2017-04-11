
Base = "../test1"
SnapBase = "snap"

BoxSize = 1.0

window, xsize= 900, ysize= 900

for num = 0, 160,1 do begin

    exts='000'
    exts=exts+strcompress(string(Num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)
    
    f= Base + "/" + Snapbase + "_"+exts 
    f= Strcompress(f, /remove_all)

    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedback=0L
    npartall=lonarr(6)	
    bytesleft= 136
    la=intarr(bytesleft/2)


    openr,1,f,/f77_unformatted
    readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,la
    print,npart,massarr
    print,time,redshift
    
    NGas=  npart(0)
    
    N=  NGas
    
    pos=fltarr(3,N)
    vel=fltarr(3,N)
    id=lonarr(N)
    readu,1, pos
    readu,1, vel
    readu,1, id
    ind=where((npart gt 0) and (massarr eq 0))
    if ind(0) ne -1 then begin
        Nm= total(npart(ind))
        mass=fltarr(Nm) ; masses for variable mass particles (usually gas+stars)
        readu,1,mass
    endif

    u=fltarr(Ngas)
    readu,1,u                   ; internal energy per unit mass
    rho=fltarr(Ngas)
    readu,1,rho                 ; comoving gas density
    
    if flag_sfr gt 0 then begin
        Nelec=fltarr(Ngas)
        readu,1 ,Nelec     ; gas electron abundance relative to hydrogen
        NH0=fltarr(Ngas)
        readu,1 ,NH0   ; neutral hydrogen abundance relative to hydrogen
    endif

    hsml=fltarr(Ngas)
    readu,1,hsml                ; smoothing length

    close,1



    f= Base + "/voronoi_mesh_"+exts 
    f= Strcompress(f, /remove_all)


    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedback=0L
    npartall=lonarr(6)	
    bytesleft= 136
    la=intarr(bytesleft/2)

    ngas = 0L
    nel  = 0L
    nedgepoints = 0L

    openr,1,f
    readu,1,ngas,nel, nedgepoints
    
    nedges = lonarr(ngas)
    nedges_offset = lonarr(ngas)
    edgelist = lonarr(nel)
    points = fltarr(2, nedgepoints)

    readu,1, nedges
    readu,1, nedges_offset
    readu,1, edgelist
    readu,1, points
    close,1





    if num eq 0 then begin
        mi = 0.0
        ma = 2.0
        loadct, 3
        tvlct, r, g, b, /get
    endif

    
    colindex= (rho - mi)/(ma-mi)*255.0
    ind = where(colindex ge 256.0)
    if ind(0) ne -1 then colindex(ind) = 255.9
    ind = where(colindex lt 0)
    if ind(0) ne -1 then colindex(ind) = 0
    colindex = byte(colindex)


    plot, [0], [0], psym=3, xrange=[0,BoxSize], yrange=[0,BoxSize], xstyle=1, ystyle=1

    for i=0, ngas-1 do begin

        x = transpose(points(0, edgelist(nedges_offset(i):nedges_offset(i)+nedges(i)-1)))
        y = transpose(points(1, edgelist(nedges_offset(i):nedges_offset(i)+nedges(i)-1)))

        polyfill, [x, x(0)], [y, y(0)], noclip=0, color=r(colindex(i))+g(colindex(i))*256L+b(colindex(i))*256L^2
        plots, [x, x(0)], [y, y(0)], noclip=0

        ind = where((x lt 0) or (x gt BoxSize) or (y lt 0) or (y gt BoxSize))
        if ind(0) ne -1 then begin
            for dx=-1,1 do begin
                for dy =-1,1 do begin
                    polyfill, [x, x(0)]+dx*BoxSize, [y, y(0)]+dy*BoxSize, noclip=0, color=r(colindex(i))+g(colindex(i))*256L+b(colindex(i))*256L^2
                    plots, [x, x(0)]+dx*BoxSize, [y, y(0)]+dy*BoxSize, noclip=0
                endfor
            endfor
        endif
    endfor

    plot, pos(0,*), pos(1,*), psym=3, xrange=[0,BoxSize], yrange=[0,BoxSize], xstyle=1, ystyle=1, /noerase


    wait, 0.5

 
endfor

end





 









