pro interpolate_gas_data, Time, Frun, Snap, NumList, TiList, N, Pos, Vel, Hsml, Mass, Temp, Id

N = 0L

ind = sort(TiList)
TiList=TiList(ind)
NumList=NumList(ind)

if (time ge min(TiList)) and (Time le max(TiList)) then begin

    indx = 0
    
    while TiList(indx+1)  lt Time do begin
        indx++
    endwhile

    print, "Indx=", indx, "  Snap=", NumList(Indx)

    load_gas_data, NumList(indx), Frun, Snap, N1, Time1, Pos1, Vel1, Hsml1, Mass1, Temp1, Id1

    load_gas_data, NumList(indx+1), Frun, Snap, N2, Time2, Pos2, Vel2, Hsml2, Mass2, Temp2, Id2

    print, Time1, N1
    print, Time2, N2

    ind1 = sort(Id1)
    ind2 = sort(Id2)

    i1 = 0L
    i2 = 0L
    N= 0L

    while (i1 lt n1) and (i2 lt n2) do begin
        if Id1(ind1(i1)) eq Id2(ind2(i2)) then begin
            N++
            i1++
            i2++
        endif else begin
            if Id1(ind1(i1)) lt Id2(ind2(i2)) then begin
                i1++
            endif else begin
                i2++
            endelse
        endelse
    endwhile

    print, "N= ", N


    Pos= fltarr(3,N)
    Vel= fltarr(3,N)
    Id = lonarr(N)
    Mass=fltarr(N)	
    Hsml=fltarr(N)
    Temp=fltarr(N)


    i1 = 0L
    i2 = 0L
    i= 0L

    t0 = Time     
    dt = Time2 - Time1
    t =  Time - Time1
    t2 = t * t                  
    t3 = t * t * t              

    if dt gt 0 then begin
        dt2inv = 1.0 / ((dt * dt) / 2) 
        dt2invb = 1.0 / (3 * dt * dt)
        linfac = t/dt 
    endif else begin
        dt2inv = 0
        dt2invb = 0
        linfac = 0
    endelse


    print, "Time1=", Time1, "  Time=", Time, "  Time2=", Time2,"  dt=", dt

    while (i1 lt n1) and (i2 lt n2) do begin
        if Id1(ind1(i1)) eq Id2(ind2(i2)) then begin

            x0 = Pos1(*,ind1(i1))
            v0 = Vel1(*,ind1(i1))
            
            dx = Pos2(*,ind2(i2)) - x0
            dv = Vel2(*,ind2(i2)) - v0
            
            c = (3 * dx - (3 * v0 + dv) * dt) * dt2inv
            f = (dv - dt * c) * dt2invb

;            ff = sqrt(f(0)^2 + f(1)^2 + f(2)^2)
;            cc = sqrt(c(0)^2 + c(1)^2 + c(2)^2)

;            if (ff * dt) lt (0.5*cc) then begin
;                Pos(*,i) = x0 + v0 * t + 0.5 * t2 * c + t3 * f
;                Vel(*,i) = v0 + t * c + 3 * t2 * f
;            endif else begin
                Pos(*,i) = x0 + dx * linfac
                Vel(*,i) = v0 + dv * linfac
;            endelse

            Hsml(i) =  Hsml1(ind1(i1)) + (Hsml2(ind2(i2)) - Hsml1(ind1(i1))) * linfac
            Mass(i) =  Mass1(ind1(i1)) + (Mass2(ind2(i2)) - Mass1(ind1(i1))) * linfac
            Temp(i) =  exp(alog(Temp1(ind1(i1))) + (alog(Temp2(ind2(i2))) - alog(Temp1(ind1(i1)))) * linfac)

            id(i) = id1(ind1(i1))


            i++
            i1++
            i2++
        endif else begin
            if Id1(ind1(i1)) lt Id2(ind2(i2)) then begin
                i1++
            endif else begin
                i2++
            endelse
        endelse
    endwhile

endif

end



