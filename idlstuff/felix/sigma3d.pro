PRO sigma3d,ppp,vvv,index,xxx,yyy,zzz,$
                    xcenter,ycenter,zcenter,bin,mindist,maxdist,x,y,count

   n          = long(maxdist-mindist)/bin
   x          = fltarr(n)
   y          = fltarr(n)
   count      = fltarr(n)
   distance3D = sqrt((ppp(xxx,index)-xcenter)^2 + (ppp(yyy,index)-ycenter)^2 + (ppp(zzz,index)-zcenter)^2)   
   vel3D      = sqrt(vvv(0,index)^2 + vvv(1,index)^2 + vvv(2,index)^2)      
   ind        = long((distance3D(*)-mindist)/bin)

   for i=0L,n-1 do begin
     iii      = where(ind eq i)
     if (min(iii) ge 0) then begin
        if (n_elements(iii) gt 2) then begin
           x(i)     = mean(distance3D(iii))
           y(i)     = sqrt(variance(vel3D(iii)))
           count(i) = n_elements(iii)
         endif	
      endif
   endfor

   cleanup,x,y,count

END
