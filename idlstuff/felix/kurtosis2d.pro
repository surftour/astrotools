PRO kurtosis2d,ppp,vvv,index,xxx,yyy,zzz,$
                    xcenter,ycenter,zcenter,bin,mindist,maxdist,x,y,count

   n          = long(maxdist-mindist)/bin
   x          = fltarr(n)
   y          = fltarr(n)
   count      = fltarr(n)
   distance2D = sqrt((ppp(xxx,index)-xcenter)^2 + (ppp(yyy,index)-ycenter)^2)
   vel2D      = vvv(zzz,*)
   ind        = long((distance2D(*)-mindist)/bin)
   
   for i=0L,n-1 do begin
     iii      = where(ind eq i)
     if (n_elements(iii) ge 2) then begin
        x(i)     = mean(distance2D(iii))
        y(i)     = kurtosis(vel2D(iii))
        count(i) = n_elements(iii)
     endif	
   endfor
   
   cleanup,x,y,count
     
END
