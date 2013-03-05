PRO sigma2d,ppp,vvv,index,xxx,yyy,zzz,$
                    xcenter,ycenter,zcenter,bin,mindist,maxdist,x,y,count

   n          = long(maxdist-mindist)/bin
   x          = fltarr(n)
   y          = fltarr(n)
   count      = fltarr(n)
   distance2D = sqrt((ppp(xxx,index)-xcenter)^2 + (ppp(yyy,index)-ycenter)^2)
   vel2D      = vvv(zzz,index)
   ind        = long((distance2D(*)-mindist)/bin)
   
   for i=0L,n-1 do begin
     iii      = where(ind eq i)
     IF (n_elements(iii) ge 2) THEN BEGIN   	
        x(i)     = mean(distance2D(iii))
        y(i)     = sqrt(variance(vel2D(iii)))
        count(i) = n_elements(iii)
     ENDIF
   endfor

   cleanup,x,y,count
     
END
