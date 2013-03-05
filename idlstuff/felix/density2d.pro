PRO density2d,ppp,vvv,mass,index,xxx,yyy,zzz,$
                    xcenter,ycenter,zcenter,bin,mindist,maxdist,x,y,count

   n          = long(maxdist-mindist)/bin
   x          = fltarr(n)
   y          = fltarr(n)
   count      = fltarr(n)
   distance2D = sqrt((ppp(xxx,index)-xcenter)^2 + (ppp(yyy,index)-ycenter)^2)
   ind        = long((distance2D(*)-mindist)/bin)
   
   for i=0L,n-1 do begin
     iii      = where(ind eq i)
     IF (min(iii) ge 0) THEN BEGIN
        x(i)     = mean(distance2D(iii))
        y(i)     = total(mass(iii)) / !Pi / ((mindist+(i+1.0)*bin)^2-(mindist+i*bin)^2) / bin^2
        count(i) = n_elements(iii)
     ENDIF
   endfor

   cleanup,x,y,count

END
