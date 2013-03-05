PRO betaprofile,ppp,vvv,index,xxx,yyy,zzz,$
                    xcenter,ycenter,zcenter,bin,maxdist,x,vrady,vphiy,vthetay,vtany
   ; Note: this is in 3D as compared to the 2D


   distance = sqrt((ppp(xxx,index)-xcenter)^2 + (ppp(yyy,index)-ycenter)^2 + (ppp(zzz,index)-zcenter)^2)      
   
   ; get the average velocity and sigma of the selected stars
   velocity         = fltarr(3,n_elements(index));
   for i=0,2 do begin
      velocity(i,*)  = vvv(i,index) - mean(vvv(i,index),/double);
   endfor   
   
   ; print,mean(velocity(0,*),/double),mean(velocity(1,*),/double),mean(velocity(2,*),/double);
   ; CHECK: OK
   
   
   vrad    = ( velocity(xxx,*)*(ppp(xxx,index)-xcenter) + $
               velocity(yyy,*)*(ppp(yyy,index)-ycenter) + $
               velocity(zzz,*)*(ppp(zzz,index)-zcenter)) / distance(*)
   
   vphi    = ( velocity(xxx,*) * (ppp(yyy,index)-ycenter) - $
               velocity(yyy,*) * (ppp(xxx,index)-xcenter) ) / sqrt((ppp(xxx,index)-xcenter)^2 + (ppp(yyy,index)-ycenter)^2)
   
   vtheta  = (- velocity(xxx,*) * (ppp(xxx,index)-xcenter) * (ppp(zzz,index)-zcenter)        - $
                velocity(yyy,*) * (ppp(yyy,index)-ycenter) * (ppp(zzz,index)-zcenter)        + $
	        velocity(zzz,*) * ((ppp(yyy,index)-ycenter)^2 + (ppp(xxx,index)-xcenter)^2)) / $
		sqrt( (ppp(xxx,index)-xcenter)^2 * (ppp(zzz,index)-zcenter)^2                + $
		      (ppp(yyy,index)-ycenter)^2 * (ppp(zzz,index)-zcenter)^2                + $
                     ((ppp(yyy,index)-ycenter)^2 + (ppp(xxx,index)-xcenter)^2)^2)

   vtan = sqrt(velocity(xxx,*)^2+velocity(yyy,*)^2+velocity(zzz,*)^2-vrad^2)
		     
   
   one     = distance
   one(*)  = 1.0   
   vrady   =                hist1d(distance,vrad*vrad    ,obin=x,binsize=bin,max=maxdist);
   vrady   = sqrt(vrady   / hist1d(distance,one          ,obin=x,binsize=bin,max=maxdist))
   
   vphiy   =                hist1d(distance,vphi*vphi    ,obin=x,binsize=bin,max=maxdist);
   vphiy   = sqrt(vphiy   / hist1d(distance,one          ,obin=x,binsize=bin,max=maxdist))
     
   vthetay =                hist1d(distance,vtheta*vtheta,obin=x,binsize=bin,max=maxdist);
   vthetay = sqrt(vthetay / hist1d(distance,one          ,obin=x,binsize=bin,max=maxdist))

   vtany   =                hist1d(distance,vtan*vtan      ,obin=x,binsize=bin,max=maxdist);
   vtany   = sqrt(vtany   / hist1d(distance,one            ,obin=x,binsize=bin,max=maxdist))

END
