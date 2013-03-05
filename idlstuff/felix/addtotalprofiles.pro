PRO addtotalprofiles,xstars,ystars,xoldstars,yoldstars,xnewstars,ynewstars,xhalo,yhalo,$
                  totalprofiles,quantity,dir,sim,snapnum,maxusesnapshotnumber,reff,interpolationpoint,mergingtime,limitmtime,havetoscale,scalefactor

if (mergingtime gt 0) then begin
      ; totalprofiles            = fltarr(number,simnumber,7,quantity,dir)
      ; 0: 3D rho, 1: 3D beta, 2: 3D sigma 3: M/r, 4: 2D rho, 5: 2D sigma, 6: 2D kurtosis	    

      ; we use the SAME scalefactor for all components of course! The scale factor is computed from all stars
      
      xgt0        = where(xstars(*) gt 0)
      
      if (havetoscale eq 1) then begin
         scalefactor = interpol(ystars(xgt0),xstars(xgt0)/reff,[interpolationpoint,interpolationpoint])
         scalefactor = scalefactor(0)
      endif else begin
         scalefactor = 1.0
      endelse	 
      
      if (quantity eq 5) then begin
         ; sigma_p
         if (snapnum lt maxusesnapshotnumber) then begin	 
	     spawn,'echo "'+string(sim)+' '+string(snapnum)+' '+string(dir)+' '+string(scalefactor)+'" >>sigmascaling.txt'
	 endif    
      endif
    
      countarray  = fltarr(n_elements(totalprofiles(*,0,0,0,0,quantity,dir)))	    

      ; for all the simulations the average sigma2d is scaled to reff, and "1" at the interpolation point.
      result               = interpol(ystars(xgt0)/scalefactor,xstars(xgt0)/reff,totalprofiles(*,0,0,0,0,quantity,dir))	    

      ; set those indices that fall outside of the real region to zero	 
      select               = where( (totalprofiles(*,0,0,0,0,quantity,dir) lt min(xstars(xgt0)/reff)) or (totalprofiles(*,0,0,0,0,quantity,dir) gt max(xstars(xgt0)/reff)))
      result(select)       = 0.0
      countarray(*)        = 1.0
      countarray(select)   = 0.0

      totalprofiles(*,sim,1,snapnum,0,quantity,dir) = result
      totalprofiles(*,sim,1,snapnum,1,quantity,dir) = countarray

      if (mergingtime lt limitmtime) then begin
	 totalprofiles(*,sim,2,snapnum,0,quantity,dir) = result
	 totalprofiles(*,sim,2,snapnum,1,quantity,dir) = countarray	    	 
      endif else begin
	 totalprofiles(*,sim,3,snapnum,0,quantity,dir) = result
	 totalprofiles(*,sim,3,snapnum,1,quantity,dir) = countarray	
      endelse

      ; old stars	    
      xgt0                 = where(xoldstars(*) gt 0)      
      result               = interpol(yoldstars(xgt0)/scalefactor,xoldstars(xgt0)/reff,totalprofiles(*,0,0,0,0,quantity,dir))	    

      ; set those indices that fall outside of the real region to zero	 
      select               = where( (totalprofiles(*,0,0,0,0,quantity,dir) lt min(xoldstars(xgt0)/reff)) or (totalprofiles(*,0,0,0,0,quantity,dir) gt max(xoldstars(xgt0)/reff)))	    
      result(select)       = 0.0
      countarray(*)        = 1.0
      countarray(select)   = 0.0

      totalprofiles(*,sim,4,snapnum,0,quantity,dir) = result
      totalprofiles(*,sim,4,snapnum,1,quantity,dir) = countarray

      ; new stars	    
      xgt0                 = where (xnewstars(*) gt 0)
      result               = interpol(ynewstars(xgt0)/scalefactor,xnewstars(xgt0)/reff,totalprofiles(*,0,0,0,0,quantity,dir))	    

      ; set those indices that fall outside of the real region to zero	 
      select               = where( (totalprofiles(*,0,0,0,0,quantity,dir) lt min(xnewstars(xgt0)/reff)) or (totalprofiles(*,0,0,0,0,quantity,dir) gt max(xnewstars(xgt0)/reff)))	    	    
      result(select)       = 0.0
      countarray(*)        = 1.0
      countarray(select)   = 0.0

      totalprofiles(*,sim,5,snapnum,0,quantity,dir) = result
      totalprofiles(*,sim,5,snapnum,1,quantity,dir) = countarray

      ; halo
      xgt0                 = where (xhalo(*) gt 0) 
      result               = interpol(yhalo(xgt0)/scalefactor,xhalo(xgt0)/reff,totalprofiles(*,0,0,0,0,quantity,dir))	    

      ; set those indices that fall outside of the real region to zero	 
      select               = where( (totalprofiles(*,0,0,0,0,quantity,dir) lt min(xhalo(xgt0)/reff)) or (totalprofiles(*,0,0,0,0,quantity,dir) gt max(xhalo(xgt0)/reff)))
      result(select)       = 0.0
      countarray(*)        = 1.0
      countarray(select)   = 0.0

      totalprofiles(*,sim,6,snapnum,0,quantity,dir) = result
      totalprofiles(*,sim,6,snapnum,1,quantity,dir) = countarray
   endif
end	 
