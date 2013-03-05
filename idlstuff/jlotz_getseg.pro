PRO GETSEG, imgfile, noisefile, segfile, positions, boxsize, thresh1, thresh2

; read in image 
  imgfile1 = imgfile + '.fits' 
  img = readfits(imgfile1)
 
  dimen = size(img, /dimensions)
  npix = dimen[0]
  seg = lonarr(npix, npix)
  seg2 = lonarr(npix, npix)

; get sigma noise
  noise = readfits(noisefile)
  sigma = robust_sigma(noise)

; smooth image

  simg = smooth(img, boxsize, /edge_truncate)

; get inital map

  max1 = simg[positions[0], positions[1]]
  max2 = simg[positions[2], positions[3]]
   dx = sqrt( (positions[0] - positions[2])^2 + (positions[1] - positions[3])^2)

  if (max1 le max2) then begin
      max = max1
  endif else begin
      max = max2
  endelse

  map = where(simg ge max)
  
; get segmap1

  ssimg = simg/sigma
  nmap = region_grow(ssimg, map, /all_neighbors, threshold=[thresh1, max(ssimg)])

  seg(nmap) = 10
  seg = sigma_filter(seg, boxsize, /ALL_PIXELS)

  seg = label_region(seg, /all_neighbors)
  nobj1 = max(seg)

  if (nobj1 gt 2) then begin
      area = fltarr(nobj1)
      for i=0, nobj1-1 do begin
          obj = where(seg eq i+1)
          area[i] = N_ELEMENTS(obj)
       endfor                   ;   

      sarea = sort(area)
      sarea = reverse(sarea)
 
      good = sarea[0:1] +1 ; 2 largest objects
      
      badobj = where(seg ne good[0] OR seg ne good[1]);
      obj1 = where(seg eq good[0])
      obj2 = where(seg eq good[1])

 
      seg(badobj) = 0
      seg(obj1) = 1
      seg(obj2) = 2
      
  endif         


; get segmap2
  
   nmap2 = where(ssimg ge thresh2)

   seg2(nmap2) = 10
   seg2 = label_region(seg2, /all_neighbors)
   nobj2 = max(seg2)
   if (nobj2 gt 2) then begin
      area = fltarr(nobj2)
      for i=0, nobj2-1 do begin
          obj = where(seg2 eq i+1)
          area[i] = N_ELEMENTS(obj)
       endfor                   ;   

      sarea = sort(area)
      sarea = reverse(sarea)
 
      good = sarea[0:1] +1 ; 2 largest objects
      
      badobj = where(seg2 ne good[0] OR seg2 ne good[1]);
      obj1 = where(seg2 eq good[0])
      obj2 = where(seg2 eq good[1])

 
      seg2(badobj) = 0
      seg2(obj1) = 1
      seg2(obj2) = 2
      
  endif         
 
; get final segmap

   if (nobj1 eq 1 AND nobj2 eq 1) then begin
       fits_write, segfile, seg2

   endif else begin
       good = where(seg gt 0)
       seg_all = seg
       seg_all(good) = 1.0

       box1 = smooth(float(seg), dx,/edge_truncate)
       box2 = smooth(float(seg_all), dx, /edge_truncate) 

       id = round(box1 /box2)
       bad = where(box2 eq 0.0)
       id(bad) = 0
       seg3 = seg

       new = where(seg2 gt 0 AND seg lt 1)
       if (new[0] ne -1) then begin
           seg3(new) = id(new)

           i = 4
           while (min(id(new)) eq 0  AND i gt 1) do begin
  
               good = where(seg gt 0)
               seg_all = seg
               seg_all(good) = 1.0

               s = npix/i
               box1 = smooth(float(seg), s,/edge_truncate)
               box2 = smooth(float(seg_all), s, /edge_truncate) 

               id = round(box1 /box2)
               bad = where(box2 eq 0.0)
               id(bad) = 0

               new = where(seg2 gt 0 AND seg lt 1 AND seg3 lt 1)
               seg3(new) = id(new)
               i = i-1          

           endwhile
       endif


       if (max(seg3) gt 1) then begin
           obj1 = where(seg3 eq 1)
           obj2 = where(seg3 eq 2)

           obj1map = seg3
           obj1map(obj2) = 0

           obj2map = seg3
           obj2map(obj1) = 0

           obj1map = label_region(obj1map)
           obj2map = label_region(obj2map)
 
           if (max(obj1map) gt 1) then begin
               area = fltarr(max(obj1map))
               for i=0, (max(obj1map)-1) do begin
                   obj = where(obj1map eq i+1)
                   area[i] = N_ELEMENTS(obj)
               endfor                  

               sarea = sort(area)
               sarea = reverse(sarea)
 
               good = sarea[0] +1 ; largest section of obj1
               bad = sarea[1:*] +1 ; other sections      
           
               obj1map_grow = dilate(obj1map, replicate(1, 5, 5))
               obj1map_grow = label_region(obj1map_grow)
 
               for i=0, N_ELEMENTS(bad)-1 do begin              
                   section = where(obj1map ne bad[i] AND obj1map_grow eq bad[i] AND seg3 ne 0)
                   if (N_ELEMENTS(section) gt 1) then begin
                       badobj = where(obj1map eq bad[i])
                       newid = round(mean(seg3(section)))
                       seg3(badobj) = newid
                   endif

               endfor
           endif
   
           if (max(obj2map) gt 1) then begin
               area = fltarr(max(obj2map))
               for i=0, (max(obj2map)-1) do begin
                   obj = where(obj2map eq i+1)
                   area[i] = N_ELEMENTS(obj)
               endfor                  
               
               sarea = sort(area)
               sarea = reverse(sarea)
 
               good = sarea[0] +1 ; largest section of obj1
               bad = sarea[1:*] +1 ; other sections      
           
               obj2map_grow = dilate(obj2map, replicate(1, 5, 5))
               obj2map_grow = label_region(obj2map_grow)
 
               for i=0, N_ELEMENTS(bad)-1 do begin              
                   section = where(obj2map ne bad[i] AND obj2map_grow eq bad[i] AND seg3 ne 0)
                   if (N_ELEMENTS(section) gt 1) then begin
                       badobj = where(obj2map eq bad[i])
                       newid = round(mean(seg3(section)))
                       seg3(badobj) = newid
                   endif
               endfor
           endif
           
           positions = round(positions)
           if ( seg3[positions[0],positions[1]] ne 1) then begin
               obj1 = where(seg3 eq 2)
               obj2 = where(seg3 eq 1)
               seg3(obj2) = 2
               seg3(obj1) = 1
           endif

       endif

       fits_write, segfile, seg3

   endelse

   
end


  
