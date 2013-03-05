PRO masstolightratio,starspos,starsmass,halopos,halomass,xxx,yyy,zzz,xcenter,ycenter,zcenter,radius,mtol

   
   distance3Dstars    = sqrt((starspos(xxx,*)-xcenter)^2 + (starspos(yyy,*)-ycenter)^2 + (starspos(zzz,*)-zcenter)^2)
   distance3Dhalo     = sqrt((halopos(xxx,*)-xcenter)^2  + (halopos(yyy,*)-ycenter)^2  + (halopos(zzz,*)-zcenter)^2)

   halovirialmass     = total(halomass (where(distance3Dhalo  lt radius)),/double)
   starsvirialmass    = total(starsmass(where(distance3Dstars lt radius)),/double)
   
   mtol = (halovirialmass+starsvirialmass)/starsvirialmass

END
