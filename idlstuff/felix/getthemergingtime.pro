PRO getthemergingtime,simulation,num,mergingradius,mergingtime
   ; read the centers text
   ; use the mergingradius to get the merging time
   ; read the snapshot header -> time
   ; return mergingtime

   readcenters,simulation+"/centers.txt",distance,time
   
   timeofmerger = interpol(time,distance,[mergingradius,mergingradius])
   timeofmerger = timeofmerger(0)

   print,timeofmerger        
   print,'timeofmerger',timeofmerger
   close,/all
   
   exts = '000'
   exts = exts+strcompress(string(num),/remove_all)
   exts = strmid(exts,strlen(exts)-3,3)
   f    = simulation+"snapshot_"+exts
   f    = strcompress(f,/remove_all)

   npart          = lonarr(6)
   massarr        = dblarr(6)
   simtime        = 0.0D
   redshift       = 0.0D
   flag_sfr       = 0L
   flag_feedback  = 0L

   bytesleft      = 256-6*4 - 6*8 - 8 - 8 - 2*4

   la             = intarr(bytesleft/2)

   print,f
   openr,1,f,/f77_unformatted

   readu,1,npart,massarr,simtime,redshift,flag_sfr,flag_feedback,la
   close,1
   
   mergingtime = simtime - timeofmerger
   print,'simtime',simtime
   print,'mergingtime',mergingtime

   ;window,2
   ;plot,time,distance,yrange=[0,20]
   ;oplot,[simtime,simtime],[-1,100],color=250
   ;oplot,[timeofmerger,timeofmerger],[-1,100],color=150
   ;oplot,[mergingtime,mergingtime],[-1,100],color=50
   
   ;a=''
   ;read,a
   ;window,0

END
