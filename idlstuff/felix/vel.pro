;  This code reads the simulations of TJ Cox and returns the plots of astro-ph/0501622

;  Felix 2005, stoehr@iap.fr

;  The data should be stored in the directory given in basename and have the same structure as the
;  directory of 
;
;      http://physics.ucsc.edu/~tj/work/felix/
;
;  In order to reproduce the plots of the paper however, only the MAJOR-MERGER simulations should be
;  in the directory basename:
;  
;  G2r/ G3/ G3b/ G3r/ Sbc201a-u4/ Sbc_hfb/ Sbc_o1/ Sbc_o2/ Sbc_o3/ Sc/
;
;  If ps = 1 then the .ps files of the paper as well as some corresponding .gif files are created
;  If ps = 0 then a directory Results is created in the current directory and filled with
;            the plots for the different simulations and outputs.

   basename  = "/home/tcox/felixtest/"

   ps           = 1
   forcecompute = 1

   charactersize = 2.5
   linethickness = 5.0
   
   !P.TICKLEN=0.04
   !P.THICK= linethickness
   !X.MARGIN = [7,1]
   !Y.MARGIN = [4,1]   
   
   maxsnapshotnumber         = 5
   maxusesnapshotnumber      = 2
        
   interpolationpoint2Dsigma = 1.0
   interpolationpointrest    = 1.0
   
   mergingradius             = 5.0 
   limitmtime                = 0.86071376
  
   maxdist                   = 48.0; kpc NECESSARY SO THAT NO ADDITIONAL BIN IS MADE
   mindist                   = 0.0;
      
   binnumber                 = 24;

   ; the stacked interpolated values (in units of reff)
   number                    = 40.0
   maxreffective             = 40.0;
   minreffective             = 0.1
      
   ; the additional factor is necessary in order to get the correct number of bins   
   bin      = (maxdist-mindist)/(binnumber+1)*1.000001 ; kpc
      
   spawn,'/bin/ls -1 '+strcompress(basename)+'>simulations.dat'
   directories = strarr(1000)
   openr,1,'simulations.dat'
   i = 0L
   stringvar = ''
   while (not eof(1)) do begin   
      readf,1,stringvar
      directories(i) = stringvar
      i = i + 1
   endwhile
   directories = directories(0:i-1)
   close,1
   simnumber = i
   print,' Found ',simnumber,' simulations'
   print,directories  
      
   setmycolors

   spawn,'rm sigmascaling.txt'
   spawn,'rm reff.txt'
   
   spawn,'mkdir Results'
  
IF ((n_elements(totalprofiles) eq 0) or (forcecompute eq 1) ) THEN BEGIN

   totalprofiles                = fltarr(number,simnumber,7,maxsnapshotnumber,2,7,3)
   totalprofiles(*,*,*,*,*,*,*)  = 0.0            
   FOR quantity = 0,6 DO BEGIN
      FOR dir=0,2 DO BEGIN    
         ; Log spaced interpolation values
         totalprofiles(*,0,0,0,0,quantity,dir)  = 10.0^(findgen(number)/double(number) * (alog10(maxreffective)-alog10(minreffective)) + alog10(minreffective))
      ENDFOR
   ENDFOR   
      
   averagemtime             = 0.0
   averagemtimecount        = 0.0

   averagevirialradius      = 0.0
   averagevirialMtoL        = 0.0   
   averagevirialradiuscount = 0.0

   averagereffective        = 0.0
   averagereffectivecount   = 0.0
   ;delvarx,averagereffectivelist

   averagescalefactorrho2d  = 0.0
   averagescalefactorsigma2d= 0.0   

   averagescalefactorrho2dcount   = 0.0
   averagescalefactorsigma2dcount = 0.0   
   
   window,xsize=400,ysize=400,retain=2;   

   FOR sim = 0,simnumber-1 DO BEGIN

      spawn,"/bin/ls -1 "+strcompress(basename)+"/"+directories(sim)+"/snap* | tr \/ \\n | grep snapshot | grep -v snapshotdata | sed 's/snapshot_//g'>snapshots.dat"
      snapshots      = read_ascii('snapshots.dat')   
      snapshots      = long(round(snapshots.field1(*)))
      snapshotnumber = n_elements(snapshots)

      print," Searching directory ",directories(sim)
      print," Found ",snapshotnumber," snapshots (",snapshots,")";

      ; this is to assure that all plots have the same range. A reff of 4 kpc is typical.
      minreffdist = mindist / 4.0
      maxreffdist = maxdist / 4.0

      for snapnum = 0,snapshotnumber-1 do begin

	 num = snapshots(snapshotnumber-1-snapnum)
	 spawn,"rm *.gif"
         
	 ; this is the starting direction 
	 dir = 0

	 pos=0L
	 vel=0L
	 simulation = basename + "/"+strcompress(directories(sim),/remove_all)+"/"            
	 readsnapshot,simulation,num,N,npart,massarr,time,redshift,pos,vel,mass

         centerpositionfile = simulation+"centerposition."+strcompress(string(num),/remove_all)+".dat"
	  	 
	 ; this reads the centerposition.txt file if it exists. Otherwise, it computes the centre and creates the file	 
	 findthecentreposition,pos,mass,num,N,center,centerpositionfile
	 getnewcenter,center,dir,num,xxx,yyy,zzz,xcenter,ycenter,zcenter
	 print,xcenter,ycenter,zcenter
 	 	  
         ;getthemergingtime,simulation,num,mergingradius,mergingtime
	 mergingtime= 1.4
	 print,simulation,mergingtime

	 selecttype,10,npart,pos,vel,mass,outstarspos   ,outstarsvel   ,outstarsmass   ,outnstars
	 selecttype,11,npart,pos,vel,mass,outoldstarspos,outoldstarsvel,outoldstarsmass,outnoldstars      
	 selecttype,4 ,npart,pos,vel,mass,outnewstarspos,outnewstarsvel,outnewstarsmass,outnnewstars
	 selecttype,1 ,npart,pos,vel,mass,outhalopos    ,outhalovel    ,outhalomass    ,outnhalo            

         virialradius,pos,mass,xxx,yyy,zzz,xcenter,ycenter,zcenter,radius
	 masstolightratio,outstarspos,outstarsmass,outhalopos,outhalomass,xxx,yyy,zzz,xcenter,ycenter,zcenter,radius,mtol
	 
	 averagevirialradius      = averagevirialradius      + radius
	 averagevirialMtoL        = averagevirialMtol        + mtol
	 averagevirialradiuscount = averagevirialradiuscount + 1.0
	 print,"Rvir, mtol at Rvir",radius,mtol

	 ; ============================================================================
	 ; 2D
	 ; ============================================================================

	 ; Stack the three projections together

	 ; 0: stars, 1: old stars, 2 new stars, 3 halo (dm)      
	 xtot         = fltarr(binnumber,4)	 
	 ytotdensity  = fltarr(binnumber,4)
	 ytotsigma    = fltarr(binnumber,4)      
	 ytotkurtosis = fltarr(binnumber,4)
	 counttot     = fltarr(binnumber,4)	             

	 print,"WHAT ABOUT THE X_AXIS? IS TAKING THE LAST ONE GOOD ENOUGH????"

	 meanreff = 0.0

	 for dir=0,2 do begin
            print,"Doing direction ",dir," ..."
            getnewcenter,center,dir,num,xxx,yyy,zzz,xcenter,ycenter,zcenter	        

	    squaredistance = (outstarspos(xxx,*)-xcenter)^2 + $
                	     (outstarspos(yyy,*)-ycenter)^2

            res = sort(squaredistance)
	    squaredistance = squaredistance(res)
	    sortedmass     = outstarsmass(res)

	    i               = 1L
	    halfmass        = 0.5*total(outstarsmass)      
	    cumulativemass  = sortedmass(0)
	    while (cumulativemass lt halfmass) do begin
               cumulativemass = cumulativemass + sortedmass(i)
               i = i + 1	 
	    endwhile
	    reff    = sqrt(squaredistance(i))
	    print,'reff:',reff
	    if (snapnum lt maxusesnapshotnumber) then begin
	       spawn,'echo '+string(reff)+' >> reff.txt'
	    endif
	    meanreff = meanreff + reff

	    image=tvrd(true=1);
	    write_bmp,"a.bmp",image,/rgb
	    spawnstring = 'convert a.bmp image.'+strcompress(string(num),/remove_all)+"."+strcompress(string(dir),/remove_all)+'.gif'
	    print,spawnstring
	    spawn,spawnstring
	    	    
	    ; DENSITY 2D PROFILES
	    ; ---------------------------------------------------------------------------------------

            density2d,outstarspos,outstarsvel,outstarsmass,lindgen(outnstars),xxx,yyy,zzz,$
                	  xcenter,ycenter,zcenter,bin,mindist,maxdist,xstars,ystars,countstars

            density2d,outoldstarspos,outoldstarsvel,outoldstarsmass,lindgen(outnoldstars),xxx,yyy,zzz,$
                	  xcenter,ycenter,zcenter,bin,mindist,maxdist,xoldstars,yoldstars,countoldstars

            density2d,outnewstarspos,outnewstarsvel,outnewstarsmass,lindgen(outnnewstars),xxx,yyy,zzz,$
                	  xcenter,ycenter,zcenter,bin,mindist,maxdist,xnewstars,ynewstars,countnewstars

	    density2d,outhalopos,outhalovel,outhalomass,lindgen(outnhalo),xxx,yyy,zzz,$
                	  xcenter,ycenter,zcenter,bin,mindist,maxdist,xhalo,yhalo,counthalo	

            addtotalprofiles,xstars,   ystars,$
	                  xoldstars,yoldstars,$
		          xnewstars,ynewstars,$
		          xhalo,    yhalo,$
                          totalprofiles,4,dir,sim,snapnum,maxusesnapshotnumber,reff,interpolationpointrest,mergingtime,limitmtime,1,scalefactor

   	    averagescalefactorrho2d      = averagescalefactorrho2d       + 1.0/scalefactor
	    averagescalefactorrho2dcount = averagescalefactorrho2dcount  + 1.0

            ytotdensity(*,0) = ytotdensity(*,0) + ystars
            ytotdensity(*,1) = ytotdensity(*,1) + yoldstars
            ytotdensity(*,2) = ytotdensity(*,2) + ynewstars	 	 
            ytotdensity(*,3) = ytotdensity(*,3) + yhalo

	    ; SIGMA 2D PROFILES
	    ; ---------------------------------------------------------------------------------------

	    sigma2d,outstarspos,outstarsvel,lindgen(outnstars),xxx,yyy,zzz,$
                	  xcenter,ycenter,zcenter,bin,mindist,maxdist,xstars,ystars,countstars

	    sigma2d,outoldstarspos,outoldstarsvel,lindgen(outnoldstars),xxx,yyy,zzz,$
                	  xcenter,ycenter,zcenter,bin,mindist,maxdist,xoldstars,yoldstars,countoldstars

	    sigma2d,outnewstarspos,outnewstarsvel,lindgen(outnnewstars),xxx,yyy,zzz,$
                	  xcenter,ycenter,zcenter,bin,mindist,maxdist,xnewstars,ynewstars,countnewstars

	    sigma2d,outhalopos,outhalovel,lindgen(outnhalo),xxx,yyy,zzz,$
                	  xcenter,ycenter,zcenter,bin,mindist,maxdist,xhalo,yhalo,counthalo
	    
            addtotalprofiles,xstars,   ystars,$
	                  xoldstars,yoldstars,$
		          xnewstars,ynewstars,$
		          xhalo,    yhalo,$
                          totalprofiles,5,dir,sim,snapnum,maxusesnapshotnumber,reff,interpolationpointrest,mergingtime,limitmtime,1,scalefactor

            ytotsigma(*,0) = ytotsigma(*,0) + ystars
	    ytotsigma(*,1) = ytotsigma(*,1) + yoldstars
	    ytotsigma(*,2) = ytotsigma(*,2) + ynewstars	 
            ytotsigma(*,3) = ytotsigma(*,3) + yhalo
      
	    ; KURTOSIS 2D PROFILES
	    ; ---------------------------------------------------------------------------------------

	    kurtosis2d,outstarspos,outstarsvel,lindgen(outnstars),xxx,yyy,zzz,$
                	  xcenter,ycenter,zcenter,bin,mindist,maxdist,xstars,ystars,countstars

	    kurtosis2d,outoldstarspos,outoldstarsvel,lindgen(outnoldstars),xxx,yyy,zzz,$
                	  xcenter,ycenter,zcenter,bin,mindist,maxdist,xoldstars,yoldstars,countoldstars

	    kurtosis2d,outnewstarspos,outnewstarsvel,lindgen(outnnewstars),xxx,yyy,zzz,$
                	  xcenter,ycenter,zcenter,bin,mindist,maxdist,xnewstars,ynewstars,countnewstars

	    kurtosis2d,outhalopos,outhalovel,lindgen(outnhalo),xxx,yyy,zzz,$
                	  xcenter,ycenter,zcenter,bin,mindist,maxdist,xhalo,yhalo,counthalo

            addtotalprofiles,xstars,   ystars,$
	                  xoldstars,yoldstars,$
		          xnewstars,ynewstars,$
		          xhalo,    yhalo,$
                          totalprofiles,6,dir,sim,snapnum,maxusesnapshotnumber,reff,interpolationpointrest,mergingtime,limitmtime,0,scalefactor

            ytotkurtosis(*,0) = ytotkurtosis(*,0) + ystars
	    ytotkurtosis(*,1) = ytotkurtosis(*,1) + yoldstars
	    ytotkurtosis(*,2) = ytotkurtosis(*,2) + ynewstars	 	 
            ytotkurtosis(*,3) = ytotkurtosis(*,3) + yhalo

            ; all 2D xvalues are the same so we can take this one
	    ; however, we have to know how many good values are used each time
            xtot(*,0) = xtot(*,0) + xstars
	    xtot(*,1) = xtot(*,1) + xoldstars
	    xtot(*,2) = xtot(*,2) + xnewstars	 
	    xtot(*,3) = xtot(*,3) + xhalo	    
            
	    countthesestars    = fltarr(n_elements(countstars))
	    counttheseoldstars = fltarr(n_elements(countoldstars))
	    countthesenewstars = fltarr(n_elements(countnewstars))
	    countthesehalo     = fltarr(n_elements(counthalo))

	    countthesestars(where(countstars gt 0))       = 1.0
	    counttheseoldstars(where(countoldstars gt 0)) = 1.0
	    countthesenewstars(where(countnewstars gt 0)) = 1.0
	    countthesehalo(where(counthalo gt 0))         = 1.0
	 	    
	    counttot(*,0) = counttot(*,0) + countthesestars
	    counttot(*,1) = counttot(*,1) + counttheseoldstars
	    counttot(*,2) = counttot(*,2) + countthesenewstars	 
	    counttot(*,3) = counttot(*,3) + countthesehalo

	 endfor

	 ; normalize the densities, sigmas and kurtosises
	 mynormalize,xtot,        counttot
         mynormalize,ytotdensity, counttot
	 mynormalize,ytotsigma,   counttot
	 mynormalize,ytotkurtosis,counttot

         ; this effective radius is also used for the 3D profile stacking
	 meanreff               = meanreff     / 3.0
         averagereffective      = averagereffective      + meanreff
	 averagereffectivecount = averagereffectivecount + 1
	 if (n_elements(averagereffectivelist) eq 0) then begin
	    averagereffectivelist    = fltarr(1)
	    averagereffectivelist(0) = meanreff
	 endif else begin
	    dummylist = averagereffectivelist
	    averagereffectivelist = fltarr(n_elements(dummylist)+1)
	    averagereffectivelist(0:n_elements(dummylist)-1) = dummylist
	    averagereffectivelist(n_elements(dummylist)) = meanreff	    
	 endelse   	 


         ; write out the data for the fiducial galaxy: all stars, new stars. Projected (all projections) v_p and 
	 ; r/meanrreff
         writeoutstars,'Sbc201a-u4',60,directories(sim),num,center,$
	                              outstarspos,outstarsvel,outnewstarspos,outnewstarsvel,$
              	                      meanreff	    	 

	 plot,[0],[0],xtitle="r [r_eff="+strcompress(string(meanreff),/remove_all)+" kpc], mtime="+strcompress(string(mergingtime),/remove_all)+"Gyr",yrange=[1e-7,1],xrange=[minreffdist,maxreffdist],xstyle=1,ytitle='rho2D(r) 1e10 Msun/kpc^2',/ylog;

	 oplot,xtot(*,0)/meanreff,ytotdensity(*,0) ,color=1;	    
	 oplot,xtot(*,1)/meanreff,ytotdensity(*,1) ,color=1,linestyle=1;	    
	 oplot,xtot(*,2)/meanreff,ytotdensity(*,2) ,color=1,linestyle=2;	                
	 oplot,xtot(*,3)/meanreff,ytotdensity(*,3) ,color=2;

	 image=tvrd(true=1);
	 write_bmp,"a.bmp",image,/rgb
	 spawnstring = 'convert a.bmp surface.gif'
	 print,spawnstring
	 spawn,spawnstring
	 plot,[0],[0],xtitle="r [r_eff="+strcompress(string(meanreff),/remove_all)+" kpc], mtime="+strcompress(string(mergingtime),/remove_all)+"Gyr",yrange=[0,300],xrange=[minreffdist,maxreffdist],xstyle=1,ytitle='sigma2D(r) [km/s]',ystyle=1

         overplotromanowsky,ytotsigma(*,0),xtot(*,0)/meanreff,interpolationpoint2Dsigma

	 oplot,xtot(*,0)/meanreff,ytotsigma(*,0) ,color=1;
	 oplot,xtot(*,1)/meanreff,ytotsigma(*,1) ,color=1,linestyle=1;
	 oplot,xtot(*,2)/meanreff,ytotsigma(*,2) ,color=1,linestyle=2;            
	 oplot,xtot(*,3)/meanreff,ytotsigma(*,3) ,color=2;

	 image=tvrd(true=1);
	 write_bmp,"a.bmp",image,/rgb
	 spawnstring = 'convert a.bmp sigma2d.gif'
	 print,spawnstring
	 spawn,spawnstring

	 averagescalefactorsigma2d      = averagescalefactorsigma2d       + 1.0/scalefactor
	 averagescalefactorsigma2dcount = averagescalefactorsigma2dcount  + 1.0

	 plot,[0],[0],xtitle="r [r_eff="+strcompress(string(meanreff),/remove_all)+" kpc], mtime="+strcompress(string(mergingtime),/remove_all)+"Gyr",yrange=[-1,2],xrange=[minreffdist,maxreffdist],xstyle=1,ytitle='kurtosis(r)';

	 oplot,xtot(*,0)/meanreff,ytotkurtosis(*,0) ,color=1;
	 oplot,xtot(*,1)/meanreff,ytotkurtosis(*,1) ,color=1,linestyle=1;
	 oplot,xtot(*,2)/meanreff,ytotkurtosis(*,2) ,color=1,linestyle=2;            
	 oplot,xtot(*,3)/meanreff,ytotkurtosis(*,3) ,color=2;
	 oplot,[-10,1000],[0,0],linestyle=1;

	 image=tvrd(true=1);
	 write_bmp,"a.bmp",image,/rgb
	 spawnstring = 'convert a.bmp kurtosis.gif'
	 print,spawnstring
	 spawn,spawnstring
	 
	 ; ============================================================================
	 ; 3D
	 ; ============================================================================
         
	 
	 ; 3D DENSITY PROFILES
	 ; ----------------------------------------------------------------------------
	 density3d,outstarspos,outstarsvel,outstarsmass,lindgen(outnstars),xxx,yyy,zzz,$
                       xcenter,ycenter,zcenter,bin,mindist,maxdist,xstars,ystars,countstars

	 density3d,outoldstarspos,outoldstarsvel,outoldstarsmass,lindgen(outnoldstars),xxx,yyy,zzz,$
                       xcenter,ycenter,zcenter,bin,mindist,maxdist,xoldstars,yoldstars,countoldstars

	 density3d,outnewstarspos,outnewstarsvel,outnewstarsmass,lindgen(outnnewstars),xxx,yyy,zzz,$
                       xcenter,ycenter,zcenter,bin,mindist,maxdist,xnewstars,ynewstars,countnewstars

	 density3d,outhalopos,outhalovel,outhalomass,lindgen(outnhalo),xxx,yyy,zzz,$
                       xcenter,ycenter,zcenter,bin,mindist,maxdist,xhalo,yhalo,counthalo

         FOR dir=0,2 DO BEGIN
            addtotalprofiles,xstars,ystars,$
	                  xoldstars,yoldstars,$
		          xnewstars,ynewstars,$
		          xhalo,yhalo,$
                          totalprofiles,0,dir,sim,snapnum,maxusesnapshotnumber,meanreff,interpolationpointrest,mergingtime,limitmtime,1,scalefactor
	 ENDFOR	       


	 plot,[0],[0],xtitle="r [kpc]",yrange=[1e-9,0.1],xrange=[mindist,maxdist],ytitle='rho(r) 1e10 Msun/kpc^3',/ylog;
	 oplot,xstars,ystars,color=1;
	 oplot,xhalo,yhalo,color=2;

	 oplot,xoldstars,yoldstars,color=1,linestyle=1;
	 oplot,xnewstars,ynewstars,color=1,linestyle=2;

	 image=tvrd(true=1);
	 write_bmp,"a.bmp",image,/rgb
	 spawnstring = 'convert a.bmp rho.gif'
	 print,spawnstring
	 spawn,spawnstring


	 savedensitystars = ystars
	 savedensityhalo  = yhalo

	 ; BETA PROFILES
	 ; ----------------------------------------------------------------------------

	 betaprofile,outstarspos,outstarsvel,lindgen(outnstars),xxx,yyy,zzz,$
        	 xcenter,ycenter,zcenter,bin,maxdist,xstars,vradystars,vphiystars,vthetaystars,vtanystars

	 betaprofile,outoldstarspos,outoldstarsvel,lindgen(outnoldstars),xxx,yyy,zzz,$
        	 xcenter,ycenter,zcenter,bin,maxdist,xoldstars,vradyoldstars,vphiyoldstars,vthetayoldstars,vtanyoldstars

	 betaprofile,outnewstarspos,outnewstarsvel,lindgen(outnnewstars),xxx,yyy,zzz,$
        	 xcenter,ycenter,zcenter,bin,maxdist,xnewstars,vradynewstars,vphiynewstars,vthetaynewstars,vtanynewstars

	 betaprofile,outhalopos,outhalovel,lindgen(outnhalo),xxx,yyy,zzz,$
        	 xcenter,ycenter,zcenter,bin,maxdist,xhalo,vradyhalo,vphiyhalo,vthetayhalo,vtanyhalo

         betastars    = 1.0-vtanystars    * vtanystars    / 2.0 / vradystars    / vradystars
	 betaoldstars = 1.0-vtanyoldstars * vtanyoldstars / 2.0 / vradyoldstars / vradyoldstars
	 betanewstars = 1.0-vtanynewstars * vtanynewstars / 2.0 / vradynewstars / vradynewstars
	 betahalo     = 1.0-vtanyhalo     * vtanyhalo     / 2.0 / vradyhalo     / vradyhalo

         FOR dir=0,2 DO BEGIN 
            addtotalprofiles,xstars,betastars,$
	                  xoldstars,betaoldstars,$
		          xnewstars,betanewstars,$
		          xhalo,betahalo,$
                          totalprofiles,1,dir,sim,snapnum,maxusesnapshotnumber,meanreff,interpolationpointrest,mergingtime,limitmtime,0,scalefactor
	 ENDFOR		  

	 plot,[0],[0],xtitle="r [kpc]",yrange=[-2,1],xrange=[mindist,maxdist],ytitle='beta(r)';

	 oplot,xstars,betastars  ,color=1;	    
	 oplot,xhalo,betahalo    ,color=2;	    

	 oplot,xoldstars,betaoldstars  ,color=1,linestyle=1;
	 oplot,xnewstars,betanewstars  ,color=1,linestyle=2;	    

	 oplot,[-10,1000],[0,0],linestyle=1;

	 image=tvrd(true=1);
	 write_bmp,"a.bmp",image,/rgb
	 spawnstring = 'convert a.bmp beta.gif'
	 print,spawnstring
	 spawn,spawnstring

	 savebetastars = betastars
	 savebetahalo  = betahalo


	 ; SIGMA PROFILES

	 sigma3d,outstarspos,outstarsvel,lindgen(outnstars),xxx,yyy,zzz,$
                       xcenter,ycenter,zcenter,bin,mindist,maxdist,xstars,ystars,countstars

	 sigma3d,outoldstarspos,outoldstarsvel,lindgen(outnoldstars),xxx,yyy,zzz,$
                       xcenter,ycenter,zcenter,bin,mindist,maxdist,xoldstars,yoldstars,countoldstars

	 sigma3d,outnewstarspos,outnewstarsvel,lindgen(outnnewstars),xxx,yyy,zzz,$
                       xcenter,ycenter,zcenter,bin,mindist,maxdist,xnewstars,ynewstars,countnewstars

	 sigma3d,outhalopos,outhalovel,lindgen(outnhalo),xxx,yyy,zzz,$
                       xcenter,ycenter,zcenter,bin,mindist,maxdist,xhalo,yhalo,counthalo

         FOR dir=0,2 DO BEGIN
            addtotalprofiles,xstars,ystars,$
	                  xoldstars,yoldstars,$
		          xnewstars,ynewstars,$
		          xhalo,yhalo,$
                          totalprofiles,2,dir,sim,snapnum,maxusesnapshotnumber,meanreff,interpolationpointrest,mergingtime,limitmtime,1,scalefactor
	 ENDFOR		  

	 plot,[0],[0],xtitle="r [kpc]",yrange=[10,250],xrange=[mindist,maxdist],ytitle='sigma3D(r) [km/s]',/ylog,ystyle=1;

	 oplot,xstars,ystars ,color=1;	    
	 oplot,xhalo,yhalo   ,color=2;	    

	 oplot,xoldstars,yoldstars ,color=1,linestyle=1;
	 oplot,xnewstars,ynewstars ,color=1,linestyle=2;

	 image=tvrd(true=1);
	 write_bmp,"a.bmp",image,/rgb
	 spawnstring = 'convert a.bmp sigma3d.gif'
	 print,spawnstring
	 spawn,spawnstring

         ; M profiles 
	 
	 m3d,outstarspos,outstarsvel,outstarsmass,lindgen(outnstars),xxx,yyy,zzz,$
                       xcenter,ycenter,zcenter,bin,mindist,maxdist,xstars,ystars,countstars

	 m3d,outoldstarspos,outoldstarsvel,outoldstarsmass,lindgen(outnoldstars),xxx,yyy,zzz,$
                       xcenter,ycenter,zcenter,bin,mindist,maxdist,xoldstars,yoldstars,countoldstars

	 m3d,outnewstarspos,outnewstarsvel,outnewstarsmass,lindgen(outnnewstars),xxx,yyy,zzz,$
                       xcenter,ycenter,zcenter,bin,mindist,maxdist,xnewstars,ynewstars,countnewstars

	 m3d,outhalopos,outhalovel,outhalomass,lindgen(outnhalo),xxx,yyy,zzz,$
                       xcenter,ycenter,zcenter,bin,mindist,maxdist,xhalo,yhalo,counthalo

         FOR dir=0,2 DO BEGIN
            addtotalprofiles,xstars,ystars,$
	                  xoldstars,yoldstars,$
		          xnewstars,ynewstars,$
		          xhalo,yhalo,$
                          totalprofiles,3,dir,sim,snapnum,maxusesnapshotnumber,meanreff,interpolationpointrest,mergingtime,limitmtime,1,scalefactor
	 ENDFOR		  

	 plot,[0],[0],xtitle="r [kpc]",yrange=[0.1,100],xrange=[mindist,maxdist],ytitle='M [1e10 Msun]',/ylog;

	 oplot,xstars,ystars ,color=1;
	 oplot,xhalo,yhalo   ,color=2;	    

	 oplot,xoldstars,yoldstars ,color=1,linestyle=1;
	 oplot,xnewstars,ynewstars ,color=1,linestyle=2;
	 oplot,xnewstars,ynewstars+yoldstars ,color=0,linestyle=1,psym=2;
	 	 
	 image=tvrd(true=1);
	 write_bmp,"a.bmp",image,/rgb
	 spawnstring = 'convert a.bmp moverr.gif'
	 print,spawnstring
	 spawn,spawnstring

	 savesigmastars = ystars
	 savesigmahalo  = yhalo

         if (mergingtime gt 0) then begin
	    averagemtime       = averagemtime      + mergingtime
            averagemtimecount  = averagemtimecount + 1.0
	 endif
        
         if (ps eq 0) then begin 
            spawn,'mkdir Results/'+strcompress(directories(sim),/remove_all)
            spawn,'mkdir Results/'+strcompress(directories(sim),/remove_all)+"/"+strcompress(string(num),/remove_all)
	    spawn,'rm a.bmp; mv *.gif Results/'+strcompress(directories(sim),/remove_all)+"/"+strcompress(string(num),/remove_all)
      	 endif      
      endfor ; this is the endloop for each snapshot!!!      	      	            
            
   ENDFOR
   averagemtime          = averagemtime           / averagemtimecount
   averagevirialradius   = averagevirialradius    / averagevirialradiuscount
   averagevirialMtoL     = averagevirialMtoL      / averagevirialradiuscount   
   averagescalefactorsigma2d = averagescalefactorsigma2d  / averagescalefactorsigma2dcount
   averagescalefactorrho2d   = averagescalefactorrho2d    / averagescalefactorrho2dcount   
   averagevirialMtoL     = averagevirialMtoL      / averagevirialradiuscount         
   averagereffective     = averagereffective      / averagereffectivecount

ENDIF

   print,'----------------------------------------------'
   print,'averagemtime       ',averagemtime   
   print,'averagevirialMtoL  ',averagevirialMtoL
   print,'averagevirialradius',averagevirialradius
   print,'averagereffective  ',averagereffective
   print,'averagescalefactorsigma2d',averagescalefactorsigma2d
   print,'averagescalefactorrho2d  ',averagescalefactorrho2d
   ;print,'averagereffectivelist, mean, sigma',mean(averagereffectivelist),sqrt(variance(averagereffectivelist))
      
   print,'Reff/Rvir          ',averagereffective/averagevirialradius

   totalprofilesaverage                = fltarr(number,7,3,7)
   totalprofilesaverage(*,*,*,*)        = 0.0

   FOR quantity=0,6 DO BEGIN
      totalprofilesaverage(*,0,0,quantity)          = totalprofiles(*,0,0,0,0,0) ; x values   
      FOR i=1L,6 DO BEGIN ; the first one are the positions
	 FOR j=0L,simnumber-1 DO BEGIN
            FOR k=0,maxusesnapshotnumber-1 DO BEGIN
	       FOR dir=0,2 DO BEGIN
                  totalprofilesaverage(*,i,0,quantity) = totalprofilesaverage(*,i,0,quantity) + totalprofiles(*,j,i,k,0,quantity,dir)
	          totalprofilesaverage(*,i,1,quantity) = totalprofilesaverage(*,i,1,quantity) + totalprofiles(*,j,i,k,1,quantity,dir)	 
	       ENDFOR
	    ENDFOR   
	 ENDFOR
	 select                                  = where(totalprofilesaverage(*,i,1,quantity) gt 0.0)
	 IF (max(select gt 0)) THEN BEGIN
  	    totalprofilesaverage(select,i,0,quantity)  = totalprofilesaverage(select,i,0,quantity)/totalprofilesaverage(select,i,1,quantity)
	 ENDIF   
      ENDFOR    
   ENDFOR   

   FOR quantity=0,6 DO BEGIN
      FOR i=1L,6 DO BEGIN ; the first one are the positions
	 FOR j=0L,simnumber-1 DO BEGIN
            FOR k=0,maxusesnapshotnumber-1 DO BEGIN
	       FOR dir=0,2 DO BEGIN
	          ; only select the snapshots that REALLY are used. Otherwise the value is zero leading to huge sigmas
                  select = where(totalprofiles(*,j,i,k,0,quantity,dir) gt 0.0)
	          IF (max(select) gt 0) THEN BEGIN
        	     totalprofilesaverage(select,i,2,quantity)   = totalprofilesaverage(select,i,2,quantity) + (totalprofiles(select,j,i,k,0,quantity,dir) - totalprofilesaverage(select,i,0,quantity))^2
	          ENDIF	    
	       ENDFOR	  
	    ENDFOR   
	 ENDFOR
         select                                 = where(totalprofilesaverage(*,i,1,quantity) gt 0.0)
	 IF (max(select gt 0)) THEN BEGIN	 	    
	    totalprofilesaverage(select,i,2,quantity)  = sqrt(totalprofilesaverage(select,i,2,quantity)/totalprofilesaverage(select,i,1,quantity))
	 ENDIF   
      ENDFOR    
   ENDFOR   

   ; maxxplot = 6.75
   maxxplot = 10.0
   minxplot = 0.39

   maxxplotsigma = 6.0
   minxplotsigma = 0.0
   plotshift     = 0.16
   
   IF (ps eq 1) THEN BEGIN
      set_plot,'PS'
      device,/tt_font,filename='sigma2dstack.eps',/color,XSIZE=20.0,YSIZE=15.0,/encapsulated
   ENDIF   
   
   ; STACKED PLOTS: FINAL 2DSigma
   ; -----------------------------------------------------------------------------------------------------
   plot,[0],[0],xrange=[minxplotsigma,maxxplotsigma],yrange=[0.4,1.6],ystyle=1,xstyle=1,xtitle="r_p/R_{eff}",ytitle="\sigma_p/\sigma_{p,eff}",$
       charsize=charactersize,xthick=linethickness,ythick=linethickness,charthick=linethickness
   ; yrange=0.35,1.5
   ;;;select = where( (totalprofilesaverage(*,0,0,5) gt minxplotsigma) and (totalprofiles(*,0,0,0,0,5) lt maxxplotsigma) and (totalprofilesaverage(*,6,0,5) gt 0))      
   select = where( (totalprofilesaverage(*,0,0,5) gt minxplotsigma+plotshift) and (totalprofilesaverage(*,6,0,5) gt 0))      
   last = max(select)

;   polyfill, [totalprofilesaverage(select,0,0,5),totalprofilesaverage(last,0,0,5),reverse(totalprofilesaverage(select,0,0,5))],$
;             [totalprofilesaverage(select,6,0,5)+totalprofilesaverage(select,6,2,5), $
;	      totalprofilesaverage(last,6,0,5)+totalprofilesaverage(last,6,2,5),$ 
;	     reverse(totalprofilesaverage(select,6,0,5)-totalprofilesaverage(select,6,2,5))],color=5,NOCLIP=0

   polyfill, [totalprofilesaverage(select,0,0,5),totalprofilesaverage(last,0,0,5),reverse(totalprofilesaverage(select,0,0,5))],$
            [totalprofilesaverage(select,5,0,5)+totalprofilesaverage(select,5,2,5), $
	      totalprofilesaverage(last,5,0,5)+totalprofilesaverage(last,5,2,5),$ 
	     reverse(totalprofilesaverage(select,5,0,5)-totalprofilesaverage(select,5,2,5))],color=4,NOCLIP=0

;  polyfill, [totalprofilesaverage(select,0,0,5),totalprofilesaverage(last,0,0,5),reverse(totalprofilesaverage(select,0,0,5))],$
;             [totalprofilesaverage(select,1,0,5)+totalprofilesaverage(select,1,2,5), $
;      totalprofilesaverage(last,1,0,5)+totalprofilesaverage(last,1,2,5),$ 
;     reverse(totalprofilesaverage(select,1,0,5)-totalprofilesaverage(select,1,2,5))],color=10


  ; oplot, [totalprofilesaverage(select,0,0,5)],[totalprofilesaverage(select,6,0,5)-totalprofilesaverage(select,6,2,5)]$
;	     ,color=5,NOCLIP=0,thick=2.0*linethickness ;color=103;== cyan

 ;  oplot, [totalprofilesaverage(select,0,0,5)],[totalprofilesaverage(select,6,0,5)+totalprofilesaverage(select,6,2,5)]$
;	     ,color=5,NOCLIP=0,thick=2.0*linethickness ;color=103;== cyan

;   overplotromanowsky,totalprofilesaverage(*,1,0,5),totalprofilesaverage(*,0,0,5),interpolationpoint2Dsigma,thick=linethickness
   ; select only the non-zero regions

   fitromanowsky,totalprofilesaverage(*,1,0,5),totalprofilesaverage(*,0,0,5),totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,1,0,5),thick=linethickness
   ; select only the non-zero regions


;;   plotmyerror,4.875,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,5,0,5),totalprofilesaverage(select,5,2,5),color=1,thick=2.0*linethickness,linestyle=2
   plotmyerror,5.125,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,1,0,5),totalprofilesaverage(select,1,2,5),color=1,thick=2.0*linethickness
   plotmyerror,5.375,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,6,0,5),totalprofilesaverage(select,6,2,5),color=2,thick=2.0*linethickness
   ;plotmyerror,5.325,totalprofilesaverage(*,6,0,5),totalprofilesaverage(pixtouse,6,2,5),color=2,thick=1.0*linethickness
   ;plotmyerror,5.450,totalprofilesaverage(*,1,0,5),totalprofilesaverage(pixtouse,1,2,5),color=1,thick=1.0*linethickness

   plot,[0],[0],xrange=[minxplotsigma,maxxplotsigma],yrange=[0.4,1.6],ystyle=1,xstyle=1,charsize=charactersize,/noerase,xthick=linethickness,ythick=linethickness

   oplot,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,6,0,5),color=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,5,0,5),color=1,linestyle=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,4,0,5),color=1,linestyle=1
   ;oplot,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,3,0,5),color=205,thick=2.0*linethickness ; late
   ;oplot,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,2,0,5),color=4,thick=2.0*linethickness ; early
   oplot,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,1,0,5),color=1,thick=2.0*linethickness
   
   ;oplot,totalprofilesaverage(select,0,0),totalprofilesaverage(select,1,0)+totalprofilesaverage(select,1,2),color=1,thick=4.0*linethickness
   ;oplot,totalprofilesaverage(select,0,0),totalprofilesaverage(select,1,0)-totalprofilesaverage(select,1,2),color=1,thick=4.0*linethickness  

   ; totalprofiles: 0 x ; 1 sigma all ; 2 sigma early ; 3 sigma late ; 4 sigma old ; 5 sigma new   6 sigmadm
   ;
   ; ,,2: 0 sigma, 1 count
   ;
   ; 0: 3D rho, 1: 3D beta, 2: 3D sigma, 3: M, 4: 2D rho, 5: 2D sigma, 6: 2D kurtosis
   
   ;quantity = 5
   ;FOR i=5,5 DO BEGIN 
   ;   FOR j=0L,simnumber-1 DO BEGIN
   ;      FOR k=0,maxusesnapshotnumber-1 DO BEGIN
   ;	    FOR dir=0,2 DO BEGIN
   ;   	       oplot,totalprofiles(*,0,0,0,0,quantity,dir),totalprofiles(*,j,i,k,0,quantity,dir),color=1
   ;	    ENDFOR 
   ;      	 ENDFOR   
   ;        ENDFOR      
   ;   ENDFOR    

   IF (ps ne 1) THEN BEGIN   
      image=tvrd(true=1);
      write_bmp,"a.bmp",image,/rgb
      spawnstring = 'convert a.bmp sigma2dstack.all.gif'
      print,spawnstring
      spawn,spawnstring
      spawn,'rm a.bmp; mv *.gif Results/'
   ENDIF ELSE BEGIN
      device,/close
      set_plot,'x'
   ENDELSE  
   
   
   IF (ps eq 1) THEN BEGIN
      set_plot,'PS'
      device,/tt_font,filename='sigma2dstack.fiducial.eps',/color,XSIZE=20.0,YSIZE=15.0,/encapsulated
   ENDIF   
   
   
; STACKED PLOTS: FINAL 2DSigma FIDUCIAL
   ; -----------------------------------------------------------------------------------------------------
   plot,[0],[0],xrange=[minxplotsigma,maxxplotsigma],yrange=[0.4,1.6],ystyle=1,xstyle=1,xtitle="r_p/R_{eff}",ytitle="\sigma_p/\sigma_{p,eff}",$
       charsize=charactersize,xthick=linethickness,ythick=linethickness,charthick=linethickness
   ; yrange=0.35,1.5
   ;;;select = where( (totalprofilesaverage(*,0,0,5) gt minxplotsigma) and (totalprofiles(*,0,0,0,0,5) lt maxxplotsigma) and (totalprofilesaverage(*,6,0,5) gt 0))      
   select = where( (totalprofilesaverage(*,0,0,5) gt minxplotsigma+plotshift) and (totalprofilesaverage(*,6,0,5) gt 0))      
   last = max(select)

;   polyfill, [totalprofilesaverage(select,0,0,5),totalprofilesaverage(last,0,0,5),reverse(totalprofilesaverage(select,0,0,5))],$
;             [totalprofilesaverage(select,6,0,5)+totalprofilesaverage(select,6,2,5), $
;	      totalprofilesaverage(last,6,0,5)+totalprofilesaverage(last,6,2,5),$ 
;	     reverse(totalprofilesaverage(select,6,0,5)-totalprofilesaverage(select,6,2,5))],color=5,NOCLIP=0

   ;polyfill, [totalprofilesaverage(select,0,0,5),totalprofilesaverage(last,0,0,5),reverse(totalprofilesaverage(select,0,0,5))],$
   ;         [totalprofilesaverage(select,5,0,5)+totalprofilesaverage(select,5,2,5), $
;	      totalprofilesaverage(last,5,0,5)+totalprofilesaverage(last,5,2,5),$ 
;	     reverse(totalprofilesaverage(select,5,0,5)-totalprofilesaverage(select,5,2,5))],color=4,NOCLIP=0

;  polyfill, [totalprofilesaverage(select,0,0,5),totalprofilesaverage(last,0,0,5),reverse(totalprofilesaverage(select,0,0,5))],$
;             [totalprofilesaverage(select,1,0,5)+totalprofilesaverage(select,1,2,5), $
;      totalprofilesaverage(last,1,0,5)+totalprofilesaverage(last,1,2,5),$ 
;     reverse(totalprofilesaverage(select,1,0,5)-totalprofilesaverage(select,1,2,5))],color=10


  ; oplot, [totalprofilesaverage(select,0,0,5)],[totalprofilesaverage(select,6,0,5)-totalprofilesaverage(select,6,2,5)]$
;	     ,color=5,NOCLIP=0,thick=2.0*linethickness ;color=103;== cyan

 ;  oplot, [totalprofilesaverage(select,0,0,5)],[totalprofilesaverage(select,6,0,5)+totalprofilesaverage(select,6,2,5)]$
;	     ,color=5,NOCLIP=0,thick=2.0*linethickness ;color=103;== cyan

;   overplotromanowsky,totalprofilesaverage(*,1,0,5),totalprofilesaverage(*,0,0,5),interpolationpoint2Dsigma,thick=linethickness
   ; select only the non-zero regions

   fitromanowsky,totalprofilesaverage(*,1,0,5),totalprofilesaverage(*,0,0,5),totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,1,0,5),thick=linethickness
   ; select only the non-zero regions


;;   plotmyerror,4.875,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,5,0,5),totalprofilesaverage(select,5,2,5),color=1,thick=2.0*linethickness,linestyle=2
   ;plotmyerror,5.125,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,1,0,5),totalprofilesaverage(select,1,2,5),color=1,thick=2.0*linethickness
   ;plotmyerror,5.375,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,6,0,5),totalprofilesaverage(select,6,2,5),color=2,thick=2.0*linethickness
   ;plotmyerror,5.325,totalprofilesaverage(*,6,0,5),totalprofilesaverage(pixtouse,6,2,5),color=2,thick=1.0*linethickness
   ;plotmyerror,5.450,totalprofilesaverage(*,1,0,5),totalprofilesaverage(pixtouse,1,2,5),color=1,thick=1.0*linethickness

   plot,[0],[0],xrange=[minxplotsigma,maxxplotsigma],yrange=[0.4,1.6],ystyle=1,xstyle=1,charsize=charactersize,/noerase,xthick=linethickness,ythick=linethickness

   ;oplot,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,6,0,5),color=2,thick=2.0*linethickness
   ;oplot,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,5,0,5),color=1,linestyle=2,thick=2.0*linethickness
   ;oplot,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,4,0,5),color=1,linestyle=1
   ;oplot,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,3,0,5),color=205,thick=2.0*linethickness ; late
   ;oplot,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,2,0,5),color=4,thick=2.0*linethickness ; early
   ;oplot,totalprofilesaverage(select,0,0,5),totalprofilesaverage(select,1,0,5),color=1,thick=2.0*linethickness
   
   ;oplot,totalprofilesaverage(select,0,0),totalprofilesaverage(select,1,0)+totalprofilesaverage(select,1,2),color=1,thick=4.0*linethickness
   ;oplot,totalprofilesaverage(select,0,0),totalprofilesaverage(select,1,0)-totalprofilesaverage(select,1,2),color=1,thick=4.0*linethickness  

   ; totalprofiles: 0 x ; 1 sigma all ; 2 sigma early ; 3 sigma late ; 4 sigma old ; 5 sigma new   6 sigmadm
   ;
   ; ,,2: 0 sigma, 1 count
   ;
   ; 0: 3D rho, 1: 3D beta, 2: 3D sigma, 3: M, 4: 2D rho, 5: 2D sigma, 6: 2D kurtosis
   
   quantity = 5
   FOR i=5,5 DO BEGIN 
      FOR j=0L,simnumber-1 DO BEGIN
         FOR k=0,maxusesnapshotnumber-1 DO BEGIN
	    FOR dir=0,2 DO BEGIN
   	       ;oplot,totalprofiles(*,0,0,0,0,quantity,dir),totalprofiles(*,j,i,k,0,quantity,dir),color=1
	    ENDFOR 
      	 ENDFOR   
        ENDFOR      
   ENDFOR    
      
   ; Find the new star line that is the lowest at about 2.5 reff and plot it as well as its all stars and edge-on counterparts
   minimum =  1e10
   maximum = -1e10
   quantity = 5
   FOR i=5,5 DO BEGIN 
      FOR j=0L,simnumber-1 DO BEGIN
         FOR k=0,maxusesnapshotnumber-1 DO BEGIN
	    FOR dir=0,2 DO BEGIN
               xgt0        = where(totalprofiles(*,0,0,0,0,quantity,dir) gt 0)
               scalefactor = interpol(totalprofiles(xgt0,j,i,k,0,quantity,dir),totalprofiles(xgt0,0,0,0,0,quantity,dir),[2.5,2.5])
	       if (minimum gt scalefactor(0)) then begin
	          minimum = scalefactor(0)
		  kmin = k
		  imin = i
		  jmin = j
	       endif

	       if (maximum lt scalefactor(0)) then begin
	          maximum = scalefactor(0)
		  kmax = k
		  imax = i
		  jmax = j
	       endif
	       
   	    ENDFOR 
      	 ENDFOR   
        ENDFOR      
   ENDFOR    
   
   FOR dir=1,2 DO BEGIN   
      xgt0        = where(totalprofiles(*,0,0,0,0,quantity,dir) gt 0.3)
      ;;oplot,totalprofiles(xgt0,0,0,0,0,quantity,dir),totalprofiles(xgt0,jmax,imax,kmax,0,quantity,dir),color=0,linestyle=4
   ENDFOR        
   FOR dir=0,1 DO BEGIN   
      xgt0        = where(totalprofiles(*,0,0,0,0,quantity,dir) gt 0.3)
      ;;oplot,totalprofiles(xgt0,0,0,0,0,quantity,dir),totalprofiles(xgt0,jmin,imin,kmin,0,quantity,dir),color=0,linestyle=2
   ENDFOR        

   ; plot the fiducial case, sim=4, output 60 (k=0)    
   j = 0
   k = 0

   FOR dir=0,2,2 DO BEGIN   
      i   = 5   
      xgt0        = where(totalprofiles(*,0,0,0,0,quantity,dir) gt 0.3)
      oplot,totalprofiles(xgt0,0,0,0,0,quantity,dir),totalprofiles(xgt0,j,i,k,0,quantity,dir),color=1,linestyle=2,thick=2.0*linethickness
      i   = 6
      xgt0        = where(totalprofiles(*,0,0,0,0,quantity,dir) gt 0.3)
      oplot,totalprofiles(xgt0,0,0,0,0,quantity,dir),totalprofiles(xgt0,j,i,k,0,quantity,dir),color=2,thick=2.0*linethickness
      i   = 1
      xgt0        = where(totalprofiles(*,0,0,0,0,quantity,dir) gt 0.3)
      oplot,totalprofiles(xgt0,0,0,0,0,quantity,dir),totalprofiles(xgt0,j,i,k,0,quantity,dir),color=1,thick=2.0*linethickness      
   ENDFOR


   
   IF (ps ne 1) THEN BEGIN   
      image=tvrd(true=1);
      write_bmp,"a.bmp",image,/rgb
      spawnstring = 'convert a.bmp sigma2dstack.fiducial.all.gif'
      print,spawnstring
      spawn,spawnstring
      spawn,'rm a.bmp; mv *.gif Results/'
   ENDIF ELSE BEGIN
      device,/close
      set_plot,'x'
   ENDELSE  
   
   IF (ps eq 1) THEN BEGIN
      set_plot,'PS'
      device,/tt_font,filename='beta3dstack.eps',/color,XSIZE=20.0,YSIZE=20.0,/encapsulated
   ENDIF   
  
   ; STACKED PLOTS: FINAL 3DBeta
   ; -----------------------------------------------------------------------------------------------------
   plot,[0],[0],xrange=[minxplot,maxxplot],yrange=[-1,1.0],ystyle=1,xstyle=1,xtitle="r/R_{eff}",ytitle="\beta",$
        charsize=charactersize,/xlog,xthick=linethickness,ythick=linethickness,charthick=linethickness
   
   select = where( (totalprofilesaverage(*,0,0,1) gt minxplot) and (totalprofilesaverage(*,0,0,1) lt maxxplot))      
   last = max(select)
   polyfill, [totalprofilesaverage(select,0,0,1),totalprofilesaverage(last,0,0,1),reverse(totalprofilesaverage(select,0,0,1))],$
             [totalprofilesaverage(select,6,0,1)+totalprofilesaverage(select,6,2,1), $
	      totalprofilesaverage(last,6,0,1)+totalprofilesaverage(last,6,2,1),$ 
	     reverse(totalprofilesaverage(select,6,0,1)-totalprofilesaverage(select,6,2,1))],color=5,NOCLIP=0 ;color=103;== cyan

   polyfill, [totalprofilesaverage(select,0,0,1),totalprofilesaverage(last,0,0,1),reverse(totalprofilesaverage(select,0,0,1))],$
             [totalprofilesaverage(select,5,0,1)+totalprofilesaverage(select,5,2,1), $
	      totalprofilesaverage(last,5,0,1)+totalprofilesaverage(last,5,2,1),$ 
	     reverse(totalprofilesaverage(select,5,0,1)-totalprofilesaverage(select,5,2,1))],color=4,NOCLIP=0
	     
   polyfill, [totalprofilesaverage(select,0,0,1),totalprofilesaverage(last,0,0,1),reverse(totalprofilesaverage(select,0,0,1))],$
             [totalprofilesaverage(select,1,0,1)+totalprofilesaverage(select,1,2,1), $
	      totalprofilesaverage(last,1,0,1)+totalprofilesaverage(last,1,2,1),$ 
	     reverse(totalprofilesaverage(select,1,0,1)-totalprofilesaverage(select,1,2,1))],color=10,NOCLIP=0 ;color=103;== cyan

;   oplot, [totalprofilesaverage(select,0,0,1),totalprofilesaverage(last,0,0,1),reverse(totalprofilesaverage(select,0,0,1))],$
;             [totalprofilesaverage(select,6,0,1)+totalprofilesaverage(select,6,2,1), $
;	      totalprofilesaverage(last,6,0,1)+totalprofilesaverage(last,6,2,1),$ 
;	     reverse(totalprofilesaverage(select,6,0,1)-totalprofilesaverage(select,6,2,1))],color=5,NOCLIP=0,thick=2.0*linethickness ;color=103;== cyan
	     
   plot,[0],[0],xrange=[minxplot,maxxplot],yrange=[-1,1.0],ystyle=1,xstyle=1,charsize=charactersize,/xlog,/noerase,xthick=linethickness,ythick=linethickness	     
   
   oplot,totalprofilesaverage(select,0,0,1),totalprofilesaverage(select,6,0,1),color=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,1),totalprofilesaverage(select,5,0,1),color=1,linestyle=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,1),totalprofilesaverage(select,4,0,1),color=1,linestyle=1
   ;oplot,totalprofilesaverage(select,0,0,1),totalprofilesaverage(select,3,0,1),color=205,thick=2.0*linethickness ; late
   ;oplot,totalprofilesaverage(select,0,0,1),totalprofilesaverage(select,2,0,1),color=4,thick=2.0*linethickness ; early
   oplot,totalprofilesaverage(select,0,0,1),totalprofilesaverage(select,1,0,1),color=1,thick=2.0*linethickness
   oplot,[minxplot,maxxplot],[0,0],linestyle=1

   ;quantity = 1
   ; FOR i=5L,5 DO BEGIN 
   ;    FOR j=0L,simnumber-1 DO BEGIN
   ;       FOR k=0,maxusesnapshotnumber-1 DO BEGIN
   ; 	    oplot,totalprofiles(*,0,0,0,0,quantity),totalprofiles(*,j,i,k,0,quantity)
   ; 	 ENDFOR   
   ;  ENDFOR      
   ;ENDFOR    
   
   IF (ps ne 1) THEN BEGIN   
      image=tvrd(true=1);
      write_bmp,"a.bmp",image,/rgb
      spawnstring = 'convert a.bmp beta3dstack.all.gif'
      print,spawnstring
      spawn,spawnstring
      spawn,'rm a.bmp; mv *.gif Results/'
   ENDIF ELSE BEGIN
      device,/close
      set_plot,'x'
   ENDELSE  


   IF (ps eq 1) THEN BEGIN
      set_plot,'PS'
      device,/tt_font,filename='sigma3dstack.eps',/color,XSIZE=20.0,YSIZE=20.0,/encapsulated
   ENDIF   

   ; STACKED PLOTS: FINAL 3DSigma
   ; -----------------------------------------------------------------------------------------------------
   plot,[0],[0],xrange=[minxplot,maxxplot],yrange=[0.5,1.7],ystyle=1,xstyle=1,xtitle="r/R_{eff}",ytitle="\sigma/\sigma_{eff}",$ 
        charsize=charactersize,/ylog,/xlog,xthick=linethickness,ythick=linethickness,charthick=linethickness
   
   select = where( (totalprofilesaverage(*,0,0,2) gt minxplot) and (totalprofilesaverage(*,0,0,2) lt maxxplot))      
   last = max(select)
   polyfill, [totalprofilesaverage(select,0,0,2),totalprofilesaverage(last,0,0,2),reverse(totalprofilesaverage(select,0,0,2))],$
             [totalprofilesaverage(select,6,0,2)+totalprofilesaverage(select,6,2,2), $
	      totalprofilesaverage(last,6,0,2)+totalprofilesaverage(last,6,2,2),$ 
	     reverse(totalprofilesaverage(select,6,0,2)-totalprofilesaverage(select,6,2,2))],color=5,NOCLIP=0 ;color=103;== cyan

   polyfill, [totalprofilesaverage(select,0,0,2),totalprofilesaverage(last,0,0,2),reverse(totalprofilesaverage(select,0,0,2))],$
             [totalprofilesaverage(select,5,0,2)+totalprofilesaverage(select,5,2,2), $
      totalprofilesaverage(last,5,0,2)+totalprofilesaverage(last,5,2,2),$ 
	     reverse(totalprofilesaverage(select,5,0,2)-totalprofilesaverage(select,5,2,2))],color=4,NOCLIP=0

   polyfill, [totalprofilesaverage(select,0,0,2),totalprofilesaverage(last,0,0,2),reverse(totalprofilesaverage(select,0,0,2))],$
             [totalprofilesaverage(select,1,0,2)+totalprofilesaverage(select,1,2,2), $
	      totalprofilesaverage(last,1,0,2)+totalprofilesaverage(last,1,2,2),$ 
	     reverse(totalprofilesaverage(select,1,0,2)-totalprofilesaverage(select,1,2,2))],color=10,NOCLIP=0 ;color=103;== cyan

;   oplot, [totalprofilesaverage(select,0,0,2),totalprofilesaverage(last,0,0,2),reverse(totalprofilesaverage(select,0,0,2))],$
;             [totalprofilesaverage(select,6,0,2)+totalprofilesaverage(select,6,2,2), $
;	      totalprofilesaverage(last,6,0,2)+totalprofilesaverage(last,6,2,2),$ 
;	     reverse(totalprofilesaverage(select,6,0,2)-totalprofilesaverage(select,6,2,2))],color=5,NOCLIP=0,thick=2.0*linethickness ;color=103;== cyan

   plot,[0],[0],xrange=[minxplot,maxxplot],yrange=[0.5,1.7],ystyle=1,xstyle=1,charsize=charactersize,/ylog,/xlog,xthick=linethickness,$
        ythick=linethickness,charthick=linethickness,/noerase
	     
   oplot,totalprofilesaverage(select,0,0,2),totalprofilesaverage(select,6,0,2),color=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,2),totalprofilesaverage(select,5,0,2),color=1,linestyle=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,2),totalprofilesaverage(select,4,0,2),color=1,linestyle=1
   ;oplot,totalprofilesaverage(select,0,0,2),totalprofilesaverage(select,3,0,2),color=205,thick=2.0*linethickness ; late
   ;oplot,totalprofilesaverage(select,0,0,2),totalprofilesaverage(select,2,0,2),color=4,thick=2.0*linethickness ; early
   oplot,totalprofilesaverage(select,0,0,2),totalprofilesaverage(select,1,0,2),color=1,thick=2.0*linethickness

   ; Plot some slopes
   xmin = 0
   xmax = 0.6
   ymax = 0.18
   textshift = 0.025
   slope = 0.2
   oplot,[10.0^xmin,10.0^xmax],[10.0^ymax,10.0^(ymax-slope*(xmax-xmin))],thick=1.0*linethickness
   xyouts,10.0^(xmax+textshift),10.0^(ymax-slope*(xmax-xmin)-0.5*textshift),"-0.2",charsize=charactersize,charthick=linethickness

   ; quantity = 2
   ; FOR i=5L,5 DO BEGIN 
   ;    FOR j=0L,simnumber-1 DO BEGIN
   ;       FOR k=0,maxusesnapshotnumber-1 DO BEGIN
   ; 	    oplot,totalprofiles(*,0,0,0,0,quantity),totalprofiles(*,j,i,k,0,quantity)
   ;	 ENDFOR   
   ;  ENDFOR      
   ; ENDFOR    


   IF (ps ne 1) THEN BEGIN   
      image=tvrd(true=1);
      write_bmp,"a.bmp",image,/rgb
      spawnstring = 'convert a.bmp sigma3dstack.all.gif'
      print,spawnstring
      spawn,spawnstring
      spawn,'rm a.bmp; mv *.gif Results/'
   ENDIF ELSE BEGIN
      device,/close
      set_plot,'x'
   ENDELSE  


   IF (ps eq 1) THEN BEGIN
      set_plot,'PS'
      device,/tt_font,filename='m3dstack.eps',/color,XSIZE=20.0,YSIZE=20.0,/encapsulated
   ENDIF   

   ; STACKED PLOTS: FINAL M
   ; -----------------------------------------------------------------------------------------------------
   plot,[0],[0],xrange=[minxplot,maxxplot],yrange=[0.1,9],ystyle=1,xstyle=1,xtitle="r/R_{eff}",ytitle="M/M_{eff}",$ 
       charsize=charactersize,/ylog,/xlog,xthick=linethickness,ythick=linethickness,charthick=linethickness
   
   select = where( (totalprofilesaverage(*,0,0,3) gt minxplot) and (totalprofilesaverage(*,0,0,3) lt maxxplot))      
   last = max(select)
   polyfill, [totalprofilesaverage(select,0,0,3),totalprofilesaverage(last,0,0,3),reverse(totalprofilesaverage(select,0,0,3))],$
             [totalprofilesaverage(select,6,0,3)+totalprofilesaverage(select,6,2,3), $
	      totalprofilesaverage(last,6,0,3)+totalprofilesaverage(last,6,2,3),$ 
	     reverse(totalprofilesaverage(select,6,0,3)-totalprofilesaverage(select,6,2,3))],color=5,NOCLIP=0 ;color=103;== cyan

   polyfill, [totalprofilesaverage(select,0,0,3),totalprofilesaverage(last,0,0,3),reverse(totalprofilesaverage(select,0,0,3))],$
             [totalprofilesaverage(select,5,0,3)+totalprofilesaverage(select,5,2,3), $
	      totalprofilesaverage(last,5,0,3)+totalprofilesaverage(last,5,2,3),$ 
	     reverse(totalprofilesaverage(select,5,0,3)-totalprofilesaverage(select,5,2,3))],color=4,NOCLIP=0
	     
   polyfill, [totalprofilesaverage(select,0,0,3),totalprofilesaverage(last,0,0,3),reverse(totalprofilesaverage(select,0,0,3))],$
             [totalprofilesaverage(select,1,0,3)+totalprofilesaverage(select,1,2,3), $
	      totalprofilesaverage(last,1,0,3)+totalprofilesaverage(last,1,2,3),$ 
	     reverse(totalprofilesaverage(select,1,0,3)-totalprofilesaverage(select,1,2,3))],color=10,NOCLIP=0 ;color=103;== cyan

;   oplot, [totalprofilesaverage(select,0,0,3),totalprofilesaverage(last,0,0,3),reverse(totalprofilesaverage(select,0,0,3))],$
;             [totalprofilesaverage(select,6,0,3)+totalprofilesaverage(select,6,2,3), $
;	      totalprofilesaverage(last,6,0,3)+totalprofilesaverage(last,6,2,3),$ 
;	     reverse(totalprofilesaverage(select,6,0,3)-totalprofilesaverage(select,6,2,3))],color=5,NOCLIP=0,thick=2.0*linethickness ;color=103;== cyan
	     
   plot,[0],[0],xrange=[minxplot,maxxplot],yrange=[0.1,9],ystyle=1,xstyle=1,charsize=charactersize,/ylog,/xlog,/noerase,xthick=linethickness,ythick=linethickness

   ; add the total V(r)/V_{eff} line
   vofr  = sqrt((totalprofilesaverage(select,6,0,3)+totalprofilesaverage(select,1,0,3))/totalprofilesaverage(select,0,0,3))
   vofre = interpol(vofr,totalprofilesaverage(select,0,0,3),[interpolationpointrest,interpolationpointrest])
   oplot,totalprofilesaverage(select,0,0,3),vofr/vofre(0),color=3,thick=2.0*linethickness
   xyouts,0.55,1.225,'V/V_{eff}',charsize=charactersize,charthick=linethickness,color=3
   
   mtotofr  = (totalprofilesaverage(select,6,0,3)+totalprofilesaverage(select,1,0,3))
   mtotofre = interpol(mtotofr,totalprofilesaverage(select,0,0,3),[interpolationpointrest,interpolationpointrest])
   ;oplot,totalprofilesaverage(select,0,0,3),mtotofr/mtotofre(0),color=9,thick=2.0*linethickness   
	     
   oplot,totalprofilesaverage(select,0,0,3),totalprofilesaverage(select,6,0,3),color=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,3),totalprofilesaverage(select,5,0,3),color=1,linestyle=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,3),totalprofilesaverage(select,4,0,3),color=1,linestyle=1
   ;oplot,totalprofilesaverage(select,0,0,3),totalprofilesaverage(select,3,0,3),color=205,thick=2.0*linethickness ; late
   ;oplot,totalprofilesaverage(select,0,0,3),totalprofilesaverage(select,2,0,3),color=4,thick=2.0*linethickness ; early
   oplot,totalprofilesaverage(select,0,0,3),totalprofilesaverage(select,1,0,3),color=1,thick=2.0*linethickness

   ; "mass to light ratio"
   mtol = mtotofr/totalprofilesaverage(select,1,0,3)
   radius = [1,3,5,8,10]
   mtolre = interpol(mtol,totalprofilesaverage(select,0,0,3),radius)
   print,'M/L     ',mtolre
   print,'at reff:',radius
   
   ; quantity = 3
   ; FOR i=5L,5 DO BEGIN 
   ;    FOR j=0L,simnumber-1 DO BEGIN
   ;       FOR k=0,maxusesnapshotnumber-1 DO BEGIN
   ; 	    oplot,totalprofiles(*,0,0,0,0,quantity),totalprofiles(*,j,i,k,0,quantity)
   ; 	 ENDFOR   
   ;  ENDFOR      
   ;ENDFOR    


   IF (ps ne 1) THEN BEGIN      
      image=tvrd(true=1);
      write_bmp,"a.bmp",image,/rgb
      spawnstring = 'convert a.bmp m3dstack.all.gif'
      print,spawnstring
      spawn,spawnstring
      spawn,'rm a.bmp; mv *.gif Results/'
   ENDIF ELSE BEGIN
      device,/close
      set_plot,'x'
   ENDELSE  


   IF (ps eq 1) THEN BEGIN
      set_plot,'PS'
      device,/tt_font,filename='rho3dstack.eps',/color,XSIZE=20.0,YSIZE=20.0,/encapsulated
   ENDIF   

   ; STACKED PLOTS: FINAL rho 3D
   ; -----------------------------------------------------------------------------------------------------
   plot,[0],[0],xrange=[minxplot,maxxplot],yrange=[0.0001,20],ystyle=1,xstyle=1,xtitle="r/R_{eff}",ytitle="\rho/\rho_{eff}",charsize=charactersize,/ylog,/xlog, $
        YTickFormat='exponent',xthick=linethickness,ythick=linethickness,charthick=linethickness
  
   select = where( (totalprofilesaverage(*,0,0,0) gt minxplot) and (totalprofilesaverage(*,0,0,0) lt maxxplot))      
   last = max(select)
   polyfill, [totalprofilesaverage(select,0,0,0),totalprofilesaverage(last,0,0,0),reverse(totalprofilesaverage(select,0,0,0))],$
             [totalprofilesaverage(select,6,0,0)+totalprofilesaverage(select,6,2,0), $
	      totalprofilesaverage(last,6,0,0)+totalprofilesaverage(last,6,2,0),$ 
	     reverse(totalprofilesaverage(select,6,0,0)-totalprofilesaverage(select,6,2,0))],color=5,NOCLIP=0 ;color=103;== cyan

   polyfill, [totalprofilesaverage(select,0,0,0),totalprofilesaverage(last,0,0,0),reverse(totalprofilesaverage(select,0,0,0))],$
             [totalprofilesaverage(select,5,0,0)+totalprofilesaverage(select,5,2,0), $
	      totalprofilesaverage(last,5,0,0)+totalprofilesaverage(last,5,2,0),$ 
	     reverse(totalprofilesaverage(select,5,0,0)-totalprofilesaverage(select,5,2,0))],color=4,NOCLIP=0

   polyfill, [totalprofilesaverage(select,0,0,0),totalprofilesaverage(last,0,0,0),reverse(totalprofilesaverage(select,0,0,0))],$
             [totalprofilesaverage(select,1,0,0)+totalprofilesaverage(select,1,2,0), $
	      totalprofilesaverage(last,1,0,0)+totalprofilesaverage(last,1,2,0),$ 
	     reverse(totalprofilesaverage(select,1,0,0)-totalprofilesaverage(select,1,2,0))],color=10 ;color=103;== cyan

;   oplot, [totalprofilesaverage(select,0,0,0),totalprofilesaverage(last,0,0,0),reverse(totalprofilesaverage(select,0,0,0))],$
;             [totalprofilesaverage(select,6,0,0)+totalprofilesaverage(select,6,2,0), $
;	      totalprofilesaverage(last,6,0,0)+totalprofilesaverage(last,6,2,0),$ 
;	     reverse(totalprofilesaverage(select,6,0,0)-totalprofilesaverage(select,6,2,0))],color=5,NOCLIP=0,thick=2.0*linethickness ;color=103;== cyan

   plot,[0],[0],xrange=[minxplot,maxxplot],yrange=[0.0001,20],ystyle=1,xstyle=1,charsize=charactersize,/ylog,/xlog, $
        YTickFormat='exponent',xthick=linethickness,ythick=linethickness,charthick=linethickness,/noerase
	     
   oplot,totalprofilesaverage(select,0,0,0),totalprofilesaverage(select,6,0,0),color=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,0),totalprofilesaverage(select,5,0,0),color=1,linestyle=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,0),totalprofilesaverage(select,4,0,0),color=1,linestyle=1
   ;oplot,totalprofilesaverage(select,0,0,0),totalprofilesaverage(select,3,0,0),color=205,thick=2.0*linethickness ; late
   ;oplot,totalprofilesaverage(select,0,0,0),totalprofilesaverage(select,2,0,0),color=4,thick=2.0*linethickness ; early
   oplot,totalprofilesaverage(select,0,0,0),totalprofilesaverage(select,1,0,0),color=1,thick=2.0*linethickness

   ; Plot some slopes
   xmin = 0.1
   xmax = 0.7
   ymax = 0.95
   textshift = 0.025
   oplot,[10.0^xmin,10.0^xmax],[10.0^ymax,10.0^(ymax-2*(xmax-xmin))],thick=1.0*linethickness
   oplot,[10.0^xmin,10.0^xmax],[10.0^ymax,10.0^(ymax-3*(xmax-xmin))],thick=1.0*linethickness   
   xyouts,10.0^(xmax+textshift),10.0^(ymax-2*(xmax-xmin)-4*textshift),'-2',charsize=charactersize,charthick=linethickness
   xyouts,10.0^(xmax+textshift),10.0^(ymax-3*(xmax-xmin)-4*textshift),'-3',charsize=charactersize,charthick=linethickness

   ; quantity = 0
   ; FOR i=5L,5 DO BEGIN 
   ;    FOR j=0L,simnumber-1 DO BEGIN
   ;       FOR k=0,maxusesnapshotnumber-1 DO BEGIN
   ; 	    oplot,totalprofiles(*,0,0,0,0,quantity),totalprofiles(*,j,i,k,0,quantity)
   ; 	 ENDFOR   
   ;  ENDFOR      
   ;ENDFOR    


   IF (ps ne 1) THEN BEGIN      
      image=tvrd(true=1);
      write_bmp,"a.bmp",image,/rgb
      spawnstring = 'convert a.bmp rho3dstack.all.gif'
   ;   spawnstring = 'convert a.bmp rho3d.new.gif'   
      print,spawnstring
      spawn,spawnstring
      spawn,'rm a.bmp; mv *.gif Results/'
   ENDIF ELSE BEGIN
      device,/close
      set_plot,'x'
   ENDELSE  


   IF (ps eq 1) THEN BEGIN
      set_plot,'PS'
      device,/tt_font,filename='rho2dstack.eps',/color,XSIZE=20.0,YSIZE=15.0,/encapsulated
   ENDIF   

   ; STACKED PLOTS: FINAL rho 2D
   ; -----------------------------------------------------------------------------------------------------
   plot,[0],[0],xrange=[minxplot,maxxplot],yrange=[0.0018,10],ystyle=1,xstyle=1,xtitle="r_p/R_{eff}",ytitle="\Sigma/\Sigma_{eff}",charsize=charactersize,/ylog,/xlog, $
        YTickFormat='exponent',xthick=linethickness,ythick=linethickness,charthick=linethickness
   
   select = where( (totalprofilesaverage(*,0,0,4) gt minxplot) and (totalprofilesaverage(*,0,0,4) lt maxxplot))      
   last = max(select)
   polyfill, [totalprofilesaverage(select,0,0,4),totalprofilesaverage(last,0,0,4),reverse(totalprofilesaverage(select,0,0,4))],$
             [totalprofilesaverage(select,6,0,4)+totalprofilesaverage(select,6,2,4), $
	      totalprofilesaverage(last,6,0,4)+totalprofilesaverage(last,6,2,4),$ 
	     reverse(totalprofilesaverage(select,6,0,4)-totalprofilesaverage(select,6,2,4))],color=5,NOCLIP=0 ;color=103;== cyan

   polyfill, [totalprofilesaverage(select,0,0,4),totalprofilesaverage(last,0,0,4),reverse(totalprofilesaverage(select,0,0,4))],$
             [totalprofilesaverage(select,5,0,4)+totalprofilesaverage(select,5,2,4), $
	      totalprofilesaverage(last,5,0,4)+totalprofilesaverage(last,5,2,4),$ 
	     reverse(totalprofilesaverage(select,5,0,4)-totalprofilesaverage(select,5,2,4))],color=4,NOCLIP=0

   polyfill, [totalprofilesaverage(select,0,0,4),totalprofilesaverage(last,0,0,4),reverse(totalprofilesaverage(select,0,0,4))],$
             [totalprofilesaverage(select,1,0,4)+totalprofilesaverage(select,1,2,4), $
	      totalprofilesaverage(last,1,0,4)+totalprofilesaverage(last,1,2,4),$ 
	     reverse(totalprofilesaverage(select,1,0,4)-totalprofilesaverage(select,1,2,4))],color=10 ;color=103;== cyan

;   oplot, [totalprofilesaverage(select,0,0,4),totalprofilesaverage(last,0,0,4),reverse(totalprofilesaverage(select,0,0,4))],$
;             [totalprofilesaverage(select,6,0,4)+totalprofilesaverage(select,6,2,4), $
;	      totalprofilesaverage(last,6,0,4)+totalprofilesaverage(last,6,2,4),$ 
;	     reverse(totalprofilesaverage(select,6,0,4)-totalprofilesaverage(select,6,2,4))],color=5,NOCLIP=0,thick=2.0*linethickness ;color=103;== cyan

    plot,[0],[0],xrange=[minxplot,maxxplot],yrange=[0.0018,10],ystyle=1,xstyle=1,charsize=charactersize,/ylog,/xlog, $
        YTickFormat='exponent',xthick=linethickness,ythick=linethickness,charthick=linethickness,/noerase
   

   n3379pe =  read_ascii("n3379_peletier.new.dat",comment_symbol="#")
   n4697pe =  read_ascii("n4697_peletier.new.dat",comment_symbol="#")
   n3379vc =  read_ascii("n3379_dVC.dat",comment_symbol="#")

   reffdata = fltarr(3)
   reffdata(0) = 55.9 ; arcsec
   reffdata(1) = 73.8 ; arcsec
   reffdata(2) = 55.9 ; arcsec   
   
   sigdata       = fltarr(3,2,200)

   sigdata(0,0,0:n_elements(n3379pe.field1(0,*))-1) = n3379pe.field1(0,*) / reffdata(0)
   sigdata(1,0,0:n_elements(n4697pe.field1(0,*))-1) = n4697pe.field1(0,*) / reffdata(1)
   sigdata(2,0,0:n_elements(n3379vc.field1(0,*))-1) = n3379vc.field1(0,*) / reffdata(2)

   sigdata(0,1,0:n_elements(n3379pe.field1(1,*))-1) = n3379pe.field1(1,*)
   sigdata(1,1,0:n_elements(n4697pe.field1(1,*))-1) = n4697pe.field1(1,*)
   sigdata(2,1,0:n_elements(n3379vc.field1(1,*))-1) = n3379vc.field1(1,*)
   
   ; normalize to the good mu_eff
   
   FOR i=0,2 DO BEGIN
      selectthis              = where(sigdata(i,0,*) gt 0)
      scaling                 = interpol(sigdata(i,1,selectthis),sigdata(i,0,selectthis),[interpolationpointrest,interpolationpointrest])
      print,'scaling sigma',10.0^(0.4*scaling(0))
      sigdata(i,1,selectthis) = 10.0^(-0.4*(sigdata(i,1,selectthis) - scaling(0)))      
   ENDFOR

   ;oplot,totalprofilesaverage(select,0,0,4),totalprofilesaverage(select,1,0,4)/10.0^(1.5),color=1   
   FOR i=2,1,-1 DO BEGIN
      colorstouse=[6,7,8]
      selectthis = where( (sigdata(i,0,*) gt 0) and (sigdata(i,0,*) gt minxplot*1.1) and (sigdata(i,0,*) lt maxxplot*0.9) )
      ;plotsym,3,/fill;      
      usersym,1.5*[-1,-1,1,1,-1],1.5*[-1,1,1,-1,-1],thick=4.0
      oplot,sigdata(i,0,selectthis),sigdata(i,1,selectthis),color=colorstouse(i),thick=2.0*linethickness,psym=7
      ;oplot,sigdata(i,0,selectthis),sigdata(i,1,selectthis)/10.0^(1.5),color=colorstouse(i),thick=2.0*linethickness      
   ENDFOR
   
   ; overplot arrow
   ;arrow,1,0.85,1,0.037,/data,color=0
   	     
   oplot,totalprofilesaverage(select,0,0,4),totalprofilesaverage(select,6,0,4),color=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,4),totalprofilesaverage(select,5,0,4),color=1,linestyle=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,4),totalprofilesaverage(select,4,0,4),color=1,linestyle=1
   ;oplot,totalprofilesaverage(select,0,0,4),totalprofilesaverage(select,3,0,4),color=205,thick=2.0*linethickness ; late
   ;oplot,totalprofilesaverage(select,0,0,4),totalprofilesaverage(select,2,0,4),color=4,thick=2.0*linethickness ; early
   oplot,totalprofilesaverage(select,0,0,4),totalprofilesaverage(select,1,0,4),color=1,thick=2.0*linethickness

   ; Plot some slopes
   xmin = -0.35;0
   xmax = 0.25;0.6
   ymax = -1.2;1.2
   textshift = 0.025
   oplot,[10.0^xmin,10.0^xmax],[10.0^ymax,10.0^(ymax-1*(xmax-xmin))],thick=1.0*linethickness
   oplot,[10.0^xmin,10.0^xmax],[10.0^ymax,10.0^(ymax-2*(xmax-xmin))],thick=1.0*linethickness
   xyouts,10.0^(xmax+textshift),10.0^(ymax-1*(xmax-xmin)-4*textshift),'-1',charsize=charactersize,charthick=linethickness
   xyouts,10.0^(xmax+textshift),10.0^(ymax-2*(xmax-xmin)-4*textshift),'-2',charsize=charactersize,charthick=linethickness

   ; totalprofiles: 0 x ; 1 sigma all ; 2 sigma early ; 3 sigma late ; 4 sigma old ; 5 sigma new   6 sigmadm
   ;
   ; ,,2: 0 sigma, 1 count
   ;
   ; 0: 3D rho, 1: 3D beta, 2: 3D sigma, 3: M, 4: 2D rho, 5: 2D sigma, 6: 2D kurtosis

   ;quantity = 4
   ;FOR i=5L,5 DO BEGIN 
   ;   FOR j=0L,simnumber-1 DO BEGIN
   ;      FOR k=0,maxusesnapshotnumber-1 DO BEGIN
;	    FOR dir=0,2 DO BEGIN
;                oplot,totalprofiles(*,0,0,0,0,quantity,dir),totalprofiles(*,j,i,k,0,quantity,dir)
;	    ENDFOR
;   	 ENDFOR   
;     ENDFOR      
;   ENDFOR    

      
   openw,1,'Sigma2d.dat'
   FOR i=0L,n_elements(select)-1 DO BEGIN
      printf,1,totalprofilesaverage(select(i),0,0,4),totalprofilesaverage(select(i),1,0,4),totalprofilesaverage(select(i),6,0,4),totalprofilesaverage(select(i),5,0,4)
   ENDFOR
   close,1
   

   IF (ps ne 1) THEN BEGIN      
      image=tvrd(true=1);
      write_bmp,"a.bmp",image,/rgb
      spawnstring = 'convert a.bmp rho2dstack.all.gif'
      print,spawnstring
      spawn,spawnstring
      spawn,'rm a.bmp; mv *.gif Results/'
   ENDIF ELSE BEGIN
      device,/close
      set_plot,'x'
   ENDELSE  

   IF (ps eq 1) THEN BEGIN
      set_plot,'PS'
      device,/tt_font,filename='kurtosis2dstack.eps',/color,XSIZE=20.0,YSIZE=15.0,/encapsulated
   ENDIF   

   ; STACKED PLOTS: FINAL kurtosis 2D
   ; -----------------------------------------------------------------------------------------------------
   plot,[0],[0],xrange=[minxplot,maxxplot],yrange=[-1,1],ystyle=1,xstyle=1,xtitle="r/R_{eff}",ytitle="\kappa",charsize=charactersize,/xlog, $
        xthick=linethickness,ythick=linethickness,charthick=linethickness
  
   select = where( (totalprofilesaverage(*,0,0,6) gt minxplot) and (totalprofilesaverage(*,0,0,6) lt maxxplot))      
   last = max(select)
   polyfill, [totalprofilesaverage(select,0,0,6),totalprofilesaverage(last,0,0,6),reverse(totalprofilesaverage(select,0,0,6))],$
             [totalprofilesaverage(select,6,0,6)+totalprofilesaverage(select,6,2,6), $
	      totalprofilesaverage(last,6,0,6)+totalprofilesaverage(last,6,2,6),$ 
	     reverse(totalprofilesaverage(select,6,0,6)-totalprofilesaverage(select,6,2,6))],color=5,NOCLIP=0 ;color=103;== cyan

   polyfill, [totalprofilesaverage(select,0,0,6),totalprofilesaverage(last,0,0,6),reverse(totalprofilesaverage(select,0,0,6))],$
             [totalprofilesaverage(select,5,0,6)+totalprofilesaverage(select,5,2,6), $
	      totalprofilesaverage(last,5,0,6)+totalprofilesaverage(last,5,2,6),$ 
	     reverse(totalprofilesaverage(select,5,0,6)-totalprofilesaverage(select,5,2,6))],color=4,NOCLIP=0

   polyfill, [totalprofilesaverage(select,0,0,6),totalprofilesaverage(last,0,0,6),reverse(totalprofilesaverage(select,0,0,6))],$
             [totalprofilesaverage(select,1,0,6)+totalprofilesaverage(select,1,2,6), $
	      totalprofilesaverage(last,1,0,6)+totalprofilesaverage(last,1,2,6),$ 
	     reverse(totalprofilesaverage(select,1,0,6)-totalprofilesaverage(select,1,2,6))],color=10,NOCLIP=0 ;color=103;== cyan

;   oplot, [totalprofilesaverage(select,0,0,6),totalprofilesaverage(last,0,0,6),reverse(totalprofilesaverage(select,0,0,6))],$
;             [totalprofilesaverage(select,6,0,6)+totalprofilesaverage(select,6,2,6), $
;	      totalprofilesaverage(last,6,0,6)+totalprofilesaverage(last,6,2,6),$ 
;	     reverse(totalprofilesaverage(select,6,0,6)-totalprofilesaverage(select,6,2,6))],color=5,NOCLIP=0,thick=2.0*linethickness ;color=103;== cyan

   plot,[0],[0],xrange=[minxplot,maxxplot],yrange=[-1,1],ystyle=1,xstyle=1,charsize=charactersize,/xlog, $
        xthick=linethickness,ythick=linethickness,charthick=linethickness,/noerase
   oplot,[0.00001,1000],[0,0],linestyle=1

   ;overplot the single kurtosises
   FOR j=0L,simnumber-1 DO BEGIN
      FOR k=0L,maxsnapshotnumber-1 DO BEGIN
       FOR dir=0,2 DO BEGIN
       ;oplot,totalprofiles(*,0,0,0,0,6,dir),totalprofiles(*,j,5,k,0,6,dir)
       ENDFOR
      ENDFOR
   ENDFOR
	     
   oplot,totalprofilesaverage(select,0,0,6),totalprofilesaverage(select,6,0,6),color=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,6),totalprofilesaverage(select,5,0,6),color=1,linestyle=2,thick=2.0*linethickness
   oplot,totalprofilesaverage(select,0,0,6),totalprofilesaverage(select,4,0,6),color=1,linestyle=1
   ;oplot,totalprofilesaverage(select,0,0,6),totalprofilesaverage(select,3,0,6),color=205,thick=2.0*linethickness ; late
   ;oplot,totalprofilesaverage(select,0,0,6),totalprofilesaverage(select,2,0,6),color=4,thick=2.0*linethickness ; early
   oplot,totalprofilesaverage(select,0,0,6),totalprofilesaverage(select,1,0,6),color=1,thick=2.0*linethickness

   IF (ps ne 1) THEN BEGIN      
      image=tvrd(true=1);
      write_bmp,"a.bmp",image,/rgb
      spawnstring = 'convert a.bmp kurtosis3dstack.all.gif'
      print,spawnstring
      spawn,spawnstring
      spawn,'rm a.bmp; mv *.gif Results/'
   ENDIF ELSE BEGIN
      device,/close
      set_plot,'x'
   ENDELSE  


   IF (ps eq 1) THEN BEGIN
      spawn,'cp *.eps Results/'   
      spawn,'latex latexfigures1 >&/dev/null'
      spawn,'latex latexfigures2 >&/dev/null'
      spawn,'latex latexfigures3 >&/dev/null'      
      spawn,'latex latexfigures4 >&/dev/null'            
      spawn,'latex latexfigures >&/dev/null'
      spawn,'dvips latexfigures >&/dev/null'      
      spawn,'dvips latexfigures1 -o dekel_etal_fig1.ps>&/dev/null'      
      spawn,'dvips latexfigures2 -o dekel_etal_fig2.ps>&/dev/null'                  
      spawn,'dvips latexfigures3 -o dekel_etal_fig1_small.ps >&/dev/null'      
      spawn,'dvips latexfigures4 -o dekel_etal_fig2_small.ps>&/dev/null'                              
      spawn,'cp latexfigures.ps Results/'
      spawn,'tar -cvf figures.tar *.eps >&/dev/null; gzip -f figures.tar;  cp figures.tar.gz Results/  '
      spawn,'rm *.gif'
      spawn,'makegifs'
      spawn,'tar -cvf dekel_etal.tar dekel_etal*.ps kurtosis2dstack.eps sigma2dstack.fiducial.eps *.gif>&/dev/null; gzip -f dekel_etal.tar;  cp dekel_etal.tar.gz Results/'
   ENDIF

end	








