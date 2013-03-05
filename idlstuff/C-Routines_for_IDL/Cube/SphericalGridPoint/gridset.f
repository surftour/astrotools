      subroutine gridset

c  set up grid:  3-D spherical grid with variable spacing in r and 
c  theta (set by exponents, rexp, texp)
c  make a bunch of arrays necessary for find_wall subroutine

c history:
c 00/09/06 (baw): write subroutine

      implicit none
      include 'tts.txt'
      include 'out.txt'
      include 'grid.txt'
      include 'opacin.txt'
      include 'dust.com'

      real*8 dr,dt,dp,rad,phi,thet,pi,r2p,densd,dense,vol
     $ ,cost,sint,pihalf,eps,tau,rming1,rming2,rming3,tiny
     $     ,rmincgs,taud,taue,mu0,tauv,taufmin,mu,thetatm
     $     ,thetmin,tauzrat,sigmu,rsonr,Tdisk,taudisk
     $     ,mudisk,c1,a1,tauw,rtop,rfac,tauRr,taumid,tauRth
     $     ,tauRos,res,thettop,maxdens,tauRmu,radamb,rbeg2
     $     ,tauave,gn,sigt,lstar,mdot,lacc,lacc2,ltot,tave,Vc
     $     ,rmine2,k,muc,mH,rcgs,mcgs,honr,h,fudge,fluxratio
     $     ,fplay,dphi

      real*4 x_central(600),y_central(600),z_central(600)
     $     ,clumpRad(600),clumpsz(600),clumpDensity(600)
     $     ,theta_clump(600),phi_clump(600),suk(600)
     $     ,mu_clump(600)

      real xyarr(200,200),x,y,xmax,dx,r,dens1,r_sort,suk_min
     $     ,suk_exp

      real chiR,ierfc,erfc
      real ran2

      integer ir,it,ip,nttmp,gflag,n,ix,iy,nr0,ntatm,nrwall,irbeg,
     $     dcount,id,ide,idd,i,idum,i1,j,N_clump
      external erfc

      pi=4.d0*datan(1.d0)
      r2p=2.d0*pi
      pihalf=pi/2.d0
      tiny=1.d-15
      rmincgs=rstar*rsol
      radamb=rmax
      gn=6.6732d-8
      sigt=5.67d-5

      if (idiskacc.eq.1) then
c     set disk accretion fraction based on alpha disk paramters
c     see eqns 5 and 6 in Whitney et al. 2003, apj, 591, 1049
         lstar=4.d0*pi*sigt*tstar**4*rmincgs**2
         Vc=dsqrt(gn*msol*massc/rmincgs)
c         print*,'alphad, pi, Vc, rho0, z1, rmincgs'
         mdot=dsqrt(18*pi**3)*alphad*Vc*rho0*(z1*rmincgs)**3
     1        /rmincgs
c     accretion in disk
         lacc=gn*(msol*massc)*mdot/2.d0/(rddust*rmincgs)
c     accretion in shock, the optically thick part from the heated
c     atmosphere  (calvet & gullbring 1998).
c     we will split it into half blackbody at Tshock, and half x-rays.
c     assumes co-Rotation radius is 5 Rstar.
         lacc2=gn*(msol*massc)*mdot*(1./rmincgs-1./(rmincgs*rtrunc))
         ltot=lacc+lacc2+lstar
         accfrac=lacc/ltot
         accfrac2=lacc2/ltot
         mdot=mdot*3600.*24.*365.25/msol
         write(12,*) 'including disk accretion luminosity'
         write(12,*) 'disk accretion rate (msol/yr) ',mdot
         write(12,*) 'fraction of luminosity from disk ',accfrac
         write(12,*) 'fraction of accretion lum. on star ',accfrac2
         write(12,*) 'total system luminosity (Lsun)',ltot/3.845d33
         print*, 'including disk accretion luminosity '
         print*, 'disk accretion rate (msol/yr) ',mdot
         print*, 'fraction of acc luminosity from disk ',accfrac
         print*, 'fraction of acc luminosity on star ',accfrac2
         print*, 'total system luminosity (Lsun)',ltot/3.845d33
c     f=calvet's filling factor, ranges from 0.01 to 0.001
c     median is 0.007
c     jon's formal solution to the grey problem
         fluxratio=0.5*lacc2/lstar/fspot
         tshock=tstar*(1.+fluxratio)**0.25
         print*,'calculated thermal shock T ',tshock
         write(12,*) 'calculated thermal shock T ',tshock
         fplay=0.5*lacc2/fspot/4.d0/pi/rmincgs**2
         print*,'log flux shock ',dlog10(fplay)
      else
         lstar=4.d0*pi*sigt*tstar**4*rmincgs**2
         mdot=0.d0
         accfrac=0.d0
         lacc=0.d0
         lacc2=0.d0
         write(12,*) 'not including disk accretion luminosity'
         print*, 'disk accretion rate (msol/yr) ',mdot
         print*, 'fraction of luminosity from disk ',accfrac         
         write(12,*) 'total system luminosity (Lsun)',lstar/3.845d33
         print*, 'total system luminosity (Lsun)',lstar/3.845d33
      endif

c     calculate an average stellar temp based on luminosity
c     of star+accretion
      ltot=lacc2+lstar
      tave=(ltot/lstar)**0.25*tstar
      tave=((lacc2+lstar-lacc)/(lstar))**0.25*tstar
      print*,'average stellar temp is now',tave
      rmine2=(tsub/tave)**(-2.1d0)
      write(12,*) 'average Tstar inc. accretion ',Tave
c     reset tstar to tave, since it is used for temperature
c     calculation (sets the luminosity scale).
      tstar=tave

      if (rmine2.gt.rmaxd.and.massd.gt.0.0) then
         print*,''
         print*,'****************'
         print*,'ERROR.  disk inner edge larger than RMAXD'
         print*,'probably because you chose a hot star and'
         print*,'dust sublimation radius is large.'
         print*,'Increase RMAXD.'
         print*,'RMIND, RMAXD in AU',rmine2/autors,rmaxd/autors
         print*,'Stopping program'
         stop
      endif

      if (irminsub.eq.1) then
         rmine=rmine2*rmine_in
         rmind=rmine2*rmind_in
         rddust=rmine2*rmind_in
         print*,'resetting Rsub to ',rmine2
         print*,'inner envelope radius ',rmine
         print*,'inner disk radius ',rmind
         write(12,*) 'updated dust destruction radius ',rmine2
         write(12,*) 'updated inner envelope radius',rmine
         write(12,*) 'updated inner disk radius',rmind
      endif

c     convert opacity to units of cm^2/gm*rstar(cgs) because distance
c     is in units if 1/rstar, and dtau=kapd*rho*ds
	do id=1,4
         kapd(id)=kappav(id)*rmincgs
      enddo
c     print*,'kappav',kappav

c     calculate some constants for the envelope density
c	TODO, modify envset for new envelope density
      call envset

c     make grid include minimum and maximum values.

      print*, ' '
      print*, 'grid setup'
      pihalf=pi/2.d0

      print*,'rddust,rmine,rmin',rddust,rmine,rmin
c     rgrid
      gflag=0
      rarr(1)=rmin
      if (rddust.lt.rmine) then
         gflag=1
         rarr(2)=rddust
         rming1=rddust
         rming2=rmine
         irbeg=2
      else if (rddust.gt.rmine) then
         gflag=2
         rarr(2)=rmine
         rming1=rmine
         rming2=rddust
         irbeg=3
      else if (rddust.eq.rmine) then
         gflag=3
         rarr(2)=rmine
         rming1=rmine
         irbeg=2
      endif
      if (massd.eq.0.) then
c     since there's empty space below rmine and rddust....
c S.C.         dr=(rmax-rming1)/(dfloat(nrg-2))**rexp
c S.C.         print*,'rexp,dr,rmin',rexp,dr/autors,rmin/autors

         do ir=3,nrg

            rarr(ir)=rming1*(exp((log(rmax/rming1))*
     $           ((dfloat(ir-2))/(nrg-2))))

c     print*,'rarr ',ir,rarr(ir)/au

         enddo
         
      else

c     lots of fine spacing in inner region of disk
         
         a1=1.d0-a
         
         mu0=z1*Rddust**b/Rddust
         tauv=taur
         tauzrat=sqrt(pihalf)*a1*mu0/((rmaxd/rddust)**a1-1.)
         tauRos=tauv*chiR(1600./11605.,1)   !setting idust2=1, disk
         tauRth=tauzrat*tauRos
c         taudisk=10.0
         taudisk=min(10.d0,tauRth)
         print*,'taudisk',taudisk
         
         print*,'tauRos,tauv,tauRth',tauRos,tauv,tauRth
         tauw=0.5
         dr=(rmax-rming1)/(dfloat(nrg-irbeg))**rexp
         rbeg2=rming1+dr
         print*,'rbeg2',rbeg2
         if (tauRos.gt.0.51) then
         rfac=(1.-(tauw*(1.-(rmaxd/rddust)**a1))/tauRos)**(1./a1)
         if (rfac.lt.1.0000001d0) then
            print*,''
            print*,'WARNING: not enough precision to resolve inner disk'
            print*,'rfac = ',rfac
c            rfac=rfac*2.d0
            rfac=1.000001d0
             print*,'resetting rfac = ',rfac
            print*,''
         endif

         mu=sqrt(2.0)*ierfc(sngl(taudisk/tauRth))
         rtop=Rddust*(1.-(taudisk*(1.-(rmaxd/rddust)**a1))/
     +        (tauRos*exp(-0.5*mu*mu))      )**(1./a1)
c     +        (tauRos)      )**(1./a1)
         
         print*,'Rddust,taudisk,rmaxd/rddust,tauRos,a1',
     1        Rddust,taudisk,rmaxd/rddust,tauRos,a1
         print*,'rfac,rtop,mu',rfac,rtop,mu
         ir=irbeg
         rarr(ir)=rddust
         if (rfac*rddust.lt.rbeg2) then
         print*,'making gaussian spacing in r, rtop',rtop
         do while ((rarr(ir).lt.rtop).and.(ir.lt.(nrg-1)))
            rarr(ir+1)=rfac*rarr(ir)
            rfac=rfac**1.1
            ir=ir+1
c            print*,'rarr',rarr(ir),rtop,ir,rfac
         end do
         if(ir.ge.nrg) then
            write(*,*) 'too many points in rgrid'
            stop
         endif
         endif
         else
	   ir=irbeg
         endif

         nrwall=ir
         print*,'nrwall,rarr(nrwall)',nrwall,rarr(nrwall)
         
         rming1=rarr(nrwall)
c S.C.         dr=(rmax-rming1)/(dfloat(nrg-nrwall))**rexp
c S.C.        print*,'rexp,dr,rmin',rexp,dr/autors,rming1

         do ir=nrwall+1,nrg
c S.C.            rarr(ir)=rming1+dr*(dfloat(ir-nrwall))**rexp

            rarr(ir)=rming1*(exp((log(rmax/rming1))*
     $         ((dfloat(ir-nrwall))/(nrg-2))))

         enddo
         
      endif

c     set rmine to a grid location
      if (gflag.eq.1) then
         call locate(rarr,nrg,rmine,ir)
         rarr(ir+1)=rmine
      endif

c     r-squared array
      do ir=1,nrg
         r2arr(ir)=rarr(ir)**2
c     print first few points of rarr
         if (ir.lt.10) then
            print*,'rarr/rmin,ir ',rarr(ir)/rmin,ir
         endif
      enddo

c     r-ave array
      do ir=1,nrg-1
         ravearr(ir)=0.5*(rarr(ir)+rarr(ir+1))
      enddo

      open(unit=15,file='rarr.dat',status='unknown')
      write(15,*) nrg,' = number of grid points in r'
      write(15,*) 'index     r/rstar  r(au)'
      do ir=1,nrg
         write(15,*) ir,(rarr(ir)),sngl(rarr(ir)/autors)
      enddo
      close(15)


      print*,'rarr(nrg),r2arr(nrg)',rarr(nrg),r2arr(nrg)

c      if (ntg.gt.2) then
c     set up a tmp array going from 0 - 90 from equ. to pole
c     make theta=0 bin 5 degrees wide (otherwise, too much noise,
c     nothing happens there anyway)
         thettop=5.*pi/180.
         if (massd.eq.0.d0) then

            ntatm=1
            thetatm=0.d0
            tmptharr(1)=0.d0
            nttmp=(ntg+1)/2

         else

c     taufmin=min(0.1,0.1*kappaf*tauv)
c	      use idust2=2, for kappaf, for disk dust properties
            taufmin=min(.001d0,0.001d0*kappaf(2)*tauv)
            
            print*,'kappaf*tauv',kappaf(2)*tauv
            
            mu=min(mu0*sqrt(2.d0*dlog(tauv*kappaf(2)/taufmin)),0.5d0)
            print*,'tauv,kappaf,taufmin',tauv,kappaf(2),taufmin
            
c     ************
c     for disk-only model, let mu=0.5 to sample entire disk height
c     at high-res
c     mu=0.5
c     ***********

            thetatm=pihalf-acos(mu)
c     test
            print*,'thetatm',thetatm*180./pi
            nttmp=(ntg+1)/2

c     *********************
c     change this for disk or envelope runs
c     for envelope with 1 degree resolution in polar region
c     except for first bin (theta=5).
            res=1.
c     for disk with nothing in the polar region
c     res = size of angle bin in poles, approximately 
c     res=10.
c     *********************

            ntatm=nttmp-(pihalf-thetatm-thettop)*90./pihalf/res

            thetmin=min(0.1d0*asin(mu0),thetatm/ntatm)
            
            tmptharr(1)=0.d0
            do it=2,ntatm
               tmptharr(it)=thetmin *
     1              (thetatm/thetmin)**((it-2.d0)/(ntatm-2.d0))
c     print*,'tmptharr ',tmptharr(it)*180.d0/pi
            enddo            
         endif
         do it=ntatm+1,nttmp-1
            tmptharr(it)=thetatm+
     1           (it-ntatm)*(pihalf-thetatm-thettop)/(nttmp-ntatm)
c            print*,'tmptharr(it)',tmptharr(it)
         end do
         tmptharr(nttmp)=pihalf
         
         thetarr(nttmp)=pihalf
         do it=2,nttmp-1
            thetarr(nttmp+1-it)=pihalf-tmptharr(it)
            thetarr(nttmp-1+it)=pihalf+tmptharr(it)
         enddo
         thetarr(1)=0.d0
         thetarr(ntg)=pi
         
c     this is not the eps used in the rest of the code!  see
c     newdisktherm for that
         eps=1.d-8
         do it=1,ntg
            if (it.eq.1) then
               thetarr(it)=0.d0
               costarr(it)=1.d0
               sintarr(it)=0.d0
               tan2arr(it)=0.d0
            else if (it.eq.ntg) then
               thetarr(it)=pi
               costarr(it)=-1.d0
               sintarr(it)=0.d0
               tan2arr(it)=0.d0
            else if (it.eq.nttmp) then
               thetarr(it)=pihalf
               costarr(it)=0.d0
               sintarr(it)=1.d0
               tan2arr(it)=-1.d0
            else
               costarr(it)=cos(thetarr(it))
               sintarr(it)=sin(thetarr(it))
               tan2arr(it)=tan(thetarr(it))**2
            endif
c            print*,'thetarr,costarr,tan2arr '
c            print*,'thetarr',thetarr(it)*180.d0/pi
         enddo
c      else
c         thetarr(1)=0.d0
c         costarr(1)=1.d0
c         sintarr(1)=0.d0
c         tan2arr(1)=0.d0
c      endif

c     theta-ave array
      do it=1,ntg-1
         thetavearr(it)=0.5*(thetarr(it)+thetarr(it+1))
      enddo

      open(unit=15,file='tharr.dat',status='unknown')
      write(15,*) ntg,' = number of grid points in theta'
      write(15,*) 
     1     'index  thet(rad)  thet(deg)    cost     tan**2(thet)'
      do it=1,ntg
         write(15,*) it,sngl(thetarr(it)),sngl(thetarr(it)*180.d0/pi),
     1        sngl(costarr(it)),sngl(tan2arr(it))
      enddo
      close(15)

      dp=r2p/dfloat(npg-1)
      do ip=1,npg
         phiarr(ip)=dp*(dfloat(ip-1))
         aarr(ip)=dsin(phiarr(ip))
         barr(ip)=-dcos(phiarr(ip))
c     print*,'phiarr, aarr, barr ',phiarr(ip),aarr(ip),barr(ip)
c     carr(ip)=
c     darr(ip)=
      enddo

c     phi-ave array
      do ip=1,npg-1
         phiavearr(ip)=0.5*(phiarr(ip)+phiarr(ip+1))
      enddo

      open(unit=15,file='phiarr.dat',status='unknown')
      write(15,*) npg
      do ip=1,npg
         write(15,*) phiarr(ip)*180.d0/pi
      enddo
      close(15)

c     set up diffusion grid
      dcount=0
      if (massd.eq.0.d0) then
         do ir =1,nrg
         do it=1,ntg
         do ip=1,npg
            diffus(ir,it,ip)=.false.
            diffdir(ir,it,ip)=0
         enddo
         enddo
         enddo
      else
         if (rddust.eq.rmin) then 
            nr0=1
         else
            nr0=3
            do it=1,ntg
               do ip=1,npg
                  diffus(1,it,ip) = .false.
                  diffdir(1,it,ip)=0
               end do
            end do
         endif
         do ir=nr0,nrg-1
c            r=0.5*(rarr(ir)+rarr(ir+1))/rddust
            r=ravearr(ir)/rddust
            sigmu=mu0*r**(b-1.)
            Rsonr=1./(rddust*r)
            Tdisk=min(1600.d0/11605.d0,max(3.d0/11605.d0,(Tstar/11605.)
     +           *(max(2./3.*(Rsonr)**3,
     +           (asin(Rsonr)-(Rsonr)*sqrt(1.-(Rsonr)**2)))/pi
c     +        +   (3.*Lacc*(Rsonr)**3*(1.-sqrt(Rsonr)))
     +           )**0.25))
            taumid=tauzrat*tauv*chiR(sngl(Tdisk),1)/r**(a-b)
            if (taumid.gt.taudisk) then
               mudisk=sigmu*sqrt(2.d0)*ierfc(sngl(taudisk/taumid))
            else
               mudisk=0.
            endif
c     print*,'ir,taudisk,taumid,mudisk',ir,taudisk,taumid,mudisk
            do it=1,ntg-1
               mu=abs(cos(0.5d0*(thetarr(it)+thetarr(it+1))))
               tauRmu=taumid*erfc(sngl(mu/(sigmu*sqrt(2.0))))
               tauRr=tauRos*exp(-0.5d0*mu*mu/(sigmu*sigmu))*
     +              (1.d0-r**a1)/(1.d0-(rmaxd/rddust)**a1)
               do ip=1,npg-1
                  if ((mu.lt.mudisk).and.(tauRr.gt.taudisk)) then
                     diffus(ir,it,ip) = .true.
                     if (tauRr.lt.tauRmu) then
                        diffdir(ir,it,ip) = -1
                     else
                        if (thetarr(it).lt.pihalf) then
                           diffdir(ir,it,ip) = -2
                        else
                           diffdir(ir,it,ip) = 2
                        endif
                     endif
                     dcount=dcount+1
c     print*,'true'
                  else
                     diffus(ir,it,ip) = .false.
                     diffdir(ir,it,ip) = 0
c     print*,'mu,mudisk,tauRr,taudisk'
c     print*,mu,mudisk,tauRr,taudisk
c     print*,'false'
                  endif
               end do
            end do
         end do
      endif

c     testing
c      do ir =1,nrg
c         do it=1,ntg
c            do ip=1,npg
c               diffus(ir,it,ip)=.false.
c            enddo
c         enddo
c      enddo
c      dcount=0
c      print*,
c     $ 'WARNING, TURNED OFF DIFFUSION, see line 378-387 of gridset.f'

ctest
c      do ir =1,nrg
c         do it=1,ntg
c            do ip=1,npg
c               if (abs(diffdir(ir,it,ip)).eq.1) then
c                  diffus(ir,it,ip)=.false.
c                  diffdir(ir,it,ip)=0
c               endif
c            enddo
c         enddo
c      enddo 
     
      print*,'number of diffusion cells in grid',dcount

      open(unit=15,file='diffuse.unf',status='unknown',
     1     form='unformatted')
      write(15) nrg-1,ntg-1,npg-1
      write(15) (ravearr(ir)/autors,ir=1,nrg-1)
      write(15) (thetavearr(it),it=1,ntg-1)
      write(15) (phiavearr(ip),ip=1,npg-1)
      write(15) (((diffus(ir,it,ip),ir=1,nrg-1),it=1,ntg-1)
     1     ,ip=1,npg-1)
      close(15)

      open(unit=15,file='diffdir.unf',status='unknown',
     1     form='unformatted')
      write(15) nrg-1,ntg-1,npg-1
      write(15) (ravearr(ir)/autors,ir=1,nrg-1)
      write(15) (thetavearr(it),it=1,ntg-1)
      write(15) (phiavearr(ip),ip=1,npg-1)
      write(15) (((diffdir(ir,it,ip),ir=1,nrg-1),it=1,ntg-1)
     1     ,ip=1,npg-1)
      close(15)

c      print*,'massd',massd
c      stop

c     calculate density in grid

      massenv=0.d0
      massdisk=0.d0
      maxdens=0.d0

      print*,'IN GRIDSET, Rmax=', rmax
      print*,'rmincgs=',rmincgs

      open (unit=12,file='my_dens.dat',status='old')
  
      do ir=1,nrg-1

c     rad=0.5d0*(rarr(ir)+rarr(ir+1))
         rad=ravearr(ir)
         dr=rarr(ir+1)-rarr(ir)
         do it=1,ntg-1
            thet=0.5d0*(thetarr(it)+thetarr(it+1))
            cost=cos(thet)
            sint=sin(thet)
            dt=thetarr(it+1)-thetarr(it)
            do ip=1,npg-1       !3-D atmosphere
               phi=0.5d0*(phiarr(ip)+phiarr(ip+1))
               dp=phiarr(ip+1)-phiarr(ip)

               read (12,*) dense

               dense=(10**dense)*1.67d-24

               print *,dense

               vol=rad**2*sint*dt*dp*dr*rmincgs**3

               massenv=massenv+dense*vol

               call densdisk(rad,sint,cost,phi,densd,idd)
               massdisk=massdisk+densd*vol
               densarr(ir,it,ip)=dense+densd
               if (densarr(ir,it,ip).gt.maxdens) 
     $              maxdens=densarr(ir,it,ip)
               massarr(ir,it,ip)=(dense+densd)*vol
               massarr(ir,it,ip)=massarr(ir,it,ip)**0.25
c     if (idd.eq.0) then
c     dustarr(ir,it,ip)=ide
c     else
               if (dense.gt.densd) then
                  dustarr(ir,it,ip)=ide
               else
                  dustarr(ir,it,ip)=idd
               endif
c     endif
            end do
         enddo
      enddo

      close(12)

      massenv=massenv/1.989d33      

      massdisk=massdisk/1.989d33
c     important note:  densarr(ir,it,ip) is the density halfway between
c     ir and ir+1, it and it+1, ip and ip+1; it is NOT the density at
c     ir,it,ip.

      print*, 'massenv (in solar units) ',massenv
      print*, 'massdisk from grid, compared to input ',massdisk,
     $     massd
      print*,'max density in disk/envelope',maxdens
      print*,'min radius where rho=rhoamb ',radamb/autors
      print*, ' '

c	TODO rescale density to input mass

c     write out rarr,thetarr,densarr at phi=0 so you can read them into
c     IDL.

c	check dustarr
	do ir=1,nrg-1
	do it=1,ntg-1
	do ip=1,npg-1
           if (rarr(ir).lt.rmine.and.rarr(ir).lt.rmind) 
     $          dustarr(ir,it,ip)=4
	   if (dustarr(ir,it,ip).eq.0) then
	     if (ir.ne.nrg.and.it.ne.ntg) then
	        print*,'oops, dustarr,r,th,phi',
     $	     rarr(ir),thetarr(it),phiarr(ip)	
            endif
            dustarr(ir,it,ip)=4    !outflow dust    
	   endif
	enddo
	enddo
	enddo
    
      open(unit=15,file='Av_thet_phi.dat',status='unknown')
c     calculate average optical depth (A_v) along all theta directions
      write(15,*) 'Av along all viewing directions'
      write(15,*) 'theta      phi       Av_env     Av_disk      Av_tot '
      tauave=0.d0
      dcount=0
      dphi=r2p/float(nph)
      do i=1,nmu
         thet=acos(u(i))
         sint=sin(thet)
         call locate(thetarr,ntg,thet,it)
         if (ntg.eq.1) it=1  
c         print*,'it',it
      do j=1,nph
         phi=float(j)*dphi-0.5*dphi
      
         if (phi.gt.r2p) print*,'oops!  phi',phi
         call locate(phiarr,npg,phi,ip)
         tau=0.d0
         taud=0.d0
         taue=0.d0
         cost=u(i)
         do ir=1,nrg-1
     
            dr=rarr(ir+1)-rarr(ir)
c     rad=0.5d0*(rarr(ir)+rarr(ir+1))
            rad=ravearr(ir)
            call densenv(rad,thet,cost,sint,pi,pihalf,phi,dense,
     &           radamb,ide,x_central,y_central,z_central,
     $        N_clump,clumpDensity,clumpsz)
            call densdisk(rad,sint,cost,phi,densd,idd)
c            dense=densarr(ir,it,ip)
            id=dustarr(ir,it,1)
            taud=taud+kapd(id)*(densd)*dr
            taue=taue+kapd(id)*(dense)*dr
         end do
         write(15,*) sngl(thet*180.d0/pi),sngl(phi*180.0/pi)
     1        ,sngl(taue*1.086)
     1        ,sngl(taud*1.086),sngl(taue+taud)*1.086
         tauave=tauave+taud+taue
         dcount=dcount+1
      end do
      end do
      tauave=tauave/dfloat(dcount)
      print*,'dcount',dcount
      print*,'tau_V average over all directions ',tauave
      write(12,*) 'A_V ave over all directions ',tauave*1.086
      close(15)

c     integrate optical depth along all the grid angle directions
      open(unit=15,file='Avall.dat',status='unknown')
      write(15,*) 'Av through each grid theta bin '
      write(15,*) 
     1 '    theta     phi     Av_env     Av_disk     Av_tot'
      do it=1,ntg-1
      do ip=1,npg-1
         tau=0.d0
         taud=0.d0
         taue=0.d0
         thet=thetavearr(it)
         if (thet.gt.pi/2.-eps.and.thet.lt.pi/2.+eps) then
c     TSC blows up at thet=90
            print*,'thet',thet*180./pi
            thet=89.999d0*pi/180.d0
         endif
         cost=cos(thet)
         sint=sin(thet)
         phi=phiarr(ip)
         do ir=1,nrg-1

            dr=rarr(ir+1)-rarr(ir)
c            rad=0.5d0*(rarr(ir)+rarr(ir+1))
            rad=ravearr(ir)

            call densenv(rad,thet,cost,sint,pi,pihalf,phi,dense,
     &           radamb,ide,x_central,y_central,z_central,
     $        N_clump,clumpDensity,clumpsz)

            call densdisk(rad,sint,cost,phi,densd,idd)
            id=dustarr(ir,it,1)
            taud=taud+kapd(id)*(densd)*dr
            taue=taue+kapd(id)*(dense)*dr
         end do
         write(15,*) sngl(thet*180.d0/pi),sngl(phi*180.d0/pi),
     1        sngl(taue)*1.086
     1        ,sngl(taud)*1.086,sngl(taue+taud)*1.086
      end do
      enddo
      close(15)

c     calculate Av along ethet,ephi directions (peeled image)
      open(unit=15,file='A_v.dat',status='unknown')
      write(15,*) 'Av along peeled direction ethet'
      write(15,*) 
     1 'theta   Av_env       Av_disk      Av_tot'
      thet=thete
      cost=cos(thet)
      sint=sin(thet)
      call locate(thetarr,ntg,thet,it)
      phi=phie
      call locate(phiarr,npg,phi,ip)
      taud=0.d0
      taue=0.d0

      do ir=1,nrg-1

         dr=rarr(ir+1)-rarr(ir)
c         rad=0.5d0*(rarr(ir)+rarr(ir+1))
         rad=ravearr(ir)
         id=dustarr(ir,it,1)

         call densenv(rad,thet,cost,sint,pi,pihalf,phi,dense,
     1	radamb,ide,x_central,y_central,z_central,
     $        N_clump,clumpDensity,clumpsz)

         call densdisk(rad,sint,cost,phi,densd,idd)
         taud=taud+kappav(id)*rmincgs*(densd)*dr
         taue=taue+kappav(id)*rmincgs*(dense)*dr
      end do
      write(15,*) sngl(thet*180.d0/pi),sngl(taue*1.086)
     1     ,sngl(taud*1.086),sngl((taue+taud)*1.086)
      write(15,*) 'average A_v along *all* directions ',tauave*1.086
      close(15)

      r=1.d0*autors
      call locate(rarr,nrg,dble(r),ir)
      dens1=densarr(ir,1,1)
      print*,'r,dens1',r,dens1

c     write out massarr
      open(unit=15,file='marr.unf',status='unknown',form='unformatted')
      write(15) nrg-1,ntg-1,npg-1
      write(15) (ravearr(ir)/autors,ir=1,nrg-1)
      write(15) (thetavearr(it),it=1,ntg-1)
      write(15) (phiavearr(ip),ip=1,npg-1)
      write(15) (((massarr(ir,it,ip),ir=1,nrg-1),it=1,ntg-1)
     1     ,ip=1,npg-1)
      close(15)

c     write out densarr
      open(unit=15,file='darr.unf',status='unknown',form='unformatted')
      write(15) nrg-1,ntg-1,npg-1
      write(15) (ravearr(ir)/autors,ir=1,nrg-1)
      write(15) (thetavearr(it),it=1,ntg-1)
      write(15) (phiavearr(ip),ip=1,npg-1)
      write(15) (((densarr(ir,it,ip),ir=1,nrg-1),it=1,ntg-1)
     1     ,ip=1,npg-1)
      close(15)
      
c	write out dustarr
      open(unit=15,file='dust.unf',status='unknown',form='unformatted')
      write(15) nrg-1,ntg-1,npg-1
      write(15) (ravearr(ir)/autors,ir=1,nrg-1)
      write(15) (thetavearr(it),it=1,ntg-1)
      write(15) (phiavearr(ip),ip=1,npg-1)
      write(15) (((dustarr(ir,it,ip),ir=1,nrg-1),it=1,ntg-1)
     1     ,ip=1,npg-1)
      close(15)

      end


      







