; general syntax:
;
; readnew,filename,varname,lable,type=type,ndim=ndim,is_present=is_present,
;                                debug=debug,status=status,parttype=parttype
;     filename  : string
;     varname   : variable
;     lable     : string (4 character) lable of block to read
;     type      : string defining type: HEAD,LONG,FLOAT,FLOATN (default=FLOAT)
;     ndim      : dimention for FLOATN (default=3)
;     is_present: array (6 int) defining the presence for different
;                 types (default=[1,0,0,0,0,0])
;     debug     : gives debugging output if set
;     status    : returns 0/1 for if block is absent/present
;
; For several known blocks, type,ndim and is_present is automaticaly
; set. Also the default allows to read any additional float without
; specifying any thing additional.
;
; examples:
; (1) reading header informations (see arround line 375 for the
;                                  content of the headear structure):
; readnew,'snap_092',h,'HEAD'
; print,h.redshift,h.npart(0),h.npart(5),'...'
; 
; (2) examin the content of a file:
; readnew,'snap_092',x,'XXXX',/debug
;
; (3) reading positions
; readnew,'snap_092',x,'POS'
;    => x is a 3 times npart_total array
;
; (4) reading only star positions:
; readnew,'snap_092',x,'POS',parttype=5
;    => x is a 3 times npart(5) array
;
; Note: 
; (a) Supports splitted file reading (just give the file name without the
;     .number extension), works also for the pattype option !
; (b) If the MASS block is read, it gives back an filled array for all 
;     particles
;
; Feel free to extend the 'IF NOT(IS_DEF(TYPE)) THEN BEGIN' branch for
; any particular new blocks ...
 
FUNCTION IS_DEF, x
aux = SIZE(x)
RETURN, aux(N_ELEMENTS(aux)-2) NE 0
END

PRO readnew,name,x,label,type=type,debug=debug,quiet=quiet,parttype=parttype,$
            status=status,ndim=ndim,is_present=is_present
@setconst
   status=-1
   IF NOT(IS_DEF(name)) THEN BEGIN
      print,'Usage: readnew,name,x,label,type=type,debug=debug'
      print,'       name  : filename'
      print,'       x     : variable containing the result'
      print,'       label : names the block identifier ("HEAD","POS ",...)'
      print,'       type  : names the data type ("FLOAT","FLOAT3","LONG","HEAD")'
      print,'       debug : give some debug information'
      return
   END

   blabel=BYTE(label+"    ")
   label=STRING(blabel(0:3))
   IF NOT IS_DEF(is_present) THEN is_present=[1,0,0,0,0,0]
   IF NOT(IS_DEF(TYPE)) THEN BEGIN
      IF strcmp(label,"HEAD")  THEN  type="HEAD"
      IF strcmp(label,"ID  ")  THEN BEGIN
         type="LONG"
         is_present=[1,1,1,1,1,1]
      END
      IF strcmp(label,"POT ")  THEN BEGIN
         type="FLOAT"
         is_present=[1,1,1,1,1,1]
      END
      IF strcmp(label,"POS ") OR strcmp(label,"VEL ") OR strcmp(label,"ACCE") THEN BEGIN
         type="FLOATN"
         ndim=3
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"TIPS") THEN BEGIN
         type="FLOATN"
         ndim=9
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"DIPS") THEN BEGIN
         type="FLOATN"
         ndim=36
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"CACO") THEN BEGIN
         type="FLOATN"
         ndim=1
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"FLDE") THEN BEGIN
         type="FLOATN"
         ndim=1
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"STDE") THEN BEGIN
         type="FLOATN"
         ndim=1
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"PSDE") THEN BEGIN
         type="FLOATN"
         ndim=1
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"ANRA") THEN BEGIN
         type="FLOATN"
         ndim=2
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"LACA") THEN BEGIN
         type="FLOATN"
         ndim=20
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"SHIN") THEN BEGIN
         type="FLOATN"
         ndim=3
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"SHOR") THEN BEGIN
         type="FLOATN"
         ndim=9
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"INDE") THEN BEGIN
         type="FLOATN"
         ndim=1
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"BFLD") OR strcmp(label,"BFSM") THEN BEGIN
         type="FLOATN"
         ndim=3
         is_present=[1,0,0,0,0,0]
      END
      IF strcmp(label,"Zs  ")  THEN BEGIN
         type="FLOATN"
         IF NOT(IS_DEF(ndim)) THEN ndim=7
         is_present=[1,0,0,0,1,0]
      END   
      IF strcmp(label,"AGE ") OR strcmp(label,"iM  ") OR $
         strcmp(label,"SLg ") THEN BEGIN
         type="FLOAT"
         is_present=[0,0,0,0,1,0]
      END
      IF NOT(IS_DEF(TYPE)) THEN type="FLOAT"
   END
   IF NOT(IS_DEF(ndim)) THEN ndim=3


   b4=bytarr(4)
   bl=0L

   get_lun,myfile
   myname=name
   numfiles=1

   ss=SIZE(FINDFILE(myname))
   IF ss(0) EQ 0 THEN BEGIN
      myname=name+'.0'
      ss=SIZE(FINDFILE(myname))
      IF ss(0) EQ 0 THEN BEGIN
         print,'Cant find file ',name,' or ',myname,' ...'
         stop
      END ELSE BEGIN
         npart=lonarr(6)
         readnew,myname,h,'HEAD',debug=debug,quiet=quiet
         npart=npart+h.npart
         WHILE ss(0) NE 0 DO BEGIN
            myname=name+'.'+strcompress(string(numfiles),/remove_all)
            ss=SIZE(FINDFILE(myname))
            IF ss(0) NE 0 THEN BEGIN
               numfiles=numfiles+1
               readnew,myname,h,'HEAD',debug=debug,quiet=quiet
               npart=npart+h.npart
            END
         END
         IF IS_DEF(debug) THEN print,'Found',numfiles,' files ...'
         IF IS_DEF(debug) THEN print,'TotNP',npart
      END
   END

   IF numfiles GT 1 THEN BEGIN
      IF strcmp(type,"HEAD") THEN BEGIN
         h.npart(*)=npart(*)
         x=h
      END ELSE BEGIN
         IF IS_DEF(parttype) THEN BEGIN
            ntot=npart(parttype)
         END ELSE BEGIN
; Keep explicite sum, as TOTAL(npart) would be a real !!!!
            ntot=npart(0)+npart(1)+npart(2)+npart(3)+npart(4)+npart(5)
         END

         IF strcmp(label,"POS ") OR strcmp(label,"ID  ") OR $
            strcmp(label,"VEL ") OR strcmp(label,"MASS") OR strcmp(label,"ACCE")$
            THEN nalloc=ntot ELSE nalloc=npart(1)
         IF strcmp(label,"TIPS") OR strcmp(label,"DIPS") OR $
            strcmp(label,"CACO") OR strcmp(label,"FLDE") OR $
            strcmp(label,"PSDE") OR strcmp(label,"ANRA") OR $ 
            strcmp(label,"LACA") OR strcmp(label,"SHOR") OR $
            strcmp(label,"INDE") OR strcmp(label,"SHIN") OR $ 
            strcmp(label,"STDE") $
            THEN nalloc=ntot
         IF strcmp(label,"AGE ") OR strcmp(label,"iM  ") OR $
            strcmp(label,"SLg ") THEN nalloc=npart(4)
         IF strcmp(label,"Zs ") THEN nalloc=npart(0)+npart(4)

         IF strcmp(type,"FLOAT")  THEN x=fltarr(nalloc) 
         IF strcmp(type,"FLOATN") THEN x=fltarr(ndim,nalloc)
         IF strcmp(type,"LONG")   THEN x=lonarr(nalloc)

         nstart=0L
         FOR i=0L,numfiles-1 DO BEGIN
            myname=name+'.'+strcompress(string(i),/remove_all)
            readnew,myname,htmp,'HEAD',debug=debug,quiet=quiet,status=status
            gof=1
            IF IS_DEF(parttype) THEN BEGIN
               IF htmp.npart(parttype) EQ 0 THEN gof=0
            END
            IF gof EQ 1 THEN BEGIN
               rfile=1
               IF strcmp(label,"MASS") THEN BEGIN
                  IF IS_DEF(parttype) THEN BEGIN
                     IF (htmp.npart(parttype) EQ 0) THEN BEGIN
                        rfile=0
                     END ELSE BEGIN
                        IF (htmp.massarr(parttype) GT 0) THEN BEGIN
                           rfile=0
                           xpart=fltarr(htmp.npart(parttype))
                           xpart(*)=htmp.massarr(parttype)
                        END
                     END
                  END ELSE BEGIN
                     nmassread=0L
                     nmasstable=0L
                     FOR im=0,5 DO BEGIN
                        IF htmp.npart(im) GT 0 THEN BEGIN
                           nmasstable=nmasstable+htmp.npart(im)
                           IF htmp.massarr(im) EQ 0 THEN nmassread=nmassread+1
                        END
                     END
                     IF nmassread EQ 0 THEN BEGIN
                        rfile=0
                        xpart=fltarr(nmasstable)
                        istart=0L
                        FOR im=0,5 DO BEGIN
                           IF htmp.npart(im) GT 0 THEN BEGIN
                              xpart(istart:istart+htmp.npart(im)-1) = htmp.massarr(im) 
                              istart=istart+htmp.npart(im)
                           END
                        END
                     END
                  END
               END
               IF rfile EQ 1 THEN $
                   readnew,myname,xpart,label,type=type,debug=debug,$
                           quiet=quiet,parttype=parttype,status=status,$
                           ndim=ndim
               ss=SIZE(xpart)
               IF strcmp(type,"FLOATN") THEN BEGIN
                  FOR ip=0L,ss(2)-1 DO BEGIN
                     FOR k=0,ndim-1 DO BEGIN
                        x(k,nstart)=xpart(k,ip)
                     END
                     nstart=nstart+1
                  END
                  IF IS_DEF(debug) THEN $
                      print,'Reading data',nstart,nstart+ss(2)-1,' ...'
               END ELSE BEGIN
                  x(nstart:nstart+ss(1)-1)=xpart(*)
                  nstart=nstart+ss(1)
                  IF IS_DEF(debug) THEN $
                      print,'Reading data',nstart,nstart+ss(1)-1,' ...'
               END
            END ELSE BEGIN
               IF IS_DEF(debug) THEN $
                   print,'Skipping file, does not contain particle of type ',parttype
            END
         END
      END
   END ELSE BEGIN


      double_flag = 0L 
      IF IS_DEF(debug) THEN print,'Testing File ',myname
      openr,myfile,myname
      readu,myfile,bl
      skip_lun,myfile, 16L + 256L - 60L 
      readu,myfile, double_flag 
      close,myfile

      IF double_flag EQ 1 THEN BEGIN
         IF IS_DEF(debug) THEN print, 'Detected double precision file'
      ENDIF

      IF bl EQ 8 THEN BEGIN
         openr,myfile,name,/f77 
         IF IS_DEF(debug) THEN print,'Open file in normal mode ..,'
      END ELSE BEGIN
         openr,myfile,name,/f77,/SWAP_ENDIAN
         IF IS_DEF(debug) THEN print,'Open file with SWAP_ENDIAN ..,'
      END

      found=0

      WHILE(NOT(EOF(myfile)) AND found EQ 0) DO BEGIN
         readu,myfile,b4,bl
         thislabel=STRING(b4)
         IF IS_DEF(debug) THEN print,thislabel," <=> ",label,bl
         IF strcmp(thislabel,'HEAD') NE 0 THEN BEGIN
            IF bl NE 264 THEN BEGIN
               bl=264L
               print,'Patching length for HEAD !!!'
            END 
         END

         IF strcmp(thislabel,label) EQ 0 THEN BEGIN
            point_lun,-myfile,pos
            point_lun,myfile,pos+bl
         END ELSE BEGIN
            found=1
         END
      ENDWHILE

      IF found EQ 0 THEN BEGIN
         IF NOT(IS_DEF(quiet)) THEN $
           print,'File ',name,' does not contain a block labled with "',label,'" !'

         xH=0.76

         IF strcmp(label,"TEMP") THEN BEGIN
            IF NOT(IS_DEF(quiet)) THEN $ 
              print,"WARNING: calculting temperature, assuming no starformation !!!"
            readnew,name,u,'U   ',quiet=quiet,parttype=parttype

            g1 = 5.0 / 3.0 - 1.

            readnew,name,xNe,'NE  ',quiet=quiet,parttype=parttype
            IF NOT(IS_DEF(xNe)) THEN xNe=1.0

            yhelium = (1. - xH) / (4. * xH)
            mu = (1. + 4. * yhelium) / (1. + 3. * yhelium + xNe)

            x = FLOAT(g1 / kcgs * mp * u * e_unit / m_unit * mu)
            status=1
         END
         IF strcmp(label,"MHOT") THEN BEGIN
            readnew,name,h,'HEAD'
            IF h.massarr(0) EQ 0.0 THEN BEGIN
               readnew,name,m,'MASS',parttype=parttype
               gmas=m(0:h.npart(0)-1)
            END ELSE BEGIN
               gmas=fltarr(h.npart(0))
               gmas(*)=h.massarr(0) 
            END
            readnew,name,mhi,'MHI ',quiet=quiet,parttype=parttype
            IF IS_DEF(mhi) THEN BEGIN
               IF NOT(IS_DEF(quiet)) THEN print,"Calculting MHOT from MHI ..."
               xfr=mhi/gmas/xH
            END ELSE BEGIN
               xfr=0.
            END
            x = (1 - xfr) * gmas
            status=1
         END
         IF strcmp(label,"RHOT") THEN BEGIN
            readnew,name,rho,'RHO ',parttype=parttype
            readnew,name,mhi,'MHI ',quiet=quiet,parttype=parttype
            IF NOT(IS_DEF(mhi)) THEN BEGIN
               IF NOT(IS_DEF(quiet)) THEN print,"Returning normal density ..."
               x=rho
            END ELSE BEGIN
               IF NOT(IS_DEF(quiet)) THEN print,"Calculting RHOT from MHI ..."
               readnew,name,h,'HEAD'
               IF h.massarr(0) EQ 0.0 THEN BEGIN
                  readnew,name,m,'MASS',parttype=parttype
                  gmas=m(0:h.npart(0)-1)
               END ELSE BEGIN
                  gmas=h.massarr(0) 
               END
               xfr=mhi/gmas/xH
               x = (1 - xfr) * rho
            END
            status=1
         END
         IF strcmp(label,"MCLD") THEN BEGIN
            readnew,name,mhi,'MHI ',quiet=quiet,parttype=parttype
            IF IS_DEF(mhi) THEN BEGIN
               IF NOT(IS_DEF(quiet)) THEN print,"Calculting MCLD from MHI ..."
               readnew,name,h,'HEAD'
               IF h.massarr(0) EQ 0.0 THEN BEGIN
                  readnew,name,m,'MASS',parttype=parttype
                  gmas=m(0:h.npart(0)-1)
               END ELSE BEGIN
                  gmas=h.massarr(0) 
               END
               xfr=mhi/gmas/xH
               x = xfr * gmas
               status=1
            END
         END
         IF strcmp(label,"RCLD") THEN BEGIN
            readnew,name,rho,'RHO ',parttype=parttype
            readnew,name,mhi,'MHI ',quiet=quiet,parttype=parttype
            IF NOT(IS_DEF(mhi)) THEN BEGIN
               IF NOT(IS_DEF(quiet)) THEN print,"Returning normal density ..."
               x=rho
            END ELSE BEGIN
               IF NOT(IS_DEF(quiet)) THEN print,"Calculting RHOT from MHI ..."
               readnew,name,h,'HEAD'
               IF h.massarr(0) EQ 0.0 THEN BEGIN
                  readnew,name,m,'MASS',parttype=parttype
                  gmas=m(0:h.npart(0)-1)
               END ELSE BEGIN
                  gmas=h.massarr(0) 
               END
               xfr=mhi/gmas/xH
               x = xfr * rho
            END
            status=1
         END
      END ELSE BEGIN
         IF strcmp(type,"HEAD")   THEN BEGIN
            npart=lonarr(6)	
            massarr=dblarr(6)
            time=0.0D
            redshift=0.0D
            flag_sfr=0L
            flag_feedback=0L
            bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4
            la=bytarr(bytesleft)
            readu,myfile,npart,massarr,time,redshift,flag_sfr,flag_feedback,la
            x = { head , npart:npart,$
                         massarr:massarr,$
                         time:time,$
                         redshift:redshift,$
                         flag_sfr:flag_sfr,$
                         flag_feedback:flag_feedback,$
                         la:la}
         END ELSE BEGIN
            status=1
            IF strcmp(label,"MASS") THEN BEGIN
               readnew,name,h,'HEAD'
               ntotm=h.npart(0)+h.npart(1)+h.npart(2)+h.npart(3)+h.npart(4)+h.npart(5)
               x=fltarr(ntotm)
               ninput=(bl-8)/4
               xr=fltarr(ninput)
               if double_flag eq 1 then begin
                 ninput=(bl-8)/8
                 IF strcmp(type,"LONG") then ninput=(bl-8)/4
                 x=dblarr(ntotm)
                 xr=dblarr(ninput)
               endif               
               readu,myfile,xr
               st=0L
               str=0L
               FOR i=0,5 DO BEGIN
                  IF h.npart(i) GT 0 THEN BEGIN
;                     print,'filling ',st,st+h.npart(i)-1
                     IF h.massarr(i) GT 0 THEN BEGIN
                        x(st:st+h.npart(i)-1)=h.massarr(i)
                        st=st+h.npart(i)
                     END ELSE BEGIN
                        x(st:st+h.npart(i)-1)=xr(str:str+h.npart(i)-1)
                        st=st+h.npart(i)
                        str=str+h.npart(i)
                     END
                  END
               END
               IF IS_DEF(parttype) THEN BEGIN
                  IF h.npart(parttype) GT 0 THEN valid=1 ELSE valid=0
              END ELSE valid=1
            END ELSE BEGIN
               IF IS_DEF(parttype) THEN BEGIN
                  readnew,name,h,'HEAD'
                  IF h.npart(parttype) GT 0 THEN valid=1 ELSE valid=0
               END ELSE valid=1

               IF strcmp(type,"FLOATN") THEN ninput=(bl-8)/4/ndim $
                                        ELSE ninput=(bl-8)/4
               if double_flag eq 1 then begin
                 IF strcmp(type,"FLOATN") THEN ninput=(bl-8)/8/ndim $
                                          ELSE ninput=(bl-8)/8
                 IF strcmp(type,"LONG") then ninput=(bl-8)/4                                          
               endif

               IF valid EQ 1 THEN BEGIN                      
                  IF strcmp(type,"FLOAT")  THEN x=fltarr(ninput) 
                  IF strcmp(type,"FLOATN") THEN x=fltarr(ndim,ninput)
                  IF strcmp(type,"LONG")   THEN x=lonarr(ninput) 
                  if double_flag eq 1 then begin
                   IF strcmp(type,"FLOAT")  THEN x=dblarr(ninput) 
                   IF strcmp(type,"FLOATN") THEN x=dblarr(ndim,ninput)
                  endif
               END
               IF valid EQ 1 THEN readu,myfile,x
            END
;            IF IS_DEF(x) THEN BEGIN
            IF valid EQ 1 THEN BEGIN
               IF IS_DEF(parttype) THEN BEGIN
                  IF ninput NE h.npart(parttype) THEN BEGIN
                     pstart=0L
                     FOR i=0,parttype-1 DO pstart=pstart+h.npart(i)*is_present(i)
                     IF NOT(IS_DEF(quiet)) THEN $
                        print,'Restricting return value to',pstart,pstart+h.npart(parttype)-1
                     IF strcmp(type,"FLOATN") THEN $
                        x=x(*,pstart:pstart+h.npart(parttype)-1) $
                     ELSE $
                        x=x(pstart:pstart+h.npart(parttype)-1) 
                  END
               END
            END ELSE BEGIN
               status=-1
               IF IS_DEF(parttype) THEN $
                  print,'File contains no "',type,'" for particle oy type',parttype,' !' $
               ELSE $
                  print,'Unknown type "',type,'" !'
            END
         END
      END

      close,myfile
   END
   free_lun,myfile
 END
   
