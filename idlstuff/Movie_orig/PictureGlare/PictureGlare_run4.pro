;call device, DECOMPOSED=0, /BYPASS_TRANSLATION, RETAIN=2 before
;openning any drawing window
;; .rn load_center_and_BHs.pro first

Object_File='../C-Routines_for_IDL/ComputeHsmlAndProject/HsmlAndProject.so'

h=0.7
omega_m=0.3
H_0 = 70                        ;km/s/mpc
t_0=13.7                        ;gyr
;tz14=0.26  ;;z0=15
tz14=0.275  ;;z=14.5
dt=0.002
tbin=0.002

UnitLength_in_cm =        3.085678d21 ;  1.0 kpc 
UnitMass_in_g    =        1.989d43 ;  1.0e10 solar masses 
UnitVelocity_in_cm_per_s = 1d5  ;  1 km/sec 
UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s 
UnitDensity_in_cgs= UnitMass_in_g/ UnitLength_in_cm^3
UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ UnitTime_in_s^2
UnitEnergy_in_cgs= UnitMass_in_g * UnitLength_in_cm^2 / UnitTime_in_s^2
GRAVITY   = 6.672d-8 
BOLTZMANN = 1.3806d-16
PROTONMASS = 1.6726e-24
G=GRAVITY/ UnitLength_in_cm^3 * UnitMass_in_g * UnitTime_in_s^2
Xh=0.76D                        ; mass fraction of hydrogen
;HubbleParam= 0.65
HubbleParam= 0.7
gamma= 5.0/3

halo=['test_2','test', 'md15_eos0.25_nowind/', 'md15_eos0.5_nowind/', $
      '/mpich2/md15_eos0.5_nowind/', 'pgdt_0306_mpi2/md15_eos0.5_nowind/' ]

halost=['md15_eos0.5_G4','md15_eos0.5_g4', 'md15_eos0.25_GM8', 'md15_eos0.5_g8', $
        'md15_eos0.5_gM8', 'md15_eos0.5_GM8' ]

din='/raid5/yxli/Runs/merger_tree/run4/ICs/1gpc_zoom1/'
dout='/raid5/yxli/Runs/merger_tree/run4/Sims/1gpc_zoom1/maps_fancy/'  

nh=n_elements(halo)

for ih=0, 0 do begin

    label=halost[ih]

    ddir='/raid5/yxli/Runs/merger_tree/run4/Sims/1gpc_zoom1/'+halo[ih]+'/'

    PicsDir = dout+halost[ih]+'/figure/'
         
    fcom=din+'combine_zoom1_0.5rd_1r200_t0_md15.dat'
    if ih eq 0 then fcom=din+'combine_zoom1_0.5rd_1r200_t0_md15_2.dat'
    readcol, fcom, gal, tb, tm, ns, $
      format='(x, a, x,x,x,x,x,x,x,x,x,x, f,f, a)', comment='#'

    ng=n_elements(gal)
    ns=fix(ns)

    if ih eq 0 then bhlist=[24771, 48747, 87321, 116734, 162614, 200198, 271103, 342757]
    if ih eq 1 then bhlist=[24771, 48747, 91920, 119355, 163039, 200537, 271776, 343123]
    if ih eq 3 then bhlist=[24771, 48747, 92021, 119392, 163044, 200607, 271540, 343293]

;    for ig=0, ng-1 do begin
    for ig=4, 4 do begin

        Frun=ddir+gal[ig]+'/'
        nsnap=fix((tm[ig]-tb[ig])/dt)
        print, 'nsnap, ns[ig]: ', nsnap, ns[ig]

        Snap = gal[ig]


        get_time_list, Frun, snap, NumList, TiList


;        PIX=1024
        PIX=768

        xpixels= long(PIX)      ; pixels in x-direction
        ypixels= long(PIX)      ; pixels in y-direction

        Len =   40.0

        seed=42L

        rep = 14

;        for num=0, nsnap do begin
        for num=16, 16 do begin

            if ig eq 0 then begin
                ord=num
            endif else begin    
                ord=fix(total(ns[0:ig-1]))+num
            endelse

;;output
            exts='000'
            exts=exts+strcompress(string(ord),/remove_all)
            exts=strmid(exts,strlen(exts)-3,3)
            fjpg=PicsDir +"pic_star_"+exts+".jpg"
            fjpg=strcompress(fjpg,/remove_all)
            feps=PicsDir +"pic_star_"+exts+".eps"
            feps=strcompress(feps,/remove_all)
;;input
            exts='000'
            exts=exts+strcompress(string(num),/remove_all)
            exts=strmid(exts,strlen(exts)-3,3)
            f=frun + snap+ "_"+exts
            f=strcompress(f,/remove_all)

            npart=lonarr(6)	
            massarr=dblarr(6)
            time=0.0D
            redshift=0.0D
            flag_sfr=0L
            flag_feedback=0L
            npartall=lonarr(6)	
            bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4
            la=intarr(bytesleft/2)

            openr,1,f,/f77_unformatted ;, /swap_endian
            readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,la
;    print,npart,massarr
            print
            print,num, time,redshift
;close,1
            N=total(npart)
            pos=fltarr(3,N)
;vel=fltarr(3,N)
            id=lonarr(N)
            readu,1,pos
            readu,1             ;,vel
            readu,1 ,id
            ind=where((npart gt 0) and (massarr eq 0))
            if ind(0) ne -1 then begin
                Nm= total(npart(ind))
                mass=fltarr(Nm)	
                readu,1,mass
            endif
            
            NGas=npart(0)
            NHalo=npart(1)
            NDisk=npart(2)
            NBulge=npart(3)
            NStars=npart(4)
            NHoles=npart(5)
            
            Nskip = NGas + NHalo + NDisk +  NBulge +  NStars
            
            
            
            u=fltarr(Ngas)
            readu,1,u
            rho=fltarr(Ngas)
            readu,1,rho
            Nelec=fltarr(Ngas)
            readu,1 ,Nelec
;  NH0=fltarr(Ngas)
            readu,1             ;,NH0
            Hsml=fltarr(Ngas)
            readu,1, Hsml
            Sfr= fltarr(Ngas)
            readu,1, Sfr
            close,1
            
            
            MeanWeight= 4.0/(3*Xh+1+4*Xh*Nelec) * PROTONMASS
            Temp = MeanWeight/BOLTZMANN * (gamma-1) * U * UnitEnergy_in_cgs/ UnitMass_in_g
            
;;YL add
            bhlist=[24771, 48747, 87321, 116734, 162614, 200198, 271103, 342757]
            
            load_center_and_BHs, f, bhlist, tz, center, center1, nbh, idbh, mbh, posbh
            print, 'cm, center1, nbh : ', center, center1, nbh
            wb=where(mbh eq max(mbh), cb)
            cbh=fltarr(3)
            cbh[0]=posbh[0,wb[0]]
            cbh[1]=posbh[1,wb[0]]
            cbh[2]=posbh[2,wb[0]]

            bhmas=total(mbh)

            load_stellar_age_met, f, mean_age, mean_smet
            print, 'mean_age, mean_smet', mean_age, mean_smet
            
            tz=tz/h
            zz=(2.0/3.0*t_0/(tz+tz14)/sqrt(omega_m))^(2.0/3.0)-1.0
            
            zstring='z='+ string(format='(f6.2)', zz)
        
            Pos(0,*)=Pos(0,*)-center1[0]
            Pos(1,*)=Pos(1,*)-center1[1]
            Pos(2,*)=Pos(2,*)-center1[2]
;            Pos(0,*)=Pos(0,*)-cbh[0]
;            Pos(1,*)=Pos(1,*)-cbh[1]
;            Pos(2,*)=Pos(2,*)-cbh[2]

;;;
            
            xc = Pos(0,Nskip)   ; position of BH
            yc = Pos(1,Nskip) 
            zc = Pos(2,Nskip)
            
            
            if Ngas gt 0 then begin
                xgas=fltarr(Ngas) &  ygas=fltarr(Ngas)  & zgas=fltarr(Ngas) & mgas=fltarr(Ngas)
                xgas(*)=pos(0,0:Ngas-1)
                ygas(*)=pos(1,0:Ngas-1)
                zgas(*)=pos(2,0:Ngas-1)
                
                if massarr(0) eq 0 then begin
                    mgas(*)=mass(0:Ngas-1)	
                endif else begin
                    mgas(*)= massarr(0)
                endelse
            endif
            
            ind =where((abs(xgas) lt len) and (abs(ygas) lt len) and (abs(zgas) lt len))
            Mgas =total(mgas(ind))
            
            
            if Nstars gt 0 then begin
                xstars=fltarr(Nstars) &  ystars=fltarr(Nstars)  & zstars=fltarr(Nstars)  & mstars=fltarr(Nstars)
                xstars(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
                ystars(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
                zstars(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
                if massarr(4) eq 0 then begin
                    skip=0L
                    for t=0,2 do begin
                        if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                            skip=skip + npart(t)
                        endif
                    endfor
                    mstars(*)=mass(0+skip:Nstars-1+skip)
                endif else begin
                    mstars(*)= massarr(4)
                endelse
            endif
            
            if Ndisk gt 0 then begin
                xdisk=fltarr(NDisk) &  ydisk=fltarr(NDisk) &  zdisk=fltarr(NDisk) & mdisk=fltarr(NDisk)
                xdisk(*)=pos(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
                ydisk(*)=pos(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
                zdisk(*)=pos(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
                if massarr(2) eq 0 then begin
                    skip=0L
                    for t=0,1 do begin
                        if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                            skip=skip + npart(t)
                        endif
                    endfor
                    mdisk(*)=mass(0+skip:Ndisk-1+skip)
                endif else begin
                    mdisk(*)= massarr(2)
                endelse
            endif
            

;;;;;;;;;;;;; start  X-Y- projection of density

            if Ndisk gt 0 then begin
                N=   long(Nstars+Ndisk) ; number of particles
                
                Xyz=  fltarr(3,N) ; generate coordinate array for particles
                Xyz(0,0:nstars-1)= xstars(*) ; fill in the x-coordinates
                Xyz(1,0:nstars-1)= ystars(*) ; ...         y-coordinates
                Xyz(2,0:nstars-1)= zstars(*) ; ...         y-coordinates
                Xyz(0,Nstars:*)= xdisk(*) ; fill in the x-coordinates
                Xyz(1,Nstars:*)= ydisk(*) ; ...         y-coordinates
                Xyz(2,Nstars:*)= zdisk(*) ; ...         y-coordinates

                Weight = fltarr(N) ; generate space for mass array
                Weight(0:nstars-1)= mstars(*) ; fill in mass values
                Weight(nstars:*)= mdisk(*) ; fill in mass values
                
            endif else begin
                N=   long(Nstars) ; number of particles
                
                Xyz=  fltarr(3,N) ; generate coordinate array for particles
                Xyz(0,0:nstars-1)= xstars(*) ; fill in the x-coordinates
                Xyz(1,0:nstars-1)= ystars(*) ; ...         y-coordinates
                Xyz(2,0:nstars-1)= zstars(*) ; ...         y-coordinates

                Weight = fltarr(N) ; generate space for mass array
                Weight(0:nstars-1)= mstars(*) ; fill in mass values
            endelse 

            Hsml = fltarr(N)                
            quantity = fltarr(N)       
            quantity(*)= 1.0    ; fill in mass values
            Axis1=   long(0)    ; select horizontal axis
            Axis2=   long(1)    ; select vertical axis
            Axis3=   long(2)    ; select projection axis
            
            xmin= float(-len)   ; left edge of map
            xmax= float( len)   ; right edge of map 
            ymin= float(-len)   ; lower
            ymax= float( len)   ; upper edge
            Zmin= float(-1.0e6) ; lower range of particles that are included
            Zmax= float(1.0e6)  ; upper range
            
            
            xpixels= long(PIX)  ; pixels in x-direction
            ypixels= long(PIX)  ; pixels in y-direction
            
            Hmax=    float( (xmax-xmin)/xpixels * 50) ; maximum used smoothing length
            BoxSize  = double(0.0)
            
            DesNgb = 32L ; for those particles that have no zero smoothing length
            
            ResultW = fltarr(YPixels, XPixels) ; projected weights (normally projected mass map)
            ResultQ = fltarr(YPixels, XPixels) ; projected `Weight' weighted quantity
                    
            
            ind =where((abs(xyz(0,*)) lt len) and (abs(xyz(1,*)) lt len) and (abs(xyz(2,*)) lt len))
            Mstars =total(Weight(ind))
            
            gasfrac =  mgas/(mgas+mstars)
            print, "gas fraction=", gasfrac
            
            S = CALL_EXTERNAL(Object_File, $
                              'findHsmlAndProject', $
                              N, $
                              Xyz, $
                              Hsml, $
                              Weight, $
                              Quantity, $
                              Xmin, $
                              Xmax, $
                              Ymin, $
                              Ymax, $
                              Zmin, $
                              Zmax, $
                              XPixels, $
                              YPixels, $
                              DesNgb, $
                              Axis1, $
                              Axis2, $
                              Axis3, $
                              Hmax,  $
                              BoxSize, $
                              ResultW, $
                              ResultQ)
            
            ResultW = transpose(ResultW)
            ResultW /= total(ResultW)

;;; make the glaring effect

            
            R = 100
            
            xs= xpixels
            ys= ypixels
            

            pic = fltarr(Xs, Ys)
            
            for bhnr=0, Nholes-1 do begin
                
                N = 200
                
                xc = (Pos(0,Nhalo+Ndisk+Ngas+NBulge+Nstars+bhnr)-xmin)/(xmax-xmin)*xs
                yc=  (Pos(1,Nhalo+Ndisk+Ngas+NBulge+Nstars+bhnr)-ymin)/(ymax-ymin)*ys
                
                print, xc, yc
                
                Nspoke = fix(80*bhmas/0.12)
                if Nspoke le 10 then Nspoke=10
                Radius= 5*bhmas/0.12

                if num gt 200 then begin
                    Nspoke=2                
                    Radius=0.01
                endif

                spoke = randomn(seed, Nspoke)
                
                
                for row = 0, ys-1 do begin
                    print, row
                    for col = 0, Xs-1 do begin
                        
                        u =  double(col - xc) 
                        v =  double(row - yc) 
                        l = sqrt (u * u + v * v)/2
                        
                        t = (atan (u, v) / (2 * !PI) + .51) * Nspoke 
                        i = floor (t) 
                        t -= i 
                        i = i mod Nspoke 
                        
                        w1 = spoke[i] * (1 - t) + spoke[(i + 1)  mod nspoke] * t 
            
                        w1 = w1 * w1 
                        
                        if l lt radius then begin
                            w1 = max(spoke^2)*(1-  l/ (radius))  + l/ (radius) * w1
                        endif
                        
                        w = 1.0 / (l + 2*radius)
                        
                        c =  w1 * w^3.8
                        
                        pic(col, row) += c
                    endfor
                endfor
            endfor
            
            

            Pic = pic /total(pic) * total(resultW)*2.0
            
            resultW *= 0.8
            
            imageTemp2 = ((+alog10(Pic)-alog10(max(resultw)/100) +3  )/4.0)*75.0
            
            ind=where(imageTemp2 lt 0) 
            if ind(0) ne -1 then ImageTemp2(ind)=0
            ind=where(ImageTemp2 gt 75.0)
            if ind(0) ne -1 then Imagetemp2(ind)=75.0
                       
            
            Map = ResultW  + Pic


;;; Now need to map the result onto a color table

            ma = 0.0004
            mi = 2.5e-08

;; now do clipping in Map

            ind=where(Map lt mi) 
            if ind(0) ne -1 then Map(ind)=mi
            ind=where(Map gt ma)
            if ind(0) ne -1 then Map(ind)=ma

;; now map the log of the field linearly onto the color table

            cols= 255
            
            ImageDens= byte((alog10(Map)-alog10(mi))/(alog10(ma/mi)) * (cols-2) +2)
            

            imageTemp = fltarr(XPixels, ypixels)
;Imagetemp(*) = (0.9-gasfrac)/0.9 * 180 
            Imagetemp(*) = (mean_age)/0.2 * 180 
            
            Imagetemp += imageTemp2
            
            
            imageTemp = fix(ImageTemp)
            
            
            cols=256
            ColPlane=fltarr(cols,cols,3)
            
            openr,1,"colplane2.dat"
            readu,1,ColPlane
            close,1
            
            JpegPic=bytarr(XPixels,YPixels,3)
            
            for i=0, XPixels-1 do begin
                for j=0, YPixels-1 do begin
                    JPegPic(i,j,0)= ColPlane(imageDens(i,j),imageTemp(i,j),0)
                    JPegPic(i,j,1)= ColPlane(imageDens(i,j),imageTemp(i,j),1)
                    JPegPic(i,j,2)= ColPlane(imageDens(i,j),imageTemp(i,j),2)
                endfor
            endfor
            
;window, xsize = xpixels, ysize=ypixels
;tv, JPegPic, true=3,0
            
;fout= "pic.jpg"
;write_jpeg, fout , Jpegpic, true=3, quality=90


            set_plot, 'z'            
;            window, xsize=Xpixels, ysize=YPixels
            device, set_resolution=[Xpixels, YPixels]
;            loadct,1   
            plot, [xmin,xmax], [ymin,ymax], /nodata, ticklen=-1, position=[0,0,1,1], $
              xstyle=1, ystyle=1
            
            xyouts, PIX*0.1, PIX*0.85, '.', /DEVICE, CHARSIZE = 1, CHARTHICK =1, col=10, font=1 
            InfoPic=TVRD(true=3) 
            InfoPic[0:5,*,*] = 0
            InfoPic[XPixels-6:XPixels-1,*] = 0
            InfoPic[*,0:5,*] = 0
            InfoPic[0,YPixels-6:YPixels-1,*] = 0
            wi = where(InfoPic gt 0)
            JPegPic[wi] = InfoPic[wi]
            write_jpeg, fjpg, JPegPic, true=3, quality=95
            
            device, /close
            
            
            set_plot,'ps'  
            device,filename=feps,/color,bits=8,/encapsulated  
            device,xsize=20.0,ysize=20.0
            
            tv, JPegPic, true=3
            
;            xyouts, PIX*1.0, PIX*20.0, zstring, /DEVICE, CHARSIZE = 4, CHARTHICK =5, col=200, font=1 
            xyouts, [xmin*0.85], [ymax*0.75], zstring,  CHARSIZE = 4, CHARTHICK =5, col=250, font=1

;            oploterror, [0], [-28], [10], [0], psym=sym(3), symsize=0.1, thick=12, $
;              col=200, errcol=200
;            xyouts, [-10], [-25], '20 h!U-1!Nkpc', charsize=3, charthick=6, col=250
;            xyouts, [-4], [-34], '3!U"!N.6', charsize=3, charthick=6, col=250
            
            device,/close   
            set_plot,'z'
            
        
        endfor

    endfor

endfor

end
        
        
        
        
        
 










 
