;call device, DECOMPOSED=0, /BYPASS_TRANSLATION, RETAIN=2 before
;openning any drawing window
;; .rn load_center_and_BHs.pro first

hubble=0.7
omega_m=0.3
H_0 = 70                        ;km/s/mpc
t_0=13.7                        ;gyr
;tz14=0.26  ;;z0=15
tz14=0.279  ;;z=14.5
dt=0.002
tbin=0.002

Object_File='/home/yxli/tools/movies/volker/C-Routines_for_IDL/Compute2dProjectionHsmlGiven/slicer.so'

halo=['test_2','test', 'md15_eos0.25_nowind/', 'md15_eos0.5_nowind/', $
      '/mpich2/md15_eos0.5_nowind/', 'pgdt_0306_mpi2/md15_eos0.5_nowind/' ]

halost=['md15_eos0.5_G4','md15_eos0.5_g4', 'md15_eos0.25_GM8', 'md15_eos0.5_g8', $
        'md15_eos0.5_gM8', 'md15_eos0.5_GM8' ]

din='/raid5/yxli/Runs/merger_tree/run4/ICs/1gpc_zoom1/'
dout='/raid5/yxli/Runs/merger_tree/run4/Sims/1gpc_zoom1/maps_fancy/test/'  

nh=n_elements(halo)

;for i=0, nh-1 do begin
;    cd, dout
;    com=string("mkdir "+ halost[i])
;    spawn, com
;endfor

;cd, "/home/yxli/tools/movies/volker/MoviePics/"

for ih=0, 0 do begin

    label=halost[ih]

    ddir='/raid5/yxli/Runs/merger_tree/run4/Sims/1gpc_zoom1/'+halo[ih]+'/'

;    PicsDir = dout+halost[ih]+'/'
    PicsDir = dout
         
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
    for ig=6, 6 do begin

        Frun=ddir+gal[ig]+'/'
        nsnap=fix((tm[ig]-tb[ig])/dt)
        print, 'nsnap, ns[ig]: ', nsnap, ns[ig]

        Snap = gal[ig]

;spawn, "mkdir "+PicsDir

        get_time_list, Frun, Snap, NumList, TiList
        
        snap_min = 0
        snap_max = nsnap
        if ih eq 0 and ig eq ng-1 then snap_max=121
        if ih eq 1 and ig eq ng-1 then snap_max=25
        if ih eq 3 and ig eq ng-1 then snap_max=132
        
;NumPix = 30
        NumPix=snap_max
        
        PIX=   768 
        
        LenBase =  40.0
        
        Len =  LenBase
        
        n0=0
;        if ig eq 0 then n0=1

;        for picnr=n0, NumPix-1 do begin
        for picnr=9, 9 do begin
            
            if ig eq 0 then begin
                ord=picnr
            endif else begin    
                ord=fix(total(ns[0:ig-1]))+picnr
            endelse
            
            exts='000'
            exts=exts+strcompress(string(ord),/remove_all)
            exts=strmid(exts,strlen(exts)-3,3)
;            f=PicsDir +"jpeg/pic_gas_"+halost[ih]+"_"+exts+".jpg"
            f=PicsDir +"pic_gas_"+halost[ih]+"_"+exts+".jpg"
            f=strcompress(f,/remove_all)
;            fps=PicsDir +"eps/pic_gas_"+halost[ih]+"_"+exts+".eps"
            fps=PicsDir +"pic_gas_"+halost[ih]+"_"+exts+".eps"
            fps=strcompress(fps,/remove_all)
            picfname= f   
            picfname_eps= fps 

;            openr,1,picfname, error=err ;, /swap_endian
;            if err eq 0 then begin
;                close,1
;                goto,skipit
;            endif
            
            openw,1,picfname,/f77_unformatted
            writeu,1,Len
            close,1
            
            print
            print, "gal[ig], PicNr, ord = ", gal[ig], Picnr, ord
            print
            
            Phi = 0.0
            
            if NumPix gt 1 then begin
                Time = TiList(snap_min) + picnr*(TiList(snap_max)-TiList(snap_min))/float(NumPix-1)
            endif else begin
                Time = TiList(snap_min) + picnr*(TiList(snap_max)-TiList(snap_min))
            endelse
            
            print,'time = ',time
            print,'frun = ',frun
            print,'snap = ',snap
            interpolate_gas_data, Time, Frun, Snap, NumList, TiList, N, Pos, Vel, Hsml, Mass, Temp
            
;;;;;load center;;;;
            ext='000'
            ext=ext+strcompress(string(picnr),/remove_all)
            ext=strmid(ext,strlen(ext)-3,3)
            name=Snap+ "_"+ext
            fname= Frun + name
            fname=strcompress(fname,/remove_all)
            
            load_center_and_BHs, fname, bhlist, tz, center, nbh, idbh, mbh, posbh
            print, 'center, nbh : ', center, nbh

            tz=tz/hubble
            zz=(2.0/3.0*t_0/(tz+tz14)/sqrt(omega_m))^(2.0/3.0)-1.0

            zstring='z='+ string(format='(f6.2)', zz)
              
;;;;;;;;;;;;; start  X-Y- projection of density

            N=   long(N)
            
            Pos2d=  fltarr(2,N) ; generate coordinate array for particles
            
            
            phi=!pi/2.0
            
            Pos(0,*)=Pos(0,*)-center[0]
            Pos(1,*)=Pos(1,*)-center[1]
            Pos(2,*)=Pos(2,*)-center[2]
            
            x = Pos(0,*)  
            z = Pos(1,*)*cos(phi) - Pos(2,*)*sin(phi)
            y = Pos(1,*)*sin(phi) + Pos(2,*)*cos(phi)
           
           bw = 40.0
           
            x = bw/(bw - z)*x
            y = bw/(bw - z)*y
            hsml = bw/(bw - z)*hsml
            Mass = (bw/(bw - z))*Mass

            Pos2d(0,*)= x       ; fill in the x-coordinates
            Pos2d(1,*)= y       ; ...         y-coordinates
            
            Hsml= float(Hsml)
            
            Mass = float(Mass)
            
            ind =where(Hsml gt 20.0)
            if ind(0) ne -1 then begin
                Mass(ind) = 0
                Hsml(ind) = 0.01
            endif
            
            ind =where(z ge bw)
            if ind(0) ne -1 then begin
                Mass(ind) = 0
                Hsml(ind) = 0.001
            endif
            

;;;rescale BH positions
            posbh(0,*)=posbh(0,*)-center[0]
            posbh(1,*)=posbh(1,*)-center[1]
            posbh(2,*)=posbh(2,*)-center[2]
            
            bx = posbh(0,*)  
            bz = posbh(1,*)*cos(phi) - posbh(2,*)*sin(phi)
            by = posbh(1,*)*sin(phi) + posbh(2,*)*cos(phi)            
           
            bx = bw/(bw - bz)*bx
            by = bw/(bw - bz)*by
      
;;;;;;;;            
            
            quantity = fltarr(N) ; generate space for mass array
            quantity(*)= Temp(*) ; fill in temperature
            
;            xmin= float(-len*4.0/3) ; left edge of map
;            xmax= float( len*4.0/3) ; right edge of map 
            xmin= float(-len)   ; left edge of map
            xmax= float( len)   ; right edge of map 
            ymin= float(-len)   ; lower
            ymax= float( len)   ; upper edge
            
;            xpixels= long(PIX)*4/3 ; pixels in x-direction
            xpixels= long(PIX)
            ypixels= long(PIX)  ; pixels in y-direction
            
            ValueXY= fltarr(ypixels,xpixels) ; holds the result
            ValueXY_Temp= fltarr(ypixels,xpixels) ; holds the result
            
            S = CALL_EXTERNAL(Object_File, $
                              'slice', $
                              N, $
                              pos2D, hsml, mass, quantity, $
                              xmin,xmax,ymin,ymax,$
                              Xpixels,Ypixels,$
                              ValueXY, ValueXY_Temp)
            
            ValueXY = transpose(ValueXY)
            ValueXY_Temp= transpose(ValueXY_Temp)
            
            
            
;;; Now need to map the result onto a color table
            
            Map= ValueXY
            
;;    ma = max(Map)
            ma = 10*3.0e-5 / (16.0/len)^2  
            mi=  ma /8000
            
            print, "max= ", ma, "min=", mi, max(map)
            
            
            
;; now do clipping in Map

            ind=where(Map lt mi) 
            if ind(0) ne -1 then Map(ind)=mi
            ind=where(Map gt ma)
            if ind(0) ne -1 then Map(ind)=ma
            
;; now map the log of the field linearly onto the color table

            cols= 255
            
            ImageDens= byte((alog10(Map)-alog10(mi))/(alog10(ma/mi)) * (cols-2) +2)
            
    
            ind=where(ValueXY eq 0)
            if ind(0) ne -1 then ValueXy(ind)=1
            TempMap= ValueXY_Temp / ValueXY ; define mass-weighted temperature map

           
            Map= TempMap
            
            ma= 3.0e6         ; select maximum value (brightest color)
            mi= 8.0e3           ; select minimum value


;; now do clipping in Map

            ind=where(Map lt mi) 
            if ind(0) ne -1 then Map(ind)=mi
            ind=where(Map gt ma)
            if ind(0) ne -1 then Map(ind)=ma
            
;; now map the log of the field linearly onto the color table

            cols= 255
            
            ImageTemp= byte((alog10(Map)-alog10(mi))/(alog10(ma/mi)) * (cols-2) +2)
 
;;now call pic_glare
            pic_glare, Frun=frun, Snap=Snap, Num=picnr, PIX=pix, Len=len, $
              bhlist=bhlist, image_glare=image_glare, image_den=image_den
           
;;now combine these two

;            ImageTemp2 = ImageTemp+image_glare
            
;            print, 'Com. max(ImageTemp)', max(ImageTemp)

;            ImageTemp2 = ImageTemp*.45 + image_glare*1.85
;            ImageTemp2 = image_glare*3.3
;            ImageTemp2 = fix(ImageTemp2)

            ImageTemp2 = image_glare
            imageDens2 = image_den

;; Load two-dimensional color table  

            cols=256
            ColPlane=fltarr(cols,cols,3)
            
            openr,1,"colplane2.dat"
            readu,1,ColPlane
            close,1
            
            JpegPic=bytarr(XPixels,YPixels,3)
            
            for i=0, XPixels-1 do begin
                for j=0, YPixels-1 do begin
                    JPegPic(i,j,0)= ColPlane(imageDens2(i,j),imageTemp2(i,j),0)
                    JPegPic(i,j,1)= ColPlane(imageDens2(i,j),imageTemp2(i,j),1)
                    JPegPic(i,j,2)= ColPlane(imageDens2(i,j),imageTemp2(i,j),2)
                endfor
            endfor

            set_plot, 'z'            
;            window, xsize=Xpixels, ysize=YPixels
            device, set_resolution=[Xpixels, YPixels]
;            loadct,1   
            plot, [xmin,xmax], [ymin,ymax], /nodata, ticklen=-1, position=[0,0,1,1], $
              xstyle=1, ystyle=1

;            tv, JPegPic, true=3, 0
            xyouts, PIX*0.1, PIX*0.85, zstring, /DEVICE, CHARSIZE = 3, CHARTHICK =3, col=200, font=1 
            xyouts, PIX*0.1, PIX*0.1, exts, /DEVICE, CHARSIZE = 3, CHARTHICK =3, col=200, font=1
;            xyouts, PIX*0.1, PIX*0.1, name, /DEVICE, CHARSIZE = 3, CHARTHICK =2, col=100, font=1 

            if nbh gt 0 then begin
                for i=0, nbh-1 do begin
                    oplot, [bx[i]], [by[i]], psym=sym(1), symsize=0.8, color=20
                endfor
            endif
            
;            wset, 0
            InfoPic=TVRD(true=3) 
            InfoPic[0:5,*,*] = 0
            InfoPic[XPixels-6:XPixels-1,*] = 0
            InfoPic[*,0:5,*] = 0
            InfoPic[0,YPixels-6:YPixels-1,*] = 0
            wi = where(InfoPic gt 0)
            JPegPic[wi] = InfoPic[wi]
            write_jpeg, picfname, JPegPic, true=3, quality=95
            
            device, /close


            set_plot,'ps'  
            device,filename=picfname_eps,/color,bits=8,/encapsulated  
            device,xsize=20.0*4./3.0,ysize=20.0

            tv, JPegPic, true=3
 
;            loadct,40   

;            !p.position=[0.0, 0.0, 0.0+20*4.0/3.0, 0.0+20.0]

;            if nbh gt 0 then begin
;                for i=0, nbh-1 do begin
;                    oplot, [bx[i]], [by[i]], psym=sym(1), symsize=0.8, color=10
;                endfor
;            endif

;            xyouts,0.1, 0.85, zstring, font=1, charsize=4, color=200   

            device,/close   
            set_plot,'z'

skipit:
            
        endfor

    endfor
    
endfor


end





 










 
