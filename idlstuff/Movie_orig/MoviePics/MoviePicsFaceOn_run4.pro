;call device, DECOMPOSED=0, /BYPASS_TRANSLATION, RETAIN=2 before
;openning any drawing window
;; .rn load_center_and_BHs.pro first

hubble=0.7
omega_m=0.3
H_0 = 70                        ;km/s/mpc
t_0=13.7                        ;gyr
tz14=0.275  ;;z0=14.5
dt=0.002
tbin=0.002

Object_File='C-Routines_for_IDL/Compute2dProjectionHsmlGiven/slicer.so'

halo=['test', 'test_2', 'md15_eos0.25_nowind/', 'md15_eos0.5_nowind/', $
      '/mpich2/md15_eos0.5_nowind/', 'pgdt_0306_mpi2/md15_eos0.5_nowind/' ]

halost=['md15_eos0.5_g4','md15_eos0.5_G4', 'md15_eos0.25_GM8', 'md15_eos0.5_g8', $
        'md15_eos0.5_gM8', 'md15_eos0.5_GM8' ]

din='/raid5/yxli/Runs/merger_tree/run4/ICs/1gpc_zoom1/'
dout='/raid5/yxli/Runs/merger_tree/run4/Sims/1gpc_zoom1/maps_fancy/'  

nh=n_elements(halo)

;for i=0, nh-1 do begin
;    cd, dout
;    com=string("mkdir "+ halost[i])
;    spawn, com
;endfor

;cd, "/home/yxli/tools/movies/volker/MoviePics/"

for ih=1, 1 do begin

    label=halost[ih]

    ddir='/raid5/yxli/Runs/merger_tree/run4/Sims/1gpc_zoom1/'+halo[ih]+'/'

    PicsDir = dout+halost[ih]+'/'
         
    fcom=din+'combine_zoom1_0.5rd_1r200_t0_md15.dat'
    if ih eq 1 then fcom=din+'combine_zoom1_0.5rd_1r200_t0_md15_2.dat'
    readcol, fcom, gal, tb, tm, ns, $
      format='(x, a, x,x,x,x,x,x,x,x,x,x, f,f, a)', comment='#'

    ng=n_elements(gal)
    ns=fix(ns)

    if ih eq 0 then bhlist=[24771, 48747, 91920, 119355, 163039, 200537, 271776, 343123]
    if ih eq 1 then bhlist=[24771, 48747, 87321, 116734, 162614, 200198, 271103, 342757]
    if ih eq 3 then bhlist=[24771, 48747, 91920, 119355, 163039, 200537, 271776, 343123]

;    for ig=0, ng-1 do begin
    for ig=ng-1, ng-1 do begin

        Frun=ddir+gal[ig]+'/'
        nsnap=fix((tm[ig]-tb[ig])/dt)
        print, 'nsnap, ns[ig]: ', nsnap, ns[ig]

        Snap = gal[ig]

;spawn, "mkdir "+PicsDir

        get_time_list, Frun, Snap, NumList, TiList
        
        snap_min = 0
        snap_max = nsnap
        if ih eq 0 and ig eq ng-1 then snap_max=121
        if ih eq 1 and ig eq ng-1 then snap_max=26
        if ih eq 3 and ig eq ng-1 then snap_max=91
        
;NumPix = 30
        NumPix=snap_max
        
        PIX=   768 
        
        LenBase =  40.0
        
        Len =  LenBase
        
;        for picnr=0, NumPix-1 do begin
        for picnr=9, 9 do begin
            
            print, 'picnr', picnr

            if ig eq 0 then begin
                ord=picnr
            endif else begin    
                ord=fix(total(ns[0:ig-1]))+picnr
            endelse
            
            exts='000'
            exts=exts+strcompress(string(ord),/remove_all)
            exts=strmid(exts,strlen(exts)-3,3)
            f=PicsDir +"pic_gas_"+halost[ih]+"_"+exts+".jpg"
            f=strcompress(f,/remove_all)
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
            
;            Time = TiList(snap_min) + picnr*(TiList(snap_max)-TiList(snap_min))/float(NumPix-1)
            
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

            if ih eq 1 then tz=tz/hubble
            zz=(2.0/3.0*t_0/(tz+tz14)/sqrt(omega_m))^(2.0/3.0)-1.0

            zstring='z ='+ string(format='(f6.2)', zz)
              
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
            
            
            
            quantity = fltarr(N) ; generate space for mass array
            quantity(*)= Temp(*) ; fill in temperature
            
            xmin= float(-len*4.0/3) ; left edge of map
            xmax= float( len*4.0/3) ; right edge of map 
            ymin= float(-len)   ; lower
            ymax= float( len)   ; upper edge
            
            xpixels= long(PIX)*4/3 ; pixels in x-direction
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
            
            
;; Load two-dimensional color table    

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
            
;;
;            window, xsize=PIX*4/3, ysize=PIX
;            tv, JPegPic, true=3,0
  
            write_jpeg, picfname, JPegPic, true=3, quality=95
            
;; now make postscript plot

            set_plot,'ps'  
            device,filename=picfname_eps,/color,bits=8,/encapsulated  
            device,xsize=20.0*4./3.0,ysize=20.0

            tv, JPegPic, true=3
 
;            loadct,40   

;            !p.position=[0.0, 0.0, 0.0+20*4.0/3.0, 0.0+20.0]

;            if nbh gt 0 then begin
;                for i=0, nbh-1 do begin
;                    oplot, [posbh[0,i]], [posbh[1,i]], psym=sym(1), symsize=0.8, color=35*i
;                endfor
;            endif

            xyouts,0.1, 0.85, zstring, font=1, charsize=4, color=200   

            device,/close   
            set_plot,'x'

skipit:
            
        endfor

    endfor
    
endfor


end





 










 
