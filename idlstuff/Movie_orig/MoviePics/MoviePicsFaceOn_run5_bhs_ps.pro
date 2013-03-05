;call device, DECOMPOSED=0, /BYPASS_TRANSLATION, RETAIN=2 before
;openning any drawing window
;; .rn load_center_and_BHs.pro first

hubble=0.7
omega_m=0.3
H_0 = 70                        ;km/s/mpc
t_0=13.7                        ;gyr
;tz14=0.26  ;;z0=15
tz14=0.275  ;;z=14.5
dt=0.002
tbin=0.002

Object_File='../C-Routines_for_IDL/Compute2dProjectionHsmlGiven/slicer.so'

halo=['md15J10E2G4','md15J10E2G4_2','md15J10E2G8','md15J10E2G8M2',$
      'md15J15E2','md10J10E2','md5J5E2', '1gpc_zoom1_1M_2/md15J10E2G4', $
      '1gpc_zoom1_1M_3/md15J10E2G4', '1gpc_zoom1_1M_4/md15J10E2G4']

halost=['t2_1M_1','t1_1M_1','md15J10G81','md15J10G81M2',$
        'md15J15G41','md10J10G42','md5J5G42', 't2_1M_2', $
        't2_1M_3', 't2_1M_4' ]

din='/raid3/yxli/Runs/merger_tree/1gpc_zoom1_hires/ICs/'
dout='/raid3/yxli/Runs/merger_tree/1gpc_zoom1_hires/movies/'

nh=n_elements(halo)

;for i=0, nh-1 do begin
;    cd, dout
;    com=string("mkdir "+ halost[i])
;    spawn, com
;endfor

;cd, "/home/yxli/tools/movies/volker/MoviePics/"

for ih=8, 8 do begin

    label=halost[ih]

    ddir='/raid3/yxli/Runs/merger_tree/1gpc_zoom1_hires/'+halo[ih]+'/'
    if ih ge 7 then ddir='/raid3/yxli/Runs/merger_tree/'+halo[ih]+'/'

    PicsDir = dout+halost[ih]+'/'
         
    fcom=din+'combine_zoom1_0.5rd_1r200_t0_ns1.dat'
    if ih eq 1 then fcom=din+'combine_zoom1_0.5rd_1r200_t0_ns2.dat'

    readcol, fcom, gal, tb, tm, ns, $
      format='(x, a, x,x,x,x,x,x,x,x,x,x, f, f, a)', comment='#'
 
    ng=n_elements(gal)
    ns=fix(ns)

    if ih eq 1 then bhlist=[24804, 45475, 84662, 121907, 235560, 426782, 792454, 1164903] 

;    if ih eq 8 then bhlist=[24804, 45475, 84559, 122001, 235246, 427309, 792300, 1164399] 
    if ih eq 8 then bhlist=[24804, 45475, 84559, 122001,  235246, 792300, 1164399]

    if ih eq 9 then bhlist=[24804, 45475, 84559, 121918, 235197, 428143, 792345, 1164475]


    for ig=0, ng-1 do begin
;    for ig=4, 4 do begin

        Frun=ddir+gal[ig]+'/'
        nsnap=fix((tm[ig]-tb[ig])/dt)
        print, 'nsnap, ns[ig]: ', nsnap, ns[ig]

        Snap = gal[ig]

;spawn, "mkdir "+PicsDir

        get_time_list, Frun, Snap, NumList, TiList
        
        snap_min = 0
        snap_max = nsnap

        if snap_max gt 300 then snap_max=300
        if ih eq 8 and ig eq ng-1 then snap_max=147
        if ih eq 9 and ig eq ng-1 then snap_max=50
     
;NumPix = 30
        NumPix=snap_max
        
        PIX=   768 
;;        PIX=   1024
        
        LenBase =  40.0
        
        Len =  LenBase
        
        n0=0
        if ig eq 0 then n0=1

        for picnr=n0, NumPix-1 do begin
;        for picnr=16, 16 do begin
            
            if ig eq 0 then begin
                ord=picnr
            endif else begin    
                ord=fix(total(ns[0:ig-1]))+picnr
            endelse
            
            exts='000'
            exts=exts+strcompress(string(ord),/remove_all)
            exts=strmid(exts,strlen(exts)-3,3)
            f=PicsDir +"pic_gas_"+exts+".jpg"
            f=strcompress(f,/remove_all)
            fps=PicsDir +"pic_gas_"+exts+".eps"
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
            
            load_center_and_BHs, fname, bhlist, tz, mcenter, g1center, nbh, idbh, mbh, posbh
            print, 'cm, center1, nbh : ', mcenter, g1center, nbh

            wb=where(mbh eq max(mbh), cb)
            cbh=fltarr(3)
            cbh[0]=posbh[0,wb[0]]
            cbh[1]=posbh[1,wb[0]]
            cbh[2]=posbh[2,wb[0]]

            center=fltarr(3)
;            center[0]=cbh[0]
;            center[1]=cbh[1]
;            center[2]=cbh[2]
            center[0]=g1center[0]
            center[1]=g1center[1]
            center[2]=g1center[2]


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
            
;            ind =where(Hsml gt 20.0)
            ind =where(Hsml gt 15.0)
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
            xmin= float(-len) ; left edge of map
            xmax= float( len) ; right edge of map 
            ymin= float(-len)   ; lower
            ymax= float( len)   ; upper edge
            
;            xpixels= long(PIX)*4/3 ; pixels in x-direction
            xpixels= long(PIX) ; pixels in x-direction
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

            set_plot, 'z'            
;            window, xsize=Xpixels, ysize=YPixels
            device, set_resolution=[Xpixels, YPixels]
;            loadct,1   
            plot, [xmin,xmax], [ymin,ymax], /nodata, ticklen=-1, position=[0,0,1,1], $
              xstyle=1, ystyle=1

;            tv, JPegPic, true=3, 0
;            xyouts, PIX*0.1, PIX*0.85, zstring, /DEVICE, CHARSIZE = 3, CHARTHICK =3, col=200, font=1 
            xyouts, PIX*0.1, PIX*0.85, '.', /DEVICE, CHARSIZE = 1, CHARTHICK =1, col=10, font=1
;            xyouts, PIX*0.1, PIX*0.1, exts, /DEVICE, CHARSIZE = 3, CHARTHICK =3, col=200, font=1
;            xyouts, PIX*0.1, PIX*0.1, name, /DEVICE, CHARSIZE = 3, CHARTHICK =2, col=100, font=1 

;            if nbh gt 0 then begin
;                for i=0, nbh-1 do begin
;                    oplot, [bx[i]], [by[i]], psym=sym(1), symsize=0.8, color=20
;                endfor
;            endif
            
;            wset, 0
            InfoPic=TVRD(true=3) 
            InfoPic[0:5,*,*] = 0
            InfoPic[XPixels-6:XPixels-1,*] = 0
            InfoPic[*,0:5,*] = 0
            InfoPic[0,YPixels-6:YPixels-1,*] = 0
            wi = where(InfoPic gt 0)
            JPegPic[wi] = InfoPic[wi]
;            write_jpeg, picfname, JPegPic, true=3, quality=95
            
            device, /close


            set_plot,'ps'  
            device,filename=picfname_eps,/color,bits=8,/encapsulated  
            device,xsize=20.0,ysize=20.0

            tv, JPegPic, true=3
 
;            oploterror, [0], [-28], [10], [0], psym=sym(3), symsize=0.1, thick=12, $
;              col=250, errcol=250
;            xyouts, [-11], [-25], '20 h!U-1!Nkpc', charsize=3, charthick=8, col=250
;            xyouts, [-4], [-34], '3!U"!N.5', charsize=3, charthick=8, col=250

            nn=n_elements(bhlist)       
            for i=0, nn-1 do begin
                w=where(idbh eq bhlist[i], c)
                if c gt 0 then begin
                    oplot, [bx[w[0]]], [by[w[0]]], psym=sym(1), symsize=0.8, color=10
                endif
            endfor

;            if nbh gt 0 then begin
;                for i=0, nbh-1 do begin
;                    oplot, [bx[i]], [by[i]], psym=sym(1), symsize=0.8, color=10
;                endfor
;            endif

;            oplot, [3.0], [5.8], psym=sym(1), symsize=0.8, color=10
;            oplot, [0.0], [9.0], psym=sym(1), symsize=0.8, color=10

            xyouts, [xmin*0.85], [ymax*0.75], zstring,  CHARSIZE = 4, CHARTHICK =5, col=250, font=1

            device,/close   
            set_plot,'z'

skipit:
            
        endfor

    endfor
    
endfor


end





 










 
