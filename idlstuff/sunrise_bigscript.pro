;
; From chayward@cfa.harvard.edu Mon Oct 29 21:46:39 2007
; Date: Sun, 28 Oct 2007 18:22:47 -0400
; From: Chris Hayward <chayward@cfa.harvard.edu>
; To: Thomas J. Cox <tcox@cfa.harvard.edu>
; Subject: script
; 
; Did you get the script?
; 
; If not, here it is again. I've changed the vulgar names to non-vulgar
; names. Maybe that will do the trick.
; 
;
;
;
;
;

;-------------------------------------------------------------

 ; plot one line for each camera of the data
pro cam_plot,x,y,n_cam, c,psym=psym,thick=thick,symsize=symsize
k = fltarr (n_elements (y))
k = y
for i = 0, n_cam- 1 do begin
    oplot, x(i:*: n_cam), k (i:*: n_cam), color = c,psym=psym, $
      thick=thick,symsize=symsize
end
end

;-------------------------------------------------------------

pro save_plot,filename
if !D.name eq 'Z' then begin
    print, "Saving image " + filename
    a=tvrd()
    tvlct,r,g,b,/get
    write_png,filename,a,r,g,b
end else if !D.name eq 'PS' then begin
    ;; if its PostScript we just need to close
    print, "Closing PostScript file"
    device,/close
end else $
  junk=get_kbrd(1)
end

;-------------------------------------------------------

pro make_2dcolor,x,y,nbinx,nbiny,weight=w
; x,y contains points to be made into a color density plot
datax=!X.CRANGE
datay=!Y.CRANGE
print,datax,datay

; check if axes are flipped, in that case we have to fix because
; hist_2d will choke
if datax[1] lt datax[0] then begin
    sx=-datax
    xp=-x
end else begin
    sx=datax
    xp=x
end
if datay[1] lt datay[0] then begin
    sy=-datay
    yp=-y
end else begin
    sy=datay
    yp=y
end


;d=hist_2d(x,y,min1=sx[0],max1=sx[1],bin1=(sx[1]-sx[0])/(nbinx-0.001),$
;          min2=sy[0],max2=sy[1],bin2=(sy [1] - sy [0])/(nbiny-0.001))
d=sunrise_hist2d(xp,yp,w,min1=sx[0],max1=sx[1],binsize1=(sx[1]-sx[0])/(nbinx-0.001),$
         min2=sy[0],max2=sy[1],binsize2=(sy [1] - sy [0])/(nbiny-0.001))

print,'Histo values ',min(d),max(d)
if min(d) eq max(d) then stop
sz=size(d)
if (sz[1] ne nbinx) or (sz[2] ne nbiny) then $
  print,'The histo array has incorrect size!'

if !D.name eq 'Z' or !D.name eq 'X' then begin
    ;; get device coordinates of plotting region
    devmin = convert_coord([datax[0],datay[0]],/data,/to_dev)
    devmax = convert_coord([datax[1],datay[1]],/data,/to_dev)
    print,'Rebinning data to ',devmin,devmax

    ;; rebin data to req # of pixels
    h=congrid(d,devmax[0]-devmin[0]+1,devmax[1]-devmin[1]+1)
    ;; display it
    bh=bytscl(max(h)-h,top=(!D.TABLE_SIZE-13))+12
    tv,bh,devmin[0],devmin[1]
end else if !D.name eq 'PS' then begin
    ;; if it's a PostScript plot we can't rebin data because the
    ;; number of pixels is huge, but we can scale the output
    bh=bytscl(max(d)-d,top=(!D.TABLE_SIZE-13))+12
;stop
                                ;tv,bh,datax[0],datay[0],xsize=sx[1]-sx[0],ysize=sy[1]-sy[0],$
    tv,bh,datax[0],datay[0],xsize=datax[1]-datax[0],ysize=datay[1]-datay[0],$
      /data
end
end
;stop

;
;-------------------------------------------------------

pro makehist,v,w,nbins,hx,hy
; get current plot range
datax=!X.CRANGE
print,datax

; check if axes are flipped, in that case we have to fix because
; hist_2d will choke
if datax[1] lt datax[0] then begin
    sx=-datax
    v=-v
end else $
  sx=datax

hy=hist1d(v,w,min=sx[0],max=sx[1],binsize=(sx[1]-sx[0])/(nbins-0.0001))
hy=[0,hy,0]
hx= (indgen(n_elements (hy)) - 0.5)/(n_elements (hy)-2)* $
  (sx[1]-sx[0]) + sx [0]
; stop
if datax[1] lt datax[0] then $
  hx=-hx
end



;-------------------------------------------------------

pro plot_stuff,index_file,PS=ps,plot_only=plot_only,bw=bw,notitle=notitle

; we use Z-buffer device to save plots, VASTLY faster than reading the
; X Windows pixmap
if keyword_set (plot_only) then set_plot, 'x' else set_plot,'z'

if (keyword_set (ps) or $
    (keyword_set (default_plot) and (!D.name eq 'PS'))) then begin
    print,"Plotting to postscript..."
    bgc = 1
    fgc = 0
    c1 = 2
    c2 = 11
    if keyword_set (ps) then begin
        set_plot,'ps'
        !P.FONT=-1
        device,filename='seds.ps',/encaps,/color,bits_per_pixel=8, $
          preview=0,pre_depth=8,pre_x=3.5,pre_y=3,/inches
        device, SET_CHARACTER_SIZE=[300,350],xsize=14,ysize=12
        if keyword_set (landscape) then device,/landscape
        !X.THICK=3
        !Y.THICK=3
        !P.THICK=2
        !X.MARGIN=[8,3]
    end
end else begin
    if !D.name eq 'X' then begin
        !P.FONT=-1
        device,true_color=24,decomposed=0,set_font='helvetica'
        device,set_character_size=[10,12]
    end
    bgc = 'ffffff'x
    fgc = 0
    c1 ='ff5500'x
    c2 ='0055ff'x
    c3 ='00ee00'x
    c4 ='aa00aa'x
    c5='8fff00'x
    c6='00ff8f'x
    c7='ee0000'x
    c8='0000ee'x
    c9='8f00ff'x
    c10='ff008f'x
    c11='000000'x
    c12='555555'x
    c13='ffff00'x
                                ;if not keyword_set(default_plot) then set_plot,'x'
end
bgc=1
c1=2
c2=3
c3=4
c4=5
c5=6
c6=7
cc = [c1,c2,c3,c4,c5,c6]

!P.MULTI=0
; load color map
if keyword_set(bw) then begin
    loadct,0                    ; bw colora map
    c1=fgc
    c2=fgc
    c3=fgc
    c4=fgc
    c5=fgc
    c6=fgc
    cc = [c1,c2,c3,c4,c5,c6]
end else loadct,3               ; red temp

colr=[0,255,255,0  ,0  ,175,143,  0,  0,  0,155,  0]
colg=[0,255,85 ,85 ,215,0  ,255,255,175,130,  0,  0]
colb=[0,255,0  ,255,0  ,175,  0,143,255,186,255,255]
tvlct,colr,colg,colb

; read fits data
if n_elements (index_file) gt 1 then begin
    print,'Loading ',index_file[0]
    d=mrdfits(index_file [0], 1)
    for i = 1, n_elements (index_file)- 1 do begin
        print,'Loading ',index_file[i]
        temp = mrdfits (index_file [i], 1)
        d = [d, temp]
    end
end else $
  d=mrdfits(index_file,1)
n=n_elements(d.snapshot_file)
n_camera=max(d.camera)+1

weight=d.dt/(max(d.camera)+1)

; and now begin plotting
;goto, magstart
;goto,there
;; time-dependent quantities

; luminosity
if !D.name eq 'PS' then device,file='l-t.eps'
yrange=[min(d.l_scat/4e26),max(d.l_bol_grid/4e26) ]
plot,d.snaptime/1e9,d.l_bol_grid/4e26,/ylog,xtitle="time/Gyr", $
  ytitle="L/L!D"+sunsymbol()+"!N", $ 
                                ;title="Luminosity vs. time", $
background = bgc, color = fgc,/nodata, $
  yrange=yrange,xtickinter=0.5
oplot,d.snaptime/1e9,d.l_bol_grid/4e26,color=c1
oplot,d.snaptime/1e9,d.l_bol_absorbed/4e26,color=c2
cam_plot,d.snaptime/1e9,d.l_scat/4e26,n_camera,c3
legend,['Bolometric','Absorbed','Not absorbed'],colors=[c1,c2,c3],textcolor=fgc,linestyle=0,/right,box=0,/bot
; b/w
;oplot,d.snaptime/1e9,d.l_bol_grid/4e26,color=fgc, linestyle=2
;oplot,d.snaptime/1e9,d.l_bol_absorbed/4e26,color=fgc,linestyle=1
;cam_plot,d.snaptime/1e9,d.l_scat/4e26,n_camera,fgc
;legend,['Bolometric','Absorbed','Not absorbed'],colors=fgc,textcolor=fgc,linestyle=[2,1,0],/right,box=0,/bot,charsize=0.8
save_plot,'l-t.png'
;stop

; attenuation
if !D.name eq 'PS' then device,file='att-t.eps'
plot, d.snaptime, d.l_bol_absorbed/d.l_bol_grid, $
  xtitle="time/yr",ytitle="Attenuation", /ylog,$
  title="Attenuation vs. time", background = bgc, color = fgc,/nodata
oplot, d.snaptime, d.l_bol_absorbed/d.l_bol_grid, color = c1
save_plot, 'att-t.png' 
; stop

; SFR
if !D.name eq 'PS' then device,file='sfr-t.eps'
plot, d.snaptime/1e9, d.sfr_tot, $
  xtitle="time/Gyr",ytitle="SFR/(M!D"+sunsymbol()+"!N/yr)", $
                                ;title="SFR vs. time", $
background = bgc, color = fgc,/nodata
oplot, d.snaptime/1e9, d.sfr_tot, color = c1
save_plot,'sfr-t.png'
; stop

; metallicity
if !D.name eq 'PS' then device,file='z-t.eps'
plot,d.snaptime,d.z_gas, $
  xtitle="time/yr",ytitle="Metallicity", $
  title="Metallicity vs. time", background = bgc, color = fgc,/nodata
oplot,d.snaptime,d.z_gas,color = c1
save_plot,'z-t.png'
; stop

; H_a and H_beta
if !D.name eq 'PS' then device,file='l_line-t.eps'
L_H_b=fltarr(n_elements(d.L_H_b))
L_H_b=d.L_H_b
yrange = [min (L_H_b(where(L_H_b gt 0))/4e26), max (d.L_H_a_ns/4e26) ]
plot,d.snaptime, d.L_H_a_ns /4e26,/ylog,xtitle="time/yr",ytitle="L/L!D"+sunsymbol()+"!N", $
  title="Line luminosity vs. time", background = bgc, color = fgc,/nodata, $
  yrange=yrange
cam_plot,d.snaptime,d.L_H_a/4e26,n_camera,c1,psym=2,syms=0.5
;cam_plot,d.snaptime,d.L_H_a/4e26,n_camera,c3
cam_plot,d.snaptime,d.L_H_b/4e26,n_camera,c4,psym=4,syms=0.5
;cam_plot,d.snaptime,d.L_H_b/4e26,n_camera,c4
oplot,d.snaptime,d.L_H_a_ns/4e26,color=c2
oplot,d.snaptime,d.L_H_b_ns/4e26,color=c3
legend,[textoidl('H\alpha'),textoidl('H\beta'),$
        textoidl('H\alpha w/o dust'),textoidl('H\beta w/o dust')],$
  psym=[2,4,-3,-3],/bottom,$
  colors=[c1,c4,c2,c3],textcolor=fgc,linestyle=0,/left,box=0
save_plot,'l_line-t.png'
;stop

; line ratio
if !D.name eq 'PS' then device,file='liner-t.eps'
yrange = [min (d.L_H_a(where(L_H_b gt 0))/L_H_b(where(L_H_b gt 0))), $
          max (d.L_H_a(where(L_H_b gt 0))/L_H_b(where(L_H_b gt 0)))]
plot,d.snaptime,d.L_H_a/d.L_H_b,/ylog,xtitle="time/yr", $
  ytitle=textoidl("L_{H\alpha}/L_{H\beta}"), $
  title="Line ratio vs. time", background = bgc, color = fgc,/nodata, $
  yrange=yrange
cam_plot,d.snaptime,d.L_H_a/d.L_H_b,n_camera,c1,psym=4,syms=0.5
;cam_plot,d.snaptime,d.L_H_a/d.L_H_b,n_camera,c1
oplot,d.snaptime,d.L_H_a_ns/d.L_H_b_ns,color=c2
legend,['w/o dust','w/dust'],box=0,psym=[-3,4],$
  colors=[c2,c1],textcolor=fgc,linestyle=0,/right
save_plot,'liner-t.png'
;stop

; H_a and H_beta equivalent width
if !D.name eq 'PS' then device,file='eqw-t.eps'
L_H_b=d.eqw_H_b
;yrange = [min (L_H_b(where(L_H_b gt 0))), max (d.eqw_H_a_ns) ]
plot,d.snaptime, d.eqw_H_a_ns ,xtitle="time/yr", $
  ytitle=textoidl("equivalent width/A"), $
  title="Line equivalent width vs. time", background = bgc, color = fgc,/nodata
                                ;yrange=yrange
cam_plot,d.snaptime,d.eqw_H_a,n_camera,c1,syms=0.5,psym=2
;cam_plot,d.snaptime,d.eqw_H_a,n_camera,c3
cam_plot,d.snaptime,d.eqw_H_b,n_camera,c4,syms=0.5,psym=4
;cam_plot,d.snaptime,d.eqw_H_b,n_camera,c4
cam_plot,d.snaptime,d.eqw_H_a_ns,n_camera,c2
cam_plot,d.snaptime,d.eqw_H_b_ns,n_camera,c3
legend,[textoidl('H\alpha'),textoidl('H\beta'),$
        textoidl('H\alpha w/o dust'),textoidl('H\beta w/o dust')], $
  /top ,psy=[2,4,-3,-3],$
  colors=[c1,c4,c2,c3],textcolor=fgc,linestyle=0,/right,box=0
save_plot,'eqw-t.png'
;stop

; UV luminosity
if !D.name eq 'PS' then device,file='luv-t.eps'
yrange = [min (d.L_1600/4e26), max (d.L_1600_ns/4e26) ]
plot,d.snaptime, d.L_1600_ns/4e26,/ylog,xtitle="time/yr", $
  ytitle=textoidl("L_{1600}/L_{"+sunsymbol()+"}"), $
  title="UV luminosity vs. time", background = bgc, color = fgc,/nodata, $
  yrange=yrange
oplot,d.snaptime,d.L_1600_ns/4e26,color=c2
cam_plot,d.snaptime,d.L_1600/4e26,n_camera,c1
legend,['w/o dust','w/dust'], $
  colors=[c2,c1],textcolor=fgc,linestyle=0,/right,box=0
save_plot,'luv-t.png'
; stop

; UV slope
if !D.name eq 'PS' then device,file='beta-t.eps'
yrange = [min ([d.betaM95_ns,d.betaH98_ns]), max ([d.betaM95,d.betaH98]) ]
plot,d.snaptime, d.betaM95_ns,xtitle="time/yr", $
  ytitle=textoidl("\beta"), $
  title="UV slope vs. time", background = bgc, color = fgc,/nodata, $
  yrange=yrange
oplot,d.snaptime,d.betaM95_ns,color=c2
cam_plot,d.snaptime,d.betaM95,n_camera,c1
oplot,d.snaptime,d.betaH98_ns,color=c4
cam_plot,d.snaptime,d.betaH98,n_camera,c3
legend,['M95 w/o dust','M95 w/dust','H98 w/o dust','H98 w/dust'], $
  colors=[c2,c1,c4,c3],textcolor=fgc,linestyle=0,/right,box=0
save_plot,'beta-t.png'
; stop

; IRX
if !D.name eq 'PS' then device,file='irx-t.eps'
irxM95=alog10(d.l_fir_iras/d.l_1600)
irxH98=alog10(d.l_fir_iras/d.l_1900)
yrange = [min ([irxM95,irxH98]), max ([irxM95,irxH98])]
plot,d.snaptime, d.betaM95_ns,xtitle="time/yr", $
  ytitle=textoidl("IRX"), $
  title="IRX vs. time", background = bgc, color = fgc,/nodata, $
  yrange=yrange
cam_plot,d.snaptime,irxM95,n_camera,c1
cam_plot,d.snaptime,irxH98,n_camera,c3
legend,['M95','H98'], $
  colors=[c1,c3],textcolor=fgc,linestyle=0,/right,box=0
save_plot,'irx-t.png'
; stop

magstart:
; magnitudes (looks for fields named AB_* (and not _NODUST))
; (nodust columns are immediately following)
f=tag_names(d)
l= 0 ; the field numbers for dusty mags
lns=0; field numbers for _nodust mags
t=''
for i = 0,n_elements(f)-1 do begin
    if (strmid(f[i],0,3) eq 'AB_') and $
      not (strmid(f[i],6,7,/reverse) eq '_NODUST') then begin
        ; find nodust
        ns=where(f eq f[i]+'_NODUST')
        print,f[i],i,ns
        sz=size(l)
        if sz[0] eq 0 then begin
            l = [i]
            lns = [ns]
            t = [f[i]]
        end else  begin
            l=[l,i]
            lns = [lns, ns]
            t = [t, f [i]]
        end

                                ;yrange[1]=min([yrange[1],d.(i)]) 
                                ;yrange[0 ]= max ([yrange[0],d.(i)])
    end
end
; now plot
nn=0
nmags_per_plot= 3

for j=0,n_elements (l) - 1,nmags_per_plot do begin
                                ; find range
    yrange=[-100.,100.]
    for i = j,j+nmags_per_plot-1 do begin ;n_elements(l)-1 do begin
        if i lt n_elements (l) then begin
            if max(d.(lns[i])) gt 100 then $
              yrange[1] = min ([yrange[1],d.(l[i])]) $
            else $
              yrange[1] = min ([yrange[1],d.(lns[i])]) 
            yrange[0] = max ([yrange[0],d.(l[i])])
        end
    end
    ; now plot
    plot_name='mag'+ strcompress(string (nn,format = '( i 2.1)'),/r) + '-t'
    if !D.name eq 'PS' then device,file=plot_name+'.eps'

    if nmags_per_plot gt 1 then begin
        ytitle="AB Magnitude"
        title="Magnitudes vs. time"
        dcolor=0
    end else begin
        ytitle=t[j]+" Magnitude"
        title=t[j]+" Magnitude vs. time"
        dcolor=1
    end
    plot,d.snaptime, d.(l[j]),xtitle="time/yr", $
      ytitle=ytitle, $
      title=title, background = bgc, color = fgc,/nodata, $
      yrange=yrange
    k = - 1
    for i = j,j+nmags_per_plot-1 do begin ;n_elements(l)-1 do begin
        if i lt n_elements (l) then begin
            cam_plot,d.snaptime,d.(l[i]),n_camera,cc[i-j] 
            cam_plot,d.snaptime,d.(lns[i]),n_camera,cc[i-j+dcolor],thick=2
            k = k+ 1
        end
    end
    if nmags_per_plot gt 1 then $
      legend,t[j:j+k], colors=cc[0:nmags_per_plot-1], $
      textcolor=fgc,linestyle=0,/right, $
      /bottom,box=0
    save_plot,plot_name+'.png'
    nn = nn+ 1
                                ;stop
end


; also make dust attenuation graphs of the magnitudes
nn=0
for j=0,n_elements (l) - 1,nmags_per_plot do begin
                                ; find range
    yrange=[-100.,100]
    for i = j,j+nmags_per_plot-1 do begin ;n_elements(l)-1 do begin
        if i lt n_elements (l) then begin
            if max(d.(l[i]+1)) lt 100 then begin
                yrange[0] = max ([yrange[0], d.(l[i]) - d.(lns[i]) ])
                yrange[1] = min ([yrange[1], d.(l[i]) - d.(lns[i]) ])
            end
        end
    end
                                ; now plot to
    plot_name='attn'+ strcompress(string (nn,format = '(i2.1)'),/r) + '-t'
    if !D.name eq 'PS' then device,file=plot_name+'.eps'

    if nmags_per_plot gt 1 then begin
        ytitle="Dust attenuation/magnitudes"
        title="Dust attenuation vs. time"
        dcolor=0
    end else begin
        ytitle="Dust attenuation in "+ t[j]+"/magnitudes"
        title="Dust attenuation in "+t[j]+" vs. time"
        dcolor=1
    end
    plot,d.snaptime, d.(l[j]),xtitle="time/yr", $
      ytitle=ytitle,$
      title=title,$
      background = bgc, color = fgc,/nodata, $
      yrange=yrange
    k = - 1
    for i = j,j+nmags_per_plot-1 do begin ;n_elements(l)-1 do begin
        if i lt n_elements (l) then begin
            cam_plot,d.snaptime,d.(l[i])-d.(lns[i]),n_camera,cc[i-j]
            k = k+ 1
        end
    end
    if nmags_per_plot gt 1 then $
      legend,t[j:j+k], colors=cc[0:nmags_per_plot-1],textcolor=fgc, $
      linestyle=0,/right, /bottom,box=0
    save_plot,plot_name+'.png'
    nn = nn+ 1
                                ;stop
end


; colors
; makeup colors
;ref = where (t eq "AB_R_SDSS")
ref = where (t eq "AB_G_SDSS")
ll=0
tt=''
colors=0
colors_nodust =0
for i = 0, n_elements (t)-1 do begin
    if (i ne ref) then begin
        sz = size (ll)
        if i lt ref then begin
            colorname=t[i]+'-'+t[ref]
            color=d.(l[i]) - d.(l[ref])
            color_nodust =d.(lns[i]) - d.(lns[ref])            
        end else begin
            colorname=t[ref]+'-'+t[i]
            color=d.(l[ref]) - d.(l[i]) 
            color_nodust =d.(lns[ref]) - d.(lns[i]) 
        end
;stop
        if sz[0] eq 0 then begin
            ll = [i]
            tt = [colorname]
            colors=[[color]]
            colors_nodust =[[color_nodust]]
        end else  begin
            ll=[ll,i]
            tt = [tt, colorname]
            colors=[[colors],[color]]
            colors_nodust =[[colors_nodust],[color_nodust]]
        end
    end
end
;stop
; now plot
;goto,there
nn=0
for j=0,n_elements (ll) - 1,nmags_per_plot do begin
                                ; find range
    yrange=[-100.,100.]
    for i = j,j+nmags_per_plot-1 do begin ;n_elements(l)-1 do begin
        if i lt n_elements (ll) then begin
            if (colors_nodust [0,i] lt 100) and $
              (colors_nodust[0,i] gt -100) then $ 
              yrange[1] = min ([yrange[1],colors_nodust [*,i]]) $
            else $
              yrange[1] = min ([yrange[1],colors [*,i]])
            yrange[0] = max ([yrange[0],colors[*,i]])
        end
    end
    plot_name ='col'+ strcompress(string (nn, format = '(i2.1)'),/r)+ '-t'
    if !D.name eq 'PS' then device,file=plot_name+'.eps'

    if nmags_per_plot gt 1 then begin
        ytitle="AB Color/magnitudes"
        title="Color vs. time"
        dcolor=0
    end else begin
        ytitle= tt[j]+" color/magnitudes"
        title=tt[j]+" color vs. time"
        dcolor=1
    end
    plot,d.snaptime, colors[*,0],xtitle="time/yr", $
      ytitle= ytitle, $
      title = title, background = bgc, color = fgc,/nodata, $
      yrange=yrange
    k=-1
    for i = j,j+nmags_per_plot-1 do begin 
        if i lt n_elements(ll) then begin
            k = k+ 1
            cam_plot,d.snaptime,colors[*,i],n_camera,cc[i-j]
            cam_plot,d.snaptime,colors_nodust [*,i],n_camera, $
              cc[i-j+dcolor],thick=2
        end
    end
    if nmags_per_plot gt 1 then $
      legend,tt[j:j+k], colors=cc[0:nmags_per_plot-1],textcolor=fgc,linestyle=0,/right,box=0
    save_plot,plot_name+'.png'
    nn = nn+ 1
                                ;stop
end


;; and now quantity-quantity plots

there:
; infrared excess-UV slope relation

readcol,'~/mhc99data.txt',mhc99beta,mhc99irx,mhc99ffir

irx=alog10(d.l_fir_iras/d.L_1600)
beta=d.betaM95
bns=d.betaM95_ns
a1600=-2.5*alog10(d.l_1600/d.l_1600_ns)

idx=where(d.l_bol_grid/4e26 lt 7e13)
irx=irx[idx]
beta = beta [idx]
bns = bns [idx]
a1600=a1600[idx]

if !D.name eq 'PS' then device,file='irx-beta.eps'
plot,[1],[1],$
;      xtitle="!7b!X(M95)",ytitle="log(L!DFIR,IRAS!N/L!D1600)", $
;      title="IRX!D1600!N-Beta  (Meurer, Heckman & Calzetti 1999)", $
background = bgc, color = fgc,/nodata, $
  xrange=[-2.5,1],yrange=[-1,4],xst=1,yst=1
make_2dcolor,beta,irx,40,40,weight=weight[idx]
plot,[1],[1],$
  xtitle=textoidl('\beta'), $
                                ;xtitle=textoidl('\beta(M95)'), $
ytitle=textoidl("log(L_{FIR,IRAS}/L_{1600})"), $
                                ;title=textoidl("IRX_{1600}-\beta  (Meurer, Heckman & Calzetti 1999)"),$
background = bgc, color = fgc,/nodata, /noerase,$
  xrange=[-2.5,1],yrange=[-1,4],xst=1,yst=1
; the MHC99 fit
bfit= indgen(100)/100.*3.0-2.2
logirxfit=alog10(10^(0.4*(4.43+1.99*bfit))-1)+.076
;oplot, d.betaM95_ns, alog10 (d.l_fir_iras/d.L_1600_ns) ,color=c2,psym=4
;oplot, d.betaM95, alog10 (d.l_fir_iras/d.L_1600) ,color=c1,psym=5
;oplot, d.betaM95, alog10 (d.l_bol_absorbed/d.L_1600) ,color=c3,psym=5
oplot,bfit, logirxfit,color= fgc
oplot,mhc99beta, mhc99irx,color= c2,ps=1
oplot,[-0.66,-0.56, -0.87, -0.09,-0.98], $
  [2.47, 2.71,2.23,3.64, 2.48], ps =4,col=c3
oplot,[-1.35, -0.64], $
  [1.13,  2.48], ps =5,col=c3
legend,['Simulations','MHC fit','MHC','LIRGs','ULIRGs'], /bottom, $
  colors=[180,fgc,c2,c3,c3],textcolor=fgc,linestyle=0,/right,box=0, $
  psy=[-3,-3,1,5,4],thick=[20,1,1,1,1]
; b/w
;oplot,bfit, logirxfit,color= fgc
;oplot,mhc99beta, mhc99irx,color= fgc,ps=1
;oplot,[-0.66,-0.56, -0.87, -0.09,-0.98], $
;      [2.47, 2.71,2.23,3.64, 2.48], ps =4,col=fgc
;oplot,[-1.35, -0.64], $
;      [1.13,  2.48], ps =5,col=fgc
;legend,['Simulations','MHC fit','MHC','(U)LIRGs'], /bottom, $
;       colors=[128,fgc,fgc,fgc],textcolor=fgc,linestyle=0,/right,box=0, $
;       psy=[-3,-3,1,4]
save_plot, 'irx-beta.png'
;stop

; MHC99 A1600 vs beta
if !D.name eq 'PS' then device,file='a1600-beta.eps'
plot,[1],[1],$
;      xtitle="!7b!X(M95)",ytitle="log(L!DFIR,IRAS!N/L!D1600)", $
;      title="IRX!D1600!N-Beta  (Meurer, Heckman & Calzetti 1999)", $
background = bgc, color = fgc,/nodata, $
  xrange=[-.5,6],yrange=[0,7],xst=1,yst=1
make_2dcolor,beta-bns,a1600,40,40,weight=weight[idx]
plot,[1],[1],$
  xtitle=textoidl('\Delta\beta'), $
                                ;xtitle=textoidl('\beta(M95)'), $
ytitle=textoidl("A_{1600}"), $
                                ;title=textoidl("A_{1600}-\beta  (Meurer, Heckman & Calzetti 1999)"),$
background = bgc, color = fgc,/nodata, /noerase,$
  xrange=[-.5,6],yrange=[0,7],xst=1,yst=1
; the MHC99 fit
bfit= indgen(100)/100.*6.0-2.2
Afit= 0+1.99*(bfit-0.5)
oplot,bfit, Afit,color= fgc
Afit= 0+8.9*(bfit+0*1.8-0.3)
oplot,bfit, Afit,color= c2
Afit= 0+3.0*(bfit+0*1.8-0.4)
oplot,bfit, Afit,color= c2,line=1
legend,['Simulations','MHC fit','D03 MW tot.','D03 MW abs.'], /bottom , $
  colors=[128,fgc,c2,c2],textcolor=fgc,linestyle=[0,0,0,1],/right,box=0       
save_plot, 'a1600-beta.png'

;stop

; SFR - line luminosity
; need to make this into a 2x2 with density plot

if !D.name eq 'PS' then device,file='sfr-line.eps'
!P.MULTI=[0,2,2,0,1]
L_H_b=d.L_H_b
L_H_a=d.L_H_a
kuken1=where(L_H_a lt 0)
kuken2=where(L_H_b lt 0)
skuken1=size(kuken1)
skuken2=size(kuken2)
if (skuken1[0] ne 0) then L_H_a[kuken1]=0
if (skuken2[0] ne 0) then L_H_b[kuken2] =0
xrange=[min(d.sfr_tot), max (d.sfr_tot)]
yrange = [min (L_H_a(where(L_H_a ge 0))/4e26), max (d.L_H_a_ns/4e26) ]
plot, [1],[1],$
  xtitle=textoidl("SFR/(M_{"+sunsymbol()+"}/yr)"), $
  ytitle=textoidl("L_{H\alpha}/L_{"+sunsymbol()+"}"), $
  title=textoidl("H\alpha w/o dust vs. SFR"), $
  background = bgc, color = fgc,/nodata, $
  yrange = yrange,xrange=xrange
make_2dcolor, d.sfr_tot, d.L_H_a_ns/4e26,40,40,weight=weight
!P.MULTI[0]=0
plot, [1],[1],$
  xtitle="",ytitle="", $
  title="", background = bgc,$
  color = fgc,/nodata, $
  yrange = yrange,xrange=xrange,/noerase
fitx=indgen(100)/100.*(xrange[1]-xrange[0])+xrange[0]
fity=fitx/7.9e-35
oplot,fitx,fity/4e26,color=fgc

;subplot 2
!P.MULTI[0]=3
plot, [1],[1],$
  xtitle=textoidl("SFR/(M_{"+sunsymbol()+"}/yr)"), $
  ytitle=textoidl("L_{H\alpha}/L_{"+sunsymbol()+"}"), $
  title=textoidl("H\alpha w/dust vs. SFR"), $
  background = bgc, color = fgc,/nodata, $
  yrange = yrange,xrange=xrange
make_2dcolor, d.sfr_tot, d.L_H_a/4e26,40,40,weight=weight
!P.MULTI[0]=3
plot, [1],[1],$
  xtitle="",ytitle="", $
  title="", background = bgc,$
  color = fgc,/nodata, $
  yrange = yrange,xrange=xrange,/noerase
oplot,fitx,fity/4e26,color=fgc
legend,['K98'],/top,/left,color=[fgc],linest=0,box=0,textc=fgc

;subplot 3 - H_beta
!P.MULTI[0]=2
yrange = [min (L_H_b(where(L_H_b gt 0))/4e26), max (d.L_H_b_ns/4e26) ]
plot, [1],[1],$
  xtitle=textoidl("SFR/(M_{"+sunsymbol()+"}/yr)"), $
  ytitle=textoidl("L_{H\beta}/L_{"+sunsymbol()+"}"), $
  title=textoidl("H\beta w/o dust vs. SFR"), $
  background = bgc, color = fgc,/nodata, $
  yrange = yrange,xrange=xrange
make_2dcolor, d.sfr_tot, d.L_H_b_ns/4e26,40,40,weight=weight
!P.MULTI[0]=2
plot, [1],[1],$
  xtitle="",ytitle="", $
  title="", background = bgc, color = fgc,/nodata, $
  yrange = yrange,xrange=xrange,/noerase

;subplot 4 - H_beta
!P.MULTI[0]=1
plot, [1],[1],$
  xtitle=textoidl("SFR/(M_{"+sunsymbol()+"}/yr)"), $
  ytitle=textoidl("L_{H\beta}/L_{"+sunsymbol()+"}"), $
  title=textoidl("H\beta w/dust vs. SFR"), $
  background = bgc, color = fgc,/nodata, $
  yrange = yrange,xrange=xrange
make_2dcolor, d.sfr_tot, d.L_H_b/4e26,40,40,weight=weight
!P.MULTI[0]=1
plot, [1],[1],$
  xtitle="",ytitle="", $
  title="", background = bgc, color = fgc,/nodata, $
  yrange = yrange,xrange=xrange,/noerase

;oplot, d.sfr_tot, d.L_H_a_ns/4e26,color=c2,psym=1
;oplot, d.sfr_tot, d.L_H_b_ns/4e26,color=c1,psym=2
;oplot, d.sfr_tot, d.L_H_b/4e26,color=c4,psym=5
;legend,['H_alpha w/o dust','H_beta w/o dust','H_alpha','H_beta', $
;        'K98 H_alpha'], $
;       colors=[c2,c1,c3,c4,fgc],textcolor=fgc,psym=[1,2,4,5,-3],/left,box=0
!P.MULTI=0
save_plot, 'sfr-line.png'
;stop


; Halpha correcting for dust using balmer lines
; there:
; absorption in mags for Halpha using MW dust
abs_halpha=2.13*alog10(d.l_h_a/d.l_h_b/2.86)
halpha_dustcor=d.l_h_a/10^(-0.4*abs_halpha)

if !D.name eq 'PS' then device,file='hadust.eps'
xrange=[min(alog10(d.l_h_a_ns/4e26)), max (alog10(d.l_h_a_ns/4e26))]
yrange=[min(alog10(d.l_h_a/4e26)), max (alog10(d.l_h_a/4e26))]
plot,[1] ,[1],$
  xtitle=textoidl("log(intrinsic L_{H\alpha}/L_{"+sunsymbol()+"})"), $
  ytitle=textoidl("log(L_{H\alpha} with dust/L_{"+sunsymbol()+"})"), $
  title=textoidl ("H\alpha luminosity without dust correction"), $
  background = bgc, color = fgc,/nodata,xrange=xrange,yrange=yrange
make_2dcolor,alog10(d.l_h_a_ns/4e26),alog10(d.l_h_a/4e26),40,40,wei=weight
plot, [1],[1],$
  xtitle="",ytitle="", $
  title="", background = bgc, color = fgc,/nodata, $
  yrange = yrange,xrange=xrange,/noerase
oplot,[1,15],[1,15],color=fgc
save_plot, 'hadust.png'

if !D.name eq 'PS' then device,file='hadustcorr.eps'
;yrange=[min(alog10(halpha_dustcor/4e26)), max (alog10(halpha_dustcor/4e26))]
plot,[1] ,[1],$
  xtitle=textoidl("log(intrinsic L_{H\alpha}/L_{"+sunsymbol()+"})"), $
  ytitle=textoidl("log(L_{H\alpha} corrected for dust/L_{"+sunsymbol()+"})"), $
  title=textoidl ("H\alpha luminosity with dust correction"), $
  background = bgc, color = fgc,/nodata,xrange=xrange,yrange=yrange
make_2dcolor,alog10(d.l_h_a_ns/4e26),alog10(halpha_dustcor/4e26),40,40,wei=weight
plot, [1],[1],$
  xtitle="",ytitle="", $
  title="", background = bgc, color = fgc,/nodata, $
  yrange = yrange,xrange=xrange,/noerase

oplot,[1,15],[1,15],color=fgc
save_plot, 'hadustcorr.png'

; SFR - infrared luminosity
if !D.name eq 'PS' then device,file='sfr-fir.eps'
xrange=[min(alog10(d.sfr_tot)), max (alog10(d.sfr_tot))]
yrange=[min(alog10(d.l_fir_iras/4e26)), max (alog10(d.l_fir_iras/4e26))]
plot,[1],[1],$
  xtitle=textoidl("log(SFR/(M_{"+sunsymbol()+"}/yr))"), $
  ytitle=textoidl("log(L_{FIR,IRAS}/L_{"+sunsymbol()+"})"), $
  title="FIR luminosity vs. SFR ", background = bgc, color = fgc,/nodata, $
  xr=xrange,yr=yrange
make_2dcolor,alog10(d.sfr_tot),alog10(d.l_fir_iras/4e26),40,40,we=weight
plot,[1],[1],$
  background = bgc, color = fgc,/nodata, $
  /noer,xr=xrange,yr=yrange
fitx=10^((indgen(100)/99.)*6-3) ;[1e-2,1e-1,1,10,100,1e3]
;fitx=indgen(100)/100.*1000
fity=fitx/4.5e-37
oplot,alog10(fitx),alog10(fity/4e26),color=c2
; to account for K98 diff imf
oplot,alog10(fitx),alog10(2.06*fity/4e26),color=c2,line=1
legend,['K98',textoidl('2\cdotK98')],colors=[c2,c2],textcolor=0, $
  linestyle=[0,1],/left,box=0
save_plot, 'sfr-fir.png'


; SFR - UV luminosity
if !D.name eq 'PS' then device,file='sfr-uv.eps'
plot,alog10(d.sfr_tot),alog10(d.l_1600_ns/4e26),$
  xtitle=textoidl("log(SFR/(M_{"+sunsymbol()+"}/yr))"), $
  ytitle=textoidl("log(L_{1600}/L_{"+sunsymbol()+"})"), $
  title="intrinsic FUV luminosity vs. SFR ", background = bgc, color = fgc, $
  /nodata,yst=16
make_2dcolor,alog10(d.sfr_tot),alog10(d.l_1600_ns/4e26),40,40,we=weight
plot,alog10(d.sfr_tot),alog10(d.l_1600_ns/4e26),$
  background = bgc, color = fgc,/nodata, $
  /noer,yst=16
save_plot, 'sfr-uv.png'

; luv dust corr based on bell+kennicutt which is based on cksb94

; UV luminosity without dust correction
if !D.name eq 'PS' then device,file='uv-dust.eps'
yrange=[min(alog10(d.l_1600/4e26)),max(alog10(d.l_1600_ns/4e26))]
plot,alog10(d.l_1600_ns/4e26),alog10(d.l_1600/4e26),$
  xtitle=textoidl("log(intrinsic L_{1600}/L_{"+sunsymbol()+"})"), $
  ytitle=textoidl("log(observed L_{1600}/L_{"+sunsymbol()+"})"), $
  title="FUV luminosity without dust correction", $
  background = bgc, color = fgc,/nodata,yr=yrange
make_2dcolor,alog10(d.l_1600_ns/4e26),alog10(d.l_1600/4e26),40,40,we=weight
plot,alog10(d.l_1600_ns/4e26),alog10(d.l_1600/4e26),$
  background = bgc, color = fgc,/nodata, $
  /noer,yr=yrange
oplot,[1,100],[1,100],col=fgc
save_plot, 'uv-dust.png'

; UV luminosity with dust correction

abs_uv=0.5*8.66*(d.betaM95+1.71)/1.88
abs_uv=beta
i=where(beta lt -1, ct)
if ct ne 0 then abs_uv[i]=4.07*beta[i]+8.09
i=where(beta ge -1, ct)
if ct ne 0 then abs_uv[i]=0.97*beta[i]+4.71
luv_dustcor=d.l_1600/10^(-0.4*abs_uv)

if !D.name eq 'PS' then device,file='uv-dustcorr.eps'
plot,alog10(d.l_1600_ns/4e26),alog10(luv_dustcor/4e26),$
  xtitle=textoidl("log(intrinsic L_{1600}/L_{"+sunsymbol()+"})"), $
  ytitle=textoidl("log(observed L_{1600}/L_{"+sunsymbol()+"})"), $
  title="FUV luminosity with dust correction", $
  background = bgc, color = fgc,/nodata,yr=yrange
make_2dcolor,alog10(d.l_1600_ns/4e26),alog10(luv_dustcor/4e26),40,40,we=weight
plot,alog10(d.l_1600_ns/4e26),alog10(luv_dustcor/4e26),$
  background = bgc, color = fgc,/nodata, $
  /noer,yr=yrange
oplot,[1,100],[1,100],col=fgc
save_plot, 'uv-dustcorr.png'

; color-magnitude plot
; take magnitude 3 (g) vs. color 2-3 (u-g)
if !D.name eq 'PS' then device,file='cmd1.eps'
yrange = [max(d.(l [3])),min(d.(l [3]))]
plot, d.(l [2]) - d.(l [3]) , d.(l [3]), $
  ytitle=t[3],xtitle=t[2] + '-' + t[3], $
  title=t[3]+' vs. '+t[2] + '-' + t[3]+' CMD', $
  background = bgc, color = fgc,/nodata,yrange=yrange
make_2dcolor, d.(l [2]) - d.(l [3]) , d.(l [3]), 40,40,weight=weight
plot, d.(l [2]) - d.(l [3]) , d.(l [3]), $
                                ;ytitle=t[3],xtitle=t[2] + '-' + t[3], $
                                ;title=t[3]+' vs. '+t[2] + '-' + t[3]+' CMD', $
background = bgc, color = fgc,/nodata,yrange=yrange,/noer
oplot, d.(lns [2]) - d.(lns [3]) , d.(lns [3]), color = c2, psym = 3
;cam_plot, d.(l [2]) - d.(l [3]) , d.(l [3]), n_camera,c2
legend,['w/o dust','w/dust'], $
  colors=[c2,c1],textcolor=fgc,psym=[2,1],/right ,box=0
save_plot, 'cmd1.png'

; take magnitude 9 (K) vs. color 4-9 (r-K)
if !D.name eq 'PS' then device,file='cmd2.eps'
yrange = [max(d.(l [9])),min(d.(l [9]))]
plot, d.(l [4]) - d.(l [9]) , d.(l [9]), $
  ytitle=t[9],xtitle=t[4] + '-' + t[9], $
  title=t[9]+' vs. '+t[4] + '-' + t[9]+' CMD', $
  background = bgc, color = fgc,/nodata,yrange=yrange
make_2dcolor, d.(l [4]) - d.(l [9]) , d.(l [9]), 40,40,weight=weight
plot, d.(l [4]) - d.(l [9]) , d.(l [9]), $
                                ;ytitle=t[9],xtitle=t[4] + '-' + t[9], $
                                ;title=t[9]+' vs. '+t[4] + '-' + t[9]+' CMD', $
background = bgc, color = fgc,/nodata,yrange=yrange,/noer
oplot, d.(lns [4]) - d.(lns [9]) , d.(lns [9]), color = c2, psym = 3
legend,['w/o dust','w/dust'], $
  colors=[c2,c1],textcolor=fgc,psym=[2,1],/right ,box=0
save_plot, 'cmd2.png'

; take magnitude 4 (r) vs. color 2-4 (u-r)
if !D.name eq 'PS' then device,file='cmd3.eps'
yrange = [max(d.(l [4])),min(d.(l [4]))]
plot, d.(l [2]) - d.(l [4]) , d.(l [4]), $
  ytitle=t[4],xtitle=t[2] + '-' + t[4], $
  title=t[4]+' vs. '+t[2] + '-' + t[4]+' CMD', $
  background = bgc, color = fgc,/nodata,yrange=yrange
make_2dcolor, d.(l [2]) - d.(l [4]) , d.(l [4]), 40,40,weight=weight
plot, d.(l [2]) - d.(l [4]) , d.(l [4]), $
  background = bgc, color = fgc,/nodata,yrange=yrange,/noer
oplot, d.(lns [2]) - d.(lns [4]) , d.(lns [4]), color = c2, psym = 3
;legend,['w/o dust','w/dust'], $
;       colors=[c2,c1],textcolor=fgc,psym=[.,1],/right ,box=0
save_plot, 'cmd3.png'

; color-color plot
; take magnitude 2-3 (u-g) vs. color 2-3 (g-r)
if !D.name eq 'PS' then device,file='ccd1.eps'
;yrange = [max(d.(l [3])),min(d.(l [3]))]
plot, d.(l [2]) - d.(l [3]) , d.(l [3])-d.(l[4]), $
  ytitle=t[3]+'-'+t[4],xtitle=t[2] + '-' + t[3], $
  title=t[3]+' vs. '+t[2] + '-' + t[3]+' CCD', $
  background = bgc, color = fgc,/nodata
make_2dcolor, d.(l [2]) - d.(l [3]) , d.(l [3])-d.(l[4]), 40,40,weight=weight
plot, d.(l [2]) - d.(l [3]) , d.(l [3])-d.(l[4]), $
                                ;ytitle=t[3],xtitle=t[2] + '-' + t[3], $
                                ;title=t[3]+' vs. '+t[2] + '-' + t[3]+' CMD', $
background = bgc, color = fgc,/nodata,/noer
oplot, d.(lns [2]) - d.(lns [3]) , d.(lns [3])-d.(lns[4]), color = c2, psym = 3
;cam_plot, d.(l [2]) - d.(l [3]) , d.(l [3]), n_camera,c2
legend,['w/o dust','w/dust'], $
  colors=[c2,c1],textcolor=fgc,psym=[2,1],/right ,box=0
save_plot, 'ccd1.png'

;there:
; for H98 plots, load H98 data so we can plot it too
readcol,'~/heckmandata.txt',h98z,junk1,junk2,junk3,h98mb,junk4,h98lir, $
  h98beta,h98luv,h98irx,h98bol,junk5,junk6,lirll
h98bol=h98bol-alog10(4e33)
lli=where(lirll eq 1)
h98c=c2

; H98 fig 2a  beta vs Z
if !D.name eq 'PS' then device,file='beta-z.eps'
plot,d.betaH98,alog10(d.z_gas/0.02)+8.93, $
  xtitle=textoidl('\beta(H98)'),ytitle='12+log(O/H) ', $
  title=textoidl('\beta vs. metallicity (H98 Fig 2a)'), $
  background = bgc, $
  color = fgc,/nodata,xrange=[-3,2],yrange=[7,9.5]
make_2dcolor,d.betaH98,alog10(d.z_gas/0.02)+8.93,40,40,weight=weight
plot,d.betaH98,alog10(d.z_gas/0.02)+8.93, $
  color = fgc,/nodata,xrange=[-3,2],yrange=[7,9.5],/noer
;oplot,d.betaH98_ns,alog10(d.z_gas/0.02)+8.93,color=c2, psym = 3
oplot,h98beta,h98z,color=h98c,psym=4
legend,['Simulations','H98'], $
  colors=[c1,h98c],textcolor=fgc,linestyle=0, $
  psym=[-3,4],/left ,box=0,thick=[10,1]
save_plot, 'beta-z.png'


; H98 fig 2b  IRX vs Z
if !D.name eq 'PS' then device,file='irx-z.eps'
plot,[1],[1],$
  xtitle=textoidl('log(L_{FIR,IRAS}/L_{1900})'), $
  ytitle='12+log(O/H) ', $
  title='IRX vs. metallicity (H98 Fig 2b)', background = bgc, $
  color = fgc,/nodata,xrange=[-1,2.5],yrange=[7,9.5]
make_2dcolor,alog10(d.L_FIR_IRAS/d.L_1900),alog10(d.z_gas/0.02)+8.93,$
  40,40,weight=weight
plot,[1],[1],$
  color = fgc,/nodata,xrange=[-1,2.5],yrange=[7,9.5],/noer
;oplot,alog10(d.L_FIR_IRAS/d.L_1900_ns),alog10(d.z_gas/0.02)+8.93,color=c2, $
;      psym = 3
oplot,h98irx,h98z,color=h98c,psym=4
plotsym,6,2
oplot,h98irx[lli],h98z[lli],color=h98c,psym=8
legend,['Simulations','H98'], $
  colors=[c1,h98c],textcolor=fgc,linestyle=0, $
  psym=[-3,4],/left ,box=0,thick=[10,1]
save_plot, 'irx-z.png'


; H98 fig 8  Z vs LUV+LIR
if !D.name eq 'PS' then device,file='z-uvplusir.eps'
plot,[1],alog10((d.L_1900+d.L_FIR_IRAS)/4e26), $
  xtitle='12+log(O/H) ', $
  ytitle=textoidl('Log((L_{1900}+L_{FIR,IRAS})/L_{'+sunsymbol()+'})'),$
  title=textoidl('Metallicity vs L_{UV}+L_{FIR} (H98 Fig 8)'), $
  background = bgc, $
  color = fgc,/nodata,xrange=[7,9.5],yrange=[7,12]
make_2dcolor,alog10(d.z_gas/0.02)+8.93, $
  alog10((d.L_1900+d.L_FIR_IRAS)/4e26), $
  40,40,weight=weight
plot,[1],alog10((d.L_1900+d.L_FIR_IRAS)/4e26), $
  color = fgc,/nodata,xrange=[7,9.5],yrange=[7,12],/noer
;oplot,alog10(d.z_gas/0.02)+8.93,alog10((d.L_1900_ns)/4e26), $
;      color=c2, psym = 3
oplot,h98z,h98bol,color=h98c,psym=4
plotsym,1,2
oplot,h98z[lli],h98bol[lli],color=h98c,psym=8
legend,['Simulations','H98'], $
  colors=[c1,h98c],textcolor=fgc,linestyle=0, $
  psym=[-3,4],/left ,box=0,thick=[10,1]
save_plot, 'z-uvplusir.png'

; H98 fig 9a  beta vs LUV+LIR
if !D.name eq 'PS' then device,file='beta-uvplusir.eps'
plot,[1],alog10((d.L_1900+d.L_FIR_IRAS)/4e26), $
  xtitle=textoidl('\beta(H98)'), $
  ytitle=textoidl('Log((L_{1900}+L_{FIR,IRAS})/L_{'+sunsymbol()+'})'),$
  title=textoidl('UV slope vs L_{UV}+L_{FIR} (H98 Fig 9a)'), $
  background = bgc, $
  color = fgc,/nodata,xrange=[-3, 2],yrange=[7,12]
make_2dcolor, d.betaH98 ,alog10((d.L_1900+d.L_FIR_IRAS)/4e26), $
  40,40,weight=weight
plot,[1],alog10((d.L_1900+d.L_FIR_IRAS)/4e26), $
  color = fgc,/nodata,xrange=[-3, 2],yrange=[7,12],/noer
;oplot, d.betaH98_ns ,alog10((d.L_1900_ns)/4e26), $
;      color=c2, psym = 3
oplot,h98beta,h98bol,color=h98c,psym=4
plotsym,1,2
oplot,h98beta[lli],h98bol[lli],color=h98c,psym=8
legend,['Simulations','H98'], $
  colors=[c1,h98c],textcolor=fgc,linestyle=0, $
  psym=[-3,4],/left ,box=0,thick=[10,1]
save_plot, 'beta-uvplusir.png'


; H98 fig 9b  IRX vs LUV+LIR
if !D.name eq 'PS' then device,file='irx-uvplusir.eps'
plot,[1],alog10((d.L_1900+d.L_FIR_IRAS)/4e26), $
  xtitle=textoidl('log(L_{FIR,IRAS}/L_{1900})'), $
  ytitle=textoidl('Log((L_{1900}+L_{FIR,IRAS})/L_{'+sunsymbol()+'})'),$
  title =textoidl('IRX vs L_{UV}+L_{FIR} (H98 Fig 9b)'), $
  background = bgc, $
  color = fgc,/nodata,xrange=[-1, 2.5],yrange=[7,12]
make_2dcolor,alog10(d.L_FIR_IRAS/d.L_1900),alog10((d.L_1900+d.L_FIR_IRAS)/4e26), $
  40,40,weight=weight
plot,[1],alog10((d.L_1900+d.L_FIR_IRAS)/4e26), $
  color = fgc,/nodata,xrange=[-1, 2.5],yrange=[7,12],/noer
;oplot, alog10(d.L_FIR_IRAS/d.L_1900_ns) ,alog10((d.L_1900_ns)/4e26), $
;      color=c2, psym = 3
oplot,h98irx,h98bol,color=h98c,psym=4
plotsym,1,2
oplot,h98irx[lli],h98bol[lli],color=h98c,psym=8
plotsym,6,2
oplot,h98irx[lli],h98bol[lli],color=h98c,psym=8
legend,['Simulations','H98'], $
  colors=[c1,h98c],textcolor=fgc,linestyle=0, $
  psym=[-3,4],/left ,box=0,thick=[10,1]
save_plot, 'irx-uvplusir.png'


; Find blue magnitude, or if not available Sloan G
ref = where (t eq "AB_B_Johnson")
if ref eq -1 then begin
    ref = where (t eq "AB_G_SDSS")
    refr = where (t eq "AB_R_SDSS")
    ;; use Fukugita et al 96 transformation to get B from SDSS bands:
;; B = 1.42g - 0.42r + 0.44
    mb=1.42*d.(l[ref]) - 0.42*d.(l[refr]) + 0.44
    mbns=1.42*d.(l[ref]+1) - 0.42*d.(l[refr]+1) + 0.44
end else $
  ;; hmmm transformation from AB to Vega should be sth like +0.16
mb=d.(l[ref])+0.16
mbns=d.(l[ref]+1)+0.16

magtitle='Johnson B'

; H98 fig 12a  beta vs. M_b
if !D.name eq 'PS' then device,file='beta-mb.eps'
;yrange = [max(d.(l[ref])),min(d.(l [ref]+1))]
yrange=[-14.,-22.5]
plot,[1],[1],$
  xtitle=textoidl('\beta(H98)'), $
  ytitle=magtitle, $
  title='UV slope vs '+magtitle+' (H98 Fig 12a)', background = bgc, $
  color = fgc,/nodata,xrange=[-3, 1.5],yrange=yrange
;make_2dcolor, d.betaH98 ,d.(l[ref]), 40,40,weight=weight
make_2dcolor, d.betaH98 ,mb, 40,40,weight=weight
plot,[1],[1],$
                                ;xtitle='!7b!X(H98)',ytitle=t[ref], $
                                ;title='UV slope vs '+t[ref]+' (H98 Fig 12a)', background = bgc, $
color = fgc,/nodata,xrange=[-3, 1.5],yrange=yrange,/noer
;oplot, d.betaH98_ns ,mbns,$
;      color=c2, psym = 3
oplot,h98beta,h98mb,color=h98c,psym=4
legend,['Simulations','H98'], $
  colors=[c1,h98c],textcolor=fgc,linestyle=0, $
  psym=[-3,1],/left ,box=0,thick=[10,1]
save_plot, 'beta-mb.png'


; H98 fig 12b  Z vs. M_b
if !D.name eq 'PS' then device,file='z-mb.eps'
plot,[1],[1],$
  xtitle='12+log(O/H) ',ytitle= magtitle, $
  title='Metallicity vs ' +magtitle + ' (H98 Fig 12b)', background = bgc, $
  color = fgc,/nodata,xrange=[7,9.5],yrange= yrange
make_2dcolor,alog10(d.z_gas/0.02)+8.93, mb, 40, 40,weight=weight
plot,[1],[1],$
  background = bgc, $
  color = fgc,/nodata,xrange=[7,9.5],yrange= yrange,/noer
;oplot,alog10(d.z_gas/0.02)+8.93, mbns, $
;      color=c2, psym = 3
oplot,h98z,h98mb,color=h98c,psym=4
legend,['Simulations','H98'], $
  colors=[c1,h98c],textcolor=fgc,linestyle=0, $
  psym=[-3,4],/left ,box=0,thick=[10,1]
save_plot, 'z-mb.png'

; H98 fig 12c  IRX vs M_b
if !D.name eq 'PS' then device,file='irx-mb.eps'
plot,[1],[1],$
  xtitle=textoidl('log(L_{FIR,IRAS}/L_{1900})'), $
  ytitle= magtitle, $
  title='IRX vs. ' +magtitle+ ' (H98 Fig 12c)', background = bgc, $
  color = fgc,/nodata,xrange=[-1,2.5],yrange= yrange 
make_2dcolor,alog10(d.L_FIR_IRAS/d.L_1900), mb, $
  40, 40,weight=weight
plot,[1],[1],$
  background = bgc, $
  color = fgc,/nodata,xrange=[-1,2.5],yrange= yrange,/noer
;oplot,alog10(d.L_FIR_IRAS/d.L_1900_ns), mbns, $
;      color = c2, psym = 3
oplot,h98irx,h98mb,color=h98c,psym=4
plotsym,6,2
oplot,h98irx[lli],h98mb[lli],color=h98c,psym=8
legend,['Simulations','H98'], $
  colors=[c1,h98c],textcolor=fgc,linestyle=0, $
  psym=[-3,4],/left ,box=0,thick=[10,1]
save_plot, 'irx-mb.png'

; H98 fig 12d  LUV+LIR vs. M_b
if !D.name eq 'PS' then device,file='uvplusir-mb.eps'
plot,[1],[1],$
  xtitle=textoidl('Log((L_{1900}+L_{FIR,IRAS})/L_{'+sunsymbol()+'})'),$
  ytitle = magtitle, $
  title=textoidl('L_{UV}+L_{FIR}')+' vs. ' +magtitle+ ' (H98 Fig  12d)', $
  background = bgc, $
  color = fgc,/nodata,xrange=[7,12],yrange= yrange 
make_2dcolor,alog10((d.L_1900+d.L_FIR_IRAS)/4e26), mb, $
  40, 40,weight=weight
plot,[1],[1],$
  color = fgc,/nodata,xrange=[7,12],yrange= yrange ,/noer
;oplot,alog10((d.L_1900_ns)/4e26), mbns, $
;      color=c2, psym = 3
oplot,h98bol,h98mb,color=h98c,psym=4
plotsym,6,2
oplot,h98bol[lli],h98mb[lli],color=h98c,psym=8
legend,['Simulations','H98'], $
  colors=[c1,h98c],textcolor=fgc,linestyle=0, $
  psym=[-3,4],/left ,box=0,thick=[10,1]
save_plot, 'uvplusir-mb.png'

; make histograms

;there:
; luminosity
xrange=[min(d.l_scat/4e26),max(d.l_bol_grid/4e26) ]
plot,[1],[1],/xlog,xtitle="L/L!D"+sunsymbol()+"", $
  ytitle='time/yr',$
  title="Luminosity", background = bgc, color = fgc,/nodata, $
  xrange=xrange
if !D.name eq 'PS' then device,file='lhist.eps'
makehist,alog10(d.l_bol_grid/4e26),weight,20,hx,hlbol
makehist,alog10(d.l_bol_absorbed/4e26),weight,20,hx,hlabs
makehist,alog10(d.l_scat/4e26),weight,20,hx,hlscat
yrange = [0, max ([hlbol, hlabs, hlscat])]
plot,[1],[1],/xlog,xtitle="L/L!D"+sunsymbol()+"", $
  ytitle='time/yr',$
  title="Luminosity histogram", background = bgc, color = fgc, $
  xrange=xrange, yrange = yrange,/nodata
oplot,10^hx, hlbol,color = c1, psym=10
oplot,10^hx,hlabs,color = c2, psym=10
oplot,10^hx, hlscat,color = c3, psym=10
legend,['Bolometric','FIR','UV/VIS'],colors=[c1,c2,c3],textcolor=fgc,linestyle=0,/right,box=0
save_plot,'lhist.png'

; attenuation
attn= d.l_bol_absorbed/d.l_bol_grid
xrange=[min(attn),max(attn)]
plot,[1],[1], $
  xtitle="Attenuation",ytitle="time/yr", /xlog,xrange=xrange,$
  title="Attenuation histogram", background = bgc, color = fgc,/nodata
makehist, alog10 (attn),weight, 20, hx, ha
if !D.name eq 'PS' then device,file='atthist.eps'
plot,10^hx,ha, $
  xtitle="Attenuation",ytitle="time/yr", /xlog,xrange=xrange,$
  title="Attenuation histogram", background = bgc, color = fgc,/nodata
oplot,10^hx,ha, color = c1,psym=10
save_plot, 'atthist.png' 
;stop

; SFR
xrange=[min(d.sfr_tot),max(d.sfr_tot)]
plot,[1],[1], xrange=xrange,/nodata
makehist, d.SFR_tot,weight, 20, hx, h
if !D.name eq 'PS' then device,file='sfrhist.eps'
plot,hx,h, $
  xtitle="SFR/(M!D"+sunsymbol()+"!N/yr)", ytitle="time/yr",xrange=xrange,$
  title="SFR histogram", background = bgc, color = fgc,/nodata
oplot,hx,h, color = c1,psym=10
save_plot,'sfrhist.png'
; stop

; metallicity
xrange=[min(d.Z_gas ),max(d.z_gas )]
plot,[1],[1], xrange=xrange,/nodata
makehist, d.z_gas,weight, 20, hx, h
if !D.name eq 'PS' then device,file='zhist.eps'
plot,hx,h, $
  xtitle="Metallicity", ytitle="time/yr",xrange=xrange,$
  title="Metallicity histogram", background = bgc, color = fgc,/nodata
oplot,hx,h, color = c1,psym=10
save_plot,'zhist.png'
; stop
;there:

; Magnitudes
nn=0
for j=0,n_elements (l) - 1,nmags_per_plot do begin
                                ; find range
    xrange=[-100.,100.]
    for i = j,j+nmags_per_plot-1 do begin ;n_elements(l)-1 do begin
        if i lt n_elements (l) then begin
            if max(d.(lns[i])) gt 100 then $
              xrange[1] = min ([xrange[1],d.(l[i])]) $
            else $
              xrange[1] = min ([xrange[1],d.(lns[i])]) 
            xrange[0] = max ([xrange[0],d.(l[i])])
        end
    end

                                ; make histograms
    plot_name='mag'+ strcompress(string (nn,format = '(i2.1)'),/r) + 'hist'
    plot,[1],[1], xrange=xrange,/nodata
    nbins=20
    h=dblarr(nbins+2, nmags_per_plot)
    hns=dblarr(nbins+2, nmags_per_plot)
    k = - 1
    for i = j,j+nmags_per_plot-1 do begin ;n_elements(l)-1 do begin
        if i lt n_elements (l) then begin
            makehist, d.(l[i]),weight, nbins, hx, hh
            h[*,i-j]=hh
            makehist, d.(lns[i]),weight, nbins, hx, hh
            hns[*,i-j]=hh
            k = k+ 1
        end
    end
    if !D.name eq 'PS' then begin
        plot_name = plot_name+ '.eps'
        device,file=plot_name
    end else $
      plot_name = plot_name+ '.png'
    
    plot,hx,h[*,0], $
      xtitle="AB Magnitude",$
      ytitle="time/yr",xrange=xrange,yrange=[0,max([h,hns])],$
      title="Magnitude histogram", background = bgc, $
      color = fgc,/nodata

    k = - 1
                                ; first plot no dust with thick lines
    for i = j,j+nmags_per_plot-1 do begin 
        if i lt n_elements(l) then begin
            oplot,hx+0.01,hns[*,i-j],color=cc[i-j],psym=10,thick=4
            k = k+ 1
        end 
    end
                                ; then with dust with thin lines
    for i = j,j+nmags_per_plot-1 do begin 
        if i lt n_elements(l) then begin 
            oplot,hx,h[*,i-j],color=cc[i-j],psym=10
        end
    end

    legend,t[j:j+k], colors=cc[0:nmags_per_plot-1], $
      textcolor=fgc,linestyle=0,/right, $
      /top,box=0
    save_plot,plot_name
    nn = nn+ 1
                                ;stop
end

; dust attenuation in Magnitudes
nn=0
for j=0,n_elements (l) - 1,nmags_per_plot do begin
                                ; find range
    xrange=[100.,-100.]
    for i = j,j+nmags_per_plot-1 do begin ;n_elements(l)-1 do begin
        if i lt n_elements (l) then begin
            if max(d.(l[i]+1)) lt 100 then begin
                xrange[1] = max ([xrange[1], d.(l[i]) - d.(lns[i]) ])
                xrange[0] = min ([xrange[0], d.(l[i]) - d.(lns[i]) ])
            end
        end
    end
    
                                ; now make histograms
    plot_name='attn'+ strcompress(string (nn,format = '(i2.1)'),/r) + 'hist'
    plot,[1],[1], xrange=xrange,/nodata
    nbins=20
    h=dblarr(nbins+2, nmags_per_plot)
    k = - 1
    for i = j,j+nmags_per_plot-1 do begin ;n_elements(l)-1 do begin
        if i lt n_elements (l) then begin
            makehist, d.(l[i])-d.(lns[i]),weight, nbins, hx, hh
            h[*,i-j]=hh
            k = k+ 1
        end
    end
    if !D.name eq 'PS' then begin
        plot_name = plot_name+ '.eps'
        device,file=plot_name
    end else $
      plot_name = plot_name+ '.png'
    
    plot,hx,h[*,0], $
      xtitle="Dust attenuation/Magnitudes",$
      ytitle="time/yr",xrange=xrange,yrange=[0,max(h)],$
      title="Dust attenuation histogram", background = bgc, $
      color = fgc,/nodata
    k = - 1
    for i = j,j+nmags_per_plot-1 do begin ;n_elements(l)-1 do begin
        if i lt n_elements (l) then begin
            oplot,hx,h[*,i-j],color=cc[i-j],psym=10
            k = k+ 1
        end
    end
    
    legend,t[j:j+k], colors=cc[0:nmags_per_plot-1],textcolor=fgc,linestyle=0,/right, $
      /top,box=0
    save_plot,plot_name
    nn = nn+ 1
                                ;stop
end

; colors
nn=0
for j=0,n_elements (ll) - 1,nmags_per_plot do begin
                                ; find range
    xrange=[100.,-100.]
    for i = j,j+nmags_per_plot-1 do begin ;n_elements(l)-1 do begin
        if i lt n_elements (ll) then begin
            if (colors_nodust [0,i] lt 100) and $
              (colors_nodust[0,i] gt -100) then $ 
              xrange[0] = min ([xrange[0],colors_nodust [*,i]]) $
            else $
              xrange[0] = min ([xrange[0],colors [*,i]])
            xrange[1] = max ([xrange[1],colors[*,i]])
        end
    end
    plot_name='col'+ strcompress(string (nn,format = '(i2.1)'),/r) + 'hist'

                                ; make histograms
    plot,[1],[1], xrange=xrange,/nodata
    nbins=20
    h=dblarr(nbins+2, nmags_per_plot)
    hns=dblarr(nbins+2, nmags_per_plot)
    k = - 1
    for i = j,j+nmags_per_plot-1 do begin ;n_elements(l)-1 do begin
        if i lt n_elements (ll) then begin
            makehist, colors[*,i],weight, nbins, hx, hh
            h[*,i-j]=hh
            makehist, colors_nodust[*,i],weight, nbins, hx, hh
            hns[*,i-j]=hh
            k = k+ 1
        end
    end
    if !D.name eq 'PS' then begin
        plot_name = plot_name+ '.eps'
        device,file=plot_name
    end else $
      plot_name = plot_name+ '.png'
    plot,hx,h[*,0], $
      xtitle="Colors/Magnitudes",$
      ytitle="time/yr",xrange=xrange,yrange=[0,max([h,hns])],$
      title="Color histogram", background = bgc, $
      color = fgc,/nodata
    
    k=-1
                                ; first plot no dust with thick lines
    for i = j,j+nmags_per_plot-1 do begin 
        if i lt n_elements(ll) then begin
            oplot,hx+0.01,hns[*,i-j],color=cc[i-j],psym=10,thick=4
            k = k+ 1
        end 
    end
                                ; then with dust with thin lines
    for i = j,j+nmags_per_plot-1 do begin 
        if i lt n_elements(ll) then begin
            oplot,hx,h[*,i-j],color=cc[i-j],psym=10
        end
    end
    legend,tt[j:j+k], colors=cc[0:nmags_per_plot-1],textcolor=fgc,linestyle=0,/right,box=0,/top
    save_plot,plot_name
    nn = nn+ 1
                                ;stop
end

;stop
end


;-------------------------------------------------------

pro plot_seds,basename,ps=ps, overplot=overplot,plot_only=plot_only, $
              cols=cols, legtxt = legtxt,save_ascii=save_ascii, $
              no_title=no_title, bw=bw

graphics_init,ps=ps,plot_only=plot_only
fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7


spawn,"ls "+basename+"_*.fits",files
;files=basename
;files='mcrx_033.fits'
if keyword_set(save_ascii) then openw,1,save_ascii

;; figure out which HDU's
start_file= files [0]
scatlam_HDU = find_HDU (start_file, "scattering_lambdas")
lambda_HDU = find_HDU (start_file, "lambda")
iq_HDU = find_HDU (start_file, "integrated_quantities")
mcrx_HDU = find_HDU (start_file, "MCRX")

;; get number of cameras
h = HEADFITS (start_file, exten= mcrx_HDU)
n_camera = FXPAR (h, 'n_camera')

;; get wavelengths
ftab_ext, start_file, ext = iq_HDU, "lambda", lambda
;; the nonscatter SED's have lambda like the images, hokey
ftab_ext, start_file, ext = lambda_HDU, "lambda", lambda_ns

ftab_ext, start_file, ext = scatlam_HDU, "entry", lambda_entries
lambda_entries = lambda_entries (1:*)
junk=where(lambda ge 5e-4)
last_lambda=junk[0]-1

if not keyword_set (cols) then $
  cols=indgen(n_camera)

if not keyword_set (overplot) then begin
    colors=intarr(n_elements(cols))+c2
end else begin
    colors=indgen(n_elements(cols))+3
    colors[0:1]=[c3,c2]
    if not keyword_set (legtxt) then begin
        legtxt='w/o dust'
        for cc=0,n_elements (cols) - 1 do $
          legtxt = [legtxt, 'camera ' + strcompress (string (cols[cc] ),/remove_all)]
    end
end

if keyword_set(bw) then begin
    print,"Setting b/w color table"
    colors[*]=fgc
    c1=fgc
    c2=fgc
    c3=fgc
    c4=fgc
    lines=[1,2,3,4,5,6,7,8,9,10,11]
end else lines=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

print,n_elements (files) 
for n = 0,n_elements (files) - 1 do begin
    fits_file = files[n]

    h = HEADFITS (fits_file, exten= iq_HDU)
    print, ' Loading ', fits_file
    d=mrdfits(fits_file,iq_HDU)
    if not keyword_set(no_title) then title_base= 'SED: '+fits_file
                                ;for c = 0,n_camera - 1 do begin
    for cc = 0,n_elements(cols)-1 do begin
        c=cols[cc]
        ;; get data from fits file
        colstring = strcompress (string (c, format = '(i 2.1)' ),/remove_all)
        colname_ns ="L_LAMBDA_NONSCATTER"+colstring
        colname="L_LAMBDA_SCATTER"+colstring 
                                ;ftab_ext, fits_file, ext = iq_HDU, colname_ns, L_lambda_ns
                                ;ftab_ext, fits_file, ext = iq_HDU, colname, L_lambda
        llnse=where(tag_names(d) eq colname_ns)
        lle=where(tag_names(d) eq colname)
        L_lambda_ns=d.(llnse)
        L_lambda=d.(lle)
        ;; plot SED
        if not (keyword_set (overplot) and c gt 0) then begin
            plot_name= 'SED_'+fits_file+"-"+ $
              string (c, format = '(i 2.2)' )
            if keyword_set (ps) then $
              device,filename=plot_name+'.eps'
            yrange = $
              [min (L_lambda[lambda_entries[1]:last_lambda]* $
                    lambda[lambda_entries[1]:last_lambda]), $
               max  ([L_lambda_ns*lambda,L_lambda*lambda])]
            if not keyword_set(no_title) then begin
                if not keyword_set (overplot) then $
                  title = title_base +" camera "+ $
                  string (c, format = '(i 2.1)' ) $
                else title = title_base 
            end else title=''

            plot,[1],[1], xrange = [1e-8,1e-3], yrange=yrange,$
              ystyle = 16 ,xtitle=textoidl('\lambda/m'), $
              ytitle=textoidl("\lambdaL_\lambda/W"), /xlog,/ylog,$
              title=title,$
              background = bgc,color=fgc,/nodata
            oplot,lambda_ns, lambda_ns*L_lambda_ns,color=c1
            if keyword_set (overplot) then begin 
                legend,legtxt, colors=[c1,colors], $
                  textcolor=fgc,linestyle=0,/right,box=0
                print 
            end else $
              legend,['w/o dust','with dust'], $ 
              colors=[c1,c2],textcolor=fgc,linestyle=[0,lines[0]],$
              /right,box=0
        end 
                                ;oplot,lambda, lambda*L_lambda_ns,color=c3
        
        oplot,lambda, lambda*L_lambda, color = colors[cc],lines=lines[cc]

                                ; save to file if we want
        if keyword_set(save_ascii) then begin
            printf,1,"# File: "+fits_file+", camera: "+colstring
            printf,1,"# lambda/m           L_lambda/W   L_lambda(no dust)/W"
            for m=0,n_elements(lambda)-1 do $
              printf,1,lambda[m],L_lambda[m],L_lambda_ns[m]
            prinfgc=0
            bgc=1
            c2=2
            c1=3
            c3=4
            c4=5
            tf,1,""
            printf,1,""
        end

        if not (keyword_set (overplot) and (cc ne n_elements(cols)-1)) then $
          save_plot,plot_name+'.png'
    end 
end

if keyword_set(save_ascii) then close,1

end


;-------------------------------------------------------

pro plot_grid,file

;set_plot,'z'
;device,pseudo_color=8,set_character_size=[9,12]
device,set_character_size=[9,12]
fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7

!P.MULTI=0
colr=[0,255,255,0  ,0  ,175,143,  0,  0,  0,155,  0]
colg=[0,255,85 ,85 ,235,0  ,255,255,175,130,  0,  0]
colb=[0,255,0  ,255,0  ,175,  0,143,255,186,255,255]
tvlct,colr,colg,colb

sz = size (file) 
if sz [0] eq 0 then file = [file]

;; figure out which HDU's
gd_HDU = find_HDU (file [0], "griddata")

; plot histogram of gas density.  We want it volume-weighted, so it's
; a little more painful than it ought to be
plot,[1],[1],/xl,/yl,xtitle='Gas density/(Msun/kpc3)',ytit='volume/kpc3', $
  /nodata, xrange = [1e-4, 1e10], yrange = [ 1e-4,1e6]
oplot,[1e7,1e7],[1e-4,1e10],linest=1

for f=0, n_elements (file) - 1 do begin
    
    a=mrdfits(file[f],gd_HDU)

    mg=a.mass_gas
    v=a.cell_volume

    nbins=100
    h=histogram(alog10(mg/v),nbins=nbins,omin=x1,omax=x2,reverse=r)
    x=indgen(nbins)/(nbins-1.0)*(x2-x1)+x1

; now do volume weighting
    hh=fltarr(nbins)
    for i=0,nbins-1 do begin
        ;; these are the cells that contributed to bin i
        cur=r[r[i]:r[i+1]-1]
        ;; add their volume to the i'th bin
        hh[i]=total(v[cur])
    end

    oplot,10^x,hh

end

stop
end


;-------------------------------------------------------

pro inclinations,file,ps=ps,plot_only=plot_only

graphics_init,ps=ps,plot_only=plot_only
fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7


a=mrdfits(file,1)

mb=a.ab_b_johnson
mbn=a.ab_b_johnson_nodust
mi=a.ab_i_cousins
min =a.ab_i_cousins_nodust
theta = a.theta

db=fltarr(11)
di=fltarr(11)
eb=fltarr(11)
ei=fltarr(11)
incl=fltarr(11)
for i=0,10 do begin
    c=where(a.camera eq i)
    db[i]=mean(mb[c]-mbn[c])
    di[i]=mean(mi[c]-min[c])
    eb[i]=stdev(mb[c]-mbn[c])
    ei[i]=stdev(mi[c]-min[c])
    incl [i] =theta[c[0]]
end

; read ferrara data
;readcol,'ferraradata.txt',f99i,f99Ab,f99Ai
f99i=[0]
f99Ab=[0]
f99Ai=[0]

if !D.name eq 'PS' then device,file='Sbc_inclinations.eps'
plot,1-cos(incl),di-di[5],xr=[0,1],yr=[-0.2,2.1],yst=3, $
  back=bgc,color=fgc,/nodata, xtitle ='1-cos(i)', $
  ytitle=textoidl('\DeltaM_{dust}') 
oplot,1-cos(incl),di,$
  color=c1
; range of observed inclinations
e=indgen(40)/39.
obsr=where(e ge 0.3 and e le 0.85)
;oplot,1-cos(incl),di+ei,line=1,color=c1
;oplot,1-cos(incl),di-ei,line=1,color=c1
oplot,e[obsr],1.37*(e[obsr])-0.45,color=c2 
oplot,e[obsr],0.8*alog10(1/(1-e[obsr])),color=c3,line=2
oplot,e[obsr],1.05*alog10(1/(1-e[obsr]))-0.1,color=c4 
; oplot,e[obsr],1.23*(e[obsr])-0.45,color=c2,line=1
; oplot,e[obsr],1.51*(e[obsr])-0.45,color=c2,line=1
oplot,1-cos(incl),db+0.2,color=c1
;oplot,1-cos(incl),db+0.2-eb,line=1,color=c1
;oplot,1-cos(incl),db+0.2+eb,line=1,color=c1
oplot,e[obsr],1.45*alog10(1/(1-e[obsr]))+0.25,color=c2,line=2

oplot,1-cos(f99i*3.14159/180),-2.5*alog10(f99Ab)+0.2,col=c3
oplot,1-cos(f99i*3.14159/180),-2.5*alog10(f99Ai),col=c3

xyouts,0.55,0.35,'I',col=fgc,chars=2
xyouts,0.55,1.45*alog10(1/0.4)+0.25,'B',chars=2,col=fgc
legend,['Simulation','Ferrara et al (99)','deVaucouleurs et al (91)',$
        'Bernstein et al (94)', $
        'Han (92)','Giovanelli et al (94)'],/left,/top, $
  col=[c1,c3,c2,c2,c3,c4],textc=fgc,box=0,$
  lines=[0,0,2,0,2,0]
save_plot,'Sbc_inclinations.png'
stop
end

;-------------------------------------------------------

pro ferrara,ps=ps,plot_only=plot_only
file='index.fits'
graphics_init,ps=ps,plot_only=plot_only
fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7


a=mrdfits(file,1)

mb=a.ab_b_johnson
mbn=a.ab_b_johnson_nodust
mi=a.ab_i_cousins
min =a.ab_i_cousins_nodust
theta = a.theta

db=fltarr(11)
di=fltarr(11)
eb=fltarr(11)
ei=fltarr(11)
incl=fltarr(11)
for i=0,10 do begin
    c=where(a.camera eq i)
    db[i]=mean(mb[c]-mbn[c])
    di[i]=mean(mi[c]-min[c])
    eb[i]=stdev(mb[c]-mbn[c])
    ei[i]=stdev(mi[c]-min[c])
    incl [i] =theta[c[0]]
end

; read ferrara data 
; files contain attn for the spec band
; sets of 8 different tau_v (0.1 0.5 1 2 5 10 20 50)
;    sets of 9 different inclinations
;       for 5 diff B/T (0, 0.1, 0.3, 0.5, 1)

;readcol,'S01_ME04-B',f99Ab
;readcol,'S01_ME04-I',f99Ai
readcol,'S04_MD20-B',f99Ab
readcol,'S04_MD20-I',f99Ai
f99i=[9.3,22.9, 30, 40, 50, 60, 70, 80, 90]*3.14/180
f99Ab=-2.5*alog10(f99Ab)
f99Ai=-2.5*alog10(f99Ai)

if !D.name eq 'PS' then device,file='f99_incl.eps'
plot,1-cos(incl),di-di[5],xr=[0,1],yr=[-0.2,2.1],yst=3, $
  back=bgc,color=fgc,/nodata, xtitle ='1-cos(i)', $
  ytitle=textoidl('\DeltaM_{dust}') ,$
  title='I band'
oplot,1-cos(incl),di,color=c1
;oplot,1-cos(incl),db,color=c1

for j=0,7 do begin
    i=72*1+j
                                ;oplot,1-cos(f99i),f99ab[i:i+72:8],col=c3
    oplot,1-cos(f99i),f99ai[i:i+72:8],col=c3
end

save_plot,'f99_incl.png'
stop
end

;-------------------------------------------------------

pro  strides,ps=ps,plot_only=plot_only

graphics_init,ps=ps,plot_only=plot_only
fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7

readcol,'stride_v194_83_sorted',x,y
readcol,'stride_v1_94_1_sorted',hx,hy

if !D.name eq 'PS' then device,file='strides.eps'
plot,hy,hx,/xl,back=bgc,/nodata,xtit='stride/bytes',col=fgc,ytit='Cumulative fraction of cells'
oplot,hy,hx,col=c1
oplot,y,x,col=c2
oplot,[4096,4096],[0,1],line=1,col=fgc
oplot,[16384,16384],[0,1],line=2,col=fgc
legend,['Default ordering','Hilbert ordering','4kb',$
       '16kb'], $
  col=[c2,c1,0,0],/right,/bot,line=[0,0,1,2],textc=fgc,box=0
save_plot,'strides.png'
end


;-------------------------------------------------------


; makes the luminosity-attenuation plots
pro lumatt,files,legtxt=legtxt,ps=ps,plot_only=plot_only,bymass=bymass,bw=bw

graphics_init,ps=ps,plot_only=plot_only

fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7
if keyword_set(bw) then begin
    c1=fgc
    c2=fgc
    c3=fgc
    c4=fgc
    c5=fgc
    c6=fgc
end

if not keyword_set(legtxt) then legtxt=files
if keyword_set (ps) then device,file='l-att.eps'

sim=['Sbc2*','Sbcm2*','Sbcp2*','G3G3*','G2G2*','G1G1*','G0G0*','Sbc11i*']
mb=[20e10,10e10,31.2e10,12.4e10,4e10,1.4e10,3.2e9,10e10]


d=mrdfits(files[0],1)
for s=0,n_elements(sim)-1 do $
  if strmatch(files[0],sim[s]) eq 1 then mass=mb[s]
print,mass

met=d.z_gas+d.z_gas_new*(d.metal_to_dust_ratio/0.4-1)

if keyword_set(bymass) then $
  plot,d.l_bol_grid/4e26/mass , $
  (d.l_bol_absorbed/d.l_bol_grid),$
  psym=4,/xl,back=1,col=0,$
  syms=0.4,xtitle=textoidl("(L/M)/(L_"+sunsymbol()+"/M_"+sunsymbol()+")"), $
  ytitle='Attenuation',$
  xra=[0.3,4],xst=1 $
else $
                                ;plots,alog10(d.l_bol_grid/4e26),met,d.l_bol_absorbed/d.l_bol_grid,ps=4,/t3d,col=0
plot,d.l_bol_grid/4e26, $
  (d.l_bol_absorbed/d.l_bol_grid),$
  psym=4,/xl,back=1,col=0,$
  syms=0.4,xtitle=textoidl("L_{bol}/L_{"+sunsymbol()+"}"), $
  ytitle=textoidl('L_{abs}/L_{bol}'),$
  xra=[3e10,1.5e12],xst=1
symbols=[4]
sizes=[0.4]
colors=[0]

for i = 1,n_elements(files)-1 do begin
    dd=mrdfits(files[i],1)
    met=dd.z_gas+dd.z_gas_new*(dd.metal_to_dust_ratio/0.4-1)
                                ;d=[d,dd]
    for s=0,n_elements(sim)-1 do $
      if strmatch(files[i],sim[s]) eq 1 then begin
        mass=mb[s]
        print,sim[s]
        if sim[s] eq 'Sbc11i*' then syms=1 else syms=0.4
    end
    print,mass
    if syms eq 1 then psym=5 else if i gt 6 then psym=5 else psym=4
    if i gt 6 then color = i-5 else color=i+1
    if keyword_set(bw) then color=fgc
    symbols=[symbols,psym]
    sizes=[sizes,syms]
    colors=[colors,color]
    if keyword_set(bymass) then $
      oplot,dd.l_bol_grid/4e26/mass ,$
      (dd.l_bol_absorbed/dd.l_bol_grid),psym=psym,color=color,syms=syms $
    else $
                                ;plots,alog10(dd.l_bol_grid/4e26),met,dd.l_bol_absorbed/dd.l_bol_grid,/t3d,col=color,psym=psym,syms=syms
    oplot,dd.l_bol_grid/4e26,$
      (dd.l_bol_absorbed/dd.l_bol_grid),psym=psym,color=color,syms=syms 
end

x=(indgen(100)/100.*3.+10.)
y=1-1./(x-9.9)/2.6
if keyword_set(bymass) then x=x-alog10(20e10)
;oplot,10^x,y,col=fgc
;stop
;legcol=indgen(n_elements(files))+1
;legcol[0]=0
if n_elements(legtxt) gt 1 then $
  legend,legtxt,col=colors,textc=0, $
  /bott,/rig,psym=symbols,syms=sizes*2,$
  box=0,spacing=1.5
;stop

if keyword_set (ps) then device,/cl
end

;-------------------------------------------------------

pro makelumattplot
fi=['Sbc201a-u4/Zsolar-imf2.35/index.fits',$
                                ;'Sbcm201-u4/Zsolar-imf2.35/index.fits',$
                                ;'Sbcp201-u4/Zsolar-imf2.35/index.fits',$
'Sbc202-u4/Zsolar-imf2.35/index.fits',$
  'Sbc203-u4/Zsolar-imf2.35/index.fits','Sbc204-u4/Zsolar-imf2.35/index.fits',$
  'Sbc205-u4/Zsolar-imf2.35/index.fits','Sbc206-u4/Zsolar-imf2.35/index.fits',$
  'Sbc217-u4/Zsolar-imf2.35/index.fits','Sbc218-u4/Zsolar-imf2.35/index.fits',$
  'Sbc219-u4/Zsolar-imf2.35/index.fits']
lumatt,fi,/ps,/bw,leg=['Sbc201'] ;,$;'Sbcm201','Sbcp201',
;'Sbc202','Sbc203','Sbc204',$
;                   'Sbc205','Sbc206','Sbc217','Sbc218','Sbc219']
end 
;-------------------------------------------------------

pro makelovermattplot
fi=['Sbc201a-u4/Zsolar-imf2.35/index.fits',$
    'Sbcm201-u4/Zsolar-imf2.35/index.fits',$
    'Sbcp201-u4/Zsolar-imf2.35/index.fits',$
    'G3G3b-u1/Zsolar-imf2.35/index.fits',$
    'G2G2-u1/Zsolar-imf2.35/index.fits',$
    'G1G1a-u1/Zsolar-imf2.35/index.fits',$
    'G0G0a-u1/Zsolar-imf2.35/index.fits',$
    'Sbc11i4-u4/Zsolar-imf2.35/index.fits']
lumatt,fi,/pl,/bymass,$
  leg=['Sbc201','Sbcm201','Sbcp201','G3G3','G2G2','G1G1',$
       'G0G0',$
       'Sbc (isol)']
end
;-------------------------------------------------------

pro makelovermattzplot
fi=['Sbc201a-u4/Zsolar-imf2.35/index.fits',$
    'Sbc201a-u4/cZsolar-imf2.35/index.fits',$
    'Sbc201a-u4/0.4Zsolar-imf2.35/index.fits',$
    'Sbc201a-u4/0.4cZsolar-imf2.35/index.fits',$
    'Sbc201a-u4/0.2Zsolar-imf2.35/index.fits',$
    'Sbc201a-u4/0.2cZsolar-imf2.35/index.fits']

for a=30,60,5 do begin
    erase,1
    scale3,xr=[10.5,12],yr=[0,0.03],zr=[0,1],ax=a,az=-90
    lumatt,fi,/pl,$
      leg=['Z!Dsolar!N',$
           'const Z!Dsolar!N',$
           '0.4 Z!Dsolar!N',$
           'const 0.4 Z!Dsolar!N',$
           '0.2 Z!Dsolar!N',$
           'const 0.2 Z!Dsolar!N']
    print,"angle",a
;       leg=['Z!D'+sunsymbol()+'!N',$
;            '0.4 Z!D'+sunsymbol()+'!N',$
;            '0.2 Z!D'+sunsymbol()+'!N']
end
end



;-------------------------------------------------------

pro cumlum,files,legtxt=legtxt,ps=ps,plot_only=plot_only
if not keyword_set(legtxt) then legtxt=files
graphics_init,ps=ps,plot_only=plot_only
fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7

d=mrdfits(files[0],1)
weight=d.dt/11
xrange=[min(d.l_scat/4e26),max(d.l_bol_grid/4e26) ]
;xrange=[min(d.l_scat/4e26),1e12 ]
plot,[1],[1],/xlog,xtitle="L/L!D"+sunsymbol()+"", $
  ytitle='time/yr',$
  title="Luminosity", background = 1, color = 0,/nodata, $
  xrange=xrange,xst=1,/yl
makehist,alog10(d.l_bol_grid/4e26),weight,40,hx,hlbol
makehist,alog10(d.l_bol_absorbed/4e26),weight,40,hx,hlabs
makehist,alog10(d.l_scat/4e26),weight,40,hx,hlscat

for i=n_elements(hx)-2,0,-1 do begin
    hlbol[i]=hlbol[i]+hlbol[i+1]
    hlabs[i]=hlabs[i]+hlabs[i+1]
    hlscat[i]=hlscat[i]+hlscat[i+1]
end
yrange = [1e7, max ([hlbol, hlabs, hlscat])]
if keyword_set (ps) then device,file='cumlum.eps'
plot,[1],[1],/xlog,xtitle="L/L!D"+sunsymbol()+"", $
  ytitle='time/yr',$
  title="Cumulative Luminosity Distribution", $
  background = 1, color = 0, $
  xrange=xrange, yrange = yrange,/nodata,/yl,xst=1
oplot,10^hx, hlbol+1,color = 2, psym=-3
oplot,10^hx,hlabs+1,color = 3, psym=-3
oplot,10^hx, hlscat+1,color = 4, psym=-3

thick=[1,1,1,1]
line=[0,0,0,0]
colors=[2,3,4,2]
for j = 1,n_elements(files)-1 do begin
    d=mrdfits(files[j],1)
    weight=d.dt/11
    
    makehist,alog10(d.l_bol_grid/4e26),weight,40,hx,hlbol
    makehist,alog10(d.l_bol_absorbed/4e26),weight,40,hx,hlabs
    makehist,alog10(d.l_scat/4e26),weight,40,hx,hlscat
    
    for i=n_elements(hx)-2,0,-1 do begin
        hlbol[i]=hlbol[i]+hlbol[i+1]
        hlabs[i]=hlabs[i]+hlabs[i+1]
        hlscat[i]=hlscat[i]+hlscat[i+1]
    end
    oplot,10^hx, hlbol+1,color = 2, psym=-3,thick=2,line=j
    oplot,10^hx,hlabs+1,color = 3, psym=-3,thick=2,line=j
    oplot,10^hx, hlscat+1,color = 4, psym=-3,thick=2,line=j
    thick=[thick,2]
    line=[line,j]
    colors=[colors,2]
end
legend,['Bolometric','FIR','UV/VIS',legtxt],colors=colors,textcolor=0,/right,box=0,thick=thick,line=line
if keyword_set (ps) then device,/cl

end

;-------------------------------------------------------

pro makecumlumplot
fi=['Sbc201a-u4/Zsolar-imf2.35/index.fits',$
    'Sbc202-u4/Zsolar-imf2.35/index.fits',$
    'Sbc203-u4/Zsolar-imf2.35/index.fits',$
    'Sbc206-u4/Zsolar-imf2.35/index.fits']
cumlum,fi,/ps,leg=['Sbc201','Sbc202','Sbc203','Sbc206']
spawn,'mv cumlum.eps cumlum1.eps'

fi=['Sbc201a-u4/Zsolar-imf2.35/index.fits',$
    'Sbc204-u4/Zsolar-imf2.35/index.fits',$
    'Sbc205-u4/Zsolar-imf2.35/index.fits']
cumlum,fi,/ps,leg=['Sbc201','Sbc204','Sbc205']
spawn,'mv cumlum.eps cumlum2.eps'

fi=['Sbc203-u4/Zsolar-imf2.35/index.fits',$
    'Sbc218-u4/Zsolar-imf2.35/index.fits',$
    'Sbc219-u4/Zsolar-imf2.35/index.fits']
cumlum,fi,/ps,leg=['Sbc203','Sbc218','Sbc219']
spawn,'mv cumlum.eps cumlum3.eps'

end

;-------------------------------------------------------

pro cumsfr,files,legtxt=legtxt,ps=ps,plot_only=plot_only
if not keyword_set(legtxt) then legtxt=files
graphics_init,ps=ps,plot_only=plot_only
fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7

d=mrdfits(files[0],1)
weight=d.dt/11
xrange=[min(d.sfr_tot),max(d.sfr_tot) ]
xrange=[1,650]
plot,[1],[1],/xlog,xtitle="L/L!D"+sunsymbol()+"", $
  ytitle='time/yr',$
  title="Luminosity", background = 1, color = 0,/nodata, $
  xrange=xrange,xst=1,/yl
makehist,alog10(d.sfr_tot),weight,40,hx,hsfr

for i=n_elements(hx)-2,0,-1 do begin
    hsfr[i]=hsfr[i]+hsfr[i+1]
end
yrange = [1e7, max (hsfr)]
if keyword_set (ps) then device,file='cumsfr.eps'
plot,[1],[1],/xlog,xtitle="SFR/M!D"+sunsymbol()+"/yr", $
  ytitle='time/yr',$
  title="Cumulative SFR Distribution", $
  background = 1, color = 0, $
  xrange=xrange, yrange = yrange,/nodata,/yl,xst=1
oplot,10^hx, hsfr+1,color = 2, psym=-3

thick=[1]
line=[00]
colors=[2]
for j = 1,n_elements(files)-1 do begin
    d=mrdfits(files[j],1)
    weight=d.dt/11

    makehist,alog10(d.sfr_tot),weight,40,hx,hsfr

    for i=n_elements(hx)-2,0,-1 do begin
       hsfr[i]=hsfr[i]+hsfr[i+1]
    end
                    ; should stop line at the time where the
                    ; simulation ends, not plot a long flat line
                    ; that's unphysical
    oplot,10^hx, hsfr+1,color = 2, psym=-3,thick=2,line=j
    thick=[thick,2]
    line=[line,j]
    colors=[colors,2]
end
legend,legtxt,colors=colors,textcolor=0,/right,box=0,thick=thick,line=line
if keyword_set (ps) then device,/cl

end

;-------------------------------------------------------

pro lummag,files,legtxt=legtxt,ps=ps,plot_only=plot_only,band=band
graphics_init,ps=ps,plot_only=plot_only
fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7

if not keyword_set(legtxt) then legtxt=files

;band=3

sim=['Sbc2*','Sbcm2*','Sbcp2*','G3G3*','G2G2*','G1G1*','Sbc11i*']
mb=[20e10,10e10,31.2e10,12.4e10,4e10,1.4e10,10e10]

d=mrdfits(files[0],1)
for s=0,n_elements(sim)-1 do $
  if strmatch(files[0],sim[s]) eq 1 then mass=mb[s]
print,mass
;plot,d.  L_bol_grid ,d.ab_r_sdss,psym=4,/xl,back=1,col=0,$
;plot,d.l_bol_grid, - (d.l_bol_absorbed-d.l_bol_grid),psym=4,/xl,back=1,col=0,$

f=tag_names(d)
l= 0
t=''
for i = 0,n_elements(f)-1 do begin
    if (strmid(f[i],0,3) eq 'AB_') and $
      not (strmid(f[i],6,7,/reverse) eq '_NODUST') then begin
        print,f[i]
        sz=size(l)
        if sz[0] eq 0 then begin
            l = [i]
            t = [f[i]]
        end else  begin
            l=[l,i]
            t = [t, f [i]]
        end
    end
end
yrange=[max(d.(l[ band])),min(d.(l[ band]))]
for i = 1,n_elements(files)-1 do begin
    dd=mrdfits(files[i],1)
    yrange[0]=max([yrange[0],dd.(l[band])])
    yrange[1]=min([yrange[1],dd.(l[band])])
end

plot_name='l-mag'+ strcompress(string (band,format = '( i 2.1)'),/r)
if !D.name eq 'PS' then device,file=plot_name+'.eps'
ytitle=t[band]+" Magnitude"

plot,d.l_bol_grid/4e26 , $
  d.(l[band]),$
  psym=4,/xl,back=1,col=0,$
  syms=0.4,xtitle="L/L!D"+sunsymbol()+"!N", $
  ytitle=ytitle,yrange=yrange,$
  xra=[4e10,2e12],xst=1
                                ;xra=[0.5,4],xst=1
symbols=[4]
sizes=[0.4]
colors=[0]
for i = 1,n_elements(files)-1 do begin
    dd=mrdfits(files[i],1)
                                ;d=[d,dd]
    for s=0,n_elements(sim)-1 do $
      if strmatch(files[i],sim[s]) eq 1 then begin
        mass=mb[s]
        print,sim[s]
        if sim[s] eq 'Sbc11i*' then syms=1 else syms=0.4
    end
    print,mass
    if syms eq 1 then psym=5 else if i gt 6 then psym=5 else psym=4
    if i gt 6 then color = i-5 else color=i+1
    symbols=[symbols,psym]
    sizes=[sizes,syms]
    colors=[colors,color]
    oplot,dd.l_bol_grid/4e26, $
      dd.(l[band]),psym=psym,color=color,syms=syms
end
legcol=indgen(n_elements(files))+1
legcol[0]=0
if n_elements(legtxt) gt 1 then $
  legend,legtxt,col=colors,textc=0, $
  /bott,/rig,psym=symbols,syms=sizes*2,$
  box=0,spacing=1.5
;stop

if keyword_set (ps) then device,/cl
end

;-------------------------------------------------------

pro makelummagplot
fi=['Sbc201a-u4/Zsolar-imf2.35/index.fits',$
                                ;'Sbcm201-u4/Zsolar-imf2.35/index.fits',$
                                ;'Sbcp201-u4/Zsolar-imf2.35/index.fits',
'Sbc202-u4/Zsolar-imf2.35/index.fits',$
  'Sbc203-u4/Zsolar-imf2.35/index.fits','Sbc204-u4/Zsolar-imf2.35/index.fits',$
  'Sbc205-u4/Zsolar-imf2.35/index.fits','Sbc206-u4/Zsolar-imf2.35/index.fits',$
  'Sbc217-u4/Zsolar-imf2.35/index.fits','Sbc218-u4/Zsolar-imf2.35/index.fits',$
  'Sbc219-u4/Zsolar-imf2.35/index.fits']
leg=['Sbc201','Sbc202','Sbc203','Sbc204',$
     'Sbc205','Sbc206','Sbc217','Sbc218','Sbc219']
lummag,fi,/ps,band=0,leg='kuk'
lummag,fi,/ps,band=2,leg='kuk'
lummag,fi,/ps,band=4,leg='kuk'
lummag,fi,/ps,band=6,leg='kuk'
lummag,fi,/ps,band=9,leg='kuk'
lummag,fi,/ps,band=13,leg=leg
end


;-------------------------------------------------------

pro lowa1600,ps=ps,plot_only=plot_only

if keyword_set (ps) then begin
    !P.FONT = 0
    set_plot,'ps'
    device,/color,/encapsulated ,bits_per_pixel=8
    device, SET_CHARACTER_SIZE=[200,300]
end else begin
    if keyword_set (plot_only) then begin
        set_plot, 'x' 
        !P.FONT=-1
        device,pseudo_color=8,set_font='helvetica'
        device,set_character_size=[10,12]
    end else begin
        set_plot,'z'
        device,set_character_size=[9,12]
    end
end

a=mrdfits('Sbc201a-u4/0.2Zsolar-imf2.35/index.fits',1)
b=mrdfits('Sbc201a-u4/0.4Zsolar-imf2.35/index.fits',1)
c=mrdfits('Sbc201a-u4/Zsolar-imf2.35/index.fits',1)

; MHC99 A1600 vs beta
if !D.name eq 'PS' then device,file='a1600-beta-lowa.eps'
plot,[1],[1],$
  xtitle=textoidl('\Delta\beta'), $
                                ;xtitle=textoidl('\beta(M95)'), $
ytitle=textoidl("A_{1600}"), $
                                ;title=textoidl("A_{1600}-\beta  (Meurer, Heckman & Calzetti 1999)"),$
background = 1, color = 0,/nodata,$
  xrange=[-2.5,-1],yrange=[0,2],xst=1,yst=1

;oplot,a.betaM95-a.betam95_ns,-2.5*alog10(a.l_1600/a.l_1600_ns),ps=4,syms=0.5,col=3
;oplot,b.betaM95-b.betam95_ns,-2.5*alog10(b.l_1600/b.l_1600_ns),ps=4,syms=0.5,col=4
;oplot,c.betaM95-c.betam95_ns,-2.5*alog10(c.l_1600/c.l_1600_ns),ps=4,syms=0.5,col=2
oplot,a.betaM95,-2.5*alog10(a.l_1600/a.l_1600_ns),ps=4,syms=0.5,col=3
oplot,b.betaM95,-2.5*alog10(b.l_1600/b.l_1600_ns),ps=4,syms=0.5,col=4
oplot,c.betaM95,-2.5*alog10(c.l_1600/c.l_1600_ns),ps=4,syms=0.5,col=2
oplot,a.betaM95_ns,-2.5*alog10(a.l_1600/a.l_1600_ns),ps=1,syms=0.5,col=3
oplot,b.betaM95_ns,-2.5*alog10(b.l_1600/b.l_1600_ns),ps=1,syms=0.5,col=4
oplot,c.betaM95_ns,-2.5*alog10(c.l_1600/c.l_1600_ns),ps=1,syms=0.5,col=2

; the MHC99 fit
bfit= indgen(100)/100.*6.0-2.2
Afit= 0+8.9*(bfit+0*1.8-0.3)
oplot,bfit, Afit,color= 0
Afit= 0+3.0*(bfit+0*1.8-0.4)
oplot,bfit, Afit,color= 0,line=1
legend,[textoidl('0.2Z_{'+sunsymbol()+'}'),textoidl('0.4Z_{'+sunsymbol()+'}'), $
        textoidl('Z_{'+sunsymbol()+'}')], /bottom , $
  colors=[3,4,2],textcolor=0,psym=[4,4,4],linestyle=[0,0,0],/right,box=0       
save_plot, 'a1600-beta-lowa.png'
stop
end

;-------------------------------------------------------

pro SbcG

files=['Sbc201a-u4/Zsolar-imf2.35/index.fits',$
       'Sbc202-u4/Zsolar-imf2.35/index.fits',$
       'Sbc203-u4/Zsolar-imf2.35/index.fits',$
       'Sbc204-u4/Zsolar-imf2.35/index.fits',$
       'Sbc205-u4/Zsolar-imf2.35/index.fits',$
       'Sbc206-u4/Zsolar-imf2.35/index.fits',$
       'Sbc217-u4/Zsolar-imf2.35/index.fits',$
       'Sbc218-u4/Zsolar-imf2.35/index.fits',$
       'Sbc219-u4/Zsolar-imf2.35/index.fits',$
       'Sbcm201-u4/0.7Zsolar-imf2.35/index.fits',$
       'Sbcp201-u4/1.1Zsolar-imf2.35/index.fits',$
       'G0G0a-u1/0.3Zsolar-imf2.35/index.fits',$
       'G1G1a-u1/0.4Zsolar-imf2.35/index.fits',$
       'G2G2-u1/0.6Zsolar-imf2.35/index.fits',$
       'G3G3b-u1/Zsolar-imf2.35/index.fits']
plot_stuff,files,/ps
end
;-------------------------------------------------------

pro Sbconly

files=['Sbc201a-u4/Zsolar-imf2.35/index.fits',$
       'Sbc202-u4/Zsolar-imf2.35/index.fits',$
       'Sbc203-u4/Zsolar-imf2.35/index.fits',$
       'Sbc204-u4/Zsolar-imf2.35/index.fits',$
       'Sbc205-u4/Zsolar-imf2.35/index.fits',$
       'Sbc206-u4/Zsolar-imf2.35/index.fits',$
       'Sbc217-u4/Zsolar-imf2.35/index.fits',$
       'Sbc218-u4/Zsolar-imf2.35/index.fits',$
       'Sbc219-u4/Zsolar-imf2.35/index.fits',$
       'Sbcm201-u4/0.7Zsolar-imf2.35/index.fits',$
       'Sbcp201-u4/1.1Zsolar-imf2.35/index.fits']
plot_stuff,files,/ps
end

;-------------------------------------------------------

pro SbcGcZ

files=['Sbc201a-u4/Zsolar-imf2.35/index.fits',$
       'Sbc202-u4/Zsolar-imf2.35/index.fits',$
       'Sbc203-u4/Zsolar-imf2.35/index.fits',$
       'Sbc204-u4/Zsolar-imf2.35/index.fits',$
       'Sbc205-u4/Zsolar-imf2.35/index.fits',$
       'Sbc206-u4/Zsolar-imf2.35/index.fits',$
       'Sbc217-u4/Zsolar-imf2.35/index.fits',$
       'Sbc218-u4/Zsolar-imf2.35/index.fits',$
       'Sbc219-u4/Zsolar-imf2.35/index.fits',$
       'Sbcm201-u4/Zsolar-imf2.35/index.fits',$
       'Sbcp201-u4/Zsolar-imf2.35/index.fits',$
       'G0G0a-u1/Zsolar-imf2.35/index.fits',$
       'G1G1a-u1/Zsolar-imf2.35/index.fits',$
       'G2G2-u1/Zsolar-imf2.35/index.fits',$
       'G3G3b-u1/Zsolar-imf2.35/index.fits']
plot_stuff,files,/ps
end
;-------------------------------------------------------

pro set1

files=[$
;'Sbc201a-u4/set1/index.fits',$
;'Sbc201a-u50/set1/index.fits',$
       'Sbc201a10x-u4/set1/index.fits',$
       'Sbc202-u4/set1/index.fits',$
       'Sbc203-u4/set1/index.fits',$
       'Sbc204-u4/set1/index.fits',$
       'Sbc205-u4/set1/index.fits',$
      ; 'Sbc206-u4/set1/index.fits',$
       'Sbc217-u4/set1/index.fits',$
     ;  'Sbc219-u4/set1/index.fits',$
       ;'Sbc219-u50/set1/index.fits',$
       'Sbcm201-u4/set1/index.fits',$
       'Sbcp201-u4/set1/index.fits'$
       ;'G0G0a-u1/set1/index.fits',$
       ;'G1G1a-u1/set1/index.fits',$
       ;'G2G2-u1/set1/index.fits',$
       ;'G3G3b-u1/set1/index.fits'$
       ]
plot_stuff,files,/pl
end

;-------------------------------------------------------

; makes plots for code test section
pro wg96plot,ps=ps,plot_only=plot_only
graphics_init,ps=ps,plot_only=plot_only

fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7
sym=[4,5,6,7,1,2]
col=[2,3,4,5,6]
;goto,here
if keyword_set (ps) then device,file='wg96-fig8.eps'
plot,[1],[1], $
  xtitle=textoidl('\tau_{eff}'),$
  ytitle=textoidl('L_{scat}/L_{tot}'),$
;     title="Luminosity", $
background = 1, color = 0,/nodata, $
  xrange=[0,10.1],yr=[1e-3,1],xst=1,/yl

; read wg96 data
readcol,'wg96data-fig8.txt',wg96ff,wg96teff,wg96ln
for i=0,2 do begin
    oplot,wg96teff[7*i:7*i+6],wg96ln[7*i:7*i+6],ps=sym[i],col=col[i]
    oplot,wg96teff[7*i:7*i+6],wg96ln[7*i:7*i+6],ps=-0,col=col[i]
end

readcol,'f0.01.txt',th,ff,k,r,N,l,s,teff,ln
oplot,teff,ln,ps=sym[3],col=col[0]
readcol,'f0.05.txt',th,ff,k,r,N,l,s,teff,ln
oplot,teff,ln,ps=sym[4],col=col[1]
readcol,'f0.1.txt',th,ff,k,r,N,l,s,teff,ln
oplot,teff,ln,ps=sym[5],col=col[2]

xyouts,.2,.045,'k=0.001',col=fgc
xyouts,9.1,.0027,'k=1',col=fgc
legend,['Sunrise ff=0.01','WG96 ff=0.01','Sunrise ff=0.05','WG96 ff=0.05','Sunrise ff=0.10','WG96 ff=0.10'],$
  col=[col[0],col[0],col[1],col[1],col[2],col[2]],psym=[sym[3],sym[0],sym[3],sym[1],sym[4],sym[2]],/top,/right,box=0,textc=fgc
save_plot,'wg96-fig8.png'



if keyword_set (ps) then device,file='wg96-fig9.eps'
plot,[1],[1], $
  xtitle=textoidl('\tau_{eff}'),$
  ytitle=textoidl('L_{scat}/L_{tot}'),$
;     title="Luminosity", $
background = 1, color = 0,/nodata, $
  xrange=[0,10.1],yr=[1e-3,1],xst=1,/yl

; read wg96 data
readcol,'wg96data-fig9.txt',wg96n,wg96teff,wg96ln
for i=0,2 do begin
    oplot,wg96teff[7*i:7*i+6],wg96ln[7*i:7*i+6],ps=sym[i],col=col[i]
    oplot,wg96teff[7*i:7*i+6],wg96ln[7*i:7*i+6],ps=-0,col=col[i]
end

readcol,'N10.txt',th,ff,k,r,N,l,s,teff,ln
oplot,teff,ln,ps=sym[3],col=col[0]
readcol,'f0.1.txt',th,ff,k,r,N,l,s,teff,ln
oplot,teff,ln,ps=sym[4],col=col[1]
readcol,'N40.txt',th,ff,k,r,N,l,s,teff,ln
oplot,teff,ln,ps=sym[5],col=col[2]

xyouts,.5,.04,'k=0.001',col=fgc
xyouts,9.1,.0027,'k=1',col=fgc
legend,['Sunrise N=10','WG96 N=10','Sunrise N=20','WG96 N=20','Sunrise N=40','WG96 N=40'],$
  col=[col[0],col[0],col[1],col[1],col[2],col[2]],psym=[sym[3],sym[0],sym[3],sym[1],sym[4],sym[2]],/top,/right,box=0,textc=fgc
save_plot,'wg96-fig9.png'

stop

;; surface brightness plots
here:
; read wg96 data, 25 per case
readcol,'wg96data-fig15.txt',wg96r,wg96sn

if keyword_set (ps) then device,file='wg96-fig15-1.eps'
plot,[1],[1], $
     xtitle=textoidl('radius'),$
     ytitle=textoidl('log(\Sigma_{scat}/L_{tot})'),$
     background = 1, color = 0,/nodata, $
     xrange=[0,2100],yr=[0,5],xst=1,ys=1

img=0
readcol,'sb/k1.txt',row,r,sn
oplot,r*2020./50,alog10(sn)+16.48,col=col[img]
oplot,wg96r[25*img:25*img+24],wg96sn[25*img:25*img+24],ps=sym[img],col=col[img]

img=1
readcol,'sb/k0.316.txt',row,r,sn
oplot,r*2020./50,alog10(sn)+16.48,col=col[img]
oplot,wg96r[25*img:25*img+24],wg96sn[25*img:25*img+24],ps=sym[img],col=col[img]

img=2
readcol,'sb/k0.1.txt',row,r,sn
oplot,r*2020./50,alog10(sn)+16.48,col=col[img]
oplot,wg96r[25*img:25*img+24],wg96sn[25*img:25*img+24],ps=sym[img],col=col[img]

img=3
readcol,'sb/k0.0316.txt',row,r,sn
oplot,r*2020./50,alog10(sn)+16.48,col=col[img]
oplot,wg96r[25*img:25*img+24],wg96sn[25*img:25*img+24],ps=sym[img],col=col[img]

legend,['k=1.000','k=0.316','k=0.100','k=0.032'], $
       col=[col[0],col[1],col[2],col[3]],psym=[sym[0],sym[1],sym[2],sym[3]],/top,/right,box=0,textc=fgc
save_plot,'wg96-fig15-1.png'

if keyword_set (ps) then device,file='wg96-fig15-2.eps'
plot,[1],[1], $
  xtitle=textoidl('radius'),$
  ytitle=textoidl('log(\Sigma_{scat}/L_{tot})'),$
;     title="Luminosity", $
background = 1, color = 0,/nodata, $
  xrange=[-50,2100],yr=[0,5],xst=1,ys=1

img=4
readcol,'sb/k0.01.txt',row,r,sn
oplot,r*2020./50,alog10(sn)+16.48,col=col[img-4]
oplot,wg96r[25*img:25*img+24],wg96sn[25*img:25*img+24],ps=sym[img-4],col=col[img-4]

img=5
readcol,'sb/k0.00316.txt',row,r,sn
oplot,r*2020./50,alog10(sn)+16.48,col=col[img-4]
oplot,wg96r[25*img:25*img+24],wg96sn[25*img:25*img+24],ps=sym[img-4],col=col[img-4]

img=6
readcol,'sb/k0.001.txt',row,r,sn
oplot,r*2020./50,alog10(sn)+16.48,col=col[img-4]
oplot,wg96r[25*img:25*img+24],wg96sn[25*img:25*img+24],ps=sym[img-4],col=col[img-4]

legend,['k=0.010','k=0.003','k=0.001'], $
  col=[col[0],col[1],col[2]],psym=[sym[0],sym[1],sym[2]],/top,/right,box=0,textc=fgc
save_plot,'wg96-fig15-2.png'



stop
end
;-------------------------------------------------------


; makes plots for code test section
pro w77plot,ps=ps,plot_only=plot_only
graphics_init,ps=ps,plot_only=plot_only

fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7

pi=4*atan(1.)

goto,p3
; read data
readcol, 'w77t12_tau1.txt', x, y
; because I don't know how to split by empty lines I'll have to do
; this
start = x [0]
starts = where (x eq x [0])
y (where (y gt 2e4*5e-6)) = 0

readcol, 'witt77results.txt', wx, wy1,wy2,wy3,wy4
wr1 = indgen (12)
wr2 = indgen (16) + 12

if keyword_set (ps) then device,file='w77-1.eps'
plot,[1],[1], $
  xtitle=textoidl('r/arcsec'),$
  ytitle=textoidl('log (S/F)/sr^{-1}'),$
  background = 1, color = 0,/nodata, $
  xrange=[-1200,1200],yr=[1.5,4.5],xst=1,yst=1

oplot,x [starts [1]:starts [2] - 1], $
  alog10(y[starts [1]:starts [2] - 1]/5e-6), color =c1
oplot,x [starts [5]:starts [6] - 1], $
  alog10(y[starts [5]:starts [6] - 1]/5e-6), color =c2
oplot,wx[wr1],wy1[wr1]+ 0.3, color = c1,ps=1
oplot,wx[wr2],wy1[wr2]+ 0.3, color = c2,ps=2

readcol, 'w77t34_tau1.txt', x, y
; because I don't know how to split by empty lines I'll have to do
; this
start = x [0]
starts = where (x eq x [0])
y (where (y gt 2e4*5e-6)) = 0

oplot,x [starts [1]:starts [2] - 1], $
  alog10(y[starts [1]:starts [2] - 1]/5e-6), color =c3
oplot,x [starts [5]:starts [6] - 1], $
  alog10(y[starts [5]:starts [6] - 1]/5e-6), color =c4
oplot,wx[wr1],wy3[wr1]+ 0.3, color = c3,ps=4
oplot,wx[wr2],wy3[wr2]+ 0.3, color = c4,ps=5

legend,['g=0.01, 10deg','g=0.01, 50deg','g=0.70, 10deg','g=0.70, 50deg'],$
  col=[c1,c2,c3,c4],psym= [1,2,4,5],/top,/left,box=0,textc=fgc,line=0
save_plot,'w77-1.png'

;;; *** plot 2 ***
p2:
; read data
readcol, 'w77t12_tau5.txt', x, y
; because I don't know how to split by empty lines I'll have to do
; this
start = x [0]
starts = where (x eq x [0])
y (where (y gt 1e5*5e-6)) = 0

if keyword_set (ps) then device,file='w77-2.eps'
plot,[1],[1], $
  xtitle=textoidl('r/arcsec'),$
  ytitle=textoidl('log (S/F)/sr^{-1}'),$
  background = 1, color = 0,/nodata, $
  xrange=[-1200,1200],yr=[1.5,5],xst=1,yst=1

oplot,x [starts [1]:starts [2] - 1], $
  alog10(y[starts [1]:starts [2] - 1]/5e-6), color =c1
oplot,x [starts [5]:starts [6] - 1], $
  alog10(y[starts [5]:starts [6] - 1]/5e-6), color =c2
oplot,wx[wr1],wy2[wr1]+ 0.3, color = c1,ps=1
oplot,wx[wr2],wy2[wr2]+ 0.3, color = c2,ps=2

readcol, 'w77t34_tau5.txt', x, y
; because I don't know how to split by empty lines I'll have to do
; this
start = x [0]
starts = where (x eq x [0])
y (where (y gt 2e4*5e-6)) = 0

oplot,x [starts [1]:starts [2] - 1], $
  alog10(y[starts [1]:starts [2] - 1]/5e-6), color =c3
oplot,x [starts [5]:starts [6] - 1], $
  alog10(y[starts [5]:starts [6] - 1]/5e-6), color =c4
oplot,wx[wr1],wy4[wr1]+ 0.3, color = c3,ps=4
oplot,wx[wr2],wy4[wr2]+ 0.3, color = c4,ps=5

legend,['g=0.01, 10deg','g=0.01, 50deg','g=0.70, 10deg','g=0.70, 50deg'],$
  col=[c1,c2,c3,c4],psym= [1,2,4,5],/top,/left,box=0,textc=fgc,line=0
save_plot,'w77-2.png'

;;; *** plot 3 ***
p3:
; read data
readcol, 'w77III_g.7_tau3_a.6.txt', x, y
; because I don't know how to split by empty lines I'll have to do
; this
start = x [0]
starts = where (x eq x [0])
y (where (y gt 1e5*5e-6)) = 0

readcol, 'witt77IIIt7.txt', wx1, wy1, wx2, wy2, wx3, wy3, wx4, wy4

if keyword_set (ps) then device,file='w77-3.eps'
plot,[1],[1], $
  xtitle=textoidl('r/arcsec'),$
  ytitle=textoidl('log (S/F)/sr^{-1}'),$
  background = 1, color = 0,/nodata, $
  xrange=[0,1200],yr=[3,6],xst=1,yst=1
oplot,x [starts [1]:starts [2] - 1], $
  alog10(y[starts [1]:starts [2] - 1]/(5e-6*exp(-1.5/cos(10*pi/180)))),$
  color =c1
oplot,x [starts [2]:starts [3] - 1], $
  alog10(y[starts [2]:starts [3] - 1]/(5e-6*exp(-1.5/cos(20*pi/180))))+.2,$
  color =c2
oplot,x [starts [3]:starts [4] - 1], $
  alog10(y[starts [3]:starts [4] - 1]/(5e-6*exp(-1.5/cos(30*pi/180))))+.4,$
  color =c3
oplot,x [starts [4]:starts [5] - 1], $
  alog10(y[starts [4]:starts [5] - 1]/(5e-6*exp(-1.5/cos(40*pi/180))))+.6,$
  color =c4
oplot,wx1,wy1+ 0.3, color = c1,ps=1
oplot,wx2,wy2+ 0.3+.2, color = c2,ps=2
oplot,wx3,wy3+ 0.3+.4, color = c3,ps=4
oplot,wx4,wy4+ 0.3+.6, color = c4,ps=5

legend,['10deg','20deg','30deg','40deg'],$
  col=[c1,c2,c3,c4],psym= [1,2,4,5],/top,/right,box=0,textc=fgc,line=0
save_plot,'w77-3.png'


;;; *** plot 4 ***
p4:
; read data
readcol, 'w77III_g.1_tau3_a.6.txt', x1, y1
start = x1 [0]
starts = where (x1 eq x1 [0])
y1 (where (y1 gt 5e5*5e-6)) = 0
readcol, 'w77III_g.3_tau3_a.6.txt', x2, y2
y2 (where (y2 gt 5e5*5e-6)) = 0
readcol, 'w77III_g.5_tau3_a.6.txt', x3, y3
y3 (where (y3 gt 5e5*5e-6)) = 0
readcol, 'w77III_g.7_tau3_a.6.txt', x4, y4
y4 (where (y4 gt 5e5*5e-6)) = 0
readcol, 'w77III_g.9_tau3_a.6.txt', x5, y5
y5 (where (y5 gt 5e5*5e-6)) = 0
;readcol, 'w77III_g.9_tau3_a.6_nf1.txt', x6, y6
;y6 (where (y6 gt 1e5*5e-6)) = 0

readcol, 'witt77IIIt10.txt', wx1, wy1, wy2, wy3, wy4, wy5

if keyword_set (ps) then device,file='w77-4.eps'
plot,[1],[1], $
  xtitle=textoidl('r/arcsec'),$
  ytitle=textoidl('log (S/F)/sr^{-1}'),$
  background = 1, color = 0,/nodata, $
  xrange=[0,1200],yr=[2,6],xst=1,yst=1

oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y1[starts [1]:starts [2] - 1]/(5e-6*exp(-1.5/cos(10*pi/180))))+.3,$
  color =c1
oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y2[starts [1]:starts [2] - 1]/(5e-6*exp(-1.5/cos(10*pi/180)))),$
  color =c2
oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y3[starts [1]:starts [2] - 1]/(5e-6*exp(-1.5/cos(10*pi/180))))-.3,$
  color =c3
oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y4[starts [1]:starts [2] - 1]/(5e-6*exp(-1.5/cos(10*pi/180))))-.7,$
  color =c4
oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y5[starts [1]:starts [2] - 1]/(5e-6*exp(-1.5/cos(10*pi/180))))-1,$
  color =c5
;oplot,x1 [starts [1]:starts [2] - 1], $
;  alog10(y6[starts [1]:starts [2] - 1]/(5e-6*exp(-1.5/cos(10*pi/180))))-1,$
;  color =c5
oplot,wx1,wy1+ 0.3+.3, color = c1,ps=1
oplot,wx1,wy2+ 0.3, color = c2,ps=2
oplot,wx1,wy3+ 0.3-.3, color = c3,ps=4
oplot,wx1,wy4+ 0.3-.7, color = c4,ps=5
oplot,wx1,wy5+ 0.3-1, color = c5,ps=6

; examine whether the large w77 bins explain rapid rise towords center
;tx=x1 [starts [1]:starts [2] - 1]
;ty=y5[starts [1]:starts [2] - 1]/(5e-6*exp(-1.5/cos(10*pi/180)))
;ty=ty[where((tx gt 50.5) and (tx lt 257.5))]
;tx=tx[where((tx gt 50.5) and (tx lt 257.5))]
;oplot,[avg(tx)],[alog10(avg(ty))]-1,color=c1,psym=2

legend,['g=0.1','g=0.3','g=0.5','g=0.7','g=0.9'],$
  col=[c1,c2,c3,c4,c5],psym= [1,2,4,5,6],/top,/right,box=0,textc=fgc,line=0
save_plot,'w77-4.png'

;;; *** plot 5 ***
p5:
; read data
readcol, 'w77III_g.7_tau1_a.2.txt', x1, y1
start = x1 [0]
starts = where (x1 eq x1 [0])
y1 (where (y1 gt 1e5*5e-6)) = 0
readcol, 'w77III_g.7_tau1_a.4.txt', x2, y2
y2 (where (y2 gt 1e5*5e-6)) = 0
readcol, 'w77III_g.7_tau1_a.6.txt', x3, y3
y3 (where (y3 gt 1e5*5e-6)) = 0
readcol, 'w77III_g.7_tau1_a.8.txt', x4, y4
y4 (where (y4 gt 1e5*5e-6)) = 0
readcol, 'w77III_g.7_tau1_a1.txt', x5, y5
y5 (where (y5 gt 1e5*5e-6)) = 0
;readcol, 'w77III_g.7_tau1_a1_nf1.txt', x6, y6
;y6 (where (y6 gt 1e5*5e-6)) = 0

readcol, 'witt77IIIt11.txt', wx1, wy1, wy2, wy3, wy4, wy5

if keyword_set (ps) then device,file='w77-5.eps'
plot,[1],[1], $
  xtitle=textoidl('r/arcsec'),$
  ytitle=textoidl('log (S/F)/sr^{-1}'),$
  background = 1, color = 0,/nodata, $
  xrange=[0,1200],yr=[2, 6 ],xst=1,yst=1

oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y1[starts [1]:starts [2] - 1]/(5e-6*exp(-.5/cos(10*pi/180)))),$
  color =c1
oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y2[starts [1]:starts [2] - 1]/(5e-6*exp(-.5/cos(10*pi/180)))),$
  color =c2
oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y3[starts [1]:starts [2] - 1]/(5e-6*exp(-.5/cos(10*pi/180)))),$
  color =c3
oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y4[starts [1]:starts [2] - 1]/(5e-6*exp(-.5/cos(10*pi/180)))),$
  color =c4
oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y5[starts [1]:starts [2] - 1]/(5e-6*exp(-.5/cos(10*pi/180)))),$
  color =c5
;oplot,x1 [starts [1]:starts [2] - 1], $
;  alog10(y6[starts [1]:starts [2] - 1]/(5e-6*exp(-.5/cos(10*pi/180)))),$
;  color =c5
oplot,wx1,wy1+ 0.3, color = c1,ps=1
oplot,wx1,wy2+ 0.3, color = c2,ps=2
oplot,wx1,wy3+ 0.3, color = c3,ps=4
oplot,wx1,wy4+ 0.3, color = c4,ps=5
oplot,wx1,wy5+ 0.3, color = c5,ps=6


legend,['a=0.2','a=0.4','a=0.6','a=0.8','a=1.0'],$
  col=[c1,c2,c3,c4,c5],psym= [1,2,4,5,6],/top,/right,box=0,textc=fgc,line=0
save_plot,'w77-5.png'

;;; *** plot 6 ***
p6:
; read data
readcol, 'w77III_g.7_tau1_a.6.txt', x1, y1
start = x1 [0]
starts = where (x1 eq x1 [0])
y1 (where (y1 gt 1e5*5e-6)) = 0
readcol, 'w77III_g.7_tau2_a.6.txt', x2, y2
y2 (where (y2 gt 1e5*5e-6)) = 0
readcol, 'w77III_g.7_tau3_a.6.txt', x3, y3
y3 (where (y3 gt 1e5*5e-6)) = 0
readcol, 'w77III_g.7_tau4_a.6.txt', x4, y4
y4 (where (y4 gt 1e5*5e-6)) = 0
readcol, 'w77III_g.7_tau5_a.6.txt', x5, y5
y5 (where (y5 gt 1e5*5e-6)) = 0

readcol, 'witt77IIIt4.txt', wx1, wy1, wy2, wy3, wy4, wy5

if keyword_set (ps) then device,file='w77-6.eps'
plot,[1],[1], $
  xtitle=textoidl('r/arcsec'),$
  ytitle=textoidl('log (S/F)/sr^{-1}'),$
  background = 1, color = 0,/nodata, $
  xrange=[-1050,1000],yr=[2.5, 6 ],xst=1,yst=1

oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y1[starts [1]:starts [2] - 1]/(5e-6*exp(-.5/cos(10*pi/180)))),$
  color =c1
oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y2[starts [1]:starts [2] - 1]/(5e-6*exp(-1/cos(10*pi/180)))),$
  color =c2
oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y3[starts [1]:starts [2] - 1]/(5e-6*exp(-1.5/cos(10*pi/180)))),$
  color =c3
oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y4[starts [1]:starts [2] - 1]/(5e-6*exp(-2/cos(10*pi/180)))),$
  color =c4
oplot,x1 [starts [1]:starts [2] - 1], $
  alog10(y5[starts [1]:starts [2] - 1]/(5e-6*exp(-2.5/cos(10*pi/180)))),$
  color =c5
oplot,wx1,wy1+ 0.3, color = c1,ps=1
oplot,wx1,wy2+ 0.3, color = c2,ps=2
oplot,wx1,wy3+ 0.3, color = c3,ps=4
oplot,wx1,wy4+ 0.3, color = c4,ps=5
oplot,wx1,wy5+ 0.3, color = c5,ps=6

legend,[textoidl('\tau=1'),textoidl('\tau=2'),textoidl('\tau=3'),textoidl('\tau=4'),textoidl('\tau=5')],$
  col=[c1,c2,c3,c4,c5],psym= [1,2,4,5,6],/top,/right,box=0,textc=fgc,line=0
save_plot,'w77-6.png'



end

;-------------------------------------------------------


; makes plots for code test section
pro polyplot,ps=ps,plot_only=plot_only
graphics_init,ps=ps,plot_only=plot_only

fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7

pi=4*atan(1.)

;;; *** plot 4 ***
p4:
; read data
readcol, 'w77III_g.1_tau3_a.6.txt', x1, y1
start = x1 [0]
starts = where (x1 eq x1 [0])
readcol, 'w77III_g.3_tau3_a.6.txt', x2, y2
readcol, 'w77III_g.5_tau3_a.6.txt', x3, y3
readcol, 'w77III_g.7_tau3_a.6.txt', x4, y4
readcol, 'w77III_g.9_tau3_a.6.txt', x5, y5

if keyword_set (ps) then device,file='polytest-1.eps'
plot,[1],[1], $
  xtitle=textoidl('r/arcsec'),$
  ytitle=textoidl('log (S_{poly}/S_{mono})'),$
  background = 1, color = 0,/nodata, $
  xrange=[0,1200],yr=[-2,2],xst=1,yst=1

; read polychromatic runs
;readcol,'poly.k1.g.1.txt',x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13
readcol,'../../sunrise.realpoly/src/poly2.txt',x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13
poly1=[[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],[p11],[p12],[p13]]
readcol,'poly.k3.g.5.txt',x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13
poly2=[[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],[p11],[p12],[p13]]
readcol,'poly.k5.g.9.txt',x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13
poly3=[[p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],[p11],[p12],[p13]]

oplot,x,x*0+1,color=fgc
oplot,x,x*0+.5,color=fgc
oplot,x,x*0,color=fgc
oplot,x,x*0-.5,color=fgc
oplot,x,x*0-1.3,color=fgc

oplot,x,alog10(poly1[*,0]/y1[starts [1]:starts [2] - 1])+1,linest=0,color=c1
oplot,x,alog10(poly1[*,1]/y2[starts [1]:starts [2] - 1])+.5,linest=0,color=c1
oplot,x,alog10(poly1[*,2]/y3[starts [1]:starts [2] - 1]),linest=0,color=c1
oplot,x,alog10(poly1[*,3]/y4[starts [1]:starts [2] - 1])-.5,linest=0,color=c1
oplot,x,alog10(poly1[*,4]/y5[starts [1]:starts [2] - 1])-1.3,linest=0,color=c1

oplot,x,alog10(poly2[*,0]/y1[starts [1]:starts [2] - 1])+1,linest=0,color=c2
oplot,x,alog10(poly2[*,1]/y2[starts [1]:starts [2] - 1])+.5,linest=0,color=c2
oplot,x,alog10(poly2[*,2]/y3[starts [1]:starts [2] - 1]),linest=0,color=c2
oplot,x,alog10(poly2[*,3]/y4[starts [1]:starts [2] - 1])-.5,linest=0,color=c2
oplot,x,alog10(poly2[*,4]/y5[starts [1]:starts [2] - 1])-1.3,linest=0,color=c2

oplot,x,alog10(poly3[*,0]/y1[starts [1]:starts [2] - 1])+1,linest=0,color=c3
oplot,x,alog10(poly3[*,1]/y2[starts [1]:starts [2] - 1])+.5,linest=0,color=c3
oplot,x,alog10(poly3[*,2]/y3[starts [1]:starts [2] - 1]),linest=0,color=c3
oplot,x,alog10(poly3[*,3]/y4[starts [1]:starts [2] - 1])-.5,linest=0,color=c3
oplot,x,alog10(poly3[*,4]/y5[starts [1]:starts [2] - 1])-1.3,linest=0,color=c3

xyouts,100,1.1,textoidl('g=0.1'),col=fgc
xyouts,100,.6,textoidl('g=0.3'),col=fgc
xyouts,100,.1,textoidl('g=0.5'),col=fgc
xyouts,100,-.4,textoidl('g=0.7'),col=fgc
xyouts,100,-1.2,textoidl('g=0.9'),col=fgc
xyouts,50,1.7,textoidl('\tau=3!ca=0.6'),col=fgc,chars=1.3
legend,[textoidl('\kappa_{ref}=1, g_{ref}=0.1'),textoidl('\kappa_{ref}=3, g_{ref}=0.5'),textoidl('\kappa_{ref}=5, g_{ref}=0.9')],$
  col=[c1,c2,c3],/top,/right,box=0,textc=fgc,line=0
save_plot,'polytest-1.png'

;;; *** plot 5 ***
p5:
; read data
readcol, 'w77III_g.7_tau1_a.2.txt', x1, y1
start = x1 [0]
starts = where (x1 eq x1 [0])
readcol, 'w77III_g.7_tau1_a.4.txt', x2, y2
readcol, 'w77III_g.7_tau1_a.6.txt', x3, y3
readcol, 'w77III_g.7_tau1_a.8.txt', x4, y4
readcol, 'w77III_g.7_tau1_a1.txt', x5, y5

if keyword_set (ps) then device,file='polytest-2.eps'
plot,[1],[1], $
  xtitle=textoidl('r/arcsec'),$
  ytitle=textoidl('log (S_{poly}/S_{mono})'),$
  background = 1, color = 0,/nodata, $
  xrange=[0,1200],yr=[-2,2 ],xst=1,yst=1

oplot,x,x*0+1,color=fgc
oplot,x,x*0+.5,color=fgc
oplot,x,x*0,color=fgc
oplot,x,x*0-.5,color=fgc
oplot,x,x*0-1.3,color=fgc

oplot,x,alog10(poly1[*,5]/y1[starts [1]:starts [2] - 1])+1,linest=0,color=c1
oplot,x,alog10(poly1[*,6]/y2[starts [1]:starts [2] - 1])+.5,linest=0,color=c1
oplot,x,alog10(poly1[*,7]/y3[starts [1]:starts [2] - 1]),linest=0,color=c1
oplot,x,alog10(poly1[*,8]/y4[starts [1]:starts [2] - 1])-.5,linest=0,color=c1
oplot,x,alog10(poly1[*,9]/y5[starts [1]:starts [2] - 1])-1.3,linest=0,color=c1

oplot,x,alog10(poly2[*,5]/y1[starts [1]:starts [2] - 1])+1,linest=0,color=c2
oplot,x,alog10(poly2[*,6]/y2[starts [1]:starts [2] - 1])+.5,linest=0,color=c2
oplot,x,alog10(poly2[*,7]/y3[starts [1]:starts [2] - 1]),linest=0,color=c2
oplot,x,alog10(poly2[*,8]/y4[starts [1]:starts [2] - 1])-.5,linest=0,color=c2
oplot,x,alog10(poly2[*,9]/y5[starts [1]:starts [2] - 1])-1.3,linest=0,color=c2

oplot,x,alog10(poly3[*,5]/y1[starts [1]:starts [2] - 1])+1,linest=0,color=c3
oplot,x,alog10(poly3[*,6]/y2[starts [1]:starts [2] - 1])+.5,linest=0,color=c3
oplot,x,alog10(poly3[*,7]/y3[starts [1]:starts [2] - 1]),linest=0,color=c3
oplot,x,alog10(poly3[*,8]/y4[starts [1]:starts [2] - 1])-.5,linest=0,color=c3
oplot,x,alog10(poly3[*,9]/y5[starts [1]:starts [2] - 1])-1.3,linest=0,color=c3

xyouts,100,1.1,textoidl('a=0.2'),col=fgc
xyouts,100,.6,textoidl('a=0.4'),col=fgc
xyouts,100,.1,textoidl('a=0.6'),col=fgc
xyouts,100,-.4,textoidl('a=0.8'),col=fgc
xyouts,100,-1.2,textoidl('a=1.0'),col=fgc
xyouts,50,1.7,textoidl('\tau=1!cg=0.7'),col=fgc,chars=1.3
legend,[textoidl('\kappa_{ref}=1, g_{ref}=0.1'),textoidl('\kappa_{ref}=3, g_{ref}=0.5'),textoidl('\kappa_{ref}=5, g_{ref}=0.9')],$
  col=[c1,c2,c3],/top,/right,box=0,textc=fgc,line=0
save_plot,'polytest-2.png'

;;; *** plot 6 ***
p6:
; read data
readcol, 'w77III_g.7_tau1_a.6.txt', x1, y1
start = x1 [0]
starts = where (x1 eq x1 [0])
readcol, 'w77III_g.7_tau2_a.6.txt', x2, y2
readcol, 'w77III_g.7_tau3_a.6.txt', x3, y3
readcol, 'w77III_g.7_tau4_a.6.txt', x4, y4
readcol, 'w77III_g.7_tau5_a.6.txt', x5, y5

if keyword_set (ps) then device,file='polytest-3.eps'
plot,[1],[1], $
  xtitle=textoidl('r/arcsec'),$
  ytitle=textoidl('log (S_{poly}/S_{mono})'),$
  background = 1, color = 0,/nodata, $
  xrange=[-1050,1000],yr=[-2,2 ],xst=1,yst=1

oplot,x,x*0+1,color=fgc
oplot,x,x*0+.5,color=fgc
oplot,x,x*0,color=fgc
oplot,x,x*0-.5,color=fgc
oplot,x,x*0-1.3,color=fgc

oplot,x,alog10(poly1[*,7]/y1[starts [1]:starts [2] - 1])+1,linest=0,color=c1
oplot,x,alog10(poly1[*,10]/y2[starts [1]:starts [2] - 1])+.5,linest=0,color=c1
oplot,x,alog10(poly1[*,3]/y3[starts [1]:starts [2] - 1]),linest=0,color=c1
oplot,x,alog10(poly1[*,11]/y4[starts [1]:starts [2] - 1])-.5,linest=0,color=c1
oplot,x,alog10(poly1[*,12]/y5[starts [1]:starts [2] - 1])-1.3,linest=0,color=c1

oplot,x,alog10(poly2[*,7]/y1[starts [1]:starts [2] - 1])+1,linest=0,color=c2
oplot,x,alog10(poly2[*,10]/y2[starts [1]:starts [2] - 1])+.5,linest=0,color=c2
oplot,x,alog10(poly2[*,3]/y3[starts [1]:starts [2] - 1]),linest=0,color=c2
oplot,x,alog10(poly2[*,11]/y4[starts [1]:starts [2] - 1])-.5,linest=0,color=c2
oplot,x,alog10(poly2[*,12]/y5[starts [1]:starts [2] - 1])-1.3,linest=0,color=c2

oplot,x,alog10(poly3[*,7]/y1[starts [1]:starts [2] - 1])+1,linest=0,color=c3
oplot,x,alog10(poly3[*,10]/y2[starts [1]:starts [2] - 1])+.5,linest=0,color=c3
oplot,x,alog10(poly3[*,3]/y3[starts [1]:starts [2] - 1]),linest=0,color=c3
oplot,x,alog10(poly3[*,11]/y4[starts [1]:starts [2] - 1])-.5,linest=0,color=c3
oplot,x,alog10(poly3[*,12]/y5[starts [1]:starts [2] - 1])-1.3,linest=0,color=c3

xyouts,-900,1.1,textoidl('\tau=1'),col=fgc
xyouts,-900,.6,textoidl('\tau=2'),col=fgc
xyouts,-900,.1,textoidl('\tau=3'),col=fgc
xyouts,-900,-.4,textoidl('\tau=4'),col=fgc
xyouts,-900,-1.1,textoidl('\tau=5'),col=fgc
xyouts,-950,1.7,textoidl('a=0.6!cg=0.7'),col=fgc,chars=1.3
legend,[textoidl('\kappa_{ref}=1, g_{ref}=0.1'),textoidl('\kappa_{ref}=3, g_{ref}=0.5'),textoidl('\kappa_{ref}=5, g_{ref}=0.9')],$
  col=[c1,c2,c3],/top,/right,box=0,textc=fgc,line=0
save_plot,'polytest-3.png'



end

;-------------------------------------------------------


; makes plots for code test section
pro wh01plot,ps=ps,plot_only=plot_only
graphics_init,ps=ps,plot_only=plot_only

fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7

pi=4*atan(1.)

readcol, 'wh01a.txt', x, yn
readcol, 'wh01a_1scat.txt', x, y1
readcol, 'wh01a_2scat.txt', x, y2

readcol, 'wh01results.txt', wx, wy0,wy1,wy2,wy3, wz0,wz1,wz2,wz3

if keyword_set (ps) then device,file='wh01-1.eps'
plot,[1],[1], $
  xtitle=textoidl('\theta/deg'),$
  ytitle=textoidl('log (L_{\infty})'),$
  background = 1, color = 0,/nodata, $
  xrange=[0, 180 ],yr=[-4,- .8],xst=1,yst=1

oplot, 180-x*180/pi, alog10(y1), color =c1,ps=2
oplot, 180-x*180/pi, alog10(y2-y1), color =c2,ps=2
oplot, 180-x*180/pi, alog10(yn-y2), color =c3,ps=2
oplot,wx, alog10 (wy0+wy1), color = c1,ps=4
oplot,wx, alog10 (wy2), color = c2,ps=4
oplot,wx, alog10 (wy3), color = c3,ps=4

legend,['WH01', 'Sunrise','0,1s','2s','>2s'], $
  col=[fgc,fgc,c1,c2,c3],psym= [4, 2,1,1,1],/top,/left ,box=0,textc=fgc,line=0
save_plot,'wh01-1.png'

;***pencil beam test

readcol, 'wh01b.txt', x, yn
readcol, 'wh01b_1scat.txt', x1, y1
readcol, 'wh01b_2scat.txt', x2, y2

if keyword_set (ps) then device,file='wh01-2.eps'
plot,[1],[1], $
  xtitle=textoidl('\theta/deg'),$
  ytitle=textoidl('log (L_{\infty})'),$
  background = 1, color = 0,/nodata, $
  xrange=[-5, 185 ],yr=[-3.5,- .6],xst=1,yst=1

;oplot, 180-x*180/pi, alog10(y), color =c1,ps=2
oplot, 180-x*180/pi, alog10(y1), color =c1,ps=2
oplot, 180-x*180/pi, alog10(y2-y1), color =c2,ps=2
oplot, 180-x*180/pi, alog10(yn-y2), color =c3,ps=2
;oplot,wx, alog10 (wz0+wz1+wz2+wz3), color = c2,ps=1
oplot,wx, alog10 (wz0+wz1), color = c1,ps=4
oplot,wx, alog10 (wz2), color = c2,ps=4
oplot,wx, alog10 (wz3), color = c3,ps=4

legend,['WH01', 'Sunrise','0,1s','2s','>2s'], $
  col=[fgc,fgc,c1,c2,c3],psym= [4, 2,1,1,1],/top,/right ,box=0,textc=fgc,line=0
save_plot,'wh01-2.png'
stop


end
;-------------------------------------------------------


; makes plots for code test section
pro scattertestplot,ps=ps,plot_only=plot_only
graphics_init,ps=ps,plot_only=plot_only

fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7

; read data

if keyword_set (ps) then device,file='scatteringtestresults.eps'
plot,[1],[1], $
  xtitle=textoidl('cos(\theta)'),$
  ytitle=textoidl('F_s/F_d'),$
  background = 1, color = 0,/nodata, $
  xrange=[-1,1],yr=[1e-7,5e-6],xst=1,/yl

readcol,'g0d50k0.001.txt',theta,theo,mcrx,mcrxerr
oplot,cos(theta),theo,col=fgc
oploterror,cos(theta),mcrx,mcrxerr,errcol=c1,psym=3,typ=1,/nohat

readcol,'g0.5d50k0.001.txt',theta,theo,mcrx,mcrxerr
oplot,cos(theta),theo,col=fgc
oploterror,cos(theta),mcrx,mcrxerr,errcol=c2,psym=3,line=0,typ=1,/nohat

legend,['g=0','g=0.5'],$
  col=[c1,c2],/top,/left,box=0,textc=fgc,line=0
save_plot,'scatteringtestresults.png'



end
;-------------------------------------------------------


; plots dust extinction curves
pro dust,ps=ps,plot_only=plot_only
graphics_init,ps=ps,plot_only=plot_only

fgc=0
bgc=1
c2=2
c1=3
c3=4
c4=5
c5=6
c6=7

readcol,'~/dust_data/kext_albedo_WD_MW_3.1_60', lambda,a,g,cext,kabs
readcol,'~/dust_data/kext_albedo_WD_SMCbar_0', lambdaSMC,aSMC,gSMC,cextSMC,kabsSMC

bfit_lambda_lower = [1.268e-7,1.309e-7,1.342e-7,1.407e-7,1.562e-7,1.677e-7,1.760e-7,1.866e-7,1.930e-7,2.400e-7]
bfit_lambda_upper = [1.284e-7,1.316e-7,1.371e-7,1.515e-7,1.583e-7,1.740e-7,1.833e-7,1.890e-7,1.950e-7,2.580e-7]

if keyword_set (ps) then device,file='dustcurve.eps'
plot,lambda*1e-6,kabs/(1-a)*2.089e-10, $
  xtitle=textoidl('\lambda/m'),$
  ytitle=textoidl('k/(kpc2/M_{'+sunsymbol()+'})'),$
  background = 1, color = 0,/nodata, $
  xrange=[7e-8,4e-7],xst=1,/yl,/xl
oplot,lambda*1e-6,kabs/(1-a)*2.089e-10, col=c1
oplot,lambdaSMC*1e-6,kabsSMC/(1-aSMC)*2.089e-10, col=c2
for i=0,n_elements(bfit_lambda_lower)-1 do $
  oplot,[bfit_lambda_lower[i],bfit_lambda_upper[i]],[1.4e-5,1.4e-5],col=fgc
legend,['MW','SMC'],col=[c1,c2],textc=fgc,box=0,line=0,/right
save_plot,'dustcurve.png'
end

;-------------------------------------------------------

pro graphics_init,ps=ps,plot_only=plot_only
if keyword_set (ps) then begin
    !P.FONT = -1
    set_plot,'ps'
    device,/encaps,/color,bits_per_pixel=8, $
      preview=0,pre_depth=8,pre_x=3.5,pre_y=3,/inches
    device, SET_CHARACTER_SIZE=[300,350],xsize=14,ysize=12
    !X.THICK=3
    !Y.THICK=3
    !P.THICK=2
    !X.MARGIN=[8,3]
end else begin
    !X.THICK=1
    !Y.THICK=1
    !P.THICK=1
    !X.MARGIN=[8,3]
    if keyword_set (plot_only) then begin
        set_plot, 'x' 
        !P.FONT=-1
        device,true_color=24,decomposed=0,set_font='helvetica'
        device,set_character_size=[10,12]
    end else begin
        set_plot,'z'
        device,set_character_size=[9,12]
    end
end
!P.MULTI=0
colr=[0,255,255,0  ,0  ,255,215,  0,155,  0,155,235]
colg=[0,255,0  ,0  ,235,0  ,215,255,155,130,  0, 95]
colb=[0,255,0  ,255,0  ,255,  0,255,155,186,215,  0]
tvlct,colr,colg,colb

end

; -------------------------------------------------------

function find_HDU, file, name
i = 0
while 1 do begin
    if (strlen (name) lt 8) then begin
        for j=1, 8-strlen (name) do name = name + ' '
    end
    h = headfits (file, exten = i)
    n =fxpar (h, 'extname')
    ;;print, '"',n,'"'
    if (strcmp (n, name,/fold_case)) then return, i else i = i+ 1
end
end



