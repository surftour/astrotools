pro lg_mollview_proj2out, planmap, Tmax, Tmin, color_bar, dx, title_display, sunits, $
              coord_out, do_rot, eul_mat, planvec, vector_scale, $
              CHARSIZE=charsize, COLT=colt, CROP=crop, GIF = gif, GRATICULE = graticule, $
              HXSIZE = hxsize, NOBAR = nobar, NOLABELS = nolabels, NOPOSITION = noposition, $
              PREVIEW = preview, $
              PS = ps, PXSIZE=pxsize, PYSIZE=pysize, ROT=rot_ang, SUBTITLE = subtitle, $
              TITLEPLOT = titleplot, XPOS = xpos, YPOS = ypos, $
              POLARIZATION=polarization, $
              PNG = png, OUTLINE = outline, $
              PROJECTION=projection, MOLLWEIDE=mollweide, GNOMIC=gnomic, CARTESIAN=cartesian, $
              ORTHOGRAPHIC=orthographic, FLIP=flip, HALF_SKY=half_sky,COORD_IN=coord_in, IGRATICULE = igraticule
;===============================================================================
;+
;  PROJ2OUT
;  ouputs on X-screen or PS file or GIF or PNG file a gnomonic or
;  mollweide or cartesian map
;
;  IN:
;    planmap, Tmax, Tmin, color_bar, dx, title_display, sunits, 
;    coord_out, do_rot, eul_mat, planvec, vector_scale
;
;  KEYWORDS:
;     CHARSIZE=charsize, COLT=colt, CROP=crop, GIF = gif, GRATICULE = graticule, HXSIZE = hxsize, $
;              NOBAR = nobar, NOLABELS = nolabels, NOPOSITION = noposition, $
;              PREVIEW = preview, PS = ps, $
;              PXSIZE=pxsize, PYSIZE=pysize, ROT = rot_ang, SUBTITLE = subtitle, $
;              TITLEPLOT = titleplot, XPOS = xpos, YPOS = ypos, $
;              POLARIZATION=polarization,$
;              PNG=png, OUTLINE = outline,$
;              PROJECTION=projection, MOLLWEIDE=mollweide, $
;              GNOMIC=gnomic, CARTESIAN=cartesian,
;              ORTHOGRAPHIC=orthographic, $
;              FLIP=flip, HALF_SKY=half_sky,COORD_IN=coord_in
;
;   for more information, see Gnomview.pro Mollview.pro
;
;   March 1999, EH
;   Nov 2000, EH, plot polarisation
;   May 2002, EH, merge gnom2out and moll2out
;   Jun 2002, EH, added the cartesian projection facility (hacked from
;       G. Giardino's pol2out)
;   Aug 2002, EH, added the orthographic projection facility
;   Jul 2002, EH, changed vector field loop index to LONG
;-
;===============================================================================

identify_projection, projtype, projection=projection, mollweide=mollweide, gnomic=gnomic, cartesian=cartesian, orthographic=orthographic
do_gnom = 0
do_moll = 0
do_cart = 0
do_orth = 0
do_fullsky = 0 ; dummy, only matters for orthview
;-------------------------------------------------

@idl_default_previewer ; defines the paper size

xsize = (size(planmap))(1)
ysize = (size(planmap))(2)

if (projtype eq 2) then begin
;  ---- Gnomonic specific definitions for the plot ----
    routine    = 'GNOMVIEW'
    proj_small = 'gnomic'
    proj_big   = 'Gnomic'
    do_gnom    = 1

;     du_dv = 1.                  ; aspect ratio
    du_dv = xsize/float(ysize)                  ; aspect ratio
    fudge = 1.00                ; 
    xc = (xsize-1)/2. & delta_x = (xsize-1 - xc)
    yc = (ysize-1)/2. & delta_y = (ysize-1 - yc)
; u and v range of the map
    umin = - dx * xc * fudge & umax = dx * xc * fudge
    vmin = - dx * yc * fudge & vmax = dx * yc * fudge
; position of the rectangle in the final window
    w_xll = 0.00 & w_xur = 1.00 & w_dx = w_xur - w_xll
    w_yll = 0.10 & w_yur = 0.90 & w_dy = w_yur - w_yll
    w_dx_dy = w_dx / w_dy       ; 1.4
; color bar, position, dimension
    cbar_dx = 1./3.
    cbar_dy = 1./70.
    cbar_xll = (1. - cbar_dx)/2.
    cbar_xur = (1. + cbar_dx)/2.
    cbar_yur = w_yll - cbar_dy
    cbar_yll = cbar_yur - cbar_dy
; polarisation color ring, position, dimension
    cring_dx = 1./15.
    cring_dy = 1./15.
    cring_xll = .025
; cring_xur = .025 + cring_dx
; cring_yur = w_yll
; cring_yll = cring_yur - cring_dy
    cring_yll = .025
; location of astro. coordinate
    x_aspos = 0.5
    y_aspos = 0.04
; pol vector scale
    vscal_x = 0.05
    vscal_y = 0.01
; location of title and subtitle
    x_title = 0.5 & y_title = 0.95
    x_subtl = 0.5 & y_subtl = 0.915
; default X dimension of hardcopy (cm)
    hxsize_def = 15.
; offset along the long axis of the page
    ;yoffset = 2  ; Europe (A4)
    yoffset = 1                 ; US (letter)
endif

if (projtype eq 1) then begin
;  ---- Mollweide specific definitions for the plot ----
    routine    = 'MOLLVIEW'
    proj_big   = 'Mollweide'
    proj_small = 'mollweide'
    do_moll    = 1

    du_dv = 2.                  ; aspect ratio
    fudge = 1.02                ; spare some space around the Mollweide egg
    xc = 0.5*(xsize-1) & delta_x = (xsize-1 - xc)
    yc = 0.5*(ysize-1) & delta_y = (ysize-1 - yc)
; x and y range of egg
    umin = - du_dv * fudge & umax = du_dv * fudge
    vmin = - fudge         & vmax =         fudge
; position of the egg in the final window
    w_xll = 0.0 & w_xur = 1.0 & w_dx = w_xur - w_xll
    w_yll = 0.1 & w_yur = 0.9 & w_dy = w_yur - w_yll
    w_dx_dy = w_dx / w_dy       ; 1./.8
; color bar, position, dimension
    cbar_dx = 1./3.
    cbar_dy = 1./70.
    cbar_xll = (1. - cbar_dx)/2.
    cbar_xur = (1. + cbar_dx)/2.
    cbar_yur = w_yll - cbar_dy
    cbar_yll = cbar_yur - cbar_dy
; polarisation color ring, position, dimension
    cring_dx = 1./10.
    cring_dy = 1./10.
    cring_xll = .025
; cring_xur = .025 + cring_dx
; cring_yur = w_yll
; cring_yll = cring_yur - cring_dy
    cring_yll = .025
; pol vector scale
    vscal_x = 0.05
    vscal_y = 0.02
; location of title and subtitle
    x_title = 0.5 & y_title = 0.95
    x_subtl = 0.5 & y_subtl = 0.905
; default X dimension of hardcopy (cm)
    hxsize_def = 26.
; offset along the long axis of the page
    ;yoffset = 2  ; Europe (A4)
    yoffset = 1                 ; US (letter)
endif

if (projtype eq 4) then begin
;  ---- Orthographic specific definitions for the plot ----
    routine    = 'ORTHVIEW'
    proj_big   = 'Orthographic'
    proj_small = 'orthographic'
    do_orth    = 1
    
    if keyword_set(half_sky) then do_fullsky = 0 else do_fullsky = 1
    if (do_fullsky) then du_dv = 2. else du_dv = 1. ; aspect ratio

    fudge = 1.02                ; spare some space around the Mollweide egg
    xc = 0.5*(xsize-1) & delta_x = (xsize-1 - xc)
    yc = 0.5*(ysize-1) & delta_y = (ysize-1 - yc)
; x and y range of disc
    umin = - du_dv * fudge & umax = du_dv * fudge
    vmin = - fudge         & vmax =         fudge
; position of the disc in the final window
    w_xll = 0.0 & w_xur = 1.0 & w_dx = w_xur - w_xll
    w_yll = 0.1 & w_yur = 0.9 & w_dy = w_yur - w_yll
    w_dx_dy = w_dx / w_dy       ; 1./.8
; color bar, position, dimension
    cbar_dx = 1./3.
    cbar_dy = 1./70.
    cbar_xll = (1. - cbar_dx)/2.
    cbar_xur = (1. + cbar_dx)/2.
    cbar_yur = w_yll - cbar_dy
    cbar_yll = cbar_yur - cbar_dy
; polarisation color ring, position, dimension
    cring_dx = 1./10.
    cring_dy = 1./10.
    cring_xll = .025
; cring_xur = .025 + cring_dx
; cring_yur = w_yll
; cring_yll = cring_yur - cring_dy
    cring_yll = .025
; pol vector scale
    vscal_x = 0.05
    vscal_y = 0.02
; location of title and subtitle
    x_title = 0.5 & y_title = 0.95
    x_subtl = 0.5 & y_subtl = 0.905
; default X dimension of hardcopy (cm)
    hxsize_def = 26.
; offset along the long axis of the page
    ;yoffset = 2  ; Europe (A4)
    yoffset = 1                 ; US (letter)
endif

if (projtype eq 3) then begin
    ;------------ cartesion projection ----------------
    routine    = 'CARTVIEW'
    proj_small = 'cartesian'
    proj_big   = 'Cartesian'
    do_cart    = 1
    
;     du_dv = 1.                  ; aspect ratio
    du_dv = xsize/float(ysize)                  ; aspect ratio
    fudge = 1.00                ; 
    xc = (xsize-1)/2. & delta_x = (xsize-1 - xc)
    yc = (ysize-1)/2. & delta_y = (ysize-1 - yc)
; u and v range of the map
    umin = - dx * xc * fudge & umax = dx * xc * fudge
    vmin = - dx * yc * fudge & vmax = dx * yc * fudge
; position of the rectangle in the final window
    w_xll = 0.00 & w_xur = 1.00 & w_dx = w_xur - w_xll
    w_yll = 0.10 & w_yur = 0.90 & w_dy = w_yur - w_yll
    w_dx_dy = w_dx / w_dy       ; 1.4
; color bar, position, dimension
    cbar_dx = 1./3.
    cbar_dy = 1./70.
    cbar_xll = (1. - cbar_dx)/2.
    cbar_xur = (1. + cbar_dx)/2.
    cbar_yur = w_yll - cbar_dy
    cbar_yll = cbar_yur - cbar_dy
; polarisation color ring, position, dimension
    cring_dx = 1./15.
    cring_dy = 1./15.
    cring_xll = .025
    cring_yll = .025
; location of astro. coordinate
    x_aspos = 0.5
    y_aspos = 0.04
; pol vector scale
    vscal_x = 0.05
    vscal_y = 0.01
; location of title and subtitle
    x_title = 0.5 & y_title = 0.95
    x_subtl = 0.5 & y_subtl = 0.915
; default X dimension of hardcopy (cm)
    hxsize_def = 15.
; offset along the long axis of the page
    yoffset = (papersize eq 'a4') ? 2 : 1
    ;yoffset = 2  ; Europe (A4)
    ;yoffset = 1                 ; US (letter)
endif
;====================================================

if defined(colt) then ct=colt else ct = 33
if undefined(polarization) then polarization=0
do_polamplitude = (polarization eq 1)
do_poldirection = (polarization eq 2)
do_polvector    = (polarization eq 3)
if defined(charsize) then charsfactor = charsize else charsfactor = 1.0


; alter the color table
; -----------------------
print,'... computing the color table ...'
common colors, r,g,b,r_curr,g_curr,b_curr
if (do_poldirection) then begin
    LOADCT, 0 , /SILENT
    ncol = 256
    one = replicate(1.,ncol)
    tvlct,[0,0,0,findgen(ncol-3)]/(ncol-3)*720,one,one,/hsv ; hue is in degrees
    tvlct,r,g,b,/get
endif else begin
    print, '... new lg_mollview_proj2out ...'
    print, 'ct= ', ct
    print, 'setting device to ps'
    SET_plot, 'ps'
    LOADCT, ct , /SILENT
endelse
; set up some specific definitions
r(0) = 0   & g(0) = 0   & b(0) = 0 ; reserve for black
r(1) = 255 & g(1) = 255 & b(1) = 255 ; reserve for white
r(2) = 175 & g(2) = 175 & b(2) = 175 ; reserve for neutral grey
TVLCT,r,g,b

; ---------------------
; open the device
; ---------------------
my_background = !p.background
my_color = !p.color
print,'... here it is.'
titlewindow = proj_big+' projection : ' + title_display
back      = REPLICATE(BYTE(!P.BACKGROUND),xsize,(ysize*cbar_dy*w_dx_dy)>1)
do_gif = keyword_set(gif)
do_png = keyword_set(png)
do_image = do_gif or do_png
if (keyword_set(ps)) then begin
    if DEFINED(hxsize) then hxsize = (hxsize > 3) < 200 else hxsize = hxsize_def
    if ((size(ps))(1) ne 7) then file_ps = 'plot_'+proj_small+'.ps' else file_ps = ps
    old_device = !d.name
    SET_plot,'ps'
    do_portrait = 0
    do_landscape = 0
    DEVICE, FILE=file_ps, /COLOR, BITS = 8 ; opens the file that will receive the PostScript plot
    if (do_gnom or (do_orth and not do_fullsky)) then begin
        do_portrait = 1
        DEVICE, /PORTRAIT,  XSIZE=hxsize, YSIZE=hxsize/du_dv*w_dx_dy, XOFFSET=4, YOFFSET=2
    endif
    if (do_moll or (do_orth and do_fullsky)) then begin
        do_landscape = 1
        DEVICE, /LANDSCAPE, XSIZE=hxsize, YSIZE=hxsize/du_dv*w_dx_dy, XOFFSET=4, YOFFSET=hxsize+yoffset
    endif
    if (do_cart) then begin
        do_landscape = 1
;         DEVICE, /LANDSCAPE, XSIZE=hxsize, YSIZE=hxsize/du_dv*w_dx_dy, XOFFSET=4, YOFFSET=hxsize+yoffset
        DEVICE, /LANDSCAPE, XSIZE=hxsize, YSIZE=hxsize/du_dv*w_dx_dy, XOFFSET=0, YOFFSET=hxsize+yoffset
    endif
    TVLCT,r,g,b
    thick_dev = 2. ; device dependent thickness factor
endif else begin
    DEVICE, PSEUDO = 8
    to_patch = ((!d.n_colors GT 256) and do_image and not keyword_set(crop))
    if (to_patch) then device, decomp = 1 else device, decomp = 0
    if (UNDEFINED(xpos) or UNDEFINED(ypos)) then begin
        WINDOW, /FREE, XSIZE = xsize, YSIZE = ysize*w_dx_dy, TITLE = titlewindow
    endif else begin
        WINDOW, /FREE, XSIZE = xsize, YSIZE = ysize*w_dx_dy, TITLE = titlewindow, XP=xpos, YP=ypos
    endelse
    TVLCT,r,g,b
    thick_dev = 1. ; device dependent thickness factor
endelse

!p.background = my_background
!p.color = my_color
; -------------------------------------------------------------
; make the plot
; -------------------------------------------------------------

plot,/nodata,[umin,umax],[vmin,vmax],pos=[w_xll,w_yll,w_xur,w_yur],XSTYLE=5,YSTYLE=5
; ---------- projection independent ------------------
; map itself
TV, planmap,w_xll,w_yll,/normal,xsize=1.

hpxv11 = 0

; the polarisation color ring
if (do_poldirection) then begin
    npring = xsize*cring_dx
    one = replicate(1.,npring)
    yy  = one # (findgen(npring) - npring/2)
    xx  = transpose(yy)
    tt  = (xx^2 + yy^2) / float(npring)^2
    mask = byte(tt lt .25 and tt gt 0.05)
    if (hpxv11) then begin
        ; angles are measured from horizontal
        angle = atan(yy,xx)  * !radeg + 180. ; angle in deg in [0,360]
    endif else begin
        ; angles are measured from vertical
        angle = atan(xx,-yy)  * !radeg + 180. ; angle in deg in [0,360]
    endelse
    color_ring = (bytscl(angle,TOP=252) + 3B) * mask + byte((1-mask)*!P.BACKGROUND); in [3,255] in ring and !p.background outside ring
    TV,color_ring,cring_xll,cring_yll,/normal,xsize = npring/float(xsize)
endif

; polarisation vector field
if (do_polvector) then begin
    dyg = 10. & dxg = dyg
    pol_rescale = float(dyg)/ysize
    xg = (lindgen(xsize/dxg)+.5) #  replicate(dxg, ysize/dyg)
    yg = (lindgen(ysize/dyg)+.5) ## replicate(dyg, xsize/dxg)
    u = umin + xg * (umax-umin) / xsize
    v = vmin + yg * (vmax-vmin) / ysize
    for i=0L, n_elements(xg)-1 do begin
        norm = planvec[xg[i],yg[i],0] * pol_rescale * (vmax-vmin)
        angle = planvec[xg[i],yg[i],1]
        if (hpxv11) then begin
            ; angles are measured from horizontal
            if (norm gt 0) then plots, u[i]+norm*cos(angle)*[-.5,.5], v[i]+norm*sin(angle)*[-.5,.5]
        endif else begin
            ; angles are measured from vertical
            if (norm gt 0) then plots, u[i]-norm*sin(angle)*[-.5,.5], v[i]+norm*cos(angle)*[-.5,.5]
        endelse
    endfor
endif

;  the color bar
if (not keyword_set(nobar) and not do_poldirection) then begin
    color_bar = BYTE(CONGRID(color_bar,xsize*cbar_dx)) # REPLICATE(1.,(ysize*cbar_dy*w_dx_dy)>1)
    back(xsize*cbar_xll,0) = color_bar
    TV, back,0,cbar_yll,/normal,xsize = 1.
endif

;  the color bar labels
if (not keyword_set(nobar) and not keyword_set(nolabels) and not do_poldirection) then begin
    format = '(g10.2)'
    if ((Tmax - Tmin) ge 50 and MAX(ABS([Tmax,Tmin])) le 1.e5) then format='(i8)'
    if ((Tmax - Tmin) ge 5  and MAX(ABS([Tmax,Tmin])) le 1.e2) then format='(f6.1)'
    strmin = STRING(Tmin,format=format)
    strmax = STRING(Tmax,format=format)
    XYOUTS, cbar_xll, cbar_yll,'!6'+STRTRIM(strmin,2)+' ',ALIGN=1.,/normal, chars=1.3*charsfactor
    XYOUTS, cbar_xur, cbar_yll,'!6 '+STRTRIM(strmax,2)+' '+sunits,ALIGN=0.,/normal, chars=1.3*charsfactor
endif

; the polarisation vector scale
if (not keyword_set(nobar)  and do_polvector) then begin
    plots, vscal_x*[1,1], vscal_y+[0, 5*pol_rescale]*w_dy, /normal
    format = '(g10.2)'
    if ((5*vector_scale) lt 1.e3 and (5*vector_scale) ge 10) then format = '(f5.1)'
    if ((5*vector_scale) lt 10 and (5*vector_scale) gt 1.e-1) then format = '(f4.2)'
    xyouts, vscal_x, vscal_y + .5*(5*pol_rescale)*w_dy, '!6  '+strtrim(string(5*vector_scale,form=format),2)+' '+sunits,ALIGN=0.,/normal, chars=1.1*charsfactor
endif

;  the title
if (not keyword_set(titleplot)) then title= '!6'+title_display else title='!6'+titleplot
XYOUTS, x_title, y_title ,title, align=0.5, /normal, chars=1.6*charsfactor

;  the subtitle
if (keyword_set(subtitle)) then begin
    XYOUTS, x_subtl, y_subtl ,'!6 '+subtitle, align=0.5, /normal, chars=1.6*.75*charsfactor
endif

; ---------- projection dependent ------------------

if (do_gnom) then begin
;  astronomical position of central point
    if (not keyword_set(noposition)) then begin
        if (undefined(rot_ang)) then rot_ang = [0.,0.,0.] else rot_ang = ([rot_ang,0,0])(0:2)
        rot_0 = STRTRIM(STRING(rot_ang(0),form='(f6.1)'),2)
        rot_1 = STRTRIM(STRING(rot_ang(1),form='(f6.1)'),2)
        XYOUTS,x_aspos,y_aspos,'('+rot_0+', '+rot_1+') '+decode_coord(coord_out),/normal,align=0.5
    endif

; ; cross in the middle of plot
; plots,0,0,ps=1,syms=5,thick=4
; phi = findgen(31)/30.*2*!pi
; x_circle = cos(phi)
; y_circle = sin(phi)
; radius = tan(1.*!dtor/2.) ; radius = fwhm/2
; xyouts,0.7*umax,-0.8*vmax,'100 GHz'
; oplot,0.92*umax+radius*x_circle,-0.8*vmax+radius*y_circle,thick=3
; radius = tan(1./1.5*!dtor/2.)
; xyouts,0.7*umax,-0.9*vmax,'150 GHz'
; oplot,0.92*umax+radius*x_circle,-0.9*vmax+radius*y_circle,thick=3

endif

grattwice=0
;  the graticule in output astrophysical coordinates
if (KEYWORD_SET(graticule)) then begin
    grattwice =1
    oplot_graticule, graticule, eul_mat, projection=proj_small, flip = flip, thick = 1.*thick_dev, color = !p.color, half_sky=half_sky, linestyle=0
endif 

;  the graticule in input coordinates
if (KEYWORD_SET(igraticule)) then begin
    lines_ig = 2*grattwice ; either 0 or 2
    oplot_graticule, igraticule, eul_mat, projection=proj_small, flip = flip, thick = 1.*thick_dev, color = !p.color, half_sky=half_sky, linestyle=lines_ig, coordsys=[coord_in,coord_out]
endif 

; outlines on the map
if (keyword_set(outline)) then begin
    outline_coord2uv, outline, coord_out, eul_mat, projection=proj_small, flip = flip, /show, thick = 3.*thick_dev, half_sky=half_sky
endif


; -----------------------------------------------
;       output the PS/GIF/PNG
; -----------------------------------------------

;  gif/png output
if do_image then begin
    if (do_gif) then begin
        if (DATATYPE(gif) ne 'STR') then file_image = 'plot_'+proj_small+'.gif' else file_image = gif
    endif else begin
        if (DATATYPE(png) ne 'STR') then file_image = 'plot_'+proj_small+'.png' else file_image = png
    endelse        
    if keyword_set(crop) then begin
        if do_gif then write_gif,file_image,planmap,r,g,b
        if do_png then write_png,file_image,planmap,r,g,b
    endif else begin
        if do_gif then write_gif,file_image,tvrd(),r,g,b
        if do_png then write_png,file_image,tvrd(),r,g,b
        if (to_patch) then begin 
            device,decomp=0 ; put back colors on X window
            tv,tvrd()
        endif
    endelse
    print,'IMAGE file is in '+file_image
    if (keyword_set(preview)) then begin
        test_preview, /crash
        if do_gif then preview_file, file_image, /gif
        if do_png then preview_file, file_image, /png
    endif
endif

if (keyword_set(ps)) then begin
    device,/close
    set_plot,old_device
    print,'PS file is in '+file_ps
    if (keyword_set(preview)) then begin
        test_preview, /crash
        preview_file, file_ps, /ps, landscape=do_landscape
    endif
endif

return
end
