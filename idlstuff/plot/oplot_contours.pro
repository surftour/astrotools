;---------------------------------------------------------------
;
;
;
;---------------------------------------------------------------
pro oplot_contours, ContourMap, x0, y0, x1, y1, $
			hmsb=hmsb, $
                        xlen=xlen, $
			pixels=pixels, $
			clr=clr, $
                        fitellipse=fitellipse, $
			fitdiskyboxy=fitdiskyboxy, $
                        pa=pa, ellip=ellip, $
			fit_a=fit_a, fit_b=fit_b, $
			fit_ellip=fit_ellip, a4diva=a4diva


        Pic= ContourMap
        if keyword_set(hmsb) then HalfMassSB= hmsb else HalfMassSB= 125

	if not keyword_set(pixels) then pixels= (size(Pic))[1]

	if not keyword_set(clr) then clr= 0

		; do some slight processing
		;    i.e., flip color scale, and trim
		; --------------------------
                idx=where(Pic eq 1)
                Pic= 256-Pic
                if idx(0) ne -1 then Pic(idx)= 1
                print, "before: min/max= ",min(Pic), max(Pic)
                print, "n > 240= ", n_elements(where(Pic gt 240))
                print, "n < 20= ", n_elements(where(Pic lt 20))
                ; smooth?
                ;width= 10      ; try other things?
                width= 2
                newMap= smooth(Pic, width, /edge_truncate)
                Pic= newMap
                print, "after: min/max= ",min(Pic), max(Pic)


		; now, determine contour levels to show
		; ---------------------------------------
                ;levels = 8
                levels = 14
                ;levels = 16
                ;step = (Max(Pic) - Min(Pic)) / levels
                ;userLevels = IndGen(levels) * step + Min(Pic)
                step = (256)/levels
                userLevels = IndGen(levels) * step
                ;userLevels = userLevels[2:7]
                ;userLevels = userLevels[10:15]
                print, "userLevels= ",userLevels
                print, '256-HalfMassSB= ', 256-HalfMassSB
		level_to_fit= 0
                ; this is 256 minus because we flip Pic around above
		;goto, skip_half_cnt_insertion
		zeroidx= 0
		;zeroidx= 2
                for i=zeroidx,n_elements(userLevels)-1 do begin
                   if userLevels[i] gt (256-HalfMassSB) then begin
                        if i eq zeroidx then nuL= userLevels[0]
                        nuL=[nuL,256-HalfMassSB,userLevels[i:levels-1]]
                        level_to_fit= i-zeroidx
                        break
                   endif else begin
                        if i eq zeroidx then nuL= [userLevels[i]] else nuL=[nuL,userLevels[i]]
                   endelse
                endfor
skip_half_cnt_insertion:
                print, "new userLevels= ",nuL
                userLevels= nuL
                contour_to_fit= -1


                ; load contours into variables (this doesn't print)
		; -----------------------------
                contour, Pic, path_xy=xy_contours, path_info=info, levels=userLevels, $
                                min_value=2, /PATH_DATA_COORDS



fload_newcolortable, 4
		; show the contours
		; ----------------------
                for i=0, (n_elements(info)-1) do begin
                        cnt_index= [indgen(info(i).N),0]
                        xycnt= xy_contours(*,info(i).offset+cnt_index)
                        n_level= info(i).level
print, "i= ", i, "  level= ",n_level, "       N= ", info(i).N, "    userLevels[i]= ", userLevels[n_level]
                        xycnt[0,*]= xycnt[0,*]*(x1-x0)/pixels + x0
                        xycnt[1,*]= xycnt[1,*]*(y1-y0)/pixels + y0
                        idx= where(xycnt[0,*] eq x0)
                        if idx(0) eq -1 then begin
                           ;if n_level eq 0 then plots, xycnt, /normal, color= 0, thick=2.0
                           ;if n_level eq 1 then plots, xycnt, /normal, color= 0, thick=2.0
                           ;if n_level eq 2 then plots, xycnt, /normal, color= 0, thick=2.0
                           ;if n_level eq 3 then plots, xycnt, /normal, color= 0, thick=2.0
                           ;if n_level eq 4 then plots, xycnt, /normal, color= 130
                           ;plots, xy_contours(*,info(i).offset+cnt_index), /normal, color= 150
                           ;plots, xycnt, /normal, color= 0, thick=2.0
                           ;plots, xycnt, /normal, color= 1, thick=2.0
                           plots, xycnt, /normal, color= clr, thick=1.0
                        endif

			; if 7 order between 10^22 cm^-2 and 10^15 and 14 levels then
			; level 13 is 10^21
			; level 11 is 10^20, etc.
			if n_level eq 13 then begin
                           plots, xycnt, /normal, color= 150, thick=4.0, linestyle= 0
			endif
			if n_level eq 11 then begin
                           plots, xycnt, /normal, color= 150, thick=3.0, linestyle= 2
			endif
			if n_level eq 9 then begin
                           plots, xycnt, /normal, color= 150, thick=2.0, linestyle= 1
			endif

                        ;if info(i).N gt 200 then contour_to_fit= i
                        if info(i).level eq level_to_fit and contour_to_fit eq -1 then begin
				contour_to_fit= i
				plots, xycnt, /normal, color= 1, thick=3.0, linestyle= 1
			endif
                endfor



        ; -------------------------
        ;
        ;  Fit Ellipse to Contour
        ;
        ;  Use the MPFit routine
        ; to fit and ellipse to one
        ; of the above contours.
        ; -------------------------
        ;fitellipse= 1
        ;fitellipse= 0
        if keyword_set(fitellipse) then begin

           ; --------------------------
           if contour_to_fit eq -1 then i= 2 else i= contour_to_fit
           ;if contour_to_fit eq -1 then i= 3 else i= contour_to_fit
           print, "fit i,N,levle= ",i,info(i).N,info(i).level
           cnt_index= [indgen(info(i).N),0]
           xy= xy_contours(*,info(i).offset+cnt_index)

           x=fltarr((size(xy))[2])
           y=fltarr((size(xy))[2])

           x(*)= xy[0,*]
           y(*)= xy[1,*]

           ; x,y are in pixel units at this point

           ; here's where we actually fit to ellipse
           params=mpfitellipse(x,y, /tilt)

           ; a and b are in normalized coordinates
           if params(0) GT params(1) then begin
                a= params(0)
                b= params(1)
                pa_rad= params(4)
           endif else begin
                a= params(1)
                b= params(0)
                ;if params(4) gt 0.0 then pa_rad= !PI*0.5 + params(4)
                ;if params(4) lt 0.0 then pa_rad= -0.5*!PI + params(4)
		; i think i'm fixing some rotation problems by commenting
		; out the above and simplifying to this to the straightforward
		; conversion below and do the fancy things elsewhere :
		pa_rad= params(4)
           endelse
           ellipticity= 1-b/a
	   ellip= ellipticity
           print, " xxxxxxxxxxxxxxxx "
           print, "  Fit Results "
           print, "Semiaxes:   a=", a, "   b=",b
           print, "ellipticity(1-b/a)= ", ellipticity

           ; convert to physics coordinates: normal= (2.0*xlen)/(x1-x0)
           ; CAREFUL: we need the image to be square to make this easy
           ; conversion
           ;a= a*2.0*xlen/(x1-x0)
           ;b= b*2.0*xlen/(x1-x0)   -> old methods

           ; new methods has x,y in pixel coordinates
           a= a*2.*xlen/pixels
           b= b*2.*xlen/pixels

	   print, "              xlen =", xlen
           print, "      (phy units) a=", a,"   b=",b
           print, " params(4) (in deg)=", params(4)*180.0/!PI
           print, "       phi (in deg)=", pa_rad*180.0/!PI
           print, "             params=", params
           print, " xxxxxxxxxxxxxxxx "

;pa_rad= 110.0
;print, "MANUALLY fixing pa_rad to ", pa_rad
;pa_rad= 110.0*!PI/180.0

           ; replot the contour we're fitting to.
;          plots, x, y, psym=-3, color=getcolor('red'), /normal

           ; plot the ellipse 
           ; -----------------
           phi = dindgen(101)*2D*!dpi/100
           ; data coord
           ;e_x = 2.*xlen/pixels * params(0)*cos(phi)
           ;e_y = 2.*xlen/pixels * params(1)*sin(phi)
           ; normal coord
           e_x = (x1-x0)*1./pixels * params(0)*cos(phi)
           e_y = (y1-y0)*1./pixels * params(1)*sin(phi)

           ; rotate to tilted frame
           ; ------------------------
           if params(4) NE 0 then begin
                e_x_prime = e_x*cos(params(4)) + e_y*sin(params(4))
                e_y_prime = -e_x*sin(params(4)) + e_y*cos(params(4))
                e_x = e_x_prime
                e_y = e_y_prime
                pa= pa_rad*180.0/!PI
           endif

	   fload_newcolortable, 4

           ; relocate to center
           ; ------------------
           ; data coord norm.
           ;e_x = 2.*xlen/pixels * params(2) + e_x
           ;e_y = 2.*xlen/pixels * params(3) + e_y
           ;oplot, e_x, e_y, color= 100, thick= 4.0
           ; normal coordinates
           e_x = x0 + (x1-x0)*1./pixels * params(2) + e_x
           e_y = y0 + (y1-y0)*1./pixels * params(3) + e_y
           ;plots, e_x, e_y, color= 150, thick= 4.0, /normal


           ; labels, if we want 'em
           ;ecclbl= strcompress(string(ellipticity),/remove_all)
           ;ecclbl= '!7e!6='+strmid(ecclbl,0,4)
           ;xyouts, x0+0.015, y1-0.13, ecclbl, /normal, size=1.4, color= 0, charthick= 3.0
           ;a= a/0.7
           ;albl= strcompress(string(a),/remove_all)
           ;;albl= '!6R!D!8a!6!N='+strmid(albl,0,4)
           ;albl= '!8a!6='+strmid(albl,0,4)
           ;xyouts, x0+0.015, y1-0.07, albl, /normal, size=1.4, color= 0, charthick= 3.0

        endif



        ; -------------------------
        ;
        ;  Fit Deviations from ellipse
        ;
        ;  Use the C routine
        ; -------------------------
	if keyword_set(fitdiskyboxy) then begin

            Thresh= 256-HalfMassSB
            SurfaceBrightness= Pic
            Scale= 2.0*xlen          ; allows the conversion to physical units

            Xfit=fltarr(360)
            YFit=fltarr(360)

            Xcon=fltarr(360)
            Ycon=fltarr(360)

            AFit=0.0
            BFit=0.0
            X0Fit=0.0
            Y0Fit=0.0
            PhiFit=0.0
            a4Fit=0.0

	    PIXELS= long((size(Pic))[1])

            S = CALL_EXTERNAL('/n/home/tcox/C-Routines_for_IDL/ComputeIsoPhotFitIDL/isophotfit.so', $
                      'fit_isophot', $
                      PIXELS, $
                      SurfaceBrightness,$
                      Thresh,$
                      Scale, $
                      Xfit, Yfit,$
                      PhiFit, AFit, BFit, X0Fit, Y0Fit, a4Fit, $
                      Xcon, Ycon)

           ; disky/boxy
	   print, " "
	   fit_a= max([AFit,BFit])
	   print, "fit_a= ", fit_a
	   fit_b= min([AFit,BFit])
	   print, "fit_b= ", fit_b
	   fit_ellip= 1-(fit_b/fit_a)
	   print, "fit_ellip= ", fit_ellip
           a4diva= a4Fit/AFit
           print, "a4/a (a4diva)= ", a4diva
	   pa= PhiFit * 180.0 / !PI
	   if AFit gt BFit then pa= pa + 90.0     ; the fit doesn't care which axis it is, just the closest one
	   print, "pa= ", pa
	   print, " "


                   ; for some reason we need axes to be reset

                   !p.position=[x0,y0,x1,y1]
                   !p.ticklen=0.03

                   ; creates axes and plot style
                   ; ---------------------------
                   plot,[0],[0], psym=3,  $
                      xrange=[-xlen,xlen], $
                      yrange=[-xlen,xlen], $
                      /NOERASE, $
                      color= 0, $
                      xcharsize=0.01, $
                      ycharsize=0.01, $
                      xthick=4.0, $
                      ythick=4.0, $
                      charthick=3.0, $
                      ;xticks=10, $
                      ;yticks=10, $
                      xtickformat='(a1)', ytickformat='(a1)', $
                      xstyle=1, ystyle=1, $
                      ;xtitle=xtit, $
                      ;ytitle=ytit, $ ;/normal, $
                      /nodata

           oplot, Xfit, Yfit, color= 150, thick= 4.0

           ; labels, if we want 'em
           ecclbl= strcompress(string(fit_ellip),/remove_all)
           ecclbl= '!7e!6='+strmid(ecclbl,0,4)
           xyouts, x0+0.015, y1-0.13, ecclbl, /normal, size=1.4, color= 0, charthick= 3.0
           ;a= fit_a/0.7
           a= fit_a
           albl= strcompress(string(a),/remove_all)
           ;albl= '!6R!D!8a!6!N='+strmid(albl,0,4)
           albl= 'a='+strmid(albl,0,4)+' kpc/h'
           xyouts, x0+0.015, y1-0.07, albl, /normal, size=1.4, color= 0, charthick= 3.0
           a4= a4diva
           albl= strcompress(string(a4),/remove_all)
           ;albl= '!6R!D!8a!6!N='+strmid(albl,0,4)
	   digs= 6
	   if a4 ge 0.0 then digs= 5
           albl= 'a4/a='+strmid(albl,0,digs)
           xyouts, x0+0.015, y1-0.19, albl, /normal, size=1.4, color= 0, charthick= 3.0



	endif





end












; =================================================================
; =================================================================

