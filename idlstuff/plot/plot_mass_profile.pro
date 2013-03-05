;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     Surface Density Plots
;   -------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------


; -------------------------------
;
;  Overplot devac profile for one
;  galaxy (and one snap).
;
; -------------------------------
pro add_one_profile, frun, snapnum, component, xmin, xmax, bins, linecolor=linecolor, msg=msg, $
					x_is_devac=x_is_devac, x_is_log=x_is_log, $
                                        fitdevac=fitdevac, $ 
                                        fitsersic=fitsersic, $ 
                                        includeerr=includeerr, $
                                        frommanyproj=frommanyproj


if not keyword_set(component) then component= "allstars"


if component eq "all" then begin


	process_and_plot_sd, frun, snapnum, "allstars", xmin, xmax, bins, linecolor=linecolor, $
                                        x_is_devac=x_is_devac, x_is_log=x_is_log, $
                                        fitdevac=fitdevac, $
                                        fitsersic=fitsersic, $ 
                                        includeerr=includeerr, $
                                        frommanyproj=frommanyproj, $
                                        h= 0.7

	x0= 0.70
	y0= 0.90


        ; gas profile
        ; --------------
        process_and_plot_sd, "loaded", 0, "gas", xmin, xmax, bins, linecolor=50, $
                                                x_is_devac=x_is_devac, $
                                                x_is_log=x_is_log

        if fload_npart(0) GT 0 then xyouts, x0, y0-0.04, 'gas', /normal, charthick= 3.0, size= 1.2, color= 50


        ; disk profile
        ; --------------
        process_and_plot_sd, "loaded", 0, "disk", xmin, xmax, bins, linecolor=100, $
                                                x_is_devac=x_is_devac, $
                                                x_is_log=x_is_log

        if fload_npart(2) GT 0 then xyouts, x0, y0-0.08, 'disk', /normal, charthick= 3.0, size= 1.2, color= 100


        ; bulge profile
        ; --------------
        process_and_plot_sd, "loaded", 0, "bulge", xmin, xmax, bins, linecolor=200, $
                                                x_is_devac=x_is_devac, $
                                                x_is_log=x_is_log

        if fload_npart(3) GT 0 then xyouts, x0, y0-0.12, 'bulge', /normal, charthick= 3.0, size= 1.2, color= 200


        ; new stars profile
        ; -------------------
        process_and_plot_sd, "loaded", 0, "newstars", xmin, xmax, bins, linecolor=150, $
                                                x_is_devac=x_is_devac, $
                                                x_is_log=x_is_log

        if fload_npart(4) GT 0 then xyouts, x0, y0-0.16, 'New Stars', /normal, charthick= 3.0, size= 1.2, color= 150


endif else begin


	process_and_plot_sd, frun, snapnum, component, xmin, xmax, bins, linecolor=linecolor, $
                                        x_is_devac=x_is_devac, x_is_log=x_is_log, $
                                        fitdevac=fitdevac, $
                                        fitsersic=fitsersic, $ 
                                        includeerr=includeerr, $
                                        frommanyproj=frommanyproj, $
                                        h= 0.7


endelse




end




;------------------------------------------------------------------------


pro plot_mass_profile, runinfo, $
			bycomp=bycomp, $
			smoothlen=smoothlen, $
			filename=filename, $
			fitsersic=fitsersic, $
			x_is_devac=x_is_devac, $
			x_is_log=x_is_log ,$
			ctab=ctab

if not keyword_set(runinfo) then begin
   print, "  "
   print, "plot_mass_profile, runinfo, smoothlen=smoothlen, filename=filename, $"
   print, "                /x_is_devac, /x_is_log  (default is linear)"
   print, "  "
   return
endif

if not keyword_set(smoothlen) then smoothlen= 0.1
if not keyword_set(filename) then filename='devac.eps'
if not keyword_set(ctab) then ctab= 4


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=ctab



; -----------------------------
; set up constants
; -----------------------------

bins = 40

yaxistitle= "!6Log !4R!6(R) !N(M!D!9n!6!N pc!E-2!N)"


; default is linear x axis (devac and log are set below)
;
xaxistitle="!6 R (kpc)"
xmax = 20.0
xmin = 0.0
ymax = 4.2
ymin = -1.0

if keyword_set(x_is_devac) then begin
	xaxistitle="!6[R (kpc)]!E1/4!N"
	xmax = 2.4
	;xmax = 4.2
	xmin = 0.5
	;ymax = 6.0
	ymax = 4.5
	;ymin = -1.0
	ymin = 0.8
	x_is_log= 0
endif

if keyword_set(x_is_log) then begin
	xaxistitle="!6Log R (kpc)"
	xmax = 1.6
	xmin = -1.0
	ymax = 4.4
	ymin = 0.8
	x_is_devac= 0
endif




; ------------------
; Plot this up
; ------------------
!p.position= [0.2, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, /nodata, $
        ;/ylog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle




; -----------------------------------------------



n_runs= (size(runinfo))[1]

for i=0, n_runs-1 do begin

        print, " "
        print, "i= ", i
        print, " "

	frun= runinfo[i].frun
	snapnum= runinfo[i].snapnum
	component= runinfo[i].component
	linecolor= runinfo[i].linecolor


	add_one_profile, frun, snapnum, component, xmin, xmax, bins, linecolor=linecolor, $
						x_is_devac=x_is_devac, $
						includeerr=includerr, $
						x_is_log=x_is_log, $
						fitsersic=fitsersic, $
						frommanyproj=frommanyproj

endfor


; -----------------------------------------------


	do_special= 0
	;do_special= 1
	if do_special eq 1 then begin

	        process_and_plot_sd, "/raid4/tcox/minor/iso_Sb", 0, "allstars", xmin, xmax, bins, linecolor=50
	        ;process_and_plot_sd, x, y, z, m, xmin, xmax, bins, linecolor=50
	        ;process_and_plot_sd, x, y, z, m, xmin, xmax, bins, linecolor=50, $
                ;                                x_is_devac=x_is_devac, $
                ;                                /includeerr, $
                ;                                x_is_log=x_is_log, $
                ;                                /fitsersic

        	xyouts, 0.60, 0.80, '!6Initial Disk Profile', /normal, charthick=3.0, size= 1.2, color= 50
	endif


        do_special= 0
        ;do_special= 1
        if do_special eq 1 then begin
                startid= 1600001 
                numpart= 170001
                x=fload_1gal_allstars_xyz('x', startid, numpart)
                y=fload_1gal_allstars_xyz('y', startid, numpart)
                z=fload_1gal_allstars_xyz('z', startid, numpart)
                m=fload_1gal_allstars_mass(startid, numpart)

                process_and_plot_sd, "special", 0, "special", m, xmin, xmax, bins, linecolor=50, $
                                                x_is_devac=x_is_devac, $
                                                x_is_log=x_is_log, $
                                                special_x= x, $
                                                special_y= y, $
                                                special_z= z, $
                                                special_m= m

        endif



; -----------------------------------------------



; extras
; --------

if keyword_set(devac) then begin
    ok = oplot_smoothing_length(smoothlen, /devac)
endif else begin
    ok = oplot_smoothing_length(smoothlen, /plog)
endelse




;  done
; -------
device, /close



end






;============================================================================




