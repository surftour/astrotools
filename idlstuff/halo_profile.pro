;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   Determine 3D Halo Profile and Velocity Structure
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------


pro halostructure, runinfo, $
			ctab=ctab, $
			filename=filename, $
			smoothlen=smoothlen


if not keyword_set(runinfo) then begin
   print, "  "
   print, "halostructure, runinfo"
   print, "  " 
   return
endif

if not keyword_set(smoothlen) then smoothlen=0.1
if not keyword_set(filename) then filename='halo.eps'
if not keyword_set(ctab) then ctab= 4

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= ctab, newxsize= 20.0, newysize= 20.0



                
; ----------------------------- 
; ----------------------------- 
;  Global Quantities
; -----------------------------
; ----------------------------- 

bins = 50
                
xaxistitle= 'Log r (kpc)'
xmax = 3.0
xmin = -1.0

n_runs= (size(runinfo))[1]



                


; ----------------------------- 
; ----------------------------- 
;  Top left: density profile
; -----------------------------
; ----------------------------- 

ymax =  1.0
ymin = -6.0

yaxistitle= 'Log !7q!3!N (M!D!9n!3!N pc!E-3!N) '

divide_by_rhocrit= 0
;divide_by_rhocrit= 1

if divide_by_rhocrit eq 1 then begin
        yaxistitle="Log !7q!3!N(r) / !7q!3!Dcrit!N"
        ymax = 8.0
        ymin = -2.0
endif


;---------------------------

!p.position= [0.11, 0.57, 0.48, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
        ;/ylog, $
        ;/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; ------------------------


for i=0, n_runs-1 do begin

        print, " "
        print, "i= ", i
        print, " "

        frun= runinfo[i].frun
        snapnum= runinfo[i].snapnum
        component= runinfo[i].component
        linecolor= runinfo[i].linecolor


	process_and_plot_rho, frun, snapnum, component, xmin, xmax, bins, linecolor= linecolor, /x_is_log, /includeerr


endfor



; -----------------
;  Plot Extras

;x0= 0.55
;y0= 0.86
;xyouts, x0, y0, fload_fid(1), /normal, charthick= 2.0, size= 1.7, color= 0
;xyouts, x0-0.23, y0+0.04, "/home/tcox/MakeNewDisk/vc/vc3c.dat", /normal, charthick= 2.0, size= 1.1, color= 0
;y0= 0.80
;xyouts, x0, y0+0.05, 'Dark Matter Profile', /normal, charthick= 2.0, size= 1.2, color= 0
;xyouts, x0, y0, 'Stellar Profile', /normal, charthick= 2.0, size= 1.2, color= 100



; smoothing length
;ok = oplot_smoothing_length(smoothlen, /plog)





; ----------------------------- 
; ----------------------------- 
;  Top right: Radial Dispersion
; -----------------------------
; ----------------------------- 

ymax = 350.0
ymin = 0.0

yaxistitle= '!7r!3!Dr!N (km s!E-1!N) '


;---------------------------

!p.position= [0.60, 0.57, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, color= 0, /nodata, /noerase, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, charthick=3.0, $
        xtitle=xaxistitle, ytitle=yaxistitle


; ------------------------


for i=0, n_runs-1 do begin

        print, " "
        print, "i= ", i
        print, " "

        frun= runinfo[i].frun
        snapnum= runinfo[i].snapnum
        component= runinfo[i].component
        linecolor= runinfo[i].linecolor

	process_and_plot_3dvel, frun, snapnum, component+'/sigr', xmin, xmax, bins, linecolor= linecolor, /x_is_log

endfor






; ------------------------------------ 
; ------------------------------------
;  Bottom left: tangential (theta, actually) dispersion
; ------------------------------------
; ------------------------------------


ymax = 350.0
ymin = 0.0

;yaxistitle= '!7r!3!Dtheta!N (km s!E-1!N) '
yaxistitle= '!7r!Dh!N!3 (km s!E-1!N) '


;---------------------------

!p.position= [0.11, 0.10, 0.48, 0.48]

plot, [0], [0], psym=-3, linestyle=0, color= 0, /nodata, /noerase, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, charthick=3.0, $
        xtitle=xaxistitle, ytitle=yaxistitle


; ------------------------


for i=0, n_runs-1 do begin

        print, " "
        print, "i= ", i
        print, " "

        frun= runinfo[i].frun
        snapnum= runinfo[i].snapnum
        component= runinfo[i].component
        linecolor= runinfo[i].linecolor


	process_and_plot_3dvel, frun, snapnum, component+'/sigtheta', xmin, xmax, bins, linecolor= linecolor, /x_is_log


endfor




; ------------------------------------ 
; ------------------------------------
;  Bottom right: anisotropy
; ------------------------------------
; ------------------------------------


ymax =  1.0
ymin = -1.0

yaxistitle= '!7b!3'


;---------------------------

!p.position= [0.60, 0.10, 0.98, 0.48]

plot, [0], [0], psym=-3, linestyle=0, color= 0, /nodata, /noerase, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, charthick=3.0, $
        xtitle=xaxistitle, ytitle=yaxistitle


; ------------------------


for i=0, n_runs-1 do begin

        print, " "
        print, "i= ", i
        print, " "

        frun= runinfo[i].frun
        snapnum= runinfo[i].snapnum
        component= runinfo[i].component
        linecolor= runinfo[i].linecolor


        process_and_plot_3dvel, frun, snapnum, component+'/beta', xmin, xmax, bins, linecolor= linecolor, /x_is_log

endfor







;--------------------------------------
;--------------------------------------

device, /close


print, " "
print, " "
print, " output to file: ", filename
print, " "
print, " "



end







;--------------------------------------------------------------------------



