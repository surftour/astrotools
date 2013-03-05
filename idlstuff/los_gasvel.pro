

;-----------------------------------------------------------------------------------

pro los_gasvel_hist, junk



if not keyword_set(junk) then begin
        print, " "
        print, " los_gasvel_hist, junk"
        print, " "
        print, " "
        return 
endif




; -------------------
filename='losgasvel.eps'



initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


; -----------------------------
; set up constants
; -----------------------------


; alignement angle
xaxistitle= "!6LOS Velocity (km sec!E-1!N)"
xmax =  1000.0
xmin = -1000.0

; number (histogram)
yaxistitle= ' '
ymax = 1.05
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.08
y0= 0.15
x1= 0.92
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
        ;/ylog, $
        ;/xlog, $ 
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
        xtitle=xaxistitle, $
        ytickformat='(a1)', $
        ;ytitle=yaxistitle, $
        /nodata   ;, /noerase






; =====================
; =====================


snapnum= 62
frun= "data/ds/vc3vc3e_2"

ok= fload_snapshot_bh(frun,snapnum)

; 
; los gas velocity
oplotit= 50
vx = fload_gas_v('x')
vy = fload_gas_v('y')
vz = fload_gas_v('z')

;
;m= fload_gas_mass(1,/HI)
rho= fload_gas_rho(1)
temp= fload_gas_temperature(1)

idx= where((rho lt 0.001) and (temp lt 1.0e5))
vx=vx(idx)
vy=vy(idx)
vz=vz(idx)

n_hat= [1.0, 0.0, 0.0]    ; along x axis
vlos= vx*n_hat[0] + vy+n_hat[1] + vz*n_hat[2]
print, "N= ", n_elements(vlos), min(vlos), max(vlos)
temp= process_histogram(vlos, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, normalization=normalization)
print, min(temp), max(temp)

n_hat= [0.0, 1.0, 0.0]    ; along y axis
vlos= vx*n_hat[0] + vy+n_hat[1] + vz*n_hat[2]
print, "N= ", n_elements(vlos), min(vlos), max(vlos)
temp= process_histogram(vlos, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, normalization=normalization)
print, min(temp), max(temp)

n_hat= [0.0, 0.0, 1.0]    ; along z axis
vlos= vx*n_hat[0] + vy+n_hat[1] + vz*n_hat[2]
print, "N= ", n_elements(vlos), min(vlos), max(vlos)
temp= process_histogram(vlos, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, normalization=normalization)
print, min(temp), max(temp)



xyouts, 0.11, 0.93, frun, /normal, size= 1.0, charthick=3.0, color= 0
xyouts, 0.13, 0.88, fload_timelbl(1,2), /normal, size= 1.0, charthick=3.0, color= 0




; -----------------------------------------------------------------------------



; -------------
;  Done
; -------------

device, /close



end








