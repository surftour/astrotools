pro minor_orbits, junk


sendto= 'ps'
filename= 'orbits.eps'

initialize_plotinfo, 1
setup_plot_stuff, sendto, filename=filename, colortable=4


; must run prodir before this, because this routime calls
; process_kenn_line, which 



; ---------------------
;  Set up Parmaters
; ---------------------

xaxistitle="Time (Gyr)"
yaxistitle="Center Separation (kpc)"

xmax = 4.0
xmin = 0.0
ymax = 250.0
ymin = 0.0


; ------------------
; Plot this up
; ------------------
!p.position= [0.20, 0.15, 0.95, 0.95]

plot, [10.0], [10.0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; ------------------

frun= "execute/G3G3b-u1"
read_centerfiles, frun, time=time, rsep=rsep
oplot, time, rsep, psym=-2, linestyle=0, color= 150
xyouts, "G3G3", x0, y0, /normal, charthick=3.0, size=1.5, color=150


; ------------------

frun= "execute/G3G2-u3"
read_centerfiles, frun, time=time, rsep=rsep
oplot, time, rsep, psym=-2, linestyle=0, color= 120
xyouts, "G3G2", x0, y0-0.05, /normal, charthick=3.0, size=1.5, color=150


; ------------------

frun= "execute/G3G1-u3"
read_centerfiles, frun, time=time, rsep=rsep
oplot, time, rsep, psym=-2, linestyle=0, color= 80
xyouts, "G3G1", x0, y0-0.10, /normal, charthick=3.0, size=1.5, color=150


; ------------------

;frun= "execute/G3G0a-u3"
;read_centerfiles, frun, time=time, rsep=rsep
;oplot, time, rsep, psym=-2, linestyle=0, color= 50
;xyouts, "G3G0", x0, y0-0.15, /normal, charthick=3.0, size=1.5, color=50


; ------------------


if sendto eq 'ps' then device, /close




end











pro read_centerfiles, filename, time=time, rsep=rsep

filename=filename+"/centers.txt"

; read this with the following
openr, 1, filename
junk= ''
lines= 0
readf, 1, junk     ; galaxy 1
print, junk
readf, 1, lines
data= fltarr(4,lines)
readf, 1, data
time= data[0,*]
gal1_xyz= data[1:3,*]
readf, 1, junk     ; galaxy 1
print, junk
readf, 1, lines
data= fltarr(4,lines)
readf, 1, data
gal2_xyz= data[1:3,*]
close, 1


x = gal1_xyz[0,*] - gal2_xyz[0,*]
y = gal1_xyz[1,*] - gal2_xyz[1,*]
z = gal1_xyz[2,*] - gal2_xyz[2,*]

rsep= sqrt(x*x + y*y + z*z)


end





