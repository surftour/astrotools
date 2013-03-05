;================================================================================
;
;
;
;
;================================================================================



pro sb99, junk

;
; read in the sb99 yield file
;
;------------------------

;sb99file="/Users/tcox/Desktop/tjtest.yield1"
sb99file="/Users/tcox/Desktop/tjtest3.yield1"

spawn, "wc "+sb99file,result
lines=long(result)
datalines=lines(0)-7

ntot= datalines

time= fltarr(datalines)
el_H= fltarr(datalines)
el_He= fltarr(datalines)
el_C= fltarr(datalines)
el_N= fltarr(datalines)
el_O= fltarr(datalines)
el_Mg= fltarr(datalines)
el_Si= fltarr(datalines)
el_S= fltarr(datalines)
el_Fe= fltarr(datalines)
chem_winds= fltarr(datalines)
chem_SN= fltarr(datalines)
chem_total= fltarr(datalines)
totmass= fltarr(datalines)

openr, 1, sb99file
print, "opening: ", sb99file
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,datalines-1 do begin
        readf, 1, junk
        ;print, junk
        tempjunk= strsplit(junk,/extract,count=count)
	time[i]= float(tempjunk(0))
	el_H[i]= float(tempjunk(1))
	el_He[i]= float(tempjunk(2))
	el_C[i]= float(tempjunk(3))
	el_N[i]= float(tempjunk(4))
	el_O[i]= float(tempjunk(5))
	el_Mg[i]= float(tempjunk(6))
	el_Si[i]= float(tempjunk(7))
	el_S[i]= float(tempjunk(8))
	el_Fe[i]= float(tempjunk(9))
	chem_winds[i]= float(tempjunk(10))
	chem_SN[i]= float(tempjunk(11))
	chem_total[i]= float(tempjunk(12))
	totmass[i]= float(tempjunk(13))
endfor

close, 1


;--------------------------------



filename='sb99_mreturn.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


; -----------------------------
; set up constants
; -----------------------------

xaxistitle= "!6Log Time (Gyr)"
;xmax = 10.0
xmax = 9.5
xmin =  4.0

yaxistitle= 'Log Mass Return Rate (M!D!9n!3!N / Yr)'
ymax = -1.0
ymin = -4.2


; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.98

!p.position= [x0, y0, x1, y1]

plot, [0], [0], psym=-3, linestyle=0, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        color= 0, xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, charthick=4.0, /nodata, /noerase, $
        xtitle=xaxistitle, ytitle=yaxistitle



time= alog10(time)
print, time[0:10]

; these are already in log
;chem_winds= alog10(chem_winds)
;chem_SN= alog10(chem_SN)
;chem_total= alog10(chem_total)


oplot, time, chem_SN, psym=-3, symsize= 1.0, color= 150, thick= 3.0
xyouts, 0.22, 0.30, 'Supernovae', size=1.3, color= 150, charthick=3.0, /normal

oplot, time, chem_winds, psym=-3, symsize= 1.0, color= 100, thick= 3.0
xyouts, 0.22, 0.25, 'Winds', size=1.3, color= 100, charthick=3.0, /normal

oplot, time, chem_total, psym=-3, symsize= 3.0, color= 0, thick= 3.0
xyouts, 0.22, 0.35, 'Total', size=1.3, color= 0, charthick=3.0, /normal


print, "Total mass returned: ", totmass(ntot-1)
print, "      mass return fraction: ", 10^(totmass(ntot-1))/1.0e+6


; done
; -------------
device, /close


end







