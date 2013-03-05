;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------
;
;  Plot something
;
;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------

pro t1, junk


if not keyword_set(junk) then begin
	print, " "
	print, " t1, junk"
	print, " "
	print, " "
	return
endif

;filename='jy_test1.eps'
;filename='jy_test2.eps'
;filename='jy_test3.eps'
filename='jy_test4.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=4
setup_plot_stuff, 'ps', filename=filename, colortable=3


;----------------------------

x0= 0.18
x1= 0.98

y0= 0.15
y1= 0.98

; 
;----------------------------


;do_t1= 1
do_t1= 0

if do_t1 eq 1 then begin

	xaxistitle = "Log Shock Factor (dD/1000)"
	xmax = 10.0
	xmin = -5.0

	yaxistitle = "Log Density Factor (!7q!6/!7q!6!Dcrit!N)"
	ymax = 6.0
	ymin = -5.0


	read_pdata, densityfactor, shockfactor, soundspeed, velocity, tsfr

	ldf= alog10(densityfactor)
	lsf= alog10(shockfactor)

	bins= 100

	contour_makegeneralpic, lsf, ldf, xmax, xmin, ymax, ymin, $
                                pixels= bins, $
                                NxNImage=NxNImage



	tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

endif



; 
;----------------------------


;do_t2= 1
do_t2= 0

if do_t2 eq 1 then begin

        xaxistitle = "Log Shock Factor (dD/1000)"
        xmax = 10.0
        xmin = -5.0

        yaxistitle = "Log Mach Number "
        ymax = 2.0
        ymin = -2.0


        read_pdata, densityfactor, shockfactor, soundspeed, velocity, tsfr

        lmach= alog10(velocity/soundspeed)
        lsf= alog10(shockfactor)

        bins= 100

        contour_makegeneralpic, lsf, lmach, xmax, xmin, ymax, ymin, $
                                pixels= bins, $
                                NxNImage=NxNImage



        tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

endif


; 
;----------------------------


;do_t3= 1
do_t3= 0

if do_t3 eq 1 then begin

        xaxistitle = "Log t!DSFR!N (Gyr)"
        xmax = 3.0
        xmin = -4.5

        yaxistitle = "Log Density Factor (!7q!6/!7q!6!Dcrit!N)"
        ymax = 6.0
        ymin = -5.0


        read_pdata, densityfactor, shockfactor, soundspeed, velocity, tsfr

	tsfr_exact= 4.5 / sqrt(densityfactor)
        ltsfr= alog10(tsfr_exact)

        ldf= alog10(densityfactor)
        ;ltsfr= alog10(tsfr)

        bins= 100

        contour_makegeneralpic, ltsfr, ldf, xmax, xmin, ymax, ymin, $
                                pixels= bins, $
                                NxNImage=NxNImage



        tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

endif


; 
;----------------------------


do_t4= 1
;do_t4= 0

if do_t4 eq 1 then begin

        ;xaxistitle = "Log t!D*!N (dD/1000)!E-0.5!N [Gyr]"
        ;xaxistitle = "Log t!D*!N (dD/1000)!E-1!N [Gyr]"
        xaxistitle = "Log t!D*!N (dD/1e7)!E-1!N [Gyr]"
        xmax = 3.0
        xmin = -4.5

        yaxistitle = "Log t!D*!N (!7q!6!Diso shock!N/!7q!6!Dcrit!N)!E-0.5!N [Gyr]"
        ymax = 6.0
        ymin = -5.0


        read_pdata, densityfactor, shockfactor, soundspeed, velocity, tsfr

	mach= velocity/soundspeed
	denfrmshck= densityfactor * mach * mach

        tsfr_exact= 4.5 / sqrt(denfrmshck)
        ltsfr= alog10(tsfr_exact)

	;shocksf= 4.5 / (shockfactor)^(1/2.)
	;shocksf= 4.5 / (shockfactor)^(1.)
	shocksf= 4.5 / (shockfactor/1.0e+4)^(1.)
	lshocksf= alog10(shocksf)

        bins= 100

        contour_makegeneralpic, lshocksf, ltsfr, xmax, xmin, ymax, ymin, $
                                pixels= bins, $
                                NxNImage=NxNImage



        tv, NxNImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

endif




;---------------------------

!p.position= [x0, y0, x1, y1]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
	;/xlog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata, /noerase





; 
;----------------------------


;x= [xmin, xmax]
;y= 4.5 / sqrt(10^x)
;oplot, alog10(y), x, psym=-3, linestyle= 1, color= 0, thick= 2.0


x= [-10, 10]
y= x
oplot, y, x, psym=-3, linestyle= 1, color= 0, thick= 2.0



;
;  done
;---------------
device, /close


end






;
;
;===================================================
pro read_pdata, densityfactor, shockfactor, soundspeed, velocity, tsfr


pdatafile= '/raid4/tcox/z5/test0/pdata.txt'


spawn, "wc "+pdatafile,result
lines=long(result)
datalines=lines(0)
densityfactor= fltarr(datalines)
shockfactor= fltarr(datalines)
soundspeed= fltarr(datalines)
velocity= fltarr(datalines)
tsfr= fltarr(datalines)

densityfactor(*)= -1
shockfactor(*)= -1
soundspeed(*)= -1
velocity(*)= -1
tsfr(*)= -1

openr, 1, pdatafile
junk=''


; read the data
for i=0L,lines(0)-1 do begin
;for i=0L,200000 do begin
	readf, 1, junk
	tempjunk= strsplit(junk,/extract,count=count)
	if count eq 8 then begin
		densityfactor(i)= float(tempjunk(2))
		shockfactor(i)=   float(tempjunk(4))
		soundspeed(i)=    float(tempjunk(5))
		velocity(i)=      float(tempjunk(6))
		tsfr(i)=          float(tempjunk(7))
	endif

	if (i mod 50000) eq 0 then print, "i= ", i
endfor

close, 1


idx= where(densityfactor gt -1)
if idx(0) ne -1 then begin
	densityfactor= densityfactor(idx)
	shockfactor= shockfactor(idx)
	soundspeed= soundspeed(idx)
	velocity= velocity(idx)
	tsfr= tsfr(idx)
endif


end



