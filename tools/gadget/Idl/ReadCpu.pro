P_save = !P

PS_PLOT = 0

;WINDOWSIZE=650
WINDOWSIZE=1050

frun="./"

f=frun+"cpu_idl.txt"

spawn,"wc "+f,result
result=strtrim(result, 1)
lines=long(result)
lines=lines(0)
columns = 34

en=fltarr(columns,LINES)

openr,1,f
readf,1,en
close,1


time = fltarr(LINES)
tot = fltarr(LINES)
treegrav = fltarr(LINES)
treebuild = fltarr(LINES)
treeupdate = fltarr(LINES)
treewalk = fltarr(LINES)
treecomm = fltarr(LINES)
treeimbal = fltarr(LINES)
pmgrav = fltarr(LINES)
sph = fltarr(LINES)
density = fltarr(LINES)
denscomm = fltarr(LINES)
densimbal = fltarr(LINES)
hydrofrc = fltarr(LINES)
hydcom = fltarr(LINES)
hydimbal = fltarr(LINES)
hmaxupdate = fltarr(LINES)
domain = fltarr(LINES)
potential = fltarr(LINES)
predict = fltarr(LINES)
kicks = fltarr(LINES)
io = fltarr(LINES)
peano = fltarr(LINES)
sfrcool = fltarr(LINES)
blackholes = fltarr(LINES)
fofsubfind = fltarr(LINES)
smoothing = fltarr(LINES)
hotngbs = fltarr(LINES)
weights_hot = fltarr(LINES)
enrich_hot = fltarr(LINES)
weights_cold =fltarr(LINES)
enrich_cold = fltarr(LINES)
cs_misc = fltarr(LINES)
misc = fltarr(LINES)


time(*) = en(0,*)
tot(*) = en(1,*)
treegrav(*)= en(2,*)
treebuild(*)= en(3,*)
treeupdate(*)= en(4,*)
treewalk(*)= en(5,*)
treecomm(*)= en(6,*)
treeimbal(*)= en(7,*)
pmgrav(*)= en(8,*)
sph(*)= en(9,*)
density(*)= en(10,*)
denscomm(*)= en(11,*)
densimbal(*)= en(12,*)
hydrofrc(*)= en(13,*)
hydcom(*)= en(14,*)
hydimbal(*)= en(15,*)
hmaxupdate(*)= en(16,*)
domain(*)= en(17,*)
potential(*)= en(18,*)
predict(*)= en(19,*)
kicks(*)= en(20,*)
io(*)= en(21,*)
peano(*)= en(22,*)
sfrcool(*)= en(23,*)
blackholes(*)= en(24,*)
fofsubfind(*)= en(25,*)
smoothing(*)= en(26,*)
hotngbs(*)= en(27,*)
weights_hot(*)= en(28,*)
enrich_hot(*)= en(29,*)
weights_cold(*)=en(30,*)
enrich_cold(*)= en(31,*)
cs_misc(*)= en(32,*)
misc(*)= en(33,*)


speed = fltarr( ceil(tot(LINES-1) / 3600) )

time_axis = fltarr(ceil(tot(LINES-1) /3600))
a_axis = fltarr(ceil(tot(LINES-1) /3600))





akt_tot = fltarr(ceil(tot(LINES-1) / 3600))
akt_treegrav = fltarr(ceil(tot(LINES-1) / 3600))
akt_treebuild = fltarr(ceil(tot(LINES-1) / 3600))
akt_treeupdate = fltarr(ceil(tot(LINES-1) / 3600))
akt_treewalk = fltarr(ceil(tot(LINES-1) / 3600))
akt_treecomm = fltarr(ceil(tot(LINES-1) / 3600))
akt_treeimbal = fltarr(ceil(tot(LINES-1) / 3600))
akt_pmgrav = fltarr(ceil(tot(LINES-1) / 3600))
akt_sph = fltarr(ceil(tot(LINES-1) / 3600))
akt_density = fltarr(ceil(tot(LINES-1) / 3600))
akt_denscomm = fltarr(ceil(tot(LINES-1) / 3600))
akt_densimbal = fltarr(ceil(tot(LINES-1) / 3600))
akt_hydrofrc = fltarr(ceil(tot(LINES-1) / 3600))
akt_hydcom = fltarr(ceil(tot(LINES-1) / 3600))
akt_hydimbal = fltarr(ceil(tot(LINES-1) / 3600))
akt_hmaxupdate = fltarr(ceil(tot(LINES-1) / 3600))
akt_domain = fltarr(ceil(tot(LINES-1) / 3600))
akt_potential = fltarr(ceil(tot(LINES-1) / 3600))
akt_predict = fltarr(ceil(tot(LINES-1) / 3600))
akt_kicks = fltarr(ceil(tot(LINES-1) / 3600))
akt_io = fltarr(ceil(tot(LINES-1) / 3600))
akt_peano = fltarr(ceil(tot(LINES-1) / 3600))
akt_sfrcool = fltarr(ceil(tot(LINES-1) / 3600))
akt_blackholes = fltarr(ceil(tot(LINES-1) / 3600))
akt_fofsubfind = fltarr(ceil(tot(LINES-1) / 3600))
akt_smoothing = fltarr(ceil(tot(LINES-1) / 3600))
akt_hotngbs = fltarr(ceil(tot(LINES-1) / 3600))
akt_weights_hot = fltarr(ceil(tot(LINES-1) / 3600))
akt_enrich_hot = fltarr(ceil(tot(LINES-1) / 3600))
akt_weights_cold =fltarr(ceil(tot(LINES-1) / 3600))
akt_enrich_cold = fltarr(ceil(tot(LINES-1) / 3600))
akt_cs_misc = fltarr(ceil(tot(LINES-1) / 3600))
akt_misc = fltarr(ceil(tot(LINES-1) / 3600))





i = 0D
beg=i

j = 0D

while i lt LINES - 1 do begin
    if tot(i) ge tot(beg) + 3600.0 then begin
        speed(j) = 3600.0 * (time(i) - time(beg)) / (tot(i) - tot(beg))
        time_axis(j) = tot(beg + round((i-beg) / 2)) / 3600.0
        a_axis(j) = time(beg + round((i-beg) / 2))
        
        akt_tot(j) = tot(i) - tot(beg)
        akt_treegrav(j) = treegrav(i) - treegrav(beg)
        akt_treebuild(j) = treebuild(i) - treebuild(beg)
        akt_treeupdate(j) = treeupdate(i) - treeupdate(beg)
        akt_treewalk(j) = treewalk(i) - treewalk(beg)
        akt_treecomm(j) = treecomm(i) - treecomm(beg)
        akt_treeimbal(j) = treeimbal(i) - treeimbal(beg)
        akt_pmgrav(j) = pmgrav(i) - pmgrav(beg)
        akt_sph(j) = sph(i) - sph(beg)
        akt_density(j) = density(i) - density(beg)
        akt_denscomm(j) = denscomm(i) - denscomm(beg)
        akt_densimbal(j) = densimbal(i) - densimbal(beg)
        akt_hydrofrc(j) = hydrofrc(i) - hydrofrc(beg)
        akt_hydcom(j) = hydcom(i) - hydcom(beg)
        akt_hydimbal(j) = hydimbal(i) - hydimbal(beg)
        akt_hmaxupdate(j) = hmaxupdate(i) -hmaxupdate(beg)
        akt_domain(j) = domain(i) - domain(beg)
        akt_potential(j) = potential(i) - potential(beg)
        akt_predict(j) = predict(i) - predict(beg)
        akt_kicks(j) = kicks(i) - kicks(beg)
        akt_io(j) = io(i) - io(beg)
        akt_peano(j) = peano(i) - peano(beg)
        akt_sfrcool(j) = sfrcool(i) - sfrcool(beg)
        akt_blackholes(j) = blackholes(i) - blackholes(beg)
        akt_fofsubfind(j) = fofsubfind(i) - fofsubfind(beg)
        akt_smoothing(j) = smoothing(i) - smoothing(beg)
        akt_hotngbs(j) = hotngbs(i) - hotngbs(beg)
        akt_weights_hot(j) = weights_hot(i) - weights_hot(beg)
        akt_enrich_hot(j) = enrich_hot(i) - enrich_hot(beg)
        akt_weights_cold(j) = weights_cold(i) - weights_cold(beg)
        akt_enrich_cold(j) = enrich_cold(i) - enrich_cold(beg)
        akt_cs_misc(j) = cs_misc(i) - cs_misc(beg)
        akt_misc(j) = misc(i) - misc(beg)

        j++
        beg = i
    endif
    i++
endwhile

ind = where(time_axis ne 0.0)

time_axis = time_axis(ind)
a_axis = a_axis(ind)
speed = speed(ind)

akt_tot = akt_tot(ind)
akt_treegrav = akt_treegrav(ind)
akt_treebuild = akt_treebuild(ind)
akt_treeupdate = akt_treeupdate(ind)
akt_treewalk = akt_treewalk(ind)
akt_treecomm = akt_treecomm(ind)
akt_treeimbal = akt_treeimbal(ind)
akt_pmgrav = akt_pmgrav(ind)
akt_sph = akt_sph(ind)
akt_density = akt_density(ind)
akt_denscomm = akt_denscomm(ind)
akt_densimbal = akt_densimbal(ind)
akt_hydrofrc = akt_hydrofrc(ind)
akt_hydcom = akt_hydcom(ind)
akt_hydimbal = akt_hydimbal(ind)
akt_hmaxupdate = akt_hmaxupdate(ind)
akt_domain = akt_domain(ind)
akt_potential = akt_potential(ind)
akt_predict = akt_predict(ind)
akt_kicks = akt_kicks(ind)
akt_io = akt_io(ind)
akt_peano = akt_peano(ind)
akt_sfrcool = akt_sfrcool(ind)
akt_blackholes = akt_blackholes(ind)
akt_fofsubfind = akt_fofsubfind(ind)
akt_smoothing = akt_smoothing(ind)
akt_hotngbs = akt_hotngbs(ind)
akt_weights_hot = akt_weights_hot(ind)
akt_enrich_hot = akt_enrich_hot(ind)
akt_weights_cold = akt_weights_cold(ind)
akt_enrich_cold = akt_enrich_cold(ind)
akt_cs_misc = akt_cs_misc(ind)
akt_misc = akt_misc(ind)




if PS_PLOT eq 0 then begin
    SET_PLOT, 'X'

    WINDOW,0, xsize=WINDOWSIZE,ysize=WINDOWSIZE, retain=2, title='progress'
    textcolor=FSC_COLOR('white')
endif else begin
    set_plot,"PS"
    device,filename='progress.eps',  /encapsulated,  /color ,bits_per_pixel=8
    device,xsize=19,ysize=19
    !p.font=0
    textcolor=FSC_COLOR('black')
endelse



plot, tot/3600.0, time, xtitle='time [hours]', ytitle='a', ystyle=1, pos=[0.11, 0.07, 0.95, 0.7],xrange=[0,max(tot)/3600.0], /xstyle

plot, time_axis, speed, pos=[0.11,0.7,0.95,0.93], xrange=[0,max(tot)/3600.0], /xstyle, /noerase,XTickformat='(A1)', yrange=[-0.0005,0.0035], /ystyle, ytitle='a/h'

Axis, XAxis=1, XTitle='Speed / Progress',xrange=[0,max(tot)/3600.0], /xstyle

oplot, [min(tot)/3600, max(tot)/3600],[0,0], color=FSC_Color('red'), linestyle=1

xyouts, 0.7,0.89, /normal, 'aktual speed: ' + string(speed(n_elements(speed)-1), format = '(F10.6)'), charsize=0.9

xyouts, 0.15, 0.65, /normal, 'aktual scale factor: ' + string(time(n_elements(time)-1), format = '(F10.6)'),charsize=0.9

xyouts, 0.15, 0.62, /normal, 'average speed (last 10 hours): ' + string(total(speed(n_elements(speed)-10:n_elements(speed)-1))/10.0,format = '(F10.6)') ,charsize = 0.9

xyouts, 0.15, 0.59, /normal, 'time remaining: ' + string((1.0-time(n_elements(time)-1))/( total(speed(n_elements(speed)-10:n_elements(speed)-1))/10.0),format = '(F6.1)') + ' hours',charsize = 0.9

if PS_PLOT eq 1 then begin
    device,/close
    set_plot,"X"
endif






if PS_PLOT eq 0 then begin
    WINDOW,1, xsize=WINDOWSIZE,ysize=WINDOWSIZE, retain=2, title='imbalance'
endif else begin
    set_plot,"PS"
    device,filename='imbalance.eps',  /encapsulated,  /color ,bits_per_pixel=8
    device,xsize=19,ysize=19
    !p.font=0
endelse

!p.multi=[0,1,2]

plot, time_axis, 100.0* akt_treeimbal/akt_tot, title='imbalance', xtitle='time[hours]', ytitle='%',xrange=[0,max(time_axis)], /xstyle, yrange=[0.0,30.0], /ystyle
oplot, [0,max(time_axis)],[100.0*treeimbal(LINES-1)/tot(LINES-1),100.0*treeimbal(LINES-1)/tot(LINES-1)], linestyle=1
oplot,  time_axis,100.0* akt_densimbal/akt_tot, color=FSC_Color('red')
oplot, [0,max(time_axis)],[100.0*densimbal(LINES-1)/tot(LINES-1),100.0*densimbal(LINES-1)/tot(LINES-1)], linestyle=1, color=FSC_Color('red')
oplot,  time_axis,100.0* akt_hydimbal/akt_tot, color=FSC_Color('green')
oplot, [0,max(time_axis)],[100.0*hydimbal(LINES-1)/tot(LINES-1),100.0*hydimbal(LINES-1)/tot(LINES-1)], linestyle=1, Color=FSC_Color('green')
legend,['tree','density','hydro'], linestyle=[0,0,0], /left, colors = [ textcolor,FSC_Color('red'),FSC_Color('green')]

plot, a_axis,100.0* akt_treeimbal/akt_tot, title='imbalance', xtitle='a', ytitle='%',xrange=[0,max(time)], /xstyle, yrange=[0.0,30.0], /ystyle
oplot, [0,max(a_axis)],[100.0*treeimbal(LINES-1)/tot(LINES-1),100.0*treeimbal(LINES-1)/tot(LINES-1)], linestyle=1
oplot,  a_axis,100.0* akt_densimbal/akt_tot, color=FSC_color('red')
oplot, [0,max(a_axis)],[100.0*densimbal(LINES-1)/tot(LINES-1),100.0*densimbal(LINES-1)/tot(LINES-1)], linestyle=1, color=FSC_Color('red')
oplot,  a_axis,100.0* akt_hydimbal/akt_tot, color=FSC_color('green')
oplot, [0,max(a_axis)],[100.0*hydimbal(LINES-1)/tot(LINES-1),100.0*hydimbal(LINES-1)/tot(LINES-1)], linestyle=1, Color=FSC_Color('green')
legend,['tree','density','hydro'], linestyle=[0,0,0], /left, colors = [ textcolor,FSC_Color('red'),FSC_Color('green')]

!P = P_save

if PS_PLOT eq 1 then begin
    device,/close
    set_plot,"X"
endif




if PS_PLOT eq 0 then begin
    WINDOW,2, xsize=WINDOWSIZE,ysize=WINDOWSIZE, retain=2, title='composition cummulative'
endif else begin
    set_plot,"PS"
    device,filename='composition_cum.eps',  /encapsulated,  /color ,bits_per_pixel=8
    device,xsize=19,ysize=19
    !p.font=0
endelse

YMARGIN_SAVE= !Y.MARGIN

!Y.MARGIN = [4,8]





xval = tot(0:LINES-1:10)/3600.0

plot_values = n_elements(xval)

xval=[xval(0), xval, xval(plot_values-1)]


value = tot(0:LINES-1:10)
plot_tot=tot(0:LINES-1:10)

plot, tot/3600.0, tot, /nodata, yrange=[0,100], /ystyle, xrange=[0,max(tot)/3600.0], /xstyle, position=[0.065,0.48, 0.97 ,0.83], xtitle='time [hours]', color=textcolor


polyfill, [tot(0)/3600.0,tot(0)/3600.0,tot(LINES-1)/3600.0,tot(LINES-1)/3600.0],[0,100.0,100.0,0], color=FSC_Color('Lime Green')

value -= treegrav(0:LINES-1:10)
value += treeimbal(0:LINES-1:10)


polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Dark Green')

value -= treeimbal(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Blue')

value -= pmgrav(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Red')

value -= sph(0:LINES-1:10)
value += densimbal(0:LINES-1:10)
value += hydimbal(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Orchid')

value -= densimbal(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Dark Red')

value -= hydimbal(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Yellow')

value -= domain(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Cyan')

value -= predict(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Brown')

value -= kicks(0:LINES-1:10)
value -= io(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Magenta')

value -= peano(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Khaki')

value -= sfrcool(0:LINES-1:10)

value -= blackholes(0:LINES-1:10)

value -= fofsubfind(0:LINES-1:10)

polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Firebrick')


plot, [0,0], [0,0], /nodata, yrange=[0,100], /ystyle, xrange=[0,max(tot)/3600.0], /xstyle, /noerase, color=FSC_Color('White'), position=[0.065,0.48, 0.97 ,0.83],XTickformat='(A1)',YTickformat='(A1)'

w = 0.02

xv = 0.065
yv = 0.976

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Lime Green') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "treegrav", /normal, color=textcolor

yv = 0.93

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Dark Green') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "treewait", /normal, color=textcolor

yv = 0.884

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Blue') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "PMgrav", /normal, color=textcolor


xv = 0.22
yv = 0.976

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Red') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "SPH", /normal, color=textcolor

yv = 0.93

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Orchid') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "densimbal", /normal, color=textcolor

yv = 0.884

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Dark Red') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "hydimbal", /normal, color=textcolor


xv = 0.375
yv = 0.976

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Yellow') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "domain", /normal, color=textcolor

yv = 0.93

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Cyan') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "predict", /normal, color=textcolor

yv = 0.884

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Brown') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "kicks/io", /normal, color=textcolor


xv = 0.53
yv = 0.976

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Magenta') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "peano", /normal, color=textcolor

yv = 0.93

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Khaki') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "sfrcool/bh/fof", /normal, color=textcolor

yv = 0.884

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Firebrick') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "misc", /normal, color=textcolor




xval = time(0:LINES-1:10)

plot_values = n_elements(xval)

xval=[xval(0), xval, xval(plot_values-1)]

value = tot(0:LINES-1:10)
plot_tot=tot(0:LINES-1:10)

!Y.MARGIN = YMARGIN_SAVE

plot, time, time, /nodata, yrange=[0,100], /ystyle, xrange=[0,max(time)], /xstyle, position=[0.065, 0.07, 0.97, 0.42], /noerase, xtitle='a', color=textcolor
polyfill, [time(0), time(0), time(LINES-1),time(LINES-1)],[0,100.0,100.0,0],color=FSC_Color('Lime Green')


value -= treegrav(0:LINES-1:10)
value += treeimbal(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Dark Green')

value -= treeimbal(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Blue')

value -= pmgrav(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Red')

value -= sph(0:LINES-1:10)
value += densimbal(0:LINES-1:10)
value += hydimbal(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Orchid')

value -= densimbal(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Dark Red')

value -= hydimbal(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Yellow')

value -= domain(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Cyan')

value -= predict(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Brown')

value -= kicks(0:LINES-1:10)
value -= io(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Magenta')

value -= peano(0:LINES-1:10)
polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Khaki')

value -= sfrcool(0:LINES-1:10)

value -= blackholes(0:LINES-1:10)

value -= fofsubfind(0:LINES-1:10)

polyfill, xval, [0, 100.0 * value/plot_tot,0], color=FSC_Color('Firebrick')

value -= misc(0:LINES-1:10)

plot, time, time, /nodata, yrange=[0,100], /ystyle, xrange=[0,max(time)], /xstyle, /noerase, color=FSC_Color('White'), position=[0.065, 0.07, 0.97, 0.42],XTickformat='(A1)',YTickformat='(A1)'

if PS_PLOT eq 1 then begin
    device,/close
    set_plot,"X"
endif








if PS_PLOT eq 0 then begin
    WINDOW,3, xsize=WINDOWSIZE,ysize=WINDOWSIZE, retain=2, title='composition differential'
endif else begin
    set_plot,"PS"
    device,filename='composition_diff.eps',  /encapsulated,  /color ,bits_per_pixel=8
    device,xsize=19,ysize=19
    !p.font=0
endelse
;xval = [time_axis(0), time_axis, max(time_axis)]
xval = fltarr(2*n_elements(time_axis)+2)
xval[0] = 0.0
xval[1] = 0.0
xval[2:n_elements(xval)-2:2]=time_axis
xval[3:n_elements(xval)-1:2]=time_axis

value = akt_tot

yval = fltarr(2*n_elements(value)+2)
yval[0] = 0.0
yval[n_elements(yval)-1] = 0.0
yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot


plot, akt_tot/3600.0, akt_tot, /nodata, yrange=[0,100], /ystyle, xrange=[0,max(time_axis)], /xstyle, position=[0.065,0.48, 0.97 ,0.83], xtitle='time [hours]', color=textcolor

polyfill, [tot(0)/3600.0,tot(0)/3600.0,tot(LINES-1)/3600.0,tot(LINES-1)/3600.0],[0,100.0,100.0,0], color=FSC_Color('Lime Green')


value -= akt_treegrav
value += akt_treeimbal

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Dark Green')


value -= akt_treeimbal

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Blue')


value -= akt_pmgrav

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Red')


value -= akt_sph
value += akt_densimbal
value += akt_hydimbal

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Orchid')


value -= akt_densimbal

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Dark Red')


value -= akt_hydimbal

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Yellow')


value -= akt_domain

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Cyan')


value -= akt_predict

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Brown')


value -= akt_kicks
value -= akt_io

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Magenta')


value -= akt_peano

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Khaki')


value -= akt_sfrcool
value -= akt_blackholes
value -= akt_fofsubfind

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Firebrick')


plot, [0,0], [0,0], /nodata, yrange=[0,100], /ystyle, xrange=[0,max(time_axis)], /xstyle, /noerase, color=FSC_Color('White'), position=[0.065,0.48, 0.97 ,0.83],XTickformat='(A1)',YTickformat='(A1)'

w = 0.02

xv = 0.065
yv = 0.976

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Lime Green') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "treegrav", /normal, color=textcolor

yv = 0.93

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Dark Green') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "treewait", /normal, color=textcolor

yv = 0.884

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Blue') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "PMgrav", /normal, color=textcolor


xv = 0.22
yv = 0.976

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Red') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "SPH", /normal, color=textcolor

yv = 0.93

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Orchid') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "densimbal", /normal, color=textcolor

yv = 0.884

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Dark Red') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "hydimbal", /normal, color=textcolor


xv = 0.375
yv = 0.976

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Yellow') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "domain", /normal, color=textcolor

yv = 0.93

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Cyan') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "predict", /normal, color=textcolor

yv = 0.884

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Brown') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "kicks/io", /normal, color=textcolor


xv = 0.53
yv = 0.976

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Magenta') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "peano", /normal, color=textcolor

yv = 0.93

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Khaki') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "sfrcool/bh/fof", /normal, color=textcolor

yv = 0.884

polyfill, [xv, xv+w, xv+w, xv], [yv, yv, yv-w, yv-w], color=FSC_Color('Firebrick') ,/normal
xyouts, xv + 1.5 * w, yv-w/2, "misc", /normal, color=textcolor




;xval = [a_axis(0), a_axis, max(a_axis)]
xval = fltarr(2*n_elements(a_axis)+2)
xval[0] = time(0)
xval[1] = time(0)
xval[2:n_elements(xval)-2:2]=a_axis
xval[3:n_elements(xval)-1:2]=a_axis


value = akt_tot

plot, time, time, /nodata, yrange=[0,100], /ystyle, xrange=[0,max(a_axis)], /xstyle, position=[0.065, 0.07, 0.97, 0.42], /noerase, xtitle='a', color=textcolor

value = akt_tot

yval = fltarr(2*n_elements(value)+2)
yval[0] = 0.0
yval[n_elements(yval)-1] = 0.0
yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, [time(0), time(0), time(LINES-1),time(LINES-1)],[0,100.0,100.0,0],color=FSC_Color('Lime Green')


value -= akt_treegrav
value += akt_treeimbal

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Dark Green')


value -= akt_treeimbal

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Blue')


value -= akt_pmgrav

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Red')


value -= akt_sph
value += akt_densimbal
value += akt_hydimbal

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Orchid')


value -= akt_densimbal

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Dark Red')


value -= akt_hydimbal

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Yellow')


value -= akt_domain

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Cyan')


value -= akt_predict

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Brown')


value -= akt_kicks
value -= akt_io

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Magenta')


value -= akt_peano

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Khaki')


value -= akt_sfrcool
value -= akt_blackholes
value -= akt_fofsubfind

yval[1:n_elements(yval)-3:2]=100.0 * value/akt_tot
yval[2:n_elements(yval)-2:2]=100.0 * value/akt_tot

polyfill, xval, yval, color=FSC_Color('Firebrick')

plot, time, time, /nodata, yrange=[0,100], /ystyle, xrange=[0,max(a_axis)], /xstyle, /noerase, color=FSC_Color('White'), position=[0.065, 0.07, 0.97, 0.42],XTickformat='(A1)',YTickformat='(A1)'

if PS_PLOT eq 1 then begin
    device,/close
    set_plot,"X"
endif

!P = P_save

end



