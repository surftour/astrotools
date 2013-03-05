
frun="../../"


f=frun+"/blackholes.txt"

spawn,"wc "+f,result
lines=long(result)
lines=lines(0)


en=dblarr(6,LINES)


openr,1,f
readf,1,en
close,1

time= en(0,*)
num_holes=en(1,*)
mass_holes=en(2,*)
mdot = en(4,*)
mass_actual = en(5,*)


 plot,time, mdot, psym=3,/ylog

 plot,time, mass_holes, psym=3

 oplot,time, mass_actual, psym=3, color=255




end






