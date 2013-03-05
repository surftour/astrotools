
frun="../../"


f=frun+"/sfr.txt"

spawn,"wc "+f,result
lines=long(result)
lines=lines(0)


en=dblarr(5,LINES)

openr,1,f
readf,1,en
close,1

ti= en(0,*)
sfr=en(2,*)


plot,ti, sfr, psym=3




end






