
pro sigz, junk

snapnum= 0
;snapnum= 10

; n= 2
;frun= 'execute/Sbc11i4-u57'
frun= '/data/tcox/Sbc11i4-u4'
;frun= 'execute/Sbc11i4-u43'

; n=1
;frun= 'execute/Sbc11i4-u56'
;frun= 'execute/Sbc11i4-u53'
;frun= 'execute/Sbc11i4-u52'

; n=0
;frun= 'execute/Sbc11i4-u55'
;frun= 'execute/Sbc11i4-u50'
;frun= 'execute/Sbc11i4-u51'


process_it, frun, snapnum

end



pro process_it, frun, snapnum


ok= fload_snapshot(frun,snapnum)


x= fload_allstars_xyz('x')
y= fload_allstars_xyz('y')
Radius= sqrt(x*x + y*y)
z= fload_allstars_xyz('z')

vz= fload_allstars_v('z')

nbins= 30
xlen= 30.0

dr= xlen/nbins

for i=0, nbins-1 do begin
	r_in= dr*i
	r_out= dr*(i+1)

	ridx= where((Radius ge r_in) and (Radius lt r_out))

	zin= z(ridx)
	vzin= vz(ridx)

	;print, 'r= ',r_in,r_out,mean(zin), sqrt(variance(zin))
	print, 'r= ',r_in,r_out,mean(vzin), sqrt(variance(vzin))
endfor


end
