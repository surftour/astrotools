pro xr_redump
nux = [16.00, 16.02, 16.04, 16.06, 16.08, 16.10, 16.12, 16.14, 16.16, 16.18, 16.20, 16.22, 16.24,$
16.26, 16.28, 16.30, 16.32, 16.34, 16.36, 16.38, 16.40, 16.42, 16.44, 16.46, 16.48, 16.50,$
16.52, 16.54, 16.56, 16.58, 16.60, 16.62, 16.64, 16.66, 16.68, 16.70, 16.72, 16.74, 16.76,$
16.78, 16.80, 16.82, 16.84, 16.86, 16.88, 16.90, 16.92, 16.94, 16.96, 16.98, 17.00, 17.02,$
17.04, 17.06, 17.08, 17.10, 17.12, 17.14, 17.16, 17.18, 17.20, 17.22, 17.24, 17.26, 17.28,$
17.30, 17.32, 17.34, 17.36, 17.38, 17.40, 17.42, 17.44, 17.46, 17.48, 17.50, 17.52, 17.54,$
17.56, 17.58, 17.60, 17.62, 17.64, 17.66, 17.68, 17.70, 17.72, 17.74, 17.76, 17.78, 17.80,$
17.82, 17.84, 17.86, 17.88, 17.90, 17.92, 17.94, 17.96, 17.98, 18.00, 18.02, 18.04, 18.06,$
18.08, 18.10, 18.12, 18.14, 18.16, 18.18, 18.20, 18.22, 18.24, 18.26, 18.28, 18.30, 18.32,$
18.34, 18.36, 18.38, 18.40, 18.42, 18.44, 18.46, 18.48, 18.50, 18.52, 18.54, 18.56, 18.58,$
18.60, 18.62, 18.64, 18.66, 18.68, 18.70, 18.72, 18.74, 18.76, 18.78, 18.80, 18.82, 18.84,$
18.86, 18.88, 18.90, 18.92, 18.94, 18.96, 18.98, 19.00, 19.02, 19.04, 19.06, 19.08, 19.10,$
19.12, 19.14, 19.16, 19.18, 19.20, 19.22, 19.24, 19.26, 19.28, 19.30, 19.32, 19.34, 19.36,$
19.38, 19.40, 19.42, 19.44, 19.46, 19.48, 19.50, 19.52, 19.54, 19.56, 19.58, 19.60, 19.62,$
19.64, 19.66, 19.68, 19.70, 19.72, 19.74, 19.76, 19.78, 19.80, 19.82, 19.84, 19.86, 19.88,$
19.90, 19.92, 19.94, 19.96, 19.98, 20.00, 20.02, 20.04, 20.06, 20.08, 20.10, 20.12, 20.14,$
20.16, 20.18, 20.20, 20.22, 20.24, 20.26, 20.28, 20.30, 20.32, 20.34, 20.36, 20.38, 20.40,$
20.42, 20.44, 20.46, 20.48, 20.50, 20.52, 20.54, 20.56, 20.58, 20.60, 20.62, 20.64, 20.66,$
20.68, 20.70, 20.72, 20.74, 20.76, 20.78, 20.80, 20.82, 20.84, 20.86, 20.88, 20.90, 20.92,$
20.94, 20.96, 20.98, 21.00, 21.02, 21.04, 21.06, 21.08, 21.10, 21.12, 21.14, 21.16, 21.18,$
21.20, 21.22, 21.24, 21.26, 21.28, 21.30, 21.32, 21.34, 21.36, 21.38, 21.40, 21.42, 21.44,$
21.46, 21.48 ]

f=load_xr_spec(nux)

	
	okbig = where(nux GT 16.38)
	nook  = where(nux LE 16.38)
	 f[nook]=10^(INTERPOL(alog10(f[okbig]),nux[okbig],nux[nook]))

	openw,1,'dumper'
	printf,1,abs(alog10(f)),FORMAT='(F7.4)'
	close,1

	plot,10^(nux)/2.418d17,alog10(f),xstyle=1,xrange=[0.01,500.0],/xlog

yy=[-2.11, -2.11, -2.11, -2.10, -2.10, -2.09, -2.09, -2.09, -2.08, -2.08, -2.07, -2.07, -2.07,$
-2.06, -2.06, -2.05, -2.05, -2.05, -2.04, -2.04, -2.03, -2.03, -2.03, -2.02, -2.02, -2.01,$
-2.01, -2.01, -2.00, -2.00, -1.99, -1.99, -1.99, -1.98, -1.98, -1.97, -1.97, -1.97, -1.96,$
-1.96, -1.95, -1.95, -1.95, -1.94, -1.94, -1.93, -1.93, -1.93, -1.92, -1.92, -1.91, -1.91,$
-1.91, -1.90, -1.90, -1.89, -1.89, -1.89, -1.88, -1.88, -1.87, -1.87, -1.87, -1.86, -1.86,$
-1.85, -1.85, -1.85, -1.84, -1.84, -1.83, -1.83, -1.83, -1.82, -1.82, -1.81, -1.81, -1.81,$
-1.80, -1.80, -1.79, -1.79, -1.78, -1.78, -1.78, -1.77, -1.77, -1.76, -1.76, -1.75, -1.75,$
-1.75, -1.74, -1.74, -1.73, -1.73, -1.72, -1.72, -1.71, -1.71, -1.70, -1.69, -1.69, -1.68,$
-1.67, -1.67, -1.66, -1.65, -1.65, -1.64, -1.63, -1.62, -1.63, -1.63, -1.62, -1.61, -1.60,$
-1.60, -1.59, -1.58, -1.57, -1.55, -1.54, -1.53, -1.51, -1.50, -1.48, -1.46, -1.45, -1.43,$
-1.41, -1.39, -1.37, -1.35, -1.33, -1.32, -1.30, -1.29, -1.28, -1.27, -1.26, -1.25, -1.24,$
-1.24, -1.23, -1.23, -1.22, -1.22, -1.22, -1.22, -1.22, -1.21, -1.21, -1.21, -1.22, -1.22,$
-1.22, -1.22, -1.22, -1.22, -1.23, -1.23, -1.23, -1.23, -1.24, -1.24, -1.25, -1.25, -1.25,$
-1.26, -1.26, -1.27, -1.28, -1.28, -1.29, -1.30, -1.30, -1.31, -1.32, -1.33, -1.34, -1.35,$
-1.36, -1.37, -1.39, -1.40, -1.41, -1.43, -1.45, -1.46, -1.47, -1.49, -1.50, -1.51, -1.53,$
-1.54, -1.56, -1.57, -1.59, -1.60, -1.62, -1.64, -1.65, -1.67, -1.69, -1.71, -1.73, -1.75,$
-1.77, -1.80, -1.82, -1.85, -1.87, -1.90, -1.93, -1.96, -1.99, -2.03, -2.06, -2.10, -2.14,$
-2.18, -2.22, -2.26, -2.31, -2.36, -2.41, -2.46, -2.51, -2.57, -2.63, -2.70, -2.76, -2.83,$
-2.91, -2.99, -3.07, -3.15, -3.24, -3.34, -3.44, -3.54, -3.65, -3.76, -3.88, -4.01, -4.14,$
-4.28, -4.42, -4.57, -4.73, -4.90, -5.07, -5.25, -5.45, -5.65, -5.86, -6.08, -6.31, -6.55,$
-6.81, -7.07, -7.35, -7.64, -7.95, -8.27, -8.61, -8.96, -9.33, -9.71,-10.12,-10.54,-10.99,$
-11.45,-11.94]
	oplot,10^(nux)/2.418d17,yy,color=250
	
	oplot,10^(nux)/2.418d17,agn_spectrum(10^(nux),12.0)-11.825,color=50
print, agn_spectrum(10^(nux),12.0)-12.
end

pro vdb
setup_plot,'idltemp.ps',/PS

bhfile= return_idl_routines_homedir(0)+'/agn_spectrum/vdbspecdat.txt'
spawn, "wc "+bhfile,result
lines=long(result)
lines=lines(0)
if lines gt 0 then bhdata= fltarr(5,lines)
openr, 1, bhfile, ERROR=err
if (err NE 0) then begin
        print, "  "
        print, "Problem: ",!ERR_STRING
        print, "  "
	close, 1
	ERR= 0
	bhmsg= "No Black Hole"
endif else begin
	print, "opening: ",bhfile
	bhdata= read_ascii(bhfile)
	lambda= bhdata.field1[0,*]
	f_lambda= bhdata.field1[1,*]
	d_f_lambda= bhdata.field1[2,*]
	close, 1
endelse
;; lambda*f_lambda = nu*f_nu
nufnu = lambda * f_lambda
nu    = 2.998d8/(lambda * 1.0d-10)

plot,lambda, f_lambda
plot,alog10(nu),nufnu,/ylog,yrange=[1.0d2,1.0d5],ystyle=1,xrange=[14.5,15.6],xstyle=1
;;plot,lambda,nufnu

f = richards_return_spectrum(alog10(nu),/RED)
	f_renorm = TOTAL(f * (nu[0]-nu[1])/nu)
	f = f * (1.0d4/f_renorm)
oplot,alog10(nu),f*20.0,COLOR=250
f = richards_return_spectrum(alog10(nu),/BLUE)
	f_renorm = TOTAL(f * (nu[0]-nu[1])/nu)
	f = f * (1.0d4/f_renorm)
oplot,alog10(nu),f*20.0,COLOR=80
f = richards_return_spectrum(alog10(nu),/ALL)
	f_renorm = TOTAL(f * (nu[0]-nu[1])/nu)
	f = f * (1.0d4/f_renorm)
oplot,alog10(nu),f*25.0,COLOR=150
f_gtr_all = f*25.0
f = marconi_return_spectrum(nu,1.0d12)
	f = 10^(f-7.6)
oplot,alog10(nu),f,COLOR=210,thick=2.0


	;; ok, continuum there well-approximated for log_nu > 14.9 -- below, too 
	;;   shallow -- this is corrected in Haziminaglou et al. 2006 -- we need to 
	;;   tilt the spectral slope appropriately to match
falt = 6.0d3 * (nu/10^(14.6))^(-0.55)
	oplot,alog10(nu),falt,color=250,thick=3.0
	;; now correct the vdb spectrum down to this where below

falt2 = 1.0d4 * (nu/10^(15.3))^(-2.0)
	oplot,alog10(nu),falt2,color=250,thick=3.0
	;; now correct the vdb spectrum down to this where below
	
	needcor = where((falt GT f_gtr_all) AND (alog10(nu) LT 15.2))
	nufnu[needcor] = nufnu[needcor] * f_gtr_all[needcor]/falt[needcor]
	oplot,alog10(nu[needcor]),nufnu[needcor],COLOR=110

	;; at the highest nu - the spectrum cuts off hard due to a mix of 
	;;   instrumental sensitivity and IGM absorption -- cut off using this 
	;;   at ~ log_nu=15.5
numax = 15.5
oplot,[numax,numax],[1.0d-10,1.0d10],THICK=2.0,COLOR=0

f_gtr_all_tmp = f_gtr_all * 0.9
	needcor = where((falt2 LT f_gtr_all_tmp) AND (alog10(nu) GT 15.0))
	nufnu[needcor] = nufnu[needcor] * f_gtr_all_tmp[needcor]/falt2[needcor]
	oplot,alog10(nu[needcor]),nufnu[needcor],COLOR=110


;; take the GTR median all spectrum to be an effective continuum
vdb_over_continuum = nufnu / f_gtr_all

plot,alog10(nu),vdb_over_continuum,/ylog,ystyle=1,yrange=[0.5,3.0]
numin = 14.575
	oplot,[numin,numin],[1.0d-10,1.0d10],THICK=2.0,COLOR=0
numax = 15.50
	oplot,[numax,numax],[1.0d-10,1.0d10],THICK=2.0,COLOR=0

nu_min=numin
nu_max=numax
;; cutoff at nu_min and nu_max
ok = where((alog10(nu) GT numin) AND (alog10(nu) LT numax))
nu_ok = nu[ok]
f_ok  = vdb_over_continuum[ok]

;; re-interpolate over an appropriate grid
d_log_nu = 0.001
nu_grid = nu_min + d_log_nu * findgen((nu_max-nu_min)/d_log_nu+1.0)
f_now = 10^(INTERPOL(alog10(vdb_over_continuum),alog10(nu),nu_grid))

;; determine the net flux over continuum in the bandpass
f_int  = TOTAL(f_now * d_log_nu * alog(10.))
f_flat = TOTAL((0.0*f_now+1.0) *  d_log_nu * alog(10.))
print, 'TOTAL flux ratio = ',f_int,f_flat,f_int/f_flat

;; renormalize to conserve the integrated luminosity of the continuum over the pass
f_renorm = f_now / (f_int/f_flat)
oplot,nu_grid,f_renorm,COLOR=110


print, n_elements(nu_grid), n_elements(f_renorm)
openw,1,'dumper'
printf,1,nu_grid,FORMAT='(F6.3)'
printf,1,f_renorm,FORMAT='(F6.3)'
close,1

nux = 0.*nu_grid
fnx = 0.*nu_grid
openr,1,'dumper'
readf,1,nux
readf,1,fnx
close,1
oplot,nux,fnx,COLOR=250

nu_BB = alog10(2.998d8/4400.0d-10)
print, 'BB = ',INTERPOL(fnx,nux,[nu_BB])

setup_plot,'idltemp.ps',/PS,/CLOSE

end


pro richdump, ALL=ALL,BLUE=BLUE,RED=RED,OPTLUM=OPTLUM,OPTDIM=OPTDIM, $
									IRLUM=IRLUM,IRDIM=IRDIM,BB=BB,SX=SX,HX=HX
setup_plot,'idltemp.ps',/PS

	OPENR,1,return_idl_routines_homedir(0)+'/agn_spectrum/richards_bol_corr_tables.dat'
	x = '1 2'
	y = STRSPLIT(x,/EXTRACT)
	i_line = 1
	i_line_max = 227
	nu = fltarr(i_line_max-1)
	f1 = nu
	f2 = nu
	f3 = nu
	f4 = nu
	f5 = nu
	f6 = nu
	f7 = nu
	while(i_line LT i_line_max) do begin
		READF,1,x
		y = STRSPLIT(x,/EXTRACT)
		nu[i_line-1] = FLOAT(y[0])
		f1[i_line-1] = FLOAT(y[1])
		f2[i_line-1] = FLOAT(y[2])
		f3[i_line-1] = FLOAT(y[3])
		f4[i_line-1] = FLOAT(y[4])
		f5[i_line-1] = FLOAT(y[5])
		f6[i_line-1] = FLOAT(y[6])
		f7[i_line-1] = FLOAT(y[7])
		i_line = i_line + 1
	endwhile
	CLOSE,1
	nuLnu = f1
		if (KEYWORD_SET(ALL)) then 		nuLnu = f1
		if (KEYWORD_SET(BLUE)) then 	nuLnu = f2
		if (KEYWORD_SET(RED)) then 		nuLnu = f3
		if (KEYWORD_SET(OPTLUM)) then 	nuLnu = f4
		if (KEYWORD_SET(OPTDIM)) then 	nuLnu = f5
		if (KEYWORD_SET(IRLUM)) then 	nuLnu = f6
		if (KEYWORD_SET(IRDIM)) then 	nuLnu = f7

	nuB = alog10(6.818d14)
	LB  = INTERPOL(nuLnu,nu,[nuB])
	print, nuB, LB

;	openw,1,'dumper'
;	printf,1,nu,FORMAT='(F5.2)'
;	printf,1,nuLnu,FORMAT='(F5.2)'
;	close,1
	
	nu_grid = 10. + 0.001*findgen((22.-10.)/0.001+1.)
	dlognu = 0.001
	nufnu = INTERPOL(nuLnu,nu,nu_grid)
	L_bol = TOTAL(10^(DOUBLE(nufnu)) * dlognu * alog(10.))
	print, ' GTR Lbol = ',alog10(L_bol)
	
	bb = INTERPOL(nufnu,nu_grid,[alog10(2.998d8/(4400.0d-10))])
	print, ' BB -> ',bb-alog10(L_bol)
	ir = INTERPOL(nufnu,nu_grid,[alog10(2.998d8/(15.0d-6))])
	print, ' IR -> ',ir-alog10(L_bol)

	keV = 2.418d17
	nsx = 0.5 + 0.001*findgen((2.0-0.5)/0.001+1.0)
		sx = INTERPOL(nufnu,nu_grid,alog10(nsx*keV))
		lsx = TOTAL(10^(DOUBLE(sx)) * (0.001/nsx))
			print, ' SX -> ',alog10(lsx)-alog10(L_bol)

	hx = INTERPOL(nufnu,nu_grid,alog10([2.998d8/(4400.0d-10)]))
	nhx = 2.0 + 0.001*findgen((10.0-2.0)/0.001+1.0)
		hx = INTERPOL(nufnu,nu_grid,alog10(nhx*keV))
		lhx = TOTAL(10^(DOUBLE(hx)) * (0.001/nhx))
			print, ' HX -> ',alog10(lhx)-alog10(L_bol)
	
	
	print, n_elements(nu)
	print, n_elements(nuLnu)

;; unnecessary 20,000 point binning
	homedir = return_idl_routines_homedir(0)+'/agn_spectrum/'
	OPENR,7,homedir+'xr.spec.dat'
	READF,7,n_to_read
	nu=fltarr(n_to_read)
	lnu=fltarr(n_to_read)
	READF,7,nu
	READF,7,lnu
	CLOSE,7
	
	nu0 = 16. + 0.0001*findgen((21.5-16.0)/0.0001)
	l0  = INTERPOL(lnu,nu,nu0)
	nu0 = 16. + 0.001*findgen((23.0-16.0)/0.001)
	l0  = INTERPOL(lnu,nu,nu0)
		nu_x_save = nu0
		l_x_save  = l0
	
	nu2  = alog10(2.0*2.418d17)
	nu10 = alog10(10.0*2.418d17)
		ok = where((nu0 GE nu2) AND (nu0 LE nu10))
		lok = l0[ok]
		ltot = TOTAL(10^(lok) * 0.0001 * alog(10.))
		print, alog10(ltot)
	
	x = 10^(interpol(l0,nu0,alog10(2.0*2.418d17)))
	
	print, 'x=',x/ltot

;	tot1 = alog10(TOTAL(10^(lnu)*

	nu0 = 16. + 0.02*findgen((21.5-16.0)/0.02)
	l0  = INTERPOL(lnu,nu,nu0)
	print, -2.0+alog10(EXP(-10^(21.48-20.1204)))

	plot,nu0,l0,ystyle=1,yrange=[-10.,0.]
		oplot,nu,lnu,color=250
		oplot,nu,-2.0+alog10(EXP(-10^(nu-20.13))),COLOR=250
	
;	openw,1,'dumper'
;	printf,1,nu0,FORMAT='(F5.2)'
;	printf,1,l0,FORMAT='(F5.2)'
;	close,1
	print, n_elements(nu0)




		OPENR,1,return_idl_routines_homedir(0)+'/agn_spectrum/sx_atten_fac.dat'
		N_NH = 0
		READF,1,N_NH
			NH = fltarr(N_NH)
			f  = fltarr(N_NH)
		READF,1,NH
		READF,1,f
		CLOSE,1



		OPENR,1,return_idl_routines_homedir(0)+'/agn_spectrum/hx_atten_fac.dat'
		N_NH = 0
		READF,1,N_NH
			NH = fltarr(N_NH)
			f  = fltarr(N_NH)
		READF,1,NH
		READF,1,f
		CLOSE,1
	bad = where(FINITE(f) NE 1,n_bad)
	if (n_bad GT 0) then begin
	good = where(FINITE(f) EQ 1)
		f[bad] = INTERPOL(f[good],NH[good],NH[bad])
	endif
;	openw,1,'dumper'
;	printf,1,NH,FORMAT='(F5.2)'
;	printf,1,f,FORMAT='(F14.8)'
;	close,1
	print, n_elements(f)



;; play w. the marconi spectrum

	;; start w. some L-something grid
	nu_grid = 10. + 0.01*findgen(1301)
	clight = 2.998d8
	nu_25 = alog10(clight/(2500.0d-10))
		l_2500_grid = 7. + 0.1*findgen(101) - nu_25
	nu_micron = alog10(clight/(1.0d-6))
	nu_13 = alog10(clight/(1300.0d-10))
	nu_12 = alog10(clight/(1200.0d-10))
	nu_5  = alog10(clight/(500.0d-10))
	nu_BB = alog10(clight/(4400.0d-10))
	nu_IR = alog10(clight/(15.0d-6))
	
	;; template left of 500 angstroms
		l_25 = 0.0
		l_nu_micron_to_13 = -0.44*(nu_grid-nu_25(0)) + l_25(0)	
		l_micron = INTERPOL(l_nu_micron_to_13,nu_grid,[nu_micron])
		l_13 = INTERPOL(l_nu_micron_to_13,nu_grid,[nu_13])
		l_nu_inf_to_micron = 2.0*(nu_grid-nu_micron(0)) + l_micron(0)
		l_nu_13_to_12 = 0.0*(nu_grid-nu_13(0)) + l_13(0)
		l_12 = l_13(0)
		l_nu_12_to_5 = -1.76*(nu_grid-nu_12(0)) + l_12(0)


		l_low = l_nu_inf_to_micron*(nu_grid LT nu_micron(0)) + $
				l_nu_micron_to_13*((nu_grid GE nu_micron(0)) AND (nu_grid LT nu_13(0))) + $
				l_nu_13_to_12*((nu_grid GE nu_13(0)) AND (nu_grid LT nu_12(0))) + $
				l_nu_12_to_5*((nu_grid GE nu_12(0)) AND (nu_grid LT nu_5(0)))
		nulnu_low = l_low + nu_grid
			nulnu_25 = INTERPOL(nulnu_low,nu_grid,[nu_25])
			nulnu_BB = INTERPOL(nulnu_low,nu_grid,[nu_BB])
			nulnu_IR = INTERPOL(nulnu_low,nu_grid,[nu_IR])
			nulnu_mu = INTERPOL(nulnu_low,nu_grid,[nu_micron])
			 fbb=nulnu_BB-nulnu_25
			 fir=nulnu_IR-nulnu_25
			 fmu=nulnu_mu-nulnu_25
				fbb=-fbb
				fir=-fir
				fmu=-fmu


		l_hig = INTERPOL(l_x_save,nu_x_save,nu_grid) - nu_grid
			;; renorm to 1.0 at 2keV
			keV=2.418d17
			nu_2k = alog10(2.0*keV)
			nu_1k = alog10(keV)
			lnorm = INTERPOL(l_hig,nu_grid,[nu_2k])
			l_hig = l_hig - lnorm(0)
		l_2keV_grid = 0.72 * (l_2500_grid+alog10(3.9)+33.0) + 4.53 - alog10(3.9)-33.0

		plot,nu_grid,l_low+nu_grid
		plot,nu_grid,l_hig+nu_grid
		plot,l_2500_grid+nu_25,l_2keV_grid+nu_2k

	lbols = 0.0*l_2500_grid
	l500 = 0.0*l_2500_grid
	l1k = 0.0*l_2500_grid
	lsxs=l1k
	lhxs=l1k
	for i=0, n_elements(l_2500_grid)-1 do begin
		l_25 = l_2500_grid[i]
		l_2k = l_2keV_grid[i]
		l_opt = l_low + l_25
		l_xry = l_hig + l_2k
		l_int = INTERPOL([l_25,l_2k],[nu_25,nu_2k],nu_grid)
		l_tot = l_opt*(nu_grid LT nu_5) + $
				l_int*((nu_grid GE nu_5) AND (nu_grid LE nu_1k)) + $
				l_xry*(nu_grid GT nu_1k)
		l_tot = l_tot + nu_grid
		 ;plot,nu_grid,l_tot-(l_25+nu_25),xrange=[12.,22.],xstyle=1,ystyle=1,yrange=[-10.,0.5]
		
		d_log_nu = 0.01
		d_ln_nu = d_log_nu * alog(10.)
		nuLnu = l_tot
			ok = where((FINITE(nuLnu) EQ 1) AND (FINITE(nuLnu,/NAN) EQ 0))
		l500[i] = INTERPOL(nuLnu,nu_grid,[alog10(2.998d8/(500.0d-10))])
		l1k[i]  = INTERPOL(nuLnu,nu_grid,[alog10(2.418d17)])
		Lbol  = TOTAL(10^(DOUBLE(nuLnu[ok])) * d_ln_nu)
		lbols[i] = alog10(Lbol)


	keV = 2.418d17
	nsx = 0.5 + 0.001*findgen((2.0-0.5)/0.001+1.0)
		sx = INTERPOL(nuLnu,nu_grid,alog10(nsx*keV))
		lsx = TOTAL(10^(DOUBLE(sx)) * (0.001/nsx))
			;print, ' SX -> ',alog10(lsx)-alog10(L_bol)

	nhx = 2.0 + 0.001*findgen((10.0-2.0)/0.001+1.0)
		hx = INTERPOL(nuLnu,nu_grid,alog10(nhx*keV))
		lhx = TOTAL(10^(DOUBLE(hx)) * (0.001/nhx))
			;print, ' HX -> ',alog10(lhx)-alog10(L_bol)
	
		lsxs[i]=alog10(lsx)
		lhxs[i]=alog10(lhx)

	endfor
	l25 = l_2500_grid + nu_25(0)
	plot,l25,lbols
	plot,lbols,l25
		xx = lbols - 12.0
		yy = -(l25-lbols)
	plot,xx+12.0,yy
		P_guess = [1.,0.5,0.2,-0.5]
		P_fitted = MPFITFUN('cubic_fitfun',xx,yy,0.0*xx+1.,P_guess,/QUIET)	
		print, P_fitted
		print, P_fitted[0]*10^fbb[0], P_fitted[1]*10^fbb[0], P_fitted[2], P_fitted[3]
		print, P_fitted[0]*10^fir[0], P_fitted[1]*10^fir[0], P_fitted[2], P_fitted[3]
		print, P_fitted[0]*10^fmu[0], P_fitted[1]*10^fmu[0], P_fitted[2], P_fitted[3]
		y0 = cubic_fitfun(xx,P_fitted)
		oplot,xx+12.0,y0,COLOR=250

	plot,l500,l1k
		P_guess = [1.,0.0,0.0,-0.5]
		;P_guess = [1.,-0.5]
		P_fitted = MPFITFUN('cubic_fitfun',l500,l1k,0.0*l500+1.,P_guess,/QUIET)	
		y0 = cubic_fitfun(l500,P_fitted)
		oplot,l500,y0,COLOR=250
		print, P_fitted




		xx = lbols - 12.0
		yy = -(lsxs-lbols)
	plot,xx+12.0,yy
		P_guess = [1.,0.5,0.2,-0.5]
		P_fitted = MPFITFUN('cubic_fitfun',xx,yy,0.0*xx+1.,P_guess,/QUIET)	
		print, P_fitted
		y0 = cubic_fitfun(xx,P_fitted)
		oplot,xx+12.0,y0,COLOR=80,linestyle=2,thick=6.0
		xx = lbols - 12.0
		yy = -(lhxs-lbols)
	plot,xx+12.0,yy
		P_guess = [1.,0.5,0.2,-0.5]
		P_fitted = MPFITFUN('cubic_fitfun',xx,yy,0.0*xx+1.,P_guess,/QUIET)	
		print, P_fitted
		y0 = cubic_fitfun(xx,P_fitted)
		oplot,xx+12.0,y0,COLOR=50,linestyle=2,thick=6.0


setup_plot,'idltemp.ps',/PS,/CLOSE
end


function cubic_fitfun, x, P
	;return, P[0]+P[1]*x+P[2]*x*x+P[3]*x*x*x
	return, alog10(P[0]*(10^(x*P[3])) + P[1]*(10^(x*P[2])))
	;return, P[0] + P[1]*x
end

