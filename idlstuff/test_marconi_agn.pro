
n_spectrum = 1000
l_bol      = 10.1334
log_nu_min = 14.0
log_nu_max = 15.0
redshift = 0.0


;x = marconi_spectrum(l_bol,redshift,n_spectrum,log_nu_min,log_nu_max)
x = fload_marconi_agn_spectrum(l_bol,redshift,n_spectrum,log_nu_min,log_nu_max)
nu   = fltarr(n_spectrum)
ll   = fltarr(n_spectrum)
l_nu = fltarr(n_spectrum)
l_l  = fltarr(n_spectrum)
nu     = x(0:n_spectrum-1)
l_nu   = x(n_spectrum:2*n_spectrum-1)

cc = 3.0e14
ll = alog10(cc/10.0^nu)
l_l = alog10(cc*(10.0^l_nu)/((10.0^ll)^2))
window,3,xsize=400,ysize=400
;plot,nu,(10^(l_nu))*(10.0^(nu)),xrange=[12.0,21.0],xstyle=1,/ylog,yrange=[10.0^(9.5),10.0^(12.0)],ystyle=1
plot,nu,(10^(l_nu)),xrange=[10.0,21.0],xstyle=1,/ylog,yrange=[10.0^(-10.),10.0^(2.0)],ystyle=1
;plot,nu,(10^(l_bol+l_nu)),xrange=[10.0,21.0],xstyle=1,/ylog,yrange=[10.0^(-10.),10.0^(2.0)],ystyle=1
;plot,nu,(10^(l_nu)/10^(l_bol)),xrange=[10.0,21.0],xstyle=1,/ylog,yrange=[10.0^(-10.),10.0^(2.0)],ystyle=1



; TJ's tests
;-------------------------------

; can i recover the total luminosity?
print, "  "
print, "can we recover the total luminosity"
for i=1,n_spectrum-1 do begin
	nu_0= 10^(nu[i-1])
	nu_1= 10^(nu[i])

	delta_nu= nu_1 - nu_0

	; l_nu
	l_nu_0= 10^(l_nu[i-1])
	l_nu_1= 10^(l_nu[i])
	avg_l_nu= 0.5*(l_nu_0 + l_nu_1)

	; nu * l_nu
	nu_l_nu_0= nu_0*10^(l_nu[i-1])
	nu_l_nu_1= nu_1*10^(l_nu[i])
	avg_nu_l_nu= 0.5*(nu_l_nu_0 + nu_l_nu_1)

	totalLum= avg_l_nu * delta_nu
endfor
print, "Total_lum= ",totalLum
total_l_nu = total(10^(l_nu))
print, "Total l_nu= ", total_l_nu
total_nu_l_nu = total(10^(l_nu)*10^(nu))
print, "Total nu_l_nu= ", total_nu_l_nu
print, "  "


; b-band
l_lambda= 4400.0
print, "lambda= ",l_lambda
nu_lambda= cc/l_lambda * 1.0e+4    ; 10^4 converts to microns
print, "nu_lambda ", nu_lambda

print, " B-band luminosity"
idx=where(10^(nu) ge nu_lambda)
print, 10^(l_nu(idx(0)))*10^(nu(idx(0)))
print, " "


end
