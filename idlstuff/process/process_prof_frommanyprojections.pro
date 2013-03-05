; --------------------------
; this is my master program which takes in
; x, y and z coordinates, plus some quantity
; and randomly projects this in many 
; different directions, averaging the results.
; it returns the average and 1 sigma 
; dispersion.
; --------------------------

; result_avg, result_1sig, xaxisvals - all need to be predefined as fltarr(bins)
;
; note that the default is to use the surface density routine,
; and if we set the /average flag, it will average instead
;
;

pro process_prof_frommanyprojections, x, y, z, quantitytoaverage, $
                                        quantitytoaverage_y=quantitytoaverage_y, $
                                        quantitytoaverage_z=quantitytoaverage_z, $
                                        center=center, $
                                        xmin, xmax, bins, $
                                        result_avg, result_1sig, xaxisvals, $
                                        x_is_log=x_is_log, x_is_devac=x_is_devac, $
                                        average=average, $
                                        y_weighting= y_weighting, $
                                        fitsersic=fitsersic, $
                                        sersic_n_avg=sersic_n_avg, $
                                        sersic_n_1sig=sersic_n_1sig, $
                                        exslt_avg=exslt_avg, $
                                        exslt_1sig= exslt_1sig,$
                                        r_e_avg=r_e_avg, $
                                        r_e_1sig=r_e_1sig, $
                                        dontoverplot=dontoverplot


;n_projections= 100    ; maybe do 1000 at some point
;n_projections= 40
n_projections= 20
seed= 154L

; 2D array which stores projected quantity for
; each projection
temp_average= fltarr(n_projections+1,bins)
sersic_n_i= fltarr(n_projections+1)
exslt_i= fltarr(n_projections+1)
r_e_i= fltarr(n_projections+1)

if keyword_set(center) then c=center else c=[0,0,0]

original_quantitytoaverage= quantitytoaverage
if keyword_set(quantitytoaverage_y) and keyword_set(quantitytoaverage_z) then begin
        qx= original_quantitytoaverage
        qy= quantitytoaverage_y
        qz= quantitytoaverage_z
endif


; ---------------------------
;  begin projection loop
; ---------------------------
for i=0,n_projections do begin

        rdphi= randomu(seed)
        rdtheta= randomu(seed)

        ; in radians
        theta= rdtheta*!PI
        phi= rdphi*2*!PI


        ; rotate
        rot_x= (x-c[0])*(cos(theta)*cos(phi)) + (y-c[1])*(cos(theta)*sin(phi)) + (z-c[2])*sin(theta)
        rot_y= -(x-c[0])*sin(theta) + (y-c[1])*cos(theta)
        ;rot_z= x*(sin(theta)*cos(phi)) - y*(sin(theta)*sin(phi)) + z*cos(theta)
        ;rot_vx= vx*(cos(theta)*cos(phi)) + vy*(cos(theta)*sin(phi)) + vz*sin(theta)
        ;rot_vy= -vx*sin(theta) + vy*cos(theta)
        ;rot_vz= vx*(sin(theta)*cos(phi)) - vy*(sin(theta)*sin(phi)) + vz*cos(theta)

        rot_r= sqrt(rot_x*rot_x + rot_y*rot_y)

        if keyword_set(quantitytoaverage_y) and keyword_set(quantitytoaverage_z) then begin
            quantitytoaverage= qx*(sin(theta)*cos(phi)) - qy*(sin(theta)*sin(phi)) + qz*cos(theta)
        endif

        if keyword_set(x_is_log) then rot_r=alog10(rot_r)
        if keyword_set(x_is_devac) then rot_r=rot_r^(0.25)

        ; calculate the profile for this projection
        rad= fltarr(bins)
        thisquant= fltarr(bins)


        if not keyword_set(average) then begin
                process_prof_sd, rot_r, quantitytoaverage, bins, xmax, xmin, $
				rad, thisquant, $
                                x_is_log=x_is_log, x_is_devac=x_is_devac, $
                                sd_1sig=mass_sd_1sig, $
                                y_weighting=y_weighting
        endif else begin
		; ***
		; not currently working
		; ***
                process_prof_avg, rot_r, quantitytoaverage, bins, xmax, xmin, $
				rad, thisquant, $
				avg_1sig=avg_1sig, $
                                y_weighting=y_weighting
        endelse

        temp_average[i,*]= thisquant

        xaxisvals= rad

        ; try to fit Sersic Profile to it
        ; --------------------------------
        if keyword_set(fitsersic) then begin

                r_tofit= rad^(4.0)          ; assumes x_is_devac is on
                sd_tofit= thisquant
                weight= mass_sd_nlog_1sig

                ;minfitradius= 0.5
                minfitradius= 1.0

                fitandoverplot_sersicprofile, r_tofit, sd_tofit, weight, $ 
                                                minfitradius=minfitradius, $
                                                sersic_n=sersic_n, $
                                                excesslight=excesslight, $
                                                r_e=r_e, $
                                                dontoverplot=dontoverplot

                sersic_n_i[i]= sersic_n
                exslt_i[i]= excesslight
                r_e_i[i]= r_e

        endif


        if (i mod 10) eq 0 then print, "i= ",i

endfor

for i=0,bins-1 do begin
        result_moment= moment(temp_average[*,i])
        result_avg[i]= result_moment[0]
        result_1sig[i]= sqrt(result_moment[1])
endfor


sersic_n_avg= mean(sersic_n_i)
sersic_n_1sig= sqrt(variance(sersic_n_i))
exslt_avg= mean(exslt_i)
exslt_1sig= sqrt(variance(exslt_i))
r_e_avg= mean(r_e_i)
r_e_1sig= sqrt(variance(r_e_i))


end


