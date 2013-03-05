
; =============================================================================
; =============================================================================


;
;  Used by plotting routine to scale the 
;  images.
; ------------------------------------------------
function ScalePic, Pic, maxvalue, minvalue, $
			var_already_inlog=var_already_inlog

    Map= Pic
    ma= maxvalue+maxvalue*0.1
    mi= minvalue-minvalue*0.1

    ind=where(Map lt mi) 
    if ind(0) ne -1 then Map(ind)=mi
    ind=where(Map gt ma) 
    if ind(0) ne -1 then Map(ind)=ma


    Cols=255 ; number of colors

    if keyword_set(var_already_inlog) then begin
        Pic=(Map-mi)/(ma-mi) * (cols-3) + 2
    endif else begin
	print, "Taking log of Map"
        Pic=(alog10(Map)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2
    endelse

    ind=where(Pic ge 256)
    if ind(0) ne -1 then Pic(ind)=255

    ind=where(Pic le 0)
    if ind(0) ne -1 then Pic(ind)=1


    invertcolors= 1
    if invertcolors eq 1 then begin
        Pic=256-Pic                          ; invert color table
        idx= where(Pic EQ 254)               ; set background to white
    endif else begin
        idx= where(Pic EQ 2)
    endelse

    if idx(0) ne -1 then Pic(idx)= 1       ; this is setting it to white
    ;if idx(0) ne -1 then Pic(idx)= 0       ; this sets background to black

   NewScaledPic= Pic

   return, NewScaledPic

end








; =============================================================================


;   Image the grid we just created


;  -------------------
;  |        |        |
;  | Vel(1) |Density |
;  |        |        |
;  |        |        |
;  --------------------
;  |        |        |
;  | Gas    |Stellar |
;  | Temp   |  Mass  |
;  |        |        |
;  --------------------



pro image_grid, junk

;------
frun='/raid4/tcox/ds/vc3vc3e_2'
;------
imgbase= 'gridimage'
imgdir= frun+'/desikagrids_img/'
;------


read_


; -------------------------------------------------------


;
Ncubed= n_elements(Grid)/3.0
N= long(Ncubed^(1/3.0))
Pic= fltarr(N,N)


zvalue= min(abs(Grid2(2,*)))
zvalue= idx(0)
;------------------------------
idx=where(Grid(2,*) eq zvalue)
Density=Density(idx)
Temp=Temp(idx)
StellarMass=StellarMass(idx)
;------------------------------
;stop
;zs= Grid(2,*)
;idx=where(zs gt zvalue)
;zvalue= min(zs(idx))
;idx=where(zs eq zvalue)
;Density=Density(idx)
;Temp=Temp(idx)
;StellarMass=StellarMass(idx)


;
; The following are already in LOG!!!!
;

; Density
; ------------
for i=0,N-1 do Pic(i,*)= Density(i*N:(i+1)*N-1)
maxtemp=max(Pic)
mintemp=min(Pic)
print, "Density: max/min ", maxtemp, mintemp
Pic= ScalePic(Pic,maxtemp,mintemp,/var_already_inlog)


; Temp
; ------------
for i=0,N-1 do Pic(i,*)= Temp(i*N:(i+1)*N-1)
maxtemp=max(Pic)
mintemp=min(Pic)
print, "Temp: max/min ", maxtemp, mintemp
maxtemp=6.5
mintemp=0.0
print, "setting pic scale to: max/min ", maxtemp, mintemp
Pic= ScalePic(Pic,maxtemp,mintemp,/var_already_inlog)



; StellarMass
; ------------
for i=0,N-1 do Pic(i,*)= StellarMass(i*N:(i+1)*N-1)
maxtemp=max(Pic)
mintemp=min(Pic)
print, "StellarMass: max/min ", maxtemp, mintemp
Pic= ScalePic(Pic,maxtemp,mintemp,/var_already_inlog)




end





; =============================================================================




