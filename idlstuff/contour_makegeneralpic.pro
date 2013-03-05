pro contour_makegeneralpic, x, y, xmax, xmin, ymax, ymin, $
                        pixels=pixels, $
                        NxNImage=NxNImage



if keyword_set(pixels) then XYPIXELS= long(pixels) else XYPIXELS= 200L


NxNImage= fltarr(XYPIXELS, XYPIXELS)
NxNImage= NxNImage*0.0 + 1.0


if not keyword_set(x) then begin
   print, "  "
   print, "PROBLEM: contour_makegeneralpic"
   print, "  "
   print, "  "
   return
endif


        

;  map x and y to 0-1 (i think)
; -------------------------------------
nbins= XYPIXELS

idx=where((x ge xmin) and (x le xmax))
if idx(0) ne -1 then begin
	xtobin_1=x(idx)
	ytobin_1=y(idx)
endif else begin
	xtobin_1=x
	ytobin_1=y
endelse
xtobin_1=xtobin_1-xmin
xtobin_1=xtobin_1*nbins/(xmax-xmin)
xtobin_1=fix(xtobin_1)

idx=where((ytobin_1 ge ymin) and (ytobin_1 le ymax))
if idx(0) ne -1 then begin
	ytobin=ytobin_1(idx)
	xtobin=xtobin_1(idx)
endif else begin
	ytobin=ytobin_1
	xtobin=xtobin_1
endelse
ytobin=ytobin-ymin
ytobin=ytobin*nbins/(ymax-ymin)
ytobin=fix(ytobin)


phist= hist_2d(xtobin, ytobin, max1=nbins-1, max2=nbins-1)


; for phase diagram
; ------------------
;phist=phist*phist           ; make stretch ^2
idx= where(phist EQ 0)
;phist[idx]= 1e8
;phist= alog10(phist)
;phist[idx]= 0
;phist= phist*236/max(phist)
;print, "max(phist)= ", max(phist)
phist= phist*236/50.0
print, " "
print, " WARNING: in contour_makegeneralpic, fixed map scale"
print, " "
phist= phist+20

; for v/sigma plot
; -------------------------
;phist=phist*phist
;;phist=phist*phist*phist
;idx= where(phist EQ 0)
;phist= phist+20
;phist= phist*235/max(phist)      ; this is std
;phist= phist+2                   ; this is std
;phist= phist+19
;phist= phist*235/800.0
;phist= phist*235/120.0
;phist= phist*235/20.0


; for shells plot
; -------------------------
;phist=alog10(phist)
;idx= where(phist EQ 0)
;phist= phist+20
;phist= phist*235/max(phist)
;phist= phist*235/10.0



; for j versus age plot
; -------------------------
;idx= where(phist le 0)
;phist(idx)= 1e+8
;phist=100*alog10(phist)
;phist(idx)= 0
;maxphist= max(phist)
;print, "Max(phist)= ", max(phist)
;print, "Min(phist)= ", min(phist)
;phist= phist+20
;phist= phist*235/maxphist



; set zero's to white
phist[idx]= 1


ind=where(phist ge 256)
if ind(0) ne -1 then phist(ind)=254

ind=where(phist le 0)
if ind(0) ne -1 then phist(ind)=1


; invert colors if we want
; ------------------------
invertcolors= 1
if invertcolors eq 1 then begin                  
        phist=256-phist                          ; invert color table
;        idx= where(phist EQ 254)               ; set background to white
endif else begin
;        idx= where(phist EQ 2)
endelse

if idx(0) ne -1 then phist(idx)= 1       ; this is setting it to white
;if idx(0) ne -1 then phist(idx)= 0       ; this sets background to black



; -----------------------------
;   set image and we're done
; -----------------------------

    NxNImage= phist


; -------------
;  Done
; -------------



end




