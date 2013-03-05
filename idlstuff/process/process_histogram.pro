;------------------------------------------------------------------------
;
;    Generate a histogram from a variable
;
;
;
;
;------------------------------------------------------------------------



function process_histogram, histvar, $
				weighting=weighting, $
				xmax=xmax, xmin=xmin, $
				levels=levels, $
				bins=bins, $
				oplotit=oplotit, $
				nonorm=nonorm, $
				mannorm=mannorm, $
				normalization=normalization


; histogram is normalized 
; to one (by default)
; -------------------------

alltts= histvar


levels= levels*1.0
step= (xmax-xmin)/(levels)
bins= (IndGen(levels+1)/(levels+1)*(xmax-xmin)) + xmin
;bins=[bins,xmax]
bins=bins+(step*0.5)


if not keyword_set(weighting) then begin
	; std histogram
	hist_alltts= histogram(alltts, binsize=step, max=xmax, min=xmin)
endif else begin
	N= n_elements(alltts)

	hist_alltts= fltarr(levels)

	for i= 0L, N-1 do begin

		thisbin= (alltts(i) - xmin) / step

		if alltts(i) lt xmin then thisbin= 0
		if alltts(i) ge xmax then thisbin= levels-1
; for density check
;if alltts(i) gt 3.0 then print, "i= ", i, alltts(i), weighting(i)

		thisbin= long(thisbin)

		hist_alltts[thisbin] = hist_alltts[thisbin] + weighting(i)
	endfor
endelse


; make sure it extends to 0
bins= [xmin,bins,xmax]
hist_alltts= [hist_alltts(0),hist_alltts,hist_alltts(levels-1)]
help, bins, hist_alltts

normalization= float(max(hist_alltts))
print, "histogram normalization= ", normalization
if keyword_set(nonorm) then normalization= 1.0
if keyword_set(mannorm) then normalization= mannorm
print, "(final) histogram normalization= ", normalization
hist_alltts= hist_alltts/normalization
hist_alltts_area = total(hist_alltts*step)

;if keyword_set(oplotit) then begin
if oplotit ge 0 then begin

	oplot, bins, hist_alltts, psym=10, color=oplotit, thick=5.0

	; fill in histogram
	; ------------------
	nbins= bins+(step*0.5)          ; make x coord
	;nbins[0]= 0.0
	nbins[0]= xmin
	nbins=[nbins,nbins]
	nbins= nbins(sort(nbins)) 
	;nbins= nbins[nbins, xmax]

	ntts= fltarr(2.*levels + 2) 
	for i=1,levels do begin
	   ntts[2.*i-1]= hist_alltts[i]
	   ntts[2.*i]= hist_alltts[i]
	endfor

	if oplotit eq 20 then begin
		polyfill, nbins, ntts, /data, color= oplotit, /fill, linestyle=0, $
					thick=3.0
	endif

	;
	;polyfill, nbins, nctts, /data, color= 50, /line_fill, linestyle=0, $
	;				thick=3.0
	;polyfill, nbins, nctts, /data, color= 50, /fill, linestyle=0, $
	;				thick=3.0
	if oplotit eq 220 then begin
		polyfill, nbins, ntts, /data, color= oplotit, /fill, linestyle=0, $
					thick=3.0
	endif

	;polyfill, nbins, ntts, /data, color= 150, /line_fill, linestyle=0, $
	;				thick=3.0, orientation=90.0
	if oplotit eq 150 then begin
		polyfill, nbins, ntts, /data, color= oplotit, /line_fill, linestyle=0, $
					thick=3.0, orientation=45.0
	endif

endif



return, hist_alltts


end









;=================================================================================







