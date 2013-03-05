;------------------------------------------------------------------------
;
;
;
;------------------------------------------------------------------------


pro process_cumulative_profile, x, y, $
				xmin, xmax, $
				bins, $
                                xs, ys, $         ; returned values
                                xlog=xlog

xs = fltarr(bins)
ys= fltarr(bins)

tempxmax= xmax
tempxmin= xmin
if keyword_set(xlog) then begin
        ;tempxmin= alog10(xmin)
        ;tempxmax= alog10(xmax)
        x= alog10(x)
endif

binsize = float((tempxmax-tempxmin))/bins

for i=0, bins-1 do begin 
    currentr= tempxmin + (i+1)*binsize
    idx= where(x le currentr)
    xs[i]= currentr 
    if idx(0) ne -1 then ys[i]= total(y(idx)) else ys[i]= 0.0
endfor



end



;==================================================================================


