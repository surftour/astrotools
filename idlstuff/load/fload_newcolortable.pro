pro fload_newcolortable, newtab

COMMON PlotInfo

;if not keyword_set(newtab) then newtab= 4

loadct, newtab
tvlct,rrr,ggg,bbb,/get
v1=[0,255]
v2=[0,255]
v3=[0,255]
tvlct,v1,v2,v3,0  ; define colors 0 and 1 as black and white

if newtab ne ct_current then ct_lastdiff= ct_current
ct_last= ct_current
ct_current= newtab

end



