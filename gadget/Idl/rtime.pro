pro rtime,name
   dummy=bytarr((13+4+9+1)*8+(6+17+6)*4+20)
; if SFR_METALS: 
;   dummy=bytarr((13+4+9+1)*8+(6+17+6)*4+20+4)
   t=dblarr(4)
   openr,1,name
   readu,1,dummy,t
   close,1
   print,t
end
