PRO cleanup,x,y,count
   ; this routine cuts the arrays such that there is no interior bin with 0 particles in it  

   ; now look for the first occurence of a zero in the count array and discard all
   ; bins with larger x
   discard = min(where(count lt 0.5))
   IF (discard gt 0) THEN BEGIN
      y(discard:n_elements(y)-1)         = 0.0
      count(discard:n_elements(count)-1) = 0.0
   ENDIF
   
end   
