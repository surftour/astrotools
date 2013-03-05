PRO mynormalize,x,count

   ; normalizing the quantity with the good count

   FOR i=0,3 DO BEGIN

      select1 = where(count(*,i) lt 0.5)
      select2 = where(count(*,i) gt 0.5)   

      ; for the undefined values set the array to zero
      if (min(select1) gt -0.5) then begin
	 x(select1,i) = 0.0
      endif   

      ; for the undefined values set the array to zero
      if (min(select2) gt -0.5) then begin
	 x(select2,i) = x(select2,i)/count(select2,i)
      endif   

   ENDFOR   

END
