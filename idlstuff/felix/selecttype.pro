
pro selecttype,type,npart,pos,vel,mass,outpos,outvel,outmass,outn

; This routine selects the "type" particles out of the total particles

   if (type eq 10) then begin
      print,'SELECTING ALL STARS!!!'
      from = total(npart(0:2-1))
      to   = total(npart(0:4))           
   endif else begin
      if (type eq 11) then begin
      print,'SELECTING OLD STARS!!!'      
         from = total(npart(0:2-1))
         to   = total(npart(0:3))
      endif else begin      
         print,'Selecting type',type
	 if (type gt 0) then begin
            from    = total(npart(0:type-1))
	 endif else begin
            from    = 0
	 endelse   
	 to      = total(npart(0:type))
      endelse   
   endelse   
   
   outpos  = pos(0:2,from:to-1)
   outvel  = vel(0:2,from:to-1)
   outmass = mass (from:to-1)
   
   outn    = long(to-from)
   
end
