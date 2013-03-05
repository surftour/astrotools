; ==================================================================

pro readsnapshot,basename,num,$
                 N,npart,massarr,time,redshift,$
		 pos,vel,mass,$
		 verbose=verbose

; ==================================================================
;
; This routine reads a snapshot of gadget including gas, stars etc.
;
; It seems that taking care of the endianness is not necessary???
;
;
; I have NO idea how the structure of the ids is. It seems wrong to
; me. At least it is for sure, that the claim that the respective
; particle types have indices that are assured to ly in the ranges
; of the respective particle numbers.
;
; ------------------------------------------------------------------

   if (keyword_set(verbose)) THEN verbose = 1 ELSE verbose = 0  
   
   close,/all

   exts = '000'
   exts = exts+strcompress(string(num),/remove_all)
   exts = strmid(exts,strlen(exts)-3,3)
   f    = basename+"snapshot_"+exts
   f    = strcompress(f,/remove_all)

   npart          = lonarr(6)
   massarr        = dblarr(6)
   time           = 0.0D
   redshift       = 0.0D
   flag_sfr       = 0L
   flag_feedback  = 0L

   bytesleft      = 256-6*4 - 6*8 - 8 - 8 - 2*4

   la             = intarr(bytesleft/2)

   print,f
   openr,1,f,/f77_unformatted

   readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,la

   N    = total(npart)

   pos  = fltarr(3,N)
   vel  = fltarr(3,N)
   id   = lonarr(N)   
   mass = fltarr(N)
   wrap = fltarr(total(npart(where(abs(massarr)  lt 1e-8))))

   if (verbose eq 1) then print,"reading the pos"
   readu,1,pos
   
   if (verbose eq 1) then print,"reading the vel"
   readu,1,vel
   
   if (verbose eq 1) then print,"read the id"
   readu,1,id  
   
   if (verbose eq 1) then print,"read the mass"   
   readu,1,wrap
   
   ; make mass array
   index  = 0L
   windex = 0L
   print,npart
   FOR i=0,4 DO BEGIN
      IF (npart(i) gt 0) then begin
	 IF (abs(massarr(i)) lt 1e-8) THEN BEGIN
            mass(index:index+npart(i)-1) = wrap(windex:windex+npart(i)-1)
	    windex = windex + npart(i)
	    index  = index  + npart(i)
	 ENDIF ELSE BEGIN
            mass(index:index+npart(i)-1) = massarr(i)
	    index = index + npart(i)         
	 ENDELSE       
      endif
   ENDFOR
   wrap=0L
   
   NGas   = npart(0)
   NHalo  = npart(1)
   NDisk  = npart(2)
   NBulge = npart(3)
   NStars = npart(4)

   if (verbose eq 1) then begin
      print,"Ngas",Ngas
      print,"NHalo",NHalo
      print,"NDisk",NDisk
      print,"NBulge",NBulge
      print,"NStars",NStars
   end
  close,1

end
