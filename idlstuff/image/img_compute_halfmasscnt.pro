pro img_compute_halfmasscnt, Map, HalfMassContourLevel, $
				Pixel_Area=Pixel_Area, $
				Total_Mass=Total_Mass, $
				FindMassFactor= FindMassFactor


if not keyword_set(FindMassFactor) then FindMassFactor= 0.5

; -----------------------------------------------------
; now find the mass contour of Total * FindMassFactor
; -----------------------------------------------------

HalfMassContourLevel= -1.0


    ; Need these to be passed into here
    ; Total_Mass
    ; Pixel_Area

    ; determine half-mass surface density
    sortedMap= Map(sort(Map))
    Nm= n_elements(sortedMap)
    print, "Nm= ", Nm

    ; Need these to be passed into here
    ; Total_Mass
    ; Pixel_Area
    Plotted_Mass= Pixel_Area * total(Map)

    print, "Pixel_Area= ", Pixel_Area
    print, "Total Mass= ", Total_Mass
    print, "Plotted_Mass= ", Pixel_Area * total(Map)
    print, "Total_Mass*FindMassFactor= ", Total_Mass*FindMassFactor


    print, " ************************************** "
    print, "     Determine HalfMassContourLevel  "
    print, " "

    if (Plotted_Mass lt (Total_Mass*FindMassFactor)) then FindMassFactor= FindMassFactor / 2.0
    if (Plotted_Mass lt (Total_Mass*FindMassFactor)) then FindMassFactor= FindMassFactor / 2.0

;----------------------------------------------------
;
;  Old brute force method
;
;    i= 0L
;    repeat begin
;	enclosedmass= Pixel_Area*total(sortedMap[Nm-1-i:Nm-1])
;;;print, "WARNING: fixing enclosedmass by hand"
;;;enclosedmass= 10.0e+12
;	i= i+1
;	if (i gt (Nm-2)) then enclosedmass=Total_Mass
;	;if (Plotted_Mass lt (Total_Mass*FindMassFactor)) then enclosedmass=Total_Mass
;	;if skiphalfmasscontour eq 1 then enclosedmass=Total_Mass
;	if Total_Mass lt 0 then break
;	if i ge n_elements(sortedMap * 0.8) then break
;    endrep until enclosedmass ge (Total_Mass*FindMassFactor)
;print, "i= "
;print, "enclosedmass= ", enclosedmass
;print, "Total_Mass*FindMassFactor= ", Total_Mass*FindMassFactor
;print, "sortedMap[Nm-1-i-1]= ", sortedMap[Nm-1-i-1]
;print, " ************************************** "
;print, "     Determine HalfMassContourLevel (now do the newer version below) "
;print, " "
;----------------------------------------------------
;
;  New bisector method
;
    i= long(Nm)-1
    di= i
    dm= 100.0
print, i, di, dm
    repeat begin
	i = i - (dm/abs(dm)) * di/2
	di= di/2
;print, "i, di= ", i, di

	enclosedmass= Pixel_Area*total(sortedMap[Nm-1-i:Nm-1])

	dm= enclosedmass - (Total_Mass*FindMassFactor)

;print, "dm, enclosedmass, tm*fmf= ", dm, enclosedmass, Total_Mass*FindMassFactor
       ; trap for problem situations
       if (i gt (Nm-2)) then enclosedmass=Total_Mass
       if Total_Mass lt 0 then break
       ;if i ge n_elements(sortedMap * 0.8) then break

	; check if we're done
	done= 0

	if di le 0 then done= 1

	if (abs(dm)/(Total_Mass*FindMassFactor)) lt 0.001 then done= 1    ; accuracy better than 0.1%

    endrep until done eq 1

;----------------------------------------------------


    if (Nm-1-i-1 ge 0) then HalfMassContourLevel= sortedMap[Nm-1-i-1] else HalfMassContourLevel=sortedMap[Nm-2]

    print, " "
    print, "  i= ", i
    print, "  enclosedmass= ", enclosedmass
    print, "  FindMassFactor= ", FindMassFactor
    print, "  HalfMassContourLevel= ",HalfMassContourLevel
    print, " "


    print, " ************************************** "

; -------------
;  Done
; -------------



end


