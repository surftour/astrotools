;
;
;Predefined IDL symbols
;
;PSYM Value       Plotting Symbol
;----------       ------------------
;1 	          Plus sign (+)
;2 	          Asterisk (*)
;3 	          Period (.)
;4 	          Diamond
;5 	          Triangle
;6 	          Square
;7 	          X
;8 	          User-defined. See USERSYM procedure.
;9 	          Undefined
;10 	          Histogram mode.
;
;
;
;
;
;
; routine to select points
; ---------------------------
pro select_thispoint, pointselection, thispsym, thiscolor, fill=fill


symsel= 3
symcolor= abs(pointselection)


;  open (or filled) square
; --------------------------
if pointselection eq 1 then begin
        symsize= 1.0
	if keyword_set(fill) then begin
           usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
	endif else begin
           usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
	endelse
        symsel= 8
        symcolor= 150
endif

;  x's
; -----
if pointselection eq 2 then begin
        symsel= 7
        symcolor= 100
endif

;  open circle
; --------------
if pointselection eq 3 then begin
        symsize= 1.0
        ;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
        symsel= 8
        symcolor= 50
endif

;  open triangle
; ----------------
if pointselection eq 4 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
        usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
        symsel= 8
        symcolor= 0
endif

;  asterisk
; ----------
if pointselection eq 5 then begin
        symsel= 2
        symcolor= 200
endif

;  filled square
; --------------------------
if pointselection eq 6 then begin
        symsize= 1.0
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 8
        symcolor= 20
endif

;  diamond
; ----------
if pointselection eq 7 then begin
        symsel= 4
        symcolor= 120
endif

;  triangle (open)
; -----------------
if pointselection eq 8 then begin
        symsel= 5
        symcolor= 180
endif

;  plus sign
; ----------
if pointselection eq 9 then begin
        symsel= 1
        symcolor= 0
endif

; black line
; ----------
if pointselection eq 11 then begin
        symsel= 3
        symcolor= 0
endif

; yellow line
; ----------
if pointselection eq 12 then begin
        symsel= 3
        symcolor= 220
endif

;  small closed circle
; -----------------------
if pointselection eq 23 then begin
        symsize= 0.01
        ;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
        symsel= 8
        symcolor= 200
endif

;  open circle
; --------------
if pointselection eq 24 then begin
        symsize= 0.01
        ;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
        symsel= 8
        symcolor= 80
endif




thispsym=symsel
thiscolor=symcolor


end




