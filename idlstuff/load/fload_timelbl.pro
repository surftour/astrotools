function fload_timelbl, h, digs_after_decimal, noteq=noteq, round_off=round_off


    COMMON GalaxyHeader

    if not keyword_set(digs_after_decimal) then digs_after_decimal= 1

    thecurrenttime = time

    ;
    ;------------
    if h gt 0 then thecurrenttime = time / h

    ;
    ;------------
    if thecurrenttime ge 10.0 then digs_after_decimal= digs_after_decimal + 1
    if thecurrenttime ge 100.0 then digs_after_decimal= digs_after_decimal + 1
    if thecurrenttime ge 1000.0 then digs_after_decimal= digs_after_decimal + 1

    ;
    ;------------
    if keyword_set(round_off) then thecurrentime= thecurrenttime + 5.0/(10^(digs_after_decimal))

    ;
    ;------------
    if keyword_set(noteq) then begin
	lblz = strcompress(string(thecurrenttime),/remove_all)
        digs= digs_after_decimal+2
        lblz = strmid(lblz,0,digs)        ; 0.x (2+digits after decimal)
    endif else begin
	lblz = '!6T='+strcompress(string(thecurrenttime),/remove_all)
	;digs= digs_after_decimal+4
	digs= digs_after_decimal+6        ; change to 6 so we can grab !6
	lblz = strmid(lblz,0,digs)        ; T=0.x (4+digits after decimal)
    endelse

    return, lblz

end





