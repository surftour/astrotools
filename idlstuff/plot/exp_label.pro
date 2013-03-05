FUNCTION exp_label, axis, index, number

  IF number EQ 0 THEN RETURN, '0' ;; Special case

  ex = String(number, '(e7.0)')
  pt = StrPos(ex, '.')

  first = StrMid(ex, pt-1, 1)
  sign = StrMid(ex, pt+2, 1)
  exponent = StrMid(ex, pt+3, 100)

  ;; Shave off leading zero in exponent

  WHILE (StrMid(exponent, 0, 1) EQ '0' and StrLen(exponent) gt 1) DO exponent = StrMid(exponent, 1, 100)

  IF first eq '1' THEN BEGIN
	IF sign EQ '-' THEN RETURN, '10!E' + sign + exponent + '!N' $
	ELSE                RETURN, '10!E' + exponent + '!N'
  ENDIF

  IF sign EQ '-' THEN   RETURN, first + 'x10!A' + sign + exponent $
     ELSE               RETURN, first + 'x10!A' + exponent
END
