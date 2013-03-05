;; uses lookup tables for IR atten + re-rad calculation to speed it up 
;;   in fitting the bolometric QLFs
;;
;; lookup for IR band luminosities L_band(L,NH)
;;   input L is a scalar, NH is a vector
;;
function lut_ir_lum_L, log_L, log_NH
	COMMON IR_TABLES_DAT
	if (log_L LT IR_TBL_L_MIN) then begin
		n = (log_L - IR_TBL_L_MIN)/IR_TBL_DL
		f1 = IR_TBL_LLb[0,*]
		f2 = IR_TBL_LLb[1,*]
		f  = alog10(f1) + alog10(f2/f1)*n
	endif
	if (log_L GT IR_TBL_L_MAX) then begin
		n = (log_L - IR_TBL_L_MAX)/IR_TBL_DL
		f1 = IR_TBL_LLb[n_elements(IR_TBL_L)-2,*]
		f2 = IR_TBL_LLb[n_elements(IR_TBL_L)-1,*]
		f  = alog10(f1) + alog10(f2/f1)*n
	endif 
	if ((log_L GE IR_TBL_L_MIN) AND (log_L LE IR_TBL_L_MAX)) then begin
		n = (log_L - IR_TBL_L_MIN)/IR_TBL_DL
		f1 = IR_TBL_LLb[n,*]
		f2 = IR_TBL_LLb[n+1.,*]
		f  = alog10(f1) + alog10(f2/f1)*(n-fix(n))
	endif
	fb = INTERPOL(f,IR_TBL_NH,log_NH)
	l_band = DOUBLE(log_L - fb)
	return, l_band
end
