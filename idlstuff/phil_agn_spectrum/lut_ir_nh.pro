;; uses lookup tables for IR atten + re-rad calculation to speed it up 
;;   in fitting the bolometric QLFs
;;
;; lookup for IR band luminosities L_band(L,NH)
;;   input L is a vector, NH is a scalar
;;
function lut_ir_NH, log_L, log_NH
	COMMON IR_TABLES_DAT
	if (log_NH LT IR_TBL_NH_MIN) then begin
		n = (log_NH - IR_TBL_NH_MIN)/IR_TBL_DNH
		f1 = IR_TBL_LLb[*,0]
		f2 = IR_TBL_LLb[*,1]
		f  = alog10(f1) + alog10(f2/f1)*n
	endif
	if (log_NH GT IR_TBL_NH_MAX) then begin
		n = (log_NH - IR_TBL_NH_MAX)/IR_TBL_DNH
		f1 = IR_TBL_LLb[*,n_elements(IR_TBL_NH)-2]
		f2 = IR_TBL_LLb[*,n_elements(IR_TBL_NH)-1]
		f  = alog10(f1) + alog10(f2/f1)*n
	endif 
	if ((log_NH GE IR_TBL_NH_MIN) AND (log_NH LE IR_TBL_NH_MAX)) then begin
		n = (log_NH - IR_TBL_NH_MIN)/IR_TBL_DNH
		f1 = IR_TBL_LLb[*,n]
		f2 = IR_TBL_LLb[*,n+1.]
		f  = alog10(f1) + alog10(f2/f1)*(n-fix(n))
	endif
	fb = INTERPOL(f,IR_TBL_L,log_L)
	l_band = DOUBLE(log_L - fb)
	return, l_band
end
