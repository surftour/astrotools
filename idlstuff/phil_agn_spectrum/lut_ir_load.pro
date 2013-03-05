;; loads lookup tables for IR atten + re-rad calculation to speed it up 
;;   in fitting the bolometric QLFs
;;
COMMON IR_TABLES_DAT,	IR_TBL_L, IR_TBL_NH, IR_TBL_LLb, IR_TBL_L_MAX, $
						IR_TBL_DL, IR_TBL_L_MIN, IR_TBL_NH_MIN, IR_TBL_NH_MAX, $
						IR_TBL_DNH

function lut_ir_load
	COMMON IR_TABLES_DAT

	homedir = return_idl_routines_homedir(0)+'/agn_spectrum/'
	OPENR,luntemp,homedir+'ir.15m.lum.dat',/get_lun
	READF,luntemp,IR_TBL_N_L
	IR_TBL_L = fltarr(IR_TBL_N_L)
	READF,luntemp,IR_TBL_L
	READF,luntemp,IR_TBL_N_NH
	IR_TBL_NH = fltarr(IR_TBL_N_NH)
	READF,luntemp,IR_TBL_NH
	IR_TBL_LLb = fltarr(IR_TBL_N_L,IR_TBL_N_NH)
	READF,luntemp,IR_TBL_LLb
	CLOSE,luntemp
	FREE_LUN,luntemp

	IR_TBL_L_MIN = IR_TBL_L[0]	
	IR_TBL_DL    = IR_TBL_L[1] - IR_TBL_L[0]
	IR_TBL_L_MAX = IR_TBL_L[IR_TBL_N_L-2]
	IR_TBL_NH_MIN= IR_TBL_NH[0]
	IR_TBL_NH_MAX= IR_TBL_NH[IR_TBL_N_NH-2]
	IR_TBL_DNH   = IR_TBL_NH[1] - IR_TBL_NH[0]
end

