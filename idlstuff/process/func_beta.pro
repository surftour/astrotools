function func_beta, x, p

	;
	; beta profile according to o'sullivan, ponman & collins '03
	;
	; p[0] = rho_0
	; p[1] = r_core
	; p[2] = 3.0 * beta - 0.5
	return, p[0]*((1.0 + ((x/p[1])^(2.0)))^(-p[2]))

	;
	; p[0] = rho_0
	; p[1] = r_core
	; p[2] = beta
	;return, p[0]*((1.0 + ((x/p[1])^(2.0)))^(-3.0*p[2]+0.5))

end


