pro methods_gas_sd, junk


   frun= "/data/tcox/Sbc201a-u4"
   sendto= 'ps'
   xlen= 80.0

   contour_gas, frun, 6, xlen, sendto, /nolabel, /pubstyle, msg='', filename='mm_gsd1.eps'
   ; contour2 actually needed to be set in contour_makeplot
   ; to make the colorbar key

   contour_gas, frun, 12, xlen, sendto, /nolabel, /pubstyle, msg=' ', filename='mm_gsd2.eps'
   contour_gas, frun, 14, xlen, sendto, /nolabel, /pubstyle, msg=' ', filename='mm_gsd3.eps'
   contour_gas, frun, 16, xlen, sendto, /nolabel, /pubstyle, msg=' ', filename='mm_gsd4.eps'
   contour_gas, frun, 18, xlen, sendto, /nolabel, /pubstyle, msg=' ', filename='mm_gsd5.eps'
   contour_gas, frun, 24, xlen, sendto, /nolabel, /pubstyle, msg=' ', filename='mm_gsd6.eps'
   contour_gas, frun, 32, xlen, sendto, /nolabel, /pubstyle, msg=' ', filename='mm_gsd7.eps'
   contour_gas, frun, 34, xlen, sendto, /nolabel, /pubstyle, msg=' ', filename='mm_gsd8.eps'
   contour_gas, frun, 36, xlen, sendto, /nolabel, /pubstyle, msg=' ', filename='mm_gsd9.eps'
   contour_gas, frun, 38, xlen, sendto, /nolabel, /pubstyle, msg=' ', filename='mm_gsd10.eps'
   contour_gas, frun, 42, xlen, sendto, /nolabel, /pubstyle, msg=' ', filename='mm_gsd11.eps'
   contour_gas, frun, 60, xlen, sendto, /nolabel, /pubstyle, msg=' ', filename='mm_gsd12.eps'




end







pro isolated_gas_sd, junk


   ;frun= "/data/tcox/Sbc11i4-u4"
   sendto= 'ps'
   xlen= 50.0
   ;contour_gas, frun, 10, xlen, sendto, /nolabel, /pubstyle, msg='n2low', filename='iso_gsd_fb1.eps'

   ;frun= "execute/Sbc11i4-u67"
   ;sendto= 'ps'
   ;xlen= 50.0
   ;contour_gas, frun, 10, xlen, sendto, /nolabel, /pubstyle, msg='n2toolow', filename='iso_gsd_fb2.eps'

   contour_gas,   'execute/Sbc11i4-u57', 10, xlen, sendto, /nolabel, /pubstyle, msg='n2low', filename='iso_gsd_2l.eps'
   contour_gas, '/data/tcox/Sbc11i4-u4', 10, xlen, sendto, /nolabel, /pubstyle, msg='n2med', filename='iso_gsd_2m.eps'
   contour_gas,   'execute/Sbc11i4-u43', 10, xlen, sendto, /nolabel, /pubstyle, msg='n2high', filename='iso_gsd_2h.eps'

   contour_gas,   'execute/Sbc11i4-u56', 10, xlen, sendto, /nolabel, /pubstyle, msg='n1low', filename='iso_gsd_1l.eps'
   contour_gas,   'execute/Sbc11i4-u53', 10, xlen, sendto, /nolabel, /pubstyle, msg='n1med', filename='iso_gsd_1m.eps'
   contour_gas,   'execute/Sbc11i4-u52', 10, xlen, sendto, /nolabel, /pubstyle, msg='n1high', filename='iso_gsd_1h.eps'

   contour_gas,   'execute/Sbc11i4-u55', 10, xlen, sendto, /nolabel, /pubstyle, msg='n0low', filename='iso_gsd_0l.eps'
   contour_gas,   'execute/Sbc11i4-u50', 10, xlen, sendto, /nolabel, /pubstyle, msg='n0med', filename='iso_gsd_0m.eps'
   contour_gas,   'execute/Sbc11i4-u51', 10, xlen, sendto, /nolabel, /pubstyle, msg='n0high', filename='iso_gsd_0h.eps'


end


