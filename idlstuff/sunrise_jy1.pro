pro plot_lir_orbits_nomf
      SET_PLOT,'PS'
      DEVICE,filename='lir_orbits_nomf.eps',/color
      loadct,13

      b5k_wagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5kv2/gtd50/fpdr0.0/agn0.1_nomf_8levels/getflux.fits',1)
      b5k_nagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5kv2/gtd50/fpdr0.0/agn0.0_nomf_8levels/getflux.fits',1)

      b5h_wagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5hv2/gtd50/fpdr0.0/agn0.1_nomf_8levels/getflux.fits',1)
      b5h_nagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5hv2/gtd50/fpdr0.0/agn0.0_nomf_8levels/getflux.fits',1)

      b5e_wagn = mrdfits('/n/data/hernquist_lab/jyounger/sunruns/z3/b5b5ev2/gtd50/fpdr0.0/agn0.1_nomf_8levels/getflux.fits',1)
      b5e_nagn = mrdfits('/n/data/hernquist_lab/jyounger/sunruns/z3/b5b5ev2/gtd50/fpdr0.0/agn0.0_nomf_8levels/getflux.fits',1)

      tref = b5k_wagn.time[where(b5k_wagn.lir_8_1000[*,0] EQ max(b5k_wagn.lir_8_1000[*,0]))]
      print,tref[0]
      b5k_wagn.time = b5k_wagn.time-tref[0]
      b5k_nagn.time = b5k_nagn.time-tref[0]

      tref = b5h_wagn.time[where(b5h_wagn.lir_8_1000[*,0] EQ max(b5h_wagn.lir_8_1000[*,0]))]
      print,tref[0]
      b5h_wagn.time = b5h_wagn.time-tref[0]
      b5h_nagn.time = b5h_nagn.time-tref[0]

      tref = b5e_wagn.time[where(b5e_wagn.lir_8_1000[*,0] EQ max(b5e_wagn.lir_8_1000[*,0]))]
      print,tref[0]
      b5e_wagn.time = b5e_wagn.time-tref[0]
      b5e_nagn.time = b5e_nagn.time-tref[0]

      plot,[-1],[-1],xrange=[-100.,100.],xstyle=1,yrange=[0.,10.],ystyle=1,$
        xtitle='Time Relative to Peak [Myr]',ytitle=TeXtoIDL('L_{IR}(8-1000\mum) [10^{12} L_{sun}]'),$
        charthick=4,charsize=1.5

      oplot,b5k_wagn.time,b5k_wagn.lir_8_1000[*,0]/1.D+12,THICK=4,COLOR=0,LINESTYLE=0
      oplot,b5k_nagn.time,b5k_nagn.lir_8_1000[*,0]/1.D+12,THICK=4,COLOR=0,LINESTYLE=2

      oplot,b5h_wagn.time,b5h_wagn.lir_8_1000[*,0]/1.D+12,THICK=4,COLOR=80,LINESTYLE=0
      oplot,b5h_nagn.time,b5h_nagn.lir_8_1000[*,0]/1.D+12,THICK=4,COLOR=80,LINESTYLE=2

      oplot,b5e_wagn.time,b5e_wagn.lir_8_1000[*,0]/1.D+12,THICK=4,COLOR=255,LINESTYLE=0
      oplot,b5e_nagn.time,b5e_nagn.lir_8_1000[*,0]/1.D+12,THICK=4,COLOR=255,LINESTYLE=2
      
      LEGEND,['b5e','b5k','b5h','w/ AGN','w/o AGN'],THICK=4,COLOR=[255,80,0,0,0],$
        LINESTYLE=0,/LEFT,BOX=0,charsize=1.3,charthick=4

      xyouts,0.4,0.2,'SF ISM; fpdr=0; z=2.5',/normal,charsize=1.3,charthick=4

      DEVICE,/CLOSE

END

pro plot_smm_orbits_nomf
      SET_PLOT,'PS'
      DEVICE,filename='smm_orbits_nomf.eps',/color
      loadct,13

      b5k_wagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5kv2/gtd50/fpdr0.0/agn0.1_nomf_8levels/getflux.fits',1)
      b5k_nagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5kv2/gtd50/fpdr0.0/agn0.0_nomf_8levels/getflux.fits',1)

      b5h_wagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5hv2/gtd50/fpdr0.0/agn0.1_nomf_8levels/getflux.fits',1)
      b5h_nagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5hv2/gtd50/fpdr0.0/agn0.0_nomf_8levels/getflux.fits',1)

      b5e_wagn = mrdfits('/n/data/hernquist_lab/jyounger/sunruns/z3/b5b5ev2/gtd50/fpdr0.0/agn0.1_nomf_8levels/getflux.fits',1)
      b5e_nagn = mrdfits('/n/data/hernquist_lab/jyounger/sunruns/z3/b5b5ev2/gtd50/fpdr0.0/agn0.0_nomf_8levels/getflux.fits',1)

      tref = b5k_wagn.time[where(b5k_wagn.s_850um[*,0] EQ max(b5k_wagn.s_850um[*,0]))]
      print,tref[0]
      b5k_wagn.time = b5k_wagn.time-tref[0]
      b5k_nagn.time = b5k_nagn.time-tref[0]

      tref = b5h_wagn.time[where(b5h_wagn.s_850um[*,0] EQ max(b5h_wagn.s_850um[*,0]))]
      print,tref[0]
      b5h_wagn.time = b5h_wagn.time-tref[0]
      b5h_nagn.time = b5h_nagn.time-tref[0]

      tref = b5e_wagn.time[where(b5e_wagn.s_850um[*,0] EQ max(b5e_wagn.s_850um[*,0]))]
      print,tref[0]
      b5e_wagn.time = b5e_wagn.time-tref[0]
      b5e_nagn.time = b5e_nagn.time-tref[0]

      plot,[-1],[-1],xrange=[-75.,75.],xstyle=1,yrange=[2.,8.],ystyle=1,$
        xtitle='Time Relative to Peak [Myr]',ytitle=TeXtoIDL('S_{850\mum} [mJy]'),$
        charthick=4,charsize=1.5

      oplot,b5k_wagn.time,b5k_wagn.S_850um[*,0],THICK=4,COLOR=0,LINESTYLE=0
      oplot,b5k_nagn.time,b5k_nagn.S_850um[*,0],THICK=4,COLOR=0,LINESTYLE=2

      oplot,b5h_wagn.time,b5h_wagn.S_850um[*,0],THICK=4,COLOR=80,LINESTYLE=0
      oplot,b5h_nagn.time,b5h_nagn.S_850um[*,0],THICK=4,COLOR=80,LINESTYLE=2

      oplot,b5e_wagn.time,b5e_wagn.S_850um[*,0],THICK=4,COLOR=255,LINESTYLE=0
      oplot,b5e_nagn.time,b5e_nagn.S_850um[*,0],THICK=4,COLOR=255,LINESTYLE=2
      
      LEGEND,['b5e','b5k','b5h','w/ AGN','w/o AGN'],THICK=4,COLOR=[255,80,0,0,0],$
        LINESTYLE=0,/LEFT,BOX=0,charsize=1.3,charthick=4

      xyouts,0.4,0.2,'SF ISM; fpdr=0; z=2.5',/normal,charsize=1.3,charthick=4

      DEVICE,/CLOSE

END

pro plot_lir_orbits
      SET_PLOT,'PS'
      DEVICE,filename='lir_orbits.eps',/color
      loadct,13

      b5k_wagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5kv2/gtd50/fpdr1.0/agn0.1_8levels/getflux.fits',1)
      b5k_nagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5kv2/gtd50/fpdr1.0/agn0.0_8levels/getflux.fits',1)

      b5h_wagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5hv2/gtd50/fpdr1.0/agn0.1_8levels/getflux.fits',1)
      b5h_nagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5hv2/gtd50/fpdr1.0/agn0.0_8levels/getflux.fits',1)

      b5e_wagn = mrdfits('/n/data/hernquist_lab/jyounger/sunruns/z3/b5b5ev2/gtd50/fpdr1.0/agn0.1_8levels/getflux.fits',1)
      b5e_nagn = mrdfits('/n/data/hernquist_lab/jyounger/sunruns/z3/b5b5ev2/gtd50/fpdr1.0/agn0.0_8levels/getflux.fits',1)

      tref = b5k_wagn.time[where(b5k_wagn.lir_8_1000[*,0] EQ max(b5k_wagn.lir_8_1000[*,0]))]
      print,tref[0]
      b5k_wagn.time = b5k_wagn.time-tref[0]
      b5k_nagn.time = b5k_nagn.time-tref[0]

      tref = b5h_wagn.time[where(b5h_wagn.lir_8_1000[*,0] EQ max(b5h_wagn.lir_8_1000[*,0]))]
      print,tref[0]
      b5h_wagn.time = b5h_wagn.time-tref[0]
      b5h_nagn.time = b5h_nagn.time-tref[0]

      tref = b5e_wagn.time[where(b5e_wagn.lir_8_1000[*,0] EQ max(b5e_wagn.lir_8_1000[*,0]))]
      print,tref[0]
      b5e_wagn.time = b5e_wagn.time-tref[0]
      b5e_nagn.time = b5e_nagn.time-tref[0]

      plot,[-1],[-1],xrange=[-100.,100.],xstyle=1,yrange=[0.,13.],ystyle=1,$
        xtitle='Time Relative to Peak [Myr]',ytitle=TeXtoIDL('L_{IR}(8-1000\mum) [10^{12} L_{sun}]'),$
        charthick=4,charsize=1.5

      oplot,b5k_wagn.time,b5k_wagn.lir_8_1000[*,0]/1.D+12,THICK=4,COLOR=0,LINESTYLE=0
      oplot,b5k_nagn.time,b5k_nagn.lir_8_1000[*,0]/1.D+12,THICK=4,COLOR=0,LINESTYLE=2

      oplot,b5h_wagn.time,b5h_wagn.lir_8_1000[*,0]/1.D+12,THICK=4,COLOR=80,LINESTYLE=0
      oplot,b5h_nagn.time,b5h_nagn.lir_8_1000[*,0]/1.D+12,THICK=4,COLOR=80,LINESTYLE=2

      oplot,b5e_wagn.time,b5e_wagn.lir_8_1000[*,0]/1.D+12,THICK=4,COLOR=255,LINESTYLE=0
      oplot,b5e_nagn.time,b5e_nagn.lir_8_1000[*,0]/1.D+12,THICK=4,COLOR=255,LINESTYLE=2
      
      LEGEND,['b5e','b5k','b5h','w/ AGN','w/o AGN'],THICK=4,COLOR=[255,80,0,0,0],$
        LINESTYLE=0,/LEFT,BOX=0,charsize=1.3,charthick=4

      xyouts,0.4,0.2,'MF ISM; fpdr=1; z=2.5',/normal,charsize=1.3,charthick=4

      DEVICE,/CLOSE

END

pro plot_smm_orbits
      SET_PLOT,'PS'
      DEVICE,filename='smm_orbits.eps',/color
      loadct,13

      b5k_wagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5kv2/gtd50/fpdr1.0/agn0.1_8levels/getflux.fits',1)
      b5k_nagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5kv2/gtd50/fpdr1.0/agn0.0_8levels/getflux.fits',1)

      b5h_wagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5hv2/gtd50/fpdr1.0/agn0.1_8levels/getflux.fits',1)
      b5h_nagn=mrdfits('/n/scratch/hernquist_lab/jyounger/sunruns/b5hv2/gtd50/fpdr1.0/agn0.0_8levels/getflux.fits',1)

      b5e_wagn = mrdfits('/n/data/hernquist_lab/jyounger/sunruns/z3/b5b5ev2/gtd50/fpdr1.0/agn0.1_8levels/getflux.fits',1)
      b5e_nagn = mrdfits('/n/data/hernquist_lab/jyounger/sunruns/z3/b5b5ev2/gtd50/fpdr1.0/agn0.0_8levels/getflux.fits',1)

      tref = b5k_wagn.time[where(b5k_wagn.s_850um[*,0] EQ max(b5k_wagn.s_850um[*,0]))]
      print,tref[0]
      b5k_wagn.time = b5k_wagn.time-tref[0]
      b5k_nagn.time = b5k_nagn.time-tref[0]

      tref = b5h_wagn.time[where(b5h_wagn.s_850um[*,0] EQ max(b5h_wagn.s_850um[*,0]))]
      print,tref[0]
      b5h_wagn.time = b5h_wagn.time-tref[0]
      b5h_nagn.time = b5h_nagn.time-tref[0]

      tref = b5e_wagn.time[where(b5e_wagn.s_850um[*,0] EQ max(b5e_wagn.s_850um[*,0]))]
      print,tref[0]
      b5e_wagn.time = b5e_wagn.time-tref[0]
      b5e_nagn.time = b5e_nagn.time-tref[0]

      plot,[-1],[-1],xrange=[-75.,75.],xstyle=1,yrange=[0.,8.],ystyle=1,$
        xtitle='Time Relative to Peak [Myr]',ytitle=TeXtoIDL('S_{850\mum} [mJy]'),$
        charthick=4,charsize=1.5

      oplot,b5k_wagn.time,b5k_wagn.S_850um[*,0],THICK=4,COLOR=0,LINESTYLE=0
      oplot,b5k_nagn.time,b5k_nagn.S_850um[*,0],THICK=4,COLOR=0,LINESTYLE=2

      oplot,b5h_wagn.time,b5h_wagn.S_850um[*,0],THICK=4,COLOR=80,LINESTYLE=0
      oplot,b5h_nagn.time,b5h_nagn.S_850um[*,0],THICK=4,COLOR=80,LINESTYLE=2

      oplot,b5e_wagn.time,b5e_wagn.S_850um[*,0],THICK=4,COLOR=255,LINESTYLE=0
      oplot,b5e_nagn.time,b5e_nagn.S_850um[*,0],THICK=4,COLOR=255,LINESTYLE=2
      
      LEGEND,['b5e','b5k','b5h','w/ AGN','w/o AGN'],THICK=4,COLOR=[255,80,0,0,0],$
        LINESTYLE=0,/LEFT,BOX=0,charsize=1.3,charthick=4

      xyouts,0.4,0.2,'MF ISM; fpdr=1; z=2.5',/normal,charsize=1.3,charthick=4

      DEVICE,/CLOSE

END

