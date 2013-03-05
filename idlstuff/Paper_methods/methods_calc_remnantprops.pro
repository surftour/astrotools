pro methods_calc_remnantprops, junk

   ; need to do this by hand - .run prodir

   ; n=2
   ; ----------
   frun="/data/tcox/Sbc201a-u4"
   sendto='ps'
   snapnum=60

   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=10.0


   frun="/data/tcox/Sbc201a-u43"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   frun="execute/Sbc201a-u57"
   ok=fload_snapshot(frun,snapnum)
   plot_Re_sigma_random, frun, sendto
   plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0


   ; n=1
   ; ----------
   ;frun="/data/tcox/Sbc201a-u40"
   frun="execute/Sbc201a-u53"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   ;frun="execute/Sbc201a-u46"
   frun="execute/Sbc201a-u52"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   frun="execute/Sbc201a-u56"
   ok=fload_snapshot(frun,snapnum)
   plot_Re_sigma_random, frun, sendto
   plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0


   ; n=0
   ; ----------
   ;frun="/data/tcox/Sbc201a-u8"
   frun="execute/Sbc201a-u50"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   ;frun="execute/Sbc201a-u45"
   frun="execute/Sbc201a-u51"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   frun="execute/Sbc201a-u55"
   ok=fload_snapshot(frun,snapnum)
   plot_Re_sigma_random, frun, sendto
   plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0


   ; collisionless
   ; ----------

   frun="/data/tcox/Sbc201a1-u4"
   ;ok=fload_snapshot(frun,61)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   frun="/data/tcox/Sbc201a2-u4"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=3.6


   ;--------------------------------
   ; resolution runs

   ; warning - this is for a strange
   ;    time, need to test against
   ;    something comparable
   frun="/data/tcox/Sbc201a10x-n4"
   ;ok=fload_snapshot(frun,99)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   frun="/data/tcox/Sbc201a4x-n4"
   ;ok=fload_snapshot(frun,600)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   frun="/data/tcox/Sbc201a2x-n4"
   ;ok=fload_snapshot(frun,600)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   frun="/data/tcox/Sbc201a-n4"
   ;ok=fload_snapshot(frun,600)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0




   ;--------------------------------
   ; c_star runs

   frun="execute/Sbc201a-u54"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   frun="execute/Sbc201a-u14"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   frun="execute/Sbc201a-u44"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   ;frun="/data/tcox/Sbc201a-u4"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto

   frun="execute/Sbc201a-u64"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0
   ;plot_Re_sigma_random, frun, sendto

   frun="execute/Sbc201a-u74"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0

   ;--------------------------------

   frun="/data/tcox/Sbc201a-u4a"
   ;ok=fload_snapshot(frun,snapnum)
   ;plot_Re_sigma_random, frun, sendto
   ;plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0




end










;===============================================



pro test_this, junk

   ; need to do this by hand - .run prodir

   frun="/data/tcox/Sbc201a-u4"
   sendto='ps'
   snapnum=60

   ok=fload_snapshot(frun,snapnum)
   plot_Re_sigma_random, frun, sendto
   plot_Re_sigma_random, frun, sendto, sigma_aperture=5.0
   plot_Re_sigma_random, frun, sendto, sigma_aperture=10.0




end



