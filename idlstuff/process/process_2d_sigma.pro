; ---------------------------------------------------------------------------
;
; given x, y and vz, plus the center and the radius with which to take
; velocities, provide with velocity dispersion.
;
; --------------------------------------------------------
function process_2d_sigma, x, y, vz, Rin, center=center


if keyword_set(center) then begin
        use_x= x-center[0]
        use_y= y-center[1]
endif else begin
        use_x= x
        use_y= y
endelse

use_rxy= sqrt(use_x*use_x + use_y*use_y)

idx= where(use_rxy LE Rin)

if idx(0) eq -1 then begin
        print, 'no particles within r=',Rin
        return, 0
endif else begin
        ;print, "n= ",n_elements(idx), " out of ", n_elements(use_rxy), " inside r=",Rin
        use_rxy= use_rxy(idx)
        use_vz= vz(idx)

        if n_elements(use_vz) eq 1 then begin
                sigma= use_vz
                return, 0
        endif

        v_moment= moment(use_vz)
        sigma= sqrt(v_moment[1])
        v_avg= v_moment[0]
        ;print, "sigma=",sigma,"   v_avg=",v_avg

        return, sigma
endelse



end




