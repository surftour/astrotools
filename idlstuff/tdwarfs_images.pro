;
; Old UPSAND procedure
;
;-------------------------
pro remnant_dwarfs_upsand, junk

	frun= "/data/tcox/Sbc201a-u4"
	xlen= 5.0
	snapnum= 60
	sendto= 'ps'
	save_center= [50.0, -55.0, -20.0]

	center_d1= save_center
	contour_gas, frun, snapnum, xlen, sendto, msg='gas', /pubstyle, filename='dwarf_g.eps', center= center_d1, /particlesonly

	center_d1= save_center
	contour_newstars, frun, snapnum, xlen, sendto, msg='new stars', /pubstyle, filename='dwarf_ns.eps', center= center_d1

	center_d1= save_center
	contour_oldstars, frun, snapnum, xlen, sendto, msg='old stars', /pubstyle, filename='dwarf_os.eps', center= center_d1

end


;
;
;   Global Images !!!Global Images !!!
;
; palantiri
; -------------------
pro make_tdwarf_images, junk

;contour_gas, "Sbc", 60, 25.0, 'ps', filename='Sbc10x25.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc", 60, 50.0, 'ps', filename='Sbc10x50.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc", 60, 100.0, 'ps', filename='Sbc10x100.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc", 60, 200.0, 'ps', filename='Sbc10x200.eps', /nolabels, /old_tj_snap, /pubstyle

;contour_gas, "Sbc10x", 60, 25.0, 'ps', filename='Sbc10x25.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc10x", 60, 50.0, 'ps', filename='Sbc10x50.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc10x", 60, 100.0, 'ps', filename='Sbc10x100.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc10x", 60, 200.0, 'ps', filename='Sbc10x200.eps', /nolabels, /old_tj_snap, /pubstyle

contour_allstars, "Sbc10x", 60, 25.0, 'ps', filename='Sbc10x25as.eps', /nolabels, /old_tj_snap, /pubstyle
contour_allstars, "Sbc10x", 60, 50.0, 'ps', filename='Sbc10x50as.eps', /nolabels, /old_tj_snap, /pubstyle
contour_allstars, "Sbc10x", 60, 100.0, 'ps', filename='Sbc10x100as.eps', /nolabels, /old_tj_snap, /pubstyle
contour_allstars, "Sbc10x", 60, 200.0, 'ps', filename='Sbc10x200as.eps', /nolabels, /old_tj_snap, /pubstyle

;contour_newstars, "Sbc10x", 60, 25.0, 'ps', filename='Sbc10x25ns.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_newstars, "Sbc10x", 60, 50.0, 'ps', filename='Sbc10x50ns.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_newstars, "Sbc10x", 60, 100.0, 'ps', filename='Sbc10x100ns.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_newstars, "Sbc10x", 60, 200.0, 'ps', filename='Sbc10x200ns.eps', /nolabels, /old_tj_snap, /pubstyle


end




;
; Old UPSAND procedure
;
;-------------------------
pro tdwarf_closeup, junk

	; 1x
	; ---
        ;frun= "Sbc"
        ;xlen= 5.0
        ;snapnum= 60
        ;save_center= [50.0, -56.0, -20.0]

	; 10x
	; ---
        frun= "Sbc10x"
        xlen= 10.0
        snapnum= 60
        ;save_center= [-22.0, -31.0, 0.0]
        ;save_center= [-32.0, 35.0, 0.0]
	;save_center= [55.7, -34.4, -3.7]   ; grp 13
	;save_center= [-29.9, -341.7, 55.4]   ; grp 6
	;save_center= [-32.5, -10.8, 3.4]   ; grp 5 - somethings wrong
	save_center= [-22.6, -30.4, -3.0]   ; grp 5 - corrected, by taking median, rather than mean
	;save_center= [58.9, -191.3, -4.4]   ; grp 3
	

	; actual particles
	; note that the point type is manually selected in contour_makeplot
        center_d1= save_center
        contour_gas, frun, snapnum, xlen, 'ps', /nolabels, /pubstyle, /old_tj_snap, filename='dwarf_g.eps', center= center_d1, /particlesonly

	; smoothed gas image
        ;center_d1= save_center
        ;contour_gas, frun, snapnum, xlen, 'ps', /nolabels, /pubstyle, /old_tj_snap, filename='dwarf_g_sm.eps', center= center_d1

	;
	; see particle info above
	;
        ;center_d1= save_center
        ;contour_newstars, frun, snapnum, xlen, 'ps', /nolabels, /pubstyle, /old_tj_snap, filename='dwarf_ns.eps', center= center_d1, /particlesonly

        ;center_d1= save_center
        ;contour_oldstars, frun, snapnum, xlen, 'ps', /nolabels, /pubstyle, /old_tj_snap, filename='dwarf_os.eps', center= center_d1, /particlesonly

        ;center_d1= save_center
        ;contour_dm, frun, snapnum, xlen, 'ps', /nolabels, /pubstyle, /old_tj_snap, filename='dwarf_dm.eps', center= center_d1, /particlesonly

end






pro tds_where_from, junk

;
;
; needs a fair amount of personal attention,
; although, all of the fiddling is done in
; one file; contour_makeplot
;
;  * change background to b/w
;  * then needed to turn on the reading of the idlists,
;    as well as adjust the files it grabs by hand
;

;contour_gas, "Sbc10x", 5, 100.0, 'ps', filename='Sbc10x50_early.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc10x", 8, 100.0, 'ps', filename='Sbc10x50_8.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc10x", 10, 100.0, 'ps', filename='Sbc10x50_10.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc10x", 13, 100.0, 'ps', filename='Sbc10x50_13.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc10x", 16, 100.0, 'ps', filename='Sbc10x50_16.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc10x", 18, 100.0, 'ps', filename='Sbc10x50_18.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc10x", 20, 200.0, 'ps', filename='Sbc10x50_20.eps', /nolabels, /old_tj_snap, /pubstyle
contour_gas, "Sbc10x", 25, 200.0, 'ps', filename='Sbc10x50_25.eps', /nolabels, /old_tj_snap, /pubstyle
contour_gas, "Sbc10x", 30, 200.0, 'ps', filename='Sbc10x50_30.eps', /nolabels, /old_tj_snap, /pubstyle
contour_gas, "Sbc10x", 36, 200.0, 'ps', filename='Sbc10x50_36.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc10x", 40, 200.0, 'ps', filename='Sbc10x50_40.eps', /nolabels, /old_tj_snap, /pubstyle
;contour_gas, "Sbc10x", 50, 200.0, 'ps', filename='Sbc10x50_50.eps', /nolabels, /old_tj_snap, /pubstyle
contour_gas, "Sbc10x", 60, 200.0, 'ps', filename='Sbc10x50_60.eps', /nolabels, /old_tj_snap, /pubstyle



end


