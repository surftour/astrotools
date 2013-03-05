

;-----------------------------------
;
;  time evolution (isolated)
;


pro doit_massprofiles, junk


gals= {frun:"data/z3/iso_b5", snapnum: 0, component:"all", linecolor: 20}
plot_mass_profile, gals, /x_is_devac, filename="sdinit_iso_b5.eps"
gals= replicate(gals, 9)
gals[*].snapnum= [0, 5, 10, 20, 30, 50, 75, 100, 125]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
gals[*].component= "allstars"
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_iso_b5.eps"


gals= {frun:"data/z3/iso_b5v2", snapnum: 0, component:"all", linecolor: 20}
plot_mass_profile, gals, /x_is_devac, filename="sdinit_iso_b5v2.eps"
gals= replicate(gals, 9)
gals[*].snapnum= [0, 5, 10, 20, 30, 50, 75, 100, 125]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
gals[*].component= "allstars"
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_iso_b5v2.eps"


gals= {frun:"data/z3/iso_b5v3", snapnum: 0, component:"all", linecolor: 20}
plot_mass_profile, gals, /x_is_devac, filename="sdinit_iso_b5v3.eps"
gals= replicate(gals, 9)
gals[*].snapnum= [0, 5, 10, 20, 30, 50, 75, 100, 125]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
gals[*].component= "allstars"
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_iso_b5v3.eps"


gals= {frun:"data/z3/iso_b5v4", snapnum: 0, component:"all", linecolor: 20}
plot_mass_profile, gals, /x_is_devac, filename="sdinit_iso_b5v4.eps"
gals= replicate(gals, 9)
gals[*].snapnum= [0, 5, 10, 20, 30, 50, 75, 100, 125]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
gals[*].component= "allstars"
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_iso_b5v4.eps"


gals= {frun:"data/z3/iso_b5v5", snapnum: 0, component:"all", linecolor: 20}
plot_mass_profile, gals, /x_is_devac, filename="sdinit_iso_b5v5.eps"
gals= replicate(gals, 9)
gals[*].snapnum= [0, 5, 10, 20, 30, 50, 75, 100, 125]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
gals[*].component= "allstars"
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_iso_b5v5.eps"


gals= {frun:"data1/z3/iso_z5test", snapnum: 0, component:"all", linecolor: 20}
plot_mass_profile, gals, /x_is_devac, filename="sdinit_iso_z5test.eps"
gals= replicate(gals, 9)
gals[*].snapnum= [0, 5, 10, 20, 30, 50, 75, 100, 125]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
gals[*].component= "allstars"
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_iso_z5test.eps"






;-----------------------------------
;
;  by component
;


plot_mass_profile, {frun:"data/z3/b5e", snapnum:75, component:"all", linecolor:0}, /x_is_devac, filename="test.eps"
;plot_sd_profile_all, "data/z3/b5e", 55, filename="sd_b5e.eps", /x_is_devac

;plot_sd_profile_all, "data/z3/b5v2e", 55, filename="sd_b5v2e.eps", /x_is_devac
;plot_sd_profile_all, "data/z3/b5v3e", 55, filename="sd_b5v3e.eps", /x_is_devac
;plot_sd_profile_all, "data/z3/b5v4e", 55, filename="sd_b5v4e.eps", /x_is_devac
;plot_sd_profile_all, "data/z3/b5v4ea", 55, filename="sd_b5v4ea.eps", /x_is_devac

;plot_sd_profile_all, "data1/z3/b5e_z5test", 55, filename="sd_b5ez5test.eps", /x_is_devac
;plot_sd_profile_all, "data1/z3/b5e_z5test_lr", 55, filename="sd_b5ez5test_lr.eps", /x_is_devac
;plot_sd_profile_all, "data1/z3/b5e_z5test_N2", 55, filename="sd_b5ez5test_N2.eps", /x_is_devac




;-----------------------------------
;
;  do a comparison
;

gals=        {frun:"data1/z3/b5e", snapnum: 0, component:"all", linecolor: 20}
gals= [gals, {frun:"data1/z3/b5v2e",        color:30,   lthick:1.0, msg:"b5ev2e",        x0:0.22, y0:0.48}]
plot_mass_profile, gals, /x_is_devac, filename="sdinit_iso_z5test.eps"
gals= replicate(gals, 9)
gals[*].snapnum= [0, 5, 10, 20, 30, 50, 75, 100, 125]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
gals[*].component= "allstars"
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_iso_z5test.eps"






;-----------------------------------
;
;  with fits
;



;plot_sd_profile, "data/z3/b5e", 55, filename="sdfit_b5e.eps", /x_is_devac

;plot_sd_profile, "data/z3/b5v2e", 55, filename="sdfit_b5v2e.eps", /x_is_devac
;plot_sd_profile, "data/z3/b5v3e", 55, filename="sdfit_b5v3e.eps", /x_is_devac
;plot_sd_profile, "data/z3/b5v4e", 55, filename="sdfit_b5v4e.eps", /x_is_devac
;plot_sd_profile, "data/z3/b5v4ea", 55, filename="sdfit_b5v4ea.eps", /x_is_devac

;plot_sd_profile, "data1/z3/b5e_z5test", 55, filename="sdfit_b5ez5test.eps", /x_is_devac
;plot_sd_profile, "data1/z3/b5e_z5test_lr", 55, filename="sdfit_b5ez5test_lr.eps", /x_is_devac
;plot_sd_profile, "data1/z3/b5e_z5test_N2", 55, filename="sdfit_b5ez5test_N2.eps", /x_is_devac

;plot_mass_profile, {frun:"data/z3/b5e", snapnum:55, component:"allstars", linecolor:15}, /x_is_devac, /fitsersic, filename="test.eps"



;-----------------------------------
;
;  time evolution (mergers)
;

gals= {frun:"data/z3/b5e", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 9)
gals[*].snapnum= [55, 65, 75, 85, 95, 105, 115, 125, 135]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_b5e.eps"

gals= {frun:"data/z3/b5v2e", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 9)
gals[*].snapnum= [55, 65, 75, 85, 95, 105, 115, 125, 135]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_b5v2e.eps"

gals= {frun:"data/z3/b5v3e", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 9)
gals[*].snapnum= [55, 65, 75, 85, 95, 105, 115, 125, 135]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_b5v3e.eps"

gals= {frun:"data/z3/b5v4ea", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 9)
gals[*].snapnum= [55, 65, 75, 85, 95, 105, 115, 125, 135]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_b5v4ea.eps"

;gals= {frun:"data1/z3/b5e_z5test", snapnum: 55, component:"allstars", linecolor: 20}
;gals= replicate(gals, 9)
;gals[*].snapnum= [55, 65, 75, 85, 95, 105, 115, 125, 135]
;gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
;plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_b5e_z5test.eps"

gals= {frun:"data1/z3/b5e_z5test_lr", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 9)
;gals[*].snapnum= [55, 65, 75, 85, 95, 105, 115, 125, 135]
gals[*].snapnum= [55, 65, 75, 85, 95, 105, 115, 125, 129]      ; the runs not quite done yet
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_b5e_z5test_lr.eps"

gals= {frun:"data1/z3/b5e_z5test_N2", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 9)
gals[*].snapnum= [55, 65, 75, 85, 95, 105, 115, 125, 135]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_b5e_z5test_N2.eps"

gals= {frun:"data1/z3/stijn_b5e_v1", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 9)
gals[*].snapnum= [55, 65, 75, 85, 95, 105, 115, 125, 135]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_stijn_b5e_v1.eps"

gals= {frun:"data1/z3/stijn_b5e_v1b", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 9)
gals[*].snapnum= [55, 65, 75, 85, 95, 105, 115, 125, 135]
gals[*].linecolor= [20, 40, 60, 80, 100, 120, 140, 160, 180]
plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_stijn_b5e_v1b.eps"


end



;-----------------------------------
;
;   mdtests
;

pro doit_mp_mdtests, junk

;
gals= {frun:"data1/mdtests/twocomp_1", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 7)
gals[*].snapnum= [0, 1, 2, 3, 5, 10, 20]
gals[*].linecolor= [10, 50, 90, 130, 170, 210, 245]
;plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_twocomp_1.eps"
;plot_mass_profile, gals, /x_is_log, ctab= 1, filename="sdevol_twocomp_1_log.eps"

gals= {frun:"data1/mdtests/twocomp_max", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 7)
gals[*].snapnum= [0, 1, 2, 3, 5, 10, 20]
gals[*].linecolor= [10, 50, 90, 130, 170, 210, 245]
;plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_twocomp_max.eps"
;plot_mass_profile, gals, /x_is_log, ctab= 1, filename="sdevol_twocomp_max_log.eps"

;
gals= {frun:"data1/mdtests/bulgehalo2_1", snapnum: 55, component:"allstars", linecolor: 20}
;gals= replicate(gals, 7)
;gals[*].snapnum= [0, 1, 2, 3, 5, 10, 20]
;gals[*].linecolor= [10, 50, 90, 130, 170, 210, 245]
gals= replicate(gals, 6)
gals[*].snapnum= [0, 1, 2, 3, 5, 10]
gals[*].linecolor= [10, 50, 90, 130, 170, 210]
;plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_bulgehalo2_1.eps"
;plot_mass_profile, gals, /x_is_log, ctab= 1, filename="sdevol_bulgehalo2_1_log.eps"

gals= {frun:"data1/mdtests/bulgehalo2_max", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 7)
gals[*].snapnum= [0, 1, 2, 3, 5, 10, 20]
gals[*].linecolor= [10, 50, 90, 130, 170, 210, 245]
;plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_bulgehalo2_max.eps"
;plot_mass_profile, gals, /x_is_log, ctab= 1, filename="sdevol_bulgehalo2_max_log.eps"

;
gals= {frun:"data1/mdtests/diskhalo_1", snapnum: 55, component:"allstars", linecolor: 20}
;gals= replicate(gals, 7)
;gals[*].snapnum= [0, 1, 2, 3, 5, 10, 20]
;gals[*].linecolor= [10, 50, 90, 130, 170, 210, 245]
gals= replicate(gals, 6)
gals[*].snapnum= [0, 1, 2, 3, 5, 10]
gals[*].linecolor= [10, 50, 90, 130, 170, 210]
;plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_diskhalo_1.eps"
plot_mass_profile, gals, /x_is_log, ctab= 1, filename="sdevol_diskhalo_1_log.eps"

gals= {frun:"data1/mdtests/diskhalo_max", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 7)
gals[*].snapnum= [0, 1, 2, 3, 5, 10, 20]
gals[*].linecolor= [10, 50, 90, 130, 170, 210, 245]
;plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_diskhalo_max.eps"
;plot_mass_profile, gals, /x_is_log, ctab= 1, filename="sdevol_diskhalo_max_log.eps"

;
gals= {frun:"data1/mdtests/threecomp_1", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 7)
gals[*].snapnum= [0, 1, 2, 3, 5, 10, 20]
gals[*].linecolor= [10, 50, 90, 130, 170, 210, 245]
;plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_threecomp_1.eps"
;plot_mass_profile, gals, /x_is_log, ctab= 1, filename="sdevol_threecomp_1_log.eps"

gals[*].component= "bulge"
plot_mass_profile, gals, /x_is_log, ctab= 1, filename="sdevol_threecomp_1_bulgelog.eps"
gals[*].component= "disk"
plot_mass_profile, gals, /x_is_log, ctab= 1, filename="sdevol_threecomp_1_disklog.eps"

gals= {frun:"data1/mdtests/threecomp_max", snapnum: 55, component:"allstars", linecolor: 20}
gals= replicate(gals, 7)
gals[*].snapnum= [0, 1, 2, 3, 5, 10, 20]
gals[*].linecolor= [10, 50, 90, 130, 170, 210, 245]
;plot_mass_profile, gals, /x_is_devac, ctab= 1, filename="sdevol_threecomp_max.eps"
;plot_mass_profile, gals, /x_is_log, ctab= 1, filename="sdevol_threecomp_max_log.eps"

gals[*].component= "bulge"
plot_mass_profile, gals, /x_is_log, ctab= 1, filename="sdevol_threecomp_max_bulgelog.eps"
gals[*].component= "disk"
plot_mass_profile, gals, /x_is_log, ctab= 1, filename="sdevol_threecomp_max_disklog.eps"



end






pro doit_sfr, junk


;------------------------
;
;
;gals=        {RUNINFO, frun:"z3/b5e",           color:0,    lthick:3.0, msg:"b5e",           x0:0.22, y0:0.52}
;gals= [gals, {RUNINFO, frun:"z3/b5v2e",        color:30,   lthick:1.0, msg:"b5ev2e",        x0:0.22, y0:0.48}]
;gals= [gals, {RUNINFO, frun:"z3/b5v3e",        color:60,   lthick:1.0, msg:"b5ev3e",        x0:0.22, y0:0.44}]
;gals= [gals, {RUNINFO, frun:"z3/b5v4e",        color:90,   lthick:1.0, msg:"b5ev4e",        x0:0.22, y0:0.40}]
;gals= [gals, {RUNINFO, frun:"z3/b5v4ea",       color:120,  lthick:1.0, msg:"b5ev4ea",       x0:0.22, y0:0.36}]
;gals= [gals, {RUNINFO, frun:"z3/b5e_z5test",    color:150,  lthick:1.0, msg:"b5e_z5test",    x0:0.22, y0:0.32}]
;gals= [gals, {RUNINFO, frun:"z3/b5e_z5test_lr", color:190,  lthick:1.0, msg:"b5e_z5test_lr", x0:0.22, y0:0.28}]
;gals= [gals, {RUNINFO, frun:"z3/b5e_z5test_N2", color:220,  lthick:1.0, msg:"b5e_z5test_N2", x0:0.22, y0:0.24}]

;plot_sfr, gals, xmax=1.5, xmin= 0.0, ymax= 4000.0, ymin= 0.1, filename="sfr_b5s.eps"




;-------------------------


;/n/scratch/hernquist_lab/chayward/gadgetruns/z3/b5e_2x_res



end




;========================================================================


pro doit_halo, junk



;.run halo_profile


;gals= {frun:"data1/z3/iso_b5", snapnum: 0, component:"halo", linecolor: 20}
;halostructure, gals, filename="halo_iso_b5.eps"




;gals=        {frun:"data1/mdtests/dmonly", snapnum: 0, component:"halo", linecolor: 150}    ; bad IC within 1 kpc
;gals= [gals, {frun:"data1/mdtests/dmonly_max", snapnum: 0, component:"halo", linecolor: 0}]    
;gals=        {frun:"data1/mdtests/dmonly_max", snapnum: 0, component:"halo", linecolor: 0}
;gals= [gals, {frun:"data1/mdtests/dmonly_2", snapnum: 0, component:"halo", linecolor: 50}]
;gals=        {frun:"data1/mdtests/junk", snapnum: 0, component:"halo", linecolor: 150} 
;gals= [gals, {frun:"data1/mdtests/junk", snapnum: 0, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/dmonly_3", snapnum: 0, component:"halo", linecolor: 100}]
;gals= [gals, {frun:"data1/mdtests/dmonly_4", snapnum: 0, component:"halo", linecolor: 200}]
;gals= [gals, {frun:"data1/mdtests/dmonly_5", snapnum: 0, component:"halo", linecolor: 200}]
;gals= [gals, {frun:"data1/mdtests/dmonly_6", snapnum: 0, component:"halo", linecolor: 100}]
;gals= [gals, {frun:"data1/mdtests/dmonly_7", snapnum: 0, component:"halo", linecolor: 050}]
;halostructure, gals, filename="halo_dmonly_init.eps"

;gals=        {frun:"data1/mdtests/dmonly_8", snapnum: 0, component:"halo", linecolor: 150} 
;gals= [gals, {frun:"data1/mdtests/dmonly_9", snapnum: 0, component:"halo", linecolor: 100}]
;halostructure, gals, filename="halo_dmonly_init_ani.eps"

;gals=        {frun:"data1/mdtests/junk", snapnum: 0, component:"halo", linecolor: 150} 
;gals= [gals, {frun:"data1/mdtests/dmonly_9", snapnum: 0, component:"halo", linecolor: 100}]
;gals= [gals, {frun:"data1/mdtests/dmonly_10", snapnum: 0, component:"halo", linecolor: 50}]
;halostructure, gals, filename="halo_dmonly_init_iso.eps"





;gals=        {frun:"data1/mdtests/dmonly_max", snapnum:  0, component:"halo", linecolor: 30}
;gals= [gals, {frun:"data1/mdtests/dmonly_max", snapnum: 10, component:"halo", linecolor: 60}]
;gals= [gals, {frun:"data1/mdtests/dmonly_max", snapnum: 20, component:"halo", linecolor: 90}]
;gals= [gals, {frun:"data1/mdtests/dmonly_max", snapnum: 30, component:"halo", linecolor: 120}]
;gals= [gals, {frun:"data1/mdtests/dmonly_max", snapnum: 40, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/dmonly_max", snapnum: 50, component:"halo", linecolor: 200}]
;halostructure, gals, ctab= 1, filename="halo_dmonly_max.eps"


;gals=        {frun:"data1/mdtests/junk", snapnum:  0, component:"halo", linecolor: 30}
;gals= [gals, {frun:"data1/mdtests/junk", snapnum: 10, component:"halo", linecolor: 60}]
;gals= [gals, {frun:"data1/mdtests/junk", snapnum: 20, component:"halo", linecolor: 90}]
;gals= [gals, {frun:"data1/mdtests/junk", snapnum: 30, component:"halo", linecolor: 120}]
;gals= [gals, {frun:"data1/mdtests/junk", snapnum: 40, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/junk", snapnum: 50, component:"halo", linecolor: 200}]
;halostructure, gals, ctab= 1, filename="halo_junk.eps"


;gals=        {frun:"data1/mdtests/dmonly_3", snapnum:  0, component:"halo", linecolor: 30}
;gals= [gals, {frun:"data1/mdtests/dmonly_3", snapnum:  2, component:"halo", linecolor: 60}]
;gals= [gals, {frun:"data1/mdtests/dmonly_3", snapnum:  4, component:"halo", linecolor: 90}]
;gals= [gals, {frun:"data1/mdtests/dmonly_3", snapnum:  6, component:"halo", linecolor: 120}]
;gals= [gals, {frun:"data1/mdtests/dmonly_3", snapnum: 10, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/dmonly_3", snapnum: 13, component:"halo", linecolor: 200}]
;halostructure, gals, ctab= 1, filename="halo_dmonly_3.eps"

;gals=        {frun:"data1/mdtests/dmonly_4", snapnum:  0, component:"halo", linecolor: 30}
;gals= [gals, {frun:"data1/mdtests/dmonly_4", snapnum:  2, component:"halo", linecolor: 60}]
;gals= [gals, {frun:"data1/mdtests/dmonly_4", snapnum:  4, component:"halo", linecolor: 90}]
;gals= [gals, {frun:"data1/mdtests/dmonly_4", snapnum:  6, component:"halo", linecolor: 120}]
;gals= [gals, {frun:"data1/mdtests/dmonly_4", snapnum: 10, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/dmonly_4", snapnum: 13, component:"halo", linecolor: 200}]
;halostructure, gals, ctab= 1, filename="halo_dmonly_4.eps"


;gals=        {frun:"data1/mdtests/dmonly_5", snapnum:  0, component:"halo", linecolor: 30}
;gals= [gals, {frun:"data1/mdtests/dmonly_5", snapnum:  2, component:"halo", linecolor: 60}]
;gals= [gals, {frun:"data1/mdtests/dmonly_5", snapnum:  4, component:"halo", linecolor: 90}]
;gals= [gals, {frun:"data1/mdtests/dmonly_5", snapnum:  6, component:"halo", linecolor: 120}]
;gals= [gals, {frun:"data1/mdtests/dmonly_5", snapnum: 10, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/dmonly_5", snapnum: 13, component:"halo", linecolor: 200}]
;halostructure, gals, ctab= 1, filename="halo_dmonly_5.eps"

;gals=        {frun:"data1/mdtests/dmonly_6", snapnum:  0, component:"halo", linecolor: 30}
;gals= [gals, {frun:"data1/mdtests/dmonly_6", snapnum:  2, component:"halo", linecolor: 60}]
;gals= [gals, {frun:"data1/mdtests/dmonly_6", snapnum:  4, component:"halo", linecolor: 90}]
;gals= [gals, {frun:"data1/mdtests/dmonly_6", snapnum:  6, component:"halo", linecolor: 120}]
;gals= [gals, {frun:"data1/mdtests/dmonly_6", snapnum: 10, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/dmonly_6", snapnum: 13, component:"halo", linecolor: 200}]
;halostructure, gals, ctab= 1, filename="halo_dmonly_6.eps"

;gals=        {frun:"data1/mdtests/dmonly_7", snapnum:  0, component:"halo", linecolor: 30}
;gals= [gals, {frun:"data1/mdtests/dmonly_7", snapnum:  2, component:"halo", linecolor: 60}]
;gals= [gals, {frun:"data1/mdtests/dmonly_7", snapnum:  4, component:"halo", linecolor: 90}]
;gals= [gals, {frun:"data1/mdtests/dmonly_7", snapnum:  6, component:"halo", linecolor: 120}]
;gals= [gals, {frun:"data1/mdtests/dmonly_7", snapnum: 10, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/dmonly_7", snapnum: 13, component:"halo", linecolor: 200}]
;halostructure, gals, ctab= 1, filename="halo_dmonly_7.eps"

;gals=        {frun:"data1/mdtests/dmonly_8", snapnum:  0, component:"halo", linecolor: 30}
;gals= [gals, {frun:"data1/mdtests/dmonly_8", snapnum:  2, component:"halo", linecolor: 60}]
;gals= [gals, {frun:"data1/mdtests/dmonly_8", snapnum:  5, component:"halo", linecolor: 90}]
;gals= [gals, {frun:"data1/mdtests/dmonly_8", snapnum: 10, component:"halo", linecolor: 120}]
;gals= [gals, {frun:"data1/mdtests/dmonly_8", snapnum: 20, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/dmonly_8", snapnum: 30, component:"halo", linecolor: 200}]
;halostructure, gals, ctab= 1, filename="halo_dmonly_8.eps"

;gals=        {frun:"data1/mdtests/dmonly_9", snapnum:  0, component:"halo", linecolor: 30}
;gals= [gals, {frun:"data1/mdtests/dmonly_9", snapnum:  2, component:"halo", linecolor: 60}]
;gals= [gals, {frun:"data1/mdtests/dmonly_9", snapnum:  4, component:"halo", linecolor: 90}]
;gals= [gals, {frun:"data1/mdtests/dmonly_9", snapnum:  6, component:"halo", linecolor: 120}]
;gals= [gals, {frun:"data1/mdtests/dmonly_9", snapnum: 10, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/dmonly_9", snapnum: 13, component:"halo", linecolor: 200}]
;halostructure, gals, ctab= 1, filename="halo_dmonly_9.eps"


gals=        {frun:"data1/mdtests/iso_1", snapnum:  0, component:"halo", linecolor: 0}
gals= [gals, {frun:"data1/mdtests/iso_1", snapnum:  2, component:"halo", linecolor: 30}]
gals= [gals, {frun:"data1/mdtests/iso_1", snapnum:  5, component:"halo", linecolor: 60}]
gals= [gals, {frun:"data1/mdtests/iso_1", snapnum: 10, component:"halo", linecolor: 90}]
gals= [gals, {frun:"data1/mdtests/iso_1", snapnum: 15, component:"halo", linecolor: 120}]
gals= [gals, {frun:"data1/mdtests/iso_1", snapnum: 20, component:"halo", linecolor: 150}]
gals= [gals, {frun:"data1/mdtests/iso_1", snapnum: 25, component:"halo", linecolor: 180}]
gals= [gals, {frun:"data1/mdtests/iso_1", snapnum: 30, component:"halo", linecolor: 210}]
halostructure, gals, ctab= 1, filename="halo_iso_1.eps"

gals=        {frun:"data1/mdtests/iso_2", snapnum:  0, component:"halo", linecolor: 0}
gals= [gals, {frun:"data1/mdtests/iso_2", snapnum:  2, component:"halo", linecolor: 30}]
gals= [gals, {frun:"data1/mdtests/iso_2", snapnum:  4, component:"halo", linecolor: 90}]
gals= [gals, {frun:"data1/mdtests/iso_2", snapnum:  6, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/iso_2", snapnum:  5, component:"halo", linecolor: 60}]
;gals= [gals, {frun:"data1/mdtests/iso_2", snapnum: 10, component:"halo", linecolor: 90}]
;gals= [gals, {frun:"data1/mdtests/iso_2", snapnum: 15, component:"halo", linecolor: 120}]
;gals= [gals, {frun:"data1/mdtests/iso_2", snapnum: 20, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/iso_2", snapnum: 25, component:"halo", linecolor: 180}]
;gals= [gals, {frun:"data1/mdtests/iso_2", snapnum: 30, component:"halo", linecolor: 210}]
halostructure, gals, ctab= 1, filename="halo_iso_2.eps"



gals=        {frun:"data1/mdtests/ani_1", snapnum:  0, component:"halo", linecolor: 0}
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum:  2, component:"halo", linecolor: 30}]
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum:  5, component:"halo", linecolor: 60}]
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum: 10, component:"halo", linecolor: 90}]
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum: 15, component:"halo", linecolor: 120}]
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum: 20, component:"halo", linecolor: 150}]
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum: 25, component:"halo", linecolor: 180}]
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum: 30, component:"halo", linecolor: 210}]
halostructure, gals, ctab= 1, filename="halo_ani_1.eps"

gals=        {frun:"data1/mdtests/ani_2", snapnum:  0, component:"halo", linecolor: 0}
gals= [gals, {frun:"data1/mdtests/ani_2", snapnum:  2, component:"halo", linecolor: 30}]
gals= [gals, {frun:"data1/mdtests/ani_2", snapnum:  4, component:"halo", linecolor: 90}]
gals= [gals, {frun:"data1/mdtests/ani_2", snapnum:  6, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/ani_2", snapnum:  5, component:"halo", linecolor: 60}]
;gals= [gals, {frun:"data1/mdtests/ani_2", snapnum: 10, component:"halo", linecolor: 90}]
;gals= [gals, {frun:"data1/mdtests/ani_2", snapnum: 15, component:"halo", linecolor: 120}]
;gals= [gals, {frun:"data1/mdtests/ani_2", snapnum: 20, component:"halo", linecolor: 150}]
;gals= [gals, {frun:"data1/mdtests/ani_2", snapnum: 25, component:"halo", linecolor: 180}]
;gals= [gals, {frun:"data1/mdtests/ani_2", snapnum: 30, component:"halo", linecolor: 210}]
halostructure, gals, ctab= 1, filename="halo_ani_2.eps"




gals=        {frun:"data1/mdtests/iso_max", snapnum:  0, component:"halo", linecolor: 0}
gals= [gals, {frun:"data1/mdtests/iso_max", snapnum:  2, component:"halo", linecolor: 30}]
gals= [gals, {frun:"data1/mdtests/iso_max", snapnum:  5, component:"halo", linecolor: 60}]
gals= [gals, {frun:"data1/mdtests/iso_max", snapnum: 10, component:"halo", linecolor: 90}]
gals= [gals, {frun:"data1/mdtests/iso_max", snapnum: 15, component:"halo", linecolor: 120}]
gals= [gals, {frun:"data1/mdtests/iso_max", snapnum: 20, component:"halo", linecolor: 150}]
gals= [gals, {frun:"data1/mdtests/iso_max", snapnum: 25, component:"halo", linecolor: 180}]
gals= [gals, {frun:"data1/mdtests/iso_max", snapnum: 30, component:"halo", linecolor: 210}]
halostructure, gals, ctab= 1, filename="halo_iso_max.eps"



end






pro doit_old, junk

;gals=        {frun:"data1/mdtests/iso",      snapnum: 0, component:"halo", linecolor: 150} 
;gals= [gals, {frun:"data1/mdtests/iso_test", snapnum: 0, component:"halo", linecolor: 100}]
;gals= [gals, {frun:"data1/mdtests/iso_H90",  snapnum: 0, component:"halo", linecolor: 50}]
;halostructure, gals, filename="halo_iso_init.eps"

;gals=        {frun:"data1/mdtests/ani_1",      snapnum: 0, component:"halo", linecolor: 150} 
;gals= [gals, {frun:"data1/mdtests/ani_1_test", snapnum: 0, component:"halo", linecolor: 100}]
;gals= [gals, {frun:"data1/mdtests/ani_H90",    snapnum: 0, component:"halo", linecolor: 50}]
;halostructure, gals, filename="halo_ani_init.eps"

gals=        {frun:"data3/iso_test2",      snapnum: 0, component:"halo", linecolor: 0} 
gals= [gals, {frun:"data1/mdtests/twocomp_1", snapnum: 0, component:"halo", linecolor: 100}]
gals= [gals, {frun:"data1/mdtests/twocomp_max",    snapnum: 0, component:"halo", linecolor: 150}]
halostructure, gals, filename="halo_iso_twocomp.eps"

end

pro doit_new, junk

gals=        {frun:"data1/mdtests/iso", snapnum: 0, component:"halo", linecolor: 140} 
gals= [gals, {frun:"data1/mdtests/iso", snapnum: 1, component:"halo", linecolor: 145}]
gals= [gals, {frun:"data1/mdtests/iso", snapnum: 2, component:"halo", linecolor: 150}]
gals= [gals, {frun:"data1/mdtests/iso", snapnum: 3, component:"halo", linecolor: 155}]
gals= [gals, {frun:"data1/mdtests/iso", snapnum: 5, component:"halo", linecolor: 160}]
gals= [gals, {frun:"data1/mdtests/iso", snapnum: 10, component:"halo", linecolor: 165}]
gals= [gals, {frun:"data1/mdtests/iso", snapnum: 20, component:"halo", linecolor: 170}]
gals= [gals, {frun:"data1/mdtests/iso_test", snapnum: 0, component:"halo", linecolor: 200}]
gals= [gals, {frun:"data1/mdtests/iso_test", snapnum: 1, component:"halo", linecolor: 205}]
gals= [gals, {frun:"data1/mdtests/iso_test", snapnum: 2, component:"halo", linecolor: 210}]
gals= [gals, {frun:"data1/mdtests/iso_test", snapnum: 3, component:"halo", linecolor: 215}]
gals= [gals, {frun:"data1/mdtests/iso_test", snapnum: 5, component:"halo", linecolor: 220}]
gals= [gals, {frun:"data1/mdtests/iso_test", snapnum: 10, component:"halo", linecolor: 225}]
gals= [gals, {frun:"data1/mdtests/iso_test", snapnum: 20, component:"halo", linecolor: 230}]
gals= [gals, {frun:"data1/mdtests/iso_H90", snapnum: 0, component:"halo", linecolor: 40}]
gals= [gals, {frun:"data1/mdtests/iso_H90", snapnum: 1, component:"halo", linecolor: 45}]
gals= [gals, {frun:"data1/mdtests/iso_H90", snapnum: 2, component:"halo", linecolor: 50}]
gals= [gals, {frun:"data1/mdtests/iso_H90", snapnum: 3, component:"halo", linecolor: 55}]
gals= [gals, {frun:"data1/mdtests/iso_H90", snapnum: 5, component:"halo", linecolor: 60}]
gals= [gals, {frun:"data1/mdtests/iso_H90", snapnum: 10, component:"halo", linecolor: 65}]
gals= [gals, {frun:"data1/mdtests/iso_H90", snapnum: 20, component:"halo", linecolor: 70}]
halostructure, gals, filename="halo_iso.eps"


gals=        {frun:"data1/mdtests/ani_1", snapnum: 0, component:"halo", linecolor: 140} 
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum: 1, component:"halo", linecolor: 145}]
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum: 2, component:"halo", linecolor: 150}]
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum: 3, component:"halo", linecolor: 155}]
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum: 5, component:"halo", linecolor: 160}]
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum: 10, component:"halo", linecolor: 165}]
gals= [gals, {frun:"data1/mdtests/ani_1", snapnum: 20, component:"halo", linecolor: 170}]
gals= [gals, {frun:"data1/mdtests/ani_1_test", snapnum: 0, component:"halo", linecolor: 200}]
gals= [gals, {frun:"data1/mdtests/ani_1_test", snapnum: 1, component:"halo", linecolor: 205}]
gals= [gals, {frun:"data1/mdtests/ani_1_test", snapnum: 2, component:"halo", linecolor: 210}]
gals= [gals, {frun:"data1/mdtests/ani_1_test", snapnum: 3, component:"halo", linecolor: 215}]
gals= [gals, {frun:"data1/mdtests/ani_1_test", snapnum: 5, component:"halo", linecolor: 220}]
gals= [gals, {frun:"data1/mdtests/ani_1_test", snapnum: 10, component:"halo", linecolor: 225}]
gals= [gals, {frun:"data1/mdtests/ani_1_test", snapnum: 20, component:"halo", linecolor: 230}]
gals= [gals, {frun:"data1/mdtests/ani_1_H90", snapnum: 0, component:"halo", linecolor: 40}]
gals= [gals, {frun:"data1/mdtests/ani_1_H90", snapnum: 1, component:"halo", linecolor: 45}]
gals= [gals, {frun:"data1/mdtests/ani_1_H90", snapnum: 2, component:"halo", linecolor: 50}]
gals= [gals, {frun:"data1/mdtests/ani_1_H90", snapnum: 3, component:"halo", linecolor: 55}]
gals= [gals, {frun:"data1/mdtests/ani_1_H90", snapnum: 5, component:"halo", linecolor: 60}]
gals= [gals, {frun:"data1/mdtests/ani_1_H90", snapnum: 10, component:"halo", linecolor: 65}]
gals= [gals, {frun:"data1/mdtests/ani_1_H90", snapnum: 20, component:"halo", linecolor: 70}]
halostructure, gals, filename="halo_ani.eps"

end







pro one_run_time, frun, filename=filename

if not keyword_set(filename) then filename= "halo_time.eps"

gals=        {frun:frun, snapnum: 0, component:"halo", linecolor: 10}
gals= [gals, {frun:frun, snapnum: 1, component:"halo", linecolor: 50}]
gals= [gals, {frun:frun, snapnum: 2, component:"halo", linecolor: 90}]
gals= [gals, {frun:frun, snapnum: 3, component:"halo", linecolor: 130}]
;gals= [gals, {frun:frun, snapnum: 4, component:"halo", linecolor: 170}]
gals= [gals, {frun:frun, snapnum: 5, component:"halo", linecolor: 170}]
gals= [gals, {frun:frun, snapnum: 10, component:"halo", linecolor: 210}]
;gals= [gals, {frun:frun, snapnum: 20, component:"halo", linecolor: 245}]
halostructure, gals, ctab= 1, filename=filename

end

