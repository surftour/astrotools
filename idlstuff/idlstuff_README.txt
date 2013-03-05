
Basic (i.e., using scripts in snapshot
directory) loading of snapshot file:


1. start IDL

2. type

	ok=fload_snapshot_bh("path_to_snapshot_file",48)

3. Now you can load the particle data via

	x= fload_allstars_xyz('x')
	vx= fload_allstars_v('x')

	etc.


===================================================


The 0_slitmap*.fits files and other kinematic data
which were provided previously where generated with
the following script:


process_kinematics, "data/ds/vc3vc3e_2", 48, xlen=17.5, /skipkinemetry


It's quite long and complicated, but at least give you a starting
point to see how things are loaded and analyzed.  




===================================================
If you use the full "idlstuff" directory,
then the following items should likely be
added to your .bashrc start-up script.




export WORK=/n/scratch/hernquist_lab/tcox
export TJHOME=$HOME

module -S load hpc/IDL-7.1
if [ -f /n/sw/IDL-7.1/idl71/bin/idl_setup.bash ]; then
        source /n/sw/IDL-7.1/idl71/bin/idl_setup.bash
fi

export IDL_PATH=$IDL_PATH:$HOME/idlstuff
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/agn_spectrum
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/attenuation
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/file
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/image
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/load
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/mpfits
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/plot
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/process
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/textoidl
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/tools
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/astron/pro/astro
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/astron/pro/fits
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/astron/pro/fits_table
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/astron/pro/fits_bintable
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/astron/pro/math
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/astron/pro/misc
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/astron/pro/structure
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/astron/pro/plot
export IDL_PATH=$IDL_PATH:$HOME/idlstuff/astron/pro/tv






