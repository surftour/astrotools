;
; Phil's routine to:
;  routine to return the IDL directory path, s.t. all routines can call 
;  it and find their data files, etc., without any trouble. 
;
; Pretty good idea, so I'll continue on with it
;
function return_idl_routines_homedir, dummy
	;;return, '/Users/phopkins/Documents/lars_galaxies/plots/idl_routines'
	;;return, '/home/phopkins/idl_routines'
	return, '/home/tcox/Tools/idlstuff'
end
