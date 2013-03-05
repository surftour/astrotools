;----------------------------------------------------------------------
pro fit_kinematic_pa_example
;
; Usage example for the routine fit_kinematic_pa.
; Needs RDFLOAT routine from IDL astro library http://idlastro.gsfc.nasa.gov/
; Michele Cappellari, Oxford, 9 October 2007

rdfloat, 'fit_kinematic_pa_example.txt', x, y, vel, SKIPLINE=3
fit_kinematic_pa, x, y, vel, pa, paErr, /DEBUG

end
;----------------------------------------------------------------------
