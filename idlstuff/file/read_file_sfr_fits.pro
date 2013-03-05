
pro read_file_sfr_fits, frun, $
                t_merg, sfr_max, t_sfr_max, $
                sfr_max_fp, t_sfr_max_fp, $
                sfr_max_sb, t_sfr_max_sb, $
                meansfr_initial, meansfr_sb, meansfr_final, $
                gSFR_max_fp, gSBtime_fp, gwidth_fp, $
                gSFR_max, gSBtime, gwidth, $
                tau, decay_const


;sfrfitsfile= '/raid4/tcox/'+frun+'/sfrfits.txt'
sfrfitsfile= frun+'/sfrfits.txt'

openr, 1, sfrfitsfile
junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_merg= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & sfr_max= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_sfr_max= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & sfr_max_fp= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_sfr_max_fp= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & sfr_max_sb= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_sfr_max_sb= float(tempjunk(1))
readf, 1, junk

readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & meansfr_initital= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & meansfr_sb= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & meansfr_final= float(tempjunk(1))
readf, 1, junk

readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gSFR_max_fp= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gSBtime_fp= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gwidth_fp= float(tempjunk(1))
readf, 1, junk

readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gSFR_max= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gSBtime= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gwidth= float(tempjunk(1))
readf, 1, junk

readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & tau= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & decay_const= float(tempjunk(1))
readf, 1, junk

close, 1


end


