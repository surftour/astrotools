Function fload_1gal_com, startid, numpart


	com = [0.0,0.0,0.0]
	lastcom = [0.0,0.0,0.0]
	dcom = 100.0


	fload_1gal_all_data, startid, numpart, allmasses, allxs, allys, allzs, allvxs, allvys, allvzs

	; find center-of-mass
	com(0) = total(allmasses*allxs)/total(allmasses)
	com(1) = total(allmasses*allys)/total(allmasses)
	com(2) = total(allmasses*allzs)/total(allmasses)

	print, "gal1 center-of-mass         ",com
	return, com

end


