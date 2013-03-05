#!/usr/bin/env python

import arrays
import os

fdin  = '/home/brant/code/idl/colors/bc03/models/Padova1994/'
fdout = '/home/brant/code/idl/colors/bc3/parse/'
fbase = 'bc2003_hr_m'
fend  = '_ssp'

model = ['chabrier','salpeter']
zm    = ['22', '32', '42', '52', '62', '72']

for m in model:
	for z in zm:
		print m

		fdata = fdin+m+'/'
		print fdata

		if (m=='chabrier'):
			fmid = '_chab'
		else:
			fmid = '_salp'

		#do Sloan ugriz first

		fext  = '.1ABmag'

		fname = fdata+fbase+z+fmid+fend+fext
		print fname

		fd = open(fname,"r")
		fl = fd.readlines()
		fd.close()

		skip = 0
		n    = 0
		for l in fl:
			if(l.split()[0][0]!='#'):
				#print l.split()[0]
				n = n+1
			else:
				skip = skip+1
		print skip
		print n

		age  = arrays.float_array(n)
		Mbol = arrays.float_array(n)
		u_AB = arrays.float_array(n)
		g_AB = arrays.float_array(n)
		r_AB = arrays.float_array(n)
		i_AB = arrays.float_array(n)
		z_AB = arrays.float_array(n)

		i = 0
		for l in fl:
			if(l.split()[0][0]!='#'):
				age[i]  = float(l.split()[0])
				Mbol[i] = float(l.split()[1])
				g_AB[i] = float(l.split()[2])
				u_AB[i] = float(l.split()[3])+g_AB[i]
				r_AB[i] = g_AB[i]-float(l.split()[4])
				i_AB[i] = g_AB[i]-float(l.split()[5])
				z_AB[i] = g_AB[i]-float(l.split()[6])
				#print age[i],Mbol[i],g_AB[i],u_AB[i]-g_AB[i],g_AB[i]-r_AB[i],g_AB[i]-i_AB[i],g_AB[i]-z_AB[i]
				#print l.split()
				i = i+1

		#do UBVK second

		fext  = '.1color'

		fname = fdata+fbase+z+fmid+fend+fext
		print fname

		fd = open(fname,"r")
		fl = fd.readlines()
		fd.close()

		skip = 0
		n    = 0
		for l in fl:
			if(l.split()[0][0]!='#'):
				#print l.split()[0]
				n = n+1
			else:
				skip = skip+1
		print skip
		print n

		U    = arrays.float_array(n)
		B    = arrays.float_array(n)
		V    = arrays.float_array(n)
		K    = arrays.float_array(n)

		i = 0
		for l in fl:
			if(l.split()[0][0]!='#'):
				U[i] = float(l.split()[2])
				B[i] = float(l.split()[3])
				V[i] = float(l.split()[4])
				K[i] = float(l.split()[5])
				#print U[i],B[i],V[i],K[i]
				#print l.split()[2],l.split()[3],l.split()[4],l.split()[5]
				i = i+1


		#do RIJH second

		fext  = '.2color'

		fname = fdata+fbase+z+fmid+fend+fext
		print fname

		fd = open(fname,"r")
		fl = fd.readlines()
		fd.close()

		skip = 0
		n    = 0
		for l in fl:
			if(l.split()[0][0]!='#'):
				#print l.split()[0]
				n = n+1
			else:
				skip = skip+1
		print skip
		print n

		R    = arrays.float_array(n)
		I    = arrays.float_array(n)
		J    = arrays.float_array(n)
		H    = arrays.float_array(n)

		i = 0
		for l in fl:
			if(l.split()[0][0]!='#'):
				R[i] = V[i]-float(l.split()[4])
				I[i] = V[i]-float(l.split()[5])
				J[i] = V[i]-float(l.split()[6])
				H[i] = J[i]-float(l.split()[9])
				#print V[i]-R[i],V[i]-I[i],V[i]-J[i],J[i]-H[i]
				#print l.split()[4],l.split()[5],l.split()[6],l.split()[9]
				i = i+1


		#do mass third

		fext  = '.4color'

		fname = fdata+fbase+z+fmid+fend+fext
		print fname

		fd = open(fname,"r")
		fl = fd.readlines()
		fd.close()

		skip = 0
		n    = 0
		for l in fl:
			if(l.split()[0][0]!='#'):
				#print l.split()[0]
				n = n+1
			else:
				skip = skip+1
		print skip
		print n

		Mstar    = arrays.float_array(n)
		Mgas     = arrays.float_array(n)

		i = 0
		for l in fl:
			if(l.split()[0][0]!='#'):
				Mstar[i] = float(l.split()[6])
				Mgas[i]  = float(l.split()[7])
				i = i+1




		#output parsed files

		fage  = open(m+'.'+z+'.'+'age'+'.txt','w')
		fmbol = open(m+'.'+z+'.'+'mbol'+'.txt','w')
		fuAB  = open(m+'.'+z+'.'+'uAB'+'.txt','w')
		fgAB  = open(m+'.'+z+'.'+'gAB'+'.txt','w')
		frAB  = open(m+'.'+z+'.'+'rAB'+'.txt','w')
		fiAB  = open(m+'.'+z+'.'+'iAB'+'.txt','w')
		fzAB  = open(m+'.'+z+'.'+'zAB'+'.txt','w')
		fU    = open(m+'.'+z+'.'+'U'+'.txt','w')
		fB    = open(m+'.'+z+'.'+'B'+'.txt','w')
		fV    = open(m+'.'+z+'.'+'V'+'.txt','w')
		fR    = open(m+'.'+z+'.'+'R'+'.txt','w')
		fI    = open(m+'.'+z+'.'+'I'+'.txt','w')
		fJ    = open(m+'.'+z+'.'+'J'+'.txt','w')
		fH    = open(m+'.'+z+'.'+'H'+'.txt','w')
		fK    = open(m+'.'+z+'.'+'K'+'.txt','w')
		fMstar    = open(m+'.'+z+'.'+'Mstar'+'.txt','w')
		fMgas     = open(m+'.'+z+'.'+'Mgas'+'.txt','w')

		for i in range(n):
			#print i,n
			fage.write(str(age[i]))
			fmbol.write(str(Mbol[i]))
			fuAB.write(str(u_AB[i]))
			fgAB.write(str(g_AB[i]))
			frAB.write(str(r_AB[i]))
			fiAB.write(str(i_AB[i]))
			fzAB.write(str(z_AB[i]))
			fU.write(str(U[i]))
			fB.write(str(B[i]))
			fV.write(str(V[i]))
			fR.write(str(R[i]))
			fI.write(str(I[i]))
			fJ.write(str(J[i]))
			fH.write(str(H[i]))
			fK.write(str(K[i]))
			fMstar.write(str(Mstar[i]))
			fMgas.write(str(Mgas[i]))
			if(z!='72'):
				fage.write(', ')
				fmbol.write(', ')
				fuAB.write(', ')
				fgAB.write(', ')
				frAB.write(', ')
				fiAB.write(', ')
				fzAB.write(', ')
				fU.write(', ')
				fB.write(', ')
				fV.write(', ')
				fR.write(', ')
				fI.write(', ')
				fJ.write(', ')
				fH.write(', ')
				fK.write(', ')
				fMstar.write(', ')
				fMgas.write(', ')
			elif(i!=(n-1)):
				fage.write(', ')
				fmbol.write(', ')
				fuAB.write(', ')
				fgAB.write(', ')
				frAB.write(', ')
				fiAB.write(', ')
				fzAB.write(', ')
				fU.write(', ')
				fB.write(', ')
				fV.write(', ')
				fR.write(', ')
				fI.write(', ')
				fJ.write(', ')
				fH.write(', ')
				fK.write(', ')
				fMstar.write(', ')
				fMgas.write(', ')

		fage.close()
		fmbol.close()
		fuAB.close()
		fgAB.close()
		frAB.close()
		fiAB.close()
		fzAB.close()
		fU.close()
		fB.close()
		fV.close()
		fR.close()
		fI.close()
		fJ.close()
		fH.close()
		fK.close()
		fMstar.close()
		fMgas.close()

fshell = open('cat_bc03.sh','w')
fshell.write('#!/bin/bash\n')

for m in model:

	fz = open('z.'+m+'.txt','w')
	fz.write('0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05')
	fz.close()

	fshell.write('rm age.'+m+'.txt\n')
	fshell.write('cat ')
	fshell.write(m+'.'+'72'+'.'+'age.txt ')
	#for z in zm:
	#	fshell.write(m+'.'+z+'.'+'age.txt ')
	fshell.write('>> age.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.age.txt\n')

	fshell.write('rm mbol.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'mbol.txt ')
	fshell.write('>> mbol.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.mbol.txt\n')

	fshell.write('rm uAB.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'uAB.txt ')
	fshell.write('>> uAB.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.uAB.txt\n')

	fshell.write('rm gAB.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'gAB.txt ')
	fshell.write('>> gAB.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.gAB.txt\n')

	fshell.write('rm rAB.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'rAB.txt ')
	fshell.write('>> rAB.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.rAB.txt\n')

	fshell.write('rm iAB.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'iAB.txt ')
	fshell.write('>> iAB.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.iAB.txt\n')

	fshell.write('rm zAB.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'zAB.txt ')
	fshell.write('>> zAB.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.zAB.txt\n')

	fshell.write('rm U.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'U.txt ')
	fshell.write('>> U.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.U.txt\n')

	fshell.write('rm B.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'B.txt ')
	fshell.write('>> B.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.B.txt\n')

	fshell.write('rm V.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'V.txt ')
	fshell.write('>> V.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.V.txt\n')

	fshell.write('rm R.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'R.txt ')
	fshell.write('>> R.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.R.txt\n')

	fshell.write('rm I.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'I.txt ')
	fshell.write('>> I.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.I.txt\n')

	fshell.write('rm J.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'J.txt ')
	fshell.write('>> J.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.J.txt\n')

	fshell.write('rm H.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'H.txt ')
	fshell.write('>> H.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.H.txt\n')

	fshell.write('rm K.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'K.txt ')
	fshell.write('>> K.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.K.txt\n')

	fshell.write('rm Mstar.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'Mstar.txt ')
	fshell.write('>> Mstar.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.Mstar.txt\n')

	fshell.write('rm Mgas.'+m+'.txt\n')
	fshell.write('cat ')
	for z in zm:
		fshell.write(m+'.'+z+'.'+'Mgas.txt ')
	fshell.write('>> Mgas.'+m+'.txt\n')
	fshell.write('rm -f '+m+'.*.Mgas.txt\n')

os.chmod('cat_bc03.sh',0777)
