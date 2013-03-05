#!/bin/bash
g77 -shared -m32 -I/usr/lib/include    -c -o extra.o extra.f 
g77 -shared -m32 -I/usr/lib/include    -c -o lenact.o lenact.f 
gcc -shared -m32 -I/usr/lib/include    -c -o udmget.o udmget.c
g77 -shared -m32 -I/usr/lib/include    -c -o xspexrav.o xspexrav_callable.f
g77 -shared -m32 -I/usr/lib/include    -c -o xspexriv.o xspexriv.f
g77 -shared -m32 -I/usr/lib/include    -c -o fgunc.o  fgunc.f
g77 -shared -m32 -I/usr/lib/include    -c -o xsbexrav.o xsbexrav.f 
g77 -shared -m32 -I/usr/lib/include    -c -o xsbexriv.o xsbexriv.f 
gcc -shared -m32 -I/usr/lib/include    -c -o main.o main.c
gcc -shared -m32 -I/usr/lib/include extra.o lenact.o udmget.o xspexrav.o xspexriv.o fgunc.o xsbexriv.o xsbexrav.o main.o -lm -lg2c -lstdc++ -L/usr/lib -o marconi_agn_spectrum
