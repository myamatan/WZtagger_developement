#!/bin/sh


for MSHIFT in 0 30 60 ; do
	for NBIN in 15 30 45 ; do
		root -b -q 'sigbkg_modeling.cc+("veryloose","MaxVVSemi_v5",'${NBIN}','${MSHIFT}')'
	done
done
	

