#! /usr/bin/env python

from MultiDarkHydra import *

snList=  n.array(glob.glob("/lustre/jcomparat/MultiDark/MD_0.4Gpc_Rockstar/snapshots/out*.list.bz2"))

box = MultiDarkSimulation(Lbox=400.0 , boxDir = "MD_0.4Gpc",snl =snList )

for ii in n.arange(len(box.snl)):
	box.writePositionCatalog(ii, vmin=65)


