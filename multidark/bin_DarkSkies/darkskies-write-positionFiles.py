from DarkSkies import *

box = DarkSkiesSimulation()

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin = 30., mmin=3*3.9*10**(12.0) )

