import ClusterScalingRelations
import numpy as n
import matplotlib.pyplot as p
import os

out_dir= '/home/comparat/data/eRoMok'

cl = ClusterScalingRelations.ClusterScalingRelations_Mantz2016()

print cl.E035, cl.E035**2.

logM500s = n.arange(13.9,15.5,0.1)

p.figure(1, (5,5))
p.plot(10**(logM500s-14), (cl.logM500_to_logMgas(logM500s, 0.35)/10**13), 'b+')
p.xlabel('M500/1e14')
p.ylabel('Mgas/1e13')
p.grid()
p.savefig(os.path.join(out_dir,'M500-Mgas.png'))
p.clf()


p.figure(1, (5,5))
p.plot(10**(logM500s-14), cl.logM500_to_kT(logM500s, 0.35), 'b+')
p.xlabel('M500/1e14')
p.ylabel('kT/keV')
p.grid()
p.savefig(os.path.join(out_dir,'M500-kT.png'))
p.clf()


p.figure(1, (5,5))
p.plot(10**(logM500s-14), (cl.logM500_to_L(logM500s, 0.35)/10**44/cl.E035**2.), 'b+')
p.xlabel('M500/1e14')
p.ylabel('L/1e44')
p.grid()
p.savefig(os.path.join(out_dir,'M500-L.png'))
p.clf()


p.figure(1, (5,5))
p.plot(10**(logM500s-14), (cl.logM500_to_Lce(logM500s, 0.35)/10**44/cl.E035**2.), 'b+')
p.xlabel('M500/1e14')
p.ylabel('Lce/1e44')
p.grid()
p.savefig(os.path.join(out_dir,'M500-Lce.png'))
p.clf()
