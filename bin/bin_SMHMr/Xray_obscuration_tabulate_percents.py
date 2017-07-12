import json                     
frac = lambda x : (x-0.4)/0.6

from scipy.interpolate import interp2d

with open("fieldsummary_comparison_obscfrac_merged_L.json") as json_file:
    json_data = json.load(json_file)

import matplotlib.pyplot as p
import numpy as n


zmin=n.array([0.50, 1.50, 2.7])
zmax=n.array([0.75, 2.10, 4.0])

#zmin=n.array([0.50, 1.50, 2.5])
#zmax=n.array([1., 2.0, 4.0])

zmean = 0.5*(zmin+zmax)

x1 = (n.transpose(json_data['0']['1'])[0]+n.transpose(json_data['0']['1'])[1])*0.5
y1 = n.ones_like(n.transpose(json_data['0']['1'])[-1])*zmean[0]
z1 = frac(n.transpose(json_data['0']['1'])[-1])
n.savetxt('fraction_zgt0_zlt10.txt', n.transpose([x1,z1]), header='LX fraction')

x2 = (n.transpose(json_data['0']['2'])[0]+n.transpose(json_data['0']['2'])[1])*0.5
y2 = n.ones_like(n.transpose(json_data['0']['2'])[-1])*zmean[1]
z2 = frac(n.transpose(json_data['0']['2'])[-1])
n.savetxt('fraction_zgt10_zlt25.txt', n.transpose([x2,z2]), header='LX fraction')

x3 = (n.transpose(json_data['0']['3'])[0]+n.transpose(json_data['0']['3'])[1])*0.5
y3 = n.ones_like(n.transpose(json_data['0']['3'])[-1])*zmean[2]
z3 = frac(n.transpose(json_data['0']['3'])[-1])
n.savetxt('fraction_zgt25_zlt50.txt', n.transpose([x3,z3]), header='LX fraction')

sys.exit()






x4 = (n.transpose(json_data['0']['1'])[0]+n.transpose(json_data['0']['1'])[1])*0.5
y4 = n.ones_like(n.transpose(json_data['0']['1'])[-1])*zmax[0]
z4 = frac(n.transpose(json_data['0']['1'])[-1])

x5 = (n.transpose(json_data['0']['2'])[0]+n.transpose(json_data['0']['2'])[1])*0.5
y5 = n.ones_like(n.transpose(json_data['0']['2'])[-1])*zmax[1]
z5 = frac(n.transpose(json_data['0']['2'])[-1])

x6 = (n.transpose(json_data['0']['3'])[0]+n.transpose(json_data['0']['3'])[1])*0.5
y6 = n.ones_like(n.transpose(json_data['0']['3'])[-1])*zmax[2]
z6 = frac(n.transpose(json_data['0']['3'])[-1])

it1 = interp2d(n.hstack((x1,x2,x3,x4,x5,x6)), n.hstack((y1,y2,y3,y4,y5,y6)), n.hstack((z1,z2,z3,z4,z5,z6)))


x1 = (n.transpose(json_data['0']['1'])[0]+n.transpose(json_data['0']['1'])[1])*0.5
y1 = n.ones_like(n.transpose(json_data['0']['1'])[-2])*zmean[0]
z1 = frac(n.transpose(json_data['0']['1'])[-2])

x2 = (n.transpose(json_data['0']['2'])[0]+n.transpose(json_data['0']['2'])[1])*0.5
y2 = n.ones_like(n.transpose(json_data['0']['2'])[-2])*zmean[1]
z2 = frac(n.transpose(json_data['0']['2'])[-2])

x3 = (n.transpose(json_data['0']['3'])[0]+n.transpose(json_data['0']['3'])[1])*0.5
y3 = n.ones_like(n.transpose(json_data['0']['3'])[-2])*zmean[2]
z3 = frac(n.transpose(json_data['0']['3'])[-2])

it2 = interp2d(n.array([x1,x2,x3]), n.array([y1,y2,y3]), n.array([z1,z2,z3]))


X,Y=n.meshgrid(n.arange(39,45,0.1),n.arange(0,5,0.1))
Z = it1(n.arange(39,45,0.1),n.arange(0,5,0.1))

p.figure(0)
# compton thin fraction '0''1,2,3' 
p.plot(n.transpose(json_data['0']['1'])[0], n.transpose(json_data['0']['1'])[-1], label='01')
p.plot(n.transpose(json_data['0']['2'])[0], n.transpose(json_data['0']['2'])[-1], label='02')
p.plot(n.transpose(json_data['0']['3'])[0], n.transpose(json_data['0']['3'])[-1], label='03')
# compton thick fraction '1''1,2,3'
p.plot(n.transpose(json_data['1']['1'])[0], n.transpose(json_data['1']['1'])[-1], label='11')
p.plot(n.transpose(json_data['1']['2'])[0], n.transpose(json_data['1']['2'])[-1], label='12')
p.plot(n.transpose(json_data['1']['3'])[0], n.transpose(json_data['1']['3'])[-1], label='13')
p.legend(frameon=False, loc=0)

p.figure(1)
p.scatter(X,Z,s=20,c=Y)
p.colorbar()
p.show()

