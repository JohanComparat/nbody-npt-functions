import time
t0 = time.time()

import os
import numpy as n
import sys
import glob

import cPickle

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as p

from scipy.interpolate import interp1d

L_box = 1000./0.6777

tracer_names = n.array(['S8_ELG', 'S8_LRG', 'S5_GAL', 'S8_QSO', 'S6_AGN', 'S5_BCG'])
marker_dict={'S5_BCG':'1', 'S5_GAL':'2', 'S6_AGN':'3', 'S8_LRG':'4', 'S8_ELG':'+', 'S8_QSO':'x'}
color_dict ={'S5_BCG':'r', 'S5_GAL':'r', 'S6_AGN':'c', 'S8_LRG':'k', 'S8_ELG':'b', 'S8_QSO':'g'}
p0 = n.array([[-1., -1.]])
points = {'S5_BCG':p0, 'S5_GAL':p0, 'S6_AGN':p0, 'S8_LRG':p0, 'S8_ELG':p0, 'S8_QSO':p0}

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

zs = n.arange(0.,4,0.001)
dc_2_z = interp1d(cosmoMD.comoving_distance(zs),zs)

import astropy.io.fits as fits

sf = fits.open(os.path.join(os.environ['MD10'],'output_MD_1.0Gpc.fits'))[1].data
plot_dir = '/afs/mpe/www/people/comparat/eRoMok/pie_plots/'
work_dir = os.path.join(os.environ['MD10'],'work_agn')


# redshift loop
#ii = 0
def get_slice(cpickle_dump_file, x_observer=0., y_observer=0., z_observer = 0., x_shift=0., y_shift=0., z_shift=0., slice_z_min=0., slice_z_max = 10., distance_min=0., distance_max = L_box):
    
    snap_selection = (sf['comoving_distance']<distance_max)&(sf['comoving_distance']>distance_min)

    z_all = sf['redshift'][snap_selection]
    z_boundaries = n.hstack((dc_2_z(distance_min), (z_all[1:]+z_all[:-1])/2., dc_2_z(distance_max)))

    for ii in range(len(z_all)):
	z_min, z_max = z_boundaries[ii], z_boundaries[ii+1]
	el = sf[ii+1]

	r_min, r_max = cosmoMD.comoving_distance(z_min).value, cosmoMD.comoving_distance(z_max).value

	position_files = n.array(glob.glob(os.path.join(work_dir, 'out_'+el['snap_name']+'_SAM_Nb_?.fits')))
	position_files.sort()

	# position file loop
	print r_min, r_max
	for index in range(len(position_files)):
	    print time.time()-t0
	    print position_files[index]
	    positions = fits.open(position_files[index])[1].data
	    tracer_files = n.array(glob.glob(os.path.join(work_dir, 'out_'+el['snap_name']+'_SAM_Nb_'+str(index)+'_4MOST_*.fits')))
	    tracer_files.sort()

	    # tracer loop
	    #path_2_tracer_file = tracer_files[0]
	    for path_2_tracer_file in tracer_files:
    		print path_2_tracer_file
		spl_bn = os.path.basename(path_2_tracer_file)[:-5].split('_')
		tracer_name = spl_bn[-2]+'_'+spl_bn[-1]
		ids = fits.open(path_2_tracer_file)[1].data['line_number']

		x_i = positions['x'][ids]/0.6777 - x_observer + x_shift
		y_i = positions['y'][ids]/0.6777 - y_observer + y_shift
		z_i = positions['z'][ids]/0.6777 - z_observer + z_shift

		shell = (x_i*x_i + y_i*y_i + z_i*z_i < r_max**2.) & (x_i*x_i + y_i*y_i + z_i*z_i > r_min**2.)
		slice = (shell) & (z_i>slice_z_min) &(z_i<slice_z_max)

		points[tracer_name] = n.vstack(( points[tracer_name], n.transpose([x_i[slice], y_i[slice]]) ))


    f=open(cpickle_dump_file, 'w')
    cPickle.dump(points,f)
    f.close()
    return points


points_1 = get_slice(os.path.join(work_dir, 'slice_1_Lbox.pkl'))
points_2 = get_slice(os.path.join(work_dir, 'slice_2_Lbox.pkl'),x_shift=L_box,distance_min=L_box, distance_max = L_box*2.)

p.figure(0, (4.5,4.5))
p.axes([0.17,0.17,0.78,0.78])
for tracer in tracer_names:
	x_pos, y_pos = points_1[tracer].T
	p.plot(x_pos, y_pos,marker=marker_dict[tracer],color=color_dict[tracer],rasterized=True,ls='None',label=tracer)

p.legend(loc=0, frameon=False)
p.xlabel('Mpc')
p.ylabel('Mpc')
p.xlim((0,L_box))
p.ylim((0,L_box))
p.savefig(os.path.join(plot_dir, 'slice_1_Lbox.png'))
p.clf()


p.figure(0, (4.5,4.5))
p.axes([0.17,0.17,0.78,0.78])
for tracer in tracer_names:
	x_pos, y_pos = points_2[tracer].T
	p.plot(x_pos, y_pos,marker=marker_dict[tracer],color=color_dict[tracer],rasterized=True,ls='None',label=tracer)

p.legend(loc=0, frameon=False)
p.xlabel('Mpc')
p.ylabel('Mpc')
p.xlim((L_box,2*L_box))
p.ylim((0.,L_box))
p.savefig(os.path.join(plot_dir, 'slice_2_Lbox.png'))
p.clf()
