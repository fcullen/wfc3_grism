import wfc3_grism

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm, Normalize

from astropy.io import ascii, fits

import glob
import os

def make_inspection_plots():
	"""
	Must be run from a reductions root directory, makes a plot folder if one doesn't exist
	then makes a plot for all the objects in the object_id_table.dat file.
	"""

	### make sure there's a plot directory
	if not os.path.exists('./PLOTS'):
		os.mkdir('./PLOTS')

	### get the object ID's:
	obj_id_data = ascii.read('object_id_table.dat')
	obj_ids = obj_id_data['idx']

	### change into plot idrectory:
	os.chdir('./PLOTS')

	spc_file = glob.glob('../DRIZZLE_*/*_opt.SPC.fits')[0]
	spc_hdu = fits.open(spc_file)

	### loop through each object and make a plot:
	for i, idx in enumerate(obj_ids):

		print "Doing object %d/%d" %(i+1, len(obj_ids))

		axe_data = spc_hdu['BEAM_%dA' %idx].data
		mask = (axe_data['LAMBDA'] >= 0.8e4) & (axe_data['LAMBDA'] <= 1.15e4)

		lam = axe_data['LAMBDA'][mask]
		flux = axe_data['FLUX'][mask]
		ferr = axe_data['FERROR'][mask]
		contam = axe_data['CONTAM'][mask]
		fcorr = flux - contam

		### set up the axes:
		fig = plt.figure(figsize=(8,6.5))

		### main specturm plot:
		ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.6])
		ax1.minorticks_on()

		_n = 1e-19
		ax1.plot(lam, fcorr / _n, color='k', lw=1.5, label=r'$\mathrm{3D-HST}$ $\mathrm{Spectrum}$')
		ax1.fill_between(lam, (fcorr - ferr) / _n, (fcorr + ferr) / _n, color='grey', alpha=0.15)

		ax1.plot(lam, contam / _n, color='green', lw=2., ls=':', label=r'$\mathrm{Contamination}$ $\mathrm{Estimate}$')

		#ax1.set_ylim(-0.25, 3.0)
		ax1.set_xlim(0.6e4, 1.15e4)

		### thumbnail plot:
		ax2 = fig.add_axes([0.1, 0.7, 0.2, 0.2])
		ax2.minorticks_on()
		ax2.set_xticklabels([])
		ax2.set_yticklabels([])

		### 2D specttrum plot:

		mef_file = glob.glob('../DRIZZLE_*/*_mef_ID%d*' %(idx))[0]
		twod = fits.open(mef_file)

		flux = twod['SCI'].data
		header = twod['SCI'].header
		cont = twod['CON'].data
		fluxcont = flux - cont

		minlam = 0.8e4
		maxlam = 1.15e4

		if minlam > header['CRVAL1']: 
			xmin = (header['CRPIX1'] + ((minlam-header['CRVAL1']) / header['CDELT1']))
		else:
			xmin = (header['CRPIX1'] - ((header['CRVAL1']-minlam) / header['CDELT1']))
		
		xmax = (header['CRPIX1'] + ((maxlam - header['CRVAL1']) / header['CDELT1']))

		fcorr2d = fluxcont[:, int(xmin):int(xmax)]
    
		ax3 = fig.add_axes([0.3, 0.7, 0.6, 0.2])
		ax3.minorticks_on()
		ax3.set_xticklabels([])
		ax3.set_yticklabels([])

		ax3.imshow(fcorr2d, cmap=cm.gist_yarg, norm=Normalize(), aspect='auto')

		ax1.set_xlabel(r'$\lambda$$/\AA$', fontsize=15)
		ax1.set_ylabel(r'$f_{\lambda}$', fontsize=15)
		ax1.legend(loc='best', frameon=False)

		ax3.set_title(r'$\mathrm{OBJECT}$ $d$' %(idx), fontsize=15)

		fig.savefig('./%s.pdf' %(idx))

	os.chdir('../')

