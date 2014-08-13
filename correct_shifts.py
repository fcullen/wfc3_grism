import figs

from drizzlepac import tweakreg

def run_tweakreg(asn_direct, verbose=True):
	"""
	Use the drizzlepac software tweakreg to compute offsets between
	the direct images.
	"""

	root = asn_direct.split('_asn.fits')[0]

	### run tweakreg in non-interactive mode:
	tweakreg.TweakReg(files=asn_direct, updatehdr=True, wcsname='TWEAK', 
					  conv_width=3.5, threshold=4.0, nsigma=1.5,
					  peakmin=0, peakmax=2000, writecat=False, clean=True,
					  verbose=verbose, runfile='%s_tweakreg.log' %(root), headerlet=False,
					  shiftfile=True, outshifts='%s_shifts.txt' %(root), 
					  outwcs='%s_wcs.fits' %(root), minobj=15, 
					  searchrad=1.0, searchunits='arcseconds', use2dhist=False,
					  see2dplot=False, separation=0.5, fitgeometry='shift', residplot='No plot',
					  nclip=3, sigma=3.0)





