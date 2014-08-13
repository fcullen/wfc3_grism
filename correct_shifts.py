import figs

from pyraf import iraf

import drizzlepac
from drizzlepac import tweakreg

yes = iraf.yes
no = iraf.no

def run_tweakreg(asn_direct):
	"""
	Use the drizzlepac software tweakreg to compute offsets between
	the direct images.
	"""

	root = asn_direct.split('_asn.fits')[0]

	### first unlearn:
	iraf.unlearn('tweakreg')
	iraf.unlearn('imagefindpars')

	### run tweakreg in non-interactive mode:
	tweakreg.TweakReg(input='*_flt.fits', updatehdr=yes, wcsname='TWEAK', 
					  conv_width=3.5, threshold=4.0, nsigma=1.5,
					  peakmin=0, peakmax=2000, writecat=no, clean=yes,
					  verbose=yes, runfile='%s_tweakreg.log' %(root), headerlet=no,
					  shiftfile=yes, outshifts='%s_shifts.txt' %(root), 
					  outwcs='%s_wcs.fits' %(root), minobj=15, 
					  searchrad=1.0, searchunits='arcseconds', use2dhist=no,
					  see2dplot=no, separation=0.5, fitgeometry='shift', residplot='No plot',
					  ncip=3, sigma=3.0)





