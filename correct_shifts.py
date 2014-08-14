import figs

from drizzlepac import tweakreg

def run_tweakreg_on_direct_exposures(asn_direct_file, verbose=True):
	"""
	Use the drizzlepac software tweakreg to compute offsets between
	the direct images.
	"""

	root = asn_direct_file.split('_asn.fits')[0]

	### run tweakreg in non-interactive mode:
	tweakreg.TweakReg(files=asn_direct_file,
			 		  updatehdr=True, 
			 		  wcsname='TWEAK', 
					  conv_width=3.5, 
					  threshold=4.0, 
					  nsigma=1.5,
					  peakmin=0, 
					  peakmax=2000, 
					  writecat=False, 
					  clean=True,
					  verbose=verbose, 
					  runfile='%s_tweakreg.log' %(root), 
					  headerlet=False,
					  shiftfile=True, 
					  outshifts='%s_shifts.txt' %(root), 
					  outwcs='%s_wcs.fits' %(root), 
					  minobj=15, 
					  searchrad=1.0, 
					  searchunits='arcseconds', 
					  use2dhist=False,
					  see2dplot=False, 
					  separation=0.5, 
					  fitgeometry='shift', 
					  residplot='No plot',
					  nclip=3, 
					  sigma=3.0)


def run_tweakreg_to_align_to_reference(asn_direct_file, reference_drz, verbose=True):
	"""
	Takes a drizzled direct image and aligns it with a reference CANDELS mosaic of the
	field.
	"""

	### get the root name:
	root = asn_direct_file.split('_asn.fits')[0]

	### run tweakreg in non-interactive mode:
	tweakreg.TweakReg(files='%s_drz_sci.fits' %(root),
					  refimage=reference_drz,
			 		  updatehdr=True, 
			 		  wcsname='TWEAK_CANDELS_ALIGN', 
					  conv_width=3.5, 
					  threshold=4.0, 
					  nsigma=1.5,
					  peakmin=0, 
					  peakmax=2000, 
					  writecat=False, 
					  clean=True,
					  verbose=verbose, 
					  runfile='%s_refalign_tweakreg.log' %(root), 
					  headerlet=False,
					  shiftfile=True, 
					  outshifts='%s_refalign_shifts.txt' %(root), 
					  outwcs='%s_refalign_wcs.fits' %(root), 
					  minobj=15, 
					  searchrad=5.0, 
					  searchunits='arcseconds', 
					  use2dhist=False,
					  see2dplot=False, 
					  separation=0.5, 
					  fitgeometry='rscale', 
					  residplot='No plot',
					  nclip=3, 
					  sigma=3.0)




