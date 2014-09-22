import figs

from astropy.io import fits
import astropy.stats as stats

import numpy as np

import os

def process_grism_images(asn_grism_file):
	"""
	Proccess the grism images in a grism pointing. The steps are:

	i)   Copy the direct image shiftfile to apply to grism image
	ii)  Run multidrizzle with native pixel scale to flag cosmic rays and do intial sky-subtraction
	iii) Blot back and make segmentation maps for each grism exposure.
	ii)  Loop through the inividual grism exposures and:
		 	-- make a segmentation map of the grism field to use in background subtraction
	        -- divide the exposure by the F140W flat-field
	        -- find best fitting "master" sky background and subtract
    iv)  Run multidrizzle to drizzle to 1/2 pixel scale
	"""

	### first make a shiftfile for the grism exposures based on the direct shiftfile
	figs.showMessage("MAKING GRISM SHIFTFILE")
	make_grism_shiftfile(asn_grism_file)

	#### drizzle them together. first pass use native pixel scale for cosmic ray rejection:
	figs.showMessage('RUNNING MULTIDRIZZLE ON GRISM IMAGES')
	figs.multidrizzle.multidrizzle_run(asn_grism_file, 
									   shiftfile='%s_final_shifts.txt' %(figs.options['ROOT_GRISM']), 
									   pixfrac=1.0, 
									   final_scale=0.128254, 
									   driz_cr=True,
									   skysub=True)

	#### blot back to original exposures (now with cosmic ray rejection and background subtraction):
	figs.showMessage('RUNNING BLOT ON DRIZZLED GRISM IMAGE')
	figs.multidrizzle.blot_run(asn_grism_file, 
							   drz_file='%s_drz.fits' %(figs.options['ROOT_GRISM']),
							   is_grism=True)

	#### generate segmentation map for each grism exposure:
	figs.showMessage('MAKING SEGEMENTAION MAPS FOR EACH GRISM EXPOSURE')
	make_grism_exposure_segmaps(asn_grism_file, sigma=1.5)

	#### copy over fresh flt files to do final sky subtraction:
	figs.utils.copy_over_fresh_flt_files(figs.options['ASN_GRISM'], from_path='../RAW')


	if figs.options['MASTER_BACKGROUND_SUBTRACTION']:
		#### do the full background subtraction on the grism exposures:
		figs.showMessage('DOING FULL GRISM SKY SUBTRACTION')
		asn = figs.utils.ASNFile(asn_grism_file)
		for grism_exposure in asn.exposures:
			 grism_sky_subtraction(grism_flt='%s_flt.fits' %(grism_exposure),
			 					   grism_segmap='%s.seg.fits' %(grism_exposure))
		### set skysub to false for the ultidrizzle tasks below:
		skysub = False
	else:
		### do sky subtraction in Multidrizzle
		skysub = True


	#### clean once more for cosmic rays and then drizzle to final resolution, skip background subtraction:
	figs.showMessage('RE-RUNNING MULTIDRIZZLE WITH BACKGROUND SUBTRACTED GRISM EXPOSURES')
	figs.multidrizzle.multidrizzle_run(asn_grism_file, 
									   shiftfile='%s_final_shifts.txt' %(figs.options['ROOT_GRISM']), 
									   pixfrac=1.0, 
									   final_scale=0.128254, 
									   driz_cr=True,
									   skysub=skysub)

	### finally drizzle to the desired resolution:
	figs.multidrizzle.multidrizzle_run(asn_grism_file, 
									   shiftfile='%s_final_shifts.txt' %(figs.options['ROOT_GRISM']), 
									   pixfrac=0.8, 
									   final_scale=0.06, 
									   driz_cr=False,
									   skysub=skysub)

def make_grism_shiftfile(asn_grism_file):
	"""
	make_grism_shiftfile(asn_direct, grism_direct)

	Make a shiftfile for grism exposures to match
	corresponding direct images
	"""

	#### Read shiftfile and ASN table
	sf = figs.correct_shifts.ShiftFile('%s_final_shifts.txt' %(figs.options['ROOT_DIRECT']))
	asn = figs.utils.ASNFile(asn_grism_file)
    
	if sf.nrows == len(asn.exposures):
		#### Assume one direct image for each grism images, so just
		#### change the image names in the shiftfile to the grism exposures
		for i,exp in enumerate(asn.exposures):
			sf.images[i] = exp+'_flt.fits'
	else:
		#### Have different number of grism and direct images.  Just use the 
		#### shifts/rotations for the first direct image
		xs = sf.xshift[0]
		ys = sf.yshift[0]
		rot = sf.rotate[0]
		scl = sf.scale[0]
		sf.images = []
		sf.xshift = []
		sf.yshift = []
		sf.rotate = []
		sf.scale = []
			for i,exp in enumerate(asn.exposures):
			sf.images.append(exp+'_flt.fits')
			sf.xshift.append(xs)
			sf.yshift.append(ys)
			sf.rotate.append(rot)
			sf.scale.append(scl)

		sf.nrows = len(asn.exposures)
        
	#### Write the new shiftfile
	sf.write('%s_final_shifts.txt' %(figs.options['ROOT_GRISM']))

def make_grism_exposure_segmaps(asn_grism_file, sigma=0.5):
	"""
	make_segmap(root='ib3701ryq_flt', sigma=1)

	Get a segmentation image for a flt file after creating its 
	BLOT SCI and WHT images.

	DETECT_THRESH = ANALYSIS_THRESH = sigma
	"""

	### get an ASN object:
	asn = figs.utils.ASNFile(asn_grism_file)

	### set the default SExtractor parameters:
	se = figs.sex.SExtractor()
	se.aXeParams()
	se.copyConvFile(grism=True)

	se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION, BACKGROUND'
	se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
	se.options['WEIGHT_TYPE'] = 'NONE' ### no weight map for segementation map
	se.options['FILTER']    = 'Y'
	se.options['BACK_TYPE']     = 'AUTO'
	se.options['BACK_FILTERSIZE']     = '2'
	se.options['FILTER_NAME'] = 'grism.conv'
	se.options['DETECT_THRESH']    = '%.1f' %sigma
	se.options['ANALYSIS_THRESH']  = '%.1f' %sigma
	se.options['MAG_ZEROPOINT'] = '%.2f' %figs.options['MAG_ZEROPOINT']
	se.overwrite = True

	for exp in asn.exposures:

		se.options['CATALOG_NAME']    = '%s.BLOT.SCI.cat' %(exp)
		se.options['CHECKIMAGE_NAME'] = '%s.seg.fits, bg.fits' %(exp)
		se.options['WEIGHT_IMAGE']    = '%s.BLOT.WHT.fits' %(exp)
	
		status = se.sextractImage('%s.BLOT.SCI.fits' %(exp))

def get_flat():
	"""
	Get the flat field for sky-subtraction
	"""

	f140w_hdu = fits.open('%suc721143i_pfl.fits' %os.environ['iref'])
	g141_hdu = fits.open('%su4m1335mi_pfl.fits' %os.environ['iref'])

	### modeify flat by the intrinsic gain variations in the G141 detector:
	flat = f140w_hdu[1].data[5:1019,5:1019] / g141_hdu[1].data[5:1019,5:1019]

	return flat

def get_column_skyvals(im_data, seg_hdu, flat=None):
	"""
	Get the median sky value across each column in an image
	using the segmentation map to mask actual sources.
	"""

	### get image and apply flat if required:
	im = im_data

	if flat is not None:
		im /= flat

	### get segmentation map:
	seg = seg_hdu[0].data

	### apply seg mask:
	seg_mask = (seg == 0)

	### get the profile:
	yprof = np.empty(im.shape[0])

	for i in range(im.shape[0]):
		### get the segemntation map of this column:
		colmask = seg_mask[:,i]
		### extract the sky pixles:
		sky_pixels = im[:,i]

		yprof[i] = np.median(sky_pixels)

	return yprof

def find_best_sky(im_hdu, seg_hdu):

	### get a list of the master-sky files:
	sky_files = figs.options['MASTER_SKIES']

	### get the image profile and plot:
	im_profile = get_column_skyvals(im_hdu[1].data, seg_hdu, flat=get_flat())

	### do the chi-squared fit:
	chi2 = 1.e10
	best_sky = None

	for sky_file in sky_files:

		### open sky image and get profile:
		sky_hdu = fits.open('%s/CONF/%s' %(figs.options['ROOT_DIR'], sky_file))
		sky_profile = get_column_skyvals(sky_hdu[0].data, seg_hdu, flat=None)

		### get the scale factor for sky to data:
		a = np.sum(sky_profile * im_profile) / np.sum(sky_profile*sky_profile)
		sky_profile *= a

		### get the chi2 value:
		chi2_model = np.sum(np.power(sky_profile - im_profile, 2))

		### check if the chi2 is better than current:
		if chi2_model < chi2:
			chi2 = chi2_model * 1.0
			best_sky = sky_file

	print "Best fitting sky file is: %s" %best_sky

	return best_sky

def grism_sky_subtraction(grism_flt, grism_segmap, stat='median'):
	"""
	Performs the sky-subtraction on each of the grism exposures.

	Have to divide out _flt images by the flat-field to be consistent
	with the Brammer et. al. 2011 master background images.
	('to separate out pixel-to-pixel variations from the more smoothly
    varying background')
	"""

	### open the grism exposure
	flt_hdu = fits.open(grism_flt)

	### open the segmap:
	seg_hdu = fits.open(grism_segmap)

	### find the best sky:
	best_fit_sky_file = find_best_sky(flt_hdu, seg_hdu)

	### open the best fitting sky:
	sky_hdu = fits.open('%s/CONF/%s' %(figs.options['ROOT_DIR'], best_fit_sky_file))

	### re-open the grism exposure:
	flt_hdu = fits.open(grism_flt)

	### divide by the flat:
	flt_hdu[1].data /= get_flat()

	### get the segmentation mass mask:
	seg = seg_hdu[0].data
	seg_mask = (seg == 0)

	### loop through each column in the images:
	for i in range(flt_hdu[1].data.shape[0]):

		### get the seg-map for that column
		colmask = seg_mask[:,i]
		### the master-sky in the source-free pixels:
		sky_in_col = sky_hdu[0].data[:,i][colmask]
		### the observed-sky in the source-free pixels:
		obs_sky = flt_hdu[1].data[:,i][colmask]

		### scale master to observed using median values:
		sky_hdu[0].data[:,i] *= ((np.median(obs_sky) / np.median(sky_in_col)))
		
		### subtract the scaled sky-value from the flt image:
		flt_hdu[1].data[:,i] -= sky_hdu[0].data[:,i]

	### get the profile along each column and subtract off the median:
	profile = get_column_skyvals(flt_hdu[1].data, seg_hdu, flat=None)
	flt_hdu[1].data -= np.median(profile)

	### multiply back out by the flat-field:
	final_image = flt_hdu[1].data * get_flat()

	### re-open the grism image:
	flt_hdu = fits.open(grism_flt, mode='update')
	flt_hdu[1].data = final_image

	bad_pixels = ~np.isfinite(flt_hdu[1].data)
	flt_hdu[1].data[bad_pixels] = 1
	flt_hdu[3].data[bad_pixels] = flt_hdu[3].data[bad_pixels] | 32

	### write the new image to the grism file:
	flt_hdu.flush()






