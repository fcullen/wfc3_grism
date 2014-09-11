import figs

from astropy.io import fits
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
	figs.showMessage('RUNNING MULTIDRIZZLE ON DIRECT IMAGES')
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
	make_grism_exposure_segmaps(asn_grism_file, sigma=0.5)

	#### copy over fresh flt files to do final sky subtraction:
	copy_over_fresh_flt_files(figs.options['ASN_GRISM'], from_path='../RAW')

	#### do the full background subtraction on the grism exposures:
	figs.showMessage('DOING FULL GRISM SKY SUBTRACTION')
	asn = figs.utils.ASNFile(asn_grism_file)
	for grism_exposure in asn.exposures:
		 grism_sky_subtraction(grism_flt='%s_flt.fits' %(grism_exposure),
		 					   grism_segmap='%s.seg.fits' %(grism_exposure))


	#### finally clean once more for cosmic rays and then drizzle to final resolution, skip background subtraction:
	# figs.showMessage('RE-RUNNING MULTIDRIZZLE WITH BACKGROUND SUBTRACTED GRISM EXPOSURES')
	# figs.multidrizzle.multidrizzle_run(asn_grism_file, 
	# 								   shiftfile='%s_final_shifts.txt' %(figs.options['ROOT_GRISM']), 
	# 								   pixfrac=1.0, 
	# 								   final_scale=0.128254, 
	# 								   driz_cr=True,
	# 								   skysub=False)

	# figs.multidrizzle.multidrizzle_run(asn_grism_file, 
	# 								   shiftfile='%s_final_shifts.txt' %(figs.options['ROOT_GRISM']), 
	# 								   pixfrac=0.8, 
	# 								   final_scale=0.06, 
	# 								   driz_cr=False,
	# 								   skysub=False)

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

def grism_sky_subtraction(grism_flt, grism_segmap):
	"""
	Performs the sky-subtraction on each of the grism exposures.
	"""

	### open the grism exposure
	flt_hdu = fits.open(grism_flt)

	### open the segmap:
	seg_hdu = fits.open(grism_segmap)

