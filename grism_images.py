import wfc3_grism

from astropy.io import fits
import astropy.stats as stats

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

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
	wfc3_grism.showMessage("MAKING GRISM SHIFTFILE")
	make_grism_shiftfile(asn_grism_file)

	#### drizzle them together. first pass use native pixel scale for cosmic ray rejection:
	#### also just use one of the shift files (x or y not really importrant here)
	wfc3_grism.showMessage('RUNNING MULTIDRIZZLE ON GRISM IMAGES')
	wfc3_grism.multidrizzle.multidrizzle_run(asn_grism_file, 
									   shiftfile='%s_final_shifts.txt' %(wfc3_grism.options['ROOT_GRISM']), 
									   pixfrac=1.0, 
									   final_scale=0.128254, 
									   driz_cr=True,
									   skysub=True,
									   updatewcs=True)

	#### blot back to original exposures (now with cosmic ray rejection and background subtraction):
	wfc3_grism.showMessage('RUNNING BLOT ON DRIZZLED GRISM IMAGE')
	wfc3_grism.multidrizzle.blot_run(asn_grism_file, 
							   drz_file='%s_drz.fits' %(wfc3_grism.options['ROOT_GRISM']),
							   is_grism=True)

	#### generate segmentation map for each grism exposure:
	wfc3_grism.showMessage('MAKING SEGEMENTAION MAPS FOR EACH GRISM EXPOSURE')
	make_grism_exposure_segmaps(asn_grism_file, sigma=5.5)

	#### copy over fresh flt files to do final sky subtraction:
	wfc3_grism.utils.copy_over_fresh_flt_files(wfc3_grism.options['ASN_GRISM'], from_path='../RAW')

	wfc3_grism.showMessage('DOING FULL GRISM SKY SUBTRACTION')

	if wfc3_grism.options['CUSTOM_MASTER_BACKGROUND_SUBTRACTION']:
		#### do the full background subtraction on the grism exposures with custom master
		#### sky backgrounds:
		asn = wfc3_grism.utils.ASNFile(asn_grism_file)
		for grism_exposure in asn.exposures:
			 grism_sky_subtraction(grism_flt='%s_flt.fits' %(grism_exposure),
			 					   grism_segmap='%s.seg.fits' %(grism_exposure),
			 					   show=True)
		### set skysub to false for the multidrizzle tasks below:
		skysub = False
	else:
		#### do the full background subtraction on the grism exposures with default master
		#### sky backgrounds:
		asn = wfc3_grism.utils.ASNFile(asn_grism_file)
		for grism_exposure in asn.exposures:
			 default_grism_sky_subtraction(grism_flt='%s_flt.fits' %(grism_exposure),
			 					   		   grism_segmap='%s.seg.fits' %(grism_exposure),
			 					   		   show=True)
		### set skysub to false for the multidrizzle tasks below:
		skysub = False


	#### clean once more for cosmic rays and then drizzle to final resolution, skip background subtraction
	#### apply the xshifts:
	wfc3_grism.showMessage('RE-RUNNING MULTIDRIZZLE WITH BACKGROUND SUBTRACTED GRISM EXPOSURES')
	wfc3_grism.multidrizzle.multidrizzle_run(asn_grism_file, 
									   shiftfile='%s_final_shifts.txt' %(wfc3_grism.options['ROOT_GRISM']), 
									   pixfrac=1.0, 
									   final_scale=0.128254, 
									   driz_cr=True,
									   skysub=skysub,
									   updatewcs=True)

	#### finally drizzle to the desired resolution
	#### apply the y-shifts:
	wfc3_grism.multidrizzle.multidrizzle_run(asn_grism_file, 
									   shiftfile='%s_final_shifts.txt' %(wfc3_grism.options['ROOT_GRISM']), 
									   pixfrac=0.8, 
									   final_scale=0.06, 
									   driz_cr=False,
									   skysub=skysub,
									   updatewcs=False)

def make_grism_shiftfile(asn_grism_file):
	"""
	make_grism_shiftfile(asn_direct, grism_direct)

	Make a shiftfile for grism exposures to match
	corresponding direct images

	Need to do both the x and y shift files for proper working
	alignment
	"""

	#### Read shiftfile and ASN table
	sf = wfc3_grism.correct_shifts.ShiftFile('%s_final_shifts.txt' %(wfc3_grism.options['ROOT_DIRECT']))
	asn = wfc3_grism.utils.ASNFile(asn_grism_file)
    
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
	sf.write('%s_final_shifts.txt' %(wfc3_grism.options['ROOT_GRISM']))

def make_grism_exposure_segmaps(asn_grism_file, sigma=0.5):
	"""
	make_segmap(root='ib3701ryq_flt', sigma=1)

	Get a segmentation image for a flt file after creating its 
	BLOT SCI and WHT images.

	DETECT_THRESH = ANALYSIS_THRESH = sigma
	"""

	### get an ASN object:
	asn = wfc3_grism.utils.ASNFile(asn_grism_file)

	### set the default SExtractor parameters:
	se = wfc3_grism.sex.SExtractor()
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
	se.options['MAG_ZEROPOINT'] = '%.2f' %wfc3_grism.options['MAG_ZEROPOINT']
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

	if wfc3_grism.options['DIRECT_FILTER'] == 'F140W':
		direct_hdu = fits.open('%suc721143i_pfl.fits' %os.environ['iref'])
	elif wfc3_grism.options['DIRECT_FILTER'] == 'F105W':
		direct_hdu = fits.open('%suc72113oi_pfl.fits' %os.environ['iref'])

	if wfc3_grism.options['GRISM_NAME'] == 'G141':
		grism_hdu = fits.open('%su4m1335mi_pfl.fits' %os.environ['iref'])
	elif wfc3_grism.options['GRISM_NAME'] == 'G102':
		grism_hdu = fits.open('%su4m1335li_pfl.fits' %os.environ['iref'])

	### modeify flat by the intrinsic gain variations in the G141 detector:
	flat = direct_hdu[1].data[5:1019,5:1019] / grism_hdu[1].data[5:1019,5:1019]

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
	sky_files = wfc3_grism.options['MASTER_SKIES']

	### get the image profile and plot:
	im_profile = get_column_skyvals(im_hdu[1].data, seg_hdu, flat=get_flat())

	### do the chi-squared fit:
	chi2 = 1.e10
	best_sky = None

	for sky_file in sky_files:

		### open sky image and get profile:
		sky_hdu = fits.open('%s/CONF/%s' %(wfc3_grism.options['ROOT_DIR'], sky_file))
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

def correct_background_residual(im_data, seg_data, grism_flt, show=True):
	"""
	Finds the peak of the backgroud histogram across all pixels,
	for the background subtracted image. 

	There normally a slight offset from zero (~ 0.005e/s).

	The final image can then be corrected for this offset
	"""

	### apply seg mask:
	seg_mask = (seg_data == 0)

	### get the background pixels:
	background = (im_data[seg_mask].flatten())
	background = background[np.isfinite(background)]

	### fit the final background:
	hist, bin_edges = np.histogram(background, bins=np.arange(-1.5, 1.5, 0.01))
	bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

	### Define model function to be used to fit to the data above:
	def gauss(x, *p):
		A, mu, sigma = p
		return A * np.exp(-(x-mu)**2/(2.*sigma**2))

	### p0 is the initial guess for the fitting coefficients
	p0 = [1., 0., 1.]
	coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)

	### get the fitted curve
	bck = np.arange(-1.5, 1.5, 0.0001)
	hist_fit = gauss(bck, *coeff)

	### find the index of fit maximum:
	idx = np.argmin(np.abs(hist_fit - max(hist_fit)))

	### get out the background redsidual
	back_res = bck[idx]

	if show():

		fig, ax = plt.subplots(figsize=(5,5))
		ax.minorticks_on()

		ax.hist(background - back_res, bins=np.arange(-1.5, 1.5, 0.01), color='k', histtype='step', lw=1.,
				label=r'$\mathrm{Final}$ $\mathrm{G141}$')

		ax.set_xlim(-1.5, 1.5)
		ax.set_yticklabels([])

		ax.set_ylabel(r'$\mathrm{N}$', fontsize=15)
		ax.set_xlabel(r'$\mathrm{e^{-}/s}$', fontsize=15)

		ax.legend(frameon=False, loc='upper left', fontsize=14)

		fig.savefig('%s_background_histogram.pdf' %(grism_flt.split('_flt.fits')[0]))

	return back_res

def grism_sky_subtraction(grism_flt, grism_segmap, stat='median', show=False):
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

	if show:

		### plt the original exposure profile:
		fig, axes = plt.subplots(figsize=(5.5,5.5), nrows=2, ncols=1)
		plt.subplots_adjust(hspace=0.1)
		for ax in axes:
			ax.set_ylabel(r'$\mathrm{e}^-$$/\mathrm{s}$', fontsize=15)
			ax.set_xlim(0, 1014)
		axes[1].set_xlabel(r'$\mathrm{xpix}$', fontsize=15)

		im_data = np.copy(flt_hdu[1].data)

		original_profile = get_column_skyvals(im_data, seg_hdu, flat=get_flat())
		xpix = np.arange(len(original_profile))

		axes[0].plot(xpix, original_profile, color='k', ls='-.')

	### find the best sky:
	best_fit_sky_file = find_best_sky(flt_hdu, seg_hdu)

	### open the best fitting sky:
	sky_hdu = fits.open('%s/CONF/%s' %(wfc3_grism.options['ROOT_DIR'], best_fit_sky_file))

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

	if show:
		im_data = np.copy(flt_hdu[1].data)
		final_profile = get_column_skyvals(im_data, seg_hdu, flat=get_flat())
		xpix = np.arange(len(original_profile))

		axes[1].plot(xpix, final_profile, color='k', ls='-')
		axes[1].set_ylim(-0.1, 0.1)

		fig.savefig('%s_background.pdf' %(grism_flt.split('_flt.fits')[0]))

	### get the final background residual to subtract:
	back_res = correct_background_residual(im_data=np.copy(flt[1].data), 
										   seg_data=seg, 
										   grism_flt=grism_flt,
										   show=True)

	### multiply back out by the flat-field:
	final_image = (flt_hdu[1].data - back_res) * get_flat()

	### re-open the grism image:
	flt_hdu = fits.open(grism_flt, mode='update')
	flt_hdu[1].data = final_image

	bad_pixels = ~np.isfinite(flt_hdu[1].data)
	flt_hdu[1].data[bad_pixels] = 1
	flt_hdu[3].data[bad_pixels] = flt_hdu[3].data[bad_pixels] | 32

	### write the new image to the grism file:
	flt_hdu.flush()

def default_grism_sky_subtraction(grism_flt, grism_segmap, stat='median', show=False):
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

	if show:

		### plt the original exposure profile:
		fig, axes = plt.subplots(figsize=(5.5,5.5), nrows=2, ncols=1)
		plt.subplots_adjust(hspace=0.1)
		for ax in axes:
			ax.set_ylabel(r'$\mathrm{e}^-$$/\mathrm{s}$', fontsize=15)
			ax.set_xlim(0, 1014)
		axes[1].set_xlabel(r'$\mathrm{xpix}$', fontsize=15)

		im_data = np.copy(flt_hdu[1].data)

		original_profile = get_column_skyvals(im_data, seg_hdu, flat=get_flat())
		xpix = np.arange(len(original_profile))

		axes[0].plot(xpix, original_profile, color='k', ls='-.')

	### open the default sky:
	if wfc3_grism.options['GRISM_NAME'] == 'G141':
		backim = 'WFC3.IR.G141.sky.V1.0.fits'
	elif wfc3_grism.options['GRISM_NAME'] == 'G102':
		backim = 'WFC3.IR.G102.sky.V1.0.fits'

	sky_hdu = fits.open('%s/CONF/%s' %(wfc3_grism.options['ROOT_DIR'], backim))

	### re-open the grism exposure:
	flt_hdu = fits.open(grism_flt)

	### divide by the flat:
	flt_hdu[1].data /= get_flat()
	sky_hdu[0].data /= get_flat()

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

	if show:
		im_data = np.copy(flt_hdu[1].data)
		final_profile = get_column_skyvals(im_data, seg_hdu, flat=get_flat())
		xpix = np.arange(len(original_profile))

		axes[1].plot(xpix, final_profile, color='k', ls='-')
		axes[1].set_ylim(-0.1, 0.1)

		fig.savefig('%s_background.pdf' %(grism_flt.split('_flt.fits')[0]))

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




