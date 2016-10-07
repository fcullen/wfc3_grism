
import wfc3_grism
import sys

from pyraf import iraf

from astropy.io import fits
from astropy import wcs

import numpy as np

def make_source_catalogue():

	### get the sextractor object with default parametesr:
	se = run_sextractor_setup()

	### first if no detection image use the direct image for anaylsis:
	if wfc3_grism.options['DETECTION_IMAGE'] is None:

		sci = '%s_SCI.fits' %(wfc3_grism.options['ROOT_DIRECT'])
		wht = '%s_WHT.fits' %(wfc3_grism.options['ROOT_DIRECT'])

		wfc3_grism.options['IMAGE_USED_FOR_CATALOGUE'] = sci

		### remove files if already exist:
		for fitsfile in [sci, wht]:
			try:
				os.remove(fitsfile)
			except:
				pass

		### separate drizzled extensions for use with SExtractor: 
		iraf.imcopy(input='%s_drz.fits[SCI]' %(wfc3_grism.options['ROOT_DIRECT']), 
				    output=sci)

		iraf.imcopy(input='%s_drz.fits[WHT]' %(wfc3_grism.options['ROOT_DIRECT']), 
				    output=wht)

		### set the correct WHT image:
		se.options['WEIGHT_IMAGE'] = wht
		### run SExtractor:
		status = se.sextractImage(sci)

	### then check if using the detection image as the analysis image also:
	elif wfc3_grism.options['USE_DETECTION_IMAGE_FOR_ANALYSIS'] == True:

		# # make sure the wcs info in header is correct:
		# f140w_hdu = fits.open('%s_drz.fits' %(wfc3_grism.options['ROOT_DIRECT']))
		# det_hdu = fits.open('%s_detection_SCI.fits' %(wfc3_grism.options['DETECTION_BAND']), 'update')

		# det_hdu[0].header = f140w_hdu[1].header
		
		# # update the image:
		# det_hdu.flush()

		# # close the image:
		# f140w_hdu.close()
		# det_hdu.close()

		## get the detection image make in wfc3_grism.direct_images.process_direct_images()
		sci = '%s_detection_SCI.fits' %(wfc3_grism.options['DETECTION_BAND'])
		wht = '%s_detection_WHT.fits' %(wfc3_grism.options['DETECTION_BAND'])

		wfc3_grism.options['IMAGE_USED_FOR_CATALOGUE'] = sci

		### set the correct WHT image:
		se.options['WEIGHT_IMAGE'] = wht
		### run SExtractor:
		status = se.sextractImage(sci)

	### final scenario is use the CANDELS as detection image but the direct image
	### as the analysis image:
	else:

		detection_im = '%s_detection_SCI.fits' %(wfc3_grism.options['DETECTION_BAND'])
		analysis_im = '%s_SCI.fits' %(wfc3_grism.options['ROOT_DIRECT'])

		wfc3_grism.options['IMAGE_USED_FOR_CATALOGUE'] = analysis_im

		### remove files if already exist:
		try:
			os.remove(analysis_im)
		except:
			pass

		### separate drizzled extensions for use with SExtractor: 
		iraf.imcopy(input='%s_drz.fits[SCI]' %(wfc3_grism.options['ROOT_DIRECT']), 
				    output=analysis_im)

		### no weight image:
		se.options['WEIGHT_TYPE'] = None
		### run SExtractor:
		status = se.sextractImage(detectionImage=detection_im,
								  analysisImage=analysis_im)

def run_sextractor_setup():

	### get the SExtractor objecT:
	se = wfc3_grism.sex.SExtractor()

	### set the output parameters required for aXe 
	se.aXeParams()
 
	### copy over .conv file:
	se.copyConvFile()

	### set up the parameters:
	se.overwrite = True
	se.options['CATALOG_NAME']    = '%s_drz.cat' %(wfc3_grism.options['ROOT_DIRECT'])
	se.options['CHECKIMAGE_NAME'] = '%s_seg.fits' %(wfc3_grism.options['ROOT_DIRECT'])
	se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
	se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'

	se.options['FILTER']    = 'Y'
	se.options['DETECT_MINAREA'] = str(wfc3_grism.options['DETECT_MINAREA'])
	se.options['DETECT_THRESH']    = str(wfc3_grism.options['DETECT_THRESH']) 
	se.options['ANALYSIS_THRESH']  = str(wfc3_grism.options['ANALYSIS_THRESH']) 
	se.options['DETECT_MINAREA'] = str(wfc3_grism.options['DETECT_MINAREA'])

	se.options['BACK_SIZE'] = str(wfc3_grism.options['BACK_SIZE'])
	se.options['DEBLEND_MINCONT'] = str(wfc3_grism.options['DEBLEND_MINCONT'])

	### if using separate detection image take the zeropoint of that detection image, else use
	### use the zeropoint of the direct image:
	if wfc3_grism.options['DETECTION_IMAGE'] and wfc3_grism.options['USE_DETECTION_IMAGE_FOR_ANALYSIS']:
		se.options['MAG_ZEROPOINT'] = str(wfc3_grism.options['DETECTION_MAG_ZEROPOINT'])
	else:
		se.options['MAG_ZEROPOINT'] = str(wfc3_grism.options['MAG_ZEROPOINT'])

	return se

def add_objects_from_premade_catalogue(extraction_image='', pre_made_catalogue=''):
	"""
	If using a pre-made catalogue this function makes sure the X-WORLD, Y-WORLD
	poistions are aligned with the reduction image. Then adds them to the 
	catalogue produced in the reduction process
	"""

	### generate a sextractor catalog object:
	sex_cat = wfc3_grism.sex.mySexCat(pre_made_catalogue)

	### get the world co-ordinates into an array:
	world_crd = np.array([[sex_cat['X_WORLD'][i], sex_cat['Y_WORLD'][i]] for i in range(len(sex_cat['X_WORLD']))], np.float_)

	### get the WCS of the extraction iamge:
	hdu = fits.open(extraction_image)
	w = wcs.WCS(hdu[0].header)

	### get the new pixel co-ordinates:
	pix_crd = w.wcs_world2pix(world_crd, 1)

	### have to ceonvert them to strings:
	x_image =  pix_crd[:,0].astype('|S8')
	y_image =  pix_crd[:,1].astype('|S8')

	### apply to the catalogue, make a new one becuase its a faff
	### using the sextractcat functions:
	outfile = open('premade_catalog.cat', 'w')

	for line in sex_cat.headerlines:
		outfile.write(line)

	for i in range(len(sex_cat['X_WORLD'])):
		outfile.write('%s\t' %(sex_cat['NUMBER'][i]) +
					  '%s\t' %(x_image[i]) +
					  '%s\t' %(y_image[i]) +
					  '%s\t' %(sex_cat['X_WORLD'][i]) +
					  '%s\t' %(sex_cat['Y_WORLD'][i]) +
					  '%s\t' %(sex_cat['A_IMAGE'][i]) +
					  '%s\t' %(sex_cat['B_IMAGE'][i]) +
					  '%s\t' %(sex_cat['THETA_IMAGE'][i]) +
					  '%s\t' %(sex_cat['A_WORLD'][i]) +
					  '%s\t' %(sex_cat['B_WORLD'][i]) +
					  '%s\t' %(sex_cat['THETA_WORLD'][i]) +
					  '%s\n' %(sex_cat['MAG_AUTO'][i]))

	outfile.close()

	### re-open and add the premade catalogue objects to end of file
	### and apply the change the MAG_AUTO column in the catalog to be MAG_{filter} for aXe
	premade_sex_cat = wfc3_grism.sex.mySexCat('premade_catalog.cat')
	sex_cat = wfc3_grism.sex.mySexCat('%s_drz.cat' %wfc3_grism.options['ROOT_DIRECT'])

	sex_cat.rowlines += premade_sex_cat.rowlines

	# aXe_filter = wfc3_grism.options['PRE_MADE_INPUT_CATALOGUE_FILTER'] 
	# sex_cat.change_MAG_AUTO_for_aXe(filter=aXe_filter)

	### finally write out the catalogue:
	sex_cat.write('%s_drz.cat' %wfc3_grism.options['ROOT_DIRECT'])
	