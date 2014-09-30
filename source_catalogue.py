import wfc3_grism
import sys

from pyraf import iraf

def make_source_catalogue():

	if wfc3_grism.options['PRE_MADE_INPUT_CATALOGUE'] is not None:

		### check there's an accompanying segmentation image:
		check_segmentation_image()

		### assumes that the MAG_AUTO column has already been changed appropriately to something like MAG_F1392W
		sexCat = wfc3_grism.sex.mySexCat(wfc3_grism.options['PRE_MADE_INPUT_CATALOGUE'])
		sexCat.write(ROOT_DIRECT+'_drz.cat')

	else:

		### get the sextractor object with default parametesr:
		se = run_sextractor_setup()

		### first if no detection image use the direct image for anaylsis:
		if wfc3_grism.options['DETECTION_IMAGE'] is None:

			sci = '%s_SCI.fits' %(wfc3_grism.options['ROOT_DIRECT'])
			wht = '%s_WHT.fits' %(wfc3_grism.options['ROOT_DIRECT'])

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

			## get the detection image make in wfc3_grism.direct_images.process_direct_images()
			sci = '%s_detection_SCI.fits' %(wfc3_grism.options['DETECTION_BAND'])
			wht = '%s_detection_WHT.fits' %(wfc3_grism.options['DETECTION_BAND'])

			### set the correct WHT image:
			se.options['WEIGHT_IMAGE'] = wht
			### run SExtractor:
			status = se.sextractImage(sci)

		### final scenario is use the CANDELS as detection image but the direct image
		### as the analysis image:
		else:

			detection_im = '%s_detection_SCI.fits' %(wfc3_grism.options['DETECTION_BAND'])
			analysis_im = '%s_SCI.fits' %(wfc3_grism.options['ROOT_DIRECT'])

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
		se.options['MAG_ZEROPOINT'] = str(wfc3_grism.options['DIRECT_MAG_ZEROPOINT'])

	return se

def check_segmentation_image():
	if wfc3_grism.options['PRE_MADE_SEGMENTATION_MAP'] is None:
		sys.exit("PRE MADE INPUT CATALOGUE DOES NOT HAVE ACCOMPANYING SEG MAP!")
	else:
		pass