import figs
import sys

from pyraf import iraf

def make_source_catalogue():

	if figs.options['PRE_MADE_INPUT_CATALOGUE'] is not None:

		### check there's an accompanying segmentation image:
		check_segmentation_image()

		### assumes that the MAG_AUTO column has already been changed appropriately to something like MAG_F1392W
		sexCat = figs.sex.mySexCat(figs.options['PRE_MADE_INPUT_CATALOGUE'])
		sexCat.write(ROOT_DIRECT+'_drz.cat')

	else:

		### get the sextractor object with default parametesr:
		se = run_sextractor_setup()

		### first if no detection image use the direct image for anaylsis:
		if figs.options['DETECTION_IMAGE'] is None:

			sci = '%s_SCI.fits' %(figs.options['ROOT_DIRECT'])
			wht = '%s_WHT.fits' %(figs.options['ROOT_DIRECT'])

			### remove files if already exist:
			for fitsfile in [sci, wht]:
				try:
					os.remove(fitsfile)
				except:
					pass

			### separate drizzled extensions for use with SExtractor: 
			iraf.imcopy(input='%s_drz.fits[SCI]' %(figs.options['ROOT_DIRECT']), 
					    output=sci)

			iraf.imcopy(input='%s_drz.fits[WHT]' %(figs.options['ROOT_DIRECT']), 
					    output=wht)

			### set the correct WHT image:
			se.options['WEIGHT_IMAGE'] = wht
			### run SExtractor:
			status = se.sextractImage(sci)

		### then check if using the detection image as the analysis image also:
		elif figs.options['USE_DETECTION_IMAGE_FOR_ANALYSIS'] == True:

			## get the detection image make in figs.direct_images.process_direct_images()
			sci = '%s_detection_SCI.fits' %(figs.options['DETECTION_BAND'])
			wht = '%s_detection_WHT.fits' %(figs.options['DETECTION_BAND'])

			### set the correct WHT image:
			se.options['WEIGHT_IMAGE'] = wht
			### run SExtractor:
			status = se.sextractImage(sci)

		### final scenario is use the CANDELS as detection image but the direct image
		### as the analysis image:
		else:

			detection_im = '%s_detection_SCI.fits' %(figs.options['DETECTION_BAND'])
			analysis_im = '%s_SCI.fits' %(figs.options['ROOT_DIRECT'])

			### remove files if already exist:
			try:
				os.remove(analysis_im)
			except:
				pass

			### separate drizzled extensions for use with SExtractor: 
			iraf.imcopy(input='%s_drz.fits[SCI]' %(figs.options['ROOT_DIRECT']), 
					    output=analysis_im)

			### no weight image:
			se.options['WEIGHT_TYPE'] = None
			### run SExtractor:
			status = se.sextractImage(detectionImage=detection_im,
									  analysisImage=analysis_im)

def run_sextractor_setup():

	### get the SExtractor objecT:
	se = figs.sex.SExtractor()

	### set the output parameters required for aXe 
	se.aXeParams()
 
	### copy over .conv file:
	se.copyConvFile()

	### set up the parameters:
	se.overwrite = True
	se.options['CATALOG_NAME']    = '%s_drz.cat' %(figs.options['ROOT_DIRECT'])
	se.options['CHECKIMAGE_NAME'] = '%s_seg.fits' %(figs.options['ROOT_DIRECT'])
	se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
	se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'

	se.options['FILTER']    = 'Y'
	se.options['DETECT_MINAREA'] = str(figs.options['DETECT_MINAREA'])

	se.options['DETECT_THRESH']    = str(figs.options['DETECT_THRESH']) 
	se.options['ANALYSIS_THRESH']  = str(figs.options['ANALYSIS_THRESH']) 
	se.options['MAG_ZEROPOINT'] = str(figs.options['MAG_ZEROPOINT'])
	se.options['DETECT_MINAREA'] = str(figs.options['DETECT_MINAREA'])

	se.options['BACK_SIZE'] = str(figs.options['BACK_SIZE'])
	se.options['DEBLEND_MINCONT'] = str(figs.options['DEBLEND_MINCONT'])

	return se

def check_segmentation_image():
	if figs.options['PRE_MADE_SEGMENTATION_MAP'] is None:
		sys.exit("PRE MADE INPUT CATALOGUE DOES NOT HAVE ACCOMPANYING SEG MAP!")
	else:
		pass