### pipeline import:
import wfc3_grism

### native imports:
import os
import shutil
import time
import glob

### pyraf/iraf imports:
import pyraf
from pyraf import iraf
from iraf import stsdas, dither, slitless, axe

### non-native python imports:
import numpy as np
from astropy.io import fits

### aXe specific imports:
import aXe2html.sexcat.sextractcat

def reduction_script(asn_grism=None, asn_direct=None):
	"""
	process_grism.reduction_script(asn_grism="ibhm39030_asn.fits", 
								   asn_direct="ibhm39040_asn.fits")

	Pipeline to process a set of grism/direct exposures.
	"""

	#### first exit if no input files are supplied:
	if asn_grism is None:
	    wfc3_grism.showMessage("No ASN grism file supplied",warn=True)
	    return False
	if asn_direct is None:
	    wfc3_grism.showMessage("No ASN driect file supplied: Using single _flt image for aligment", warn=True)
	    if not wfc3_grism.options['SINGLE_FLT_DIRECT']:
	    	wfc3_grism.showMessage("No direct image for alignment. Shutdown", warn=True)
	    	return False

	### make holders for the original ASN filenames used for copying
	### over the original files:
	if asn_direct:
		wfc3_grism.options['ASN_DIRECT'] = asn_direct
	wfc3_grism.options['ASN_GRISM'] = asn_grism

	wfc3_grism.showMessage('STAGE I: PREPARING ENVIRONMENT, DIRECTORS AND FILES FOR REDUCTION')

	### make sure we're in the root directory of the reduction:
	os.chdir(wfc3_grism.options['ROOT_DIR'])

	### clean up and files there from previous reductions:
	cleanup_reduction_directories()

	#### timer to check how long the aXe process is taking:
	start_time = time.time()

	### set up the aXe environment:
	set_aXe_environment(grating=wfc3_grism.options['GRISM_NAME'])

	#### check that we're in the home directory of a 3D-HST field
	check_3dhst_environment(makeDirs=True)

	### make sure the raw flt files are in the ./RAW directory
	check_flt_files_in_raw()

	### finally get all the config files in place:
	copy_over_config_files()

	### show the default rediction options:
	wfc3_grism.showOptions(outfile=None)

	#### change to the data directory:
	os.chdir('./DATA')

	#### copy over the ./RAW files into the DATA directory:
	if asn_direct:
		wfc3_grism.utils.copy_over_fresh_flt_files(asn_filename=wfc3_grism.options['ASN_DIRECT'], from_path='../RAW')
	else:
		shutil.copy('../RAW/%s' %(wfc3_grism.options['SINGLE_FLT_DIRECT']), './')

	wfc3_grism.utils.copy_over_fresh_flt_files(asn_filename=wfc3_grism.options['ASN_GRISM'], from_path='../RAW')

	### adjust the target names for the grism and direct files
	if asn_direct:
		asn_direct_file = wfc3_grism.utils.make_targname_asn(asn_direct)
	asn_grism_file = wfc3_grism.utils.make_targname_asn(asn_grism)

	### get rootnames of the grism and direct files:
	wfc3_grism.options['ROOT_GRISM'] = asn_grism_file.split('_asn.fits')[0]
	if asn_direct:
		wfc3_grism.options['ROOT_DIRECT'] = asn_direct_file.split('_asn.fits')[0]
	else:
		wfc3_grism.options['ROOT_DIRECT'] = '%s-%s' %(wfc3_grism.options['ROOT_GRISM'][:-5], wfc3_grism.options['SINGLE_FLT_FILTER'])
	
	### process the direct exposures:
	wfc3_grism.showMessage('STAGE II: PREPARE DIRECT IMAGES')
	if asn_direct:
		wfc3_grism.direct_images.process_direct_images(asn_direct_file)
	else:
		wfc3_grism.direct_images.process_single_direct_image(wfc3_grism.options['SINGLE_FLT_DIRECT'])

	### now process the grism exposures:
	wfc3_grism.showMessage('STAGE III: PREPARE GRISM IMAGES')
	wfc3_grism.grism_images.process_grism_images(asn_grism_file)

	### now make the source catalog for extraction:
	wfc3_grism.showMessage('STAGE IV: MAKING SOURCE CATALOGUE')
	wfc3_grism.source_catalogue.make_source_catalogue()

	### now check if there is a pre-made input direct catalogue,
	### since the segementaion immge for the fluxcude has been made
	### can just substitute the pre-made catalogue for the one made above:
	if wfc3_grism.options['PRE_MADE_INPUT_CATALOGUE']:
		wfc3_grism.source_catalogue.add_objects_from_premade_catalogue(extraction_image=wfc3_grism.options['IMAGE_USED_FOR_CATALOGUE'],
										                               pre_made_catalogue=wfc3_grism.options['PRE_MADE_INPUT_CATALOGUE'])

	### add the changes to sex file:
	sex_cat = wfc3_grism.sex.mySexCat('%s_drz.cat' %wfc3_grism.options['ROOT_DIRECT'])
	
	if wfc3_grism.options['DETECTION_BAND'] is None:
		aXe_filter = wfc3_grism.options['DIRECT_FILTER']
	else:
		aXe_filter = wfc3_grism.options['DETECTION_BAND']

	sex_cat.change_MAG_AUTO_for_aXe(filter=aXe_filter)

	### set up final config parameters for the grism run:
	wfc3_grism.showMessage('STAGE V: FINAL PREPARATION: TUNING CONFIGURATION FILES, MAKING AXE LIST')
	set_confguration_parameters()
	wfc3_grism.utils.make_aXe_lis(asn_grism_file, asn_direct_file)

	### now start the first aXe routine: iolprep:
	wfc3_grism.showMessage('STAGE VI: RUNNING AXE.IOLPREP')    
	wfc3_grism.utils.iraf_flpr()

	iraf.iolprep(mdrizzle_ima='%s_drz.fits' %wfc3_grism.options['ROOT_DIRECT'],
				 input_cat='%s_drz.cat' %wfc3_grism.options['ROOT_DIRECT'], 
				 dimension_in=wfc3_grism.options["AXE_EDGES"])

	### start the aXe routine: axeprep:
	wfc3_grism.showMessage('STAGE VII: RUNNING AXE.AXEPREP')    
	### cange the root directort:
	os.chdir(wfc3_grism.options['ROOT_DIR'])
	### check first whether using the 3D-HST master skies or not:

	# if wfc3_grism.options['GRISM_NAME'] == 'G141':
	# 	backim = 'WFC3.IR.G141.sky.V1.0.fits'
	# elif wfc3_grism.options['GRISM_NAME'] == 'G102':
	# 	backim = 'WFC3.IR.G102.sky.V1.0.fits'

	if wfc3_grism.options['CUSTOM_MASTER_BACKGROUND_SUBTRACTION'] == False:
		backgr=True
	else:
		backgr=False

	### now run the routine:
	wfc3_grism.utils.iraf_flpr()
	iraf.axeprep(inlist='%s_prep.lis' %(wfc3_grism.options['ROOT_GRISM']),
				 configs=wfc3_grism.options['FINAL_AXE_CONFIG'],
				 backgr=backgr, 
				 backims=wfc3_grism.options['SKY_BACKGROUND'], 
				 mfwhm=3.0,
				 norm=False)

	### now set up files for the fluxcube:
	if wfc3_grism.options['APPLY_FLUXCUBE_MODEL'] == True:
		os.chdir('./DATA')
		wfc3_grism.showMessage('STAGE VIII: SETTING UP FLUXCUBE')
		wfc3_grism.contamination.setup_fluxcube()
		os.chdir(wfc3_grism.options['ROOT_DIR'])
		cont_model="fluxcube"
	else:
		cont_model="gauss"

	### start the aXe routine: axecore:
	wfc3_grism.showMessage('STAGE IX: RUNNING AXE.AXECORE')

	wfc3_grism.utils.iraf_flpr() 

	iraf.axecore(inlist='%s_prep.lis' %(wfc3_grism.options['ROOT_GRISM']), 
				 configs=wfc3_grism.options['FINAL_AXE_CONFIG'],
				 back=False,
				 extrfwhm=4.0, 
				 drzfwhm=3.0, 
				 backfwhm=4.0,
				 slitless_geom=wfc3_grism.options['FULL_EXTRACTION_GEOMETRY'], 
				 orient=wfc3_grism.options['FULL_EXTRACTION_GEOMETRY'], 
				 exclude=False, 
	 			 lambda_mark=wfc3_grism.options['FILTWAVE'], 
				 cont_model=cont_model, 
				 model_scale=4.0, 
				 lambda_psf=wfc3_grism.options['FILTWAVE'],
				 inter_type="linear", 
				 np=10, 
				 interp=0, 
				 smooth_lengt=0,
				 smooth_fwhm=0.0,
				 spectr=False, 
				 adj_sens=wfc3_grism.options['AXE_ADJ_SENS'],
				 weights=True,
				 sampling="drizzle")   
    
	### start the aXe routine: drzprep
	wfc3_grism.showMessage('STAGE X: RUNNING AXE.DRZPREP')

	### set drizzle path here to avoid axe.axedrizzle() crashing
	### becuase of too long file names:
	os.environ['AXE_DRIZZLE_PATH'] = ('./DRIZZLE_%s' %wfc3_grism.options['GRISM_NAME'])
	wfc3_grism.options['AXE_DRIZZLE_PATH'] = ('%s/DRIZZLE_%s' %(wfc3_grism.options['ROOT_DIR'], wfc3_grism.options['GRISM_NAME']))

	wfc3_grism.utils.iraf_flpr() 

	iraf.drzprep(inlist='%s_prep.lis' %(wfc3_grism.options['ROOT_GRISM']), 
				 configs=wfc3_grism.options['FINAL_AXE_CONFIG'],
				 opt_extr=wfc3_grism.options['AXE_OPT_EXTR'], 
				 back=False)

	### start the aXe routine: axedrizzle
	wfc3_grism.showMessage('STAGE X: RUNNING AXE.AXEDRIZZLE')

	wfc3_grism.utils.iraf_flpr() 

	iraf.axedrizzle(inlist='%s_prep.lis' %(wfc3_grism.options['ROOT_GRISM']),
					configs=wfc3_grism.options['FINAL_AXE_CONFIG'],
					infwhm=4.0,
					outfwhm=3.0, 
					back=False,
					makespc=True,
					opt_extr=wfc3_grism.options['AXE_OPT_EXTR'],
					adj_sens=wfc3_grism.options['AXE_ADJ_SENS'], 
					driz_separate=False)

	### make a drizzled contamination image to test accuracy of contamination model:
	wfc3_grism.showMessage('STAGE XI: MAKING DRIZZLED CONTAMINATION IMAGE')
	wfc3_grism.multidrizzle.make_drizzled_contamination_image(asn_grism_file)

	### finally make an object table:
	wfc3_grism.utils.make_object_id_table()

	### change back to root directory 
	os.chdir(wfc3_grism.options['ROOT_DIR'])

	### make a file containing all the options:
	wfc3_grism.showOptions(outfile='./reduction_parameters.txt')

	### finally remove any leftover files:
	final_cleanup()

	# #### get end time
	end_time = time.time()
	total_time_minutes = (end_time - start_time) / 60.
	wfc3_grism.showMessage('FINISHED!! REDUCTION HAS TAKEN %.1f MINUTES' %(total_time_minutes))

def set_aXe_environment(grating='G141'):
	"""
	set_aXe_environment(grating='G141')
	    
	Setup aXe environment variables:

	AXE_IMAGE_PATH   = ./DATA
	AXE_CONFIG_PATH  = ./CONF
	AXE_OUTPUT_PATH  = ./OUTPUT_G141
	AXE_DRIZZLE_PATH = ./DRIZZLE_G141
	    
	CONF can be symlinked from e.g. /research/HST/GRISM/CONF
	"""

	grating = grating.upper()

	os.environ['AXE_IMAGE_PATH'] = wfc3_grism.options['ROOT_DIR'] + '/DATA/'
	print '--> variable AXE_IMAGE_PATH   set to "./DATA"'

	os.environ['AXE_CONFIG_PATH'] = './CONF/'
	print '--> variable AXE_CONFIG_PATH  set to "./CONF/"'
	 
	os.environ['AXE_OUTPUT_PATH'] = wfc3_grism.options['ROOT_DIR'] +  '/OUTPUT_' + grating + '/'
	print '--> variable AXE_OUTPUT_PATH  set to "./OUTPUT_' + grating + '/"'

	os.environ['AXE_DRIZZLE_PATH'] = wfc3_grism.options['ROOT_DIR'] + '/DRIZZLE_' + grating + '/'
	print '--> variable AXE_DRIZZLE_PATH set to' + '"./DRIZZLE_' + grating + '/"'

def check_3dhst_environment(makeDirs=True):
	"""
	check_3dhst_environment(makeDirs=False)
	    
	Check that all of the expected directories exist for 
	3D-HST data reduction.
	    
	If makeDirs is True, then mkdir any that isn't found in ./
	"""  

	directories = ['DATA',
				   'RAW',
				   'OUTPUT_%s' %(wfc3_grism.options['GRISM_NAME']),
				   'DRIZZLE_%s' %(wfc3_grism.options['GRISM_NAME']),
				   'CONF']

	for dir in directories:
		if not os.path.exists(dir):
			if makeDirs:
				os.mkdir(dir)
			else:
				raise IOError('Directory %s doesn\'t exist in %s.'
							  %(dir,os.getcwd()))

def cleanup_reduction_directories():

	directories = ['DATA',
				   'RAW',
				   'OUTPUT_%s' %(wfc3_grism.options['GRISM_NAME']),
				   'DRIZZLE_%s' %(wfc3_grism.options['GRISM_NAME']),
				   'CONF']

	for dir in directories:
		if os.path.exists(dir):
			shutil.rmtree(dir)
		else:
			pass

	files = glob.glob('*.lis')
	files.append('object_id_table.dat')
	files.append('reduction_parameters.txt')

	for file in files:
		if os.path.exists(file):
			os.remove(file)
		else:
			pass

def check_flt_files_in_raw():
	"""
	Checks to see if the raw flt files are in the root directory and
	if they are moves them into raw:
	"""

	### check for the raw files:
	raw_files = glob.glob('*_flt.fits')
	raw_files.extend(glob.glob('*_asn.fits'))

	if len(raw_files) == 0:
		pass
	else:
		for _file in raw_files:
			shutil.move(_file, "./RAW")

def copy_over_config_files():

	config_files = glob.glob('%s/*.fits' %(wfc3_grism.options['GENERAL_CONFIG_FILE_DIRECTORY']))
	config_files.extend(glob.glob('%s/*.conf' %(wfc3_grism.options['GENERAL_CONFIG_FILE_DIRECTORY'])))

	### set the reduction configuration directory here:
	wfc3_grism.options['REDUCTION_CONFIG_FILE_DIRECTORY'] = '%s/CONF' %(wfc3_grism.options['ROOT_DIR'])

	for cfile in config_files:
		shutil.copy(cfile, wfc3_grism.options['REDUCTION_CONFIG_FILE_DIRECTORY'])

def set_confguration_parameters():

	#### Initialize parameters, update the config file in CONF
	conf = wfc3_grism.utils.ConfFile(wfc3_grism.options['CONFIG_FILE'])

	### need to scale the 0th order sensitivity curve
	### if conf.params['SENSITIVITY_B'] == 'wfc3_abscal_IRg141_0th_sens.fits'.
	### with default the 0th order contamination is being underestimated by the pipeliene
	zeroth_list = ['wfc3_abscal_IRg141_0th_sens.fits', 'WFC3.IR.G141.0th.sens.1.fits']

	if conf.params['SENSITIVITY_B'] in zeroth_list:

		zeroth_file = fits.open('%s/%s' %(conf.path, conf.params['SENSITIVITY_B']))
		zeroth_data = zeroth_file[1].data
		sens = zeroth_data.field('SENSITIVITY')
		err = zeroth_data.field('ERROR')
		scale_factor = 3.6
		sens *= scale_factor
		err *= scale_factor
		zeroth_file.writeto('%s/WFC3_G141_0th_SCALED.fits' %(conf.path), clobber=True)
		conf.params['SENSITIVITY_B'] = 'WFC3_G141_0th_SCALED.fits'

	##### Parameters for aXe
	conf.params['DRZROOT'] = wfc3_grism.options['ROOT_GRISM']
	conf.params['DRZRESOLA'] = wfc3_grism.options['DRZRESOLA']
	conf.params['DRZSCALE'] = wfc3_grism.options['DRZSCALE']
	conf.params['DRZPFRAC'] = wfc3_grism.options['PIXFRAC']
    
	#### Parameters for BEAM order extraction, all higher orders set to 10
	#### to ensure they are not extracted.
	conf.params['MMAG_EXTRACT_A'] = str(wfc3_grism.options['LIMITING_MAGNITUDE'])
	conf.params['MMAG_EXTRACT_B'] = str(10.0)
	conf.params['MMAG_EXTRACT_C'] = str(10.0)
	conf.params['MMAG_EXTRACT_D'] = str(10.0)
	conf.params['MMAG_EXTRACT_E'] = str(10.0)

	#### Contamination estimation extarction limits, set to ratio 1:10 for a 
	#### object with limiting magnitude defined by options['LIMITING_MAGNITUDE'].
	#### The equations are described in aXe handbook page 44.

	contam_mag = -2.5*np.log(wfc3_grism.options['LIMITING_CONTAM'])

	MagMarkA = wfc3_grism.options['LIMITING_MAGNITUDE']+ contam_mag
	MagMarkB = wfc3_grism.options['LIMITING_MAGNITUDE']+ contam_mag - 0.5
	MagMarkC = wfc3_grism.options['LIMITING_MAGNITUDE']+ contam_mag - 3.1
	MagMarkD = wfc3_grism.options['LIMITING_MAGNITUDE']+ contam_mag - 5.9
	MagMarkE = wfc3_grism.options['LIMITING_MAGNITUDE']+ contam_mag - 4.5

	conf.params['MMAG_MARK_A'] = str(MagMarkA)
	conf.params['MMAG_MARK_B'] = str(MagMarkB)
	conf.params['MMAG_MARK_C'] = str(MagMarkC)
	conf.params['MMAG_MARK_D'] = str(MagMarkD)
	conf.params['MMAG_MARK_E'] = str(MagMarkE)

	## Try expanding the SMFACTOR to account for different pixel scales
	## in the sensitivity smoothing.  Bug in aXe???
	if wfc3_grism.options['GRISM_NAME'] == 'G141':
		conf.params['SMFACTOR'] = '%.3f' %(0.128254/np.float(wfc3_grism.options['DRZSCALE']))

	conf.params['DQMASK'] = np.str(np.int(conf.params['DQMASK'].split()[0]) | 4096 | 2048)

	#### Workaround to get 0th order contam. in the right place for the fluxcube
	if wfc3_grism.options['CONFIG_FILE'] == 'WFC3.IR.G141.V1.0.conf':
		conf.params['BEAMB'] = '-220 220'    
    
        
	conf.writeto('%s_full.conf' %(wfc3_grism.options['ROOT_GRISM']))

	wfc3_grism.options['FINAL_AXE_CONFIG'] = '%s_full.conf' %(wfc3_grism.options['ROOT_GRISM'])

def final_cleanup():

	os.chdir('%s/DATA' %(wfc3_grism.options['ROOT_DIR']))

	### remove lefover files:
	rmfiles = ['bg.fits', 'default.conv', 'default.nnw', 'grism.conv', 'sex_stderr', 'wfc3_grism_auto.param', 'wfc3_grism_auto.sex']

	for rfile in rmfiles:
		try:
			os.remove(rfile)
		except:
			pass

	os.chdir('%s' %(wfc3_grism.options['ROOT_DIR']))
