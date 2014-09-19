### pipeline import:
import figs

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

def reduction_script(asn_grism=None, asn_direct=None, test_run=False):
	"""
	process_grism.reduction_script(asn_grism="ibhm39030_asn.fits", 
								   asn_direct="ibhm39040_asn.fits")

	Pipeline to process a set of grism/direct exposures.
	"""

	#### first exit if no input files are supplied:
	if not asn_grism:
	    figs.showMessage("No ASN grism file supplied",warn=True)
	    return False
	if not asn_direct:
	    figs.showMessage("No ASN driect file supplied", warn=True)
	    return False

	### make holders for the original ASN filenames used for copying
	### over the original files:
	figs.options['ASN_DIRECT'] = asn_direct
	figs.options['ASN_GRISM'] = asn_grism

	figs.showMessage('STAGE I: PREPARING ENVIRONMENT, DIRECTORS AND FILES FOR REDUCTION')

	### make sure we're in the root directory of the reduction:
	os.chdir(figs.options['ROOT_DIR'])

	#### timer to check how long the aXe process is taking:
	start_time = time.time()

	### set up the aXe environment:
	set_aXe_environment(grating=figs.options['GRISM_NAME'])

	#### check that we're in the home directory of a 3D-HST field
	check_3dhst_environment(makeDirs=True)

	### make sure the raw flt files are in the ./RAW directory
	check_flt_files_in_raw()

	### finally get all the config files in place:
	copy_over_config_files()

	#### change to the data directory:
	os.chdir('./DATA')

	#### copy over the ./RAW files into the DATA directory:
	figs.utils.copy_over_fresh_flt_files(asn_filename=figs.options['ASN_DIRECT'], from_path='../RAW')
	figs.utils.copy_over_fresh_flt_files(asn_filename=figs.options['ASN_GRISM'], from_path='../RAW')

	#### adjust the target names for the grism and direct files
	asn_direct_file = figs.utils.make_targname_asn(asn_direct)
	asn_grism_file = figs.utils.make_targname_asn(asn_grism)

	### get rootnames of the grism and direct files:
	figs.options['ROOT_DIRECT'] = asn_direct_file.split('_asn.fits')[0]
	figs.options['ROOT_GRISM'] = asn_grism_file.split('_asn.fits')[0]

	### process the direct exposures:
	figs.showMessage('STAGE II: PREPARE DIRECT IMAGES')
	figs.direct_images.process_direct_images(asn_direct_file)

	### now process the grism exposures:
	figs.showMessage('STAGE III: PREPARE GRISM IMAGES')
	figs.grism_images.process_grism_images(asn_grism_file)

	### now make the source catalog for extraction:
	figs.showMessage('STAGE IV: MAKING SOURCE CATALOGUE')
	figs.source_catalogue.make_source_catalogue()

	### now set up files for the fluxcuve:
	figs.showMessage('STAGE V: SETTING UP FLUXCUBE')
	figs.contamination.setup_fluxcube()

	### set up final config parameters for the grism run:
	figs.showMessage('STAGE VI: FINAL PREPARATION: TUNING CONFIGURATION FILES, MAKING AXE LIST')
	set_confguration_parameters()
	figs.utils.make_aXe_lis(asn_grism_file, asn_direct_file)

	### now start the first aXe routine: iolprep:
	figs.showMessage('STAGE VII: RUNNING AXE.IOLPREP')    
	figs.utils.iraf_flpr()

	iraf.iolprep(mdrizzle_ima='%s_drz.fits' %figs.options['ROOT_DIRECT'],
				 input_cat='%s_drz.cat' %figs.options['ROOT_DIRECT'], 
				 dimension_in=figs.options["AXE_EDGES"])

	### start the aXe routine: axeprep:
	figs.showMessage('STAGE VIII: RUNNING AXE.AXEPREP')    
	### cange the root directort:
	os.chdir(figs.options['ROOT_DIR'])
	### check first whether using the 3D-HST master skies or not:
	if figs.options['MASTER_BACKGROUND_SUBTRACTION'] == True:
		backgr = False
		backim = ''
	else:
		backgr = True
		backim = 'WFC3.IR.G141.sky.V1.0.fits'
	### now run the routine:
	figs.utils.iraf_flpr()
	iraf.axeprep(inlist='%s_prep.lis' %(figs.options['ROOT_GRISM']),
				 configs=figs.options['FINAL_AXE_CONFIG'],
				 backgr=backgr, 
				 backims=backim, 
				 mfwhm=3.0,
				 norm=False)

	### start the aXe routine: axecore:
	figs.showMessage('STAGE IX: RUNNING AXE.AXECORE')

	figs.utils.iraf_flpr() 

	iraf.axecore(inlist='%s_prep.lis' %(figs.options['ROOT_GRISM']), 
				 configs=figs.options['FINAL_AXE_CONFIG'],
				 back=False,
				 extrfwhm=4.0, 
				 drzfwhm=3.0, 
				 backfwhm=4.0,
				 slitless_geom=figs.options['FULL_EXTRACTION_GEOMETRY'], 
				 orient=figs.options['FULL_EXTRACTION_GEOMETRY'], 
				 exclude=False, 
	 			 lambda_mark=figs.options['FILTWAVE'], 
				 cont_model="fluxcube", 
				 model_scale=4.0, 
				 lambda_psf=figs.options['FILTWAVE'],
				 inter_type="linear", 
				 np=10, 
				 interp=0, 
				 smooth_lengt=0,
				 smooth_fwhm=0.0,
				 spectr=False, 
				 adj_sens=figs.options['AXE_ADJ_SENS'],
				 weights=True,
				 sampling="drizzle")   
                 
	### change back to root directory 
	os.chdir(figs.options['ROOT_DIR'])

	#### get end time
	end_time = time.time()
	total_time_minutes = (end_time - start_time) / 60.
	figs.showMessage('FINISHED!! REDUCTION HAS TAKEN %.1f MINUTES' %(total_time_minutes))

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

	os.environ['AXE_IMAGE_PATH'] = figs.options['ROOT_DIR'] + '/DATA/'
	print '--> variable AXE_IMAGE_PATH   set to "./DATA"'

	os.environ['AXE_CONFIG_PATH'] = './CONF/'
	print '--> variable AXE_CONFIG_PATH  set to "./CONF/"'
	 
	os.environ['AXE_OUTPUT_PATH'] = figs.options['ROOT_DIR'] +  '/OUTPUT_' + grating + '/'
	print '--> variable AXE_OUTPUT_PATH  set to "./OUTPUT_' + grating + '/"'

	os.environ['AXE_DRIZZLE_PATH'] = figs.options['ROOT_DIR'] + '/DRIZZLE_' + grating + '/'
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
				   'OUTPUT_%s' %(figs.options['GRISM_NAME']),
				   'DRIZZLE_%s' %(figs.options['GRISM_NAME']),
				   'CONF', 
				   'OTHERBANDS']

	for dir in directories:
		if not os.path.exists(dir):
			if makeDirs:
				os.mkdir(dir)
			else:
				raise IOError('Directory %s doesn\'t exist in %s.'
							  %(dir,os.getcwd()))

	if figs.options['SKY_BACKGROUND'] is not None:
		if not os.path.exists('CONF/'+figs.options['SKY_BACKGROUND']):
			raise IOError("options['SKY_BACKGROUND'] doesn't exist:" +   
						  "CONF/"+figs.options['SKY_BACKGROUND'])

def cleanup_reduction_directories():

	directories = ['DATA',
				   'RAW',
				   'OUTPUT_%s' %(figs.options['GRISM_NAME']),
				   'DRIZZLE_%s' %(figs.options['GRISM_NAME']),
				   'CONF', 
				   'OTHERBANDS']

	for dir in directories:
		if os.path.exists(dir):
			shutil.rmtree(dir)
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

	config_files = glob.glob('%s/*.fits' %(figs.options['GENERAL_CONFIG_FILE_DIRECTORY']))
	config_files.extend(glob.glob('%s/*.conf' %(figs.options['GENERAL_CONFIG_FILE_DIRECTORY'])))

	for cfile in config_files:
		shutil.copy(cfile, figs.options['REDUCTION_CONFIG_FILE_DIRECTORY'])

def set_confguration_parameters():

	#### Initialize parameters, update the config file in CONF
	conf = figs.utils.ConfFile(figs.options['CONFIG_FILE'])

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
	conf.params['DRZROOT'] = figs.options['ROOT_GRISM']
	conf.params['DRZRESOLA'] = figs.options['DRZRESOLA']
	conf.params['DRZSCALE'] = figs.options['DRZSCALE']
	conf.params['DRZPFRAC'] = figs.options['PIXFRAC']
    
	#### Parameters for BEAM order extraction, all higher orders set to 10
	#### to ensure they are not extracted.
	conf.params['MMAG_EXTRACT_A'] = str(figs.options['LIMITING_MAGNITUDE'])
	conf.params['MMAG_EXTRACT_B'] = str(10.0)
	conf.params['MMAG_EXTRACT_C'] = str(10.0)
	conf.params['MMAG_EXTRACT_D'] = str(10.0)
	conf.params['MMAG_EXTRACT_E'] = str(10.0)

	#### Contamination estimation extarction limits, set to ratio 1:10 for a 
	#### object with limiting magnitude defined by options['LIMITING_MAGNITUDE'].
	#### The equations are described in aXe handbook page 44.

	contam_mag = -2.5*np.log(figs.options['LIMITING_CONTAM'])

	MagMarkA = figs.options['LIMITING_MAGNITUDE']+ contam_mag
	MagMarkB = figs.options['LIMITING_MAGNITUDE']+ contam_mag - 0.5
	MagMarkC = figs.options['LIMITING_MAGNITUDE']+ contam_mag - 3.1
	MagMarkD = figs.options['LIMITING_MAGNITUDE']+ contam_mag - 5.9
	MagMarkE = figs.options['LIMITING_MAGNITUDE']+ contam_mag - 4.5

	conf.params['MMAG_MARK_A'] = str(MagMarkA)
	conf.params['MMAG_MARK_B'] = str(MagMarkB)
	conf.params['MMAG_MARK_C'] = str(MagMarkC)
	conf.params['MMAG_MARK_D'] = str(MagMarkD)
	conf.params['MMAG_MARK_E'] = str(MagMarkE)

	## Try expanding the SMFACTOR to account for different pixel scales
	## in the sensitivity smoothing.  Bug in aXe???
	if figs.options['GRISM_NAME'] == 'G141':
		conf.params['SMFACTOR'] = '%.3f' %(0.128254/np.float(figs.options['DRZSCALE']))

	conf.params['DQMASK'] = np.str(np.int(conf.params['DQMASK'].split()[0]) | 4096 | 2048)

	#### Workaround to get 0th order contam. in the right place for the fluxcube
	if figs.options['CONFIG_FILE'] == 'WFC3.IR.G141.V1.0.conf':
		conf.params['BEAMB'] = '-220 220'    
    
        
	conf.writeto('%s_full.conf' %(figs.options['ROOT_GRISM']))

	figs.options['FINAL_AXE_CONFIG'] = '%s_full.conf' %(figs.options['ROOT_GRISM'])

	