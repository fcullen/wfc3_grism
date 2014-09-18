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

	# ### now process the grism exposures:
	# figs.showMessage('STAGE III: PREPARE GRISM IMAGES')
	# figs.grism_images.process_grism_images(asn_grism_file)

	# ### now make the source catalog for extraction:
	# figs.showMessage('STAGE IV: MAKING SOURCE CATALOGUE')
	# figs.source_catalogue.make_source_catalogue()

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

	config_files = glob.glob('%s/*.fits' %(figs.options['CONFIG_FILE_DIRECTORY']))
	config_files.extend(glob.glob('%s/*.conf' %(figs.options['CONFIG_FILE_DIRECTORY'])))

	for cfile in config_files:
		shutil.copy(cfile, '%s/CONF/' %(figs.options['ROOT_DIR']))


	