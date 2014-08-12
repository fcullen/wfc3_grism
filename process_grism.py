import figs

import os
import pyfits
import pyraf
from pyraf import iraf
from iraf import stsdas, dither, slitless, axe

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

    directories = ['DATA','RAW','OUTPUT_'+figs.options['GRISM_NAME'],
                   'DRIZZLE_'+figs.options['GRISM_NAME'],'CONF',
                   'fitting', 'OTHERBANDS', 'plots', 'Thumbnails']

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

def reduction_script(asn_grism=None, asn_direct=None):
	"""
	process_grism.reduction_script(asn_grism="ibhm39030_asn.fits", 
								   asn_direct="ibhm39040_asn.fits")

	Pipeline to process a set of grism/direct exposures.
	"""

	import shutil
	import glob
	import numpy as np
	import aXe2html.sexcat.sextractcat
	import time

	#### first exit if no input files are supplied:
	if not asn_grism:
	    figs.showMessage("No ASN grism file supplied",warn=True)
	    return False
	if not asn_direct:
	    figs.showMessage("No ASN driect file supplied", warn=True)
	    return False

	#### check how long the aXe process is taking:
	start_time = time.time()

	### set up the aXe environment:
	set_aXe_environment(grating=figs.options['GRATING'])

	#### check that we're in the home directory of a 3D-HST field
	check_3dhst_environment(makeDirs=True)

	#### change to the data directory:
	os.chdir('./DATA')

	figs.showMessage('PREPARING DIRECT IMAGES')

	#### adjust the target names for the grism and diect files
	asn_direct_file = figs.utils.make_targname_asn(asn_direct)
	asn_grism_file = figs.utils.make_targname_asn(asn_grism)

	#### copy over the ./RAW files into the DATA directory:
	figs.process_direct_images.copy_over_fresh_flt_files(asn_filename=asn_direct_file,
													     from_path='%s/RAW' %(figs.options['ROOT_DIR']))

	#### now align the direct images to a referecnce file
	asn = figs.utils.ASNFile(asn_direct_file)
	flt = '%s_flt.fits' %(asn.exposures[0])
	figs.process_direct_images.align_raw_flt_to_reference(raw_flt=flt,
														  reference_image=figs.options['ALIGN_IMAGE'])




	