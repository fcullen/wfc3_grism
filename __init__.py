"""
A pipeline for the reduction fo the WFC3 3D-HST grism data

Based on: https://threedhst.googlecode.com/svn/threedhst/
"""

### utilities:
import terminal_controller
import utils

### wrappers:
import multidrizzle
import sex

### pipeline tasks:
import direct_images
import grism_images
import correct_shifts
import source_catalogue

### the main program:
import pipeline

options = {}

def showMessage(msg, warn=False):
    """
	showMessage(msg)
    
    Print a system message formatted like:
    
    ***********************************
    *  FIGS.`module`.`function`  *
    ***********************************
    
    `msg`
    
    ***********************************
    """       
    import os
    import sys
    import inspect
    import time
    
    calling_function_name = sys._getframe(1).f_code.co_name    
    module_name =  os.path.basename(inspect.stack()[1][1]).split('.py')[0]
    
    term = terminal_controller.TerminalController()
    
    char = '='
    substr = 'FIGS.'+module_name+'.'+calling_function_name
    
    substr = char+'  '+substr+'  '+char
    
    NL = len(substr)
    topbar = char*NL
    
    t0 = time.localtime()
    theDate = '%0d/%0d/%0d %0d:%0d' %(t0[0],t0[1],t0[2],t0[3],t0[4])
    N2 = (NL-2-len(theDate))/2
    botbar = char*N2+' '+theDate+' '+char*(NL-N2-2-len(theDate))
    
    if warn:
        text_color = term.WHITE
        bg_color = term.BG_RED
    else:
        text_color = term.BLUE
        bg_color = term.BG_WHITE
        
    print (bg_color+text_color+term.BOLD+'\n'+topbar+
           '\n'+substr+'\n'+topbar+'\n\n'+term.NORMAL+
           msg+'\n\n'+
           bg_color+text_color+term.BOLD+botbar+'\n'+term.NORMAL)

def defaultOptions():
    """
	defaultOptions()
    
    Set FIGS default options.
    
    To see the defaults, run
    
    >>> figs.defaultOptions()
    >>> figs.showOptions()
    """

    showMessage('Initializing FIGS parameters')
    
    #### delete all keywords and reset
    for key in options.keys():
        pop = options.popitem()
    
    ### the root directory for a 3D-HST reduction:
    options['ROOT_DIR'] = '/disk1/fc/FIGS/tests/pipeline_test'

    ### the name of the wfc3 grism grating being used:
    options['GRISM_NAME'] = 'G141'
    options['MAG_ZEROPOINT'] = 26.46

    ### an image used to align the raw 3dhst data:
    options['ALIGN_IMAGE'] = '/disk1/fc/FIGS/tests/candels_mosaics/gs_all_candels_ers_udf_f160w_v0.5_drz.fits'

    ### directory where all config files are kept:
    options['CONFIG_FILE_DIRECTORY'] = '/disk1/fc/FIGS/tests/conf_files'

    ### configuration options
    options['CONFIG_FILE'] = 'WFC3.IR.G141.V2.0.conf'
    options['SKY_BACKGROUND'] = None
    options['DRZRESOLA'] = '46.5'
    options['DRZSCALE'] = '0.128254'

    ### location of the 3D-HST master sky backgrounds:
    options['MASTER_BACKGROUND_SUBTRACTION'] = True
    options['MASTER_SKIES'] = ['sky_cosmos.fits', 'sky_goodsn_hi.fits', 'sky_goodsn_lo.fits', 'sky_goodsn_vhi.fits']

    ### option to input a pre-made source catalog to the reduction pipeline:
    options['PRE_MADE_INPUT_CATALOGUE'] = None
    options['PRE_MADE_SEGMENTATION_MAP'] = None

    ### the detection image and band for making source catalogue:
    options['DETECTION_IMAGE'] = '/disk1/fc/FIGS/tests/candels_mosaics/gs_all_candels_ers_udf_f160w_v0.5_drz.fits'
    options['DETECTION_WHT'] = '/disk1/fc/FIGS/tests/candels_mosaics/gs_all_candels_ers_udf_f160w_v0.5_wht.fits'
    optiosn['DETECTION_BAND'] = 'F160W'
    options['DETECTION_MAG_ZEROPOINT'] = 25.95
    options['DETECTION_FILTWAVE'] = 1535.
    options['USE_DETECTION_IMAGE_FOR_ANALYSIS'] = True

    ### general source detection parameters:
    options['DETECT_THRESH'] = 5     ## Default 2.5
    options['ANALYSIS_THRESH']  = 5  ## Default 2.5
    options['GRISM_NAME'] = 'G141'
    options['MAG_ZEROPOINT'] = 26.46
    options['FILTWAVE'] = 1392.
    options['DETECT_MINAREA'] = 5 ## Default 2.5
    options['BACK_SIZE'] = 256.
    options['DEBLEND_MINCONT'] = 0.01

# set the default options    
defaultOptions()

def showOptions(outfile=None):
    """
    printOptions()
    
    Show the current FIGS option set.
    """
    import time
    
    print '\n'

    if outfile is None:
        for key in options.keys():
            print '%s = %s' %(key,str(options[key]))
    else:
        fp = open(outfile,"w")
        fp.write('#######################################\n')
        fp.write('###                                 ###\n')
        fp.write('###    figs   %s      ###\n' %__version__)
        fp.write('###                                 ###\n')
        fp.write('###    %s     ###\n' %time.asctime())
        fp.write('###                                 ###\n')
        fp.write('#######################################\n')
        for key in options.keys():
            fp.write('%s = %s\n' %(key,str(options[key])))
        fp.close()

    print '\n'