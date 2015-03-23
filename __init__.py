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
import contamination

### the main program:
import pipeline

### post-production:
import plotting

options = {}

def defaultOptions():
    """
	defaultOptions()
    
    Set default options.
    
    To see the defaults, run
    
    >>> wfc3_grism.defaultOptions()
    >>> wfc3_grism.showOptions()
    """

    #### delete all keywords and reset
    for key in options.keys():
        pop = options.popitem()

    ### the observation program the reduction is being used for,
    ### currently either 3D-HST, FRONTIER_FIELDS, FIGS
    options['OBS'] = ''
    
    ### the root directory for a 3D-HST reduction:
    options['ROOT_DIR'] = '/disk1/fc/FIGS/tests/pipeline_test'

    ### the name of the wfc3 grism grating being used:
    options['GRISM_NAME'] = 'G141'
    options['DIRECT_FILTER'] = 'F140W'
    options['MAG_ZEROPOINT'] = 26.46

    ### an image used to align the raw 3dhst data:
    options['ALIGN_IMAGE'] = '/disk1/fc/FIGS/tests/candels_mosaics/gs_all_candels_ers_udf_f160w_v0.5_drz.fits'

    ### provide single flt image if not multiple images to drizzle etc.
    options['SINGLE_FLT_DIRECT'] = None
    options['SINGLE_FLT_FILTER'] = 'F125W'

    ### directory where all config files are kept and unchanged, copy config files from
    ### here for a grism reduction:
    options['GENERAL_CONFIG_FILE_DIRECTORY'] = '/disk1/fc/grism_pipeline/g141_conf'

    ### configuration options
    options['CONFIG_FILE'] = 'WFC3.IR.G141.V2.0.conf'
    options['SKY_BACKGROUND'] = None
    options['DRZRESOLA'] = '46.5'

    ### parameters to input for final multidrizzle:
    options['DRZSCALE'] = '0.06'
    options['PIXFRAC'] = '0.8'

    ### location of the 3D-HST master sky backgrounds and wheter to use them.
    ### if options['MASTER_BACKGROUND_SUBTRACTION'] = False use the default aXe
    ### background subtraction
    options['CUSTOM_MASTER_BACKGROUND_SUBTRACTION'] = True
    options['MASTER_SKIES'] = ['sky_cosmos.fits', 'sky_goodsn_hi.fits', 'sky_goodsn_lo.fits', 'sky_goodsn_vhi.fits']

    ### the sigma used for estimating segmentation map used in background
    ### subtractions, will generalyl chage for different reductions (e.g. 3D-HST, FIGS)
    ### so is worth optimizing.
    options['BACKGROUND_SEGMAP_SIGMA'] = 1.0

    ### option to input a pre-made source catalog to the reduction pipeline:
    options['PRE_MADE_INPUT_CATALOGUE'] = None
    options['PRE_MADE_INPUT_CATALOGUE_FILTER'] = 'F1537W'
    options['PRE_MADE_SEGMENTATION_MAP'] = None

    ### the detection image and band for making source catalogue:
    options['DETECTION_IMAGE'] = '/disk1/fc/FIGS/tests/candels_mosaics/gs_all_candels_ers_udf_f160w_v0.5_drz.fits'
    options['DETECTION_WHT'] = '/disk1/fc/FIGS/tests/candels_mosaics/gs_all_candels_ers_udf_f160w_v0.5_wht.fits'
    options['DETECTION_BAND'] = 'F160W'
    options['DETECTION_MAG_ZEROPOINT'] = 25.95
    options['DETECTION_FILTWAVE'] = 1535.
    options['USE_DETECTION_IMAGE_FOR_ANALYSIS'] = True

    ### general source detection parameters:
    options['DETECT_THRESH'] = 20     ## Default 2.5
    options['ANALYSIS_THRESH']  = 20  ## Default 2.5
    options['FILTWAVE'] = 1392.
    options['DETECT_MINAREA'] = 5 ## Default 2.5
    options['BACK_SIZE'] = 256.
    options['DEBLEND_MINCONT'] = 0.01

    ### bands for the contamination estimation:
    options['APPLY_FLUXCUBE_MODEL'] = True
    options['FLUXCUBE_FILTERS_DIR'] = '/disk1/fc/FIGS/tests/candels_mosaics'
    options['FLUXCUBE_FILTERS'] = [['/disk1/fc/FIGS/tests/candels_mosaics/gs_all_candels_ers_udf_f105w_v0.5_drz.fits', 'F125W', 1248.6, 26.23],
                                   ['/disk1/fc/FIGS/tests/candels_mosaics/gs_all_candels_ers_udf_f125w_v0.5_drz.fits', 'F160W', 1536.9, 25.95],
                                   ['/disk1/fc/FIGS/tests/candels_mosaics/gs_all_candels_ers_udf_f160w_v0.5_drz.fits', 'F105W', 1055.2, 26.27]]

    #### Limiting (direct) magnitude for objects run through the grism 
    #### reduction.  This is useful for cases where you want a very low 
    #### DETECT_THRESH to get good segmentation images but don't want to
    #### extract spectra for all of the faint sources.
    options['LIMITING_MAGNITUDE'] = 26.

    #### This determines what fraction of the flux must the contamination
    #### exceed to be considered in the reduction. e.g if set to 0.1 then
    #### if contamination is less than 10% of measured flux it will not be
    #### considered.
    options['LIMITING_CONTAM'] = 0.01

    #### Axe extraction geometry to get objects outside the edges
    #### of the grism images.  Currently doesn't work correctly and 
    #### the model is offset w.r.t. the observed spectra. If you want
    #### these objects need to add a border around the edge of the 
    #### direct image:
    options['AXE_EDGES'] = "300,300,0,0"
    options['ADD_IMAGE_BORDER'] = True

    #### aXe extraction geometry
    #### currently set slitless_geom=NO, orient=NO in aXecore
    #### to get the 2D spectra to line up with the orientation
    #### of the direct thumbnail.
    options['FULL_EXTRACTION_GEOMETRY'] = False
    #### aXe adjust sensitivity - convolve grism throughput with source profile
    options['AXE_ADJ_SENS'] = True
    #### aXe extract with "optimal weights"
    options['AXE_OPT_EXTR'] = True
    
# set the default options    
defaultOptions()

def showMessage(msg, warn=False):
    """
    showMessage(msg)
    
    Print a system message formatted like:
    
    ***********************************
    *  WFC3_GRISM.%s.%s.%s`module`.`function`  *
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
    substr = 'WFC3_GRISM.%s.%s.%s' %(options['OBS'], module_name, calling_function_name)
    
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

showMessage('PARAMETERS INITIALIZED')

def showOptions(outfile=None):
    """
    printOptions()
    
    Show the current option set.
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
        fp.write('###       %s     ###\n' %(options['OBS']))
        fp.write('###                                 ###\n')
        fp.write('###    %s     ###\n' %time.asctime())
        fp.write('###                                 ###\n')
        fp.write('#######################################\n')
        for key in options.keys():
            fp.write('%s = %s\n' %(key,str(options[key])))
        fp.close()

    print '\n'