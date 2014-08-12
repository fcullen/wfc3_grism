"""
A pipeline for the reduction fo the WFC3 3D-HST grism data

Based on: https://threedhst.googlecode.com/svn/threedhst/
"""

import terminal_controller
import utils
import shifts
import astrodrizzle
import process_grism
import regions
import sex
import process_direct_images

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
    
    term = TerminalController.TerminalController()
    
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
    
    # delete all keywords and reset
    for key in options.keys():
        pop = options.popitem()
    
    # the root directory for a 3D-HST reduction:
    options['ROOT_DIR'] = '/disk1/fc/FIGS'

    # the name of the wfc3 grism grating being used:
    options['GRATING'] = 'G141'

    # an image used to align the raw 3dhst data:
    options['ALIGN_IMAGE'] = None

    # configuration options
    options['CONFIG_FILE'] = 'WFC3.IR.G141.V2.0.conf'
    options['SKY_BACKGROUND'] = None
    options['DRZRESOLA'] = '46.5'
    options['DRZSCALE'] = '0.128254'

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