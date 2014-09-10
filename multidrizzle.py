import figs

from pyraf import iraf

import os
import glob

def direct_multidrizzle_run(asn_direct, shiftfile):
    """
    Performs an multidrizzle run on a set of input direct images.

    Some CANDELS defaults:

    driz_cr_snr='6 3.0',
    driz_cr_scale='1.6 0.7'
    """

    ### get the root name:
    root = figs.options['ROOT_DIRECT']

    ### iraf flpr()
    figs.utils.iraf_flpr()

    ### unlearn the routine:
    iraf.unlearn('multidrizzle')

    iraf.multidrizzle(input=asn_direct,
                      shiftfile=shiftfile,
                      output ='', 
                      skysub =True, 
                      updatewcs =True,
                      static=True,
                      static_sig=4.0,
                      driz_separate=True, 
                      driz_sep_kernel='turbo',
                      median=True, 
                      blot=True, 
                      driz_cr=True,
                      driz_cr_snr='6 3.0',
                      driz_cr_scale='1.6 0.7',
                      final_scale = 0.06, 
                      final_pixfrac = 0.8,
                      driz_combine=True,
                      final_wht_type='IVM',
                      clean=False)

    ### clean up:
    clean_multidrizzle_output()

def clean_multidrizzle_output():
    """
    Addidtional step for cleaning up multidrizzle outputs
    
    Removes *single_[sci/wht].fits, *sci1_blt.fits, *flt*mask1.fits, *coeffs1.dat
    """

    rmfiles = glob.glob('*single_???.fits')
    rmfiles.extend(glob.glob('*sci[12]_blt.fits'))
    rmfiles.extend(glob.glob('*flt*mask[12].fits'))
    rmfiles.extend(glob.glob('*_med.fits'))

    if len(rmfiles) > 0:
        for rmfile in rmfiles:
            os.remove(rmfile)
    else:
        pass













