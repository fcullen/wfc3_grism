import figs
import os

import numpy as np
from astropy.io import fits

def process_direct_images(asn_direct_file):
    """
    Proccess the direct images in a grism pointing. The steps are:

    i)   Correct pointing error shifts with tweakshifts
    ii)  Run multidrizzle with native pixel scale to flag cosmic rays
    iii) Run multidrizzle to drizzle to 1/2 pixel scale
    iv)  Resgister to CANDELS mosiac to get rough area of overlap
    v)   Align drizzled image to mosaic using figs.correct_shifts.align_direct_to_reference()
    vi)  Copy over fresh _flt files and re-run multidrizzle with new shifts
    """

    #### first get the shifts between the individual direct exposures:
    figs.showMessage('RUNNING TWEAKSHIFT ON DIRECT IMAGES')
    figs.correct_shifts.run_tweakshifts_on_direct_exposures(asn_direct_file, 
                                                            verbose=True)

    #### drizzle them together:
    figs.showMessage('RUNNING MULTIDRIZZLE ON DIRECT IMAGES')

    ### first pass use native pixel scale for cosmic ray rejection:
    figs.multidrizzle.multidrizzle_run(asn_direct_file, 
                                       shiftfile='%s_initial_shifts.txt' %(figs.options['ROOT_DIRECT']), 
                                       pixfrac=1.0, 
                                       final_scale=0.128254, 
                                       driz_cr=True,
                                       skysub=True)

    ### now sample at 1/2 pixel scale for registering to CANDELS images:
    figs.multidrizzle.multidrizzle_run(asn_direct_file, 
                                       shiftfile='%s_initial_shifts.txt' %(figs.options['ROOT_DIRECT']), 
                                       pixfrac=0.8, 
                                       final_scale=0.06, 
                                       driz_cr=False,
                                       skysub=True)

    ### cut out a region of the CANDEL image to use to align the grism exposures:
    figs.showMessage('CUTTING OUT CANDELS REGION TO ALIGN')
    figs.correct_shifts.run_sregister_for_align_image(asn_direct_file, mosiac_drz=figs.options['ALIGN_IMAGE'])

    ### align dirzzled image to reference CANDELS mosaic:
    figs.showMessage('ALIGNING DIRECT IMAGE TO CANDELS')
    figs.correct_shifts.align_direct_to_reference(asn_direct_file, verbose=True)

    ### copy over the new .flt files into the DATA directory to apply new shfits:
    figs.showMessage('RE-RUNNING MULTIDRIZZLE WITH THE NEW SHIFTS')
    figs.utils.copy_over_fresh_flt_files(asn_filename=figs.options['ASN_DIRECT'], from_path='../RAW')
    figs.multidrizzle.multidrizzle_run(asn_direct_file, 
                                       shiftfile='%s_final_shifts.txt' %(figs.options['ROOT_DIRECT']),
                                       pixfrac=0.8, 
                                       final_scale=0.06, 
                                       driz_cr=False,
                                       skysub=True)

    ### now run sregister on CANDELS mosaic if using one as a detection image:
    if figs.options['DETECTION_IMAGE'] is not None:
      figs.showMessage('ADDING BORDER TO DRZ AND MAKING FINCAL CANDELS CUTOUT')
      add_border_to_drz(drz_image='%s_drz.fits' %(figs.options['ROOT_DIRECT']), extension=300.)
      figs.correct_shifts.run_sregister_for_detection_image(asn_direct_file, 
                                                            sci_image=figs.options['DETECTION_IMAGE'],
                                                            wht_image=figs.options['DETECTION_WHT'] )

    ### at this point can perform background subtraction etc. on the direct image but
    ### don't need to do this if using different detection image, it is enough that the
    ### correct WCS solution has been found to align the image to CANDELS. The CANDELS images
    ### will now be used as the detection image for the pipeline.

def add_border_to_drz(drz_image, extension=300.):
    """
    Adds a border to the drizzled image so when registering the 
    CANDELS mosaics cuts out a larger section which is needed to can
    estimate contamination from objects outside the WFC3 Grism field of view.
    """
    
    hdulist = fits.open(drz_image)
    
    for ext in 'SCI', 'WHT':

        data = hdulist[ext].data
        hdr = hdulist[ext].header
        
        newx = hdr['NAXIS1']*0.5+extension
        newy = hdr['NAXIS2']*0.5+extension
        
        hdr.update('CRPIX1', newx)
        hdr.update('CRPIX2', newy)
        
        largeim = np.zeros((hdr['NAXIS2']+2*extension, hdr['NAXIS1']+2*extension))
        largeim[int(extension):-int(extension),int(extension):-int(extension)]=data[:,:]
        
        fits.update(drz_image, largeim, hdr, ext)


























