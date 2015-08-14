import wfc3_grism
import os

import numpy as np

import matplotlib.pyplot as plt

from astropy.io import fits
import astropy.stats as stats

import shutil

def process_single_direct_image(direct_flt):

    ### re-name the _flt file as a fake _drz file to use
    ### in the correct_shifts.align_direct_to_reference() routine:
    shutil.move(direct_flt, '%s_drz.fits' %(wfc3_grism.options['ROOT_DIRECT']))

    ### cut out a region of the CANDEL image to use to align the grism exposures:
    wfc3_grism.showMessage('CUTTING OUT CANDELS REGION TO ALIGN')
    wfc3_grism.correct_shifts.run_sregister_for_align_image(mosiac_drz=wfc3_grism.options['ALIGN_IMAGE'])

    ### align to the reference mosaic:
    wfc3_grism.showMessage('ALIGNING DIRECT IMAGE TO CANDELS')
    wfc3_grism.correct_shifts.align_direct_to_reference(verbose=True, n_iter=5, drizzled_image=False)

def process_direct_images(asn_direct_file):
    """
    Proccess the direct images in a grism pointing. The steps are:

    i)   Correct pointing error shifts with tweakshifts
    ii)  Run multidrizzle with native pixel scale to flag cosmic rays
    iii) Run multidrizzle to drizzle to 1/2 pixel scale
    iv)  Resgister to CANDELS mosiac to get rough area of overlap
    v)   Align drizzled image to mosaic using wfc3_grism.correct_shifts.align_direct_to_reference()
    vi)  Copy over fresh _flt files and re-run multidrizzle with new shifts
    """

    #### first get the shifts between the individual direct exposures:
    wfc3_grism.showMessage('RUNNING TWEAKSHIFT ON DIRECT IMAGES')
    wfc3_grism.correct_shifts.run_tweakshifts_on_direct_exposures(asn_direct_file, 
                                                            verbose=True)

    #### drizzle them together:
    wfc3_grism.showMessage('RUNNING MULTIDRIZZLE ON DIRECT IMAGES')

    ### first pass use native pixel scale for cosmic ray rejection:
    wfc3_grism.multidrizzle.multidrizzle_run(asn_direct_file, 
                                       shiftfile='%s_initial_shifts.txt' %(wfc3_grism.options['ROOT_DIRECT']), 
                                       pixfrac=1.0, 
                                       final_scale=wfc3_grism.options['INSTRUMENT_PIXEL_SCALE'], 
                                       driz_cr=True,
                                       skysub=True,
                                       blot_back=True)

    ### now sample at 1/2 pixel scale for registering to CANDELS images:
    wfc3_grism.multidrizzle.multidrizzle_run(asn_direct_file, 
                                       shiftfile='%s_initial_shifts.txt' %(wfc3_grism.options['ROOT_DIRECT']), 
                                       pixfrac=wfc3_grism.options['PIXFRAC'], 
                                       final_scale=wfc3_grism.options['ALIGN_IMAGE_PIXEL_SCALE'], 
                                       driz_cr=False,
                                       skysub=True,
                                       blot_back=False)

    shutil.copy('%s_drz.fits' %(wfc3_grism.options['ROOT_DIRECT']), 'INITIAL_SHIFTED_DRZ.fits')

    ### cut out a region of the CANDEL image to use to align the grism exposures:
    wfc3_grism.showMessage('CUTTING OUT CANDELS REGION TO ALIGN')
    wfc3_grism.correct_shifts.run_sregister_for_align_image(mosiac_drz=wfc3_grism.options['ALIGN_IMAGE'])

    ### align dirzzled image to reference CANDELS mosaic:
    wfc3_grism.showMessage('ALIGNING DIRECT IMAGE TO CANDELS')
    wfc3_grism.correct_shifts.align_direct_to_reference(verbose=True)

    ### copy over the new .flt files into the DATA directory to apply new shfits:
    wfc3_grism.showMessage('RE-RUNNING MULTIDRIZZLE WITH THE NEW SHIFTS FOR COSMIC RAY REJECTION')
    wfc3_grism.utils.copy_over_fresh_flt_files(asn_filename=wfc3_grism.options['ASN_DIRECT'], from_path='../RAW')

    ### first pass cosmic-ray rejection:
    wfc3_grism.multidrizzle.multidrizzle_run(asn_direct_file, 
                                       shiftfile='%s_initial_shifts.txt' %(wfc3_grism.options['ROOT_DIRECT']), 
                                       pixfrac=1.0, 
                                       final_scale=wfc3_grism.options['INSTRUMENT_PIXEL_SCALE'], 
                                       driz_cr=True,
                                       skysub=False,
                                       updatewcs=False,
                                       blot_back=True)

    ### blot back to original exposures (now with cosmic ray rejection and background subtraction):
    wfc3_grism.showMessage('RUNNING BLOT ON DRIZZLED DIRECT IMAGE')
    wfc3_grism.multidrizzle.blot_run(asn_direct_file, 
                               drz_file='%s_drz.fits' %(wfc3_grism.options['ROOT_DIRECT']),
                               is_grism=False)

    ### generate segmentation map for each grism exposure:
    wfc3_grism.showMessage('MAKING SEGEMENTAION MAPS FOR EACH DIRECT EXPOSURE')
    make_direct_exposure_segmaps(asn_direct_file, sigma=1.0)

    ### do the background subtraction:
    wfc3_grism.showMessage('DOING FULL DIRECT IMAGE SKY SUBTRACTION')
    asn = wfc3_grism.utils.ASNFile(asn_direct_file)
    for direct_exposure in asn.exposures:
        direct_image_sky_subtraction(flt='%s_flt.fits' %(direct_exposure),
                                     segmap='%s.seg.fits' %(direct_exposure),
                                     show=True)

    ### final drizzle to 1/2 pixel resolution with the background subtraction:
    wfc3_grism.showMessage('FINAL MULTIDRIZZLE WITH THE NEW SHIFTS + BACKGROUND SUBTRACTION')
    wfc3_grism.multidrizzle.multidrizzle_run(asn_direct_file, 
                                       shiftfile='%s_final_shifts.txt' %(wfc3_grism.options['ROOT_DIRECT']),
                                       pixfrac=wfc3_grism.options['PIXFRAC'], 
                                       final_scale=wfc3_grism.options['FINAL_DRIZZLE_PIXEL_SCALE'], 
                                       driz_cr=False,
                                       skysub=False,
                                       updatewcs=True,
                                       blot_back=False)


    ### now run sregister on CANDELS mosaic if using one as a detection image:
    if wfc3_grism.options['ADD_IMAGE_BORDER'] == True:
        wfc3_grism.showMessage('ADDING BORDER TO DRZ AND MAKING FINCAL CANDELS CUTOUT')
        add_border_to_drz(drz_image='%s_drz.fits' %(wfc3_grism.options['ROOT_DIRECT']), extension=500.)

    if wfc3_grism.options['DETECTION_IMAGE']:
        wfc3_grism.correct_shifts.run_sregister_for_detection_image(asn_direct_file, 
                                                                    sci_image=wfc3_grism.options['DETECTION_IMAGE'],
                                                                    wht_image=wfc3_grism.options['DETECTION_WHT'] )

def add_border_to_drz(drz_image, extension=300.):
    """
    Adds a border to the drizzled image so when registering the 
    CANDELS mosaics cuts out a larger section which is needed to can
    estimate contamination from objects outside the WFC3 Grism field of view.
    """
    
    ### open the fits file:
    drz_hdu = fits.open(drz_image, mode='update')
    
    ### loop through the two extensions:
    for ext in 'SCI', 'WHT':

        ### set the new x and y dimensions of the image
        newx = drz_hdu[ext].header['NAXIS1']*0.5+extension
        newy = drz_hdu[ext].header['NAXIS2']*0.5+extension
        
        ### adjust the centre pixel value:
        drz_hdu[ext].header.update('CRPIX1', newx)
        drz_hdu[ext].header.update('CRPIX2', newy)
        
        ### make the new image:
        largeim = np.zeros((drz_hdu[ext].header['NAXIS2']+2*extension, drz_hdu[ext].header['NAXIS1']+2*extension))
        largeim[int(extension):-int(extension),int(extension):-int(extension)]=drz_hdu[ext].data[:,:]
        
        ### assign to the fits data
        drz_hdu[ext].data = largeim

    ### write the new image to file:
    drz_hdu.flush()

def make_direct_exposure_segmaps(asn_direct, sigma=1.5):
    """
    make_segmap(root='ib3701ryq_flt', sigma=1)

    Get a segmentation image for a flt file after creating its 
    BLOT SCI and WHT images.

    DETECT_THRESH = ANALYSIS_THRESH = sigma
    """

    ### get an ASN object:
    asn = wfc3_grism.utils.ASNFile(asn_direct)

    ### set the default SExtractor parameters:
    se = wfc3_grism.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile(grism=False)

    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION, BACKGROUND'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_TYPE'] = 'NONE' ### no weight map for segementation map
    se.options['FILTER']    = 'Y'
    se.options['BACK_TYPE']     = 'AUTO'
    se.options['BACK_FILTERSIZE']     = '2'
    se.options['FILTER_NAME'] = 'default.conv'
    se.options['DETECT_THRESH']    = '%.1f' %sigma
    se.options['ANALYSIS_THRESH']  = '%.1f' %sigma
    se.options['MAG_ZEROPOINT'] = '%.2f' %wfc3_grism.options['MAG_ZEROPOINT']
    se.options['DETECT_MINAREA'] = '1'
    se.overwrite = True

    for exp in asn.exposures:

        se.options['CATALOG_NAME']    = '%s.BLOT.SCI.cat' %(exp)
        se.options['CHECKIMAGE_NAME'] = '%s.seg.fits, bg.fits' %(exp)
        se.options['WEIGHT_IMAGE']    = '%s.BLOT.WHT.fits' %(exp)
    
        status = se.sextractImage('%s.BLOT.SCI.fits' %(exp))

def direct_image_sky_subtraction(flt, segmap, sig_clip=3.0, show=True):
    """
    Subtract sky from direct exposure procedure is:

    - mask sources with segmentation map
    - get sky distribution from remaining pixels
    - iterativly sigma-clip this distribution to remove outliers
    - taken median of clipped sky distribution as the sky value across
      the image
    """
    
    ### open the fits files:
    flt_hdu = fits.open(flt)
    seg_hdu = fits.open(segmap)

    ### get the mask for extracting sky values:
    mask_sky = (seg_hdu[0].data == 0)

    ### make a copy of the image data to work on:
    im_data = np.copy(flt_hdu[1].data)

    if show:
        fig, ax = plt.subplots(figsize=(5.5, 5.5))
        ax.hist(im_data[mask_sky].flatten(), bins=np.arange(-2.0, 2.0, 0.02), 
                histtype='step', color='grey', lw=1.0)

    ### get a 1D array of the sky pixels:
    sky_pixels = np.copy(im_data[mask_sky].flatten())

    ### use sigma clipping to get the true sky values:
    true_sky = wfc3_grism.utils.sigma_clip(sky_pixels, sig=3, iters=None)
 
    ### find the median of this array:
    med_sky = np.median(true_sky.data[np.logical_not(true_sky.mask)])

    ### re-open the fits file, write new data to it:
    flt_hdu = fits.open(flt, mode='update')
    flt_hdu[1].data = im_data - med_sky

    flt_hdu.flush()

    if show:
        fig, ax = plt.subplots(figsize=(5.5, 5.5))
        final_data = im_data - med_sky
        ax.hist(final_data[mask_sky].flatten(), bins=np.arange(-2.0, 2.0, 0.02), 
                histtype='step', color='red', lw=1.0)
        fig.savefig('%s_background_historgram.pdf' %(flt.split('_flt.fits')[0]))




















