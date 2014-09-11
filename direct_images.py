import figs
import os

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
    figs.correct_shifts.run_sregister_to_cutout_CANDELS_region(asn_direct_file, mosiac_drz=figs.options['ALIGN_IMAGE'])

    ### align dirzzled image to reference CANDELS mosaic:
    figs.showMessage('ALIGNING DIRECT IMAGE TO CANDELS')
    figs.correct_shifts.align_direct_to_reference(asn_direct_file, verbose=True)

    #### copy over the new .flt files into the DATA directory to apply new shfits:
    figs.showMessage('RE-RUNNING MULTIDRIZZLE WITH THE NEW SHIFTS')
    figs.utils.copy_over_fresh_flt_files(asn_filename=figs.options['ASN_DIRECT'], from_path='../RAW')
    figs.multidrizzle.multidrizzle_run(asn_direct_file, 
                                       shiftfile='%s_final_shifts.txt' %(figs.options['ROOT_DIRECT']),
                                       pixfrac=0.8, 
                                       final_scale=0.06, 
                                       driz_cr=False,
                                       skysub=True)

    ### at this point can perform background subtraction etc. on the direct image but
    ### don't need to do this if using different detection image, it is enough that the
    ### correct WCS solution has been found to align the image to CANDELS