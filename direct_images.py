import figs

import pyfits
import os

def process_direct_images(asn_direct_file):

    #### first get the shifts between the individual direct exposures:
    print "\n Done!"
    figs.showMessage('RUNNING TWEAKSHIFT ON DIRECT IMAGES')
    figs.correct_shifts.run_tweakshifts_on_direct_exposures(asn_direct_file, verbose=True)

    # #### drizzle them together:
    print "\n Done!"
    figs.showMessage('RUNNING MULTIDRIZZLE ON DIRECT IMAGES')
    figs.multidrizzle.direct_multidrizzle_run(asn_direct_file, shiftfile='%s_shifts.txt' %(figs.options['ROOT_DIRECT']))

    ### cut out a region of the CANDEL image to use to align the grism exposures:
    print "\n Done!"
    figs.showMessage('CUTTING OUT CANDELS REGION TO ALIGN')
    figs.correct_shifts.run_sregister_to_cutout_CANDELS_region(asn_direct_file, mosiac_drz=figs.options['ALIGN_IMAGE'])

    # #### align dirzzled image to reference CANDELS mosaic:
    print "\n Done!"
    figs.showMessage('RUNNING TWEAKREG TO ALIGN DIRECT IMAGE TO CANDELS')
    figs.correct_shifts.align_direct_to_reference(asn_direct_file, verbose=True)

    #### copy over the new .flt files into the DATA directory to apply new shfits:
    print "\n Done!"
    figs.showMessage('RE-RUNNING MULTIDRIZZLE WITH THE NEW SHIFTS')
    copy_over_fresh_flt_files(asn_filename=figs.options['ASN_DIRECT'], from_path='../RAW')
    figs.multidrizzle.direct_multidrizzle_run(asn_direct_file, shiftfile='%s_reference_image_shifts.txt' %(figs.options['ROOT_DIRECT']))

    ### get the new cut-out of the CANDELS region, should now be aligned
    figs.correct_shifts.run_sregister_to_cutout_CANDELS_region(asn_direct_file, mosiac_drz=figs.options['ALIGN_IMAGE'])

def copy_over_fresh_flt_files(asn_filename, from_path='../RAW'):
	"""
	Copys over the raw _flt.fits files into the DATA directory, applies a best flat frame
	and updates to the most recent ICDTAB.
	"""

	# make an ASN object from the asn file
	ASN  = figs.utils.ASNFile(file='%s/%s' %(from_path, asn_filename))

	# extract the list of exposrues and loop through them:
	explist = []
	explist.extend(ASN.exposures)

	for exp in explist:

		# first find the file and open it as a pyfits object
		fits_file = figs.utils.find_fits_gz('%s/%s_flt.fits' %(from_path, exp))
		fi = pyfits.open(fits_file)

		# remove the current copy if one alrady exists:
		try:
			os.remove('./%s_flt.fits' %(exp))
		except:
			pass

		# write the fits file to the current ('/DATA') directory
		fi.writeto('./%s_flt.fits' %(exp), clobber=True)

		# now see that the flat-flied applied to image is the best available and apply
		# better flat if one is avaiable, can comment this out if just want to stick
		# with the original flat field.
		apply_best_flat('%s_flt.fits' %(exp), verbose=True)

def find_best_flat(flt_fits, verbose=True):
    """
    Find the most recent PFL file in $IREF for the filter used for the 
    provided FLT image.  Doesn't do any special check on USEAFTER date, just
    looks for the most-recently modified file. 
    """
    import glob
    import os.path
    import time
    
    IREF = os.environ["iref"]

    the_filter = pyfits.getheader(flt_fits,0).get('FILTER')
    
    pfls = glob.glob(IREF+'*pfl.fits')
    latest = 0
    best_pfl = None
    
    for pfl in pfls:
        head = pyfits.getheader(pfl)
        if head.get('FILTER') != the_filter:
            continue    
        
        this_created = os.path.getmtime(pfl)
        if this_created > latest:
            best_pfl = pfl
            latest = this_created
            
        if verbose:
            print '%s %s %s' %(pfl, the_filter, time.ctime(latest))
    
    return best_pfl #, the_filter, time.ctime(latest)

def apply_best_flat(fits_file, verbose=False):
    """
    Check that the flat used in the pipeline calibration is the 
    best available.  If not, multiply by the flat used and divide
    by the better flat.
    
    Input fits_file can either be an ASN list or an individual FLT file
    """

    fits_list = [fits_file]
    
    if fits_file.find('_asn.fits') > 0:
        asn = figs.utils.ASNFile(fits_file)
        fits_list = []
        for exp in asn.exposures:
            fits_list.append(exp+'_flt.fits')
    
    for file in fits_list:

        im = pyfits.open(file, 'update')

        USED_PFL = im[0].header['PFLTFILE'].split('$')[1]
        BEST_PFL = find_best_flat(file, verbose=False)

        if BEST_PFL is None:
        	BEST_PFL = USED_PFL

        IREF = os.environ["iref"]+"/"
            
        MSG = 'PFLAT, %s: Used= %s, Best= %s' %(file, USED_PFL, BEST_PFL)
        
        if BEST_PFL is None:
            figs.showMessage("No PFL file found! (NEED %s)" %(USED_PFL), warn=True)
                    
        BEST_PFL = os.path.basename(BEST_PFL)
                
        if USED_PFL != BEST_PFL:
            MSG += ' *'
            used = pyfits.open(IREF+USED_PFL)
            best = pyfits.open(IREF+BEST_PFL)
            
            im[1].data *= (used[1].data/best[1].data)[5:-5,5:-5]
            im[0].header.update('PFLTFILE', 'iref$'+BEST_PFL)
            im.flush()
            
        if verbose:
            print MSG
