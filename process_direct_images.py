import figs

import pyfits
import os

def copy_over_fresh_flt_files(asn_filename, from_path='../RAW'):
	"""
	Copys over the raw _flt.fits files into the DATA directory, applies a best flat frame
	and updates to the most recent ICDTAB.
	"""

	# make an ASN object from the asn file
	ASN  = figs.utils.ASNFile(file=asn_filename)

	# extract the list of exposrues and loop through them:
	explist = []
	explist.extend(ASN.exposures)

	print explist

	for exp in explist:

		print exp

		# first find the file and open it as a pyfits object
		fits_file = figs.utils.find_fits_gz(from_path+'/'+exp+'_flt.fits')
		fi = pyfits.open(fits_file)

		# remove the current copy if one alrady exists:
		try:
			os.remove('./'+exp+'_flt.fits')
		except:
			pass

		# write the fits file to the current ('/DATA') directory
		fi.writeto('./'+exp+'_flt.fits', clobber=True)

		# now see that the flat-flied applied to image is the best available and apply
		# better flat if one is avaiable, can comment this out if just want to stick
		# with the original flat field.
		apply_best_flat(exp+'_flt.fits', verbose=True)

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

def align_raw_flt_to_reference(raw_flt, reference_image):

	# get out the extension of the flt file:
	root_name = raw_flt.split('_flt.fits')[0]

	#### Use swarp to combine the alignment images to the same image 
	#### dimensions as the direct mosaic
	try:
		os.remove(root_name+'_align.fits')
	except:
		pass

	swarp_to_image(input_image=raw_flt,
				   ref_image=reference_image,
				   output=root_name+'_align.fits', image_extension=1,
				   ref_extension=0)

def swarp_to_image(input_image=None, ref_image=None,output=None,
                   image_extension=0, ref_extension=0):
    """
	swarp_to_image(input=None,matchImage=None,output=None,
                 match_extension=0)
    
    SWarp input image to same size/scale as matchIMage.
    Output default is input_root+'.match.fits'
    
    Input is a list of images
    """
    from figs.sex import SWarp

    if not ref_image:
        return False
    
    if not output:
        output = os.path.basename(input).split('.fits')[0]+'.match.fits'
    
    #### initialize SWarp
    sw = SWarp()
    sw._aXeDefaults()
    sw.overwrite = True
    
    #### Get first guess coordinates
    sw.swarpMatchImage(input_image, extension=image_extension)
    status = sw.swarpImage(input_image + '[%d]' %image_extension, mode='wait')
    os.remove('coadd.fits')
    os.remove('coadd.weight.fits')
    
    #### Recenter
    sw.swarpRecenter()
    
    #### Make the final output image
    sw.options['IMAGEOUT_NAME'] = output
    #sw.options['WEIGHTOUT_NAME'] = base+'.match.weight.fits'
    
    ref_image += '[%d]' %ref_extension
    status = sw.swarpImage(ref_image, mode='direct')
    os.remove('coadd.weight.fits')