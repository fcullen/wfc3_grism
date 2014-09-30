import wfc3_grism

from pyraf import iraf

import os
import glob

import numpy as np

from astropy.io import fits

def multidrizzle_run(asn_file, shiftfile, pixfrac=0.8, final_scale=0.06, driz_cr=False, skysub=True):
    """
    Performs an multidrizzle run on a set of input direct images.

    Some CANDELS defaults:

    driz_cr_snr='6 3.0',
    driz_cr_scale='1.6 0.7'
    """

    ### iraf flpr()
    wfc3_grism.utils.iraf_flpr()

    ### unlearn the routine:
    iraf.unlearn('multidrizzle')

    iraf.multidrizzle(input=asn_file,
                      shiftfile=shiftfile,
                      output ='', 
                      skysub = skysub, 
                      updatewcs =True,
                      static=True,
                      static_sig=4.0,
                      driz_separate=True, 
                      driz_sep_kernel='turbo',
                      median=True, 
                      blot=True, 
                      driz_cr=driz_cr,
                      driz_cr_snr='6 3.0',
                      driz_cr_scale='1.6 0.7',
                      final_scale = final_scale, 
                      final_pixfrac = pixfrac,
                      driz_combine=True,
                      final_wht_type='IVM',
                      clean=False)

    ### clean up:
    clean_multidrizzle_output()

def blot_run(asn_file, drz_file, is_grism=True):

    ### iraf flpr()
    wfc3_grism.utils.iraf_flpr()

    ### unlearn the routine:
    iraf.unlearn('blot')

    ### get an ASN object:
    asn = wfc3_grism.utils.ASNFile(asn_file)

    ### get the multidrizzle run file:
    if is_grism:
        runfile = '%s.run' %(wfc3_grism.options['ROOT_GRISM'])
    else:
        runfile = '%s.run' %(wfc3_grism.options['ROOT_DIRECT'])

    ### find the shifts for each exposure which you need as input to blot:
    xsh = []
    ysh = []
    rot = []

    for line in open(runfile, 'r'):
        if line.startswith('drizzle'):
            params = line.split()
            for param in params:
                if param.startswith('xsh'):
                    xsh.append(np.float(param.split('=')[1]))
                if param.startswith('ysh'):
                    ysh.append(np.float(param.split('=')[1]))
                if param.startswith('rot'):
                    rot.append(np.float(param.split('=')[1]))

    ### blot is seg faulting when using the _drz.fits files so
    ### try copy out the ['SCI'] extension and use that instead:
    drz_hdu = fits.open(drz_file)
    sci_data = drz_hdu['SCI'].data
    sci_header = drz_hdu['SCI'].header

    ### need to include exposure time in header for blot:
    exptime = drz_hdu[0].header['EXPTIME']
    sci_header.set('EXPTIME', exptime)

    new_hdu = fits.PrimaryHDU(data=sci_data, header=sci_header)
    new_hdu.writeto('%s_drz_sci.fits' %(wfc3_grism.options['ROOT_GRISM']), clobber=True)

    blot_infile = '%s_drz_sci.fits' %(wfc3_grism.options['ROOT_GRISM'])

    ### loop through each exposure and blot back:
    for i, exp in enumerate(asn.exposures):

        iraf.blot(data=blot_infile,
                outdata='%s.BLOT.SCI.fits' %(exp), 
                scale=1.0,
                coeffs='%s_flt_coeffs1.dat' %(exp), 
                xsh=xsh[i], 
                ysh=ysh[i], 
                rot=rot[i], 
                outnx=1014, 
                outny=1014, 
                align='center', 
                shft_un='input', 
                shft_fr='input', 
                in_un='cps', 
                out_un='cps', 
                interpol='poly5', 
                sinscl='1.0', 
                expout=exptime,
                expkey='EXPTIME',
                fillval=0.0)

    ### cleanup:
    rmfiles = [blot_infile]

    for rfile in rmfiles:
        try:
            os.remove(rfile)
        except:
            pass

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

def make_drizzled_contamination_image(asn_grism_file):
    """
    Drizzle the output contamination images to compare with the drizzled real image
    """

    ### get starting directory, if not in data directory, go there!
    cwd = os.get_cwd()
    if cwd != wfc3_grism.options['ROOT_DIR']:
        os.chdir('%s/DATA' %wfc3_grism.optionsp['ROOT_DIR'])

    ### get the asn girsm object:
    asn_grism = wfc3_grism.utils.ASNFile(asn_grism_file)

    for exp in asn_grism.exposures:
        ### copy OUTPUT_G141/*_CONT.fits into existing grism FLT files in DATA directory
        flt = fits.open('%s_flt.fits' %exp,'update')
        cont = fits.open('../OUTPUT_%s/%s_flt_2.CONT.fits' %(wfc3_grism.options['GRISM_NAME'],
                                                             exp))
        flt[1].data  = cont[1].data
        
        flt.flush()

    ### make a separate CONT asn file:
    asn_grism.product = '%s_CONT' %(wfc3_grism.options['ROOT_GRISM'])
    asn_grism.write('%s_CONT_asn.fits'%(wfc3_grism.options['ROOT_GRISM']))

    ### open it as the asn_file for Multidrizzle:
    asn_file = wfc3_grism.utils.ASNFile('%s_CONT_asn.fits' %(wfc3_grism.options['ROOT_GRISM']))

    ### iraf flpr()
    wfc3_grism.utils.iraf_flpr()

    ### unlearn the routine:
    iraf.unlearn('multidrizzle')

    ### first apply x_shift
    iraf.multidrizzle(input=asn_file,
                      shiftfile='%s_final_shifts_xshift' %(wfc3_grism.options['ROOT_GRISM'),
                      output ='', 
                      skysub = False, 
                      updatewcs =True,
                      static=True,
                      static_sig=4.0,
                      driz_separate=True, 
                      driz_sep_kernel='turbo',
                      median=True, 
                      blot=True, 
                      driz_cr=False,
                      driz_cr_snr='6 3.0',
                      driz_cr_scale='1.6 0.7',
                      final_scale = wfc3_grism.options['DRZSCALE'], 
                      final_pixfrac = wfc3_grism.options['PIXFRAC'],
                      driz_combine=True,
                      final_wht_type='IVM',
                      clean=False)

    ### iraf flpr()
    wfc3_grism.utils.iraf_flpr()

    ### unlearn the routine:
    iraf.unlearn('multidrizzle')

    ### now apply y_shift
    iraf.multidrizzle(input=asn_file,
                      shiftfile='%s_final_shifts_yshift' %(wfc3_grism.options['ROOT_GRISM'),
                      output ='', 
                      skysub = False, 
                      updatewcs =True,
                      static=True,
                      static_sig=4.0,
                      driz_separate=True, 
                      driz_sep_kernel='turbo',
                      median=True, 
                      blot=True, 
                      driz_cr=False,
                      driz_cr_snr='6 3.0',
                      driz_cr_scale='1.6 0.7',
                      final_scale = wfc3_grism.options['DRZSCALE'], 
                      final_pixfrac = wfc3_grism.options['PIXFRAC'],
                      driz_combine=True,
                      final_wht_type='IVM',
                      clean=False)

    ### clean up:
    clean_multidrizzle_output()

    ### copy over the raw grism flt files which were overwritten:
     for exp in asn_grism.exposures:
        shutil.copy('../RAW/%s_flt.fits' %exp, './')

    ### change back to current directory:
    os.chdir(cwd)




