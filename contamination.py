import wfc3_grism

import os
import glob

import numpy as np
import matplotlib.pyplot as plt

from pyraf import iraf

from astropy.io import fits

def register_fluxcube_images(fluxcube_filters):
    """
    Runs sregister assuming the images are already aligned with a larger mosiac

    e.g.:

    	ROOT_DIRECT = 'GOODS-S-34-F140W'
    	other_bands = [['/disk1/fc/research/HST/CANDELS/gsa_6epoch_wfc3_f125w_060mas_v0.2_drz.fits', 'F125W', 1248.6, 26.25],
                   	  ['/disk1/fc/research/HST/CANDELS/gsa_6epoch_wfc3_f160w_060mas_v0.2_drz.fits', 'F160W', 1536.9, 25.96],
                      ['/disk1/fc/research/HST/CANDELS/gs_acs_old_f850lp_060mas_v2_drz.fits', 'F850LP', 903.3, 24.84]]
    	DIRECT_MOSAIC = 'GOODS-S-34-F140W_drz.fits'

    """

    ### get the root name:
    root = wfc3_grism.options['ROOT_DIRECT']
    
    ### make the SExtractor files for the direct image and mosiac:
    se = wfc3_grism.sex.SExtractor()
    ### Set the output parameters required for image registration 
    ### (stored in [wfc3_grism source]/data/imreg.param) 
    se.imregParams()
    se.copyConvFile()
    se.overwrite = True

    for ib, filt in enumerate(fluxcube_filters):
              
        wfc3_grism.showMessage('REGISTERING: %s' %(filt[0]))
        
        #### First find its corresponding weight file:
        wht_image = glob.glob('%s/*%s*_wht*.fits' %(wfc3_grism.options['FLUXCUBE_FILTERS_DIR'], filt[1].lower()))[0]
        print 'WEIGHT IMAGE: %s' %(wht_image)
            
        ### iraf flpr()
        wfc3_grism.utils.iraf_flpr()

        #### sregister is being funny about trying to delete files IT creates ITSELF, so
        #### added a bit of a fudge to skip this error as it doesn't affect
        #### the creation of the registered mosiac image       
        iraf.unlearn('sregister')
        
        #### First register the _drz.fits mosiac file
        try:
            iraf.sregister(input=filt[0],
            			   reference='%s_drz.fits[SCI]' %(root),
                           output='%s_drz_reg.fits' %(filt[1]),
                           verbose=False)
        except:
            pass

        iraf.imcopy(input='%s_drz_reg.fits' %(filt[1]), 
                    output='%s_drz.fits[SCI,1,append]' %(filt[1]))
        
        ### iraf flpr()
        wfc3_grism.utils.iraf_flpr()

        #### sregister is being funny about trying to delete files IT creates ITSELF, so
        #### added a bit of a fudge to skip this error as it doesn't affect
        #### the creation of the registered mosiac image       
        iraf.unlearn('sregister')

        try:
            iraf.sregister(input=wht_image,
            			   reference='%s_drz.fits[SCI]' %(root),
                           output='%s_wht_reg.fits' %(filt[1]),
                           verbose=False)
        except:
            pass
        
        iraf.imcopy(input='%s_wht_reg.fits' %(filt[1]), 
                    output='%s_drz.fits[WHT,1,append]' %(filt[1]))     
           
        ### iraf flpr()            
        wfc3_grism.utils.iraf_flpr()
        
        ### plot histogram of object positions to check that registration has worked properly:
        ### set sextractor parameters: 
        se.options['DETECT_THRESH']    = '3.'
        se.options['ANALYSIS_THRESH']  = '3.'
        
        ### sort out the calalog for the registered image:
        se.options['CATALOG_NAME'] = 'reg.cat'
        status = se.sextractImage('%s_drz_reg.fits' %(filt[1]))
        reg_sexCat = wfc3_grism.sex.mySexCat('reg.cat')
        
        ### make catalog for f140w image, need to copy out the science extension to
        ### use with sextractor:
        se.options['CATALOG_NAME'] = 'direct.cat'

        if os.path.exists('./%s_SCI.fits' %(root)):
            pass
        else:
            iraf.imcopy(input='%s_drz.fits[SCI]' %(root), output='%s_SCI.fits' %(root))

        status = se.sextractImage('%s_SCI.fits' %(root))
        dir_sexCat = wfc3_grism.sex.mySexCat('direct.cat')
        
        separation = np.empty(len(dir_sexCat['X_IMAGE']))
        
        for i in range(len(dir_sexCat['X_IMAGE'])):
    
            ra_dir = dir_sexCat['ALPHA_J2000'][i]
            dec_dir = dir_sexCat['DELTA_J2000'][i]
    
            diff_ra = np.subtract(reg_sexCat['ALPHA_J2000'],ra_dir)
            diff_dec = np.subtract(reg_sexCat['DELTA_J2000'],dec_dir)
    
            sep = np.sqrt(np.add(np.square(diff_ra),np.square(diff_dec)))*3600.
        
            separation[i] = min(sep)
        
        print len(separation)
        fig, ax = plt.subplots(figsize=(5,5))
        ax.minorticks_on()
        ax.hist(separation, bins=np.arange(-5.0, 5.0, 0.05), histtype='step', color='k')
        ax.set_xlabel(r'$\mathrm{Separation}$ $/$ $\mathrm{arcsec}$', fontsize=14)
        ax.set_ylabel(r'$\mathrm{N}$', fontsize=14)
        fig.savefig('./%s_separation.pdf' %(filt[1]))
            
        ### these are the files IRAF can't find to remove... stupid IRAF
        tmps = glob.glob('tmps*')

        ### remove files that know are in directory:
        os.remove('./reg.cat')
        os.remove('./direct.cat')
        os.remove('./%s_drz_reg.fits' %(filt[1]))
        os.remove('./%s_wht_reg.fits' %(filt[1]))
        
        try:            
            for file in tmps:
                os.remove(file)
        except:
            pass

def make_zeropoint_file(fluxcube_filters):
    """
    Create the zeropoints.lis file needed to input to the fluxcube routine
    """

    ### open the file:
    outfile = open('zeropoints.lis','w')

    ### generate the lines:
    lines=[]
    for ib, band in enumerate(fluxcube_filters):
        lines.append('%s_drz.fits, %6.1f, %5.2f\n' %(band[1], band[2], band[3]))
   
    ### write lines to file and close:
    outfile.writelines(lines)
    outfile.close()

def setup_fluxcube():
    """
    Sets up all things needed for the aXe fluxcube task:
    -> registers the CANDELS mosaic images in each filter
    -> make the file with zeropoints etc ..
    """

    # get the fluxcube filters
    fluxcube_filters = wfc3_grism.options['FLUXCUBE_FILTERS']

    # register the images to the pointing:
    register_fluxcube_images(fluxcube_filters)

    # make the zeropoint file ('zeropoints.lis'):
    make_zeropoint_file(fluxcube_filters)

    # iraf.flpr:
    wfc3_grism.utils.iraf_flpr()
    wfc3_grism.utils.iraf_flpr()

    # delete descriptors in _seg files
    # or else you run into trouble with the blot call
    # in iraf.fcubeprep
    print 'Deleting descriptors in %s_seg.fits' %wfc3_grism.options['ROOT_DIRECT']
    wfc3_grism.utils.delete_multidrizzle_descriptions('%s_seg.fits' %wfc3_grism.options['ROOT_DIRECT'],
                                                       hdu_ext=0)

    for filt in fluxcube_filters:
        print 'Deleting descriptors in %s_drz.fits' %filt[1]
        wfc3_grism.utils.delete_multidrizzle_descriptions('%s_drz.fits' %filt[1], hdu_ext=1)

    # run fcube prep:
    iraf.fcubeprep(grism_image = '%s_drz.fits' %wfc3_grism.options['ROOT_GRISM'],
                   segm_image = '%s_seg.fits' %wfc3_grism.options['ROOT_DIRECT'],
                   filter_info = 'zeropoints.lis', 
                   AB_zero = True, 
                   dimension_info = wfc3_grism.options["AXE_EDGES"], 
                   interpol="poly5")


    