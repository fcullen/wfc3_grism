import figs

import os
import glob

import numpy as np
import matplotlib.pyplot as plt

from pyraf import iraf

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
    root = figs.options['ROOT_DIRECT']
    
    ### make the SExtractor files for the direct image and mosiac:
    se = figs.sex.SExtractor()
    ### Set the output parameters required for image registration 
    ### (stored in [figs source]/data/imreg.param) 
    se.imregParams()
    se.copyConvFile()
    se.overwrite = True

    for ib, filt in enumerate(fluxcube_filters):
              
        figs.showMessage('REGISTERING: %s' %(filt[0]))
        
        #### First find its corresponding weight file:
        wht_image = glob.glob('%s/*%s*_wht.fits' %(figs.options['FLUXCUBE_FILTERS_DIR'], filt[1].lower()))[0]
        print 'WEIGHT IMAGE: %s' %(wht_image)
            
        ### iraf flpr()
        figs.utils.iraf_flpr()

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
        figs.utils.iraf_flpr()

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
        figs.utils.iraf_flpr()
        
        ### plot histogram of object positions to check that registration has worked properly:
        ### set sextractor parameters: 
        se.options['DETECT_THRESH']    = '5.'
        se.options['ANALYSIS_THRESH']  = '5.'
        
        ### sort out the calalog for the registered image:
        se.options['CATALOG_NAME'] = 'reg.cat'
        status = se.sextractImage('%s_drz_reg.fits' %(filt[1]))
        reg_sexCat = figs.sex.mySexCat('reg.cat')
        
        ### make catalog for f140w image, need to copy out the science extension to
        ### use with sextractor:
        se.options['CATALOG_NAME'] = 'direct.cat'

        if os.path.exists('./%s_SCI.fits' %(root)):
            pass
        else:
            iraf.imcopy(input='%s_drz.fits[SCI]' %(root), output='%s_SCI.fits' %(root))

        status = se.sextractImage('%s_SCI.fits' %(root))
        dir_sexCat = figs.sex.mySexCat('direct.cat')
        
        separation = np.empty(len(dir_sexCat['X_IMAGE']))
        
        for i in range(len(dir_sexCat['X_IMAGE'])):
    
            ra_dir = dir_sexCat['ALPHA_J2000'][i]
            dec_dir = dir_sexCat['DELTA_J2000'][i]
    
            diff_ra = np.subtract(reg_sexCat['ALPHA_J2000'],ra_dir)
            diff_dec = np.subtract(reg_sexCat['DELTA_J2000'],dec_dir)
    
            sep = np.sqrt(np.add(np.square(diff_ra),np.square(diff_dec)))*3600.
        
            separation[i] = min(sep)
        
        fig, ax = plt.subplots(figsize=(5,5))
        ax.minorticks_on()
        ax.hist(separation, bins=np.arange(-1.0, 1.0, 0.05), histtype='step', color='k')
        ax.set_xlabel(r'$\mathrm{Separation}$ $/$ $\mathrm{arcsec}$', fontsize=14)
        ax.set_ylabel(r'$\mathrm{N}$', fontsize=14)
        fig.savefig('./%s_separation.pdf' %(filt[1]))
            
        #### These are the files IRAF can't find to remove... stupid IRAF
        tmps = glob.glob('tmps*')
        
        try:
            os.remove('./temp.db')
            os.remove('./direct_tmp.fits')
            os.remove('./mosaic_imreg.cat')
            os.remove('./direct_imreg.cat')
            os.remove('./%s_drz_reg.fits' %(filt[1]))
            os.remove('./%s_wht_reg.fits' %(filt[1]))
            os.remove('./%s_SCI.fits' %(root))
            os.remove('./sregister.db')
            os.remove('./reg.cat')
            os.remove('./direct.cat')
            for file in tmps:
                os.remove(file)
        except:
            pass

def setup_fluxcube():
    """
    Sets up all things needed for the aXe fluxcube task:
    -> registers the CANDELS mosaic images in each filter
    -> make the file with zeropoints etc ..
    """

    ### get the fluxcube bans
    fluxcube_filters = figs.options['FLUXCUBE_FILTERS']
    register_fluxcube_images(fluxcube_filters)

    