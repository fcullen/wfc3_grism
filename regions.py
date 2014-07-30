import os
import string
import time

import pyfits
import numpy as np

import figs

def asn_region(asn_file, path_to_flt='./'):
    """
    asn_region(asn_file)
    
    Create a DS9 region file for the exposures defined in an ASN file.
    
    """
    ##### Output file
    output_file = asn_file.split('.fits')[0]+'.pointing.reg'
    fp = open(output_file,'w')
    fp.write('fk5\n') ### WCS coordinates
    
    ##### Read ASN file
    asn = figs.utils.ASNFile(asn_file)
    NEXP = len(asn.exposures)
    RAcenters  = np.zeros(NEXP)
    DECcenters = np.zeros(NEXP)

    ##### Loop through exposures and get footprints
    for i, exp_root in enumerate(asn.exposures):
        flt_file =figs.utils.find_fits_gz(path_to_flt + '/' + exp_root.lower()+'_flt.fits', hard_break = True)
        
        #head = pyfits.getheader(exp_root.lower()+'_flt.fits')
        head = pyfits.getheader(flt_file)
        if head.get('INSTRUME') == 'ACS':
            extensions=[1,4]
        else:
            extensions=[1]
        
        for ext in extensions:
            regX, regY = wcs_polygon(flt_file,extension=ext)
            line = "polygon(%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f)"  %(regX[0],regY[0],regX[1],regY[1],regX[2],regY[2],regX[3],regY[3])

            RAcenters[i] += np.mean(regX)/len(extensions)
            DECcenters[i] += np.mean(regY)/len(extensions)
            
            fp.write(line+' # color=magenta\n')
        
    ##### Text label with ASN filename
    fp.write('# text(%10.6f,%10.6f) text={%s} font="Helvetica 14 normal" color=magenta\n' \
        %(np.mean(RAcenters),np.mean(DECcenters),
          asn_file.split('_asn.fits')[0]))
    fp.close()

    figs.showMessage('Create region file, %s.' %output_file)

def wcs_polygon(fits_file, extension=1, use_pywcs=False):
    """    
    X, Y = wcs_polygon(fits_file, extension=1)
    
    Calculate a DS9/region polygon from WCS header keywords.  
        
    Will try to use pywcs.WCS.calcFootprint if pywcs is installed.  Otherwise
    will compute from header directly.
    
    """
    ##### Open the FITS file
    hdulist = pyfits.open(fits_file) 
    ##### Get the header
    try:
        sci = hdulist[extension].header
    except IndexError:
        print 'ERROR 3D-HST/wcs_polygon:\n'+\
              'Extension #%d out of range in %s' %(extension, fits_file)
        raise
    
    #### Try to use pywcs if it is installed
    pywcs_exists = True
    try:
        import pywcs
    except:
        pywcs_exists = False   
    
    if pywcs_exists & use_pywcs:
        wcs = pywcs.WCS(sci)
        footprint = wcs.calcFootprint()
        regX = footprint[:,0]    
        regY = footprint[:,1]    
        return regX, regY
    
    #### Do it by hand if no pywcs    
    NAXIS = [sci['NAXIS1'],sci['NAXIS2']]
    CRPIX = [sci['CRPIX1'],sci['CRPIX2']]
    CRVAL = [sci['CRVAL1'],sci['CRVAL2']]
    
    CD1_1 = sci['CD1_1']
    CD2_2 = sci['CD2_2']
    try:
        CD1_2 = sci['CD1_2']
        CD2_1 = sci['CD2_1']
    except:
        CD1_2 = 0.
        CD2_1 = 0.
        
    cosDec = np.cos(CRVAL[1]/180*np.pi)
    ##### Make region polygon from WCS keywords
    regX = CRVAL[0] + \
            ((np.array([0,NAXIS[0],NAXIS[0],0])-CRPIX[0])*CD1_1 +                        
             (np.array([0,0,NAXIS[1],NAXIS[1]])-CRPIX[1])*CD1_2) / cosDec
    
    regY = CRVAL[1] + \
            ((np.array([0,NAXIS[0],NAXIS[0],0])-CRPIX[0])*CD2_1 +         
             (np.array([0,0,NAXIS[1],NAXIS[1]])-CRPIX[1])*CD2_2)
             
    return regX, regY