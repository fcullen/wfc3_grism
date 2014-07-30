import figs

from pyraf import iraf
from iraf import stsdas,dither,slitless,axe

import pyfits

import drizzlepac
from drizzlepac import astrodrizzle

import os

import numpy as np

INDEF = iraf.INDEF
no = iraf.no
yes = iraf.yes

class MultidrizzleRun():
    """
    MultidrizzleRun(root='IB3728050')
    
    Read a .run file output from MultiDrizzle.
    
    Get list of flt files and their shifts as used by multidrizzle.
    """
    def __init__(self, root='IB3728050'):
        import numpy as np
        
        runfile = root+'.run'
        self.root = root
        
        self.flt = []
        self.xsh = []
        self.ysh = []
        self.rot = []
        self.scl = 1.
        self.exptime = []
        
        for line in open(runfile,'r'):
            if line.startswith('drizzle.scale'):
                self.scl = line.split()[2]
            if line.startswith('drizzle '):
                spl = line.split()
                self.flt.append(spl[1].split('.fits')[0])
                self.exptime.append(-1)
                for tag in spl:
                    if tag.startswith('xsh'):
                        self.xsh.append(np.float(tag.split('=')[1]))
                    if tag.startswith('ysh'):
                        self.ysh.append(np.float(tag.split('=')[1]))
                    if tag.startswith('rot'):
                        self.rot.append(np.float(tag.split('=')[1]))
        
        self.count = len(self.flt)
        
    def blot_back(self, ii=0, SCI=True, WHT=True, copy_new=True, shape = None, ACS_CHIP=None):
        """
        blot_back(self, ii=0, SCI=True, WHT=True, copy_new=True)
    
        Blot the output DRZ file back to exposure #ii pixels.
        
        if SCI is True:
            blot science extension to FLT+'.BLOT.SCI.fits'

        if WHT is True:
            blot weight extension to FLT+'.BLOT.WHT.fits'
        
        if copy_new is True:
            imcopy SCI and WHT extensions of DRZ image to separate files.
        
        """
        #flt_orig = pyfits.open('../RAW/'+self.flt[ii]+'.fits.gz')
        figs.process_grism.flprMulti()
        
        if ACS_CHIP is None:
            ACS = False
            coeffs_ext = '_coeffs1.dat'
            EXT = 1
        else:
            ACS = True
            coeffs_ext = '_coeffs%0d.dat' %(3-ACS_CHIP)
            if ACS_CHIP == 1:
                EXT = 1
            else:
                EXT = 4
                
        ## Copy the coeffs1.dat file if necessary
        if not os.path.exists(self.flt[ii]+coeffs_ext):
            coeffs = figs.utils.get_package_data('wfc3_coeffs1.dat')
            fp = open(self.flt[ii]+'_coeffs1.dat','w')
            fp.writelines(coeffs)
            fp.close()
            
        if self.exptime[ii] < 0:
            try:
                flt_orig = pyfits.open(self.flt[ii]+'.fits')
                exptime = flt_orig[0].header.get('EXPTIME')
                if not ACS:
                    filter = flt_orig[0].header.get('FILTER').strip()
                else:
                    filter = (flt_orig[0].header.get('FILTER1').strip(),flt_orig[0].header.get('FILTER2').strip())
                
                flt_orig.close()
            except:
                exptime = 1.
                filter='INDEF'
        else:
            exptime = self.exptime[ii]
        
        if shape is None:  
            try:
                inNX = flt_orig[EXT].header.get('NAXIS1')
                inNY = flt_orig[EXT].header.get('NAXIS2')
            except:
                inNX = 1014
                inNY = 1014
            
            shape = (inNX, inNY)
        else:
            inNX, inNY = shape
        
        #### Need to update reference position of coeffs file
        #### for an output shape different than 1014, 1014
        coeffs = self.flt[ii]+coeffs_ext #'_coeffs1.dat'
        fp = open(coeffs)
        coeffs_lines = fp.readlines()
        fp.close()
        
        if shape != (1014, 1014):
            for i, line in enumerate(coeffs_lines):
                if line.strip().startswith('refpix'):
                    ### Default to center pixel
                    coeffs_lines[i] = 'refpix %9.3f %9.3f\n' %(inNX*1./2, inNY*1./2)
        
        fp = open('tmp'+coeffs_ext,'w')
        fp.writelines(coeffs_lines)
        fp.close()
        
                
        #iraf.delete(self.flt[ii]+'.BLOT.*.fits')
        files = glob.glob(self.flt[ii]+'.BLOT.*.fits')
        for file in files:
            os.remove(file)
            
        if copy_new:
            iraf.delete('drz_*.fits')
            # iraf.imcopy(self.root+'_drz.fits[1]','drz_sci.fits')
            # iraf.imcopy(self.root+'_drz.fits[2]','drz_wht.fits')
            
            ### NEED TO STRIP FITS HEADER
            im_drz = pyfits.open(self.root+'_drz.fits')
            sci = im_drz[1].data            
            s_hdu = pyfits.PrimaryHDU(sci)
            s_list = pyfits.HDUList([s_hdu])
            copy_keys = ['CTYPE1','CTYPE2','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2','LTM1_1','LTM2_2']
            s_list[0].header.update('EXPTIME',im_drz[0].header.get('EXPTIME'))
            s_list[0].header.update('CDELT1',im_drz[1].header.get('CD1_1'))
            s_list[0].header.update('CDELT2',im_drz[1].header.get('CD2_2'))
            for key in copy_keys:
                s_list[0].header.update(key, im_drz[1].header.get(key))
            s_list.writeto('drz_sci.fits', clobber=True)
            
            wht = im_drz[2].data
            w_hdu = pyfits.PrimaryHDU(wht)
            w_list = pyfits.HDUList([w_hdu])
            copy_keys = ['CTYPE1','CTYPE2','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2','LTM1_1','LTM2_2']
            w_list[0].header.update('EXPTIME',im_drz[0].header.get('EXPTIME'))
            w_list[0].header.update('CDELT1',im_drz[1].header.get('CD1_1'))
            w_list[0].header.update('CDELT2',im_drz[1].header.get('CD2_2'))
            for key in copy_keys:
                w_list[0].header.update(key, im_drz[1].header.get(key))
            w_list.writeto('drz_wht.fits', clobber=True)
            
        if SCI:
            iraf.blot(data='drz_sci.fits',
                outdata=self.flt[ii]+'.BLOT.SCI.fits', scale=self.scl,
                coeffs='tmp'+coeffs_ext, xsh=self.xsh[ii], 
                ysh=self.ysh[ii], 
                rot=self.rot[ii], outnx=inNX, outny=inNY, align='center', 
                shft_un='input', shft_fr='input', in_un='cps', out_un='cps', 
                interpol='poly5', sinscl='1.0', expout=exptime, 
                expkey='EXPTIME',fillval=0.0)
        
        if WHT:
            iraf.blot(data='drz_wht.fits',
                outdata=self.flt[ii]+'.BLOT.WHT.fits', scale=self.scl,
                coeffs='tmp'+coeffs_ext, xsh=self.xsh[ii], 
                ysh=self.ysh[ii], 
                rot=self.rot[ii], outnx=inNX, outny=inNY, align='center', 
                shft_un='input', shft_fr='input', in_un='cps', out_un='cps', 
                interpol='poly5', sinscl='1.0', expout=exptime, 
                expkey='EXPTIME',fillval=0.0)

def startMultidrizzle(root='ib3727050_asn.fits', use_shiftfile = True,
    skysub=True, updatewcs=True, driz_cr=True, median=True,
    final_scale=0.06, pixfrac=0.8, clean=True,
    final_outnx='', final_outny='', final_rot=0., ra='', dec='', 
    refimage='', unlearn=True, use_mdz_defaults=True, ivar_weights=True, build_drz=True):
    """
	startMultidrizzle(root='ib3727050_asn.fits', use_shiftfile = True,
                  skysub=True, final_scale=0.06, updatewcs=True, driz_cr=True,
                  median=True, final_scale=0.06, pixfrac=0.8, 
                  final_outnx='', final_outny='', final_rot=0., ra='', dec='',
                  refimage='', unlearn=True)
    
    Run multidrizzle on an input asn table.
    
    if `use_shiftfile` is True:
        use a root+'_shifts.txt' shiftfile.
    else: 
        no shiftfile
    
    if skysub is True:
        Run multidrizzle WITH sky subtraction
    else:
        "        "       WITHOUT   "         "
        
    final_scale: Final pixel scale of output image (arcsec)
    """

    asn_direct_file = root
    
    asn_direct = figs.utils.ASNFile(file=asn_direct_file)
    ROOT_DIRECT = asn_direct_file.split('_asn')[0]
    
    if use_shiftfile:
        shiftfile=ROOT_DIRECT+'_shifts.txt'
    else:
        shiftfile=''
    
    if skysub:
        skysub=yes
    else:
        skysub=no
    
    if updatewcs:
        updatewcs=yes
    else:
        updatewcs=no
    
    if driz_cr:
        driz_cr=yes
    else:
        driz_cr=no

    if median:
        median=yes
    else:
        median=no
    
    if unlearn:
        iraf.unlearn('multidrizzle')
    
    #### Set default parameters from pipeline mdz file
    if use_mdz_defaults:
        #### Read the first FLT image and read its MDRIZTAB file
        flt = pyfits.open(asn_direct.exposures[0]+'_flt.fits')
        
        #### Get the filter string, WFC3 or ACS
        if flt[0].header['INSTRUME'] == 'WFC3':
            filter=flt[0].header['FILTER']
            REF = 'iref'
        else:
            filter=(flt[0].header['FILTER1']+','+flt[0].header['FILTER2']).strip()
            REF = 'jref'
        
        mdz = pyfits.open(flt[0].header['MDRIZTAB'].replace(REF+'$',os.getenv(REF)+'/'))[1].data
        
        #### Force direct filter because parameters are a bit strange for grisms
        if filter.startswith('G1'):
            filter='F140W'
        
        if filter.startswith('G8'):
            filter='F814W'
        
        #### find 
        idx = np.where(mdz.field('filter') == filter)[0]
        if len(idx) == 0:
            filter='ANY'
            idx = np.where(mdz.field('filter') == filter)[0]
        
        #### Find right column for given "numimages" = len(exposures)  
        use = idx[0]
        for i in idx[1:]:
            if len(asn_direct.exposures) >= mdz.field('numimages')[i]:
                use = i
        
        #### Now set all of the parameters
        for param in mdz.names:
            try:
                value = mdz.field(param)[use]
                if (not np.isfinite(value)) | (value < -1000):
                    value = iraf.INDEF
                #
                iraf.dither.multidrizzle.setParam(param, value)
            except:
                #### some columns in the MDZ file aren't actually parameters, skip
                pass
        
        #### Don't use these default values from the mdrz file
        iraf.dither.multidrizzle.setParam('crbit','')
        iraf.dither.multidrizzle.setParam('combine_type','minmed')
        iraf.dither.multidrizzle.setParam('mdriztab',iraf.no)
        iraf.dither.multidrizzle.setParam('context',iraf.no)
        iraf.dither.multidrizzle.setParam('clean',iraf.no)
        iraf.dither.multidrizzle.setParam('ra','')
        iraf.dither.multidrizzle.setParam('dec','')
        iraf.dither.multidrizzle.setParam('dec','')
        iraf.dither.multidrizzle.setParam('runfile','')
        # iraf.dither.multidrizzle.setParam('driz_cr_snr','3.5 3.0')
        
    ### Set CR SNR parameter following candels
    if flt[0].header['INSTRUME'] == 'WFC3':
        iraf.dither.multidrizzle.driz_cr_snr = '6 3.0'
        iraf.dither.multidrizzle.driz_cr_scale = '1.6 0.7'
        
    if ivar_weights:
        #### Generate inverse variance weight map, will need flat + dark images
        #### in the iref or jref directories
        iraf.dither.multidrizzle.setParam('final_wht_type','IVM')
    
    if build_drz:
        build=iraf.yes
    else:
        build=iraf.no
          
    #### Run Multidrizzle
    # iraf.multidrizzle(input=asn_direct_file, \
    #    shiftfile=shiftfile, \
    #    output = '', skysub = skysub, updatewcs = updatewcs, driz_cr=driz_cr,
    #    final_scale = final_scale, final_pixfrac = pixfrac, median=median, 
    #    blot=median, driz_separate=median, static=median,
    #    driz_sep_outnx = final_outnx, driz_sep_outny = final_outny, 
    #    final_outnx=final_outnx, final_outny=final_outny, 
    #    final_rot=final_rot, ra=ra, dec=dec, refimage=refimage, build=build)
    
    astrodrizzle.AstroDrizzle(input=asn_direct_file, \
       shiftfile=shiftfile, \
       output = '', skysub = skysub, updatewcs = updatewcs, driz_cr=driz_cr,
       final_scale = final_scale, final_pixfrac = pixfrac, median=median, 
       blot=median, driz_separate=median, static=median,
       driz_sep_outnx = final_outnx, driz_sep_outny = final_outny, 
       final_outnx=final_outnx, final_outny=final_outny, 
       final_rot=final_rot, ra=ra, dec=dec, refimage=refimage, build=build)

    #### Delete created files    
    if clean is True:
        figs.process_grism.cleanMultidrizzleOutput()

