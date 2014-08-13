import os
import string
import time

import pyfits
import numpy as np

class ASNFile(object):
    """
    ASNFile()
            
    Class for handling ASN fits files.
        
    >>> asn = ASNFile(file='ib3701050_asn.fits', grow=False)
    >>> asn.exposures
    ['ib3701ryq', 'ib3701sbq', 'ib3701sdq', 'ib3701sqq']
    >>> asn.product
    'IB3701050'
    
    If grow=True, allow file rootnames to be 20 characters rather than 14.
    """

    def __init__(self, file=None, grow=True):
        self.file = file
        self.exposures = []
        self.product = None
        if file:
            self._read_asn_file(grow=grow)

    def _read_asn_file(self, grow=True):
        """
        _read_asn_file(self)
        
        Read an ASN FITS file (self.file).
        """
        import numpy as np
        from warnings import warn
        
        self.in_fits = pyfits.open(self.file)
        
        data = self.in_fits[1].data
        
        if grow:
            #### Allow more characters in the MEMNAME column
            memname = pyfits.Column(name='MEMNAME', format='40A', array=self.in_fits[1].columns[0].array.astype('S40'), disp='A40')
            memtype = self.in_fits[1].columns[1]
            memprsnt = self.in_fits[1].columns[2]
            coldefs = pyfits.ColDefs([memname, memtype, memprsnt])
            hdu = pyfits.new_table(coldefs)
            hdu.header = self.in_fits[1].header
            hdu.header['TFORM1'] = '40A'
            hdu.header['TDISP1'] = 'A40'
            hdu.header['NAXIS1'] += 26
            self.in_fits[1] = hdu        
    
        data = self.in_fits[1].data
        
        self.header = self.in_fits[0].header
        
        names = data.field('MEMNAME')
        types = data.field('MEMTYPE')
        
        ##### Exposures
        exp_idx  = np.where(types == 'EXP-DTH')
        if exp_idx[0].shape[0] == 0:
            warn('ASN file %s has no EXP-DTH items')
        else:
            self.exposures = []
            for exp in names[exp_idx]:
                self.exposures.append(exp.lower())
        
        ##### Products
        prod_idx = np.where(types == 'PROD-DTH')
        if prod_idx[0].shape[0] != 1:
            warn('ASN file %s has N != 1 PROD-DTH items' %self.file )
            self.product = None
        else:
            self.product = names[prod_idx[0]][0].upper()
    

    def write(self, out_file=None, clobber=True):
        """
        write(self,out_file=None, clobber=True)

        out_file='self' writes to `self.file`.
        """
        if not out_file:
            print "USAGE:: asn.write(out_file='output_asn.fits')"
        else:
            if out_file == 'self':
                out_file = self.file
            
            nexp  = self.exposures.__len__()
            if self.product:
                nprod = 1
            else:
                nprod = 0
            nrows = nexp + nprod
            #### Primary HDU
            hdu = self.in_fits[0].copy()
            #### BinTable HDU
            tbhdu = pyfits.new_table(self.in_fits[1].columns, nrows=nrows, fill=True)
            for i in range(nexp):
                tbhdu.data[i] = (self.exposures[i].upper(), 'EXP-DTH', True)
            if nprod > 0:
                tbhdu.data[i+1] = (self.product, 'PROD-DTH', True)
            
            tbhdu.header = self.in_fits[1].header.copy()
            tbhdu.header.update('ASN_ID',out_file.split('_asn.fits')[0])
            tbhdu.header.update('ASN_TAB',out_file)
            #### Create HDUList and write it to output file
            self.out_fits = pyfits.HDUList([hdu,tbhdu])
            if 'EXTEND' not in hdu.header.keys():
                hdu.header.update('EXTEND', True, after='NAXIS')
                
            self.out_fits.writeto(out_file, clobber=clobber)
    
    def showContents(self):
        """
        showContents()
        
        >>> x = ASNFile(file='ib3702060_asn.fits')
        >>> x.showContents()
        1   ib3703uxq    EXP-DTH      yes
        2   ib3703vaq    EXP-DTH      yes
        3   ib3703vcq    EXP-DTH      yes
        4   ib3703vsq    EXP-DTH      yes
        5   IB3703050    PROD-DTH     yes
        """
        if self.exposures.__len__() > 0:
            for i,exp in enumerate(self.exposures):
                print '%5d   %s    EXP-DTH      yes' %(i+1,exp)
            print '%5d   %s    PROD-DTH     yes' %(i+2,self.product)
    
    def append(self, new):
        """
        append(self, new)
        
        `new` must be an instance of ASNFile.
            
        `new.exposures` are added to the `self.exposures` list.
        """

        from warnings import warn
        if not isinstance(new,self.__class__):
            warn("argument is not an instance of ASNFile()")
        else:
            self.exposures.extend(new.exposures)

def find_fits_gz(fits_file, hard_break = True):
    """
    foo = find_fits_gz(fits_file, hard_break = True)
    
    With ``fits_file`` being some filename with an extension ``.fits``
    (ib3713wvq_flt.fits), check to see if the file itself or its gzipped
    counterpart (fits_file+'.gz') exists.
    
    If neither is found, 
        hard_break = True  : raise an IOError
        hard_break = False : return None
    
    """
    import os
    if os.path.exists(fits_file):
        return fits_file
    if os.path.exists(fits_file+'.gz'):
        return fits_file+'.gz'
    #### File not found.  Either raise an error or return None
    if hard_break:
        raise IOError('File %s[.gz] not found in %s' 
                              %(fits_file, os.getcwd()))
    else:
        return None

def make_targname_asn(asn_file, newfile=True, use_filtname=True, path_to_flt='../RAW/'):
    """
    Take an ASN file like 'ibhm51030_asn.fits' and turn it into 
    'COSMOS-3-F140W_asn.fits'
    """
    asn = ASNFile(path_to_flt+asn_file)
    
    flt_file = find_fits_gz(path_to_flt+'/'+asn.exposures[0]+'_flt.fits')
    im = pyfits.open(flt_file)
    
    instrum = im[0].header['INSTRUME']
    if instrum == 'ACS':
        filter=im[0].header['FILTER1']
    else:
        filter=im[0].header['FILTER']
        
    if filter.startswith('F'):
        type='D'
    else:
        type='G'
    
    if use_filtname:
        type=filter
        
    target = im[0].header['TARGNAME']
    target = target.replace('SOUTH','S')
    target = target.replace('GNGRISM','GOODS-N-')
    target = target.replace('GEORGE','GEORGE-')
    
    if target == 'MARSHALL':
        #### Add the pointing number, 1-6
        ID = asn.exposures[0][5]
        target+='-'+ID

    if target == 'PRIMO':
        #### Add the date, like "1026"
        date = ''.join(im[0].header['DATE-OBS'].split('-')[1:])
        target+='-'+date
    
    if target.startswith('GEORGE'):
        #### Add the date, like "1026"
        hour = np.int(im[0].header['TIME-OBS'].split(':')[0])
        if hour > 12:
            target='GEORGE-2'
        else:
            target='GEORGE-1'
    # 
    # if target.startswith('GOODS-N'):
    #     #### Some Visits were redone
    #     date = im[0].header['DATE-OBS']
    #     if date > '2011-04-01':
    #         target = target.replace('GOODS-N','GOODS-N2')
    
    product = target+'-'+type
    asn.product = product
    if newfile:
        asn.write(product+'_asn.fits', clobber=True)
    return product+'_asn.fits'

def get_package_data(dataname):
    """
    (taken from astropysics.io)
    Use this function to load data files distributed with astropysics in the 
    astropysics/data directory
    
    `dataname` is the file name of a file in the *threedhst/data* directory, and
    a string with the contents of the file will be returned
    """
    try:
        ### Find the data directory in the root
        ### directory of the threedhst package
        from . import __name__ as rootname
        from . import __file__ as rootfile
        from pkgutil import get_loader
        from os.path import dirname
        path = dirname(rootfile)+'/data/'+dataname
        return get_loader(rootname).get_data(path)
    except:
        ### Hardwired  in case relative import doesn't work
        fp = open('/research/HST/GRISM/3DHST/progs/threedhst/data/'+dataname)
        output = fp.read()
        fp.close()
        return output