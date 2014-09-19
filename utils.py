import figs

import os
import string
import time

from astropy.io import fits
import numpy as np

from pyraf import iraf

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
        from warnings import warn
        
        self.in_fits = fits.open(self.file)
        
        data = self.in_fits[1].data
        
        if grow:
            #### Allow more characters in the MEMNAME column
            memname = fits.Column(name='MEMNAME', format='40A', array=self.in_fits[1].columns[0].array.astype('S40'), disp='A40')
            memtype = self.in_fits[1].columns[1]
            memprsnt = self.in_fits[1].columns[2]
            coldefs = fits.ColDefs([memname, memtype, memprsnt])
            hdu = fits.new_table(coldefs)
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
            tbhdu = fits.new_table(self.in_fits[1].columns, nrows=nrows, fill=True)
            for i in range(nexp):
                tbhdu.data[i] = (self.exposures[i].upper(), 'EXP-DTH', True)
            if nprod > 0:
                tbhdu.data[i+1] = (self.product, 'PROD-DTH', True)
            
            tbhdu.header = self.in_fits[1].header.copy()
            tbhdu.header.update('ASN_ID',out_file.split('_asn.fits')[0])
            tbhdu.header.update('ASN_TAB',out_file)
            #### Create HDUList and write it to output file
            self.out_fits = fits.HDUList([hdu,tbhdu])
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

class ConfFile(object):
    """
    Conf(infile='WFC3.IR.G141.V1.0.conf')
    
    Read an aXe configuration file for easy manipulation of the parameters.
    """            
    def _getPath(self):
        """
        _getPath()
        
        Figure out if we're in the root directory or in DATA
        """

        self.path = figs.options['REDUCTION_CONFIG_FILE_DIRECTORY']
        
    def _processLines(self):
        """
        _processLines()
        
        Read the lines, extracting PARAMETER VALUE pairs
        """
        self.nlines = len(self.lines)
        self.params = {}
        self._pline = {}
        for i,line in enumerate(self.lines):
            if (line[0] is not '#') & (line.strip() is not  ''):
                spl = line.split()
                self.params[spl[0]] = ' '.join(spl[1:])
                self._pline[spl[0]] = i
        self.nkeys = self.params.keys().__len__()
        
    def _readLines(self):
        """
        _readlines()
        """
        #self._getPath()
        fp = open('%s/%s' %(self.path, self.infile),'r')
        self.lines = fp.readlines()
        fp.close()
        
    def __init__(self, infile='WFC3.IR.G141.V1.0.conf', path=None):
        self.infile = infile
        if path is None:
            self._getPath()
        else:
            if not path.startswith('/'):
                path = os.getcwd()+'/'+path
            self.path = path
        self._readLines()
        self._processLines()
        
    def _assignPars(self):
        """
        _assignPars()
        
        Apply parameter values to self.lines
        """
        import numpy as np
        for key in self.params.keys():
            param = self.params[key]
            if type(param) is not type('xxx'):
                param = np.str(param)
            
            #### New parameter?
            if self._pline.has_key(key):
                self.lines[self._pline[key]] = key + ' ' + param +'\n'
            else:
                self.lines.append(key + ' ' + param +'\n')
        
        self._processLines()
    
    def writeto(self, output='tmp.conf'):
        """
        writeto(output='tmp.conf')
        """ 
        #self._getPath()
        self._assignPars()
        fp = open('%s/%s' %(self.path, output),'w')
        fp.writelines(self.lines)
        fp.close()

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
    im = fits.open(flt_file)
    
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

def iraf_flpr():

    iraf.flpr()
    iraf.flpr()
    iraf.flpr()

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

        # first find the file and open it as a fits object
        fits_file = figs.utils.find_fits_gz('%s/%s_flt.fits' %(from_path, exp))
        fi = fits.open(fits_file)

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

    the_filter = fits.getheader(flt_fits,0).get('FILTER')
    
    pfls = glob.glob(IREF+'*pfl.fits')
    latest = 0
    best_pfl = None
    
    for pfl in pfls:
        head = fits.getheader(pfl)
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

        im = fits.open(file, 'update')

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
            used = fits.open(IREF+USED_PFL)
            best = fits.open(IREF+BEST_PFL)
            
            im[1].data *= (used[1].data/best[1].data)[5:-5,5:-5]
            im[0].header.update('PFLTFILE', 'iref$'+BEST_PFL)
            im.flush()
            
        if verbose:
            print MSG

def sigma_clip(data, sig=3, iters=1, cenfunc=np.ma.median, varfunc=np.var, axis=None, copy=True):
    """
    Froms astropy, had to copy over because using pyraf on ROE machines the astropy version
    is well behind.

    Perform sigma-clipping on the provided data.

    This performs the sigma clipping algorithm - i.e. the data will be iterated
    over, each time rejecting points that are more than a specified number of
    standard deviations discrepant.

    .. note::
        `scipy.stats.sigmaclip
        <http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sigmaclip.html>`_
        provides a subset of the functionality in this function.

    Parameters
    ----------
    data : array-like
        The data to be sigma-clipped (any shape).
    sig : float
        The number of standard deviations (*not* variances) to use as the
        clipping limit.
    iters : int or `None`
        The number of iterations to perform clipping for, or `None` to clip
        until convergence is achieved (i.e. continue until the last
        iteration clips nothing).
    cenfunc : callable
        The technique to compute the center for the clipping. Must be a
        callable that takes in a masked array and outputs the central value.
        Defaults to the median (numpy.median).
    varfunc : callable
        The technique to compute the standard deviation about the center. Must
        be a callable that takes in a masked array and outputs a width
        estimator::

             deviation**2 > sig**2 * varfunc(deviation)

        Defaults to the variance (numpy.var).

    axis : int or `None`
        If not `None`, clip along the given axis.  For this case, axis=int will
        be passed on to cenfunc and varfunc, which are expected to return an
        array with the axis dimension removed (like the numpy functions).
        If `None`, clip over all values.  Defaults to `None`.
    copy : bool
        If `True`, the data array will be copied.  If `False`, the masked array
        data will contain the same array as ``data``.  Defaults to `True`.

    Returns
    -------
    filtered_data : `numpy.ma.MaskedArray`
        A masked array with the same shape as ``data`` input, where the points
        rejected by the algorithm have been masked.

    Notes
    -----
     1. The routine works by calculating::

            deviation = data - cenfunc(data [,axis=int])

        and then setting a mask for points outside the range::

            data.mask = deviation**2 > sig**2 * varfunc(deviation)

        It will iterate a given number of times, or until no further points are
        rejected.

     2. Most numpy functions deal well with masked arrays, but if one would
        like to have an array with just the good (or bad) values, one can use::

            good_only = filtered_data.data[~filtered_data.mask]
            bad_only = filtered_data.data[filtered_data.mask]

        However, for multidimensional data, this flattens the array, which may
        not be what one wants (especially is filtering was done along an axis).

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    a masked array in which all points that are more than 2 *sample* standard
    deviation from the median are masked::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> filtered_data = sigma_clip(randvar, 2, 1)

    This will clipping on a similar distribution, but for 3 sigma relative to
    the sample *mean*, will clip until converged, and does not copy the data::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> from numpy import mean
        >>> randvar = randn(10000)
        >>> filtered_data = sigma_clip(randvar, 3, None, mean, copy=False)

    This will clip along one axis on a similar distribution with bad points
    inserted::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import normal
        >>> from numpy import arange, diag, ones
        >>> data = arange(5)+normal(0.,0.05,(5,5))+diag(ones(5))
        >>> filtered_data = sigma_clip(data, axis=0, sig=2.3)

    Note that along the other axis, no points would be masked, as the variance
    is higher.

    """

    if axis is not None:
        cenfunc_in = cenfunc
        varfunc_in = varfunc
        cenfunc = lambda d: np.expand_dims(cenfunc_in(d, axis=axis), axis=axis)
        varfunc = lambda d: np.expand_dims(varfunc_in(d, axis=axis), axis=axis)

    filtered_data = np.ma.array(data, copy=copy)

    if iters is None:
        i = -1
        lastrej = filtered_data.count() + 1
        while(filtered_data.count() != lastrej):
            i += 1
            lastrej = filtered_data.count()
            do = filtered_data - cenfunc(filtered_data)
            filtered_data.mask |= do * do > varfunc(filtered_data) * sig ** 2
        iters = i + 1
    else:
        for i in range(iters):
            do = filtered_data - cenfunc(filtered_data)
            filtered_data.mask |= do * do > varfunc(filtered_data) * sig ** 2

    return filtered_data

def make_aXe_lis(asn_grism_file, asn_direct_file):
    """
    status = make_aXe_lis(asn_grism_file, asn_direct_file)
    
    Make "inlist" file for aXe routines, with format
    
        grismA_flt.fits directA_flt_1.cat directA_flt.fits 0.0
        grismB_flt.fits directB_flt_1.cat directB_flt.fits 0.0
            
    """
    
    ### get the ASN files:
    asn_direct = ASNFile(asn_direct_file)
    asn_grism = ASNFile(asn_grism_file)

    ### get number of exposures
    n_exposures = len(asn_grism.exposures)

    ### create and open the output file, need to be saved in the root directory:
    outfile = '%s/%s_prep.lis' %(figs.options['ROOT_DIR'], asn_grism_file.split('_asn.fits')[0])
    fp = open(outfile,'w')
    
    for i in range(n_exposures):
        fp.write('%s_flt.fits ' %asn_grism.exposures[i] +
                 '%s_flt_1.cat ' %asn_direct.exposures[i] +
                 '%s_flt.fits 0.0\n' %asn_direct.exposures[i])

    fp.close()

def make_object_id_table():
    """
    assumes run from the ROOT_DIR directory
    """
    
    ### get out relevant parameters
    path_spec = threedhst.options['DRIZZLE_PATH']
    root_grism = figs.options['ROOT_GRISM']

    sexcat = './DATA/%s_drz.cat' %(figs.options['ROOT_DIRECT'])

    ### open the output file:
    file_out = open('./object_id_table.dat', 'w')
    file_out.write('#  idx  ra  dec')
    
    spc_hdu = fits.open('%s/%s_2_opt.SPC.fits' %(path_spec, root_grism))
    
    cat = threedhst.sex.mySexCat(sexcat)
    
    for i in range(len(spc_hdu)-1):
        
        ID = int(spc_hdu[i+1].header['OBJECTID'])
        
        num_list = cat['NUMBER']
        ra_list = cat['X_WORLD']
        dec_list = cat['Y_WORLD']
        
        for j in range(len(num_list)):
            if(num_list[j] == ID):
                ra = ra_list[j]
                dec = dec_list[j]
                file_out.write('%-15d%-15.7f%15.7f\n' %(ID, ra, dec))
    
    file_out.close()