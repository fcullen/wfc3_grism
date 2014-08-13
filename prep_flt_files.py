import figs

from pyraf import iraf
from iraf import stsdas,dither,slitless,axe

from figs.multidrizzle import startMultidrizzle, MultidrizzleRun

import pyfits

import os

INDEF = iraf.INDEF
no = iraf.no
yes = iraf.yes

class fit_2D_background():
    
    def __init__(self, ORDER=-1, x0=None, y0=None, DQMAX=10, IMAGES=['/research/HST/GRISM/CONF/G141_sky_cleaned.fits']):
        """
        __init__(self, ORDER=-1, x0=None, y0=None, DQMAX=10,
         IMAGES=['/research/HST/GRISM/CONF/G141_sky_cleaned.fits'])
    
	    ORDER: Polynomial order of the fit, e.g.
	        ORDER=2 - 0th order, x, y, x**2, y**2, x*y
	        ORDER=3 - 0th order, x, y, x**2, y**2, x*y, x**2 y, x y**2, x**3, y**3
	           
	    x0, y0: reference coordinate for polynomical fit.  Defaults to image center
	    
	    DQMAX: Maximum value in FLT[DQ] extension considered OK
	    
	    IMAGES: Include actual images in the model that is fit.  This can be multiple files, such as the various background images for different sky levels.
    
        """
        self.ORDER = ORDER
        self.x0 = x0
        self.y0 = y0
        self.DQMAX = DQMAX
        self.IMAGES = IMAGES
        self.setup_matrices()
            
    def setup_matrices(self):
        """
		setup_matrices()
    
    	Setup self.A matrix for polynomial fit.
        """
        NX = 1014
        NY = 1014
        
        #### Image matrix indices
        xi,yi = np.indices((NX,NY))
        
        #### Default reference position is image center
        if self.x0 is None:
            self.x0 = NX/2.
        if self.y0 is None:
            self.y0 = NY/2.
                
        xi = (xi-self.x0*1.)/NX
        yi = (yi-self.y0*1.)/NY
        
        NPARAM  = np.sum(np.arange(self.ORDER+2)) #+1 #+1
        NPARAM += len(self.IMAGES)
        self.NPARAM = NPARAM
        
        self.A = np.zeros((NPARAM,NX,NY))
        
        #### Read images to add to the "model"
        count=0
        for img in self.IMAGES:
            hdu = pyfits.open(img)
            self.A[count,:,: ] = hdu[0].data
            hdu.close()
            count += 1
        
        #### zeroth order, flat background
        if self.ORDER >= 0:
            self.A[count,:,:] += 1
            count+=1
        
        for pow in range(1,self.ORDER+1):
            pi = pow-1

            #### Cross terms
            while (pi > pow/2.):
                self.A[count,:,:] = xi**pi*yi**(pow-pi)
                print 'A[%d,:,:] = xi**%d*yi**%d' %(count,pi,pow-pi)
                count+=1
                self.A[count,:,:] = xi**(pow-pi)*yi**pi
                print 'A[%d,:,:] = xi**%d*yi**%d' %(count,pow-pi,pi)
                count+=1
                pi-=1
            
            #### x**pow/2 * y**pow/2 term
            if (pow/2. == np.int(pow/2.)):
                self.A[count,:,:] = xi**(pow/2)*yi**(pow/2)
                print 'A[%d,:,:] = xi**%d*yi**%d' %(count,pow/2,pow/2)
                count+=1
            
            #### x**pow, y**pow terms
            print 'A[%d,:,:] = xi**%d' %(count,pow)
            self.A[count,:,:] = xi**pow
            count+=1
            print 'A[%d,:,:] = yi**%d' %(count,pow)
            self.A[count,:,:] = yi**pow
            count+=1
        
        # #### Oth order for `grism` is True is G141 median image
        # #medF = pyfits.open('../CONF/WFC3.IR.G141.sky.V1.0.fits') # from aXe web
        # # cleaned of value=0 pixels
        # medF = pyfits.open(GRISM_SKY)
        # med_g141 = medF[0].data
        # medF.close()
        # 
        # self.B = self.A*1.
        # self.B[0,:,:] = med_g141*1.
        # self.NPARAM = NPARAM
    
    def fit_image(self, root, A=None, overwrite=False, show=True, save_fit=False):
        """
		fit_image(self, root, A=None, overwrite=False, show=True)
    
	    Fit and optionally subtract the background from a FLT image.
	    
	    `root` is like "ib3728d2q"
	    
	    `A` is a matrix computed by `setup_matrices` or `__init__`.
	    
	    if `overwrite` is True:
	        Write the background-subtracted image to the original FLT file
	        
	    if `show` is True:
	        Show a plot of the original, background, and bg-subtracted images.
        """
        import os
        import glob
        
        if A is None:
            A = self.A
        
        #### Read the FLT file and get dimensions 
        #### (should always be 1014x1014 for WFC3)
        fi = pyfits.open(root+'_flt.fits',mode='update')
        IMG = fi[1].data
        DQ = fi[3].data
        NX,NY = IMG.shape
        
        #### array indices for each pixel position
        xi,yi = np.indices((NX,NY))
        xi = (xi-NX)/2./NX
        yi = (yi-NY)/2./NY
        
        #### If a segmentation image is available, read it
        #### for use as an object mask.
        segfile = glob.glob(root+'_flt.seg.fits*')
        if len(segfile) == 0:
            seg = IMG*0.
        else:
            print 'Segmentation image: %s' %segfile[0]
            fis = pyfits.open(segfile[0])
            seg = fis[0].data
        
        #### Apply segmentation mask, also mask out extreme IMG values 
        #### and any pixel with DQ flag > self.DQMAX
        
        q = np.where((seg == 0) & (IMG > -1) & (IMG < 4) & (DQ < self.DQMAX)) 
        qb = np.where((seg > 0) | (IMG < -1) | (IMG > 4) | (DQ >= self.DQMAX))
        IMGb = IMG*1.
        IMGb[qb] = np.nan
        
        #### Apply mask to IMG and fit matrices
        Aq = np.transpose(A[:,q[0],q[1]])
        IMGq = IMG[q[0],q[1]]
        
        #### Get fit parameters with least-sq. fit.
        p, resid, rank, s = scipy.linalg.lstsq(Aq,IMGq)

        print p
        
        #### Create the bg fit image from the fit parameters
        IMGout = IMG*0.
        for i in range(self.NPARAM):
            IMGout += A[i,:,:]*p[i]
        print 'Done'
        
        #### Save fit parameters to an ASCII file
        fp = open(root+'_flt.polybg','w')
        for pi in p:
            fp.write('%13.5e\n' %pi)
        fp.close()
        
        #### Show the results, note that subsequent 
        #### plots aren't cleared from memory, so memory 
        #### fills up quickly with repeated calls with show=True.
        if show:
            dshow = 0.3
            plt.figure()
            plt.subplot(221)
            plt.imshow(IMGb,vmin=np.median(IMGb)-dshow,
                vmax=np.median(IMGb)+dshow)
            plt.subplot(222)
            plt.imshow(IMG-IMGout,vmin=0-dshow,vmax=dshow)
            plt.subplot(223)
            plt.imshow(IMGout,vmin=np.median(IMGb)-dshow,
                vmax=np.median(IMGb)+dshow)
            plt.subplot(224)
            plt.imshow(DQ,vmin=0,vmax=10)
        
        if save_fit:
            hdu0 = pyfits.PrimaryHDU(header=fi[0].header)
            hdu1 = pyfits.ImageHDU(data=IMG, header=fi[1].header)
            hdu2 = pyfits.ImageHDU(data=IMGb, header=fi[1].header)
            hdu3 = pyfits.ImageHDU(data=IMGout, header=fi[1].header)
            save_im = pyfits.HDUList([hdu0,hdu1,hdu2,hdu3])
            #save_im[0].data = IMGout
            id=len(glob.glob(root+'_flt.BG*fits'))
            save_im.writeto(root+'_flt.BG_%02d.fits' %(id+1), clobber=True)
        
        #### Subtract fitted background, write
        #### bg-subtracted image to original FLT file if `overwrite` is True
        FIX = IMG-IMGout
        if overwrite:
            print 'Overwrite: '+root
            fi[1].data = FIX
            fi.flush()
        
        #### Save images to self.
        self.FIX = FIX
        self.IMG = IMG
        self.MODEL = IMGout
        fi.close()

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

def process_grism_direct_pair(asn_direct_file='ibhj34030_asn.fits',
	asn_grism_file='ibhj34040_asn.fits',
	adjust_targname=True,
	align_image='/disk1/fc/research/HST/CANDELS/RAW/gsa_6epoch_wfc3_f125w_060mas_v0.2_drz.fits',
	align_extension=0,
	align_geometry=True,
	path_to_raw='../RAW',
	images=['../CONF/G141_sky_cleaned.fits',
		'../CONF/G141wHI_fixed_sky.fits',
		'../CONF/G141wLO_fixed_sky.fits'],
	skip_grism=False,
	skip_direct=False,
	get_shift=True,
	tweakshifts_only=False,
	direct_higher_order=2,
	grism_higher_order=1,
	save_fit=False,
	second_pass=True,
	overall=True,
	sky_images=['sky_cosmos.fits', 'sky_goodsn_lo.fits', 
			'sky_goodsn_hi.fits', 'sky_goodsn_vhi.fits']):


	#### first adjust the names to make them sesnible
	#### e.g. 'ibhm51030_asn.fits' --> 'COSMOS-3-F140W_asn.fits'
	if adjust_targname:
		asn_direct_file = make_targname_asn(asn_direct_file)
		asn_grism_file = make_targname_asn(asn_grism_file)

	print 'DIRECT: %s, GRISM: %s\n' %(asn_direct_file, asn_grism_file)

	#### now process the direct images:
	if not skip_direct:

		# copy over the files and check best flat-frame etc..
		fresh_flt_files(asn_direct_file, from_path=path_to_raw)

		#  make region files for the pointing
		if not os.path.exists(asn_direct_file.replace('fits','pointing.reg')):
		    figs.regions.asn_region(asn_direct_file)

		figs.prep_flt_files.prep_flt(asn_file=asn_direct_file,
										get_shift=get_shift, first_run=True,
										bg_only=False, bg_skip=False, redo_background=True,
										align_image=align_image, 
										align_ext=align_extension,
										skip_drz=False, final_scale=0.06, pixfrac=0.8,
										images=[],
										align_geometry=align_geometry, clean=True,
										initial_order=0, save_fit=save_fit,
										tweakshifts_only=tweakshifts_only, make_persistence_mask=True)

		if DIRECT_HIGHER_ORDER > 0:
			figs.prep_flt_files.prep_flt(asn_file=asn_direct_file,
											get_shift=False, first_run=False,
											bg_only=False, bg_skip=False, redo_background=False,
											skip_drz=False, final_scale=0.06, pixfrac=0.8,
											IMAGES=[], clean=True,
											initial_order=DIRECT_HIGHER_ORDER, save_fit=save_fit)

def prep_flt(asn_file=None, get_shift=True, bg_only=False, bg_skip=False,
	first_run=True, redo_background=True,
	align_image='../ACS/h_nz_sect*img.fits', align_ext=0, 
	skip_drz=False, final_scale=0.06, pixfrac=0.8,
	images=['/research/HST/GRISM/CONF/G141_fixed_sky.fits'],
	align_geometry='shift', clean=True,
	initial_order=-1, save_fit=False,
	tweakshifts_only=False,
	oned_background=True, make_persistence_mask=False):
    """
		prep_flt(asn_file=None, get_shift=True, bg_only=False,
	            redo_background=True)

	    
	    Subtract background and align WCS of direct/grism FLT files.
	    
	    1) Apply the DQ masks defined in the *mask.reg files, as created
	       by figs.dq
	    
	    2) First pass on background subtraction 
	        
	        o 0th order is median background image, e.g. G141 sky 
	        
	        o 1-nth order is polynomial fit with x-y cross terms.
	    
	        o [if `bg_only` is True then return]
	        
	    3) Run tweakshifts  [if `get_shift` is True & `grism` is False]
	    
	    4) Run Multidrizzle with first guess tweakshifts
	    
	    5) Get object mask for FLT files:
	    
	        o Blot DRZ output to all the individual FLT frames
	        
	        o Run SExtractor on blot images for object (segmentation)
	          mask  
	    
	    6) Use figs routines to align the F140W direct image to
	       ACS reference.   [if `get_shift` is True & `grism` is False]
	       
	    7) Redo background subtraction with improved object mask defined
	       by the segmentation images [if redo_background is True]
	       
	    8) Subtract the collapsed background to fix residuals from the grism sky fit
	    
	    9) Run multidrizzle again [if redo_background is True]
    """

    import os
    import glob
        
    if bg_skip:
        bg_only=False
        redo_background=False
    
    # make and asn object from the asn file supplied:
    ROOT_DIRECT = asn_file.split('_asn.fits')[0]
    asn = figs.utils.ASNFile(asn_file)
    
    # first guess at shifts
    if get_shift:

        figs.showMessage("Running tweakshifts first pass ...", warn=False)

        figs.shifts.run_tweakshifts(asn_file, verbose=True)
        figs.shifts.checkShiftfile(asn_file)

        # corrects for the faults in pre-April 2011 ICDTAB files
        # if necessary. See the function docstring for details:
        figs.shifts.default_rotation(asn_file, path_to_flt='./')
        
        # now update the shift files to align the direct images to a
        # reference image:
        if not tweakshifts_only:

        	if not align_image:
        		figs.showMessage("No align image file supplied", warn=True)
        		return False

        for ig, geom in enumerate(align_geometry.split(',')):

            first = ig == 0

            startMultidrizzle(asn_file, use_shiftfile=True, skysub=True,
            final_scale=final_scale, pixfrac=pixfrac, driz_cr=first,
            updatewcs=0, clean=True, median=first)

            figs.shifts.refine_shifts(ROOT_DIRECT=ROOT_DIRECT,
            ALIGN_IMAGE=align_image, 
            ALIGN_EXTENSION = align_ext,
            fitgeometry=geom.strip(), clean=clean)
                   
        # need fresh FLT files now
        fresh_flt_files(asn_file)
             
    # first pass background subtraction to individual exposires
    # fits to background without a segmentation map
    if not bg_skip:
        # set up matrix for fitting
        fit = fit_2D_background(ORDER=initial_order,
                                IMAGES=images)
        
        for exp in asn.exposures:
            fit.fit_image(exp, A=fit.A, show=False, overwrite=True,
                          save_fit=save_fit)
    
    # stop here if only want background subtraction
    if bg_only:
        return
    
    # not make a the drizzled file:        
    if not skip_drz:
        if first_run:
            # first MDRZ run needs native pixels to avoid flagging stars as CRs
            startMultidrizzle(asn_file, use_shiftfile=True, 
                skysub=True, final_scale=0.128254, pixfrac=1.0,
                driz_cr=True, median=True, updatewcs=True)
        
        if (final_scale != 0.128254) | (pixfrac != 1.0) | (not first_run):
            startMultidrizzle(asn_file, use_shiftfile=True, 
                skysub=bg_skip, final_scale=final_scale, pixfrac=pixfrac,
                driz_cr=False, median=False, updatewcs=False)
     
    #### blot combined images back to reference frame and make a 
    #### segmentation mask
    run = MultidrizzleRun((asn_file.split('_asn.fits')[0]).upper())
    
    #### ACS has entries in run file for each of two WFC chips
    flt = pyfits.open(asn.exposures[0]+'_flt.fits')
    inst = flt[0].header.get('INSTRUME').strip()
    if inst is 'ACS':
        skip=2
    else:
        skip=1
        
    if redo_background:
        for i,exp in enumerate(asn.exposures):
            run.blot_back(ii=i*skip, copy_new=(i is 0))
            make_segmap(run.flt[i])
        
        ### Flag bright stars in segmentation map
        asn_mask = asn_file+'.mask.reg'
        if os.path.exists(asn_mask):
            figs.showMessage("Apply ASN mask: %s" %(asn_mask))
            for flt in asn.exposures:
                seg = flt+'_flt.seg.fits'
                figs.regions.apply_dq_mask(seg, extension=0, 
                    mask_file = asn_mask)
                            
    #### Run BG subtraction with improved mask and run multidrizzle again
    if redo_background:
        fit = fit_2D_background(ORDER=initial_order, IMAGES=IMAGES)
        for exp in asn.exposures:
            #### 2-D background, fit
            fit.fit_image(exp, A=fit.A, show=False, overwrite=True, 
                          save_fit=save_fit)
            
            #### 1-D background, measured
            if oned_background:
                print '\n Extract 1D background \n'
                test = oned_grism_background_subtract(exp, nbin=26, savefig=True, verbose=False)
            
        startMultidrizzle(asn_file, use_shiftfile=True, skysub=False,
            final_scale=final_scale, pixfrac=pixfrac, driz_cr=False,
            updatewcs=False, median=False, clean=clean)
    
    if clean:
        files=glob.glob('*BLOT*')
        files.extend(glob.glob('drz_???.fits'))
        for file in files: 
            #print 'rm '+file
            os.remove(file)