import figs

from drizzlepac import astrodrizzle

def direct_astrodrizzle_run(asn_direct):
    """
    Performs an astrodrizzle run on a set of input direct images

    Some CANDELS defaults:

    driz_cr_snr='6 3.0',
    driz_cr_scale='1.6 0.7'
    """

    root = asn_direct.split('_asn.fits')[0]

    astrodrizzle.AstroDrizzle(input=asn_direct, 
    						  output='', 
    						  runfile='%s_astrodrizzle.log' %(root),
    						  restore=False, 
    						  preserve=True, 
    						  overwrite=False, 
    						  clean=True, 
    						  static=True,
    						  static_sig=4.0, 
    						  skysub=True, 
    						  driz_separate=True, 
    						  driz_sep_kernel='turbo',
    						  driz_sep_wcs=False, 
    						  median=True,
    						  blot=True,
    						  driz_cr=True,
    						  driz_cr_snr='6 3.0',
    						  driz_cr_scale='1.6 0.7',
    						  driz_combine=True,
    						  final_wht_type='IVM',
    						  final_kernel='turbo',
    						  final_pixfrac=0.8,
    						  final_scale=0.06,
    						  final_rot=0.0,
    						  final_wcs=True)


