import figs

from pyraf import iraf

import numpy as np

from astropy.io import fits

import os
import glob

import shutil

def run_sregister_for_align_image(mosiac_drz):
	"""
	Run sregister on the drizzled direct mosaic and the CANDELS mosaic to
	cut out an area to use for aligning the images with tweakreg()
	"""

	### get the root name:
	root = figs.options['ROOT_DIRECT']

	### iraf flpr()
	figs.utils.iraf_flpr()

	### first unlearn the routine:
	iraf.unlearn('sregister')

	### remove previous versions if they exist:
	try:
		os.remove('%s_align_reference.fits' %(root))
	except:
		pass

	### run sregister, put in a try-except because
	### sregister crashes once it has completed when attempting
	### to delete files it creates. By ignoring this error there is
	### no change to the actual output:
	try:
		iraf.sregister(input=mosiac_drz,
					   reference='%s_drz.fits[SCI]' %(root),
					   output='%s_align_reference.fits' %(root),
					   verbose=False)
	except:
		pass


	### cleanup files:
	tmps_files = glob.glob('tmps*')
	for tmps in tmps_files:
		os.remove(tmps)

	try:
		os.remove('sregister.db')
	except:
		pass

def run_sregister_for_detection_image(asn_direct_file, sci_image, wht_image):
	"""
	Run sregister on the drizzled direct mosaic and the CANDELS mosaic to
	cut out an area to use for creating source catalogue.
	"""

	### get root name for the reference image:
	root = figs.options['ROOT_DIRECT']

	images = [sci_image, wht_image]
	extensions = ['SCI', 'WHT']

	for im, ext in zip(images, extensions):

		### iraf flpr()
		figs.utils.iraf_flpr()

		### first unlearn the routine:
		iraf.unlearn('sregister')

		### remove previous versions if they exist:
		try:
			os.remove('%s_detection_%s.fits' %(figs.options['DETECTION_BAND'], ext))
		except:
			pass

		try:
			os.remove('%s_detection_%s.fits' %(figs.options['DETECTION_BAND'], ext))
		except:
			pass

		### run sregister, put in a try-except because
		### sregister crashes once it has completed when attempting
		### to delete files it creates. By ignoring this error there is
		### no change to the actual output:
		iraf.sregister(input=im,
					   reference='%s_drz.fits[SCI]' %(root),
					   output='%s_detection_%s.fits' %(figs.options['DETECTION_BAND'], ext),
					   verbose=False)



	### cleanup files:
	tmps_files = glob.glob('tmps*')
	for tmps in tmps_files:
		os.remove(tmps)

	try:
		os.remove('sregister.db')
	except:
		pass

def run_tweakshifts_on_direct_exposures(asn_direct_file, verbose=True):
	"""
	Use the drizzlepac software tweakreg to compute offsets between
	the direct images.
	"""

	root = figs.options['ROOT_DIRECT']

	### iraf flpr()
	figs.utils.iraf_flpr()

	### unlearn the routine:
	iraf.unlearn('tweakshifts')

	### run tweakreg in non-interactive mode:
	iraf.tweakshifts(input= asn_direct_file, 
							shiftfile='',
							reference='%s_tweak.fits' %root,
							output = '%s_initial_shifts.txt' %root, 
							findmode = 'catalog',
							gencatalog = 'daofind', 
							sextractpars = '', 
							undistort = True, 
							computesig = True, 
							idckey = 'idctab',
							clean = True, 
							verbose = False, 
							catfile = '', 
							xcol = 1, 
							ycol = 2,
							fluxcol = 3, 
							fluxmax = iraf.INDEF, 
							fluxmin = iraf.INDEF, 
							fluxunits = 'counts',
							nbright = iraf.INDEF, 
							refcat = '', 
							refxcol = 1, 
							refycol = 2, 
							rfluxcol = 3,
							rfluxmax = iraf.INDEF, 
							rfluxmin = iraf.INDEF, 
							rfluxunits = 'counts',
							refnbright = iraf.INDEF, 
							minobj = 15, 
							nmatch = 30,
							matching = 'tolerance',
							xyxin = iraf.INDEF, 
							xyyin = iraf.INDEF, 
							tolerance = 4.0, 
							fwhmpsf = 1.5,
							sigma = 0.0, 
							datamin = iraf.INDEF, 
							datamax = iraf.INDEF, 
							threshold = 4.0,
							nsigma = 1.5, 
							fitgeometry = 'shift', 
							function = 'polynomial',
							maxiter = 3, 
							reject = 3.0, 
							crossref = '', 
							margin = 50, 
							tapersz = 50,
							pad = False, 
							fwhm = 7.0, 
							ellip = 0.05, 
							pa = 45.0, 
							fitbox = 7,
							Stdout=1)

def align_direct_to_reference(verbose=True, n_iter=20, drizzled_image=True):
	"""
	Use iraf software geomap to get shift solutions between images.

	n_iter = number of times to iterate over the alignment routines
	drizzled_image = boolean to distinguish between a drizzled image used for alignment,
	                 in which case the shifts must be transferred back to oiginal flt.
	"""

	### get the root name:
	root = figs.options['ROOT_DIRECT']

	### get the alignment image which has been produced
	### with run_sregister_to_cutout_CANDELS_region()
	align_image = '%s_align_reference.fits' %(root)

	### now run SExtractor on the direct and refereance
	### images to build up 2 catalogs:
	se = figs.sex.SExtractor()
	se.aXeParams()
	se.copyConvFile()
	se.overwrite = True
	se.options['CHECKIMAGE_TYPE'] = 'NONE'
	se.options['FILTER']    = 'Y'
	se.options['DETECT_THRESH']    = '50' 
	se.options['ANALYSIS_THRESH']  = '50' 
	se.options['MAG_ZEROPOINT'] = '%.2f' %(figs.options['MAG_ZEROPOINT'])
	se.options['DETECT_MINAREA'] = '100'

	### generate the direct image catalog:
	se.options['CATALOG_NAME']    = 'direct.cat'
	iraf.imcopy('%s_drz.fits[SCI]' %(root), "SCI.fits", verbose=False)

	## if not using drizzled image for alignment done't include
	### a weight image:
	if drizzled_image:
		se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
		se.options['WEIGHT_IMAGE']    = 'WHT.fits'
		iraf.imcopy('%s_drz.fits[WHT]' %(root), "WHT.fits", verbose=False)
	else:
		se.options['WEIGHT_TYPE']     = 'NONE'
		se.options['WEIGHT_IMAGE']    = 'WHT.fits'

	status = se.sextractImage('SCI.fits')

	### generate the alignment image catalog:
	se.options['CATALOG_NAME']    = 'align.cat'
	se.options['WEIGHT_TYPE']     = 'NONE'
	se.options['WEIGHT_IMAGE']    = 'WHT.fits'
	status = se.sextractImage(align_image)

	### Read the catalogs
	direct_cat = figs.sex.mySexCat('direct.cat')
	align_cat = figs.sex.mySexCat('align.cat')

	### initialize x,y shift parameters so can be 
	### updated with each iteration, the x,y shifts 
	### will converge twoard an optimal value over the
	### iterations:
	xshift = 0
	yshift = 0

	### empty line in the output to make things clearer:
	print ""

	### now loop through the process:
	iteration = 0
	while(iteration < n_iter):

		print "Running matching algorithm on iteration #%d" %(iteration)

		### Get x,y coordinates of detected objects direct image
		fp = open('direct.xy','w')
		for i in range(len(direct_cat.X_IMAGE)):
			fp.write('%-15s%-15s\n' %(direct_cat.X_IMAGE[i],direct_cat.Y_IMAGE[i]))
		fp.close()

		### Get x,y coordinates of detected objects alignment image
		fp = open('align.xy','w')
		for i in range(len(align_cat.X_IMAGE)):
			fp.write('%-15s%-15s\n' %(np.float(align_cat.X_IMAGE[i])+xshift, np.float(align_cat.Y_IMAGE[i])+yshift))
		fp.close()

		### iraf flpr()
		figs.utils.iraf_flpr()

		### remove previous solution:
		if iteration > 0:
			os.remove('align.match')

		### get the alignment catalog:
		status = iraf.xyxymatch(input="direct.xy", 
								reference="align.xy",
								output="align.match",
								tolerance=figs.options['ALIGN_TOLERANCE'], 
								separation=0, 
								verbose=True, 
								Stdout=1)

		### iraf flpr()
		figs.utils.iraf_flpr()

		### now get shifts with geomap:
		if iteration > 0:
			os.remove('align.map')

		iraf.geomap(input="align.match", 
					database="align.map",
					fitgeometry="rscale", 
					interactive=False, 
					xmin=iraf.INDEF, 
					xmax=iraf.INDEF, 
					ymin=iraf.INDEF, 
					ymax=iraf.INDEF,
					maxiter = 10, 
					reject = 2.0, 
					Stdout=1)

		### get the output from geomap:
		fp = open("align.map", "r")
		for line in fp.readlines():
			spl = line.split()
			if spl[0].startswith('xshift'):
				xshift += float(spl[1])    
			if spl[0].startswith('yshift'):
				yshift += float(spl[1])    
			if spl[0].startswith('xrotation'):
				rot = float(spl[1])    
			if spl[0].startswith('xmag'):
				scale = float(spl[1])    
			if spl[0].startswith('xrms'):
				xrms = float(spl[1])    
			if spl[0].startswith('yrms'):
				yrms = float(spl[1])    
		fp.close()

		### update iteration counter:
		iteration += 1

		print 'Shift iteration #%d, xshift=%f, yshift=%f, rot=%f, scl=%f (rms: %5.2f,%5.2f)' %(iteration, xshift, yshift, rot, scale, xrms, yrms)

	### copy the final align.map for posterity:
	shutil.copy('align.map', '%s_align.map' %(root))

	## cleanup from the alignment process:
	remvfiles = ['SCI.fits','WHT.fits','align.cat',
				'align.map','align.match','align.reg','align.xy',
				'direct.cat','direct.reg','direct.xy',
				'drz_sci.fits','drz_wht.fits','bg.fits', 'imxymatch.1', 'sex_stderr',
				'figs_auto.sex', 'figs_auto.param', 'default.nnw', 'default.conv']

	for file in remvfiles:
		try:
			os.remove(file)
		except:
			pass

	if drizzled_image:
		#### shifts measured in drz frame.  Translate to the flt frame:
		drz = fits.open('%s_drz.fits' %(root))

		#### Get reference angle from first image in the ASN file:
		asn = figs.utils.ASNFile('%s_asn.fits' %(root))
		alpha = (180. - fits.getheader(asn.exposures[0]+'_flt.fits', 1)['PA_APER']) / (360. * 2 * np.pi)

		### Get the drizzle scale from the MultiDrizzle '.run' file:
		for line in open('%s.run' %(root),'r'):
			if line.startswith('drizzle.scale'):
				drizzle_scale = line.split()[2]

		print drizzle_scale

		### get the equivalent shifts in the FLT frames:
		xsh = (xshift*np.cos(alpha) - yshift*np.sin(alpha))*np.float(drizzle_scale)
		ysh = (xshift*np.sin(alpha) + yshift*np.cos(alpha))*np.float(drizzle_scale)

		print 'Final shift:', xsh, ysh, drz[1].header['PA_APER']
	else:
		xsh = xshift
		ysh = yshift

	fp = open('%s_align.info' %(root),'w')
	fp.write('%s %8.3f %8.3f %8.3f\n' %(align_image, xsh, ysh, rot)) 
	fp.close()

	#### Read the shiftfile
	if drizzled_image:
		shiftF = ShiftFile('%s_initial_shifts.txt' %(root))

		#### Apply the alignment shifts to the shiftfile
		#### NB there's a bug in multidrizzle so that sometimes it only applies
		#### either the x or y shift but not both.
		#### Therefore make 2 shiftfiles, one with only x-shift, one with only
		#### y-shift. And run multidrizlle twice for both shifts to be applied!
		shiftF.xshift = list(np.array(shiftF.xshift)-np.array(shiftF.xshift))
		shiftF.yshift = list(np.array(shiftF.yshift)-ysh)
		shiftF.rotate = list((np.array(shiftF.rotate)+rot) % 360)
		shiftF.scale = list(np.array(shiftF.scale)*scale)

		shiftF.write('%s_final_shifts_yshift.txt' %(root))

		#### Apply the alignment shifts to the shiftfile
		shiftF.xshift = list(np.array(shiftF.xshift)-xsh)
		shiftF.yshift = list(np.array(shiftF.yshift)-np.array(shiftF.yshift))
		shiftF.rotate = list((np.array(shiftF.rotate)+rot) % 360)
		shiftF.scale = list(np.array(shiftF.scale)*scale)

		shiftF.write('%s_final_shifts_xshift.txt' %(root))

	else:
		### use the default shift file in figs data folder:
		shiftF = ShiftFile('/disk1/fc/FIGS/figs/data/default_shift_file.txt')

		### add the reference image as needed by multidrizzle:
		shiftF.headerlines[1] = '# refimage: %s \n' %(align_image)

		#### Apply the alignment shifts to the shiftfile
		shiftF.xshift = [np.array(xsh)]
		shiftF.yshift = [np.array(ysh)]
		shiftF.rotate = [np.array(rot) % 360]
		shiftF.scale = [np.array(scale)]

		shiftF.write('%s_final_shifts.txt' %(root))

class ShiftFile():
	"""
	ShiftFile(infile)

	Read and manipulate a shiftfile produced by, e.g., tweakshifts

	Based on aXe2html.sextractcat
	"""
	def __init__(self, filename):
		linelist = self.opencat(filename)
		self.headerlines = self.extractheader(linelist)
		rowlines    = self.extractrows(linelist)
		self.processrows(rowlines)
		self.nrows = len(rowlines)
		self.filename = filename

	def opencat(self, filename):
		"""
		Input:
		filename - the name of the sextractor ascii catalog

		Return:
		linelist - a list wit all lines of the sextractor
		catalog

		Description:
		Opens the file and stores all the lines in a list.

		"""
		listfile = open(filename,'r')
		linelist = listfile.readlines()
		listfile.close()

		return linelist

	def extractheader(self, linelist):
		"""
		Input:
		linelist - all lines of a sextractor catalog

		Return:
		headerlines - the lines which contain header
		information

		Description:
		Extracts the header lines from the list of lines
		Header lines have a '#' as the first digits and are
		longer than two characters. This allows for small
		cosmetical changes to the original Sextractor tables
		"""
		headerlines = []
		for index in range(len(linelist)):
			oneline = linelist[index]
			if (oneline[0] == "#") and (len(oneline)) > 2:
				headerlines.append(oneline)
		return headerlines

	def extractrows(self, linelist):
		"""
		Input:
		linelist - all lines of a sextractor catalog

		Return:
		rowlines - the content lines of a sextractor
		catalog

		Description:
		Extracts the content lines from the list of lines.
		Content lines have a all lines except the ones which start with a 
		'#' and are longer than one character. This allows for
		e.g. whitespaces as the last 'rows', which often remain 
		after editing an ascii-file

		"""
		rowlines = []
		for index in range(len(linelist)):
			oneline = linelist[index]
			if oneline[0] != "#" and len(oneline) > 1:
				rowlines.append(oneline)        
		return rowlines

	def processrows(self, rowlines):
		"""
		processrows(self, rowlines)

		Read: image xshift yshift rotate scale from shiftfile.

		"""
		nlines = len(rowlines)
		self.images = []
		self.xshift = []
		self.yshift = []
		self.rotate = []
		self.scale = []
		for i in range(nlines):
			line = rowlines[i].split()
			self.images.append(line[0])
			self.xshift.append(float(line[1]))
			self.yshift.append(float(line[2]))
			if len(line) > 3:
				self.rotate.append(float(line[3]))
			else:
				self.rotate.append(0.)
				self.scale.append(1.0)
			if len(line) > 4:
				self.scale.append(float(line[4]))

	def write(self, outfile):
		"""write(outfile)"""
		fp = open(outfile,'w')
		fp.writelines(self.headerlines)
		for i in range(self.nrows):
			line = '%s    %.4f  %4f    %.3f    %.3f\n' %(self.images[i],
														 self.xshift[i], 
														 self.yshift[i], 
														 self.rotate[i], 
														 self.scale[i])
			fp.write(line)
		fp.close()

	def pop(self,idx):
		"""
		pop(self,idx)
		"""
		out = self.images.pop(idx) 
		out = self.xshift.pop(idx) 
		out = self.yshift.pop(idx) 
		out = self.rotate.pop(idx) 
		out = self.scale.pop(idx) 
		self.nrows -= 1

	def append(self,image, xshift=0., yshift=0., rotate=0.0, scale=1.0):
		"""
		append(self,image, xshift=0., yxhift=0.,
		   rotate=0.0, scale=1.0)
		"""
		self.images.extend([image])
		self.xshift.extend([xshift])
		self.yshift.extend([yshift])
		self.rotate.extend([rotate])
		self.scale.extend([scale])
		self.nrows += 1




	

 

