# -*- coding: utf-8 -*-

# Calculate the completeness as a function of magnitude for an astronomical image.
#	-	Carl J. Mundy, July 2015
# Based on scripts by Kenneth J. Duncan, Dec 2012.

import os, numpy as np, glob, shutil, sys, subprocess
from astropy.table import Table, vstack
from astropy.io import fits
from pyraf.iraf import noao, images
from pyraf.iraf import artdata as art
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator

###########################################
# QUICK FLAGS
doMocks			= False
doSex			= False
doCat			= False
doStats			= True
###########################################
# INPUT IMAGE PARAMETERS
inputImage 		= '/home/ppxcjm/data/video/images/CFHTLS-D1_Ks_maxseeing0p90_2015-03-25_cfht.fits'#
inputImageZP 	= 30.0
inputImageGain 	= 1207.06
inputImageExp 	= 1.0
outputPath 		= '/home/ppxcjm/data/video/completeness/output/'
splitRatio 		= 3. # Images < 500MB seem to be okay.
###########################################
# PSF PARAMETERS
psfTablePath 	= './CFHTLS/psf_list_CFHTLS_Ks.txt'
psfPath 		= './CFHTLS/PSF_Ks/'
psfPref 		= 'psf_Ks_CFHTLS_*' ## PSF files are psf_Ks_CFHTLS_XXX_YYY.fits
psfSize 		= 75/2. #pixels #radius of *file* in pixels
###########################################
# SIMULATED GALAXY PROPERTIES
numSims 		= 2
gallist_out 	= 'test.txt'
outPrefix 		= 'video_'
ngals 			= 7000
poisson			= False
###########################################
# SIMULATED GALAXY DISTRIBUTIONS
sersic_dist 	= np.clip(0.5*np.round(np.random.lognormal(0.7, 0.6, size=ngals)*2.), 0.5, 10.) # crop the distribution before science
axis_dist 		= np.clip(np.random.normal(0.8, 0.2, size=ngals), 0.1, 1.)
mag_dist 		= np.clip(np.random.lognormal(3.1, 0.1, size=ngals), 17., 27.)
pa_dist			= np.random.rand(ngals)*360.
radii_dist 		= np.clip(np.random.gamma(1.5, 1.5, size=ngals)*0.723, 0.5, 1000.)
###########################################
# SOURCE EXTRACTOR PARAMETERS
sexCmd			= 'sex'
sexConfig		= '/home/ppxcjm/data/video/sextractor_files/video_cfht.sex'
sexOutFile		= 'test.cat'
###########################################
# TOPCAT PARAMETERS
topCmd			= 'topcat -stilts'
topMaxSep		= 20 #pixels
topSave			= 'matched_objects.fits'
###########################################
# STATISTICS PARAMETERS
statMagBins		= np.arange(17, 27, 0.25)
statRMax		= 5 #pixels
statMagMax		= 0.3 #mag
###########################################

# Generate a table of x, y, psf_file
PSFfiles = glob.glob(psfPath+psfPref)
PSF_x = [ int(i[:-5].split('_')[-2]) for i in PSFfiles]
PSF_y = [ int(i[:-5].split('_')[-1]) for i in PSFfiles]
PSF_file = [ os.path.split(i)[-1] for i in PSFfiles]
psfTable = Table()
psfTable['x'], psfTable['y'] = PSF_x, PSF_y
psfTable['file'] = PSF_file

# Read in the header of our image
header = fits.getheader(inputImage, 0)
NAXIS1, NAXIS2 = header['NAXIS1'], header['NAXIS2']

if doMocks:
	# Loop through simulations
	print 'Creating mock catalogues...'
	for s in range(numSims):
		# Output image for this simulation
		outputImage = outputPath + outPrefix + str(s+1) + '.fits'

		# Generate random (x,y) co-ordinates to put the mock galaxies in
		header = fits.getheader(inputImage, 0)
		NAXIS1, NAXIS2 = header['NAXIS1'], header['NAXIS2']
		ra = np.random.randint(100, NAXIS1-100, ngals)
		dec = np.random.randint(100, NAXIS2-100, ngals)

		# Create the table of random data for this run
		simTable = Table()
		simTable['x'] = ra
		simTable['y'] = dec
		simTable['mag'] = mag_dist
		simTable['sersic'] = sersic_dist
		simTable['radii'] = radii_dist
		simTable['axis_ratio'] = axis_dist
		simTable['pa'] = pa_dist
		simTable.write(outputPath+'simulated_galaxies.{0}.tab'.format(s), format='ascii.commented_header')

		# Write out a param file for IRAF's mkobject to understand
		with open(gallist_out+'.modx', 'w') as outfile:
			for i in range(ngals):
				outfile.write('{0:>9.1f} {1:>9.1f} {2:>9.3f} {3:>s} {4:>9.3f} {5:>9.3f} {6:>9.1f} {7}'
					.format(simTable['x'][i],
							simTable['y'][i],
							simTable['mag'][i],
							'sersic'+str(simTable['sersic'][i]),
							simTable['radii'][i],
							simTable['axis_ratio'][i],
							simTable['pa'][i],
							'\n'))

		# Split the image up here if there are problems.
		split = int(np.round(NAXIS2 / splitRatio))
		outfilelist = []
		for si in range(int(splitRatio)):
			# Sort out the splitting up of the image
			outfile = outputImage+'.split{0}.fits'.format(si)
			outfilelist.append(outfile)
			infile = inputImage+'[1:{0},{2}:{1}]'.format(NAXIS1, split*(si+1)+1, (split*si)+1)
			
			# Copy the split files to the output directory.
			if os.path.isfile(outfile):
				os.remove(outfile)
			images.imcopy(input=infile, output=outfile, Stdout=1)

			# If there was no spatial PSF dependence, we could run mkobject here. But bum.
			# Choose a random PSF file
			psfidx = np.random.randint(0, len(psfTable)-1)
			psfImage = os.path.split(psfTable['file'][psfidx])
			psfFile = psfPath + psfImage[-1]
			psfSize = psfSize # Can change this later to get it from file

			# Make a new image with simulated galaxies in it
			newImage = art.mkobject(outfile,
									objects = gallist_out+'.modx',
									yoffset = -1*(split*si),
									star = psfFile,
									radius = psfSize,
									magzero = inputImageZP,
									exptime = inputImageExp,
									comments = False,
									gain = inputImageGain,
									poisson = poisson)

		# Now we want to 
		inputfiles = '\n'.join([filename for filename in outfilelist]) + '\n'
		inputstring = outputPath+'irafinput.txt'
		txt_file = open(inputstring, 'w')
		txt_file.write(inputfiles)
		txt_file.close()

		if os.path.isfile(outputImage):
			os.remove(outputImage)
		sub_stack = images.imcombine(input="@"+inputstring, output = outputImage,
									combine = "median", offset="wcs", Stdout=1)
		print('Simulated image {0} saved to {1}'.format(s, outputImage))

		# Clean up
		os.remove(inputstring)
		os.remove(gallist_out+'.modx')
		for si in range(int(splitRatio)):
			os.remove(outputImage+'.split{0}.fits'.format(si))

if doSex:
	# Now we want to run SExtractor on the new images
	print 'Running SExtractor on mock images...'
	for s in range(numSims):
		outputImage = outputPath + outPrefix + str(s+1) + '.fits'
		cmd = '{0} {1} -c {2}'.format(sexCmd, outputImage, sexConfig)
		cmdOutput = subprocess.call(cmd, stdout = open('sextractor.stdout', 'w'),
									stderr = open('sextractor.stderr', 'w'),
									shell = True)

		if cmdOutput is not 0:
			print 'SExtractor return code not 0. Check the terminal and error logs.'
			sys.exit()

		# Copy and rename the output file
		shutil.move(sexOutFile, outputPath+outPrefix+'{0}.sexed'.format(s))
		print 'Run SExtractor on mock image {0}'.format(s)


if doCat:
	print 'Matching catalogues with TOPCAT...'
	tmpCatPath = outputPath+'matched_tmp.fits'
	masterMatches = Table()
	for s in range(numSims):
		# Now match with TOPCAT
		cmd = '{0} tmatch2 in1={1} ifmt1=fits in2={2} ifmt2=ascii join=all2\
					matcher=2d values1="X_IMAGE Y_IMAGE" values2="x y" out={3} \
					params={4:1.1f} omode=out'.format(	topCmd, 
										outputPath+outPrefix+'{0}.sexed'.format(s),
										outputPath+'simulated_galaxies.{0}.tab'.format(s),
										tmpCatPath,
										topMaxSep)

		cmdOutput = subprocess.call(cmd, stdout = open('topcat.stdout', 'w'),
									stderr = open('topcat.stderr', 'w'),
									shell = True)

		if cmdOutput is not 0:
			print 'TOPCAT return code is {0}. Check the terminal and error logs.'.format(cmdOutput)
			sys.exit()

		# Read in the matched catalogue
		matches = Table.read(tmpCatPath)
		masterMatches = vstack([masterMatches, matches])

	# Save the master matched catalogue
	if os.path.isfile(topSave):
		os.remove(topSave)
	masterMatches.write(topSave)
	print 'All matches saved to {0}'.format(topSave)
		
if doStats:
	# Load in the master table
	print 'Performing completeness calculations...'
	masterMatches = Table.read(topSave)
	
	# Set the conditions mask
	magDiff = np.abs(masterMatches['MAG_AUTO'] - masterMatches['mag'])
	masterMask = ((masterMatches['Separation'] >= 0.) * (masterMatches['Separation'] <= statRMax)
				* (magDiff <= statMagMax))
	
	dm = np.clip(np.diff(statMagBins)[0] / 2., 0.1, 0.5)
	binResults = []
	realResults = []
	badFrac = []
	for mag in statMagBins:
		binMask = ((masterMatches['mag'] >= (mag-dm)) * (masterMatches['mag'] < (mag+dm)))
		realResults.append(binMask.sum())
		binObjs = (binMask * masterMask)
		binResults.append(binObjs.sum())
		bad = float(np.invert(np.isfinite(masterMatches['MAG_AUTO'][binMask])).sum())
		badFrac.append( bad/binMask.sum(dtype=float) )

	# Set up PDF for plots
	with PdfPages('look_at_these_plots.pdf') as pdf:
		# Plot the completeness
		Fig, ax = plt.subplots(1,1,figsize=(6,4))
		ax.plot(statMagBins, np.array(binResults, dtype=float)/np.array(realResults, dtype=float),
			lw=3, color='royalblue', label='matched')
		plt.plot(statMagBins, 1.-np.array(badFrac), '-k', alpha=0.8, lw=2, label='detections')
		plt.legend(loc='best', fontsize=8).draw_frame(False)
		plt.plot([10,30], [0.95, 0.95], ':k', lw=1)
		plt.plot([10,30], [0.9, 0.9], ':k', lw=1)
		plt.plot([10,30], [0.8, 0.8], ':k', lw=1)
		plt.plot([10,30], [0.5, 0.5], ':k', lw=1)
		plt.xlim(np.min(statMagBins), np.max(statMagBins))
		plt.ylim(-0.05, 1.05)
		plt.xlabel(r'$K_s$ magnitude (AB)')
		plt.ylabel('completeness')
		plt.text(18, 0.3, r'N = ' + str(numSims) + '\ndmag = '+str(statMagMax) +'\nr_max = '+str(statRMax), fontsize=10)
		plt.text(25.5, 0.94, '95%', fontsize=7, backgroundcolor='w')
		plt.text(25.5, 0.89, '90%', fontsize=7, backgroundcolor='w')
		plt.text(25.5, 0.79, '80%', fontsize=7, backgroundcolor='w')
		plt.text(25.5, 0.49, '50%', fontsize=7, backgroundcolor='w')
		majTicks = MultipleLocator(1.)
		minTicks = MultipleLocator(.5)
		ax.xaxis.set_major_locator(majTicks)
		ax.xaxis.set_minor_locator(minTicks)
		plt.tight_layout()
		pdf.savefig()
		plt.close()

		# Plot the magnitude difference
		plt.figure(figsize=(6,4))
		n, bins, patches = plt.hist(masterMatches['MAG_AUTO'] - masterMatches['mag'], bins=np.arange(-3., 3., 0.1), histtype='step',
			 						lw=3, color='royalblue')
		plt.plot([0,0],[0,n.max()], ':k', lw=1.5)
		plt.plot([statMagMax, statMagMax], [0, n.max()], ':k', lw=1.5)
		plt.plot([-statMagMax, -statMagMax], [0, n.max()], ':k', lw=1.5)
		plt.xlabel('(extracted - original) magnitude')
		plt.ylabel('counts')
		plt.tight_layout()
		pdf.savefig()
		plt.close()

		# Plot the magnitude difference as a function of magnitude
		plt.figure(figsize=(6,4))
		plt.plot(masterMatches['mag'], masterMatches['MAG_AUTO'], '.k', ms=1.3, rasterized=True, alpha=0.7)
		plt.plot([0,30],[0,30], '-r', lw=1, alpha=0.7)
		plt.ylim(18, 27), plt.xlim(18,27)
		plt.ylabel('output magnitude')
		plt.xlabel('input magnitude')
		plt.tight_layout()
		pdf.savefig()
		plt.close()

		# Plot the magnitude distribution in and out
		plt.figure(figsize=(6,4))
		plt.hist(masterMatches['mag'], bins=statMagBins, histtype='step', lw=1.5, 
			label='input magnitudes')
		plt.hist(masterMatches['MAG_AUTO'], bins=statMagBins, histtype='step', lw=1.5,
			label='output magnitudes')
		plt.hist(masterMatches['mag'][masterMask], bins=statMagBins, histtype='step', lw=1.5,
			label='matched input mags')
		plt.legend(loc='best', fontsize=9).draw_frame(False)
		plt.xlabel('K-band magnitude'), plt.ylabel('counts')
		plt.tight_layout()
		pdf.savefig()
		plt.close()

		# Plot the separations
		plt.figure(figsize=(6,4))
		plt.plot(masterMatches['mag'], np.log10(masterMatches['Separation']), '.k', ms=1.3, rasterized=True, alpha=0.7)
		# plt.ylim(0,5)
		plt.xlabel('input magnitude'), plt.ylabel('log separation (pixels)')
		plt.tight_layout()
		pdf.savefig()
		plt.close()

		# Plot the separation as a function of mag difference
		plt.figure(figsize=(6,4))
		plt.plot((masterMatches['mag']-masterMatches['MAG_AUTO']), np.log10(masterMatches['Separation']), '.k', ms=1.3, rasterized=True, alpha=0.7)
		plt.xlim(-2,2)
		plt.xlabel('magnitude difference'), plt.ylabel('log separation (pixels)')
		plt.tight_layout()
		pdf.savefig()
		plt.close()
