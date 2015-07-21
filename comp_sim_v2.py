# -*- coding: utf-8 -*-

# For VIDEO, we have PSFs as a function of position.
# Go through each galaxy separately and add it to the image.

import os, numpy as np, glob, shutil, sys, subprocess
from astropy.table import Table, vstack
from astropy.io import fits
from pyraf.iraf import noao, images
from pyraf.iraf import artdata as art
import matplotlib.pyplot as plt

###########################################
# QUICK FLAGS
doMocks			= True
doSex			= True
doCat			= True
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
numSims 		= 1
gallist_out 	= 'test.txt'
outPrefix 		= 'video_'
ngals 			= 7000
###########################################
# SIMULATED GALAXY DISTRIBUTIONS
sersic_dist 	= 0.5*np.round(np.random.lognormal(0.4, 0.3, size=ngals)*2.) # crop the distribution before science
axis_dist 		= np.clip(np.random.normal(0.7, 0.2, size=ngals), 0.1, 1.)
mag_dist 		= np.clip(np.random.lognormal(3.1, 0.1, size=ngals), 17., 27.)
pa_dist			= np.random.rand(ngals)*360.
radii_dist 		= np.clip(np.random.lognormal(2., 0.2, size=ngals), 2.5, 500.)
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
statMagBins		= np.arange(17, 27.2, 0.2)
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
		ra = np.random.randint(50, NAXIS1-50, ngals)
		dec = np.random.randint(50, NAXIS2-50, ngals)

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
			psfImage = os.path.split(psfTable['file'][5000])
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
									poisson = True)

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
	for mag in statMagBins:
		binMask = ((masterMatches['mag'] >= (mag-dm)) * (masterMatches['mag'] < (mag+dm)))
		realResults.append(binMask.sum())
		binObjs = (binMask * masterMask)
		binResults.append(binObjs.sum())


	plt.figure(figsize=(6,4))
	plt.plot(statMagBins, np.array(binResults, dtype=float)/np.array(realResults, dtype=float),
		lw=3, color='royalblue')
	plt.plot([10,30], [0.95, 0.95], ':k', lw=1)
	plt.plot([10,30], [0.9, 0.9], ':k', lw=1)
	plt.plot([10,30], [0.8, 0.8], ':k', lw=1)
	plt.plot([10,30], [0.5, 0.5], ':k', lw=1)
	plt.xlim(np.min(statMagBins), np.max(statMagBins))
	plt.xlabel(r'$K_s$ magnitude (AB)')
	plt.ylabel('completeness')
	plt.tight_layout()
	plt.show()