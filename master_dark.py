import numpy as np
import pyfits 
import glob 
import string as s

darkpath = '/Users/Briana/Documents/Grad_School/LRS2/Casey_galCluster/20160909/camra0009920/'
biaspath = '/Users/Briana/Documents/Grad_School/LRS2/Casey_galCluster/20160909/camra0000007/'

amps  = ['LL','LU','RL','RU']

#for getting a regions using header coordinates (returns array region from coordinates)
def get_region(hdr_name,data,header,tag):
	reg = header[hdr_name]
    #print str(tag)+' region aquired: '+str(reg)
	h = s.split(s.split(s.split(reg,'[')[1],']')[0],',')
	hx = s.split(h[0],':')
	hy = s.split(h[1],':')   
	region = data[(int(hy[0])-1):(int(hy[1])-1),(int(hx[0])-1):(int(hx[1])-1)]   
	return region

#iterate through the amps 
for a in amps:
	darks  = glob.glob(darkpath + 'exp*/camra/*' + a + '*.fits')
	biases = glob.glob(biaspath + 'exp*/camra/*' + a + '*.fits')

	#for each amp build a masterbias 
	bias_list = []
	for b in biases:
		im = pyfits.open(b)
		dat = im[0].data
		hdr = im[0].header

		dat_region = get_region('TRIMSEC',dat,hdr,'data')
		ref_region = get_region('BIASSEC',dat,hdr,'reference')

    	#subtract overscan from each bias frames
		ADC_off = np.average(ref_region)
		dat_subOff = np.subtract(dat_region, ADC_off)

		bias_list.append(dat_subOff)

    #stack all of the baises and take a median frame called masterbias for that amp
	amp_biases = np.dstack(bias_list)
	master_bias = np.median(amp_biases,axis=2)

    #for each amp build a masterdark 
	dark_list = []
	for d in darks:
		im = pyfits.open(d)
		dat = im[0].data
		hdr = im[0].header

		dat_region = get_region('TRIMSEC',dat,hdr,'data')
		ref_region = get_region('BIASSEC',dat,hdr,'reference')

    	#subtract overscan from each of the dark frames
		ADC_off = np.average(ref_region)
		dat_subOff = np.subtract(dat_region, ADC_off)

    	#subtract the amp's masterbias from each dark frame 
		dat_subBias = np.subtract(dat_subOff,master_bias)

		dark_list.append(dat_subBias)

    #stack all of the dark frames and take a median frame called the masterdark for that amp
	amp_darks = np.dstack(dark_list)
	master_frame = np.median(amp_darks,axis=2)

	#save the masterdark for that amp to a fits file 
	#save the header of that last dark frame 
	pyfits.writeto(darkpath+'master_dark_'+a+'.fits', data = master_frame, header = hdr, clobber=True)

	print 'AMP: ' + a 
	print 'median: ' +str(np.median(master_frame))
	print 'std   : '+str(np.std(master_frame))
	print '\n'


