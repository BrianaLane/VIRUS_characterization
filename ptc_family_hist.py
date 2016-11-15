import pyfits 
from numpy import * 
import numpy.ma as ma
from pylab import * 
from string import *
from astropy.io import fits
import glob 
import scipy.stats as stats
import matplotlib.mlab as mlab

#------------------------#
# User Defined Variables #
#------------------------#

#Chip to work on (LL,LU,RL,RU)
chip = 'RU'
#Camera number 
cam = '501'
#Bin size for steps in signal level
bin_size = 200
#Number of points to cut on shot noise curve on beginning and end
shot_Begin = 0
shot_End = 0
#Number of sigma to clip away from shot noise fit
shot_sigma = 1

#--------------------#
# User Paths to Data #
#--------------------#

#path the the ptc data
ptc_path = '/Users/Briana/Desktop/VIRUS/Cam'+cam+'/ptc/exp*/camra/'

# give it a path the the set of biases if you want it to calculate RN
bias_path = '/Users/Briana/Desktop/VIRUS/Cam'+cam+'/bias/exp*/camra/'

#-------------------------#
# Retrive and Define Data #
#-------------------------#

#makes a list of all of the ptc data for the specified chip
all_list = glob.glob(str(ptc_path)+'*'+str(chip)+'*.fits')

#makes a list of all of the odd number exposures 
im_list = all_list[0::2]
#makes a list of all of the correpsonding exposure time even number exposures
twin_list = all_list[1::2]

S_DN_list            = []
total_noise_list     = []
Read_Shot_Noise_list = []

#for getting a regions using header coordinates (returns array region from coordinates)
def get_region(hdr_name,data,header,tag):
    reg = header[hdr_name]
    #print str(tag)+' region aquired: '+str(reg)
    h = split(split(split(reg,'[')[1],']')[0],',')
    hx = split(h[0],':')
    hy = split(h[1],':')   
    region = data[(int(hy[0])-1):(int(hy[1])-1),(int(hx[0])-1):(int(hx[1])-1)]   
    return region

#----------------------------#
# Bin Pixels in Median Image #
#----------------------------#

mask_list = []

#Iterate through all of the images and save the data array to a list
dat_list        = []
dat_subOff_list = []
dat_t_list      = []

for i in range(len(im_list)):
    im  = pyfits.open(im_list[i])
    hdr = im[0].header
    dat = im[0].data 

    #save the region of the chip with data (dat_region), and reference pixels (ref_region)
    dat_region = get_region('TRIMSEC',dat,hdr,'data')
    ref_region = get_region('BIASSEC',dat,hdr,'reference')

    #Find ADC offset and subtract it from the data 
    ADC_off = average(ref_region)
    dat_subOff = subtract(dat_region, ADC_off)

    im_t  = pyfits.open(twin_list[i])
    hdr_t = im_t[0].header
    dat_t = im_t[0].data 

    #Twin region for data analysis 
    twin_region = get_region('TRIMSEC',dat_t,hdr_t,'twin data')

    dat_list.append(dat_region)
    dat_subOff_list.append(dat_subOff)
    dat_t_list.append(twin_region)

#from all of the images stack them and build a median image 
master_image  = vstack(dat_list)
master_subOff = vstack(dat_subOff_list)
master_twin   = vstack(dat_t_list)

#------------------------------#
# Build the Masks for Each Bin #
#------------------------------#

#build a list of the signal bins from the min to max signal in the median image
#use defines size in signal of each bin
#max number is non inclusive so add an extra bin size to includes up to max in array 
bin_space = arange(0, 40000, bin_size)

print str(len(bin_space))+' bins used with a bin size of '+str(bin_size)+' DN'
print '\n'

#build mask for each bin based on the median image
cnt = 0
for b in range(len(bin_space)-1):
    bin_mask = ma.masked_outside(master_subOff, bin_space[b], bin_space[b+1]-1)
    print 'bin avg: ' + str(average(bin_mask.compressed())) + ', # in bin: ' + str(bin_mask.count())
    mask_list.append(bin_mask.mask)

    cnt = cnt + bin_mask.count()

    #-----------------------------#
    # Process Images + Remove FPN #
    #-----------------------------#

    #mask each array with the bin mask 
    image_mask  = ma.masked_array(master_image, mask = bin_mask.mask).compressed()
    twin_mask   = ma.masked_array(master_twin,  mask = bin_mask.mask).compressed()
    subOff_mask = bin_mask.compressed()

    #Find average signal value S(DN) (equation 5.1)
    S_DN = average(subOff_mask)
    S_DN_list.append(S_DN)
    
    #Find Total Noise 
    total_noise = std(subOff_mask)
    total_noise_list.append(total_noise)

    #Remove FPN 
    FPN_diff_im = image_mask - twin_mask 

    #Find Read Noise + Shot Noise (DN)
    Read_Shot_Noise = std(FPN_diff_im)/sqrt(2.0)
    Read_Shot_Noise_list.append(Read_Shot_Noise)

print str(100-(((master_subOff.size - cnt)/float(master_subOff.size))*100)) + " percent of points used"
print '\n'

#------------------------#
# Noise Source Functions #
#------------------------#

#this function finds the read noise by taking the stdv of differenced consecutive bias frames and taking the mean of those stdv's 
#Returns single value for the read noise
def find_RN_bias():

    bias_chip_list = glob.glob(str(bias_path)+'*'+str(chip)+'*.fits')

    im1 = pyfits.open(bias_chip_list[0])
    dat1 =  im1[0].data

    stack = dat1
    for i in range(len(bias_chip_list)-1):

        im = pyfits.open(bias_chip_list[i+1])
        dat = im[0].data

        array1 = dstack((stack,dat))
        stack = array1

    readnoise_im = std(stack, axis=2)
    #print std(readnoise_im),average(readnoise_im), average(readnoise_im)+(1*std(readnoise_im)),average(readnoise_im)-(1*std(readnoise_im))
    rdnoiseclipped   = stats.sigmaclip(readnoise_im)[0]
    #print stats.sigmaclip(readnoise_im)[1], stats.sigmaclip(readnoise_im)[2]

    # the histogram of the data
    n, bins, patches = plt.hist(rdnoiseclipped, 50, facecolor='red', alpha=0.75)
    #plt.hist(rdnoiseclipped, 50, normed=1, facecolor='blue', alpha=0.75)
    # add a 'best fit' line
    y = mlab.normpdf( bins, average(rdnoiseclipped), std(rdnoiseclipped))
    l = plt.plot(bins, y, 'r--', linewidth=1)

    plt.xlabel('Standard Deviation (DN)')
    plt.ylabel('Number of Pixels')
    plt.title(r'$\mathrm{Histogram\ of\ Standard\ Deviation\ Bias\ Frame}$')
    #plt.axis([0, 6, 0, 1])
    plt.grid(True)

    #plt.show()

    avg_chip_RN = mean(rdnoiseclipped)

    print 'Read Noise from biase frames (DN): '+str(avg_chip_RN)
    print '\n'

    return avg_chip_RN
 
#Uses the equation  5.11 from Janesick to calculate the shot noise given the read and read+shot noise. 
#Returns array of shot noise values   
def find_shot_noise(read_shot_list, read):
    shot_noise_list = []
    for rs in range(len(read_shot_list)):
        shot_noise = sqrt(pow((read_shot_list[rs]),2)-pow(read,2))
        shot_noise_list.append(shot_noise)
    return shot_noise_list

#-----------------#
# Find Read Noise #
#-----------------#

#Find Read Noise from biases 
read_noise_bias = find_RN_bias()

#-----------------#
# Find Shot Noise #
#-----------------#

#Finds Shot Noise List
shot_noise_list = find_shot_noise(Read_Shot_Noise_list, read_noise_bias)
#get rid of nan values
shot_noise_list = nan_to_num(shot_noise_list)

#choose start and end points to exclude saturated points at the end and goofy points at the beginning. 
#defined by user
st = shot_Begin
ed = len(shot_noise_list)-shot_End

S_DN_clipped = S_DN_list[st:ed]
shot_noise_clipped = shot_noise_list[st:ed]

#This is fitting an initial function to the data that will be subtracted to eliminate points far from fit
#fit a line to the log of the shot noise to get a linear fit in log space 
p = polyfit(log10(S_DN_clipped),log10(shot_noise_clipped),1)
Shot_fit_log = polyval(p,log10(S_DN_clipped))
#have to convert the log data back into linear space to be plotted correctly
Shot_fit=pow(10,Shot_fit_log)

#find data points far from initial fit  
data_fit_diff = subtract(array(shot_noise_clipped),array(Shot_fit))
avg_diff = average(data_fit_diff)
std_diff = std(data_fit_diff)

#Sigma clip points far from intial fit (defined by user)
shot_mask = ma.masked_outside(data_fit_diff,avg_diff-(shot_sigma*std_diff), avg_diff+(shot_sigma*std_diff))

#mask sigma clipped points from signal and shot lists
shot_noise_mask1 = ma.masked_array(shot_noise_clipped, mask = shot_mask.mask).compressed()
S_DN_mask1 = ma.masked_array(S_DN_clipped, mask = shot_mask.mask).compressed()

#mask extremely low signal value points in case there are still some left over
S_DN_mask = ma.masked_less(S_DN_mask1,5)

#apply all masks to both signal and shot lists 
shot_noise_masked = ma.masked_array(shot_noise_mask1, mask = S_DN_mask.mask).compressed()
S_DN_masked = ma.masked_array(S_DN_mask1, mask = S_DN_mask.mask).compressed()

#fit a new shot noise curve to the masked data 
p = polyfit(log10(S_DN_masked),log10(shot_noise_masked),1)
Shot_fit_log = polyval(p,log10(S_DN_masked))
#have to convert the log data back into linear space to be plotted correctly 
Shot_fit=pow(10,Shot_fit_log)

#Find the Values of K(ADC) from shot noise 
#find slope of the shot noise curve
#these used only the masked points to caluclate the gain and slope
slope_shot_plot = round(p[0],2)
K_ADC_lis = divide(S_DN_masked,pow(shot_noise_masked,2)) #this is equation 4.20 from Janesick and assumes a quantum_yield=1
K_ADC = average(K_ADC_lis)
print 'K(ADC): '+str(K_ADC)
print 'slope of shot curve: '+str(slope_shot_plot)
print '\n'

#---------------#
# Print Results #
#---------------#

# #just prints a list of the results 
# print str('Image  '), str('Total     '), str('Read_Shot           '), str('Shot              '), str('S(DN) ')
# for i in range(len(total_noise_list)):
#     if i%10:
#         print str(i)+'      '+str(total_noise_list[i])+'      '+str(Read_Shot_Noise_list[i])+'      '+str(shot_noise_list[i])+'      '+str(S_DN_list[i])

#------------------#
# Plot PTC Results #
#------------------#

#Plot PTC Family          
plot(S_DN_list, total_noise_list, marker = '*',mfc = 'None',mec = 'r',mew =5, ls = 'None',label='Total')
loglog(S_DN_list, Read_Shot_Noise_list, marker = 'o', mfc = 'None',mec = 'blue',mew =2,ls='None',label='Read+Shot')

#reference arrays 
x = logspace(0,5,40)
x2 = logspace(-2,5,40)
half_x = sqrt(x2)*10

#plot read noise
read_noise_line = ones(len(x2))*read_noise_bias
loglog(x2,read_noise_line, linewidth = 3, color='indigo', label = 'Read_Noise: '+str(round(read_noise_bias,2))+'(DN)')

#plot shot noise and its fit 
loglog(S_DN_list,shot_noise_list, marker = '^',mfc='None',mec='green',mew=3,ls='None',label='Shot Noise')
loglog(S_DN_masked,shot_noise_masked, marker = '^',mfc='None',mec='black',mew=3,ls='None',label='Shot Noise Mask')
loglog(S_DN_masked,Shot_fit, linewidth = 2, color='darkorange',label = 'Shot Noise Fit: slope '+str(slope_shot_plot))

text(4,1000,'K(ADC)(e-/DN) = '+str(round(K_ADC,2)),color='navy',size='x-large', bbox=dict(facecolor='none', edgecolor='navy', pad=10.0))
text(4,400,'Read Noise = '+str(round(K_ADC*read_noise_bias,2))+' e-',color='indigo',size='x-large', bbox=dict(facecolor='none', edgecolor='indigo', pad=10.0))

grid(which='both')
xlabel('Signal (DN)',fontsize='large')
ylabel('Noise (DN)',fontsize='large')
ylim(1,10000)
xlim(0,100000)
title('PTC Set for Cam'+cam+' '+str(chip),fontsize='x-large')
legend(loc=2)

show()


    