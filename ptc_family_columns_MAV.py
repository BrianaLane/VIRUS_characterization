import pyfits 
from numpy import * 
from pylab import * 
from string import *
from astropy.io import fits
import glob 
import scipy.stats as stats
import matplotlib.mlab as mlab
import argparse as ap
import os.path as op

#------------------------#
# User Defined Variables #
#------------------------#

def parse_args(argv=None):
    """Parse the command line arguments
    Parameters
    ----------
    argv : list of string
        arguments to parse; if ``None``, ``sys.argv`` is used
    Returns
    -------
    Namespace
        parsed arguments
    """
    description = "PTC curve"
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.ArgumentDefaultsHelpFormatter)
                        
    parser.add_argument("-ptc", nargs='?', type=str, help='''Location of the ptc files.  For example: 
                        "/work/03946/hetdex/maverick/virus_lab/20160321/camra0000020"''')

    parser.add_argument("-zro", nargs='?', type=str, help='''Location of the bias (or zro) files.  For example: 
                    "/work/03946/hetdex/maverick/virus_lab/20160321/camra0000020"''')

    parser.add_argument("-cam", type=str, help="Camera or SPECID", 
                    default=None)
    
    parser.add_argument("-shot_beg", type=int, help="Point to begin using shot noise values in fit", 
                    default=200)

    parser.add_argument("-shot_end", type=int, help="Point to stop using shot noise values in fit", 
                    default=300)

    parser.add_argument("-size_cols", type=int, 
                        help='''pixel size of columns to divide images by to make K(ADC) measurement''',
                        default=50)

    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.ptc is None:
        msg = 'The ptc folder location was not provided'
        parser.error(msg)
    if args.zro is None:
        msg = 'The zro folder location was not provided'
        parser.error(msg)
    if args.cam is None:
        msg = 'A camera number was not provided'
        parser.error(msg)
    else:
        try:
            args.cam = "%03d" %(int(args.cam))
        except ValueError:
            msg = 'The camera number was not proper format.  Try "008" or "8".'
            parser.error(msg)      
           
    return args

args = parse_args()

#---------------#
# Set Variables #
#---------------#

#Chip to work on (LL,LU,RL,RU)
amp_list = ['LU','RU','LL','RL']
#Camera number 
cam = str(args.cam)
#Column width to use as ptc region 
col_size = args.size_cols
#Number of points to cut on shot noise curve on beginning and end
shot_Begin = args.shot_beg
shot_End = args.shot_end
#Number of sigma to clip away from shot noise fit
shot_sigma = 0.3

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(13,10))
ax_DICT = { 'LL':ax1, 'RL':ax2, 'LU':ax3, 'RU':ax4 }

#------------------------#
# Noise Source Functions #
#------------------------#

#for getting a regions using header coordinates (returns array region from coordinates)
def get_region(hdr_name,data,header,tag):
    reg = header[hdr_name]
    #print str(tag)+' region aquired: '+str(reg)
    h = split(split(split(reg,'[')[1],']')[0],',')
    hx = split(h[0],':')
    hy = split(h[1],':')   
    region = data[(int(hy[0])-1):(int(hy[1])-1),(int(hx[0])-1):(int(hx[1])-1)]   
    return region

#this function finds the read noise by taking the stdv of differenced consecutive bias frames and taking the mean of those stdv's 
#Returns single value for the read noise
def find_RN_bias(bias_chip_list):

    im1 = pyfits.open(bias_chip_list[0])
    dat1 =  im1[0].data

    stack = dat1
    for i in range(len(bias_chip_list)-1):

        im = pyfits.open(bias_chip_list[i+1])
        dat = im[0].data

        array1 = dstack((stack,dat))
        stack = array1

    readnoise_im = std(stack, axis=2)
    rdnoiseclipped   = stats.sigmaclip(readnoise_im)[0]

    avg_chip_RN = mean(rdnoiseclipped)

    print 'Read Noise from biase frames (DN): '+str(avg_chip_RN)

    return avg_chip_RN

#Uses the equation  5.11 from Janesick to calculate the shot noise given the read and read+shot noise. 
#Returns array of shot noise values   
def find_shot_noise(read_shot_list, read):
    shot_noise_list = []
    for rs in range(len(read_shot_list)):
        shot_noise = sqrt(pow((read_shot_list[rs]),2)-pow(read,2))
        shot_noise_list.append(shot_noise)
    return shot_noise_list

#----------------------#
# Iterate through amps #
#----------------------#
for amp in amp_list:
    print '******'
    print '* '+amp+' *'
    print '******'
    #-------------------------#
    # Retrive and Define Data #
    #-------------------------#

    #makes a list of all of the ptc data for the specified chip
    ptc_list = glob.glob(op.join(args.ptc,'exp*/camra/*'+amp+'*.fits' ))
    bias_list = glob.glob(op.join(args.zro,'exp*/camra/*'+amp+'*.fits' ))

    #makes a list of all of the odd number exposures 
    im_list = ptc_list[0::2]
    #makes a list of all of the correpsonding exposure time even number exposures
    twin_list = ptc_list[1::2]

    exp_time_list = []
    S_DN_list = []
    total_noise_list = []
    Read_Shot_Noise_list = []

    #coordinate for ptc region (region without reference pixels)
    x_lis = arange(0,2100,col_size)

    #-----------------------------#
    # Process Images + Remove FPN #
    #-----------------------------#

    #iterates through all of the ptc images   
    for i in range(len(im_list)): 

        #open data and header from .fits file
        im = pyfits.open(im_list[i])
        hdr = im[0].header
        dat =  im[0].data

        #get exposure time from header and add it to the exposure time list     
        exp_time = hdr['EXPTIME']
        exp_time_list.append(int(exp_time))  

        #save the region of the chip with data (dat_region) 
        dat_region = get_region('TRIMSEC',dat,hdr,'data')

        #Save the region of the chip that is reference pixels (ref_region)
        ref_region = get_region('BIASSEC',dat,hdr,'reference')

        #Twin image header and data for FPN subtraction
        #print 'Twin File Name: '+str(twin_list[i])
        twin = pyfits.open(twin_list[i])
        hdr_twin = twin[0].header
        dat_twin =  twin[0].data  

        exp_time1 = hdr_twin['EXPTIME']
        #print 'Exposure Time: '+str(exp_time1)

        ADC_off = average(ref_region)

        #Twin region for data analysis 
        twin_region = get_region('TRIMSEC',dat_twin,hdr_twin,'twin data')

        for j in range(len(x_lis)-1):

            #coordinate for ptc region (region without reference pixels)
            x1 = x_lis[j]
            x2 = x_lis[j+1]
            y1 = 100
            y2 = 950

            #save the region for ptc analysis from dat_region (ptc_region)
            ptc_region = dat_region[y1:y2,x1:x2]
            
            #Find ADC offset and subtract it from the data 
            ptc_subOff = subtract(ptc_region, ADC_off)

            #Find average signal value S(DN) (equation 5.1)
            S_DN = average(ptc_subOff)
            S_DN_list.append(S_DN)
            
            #Find Total Noise 
            total_noise = std(ptc_subOff)
            total_noise_list.append(total_noise)

            twin_ptc_region = twin_region[y1:y2,x1:x2]
         
            #Remove FPN 
            FPN_diff_im = ptc_region - twin_ptc_region 
            
            #Find Read Noise + Shot Noise (DN)
            Read_Shot_Noise = std(FPN_diff_im)/sqrt(2.0)
            Read_Shot_Noise_list.append(Read_Shot_Noise)
            #print 'S: '+str(S_DN)+'  Tot: '+str(total_noise)+'  FP: '+str(average(FPN_diff_im))+'  Read+Shot: '+str(Read_Shot_Noise)

    #-----------------#
    # Find Read Noise #
    #-----------------#

    #Find Read Noise from biases 
    read_noise_bias = find_RN_bias(bias_list)

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

    #mask points that show with signal below 0 in case there are goofy points
    SDN_mask2 = ma.masked_less(S_DN_clipped,0)

    #apply mask to shot and signal lists
    shot_noise_maskz = ma.masked_array(shot_noise_clipped, mask = SDN_mask2.mask).compressed()
    S_DN_maskz = ma.masked_array(S_DN_clipped, mask = SDN_mask2.mask).compressed()

    #This is fitting an initial function to the data that will be subtracted to eliminate points far from fit
    #fit a line to the log of the shot noise to get a linear fit in log space 
    p = polyfit(log10(S_DN_maskz),log10(shot_noise_maskz),1)
    Shot_fit_log = polyval(p,log10(S_DN_maskz))
    #have to convert the log data back into linear space to be plotted correctly
    Shot_fit=pow(10,Shot_fit_log)

    #find data points far from initial fit  
    data_fit_diff = subtract(array(shot_noise_maskz),array(Shot_fit))
    avg_diff = average(data_fit_diff)
    std_diff = std(data_fit_diff)

    #Sigma clip points far from intial fit (defined by user)
    shot_mask = ma.masked_outside(data_fit_diff,avg_diff-(shot_sigma*std_diff), avg_diff+(shot_sigma*std_diff))

    #mask sigma clipped points from signal and shot lists
    shot_noise_masked1 = ma.masked_array(shot_noise_maskz, mask = shot_mask.mask).compressed()
    S_DN_masked1 = ma.masked_array(S_DN_maskz, mask = shot_mask.mask).compressed()

    #mask extremely low signal value points in case there are still some left over
    S_DN_mask = ma.masked_less(S_DN_masked1,5)

    #apply all masks to both signal and shot lists 
    shot_noise_masked = ma.masked_array(shot_noise_masked1, mask = S_DN_mask.mask).compressed()
    S_DN_masked = ma.masked_array(S_DN_masked1, mask = S_DN_mask.mask).compressed()

    #fit a new shot noise curve to the masked data 
    p = polyfit(log10(S_DN_masked),log10(shot_noise_masked),1)
    Shot_fit_log = polyval(p,log10(S_DN_masked))
    #have to convert the log data back into linear space to be plotted correctly 
    Shot_fit=pow(10,Shot_fit_log)

    #----------------------------------#
    # Find Gain Given Shot Noise Curve #
    #----------------------------------#

    #Find the Values of K(ADC) from shot noise 
    #find slope of the shot noise curve
    #these used only the masked points to caluclate the gain and slope
    slope_shot_plot = round(p[0],2)
    K_ADC_lis = divide(S_DN_masked,pow(shot_noise_masked,2)) #this is equation 4.20 from Janesick and assumes a quantum_yield=1
    K_ADC = average(K_ADC_lis)
    print 'K(ADC)(e-/DN) : '+str(K_ADC)
    print 'Std of K(ADC): '+str(std(K_ADC_lis))
    print 'slope of shot curve: '+str(slope_shot_plot)
    print '\n'

    #---------------#
    # Print Results #
    #---------------#

    #just prints a list of the results 
    #print str('Image  '), str('Total     '), str('Read_Shot           '), str('Shot              '), str('S(DN) ')
    #for i in range(len(total_noise_list)):
    #    print str(i)+'      '+str(total_noise_list[i])+'      '+str(Read_Shot_Noise_list[i])+'      '+str(shot_noise_list[i])+'      '+str(S_DN_list[i])

    #------------------#
    # Plot PTC Results #
    #------------------#

    #Plot PTC Family          
    ax_DICT[amp].plot(S_DN_list, total_noise_list, marker = '*',mfc = 'None',mec = 'r',mew =5, ls = 'None',label='Total')
    ax_DICT[amp].loglog(S_DN_list, Read_Shot_Noise_list, marker = 'o', mfc = 'None',mec = 'blue',mew =2,ls='None',label='Read+Shot')

    #reference arrays 
    x = logspace(0,5,40)
    x2 = logspace(-2,5,40)
    half_x = sqrt(x2)*10

    #plot read noise
    read_noise_line = ones(len(x2))*read_noise_bias
    ax_DICT[amp].loglog(x2,read_noise_line, linewidth = 3, color='indigo', label = 'Read_Noise: '+str(round(read_noise_bias,2))+'(DN)')

    #plot shot noise and its fit 
    ax_DICT[amp].loglog(S_DN_list,shot_noise_list, marker = '^',mfc='None',mec='green',mew=3,ls='None',label='Shot Noise')
    ax_DICT[amp].loglog(S_DN_masked,shot_noise_masked, marker = '^',mfc='None',mec='black',mew=3,ls='None',label='Shot Noise Mask')
    ax_DICT[amp].loglog(S_DN_masked,Shot_fit, linewidth = 2, color='darkorange',label = 'Shot Noise Fit: slope '+str(slope_shot_plot))

    ax_DICT[amp].text(0.05,2500,'K(ADC)(e-/DN) = '+str(round(K_ADC,2)),color='navy',size='medium', bbox=dict(facecolor='none', edgecolor='navy', pad=10.0))
    ax_DICT[amp].text(0.05,400,'Read Noise = '+str(round(K_ADC*read_noise_bias,2))+' e-',color='indigo',size='medium', bbox=dict(facecolor='none', edgecolor='indigo', pad=10.0))

    ax_DICT[amp].grid(which='both')
    ax_DICT[amp].set_ylim(1,10000)
    ax_DICT[amp].set_xlim(0,100000)
    ax_DICT[amp].set_title('Amp '+amp,fontsize='x-large')

    box = ax_DICT[amp].get_position()
    ax_DICT[amp].set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

labels = [line.get_label() for line in ax1.lines]
f.legend(ax1.lines, labels,  loc = 'lower right', ncol=3, labelspacing=0., borderaxespad=2, bbox_to_anchor=(0.8, -0.03))
#f.legend(ax1.lines, labels, loc = 'upper right', bbox_to_anchor = (0,-0.1,1,1), bbox_transform = plt.gcf().transFigure )

f.text(0.5, 0.08, 'Signal (DN)',fontsize='large', ha='center')
f.text(0.04, 0.5, 'Noise (DN)',fontsize='large', va='center', rotation='vertical')
f.suptitle('PTC Set for Cam'+cam,fontsize='x-large')
f.savefig(op.join(args.ptc,'..','PTC_Cam'+cam+'.png'))
