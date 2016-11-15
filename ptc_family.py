import pyfits 
from numpy import * 
from pylab import * 
from string import *
import glob 
import scipy.stats as stats

#Choose the chip to work on (LL,LU,RL,RU)
chip = 'LL'
cam = '31'

#coordinate for ptc region (region without reference pixels)
x1 = 300
x2 = 1700
y1 = 100
y2 = 1000

#path the the ptc data
#ptc_path = '/Volumes/USB_BRI/Charactization/Cam'+cam+'/ptc/exp*/camra/'
ptc_path = '/Users/Briana/Desktop/VIRUS/Cam'+cam+'/ptc/exp*/camra/'
#ptc_path = '/Users/Briana/Documents/Grad_School/LRS2/LRS2_red/lrs20002015/exp*/lrs2/'

# give it a path the the set of biases if you want it to calculate RN
#bias_path = '/Volumes/USB_BRI/Charactization/Cam'+cam+'/bias/exp*/camra/'
bias_path = '/Users/Briana/Desktop/VIRUS/Cam'+cam+'/bias/exp*/camra/'
#bias_path = '/Users/Briana/Documents/Grad_School/LRS2/LRS2_red/lrs20002010/exp*/lrs2'

#makes a list of all of the ptc data for the specified chip
all_list = glob.glob(str(ptc_path)+'*'+str(chip)+'*.fits')

#makes a list of all of the odd number exposures 
im_list = all_list[0::2]
#makes a list of all of the correpsonding exposure time even number exposures
twin_list = all_list[1::2]

exp_time_list = []
S_DN_list = []
total_noise_list = []
Read_Shot_Noise_list = []

#for getting a regions using header coordinates (returns array region from coordinates)
def get_region(hdr_name,data,header,tag):
    reg = header[hdr_name]
    print str(tag)+' region aquired: '+str(reg)
    h = split(split(split(reg,'[')[1],']')[0],',')
    hx = split(h[0],':')
    hy = split(h[1],':')   
    region = data[(int(hy[0])-1):(int(hy[1])-1),(int(hx[0])-1):(int(hx[1])-1)]   
    return region

def find_ADC_offset(overscan_reg):
    return overscan_reg

#iterates through all of the ptc images   
for i in range(len(im_list)):
    
    print 'Image: '+str(i+1)
    print 'File Name: '+str(im_list[i])
    
    #open data and header from .fits file
    im = pyfits.open(im_list[i])
    hdr = im[0].header
    dat =  im[0].data
    
    #get exposure time from header and add it to the exposure time list     
    exp_time = hdr['EXPTIME']
    exp_time_list.append(int(exp_time))
    print 'Exposure Time: '+str(exp_time)    
    
    #save the region of the chip with data (dat_region) 
    dat_region = get_region('TRIMSEC',dat,hdr,'data')
    print 'Average dat: '+str(average(dat_region))
    print ' Shape: '+str(dat_region.shape)

    #Save the region of the chip that is reference pixels (ref_region)
    ref_region = get_region('BIASSEC',dat,hdr,'reference')
    print 'Average ref: '+str(average(ref_region))
    print ' Shape: '+str(ref_region.shape)
    
    #save the region for ptc analysis from dat_region (ptc_region)
    ptc_region = dat_region[y1:y2,x1:x2]
    print 'Average ptc: '+str(average(ptc_region))
    print ' Shape: '+str(ptc_region.shape)
    
    #Find ADC offset and subtract it from the data 
    #avg_ref = average(ref_region, axis=1)
    #x = arange(1,len(avg_ref)+1,1)
    #funct = polyfit(x,avg_ref,3)
    #ref_fit = polyval(funct,x)
    #fit_M = matrix(ref_fit).T
    #ADC_off_reg = subtract(ptc_region,fit_M)
    #ADC_off = average(fit_M)
    #print 'ADC_off: '+str(ADC_off)
    ADC_off = average(ref_region)
    ptc_subOff = subtract(ptc_region, ADC_off)
    print 'ADC_off: '+str(ADC_off)

    #Find average signal value S(DN) (equation 5.1)
    S_DN = average(ptc_subOff)
    S_DN_list.append(S_DN)
    print 'Average Signal Value S(DN): '+str(S_DN)
    
    #Find Total Noise 
    total_noise = std(ptc_subOff)
    total_noise_list.append(total_noise)
    print 'Total Noise: '+str(total_noise)
    
    #Twin image header and data for FPN subtraction
    print 'Twin File Name: '+str(twin_list[i])
    twin = pyfits.open(twin_list[i])
    hdr_twin = twin[0].header
    dat_twin =  twin[0].data  

    exp_time1 = hdr_twin['EXPTIME']
    print 'Exposure Time: '+str(exp_time1)
    print shape(dat_twin)
    
    #Twin region for data analysis 
    #DETSEC
    twin_region = get_region('TRIMSEC',dat_twin,hdr_twin,'twin data')
    print shape(twin_region)
    twin_ptc_region = twin_region[y1:y2,x1:x2]
    print shape(twin_ptc_region)
    print 'Twin Average Signal: '+str(average(twin_ptc_region))
 
    #Remove FPN 
    FPN_diff_im = ptc_region - twin_ptc_region 
    print 'Average of FPN difference image: '+str(average(FPN_diff_im))
    
    #Find Read Noise + Shot Noise (DN)
    Read_Shot_Noise = std(FPN_diff_im)/sqrt(2.0)
    Read_Shot_Noise_list.append(Read_Shot_Noise)
    print 'Read+Shot Noise: '+str(Read_Shot_Noise)
    
    print ''

#finds the read noise from the total_noise_list. have to give it the start and end values to average over. 
#Returns single value for the read noise 
def find_read_noise(avg_sig_list,st,end):
    read_noise = average(avg_sig_list[st:end+1])
    print "Read Noise (DN): "+str(read_noise)
    return read_noise

#this function finds the read noise by taking the stdv of differenced consecutive bias frames and taking the mean of those stdv's 
#Returns single value for the read noise
def find_RN_bias_sub():

    bias_chip_list = glob.glob(str(bias_path)+'*'+str(chip)+'*.fits')

    RN_list = []
    for i in range(len(bias_chip_list)/2):

        #image in odd list
        im1 = pyfits.open(bias_chip_list[0::2][i])
        dat1 =  im1[0].data

        #image in even list
        im2 = pyfits.open(bias_chip_list[1::2][i])
        dat2 =  im2[0].data

        #subtracts consecutive images
        sub_im = subtract(dat1,dat2)

        avg = mean(sub_im)
        variance = var(sub_im)
        sigma = sqrt(variance)
        peak = avg
        sig_clip_h = peak+(4*sigma)
        sig_clip_l = peak-(4*sigma)
        sub_im_clipped = ma.masked_outside(sub_im, sig_clip_l,sig_clip_h)

        avg_clip = mean(sub_im_clipped)
        std_sub = sqrt(average(pow(subtract(sub_im_clipped,avg_clip),2)))
        RN_diff = std_sub/sqrt(2)
        RN_list.append(RN_diff)

    avg_chip_RN = average(RN_list)
    print 'Read Noise from biase frames (DN): '+str(avg_chip_RN)
    print 'Error on RN from bias_frames: '+str(std(RN_list))
    print '\n'

    return avg_chip_RN

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
    rdnoiseclipped   = stats.sigmaclip (readnoise_im)[0]
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

def find_FPN(read_shot,tot_noise):
    rs_sq = square(read_shot)
    tot_noise_sq = square(tot_noise)
    FPN_list = sqrt(tot_noise_sq - rs_sq)    
    return FPN_list

#This funciton finds the full well by finding the first value when the signal is less than the previous value
#returns the index of the full well in the signal list if there is a drop/full well detected
#returns 0 if full well is not found or values continuously increase
def find_full_well(sig_list):
    for i in range(len(sig_list)-1):
        if float(sig_list[i+1]-sig_list[i]) < 0:
            return i  

#Find FPN 
#FPN_list = find_FPN(Read_Shot_Noise_list,total_noise_list)

#Find Read Noise 
read_noise = find_read_noise(total_noise_list,0,0)

#Find Read Noise from biases 
read_noise_bias = find_RN_bias()

#Finds Shot Noise List
shot_noise_list = find_shot_noise(Read_Shot_Noise_list, read_noise_bias)
shot_noise_list = nan_to_num(shot_noise_list)

#choose start and end points to exclude saturated points at the end and goofy points at the beginning. 
st = 4
ed = len(shot_noise_list)-3
#fit a line to the log of the shot noise to get a linear fit in log space 
p = polyfit(log10(S_DN_list)[st:ed],log10(shot_noise_list)[st:ed],1)
Shot_fit_log = polyval(p,log10(S_DN_list))
Shot_fit=[]
#have to convert the log data back into linear space to be plotted correctly 
for i in range(len(Shot_fit_log)):
    Shot_fit_val = pow(10,Shot_fit_log[i])
    Shot_fit.append(Shot_fit_val)

#Find the Values of K(ADC) from shot noise 
#find slope of the shot noise curve
slope_shot_plot = round(p[0],2)
K_ADC_lis = divide(S_DN_list[st:ed],pow(shot_noise_list[st:ed],2)) #this is equation 4.20 from Janesick and assumes a quantum_yield=1
K_ADC = average(K_ADC_lis)
print 'K(ADC): '+str(K_ADC)
print 'slope from shot curve: '+str(slope_shot_plot)

#just prints a list of the results 
print str('Image  '), str('Total     '), str('Read_Shot           '), str('Shot              '), str('S(DN) '), str( 'Exp_Time')
for i in range(len(total_noise_list)):
    print str(i)+'      '+str(total_noise_list[i])+'      '+str(Read_Shot_Noise_list[i])+'      '+str(shot_noise_list[i])+'      '+str(S_DN_list[i])+'      '+str(exp_time_list[i])

#Find Full Well using the find_full_well function and print the result
FW = find_full_well(total_noise_list)
if FW == None:
    print 'Full Well Not Reached!'
else:
    print 'Full Well: '+str(S_DN_list[FW])+' (DN)'


#Plot PTC Family          
plot(S_DN_list, total_noise_list, marker = '*',mfc = 'None',mec = 'r',mew =5, ls = 'None',label='Total')
loglog(S_DN_list, Read_Shot_Noise_list, marker = 'o', mfc = 'None',mec = 'blue',mew =2,ls='None',label='Read+Shot')
#loglog(S_DN_list,FPN_list, marker = 'd',mfc='None',mec = 'g', mew =2,ls = 'None',label='FPN')
#reference arrays 
x = logspace(0,5,40)
x2 = logspace(-2,5,40)
half_x = sqrt(x2)*10
#plot read noise
read_noise_line = ones(len(x2))*read_noise_bias
loglog(x2,read_noise_line, linewidth =3, color ='indigo', label = 'Read_Noise: '+str(round(read_noise_bias,2))+'(DN)')
#plot shot noise and its fit 
loglog(S_DN_list,shot_noise_list, marker = '^',mfc='None',mec='green',mew=3,ls='None',label='Shot Noise')
loglog(S_DN_list,Shot_fit, linewidth=2, color='darkorange',label = 'Shot Noise Fit: slope '+str(slope_shot_plot))
#plot reference lines
#loglog(x,x,label = 'Slope=1 Line')
#loglog(x2,half_x,label = 'Slope=1/2 Line')
text(4,1000,'K(ADC)(e-/DN) = '+str(round(K_ADC,3)),color='navy',size='x-large', bbox=dict(facecolor='none', edgecolor='navy', pad=10.0))
text(4,300,'Read Noise = '+str(round(K_ADC*read_noise_bias,2))+' e-',color='indigo',size='x-large', bbox=dict(facecolor='none', edgecolor='indigo', pad=10.0))

grid(which='both')
xlabel('Signal (DN)',fontsize='large')
ylabel('Noise (DN)',fontsize='large')
ylim(1,10000)
xlim(0,150000)
title('PTC Set for Cam'+cam+' '+str(chip),fontsize='x-large')
legend(loc=2)

show()

#Variance PTC
#var_list = square(Read_Shot_Noise_list)

#v = polyfit(S_DN_list,var_list,1)
#var_fit = polyval(v,S_DN_list)

#print 'Variance Read Noise: '+str(v[1])
#print 'Variance K(ADC): '+str(v[0])

#plot(S_DN_list,var_list, marker = 'x',ls='None')
#plot(S_DN_list,var_fit)

#xlabel('Signal (DN)')
#ylabel('Noise (DN)')
#title('PTC Variance for '+str(chip))

#show()


    