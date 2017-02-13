# -*- coding: utf-8 -*-
"""
Created on Friday April  8 02:10:02 2016

Reduces LRS2 data for either the blue or red channels
This script reads user parameters from lrs2_config.py 

@author: gregz, brianaindahl
"""
from __future__ import print_function

import numpy as np
from scipy import interpolate
import pyfits
import glob
from datetime import datetime
import shutil
import sys
import os
import copy
import os.path as op
from os import environ
import re
import string 
import textwrap
import argparse as ap

########################
# Parse User Arguments #
########################

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
    description = textwrap.dedent('''Science Reduction Script - 
    
                     This script produces the following output:
                         Fe*.fits for fiber extracted spectra
                     
                     The script places the files in the given output directory.
                     
                     ''')
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.RawTextHelpFormatter)
            
    parser.add_argument("--cam", nargs='?', type=str, 
            help='''SPECID or camra number [REQUIRED] 
            Ex: "020".''', default = None)

    parser.add_argument("--output", nargs='?', type=str, 
            help='''Output Directory [REQUIRED] 
            Ex: "/Users/Documents/shela_z4"''', default=None)                       

    parser.add_argument("--date_fold", nargs='?', type=str, 
            help='''Date folder in root directory containing data [REQUIRED] 
            Ex: 20160909''', default=None)

    parser.add_argument("--flt", nargs='?', type=str, 
            help='''Directory containing LDLS flats [REQUIRED] 
            Ex. camra0000020''', default=None)

    parser.add_argument("--cmp", nargs='?', type=str, 
            help='''Directory containing HgCd comps [REQUIRED] 
            Ex. camra0000021''', default=None)

    parser.add_argument("--zro", nargs='?', type=str, 
            help='''Directory containing biases [REQUIRED] 
            Ex. camra0000020''', default=None)

    parser.add_argument('-step', help='''Choose what reduction step to start with 
            Default: basic''', dest='step', choices=['basic','deformer','fibextract'])
    parser.set_defaults(step='basic')

    parser.add_argument('-dont_clean', help='''Will NOT remove intermediate files 
            Default: False''', dest='dont_clean', action='store_true')
    parser.set_defaults(dont_clean=False)
                        
    parser.add_argument("--rootdir", nargs='?', type=str, 
            help='''Root Directory
            Default: \"/work/03946/hetdex/maverick\"''', 
            default="/work/03946/hetdex/maverick")

    parser.add_argument("--configdir", nargs='?', type=str, 
            help='''Config Directory
            Default: \"/work/03946/hetdex/maverick\"''', 
            default="/work/03946/hetdex/maverick")
                          
    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.cam is None:
        msg = 'A camera number was not provided'
        parser.error(msg)
    else:
        try:
            args.cam = "%03d" %(int(args.cam))
        except ValueError:
            msg = 'The camera number was not proper format.  Try "008" or "8".'
            parser.error(msg)      
        
    if args.output is None:
        msg = 'No output directory was provided'
        parser.error(msg)    
    if args.date_fold is None:
        msg = 'No date folder was provided'
        parser.error(msg) 
    if args.flt is None:
        msg = 'No flat folder was provided'
        parser.error(msg) 
    if args.cmp is None:
        msg = 'No comp folder was provided'
        parser.error(msg)
    if args.zro is None:
        msg = 'No bias folder was provided'
        parser.error(msg)   
                
    return args  

#########################################
# Read User Arguments and Set Variables #
#########################################

args = parse_args()

ucam = args.cam
LAMP_DICT = {0:'HgCd'}

#path and name of the folder reduction is run - folder created by script
redux_dir       = args.output    
#path to date folder containing raw data 
date_folder     = op.join(args.rootdir,args.date_fold)        
#path to virus config folder
configdir       = args.configdir 

if   args.step == 'deformer':
    basic           = False      
    run_deformer    = True      
    fiberextract    = True
elif args.step == 'fibextract':
    basic           = False      
    run_deformer    = False      
    fiberextract    = True
else:
    basic           = True      #run basic reduction (overscan + bias subtract, ccd combine, build mastertrace + masterarc)
    run_deformer    = True      #run deformer to map spectral traces and build wavelength solution for fiber extraction
    fiberextract    = True      #extract spectra and save into fits file - Need to have run deformer

if dont_clean:
    CLEAN_AFTER_DONE = False
else:
    CLEAN_AFTER_DONE = True      #If true it will delete intermediate reduction files for the calibration data  

##################################################
# Setting CUREBIN and check LRS2 defined in CURE #
##################################################

#setting CUREBIN
CUREBIN = None
if not CUREBIN:
    CUREBIN = environ.get('CUREBIN')
if not CUREBIN:
    sys.exit("Please set CUREBIN as environment variable")

#checking that LRS2 is defined in specconf.h 
cureversion = os.popen(op.join(CUREBIN, 'cureversion')).readlines()
spec_define = cureversion[4].split(' ')[1]
instrument = spec_define.rstrip('\n')

if instrument == 'VIRUS_HET':
    print ('CURE is set for VIRUS reduction')
else: 
    print ('You need to update specconf.h in CURE to define VIRUS_HET')
    sys.exit('Right now specconf.h defines '+instrument)

###################################
# Defining which functions to run #
###################################

#if basic reduction is run need to specify specific routines to run 
# divide pixel flat and masterdark are not being used now
if basic:
    masterarc       = True  
    mastertrace     = True 
else:
    masterarc       = False  
    mastertrace     = False

# This makes sure that the redux folder is only overwritten if the user chooses to run basic reduction
# If you user only wants to run deformer, skysubtract, fiberextract, or mkcube it used the data in redux 
if basic:
    all_copy = True
    RESTART_FROM_SCRATCH = True
else:
    all_copy = False
    RESTART_FROM_SCRATCH = False

#sets the initial action base to start with. CURE adds bases to data when a change is made 
initial_base     = ''

########################
# Specifying CURE opts #
########################

#specify opts for CURE routines used in this script
meanfitsopts    = "--new_error -k 3"
headfitsopts    = "-m -k EXPTIME -v 1 -V"
arcopts         = "--maxmem 1024 -s -t -m -k 2.8"
traceopts       = "--maxmem 1024 -s -t -m -k 2.8"
deformeropts    = "-s 6 -C 10 --debug --dump_psf_data"
fibextractopts  = "-P"

#########################
# Defining data folders #
#########################

#reads folder names from config file
zro_file_loc = [op.join(date_folder,args.zro)]
flt_file_loc = [op.join(date_folder,args.flt)]
cmp_file_loc = [op.join(date_folder,args.cmp)]


#specify folders where data types are stored
zro_dir  = "zro"
flt_dir  = "flt"
cmp_dir  = "cmp"

##########################
# Building Spec Libaries #
##########################

#dictionary of data type folders 
file_loc_dir = [ zro_file_loc, cmp_file_loc, flt_file_loc ] # Should match order of dictionary, otherwise there will be mismatches
DIR_DICT     = {    0:zro_dir,  1:cmp_dir,    2:flt_dir } # Dictionary for looping through directory names

#Detector Amps and spectrograph side lists 
SPEC = ["LL","LU","RL","RU"]
SPECBIG = ["L","R"]

#############################
# Define config directories #
#############################

#specifies directories for lines and mapping/cen files for LRS2 in the LRS2 config directory 
linesdir    = configdir + '/lines_files'
mapdir      = configdir + '/mapping_files/'

##################
# Define Classes #
##################

'''VirusFrame class.  This is critical to call cure commands from 
"call_cure.py", however this is a local definition due to the ever changing
nature of data to be read in.  Some have header key words set, others do not.
When finished, this will be a global class for all reductions, but for now,
I will just use it here, locally.'''
class VirusFrame:
    def __init__ ( self, filename = None, dir_dict = None):
        '''
        Initializing a VirusFrame for a given filename.
        This includes reading the header and setting reduction keywords.
        From this, the file with have attributes that can be tracked.
        '''
        
        if filename is None:
            print("No filename was given for VirusFrame to be initialized.")
            return None
        if filename is None:
            print("No file type was given for VirusFrame to be initialized.")
            return None
        else:
            ######## OPEN AND KEEP FILES IN SOME SMART WAY ########
            self.filename             = filename
            self.origloc              = op.dirname(self.filename)
            self.basename, tmp1, tmp2 = op.basename(self.filename ).split('_')
            self.type                 = dir_dict
            self.ifuslot              = tmp1[0:3]
            self.time                 = self.basename.split('T')[1]
            self.hr                   = self.basename.split('T')[1][0:2]
            self.min                  = self.basename.split('T')[1][2:4]
            self.sec                  = self.basename.split('T')[1][4:8]
            self.clean                = CLEAN_AFTER_DONE

            
            ###### READ IN THE HEADER AND NECESSARY KEYWORDS ######
            self.trimsec    = {}
            self.biassec    = {}
            self.actionbase = {}
            for amp in SPEC:    
                rootname          = op.join(self.origloc, self.basename + '_' + self.ifuslot + amp + '_' + self.type + '.fits' )
                hdulist = pyfits.open(rootname)
                trimstr = re.split('[\[ \] ]',hdulist[0].header['TRIMSEC'])[1]
                biasstr = re.split('[\[ \] ]',hdulist[0].header['BIASSEC'])[1]
                self.trimsec[amp] = "\"" + trimstr + "\""
                self.biassec[amp] = "\"" + biasstr + "\""
                self.actionbase[amp] = ''            

            self.actionbase["L"] = initial_base  
            self.actionbase["R"] = initial_base  
            
            # Use mapping because if headers don't include SPECID
            if usemapping:
                self.specid = IFUSLOT_DICT[self.ifuslot][0]
            else:
                self.specid = hdulist[0].header['SPECID']
                
            self.object      = hdulist[0].header['OBJECT']
            self.orggain     = hdulist[0].header['GAIN']
            self.orgrdnoise  = hdulist[0].header['RDNOISE']
            self.exptime     = hdulist[0].header['EXPTIME']
                    


            #print (self.type, self.object, self.cal_side)

    def addbase(self, action, amp = None, side = None):
        if self.clean:
            if amp is not None:
                filename   = op.join ( self.origloc,        self.actionbase[amp] + self.basename + '_' + self.ifuslot + amp + '_' + self.type + '.fits')
                filename_e = op.join ( self.origloc, 'e.' + self.actionbase[amp] + self.basename + '_' + self.ifuslot + amp + '_' + self.type + '.fits')            
                self.actionbase[amp] = action + self.actionbase[amp]
                if op.exists  ( filename ):
                    os.remove ( filename )
                if op.exists  ( filename_e ):
                    os.remove ( filename_e )
                    
            if side is not None:
                filename   = op.join ( self.origloc,        self.actionbase[side] + self.basename + '_' + self.ifuslot + '_' + self.type + '_' + side + '.fits')
                filename_e = op.join ( self.origloc, 'e.' + self.actionbase[side] + self.basename + '_' + self.ifuslot + '_' + self.type + '_' + side + '.fits')            
                self.actionbase[side] = action + self.actionbase[side]
                if op.exists  ( filename ):
                    os.remove ( filename )
                if op.exists  ( filename_e ):
                    os.remove ( filename_e )                           
        else:
            if amp is not None:
                self.actionbase[amp] = action + self.actionbase[amp]
                
            if side is not None:
                self.actionbase[side] = action + self.actionbase[side]


#####################################################
# Define For Running Data Reduction Steps Functions #
#####################################################

def run_cure_command(command, suppress_output=0, debug=1):
    '''
       Run a cure command
       
       Parameter suppress_output is used to quiet down. (default)
       Return the value of os.system(command)
    '''
    # store info in log?

    if debug:
        print('Running: \'%s\'' % command)
    if not suppress_output:
        return os.system(op.join(CUREBIN, command) +' 1>>output.log  2>> error.log')
    else:
        return os.system(op.join(CUREBIN, command))
        
        
def mkerrorframe(frames, amp):
    
    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]

    command = 'mkerrorframe' 
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
    
    run_cure_command( command, 0 )
        
    return command
    
        
def subtractoverscan(biassec, frames, amp):
    
    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]

    command = 'subtractfits -s -a -k 2.8 -t -o %s -z' % (biassec)

    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('s',amp=amp) for f in frames] 

    return command
    
    
def subtractbias(frames, masterbiasname, amp):
        
    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]

    command = 'subtractfits -f %s' % ( masterbiasname ) 
        
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('s',amp=amp) for f in frames] 

    return command
    

def meanbiasfits(amp, specid, dest_dir, basename , opts, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]

    mastername = op.join ( dest_dir , basename + '_' + specid + '_' + amp + '.fits')
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('master',amp=amp) for f in frames]

    return command

def extractfits(trimsec, frames, amp):

    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]
    
    command = 'extractfits -r %s' % (trimsec)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command ( command, 0 )
    
    [f.addbase ('e',amp=amp) for f in frames] 
    
    return command


def ccdcombine(frames):
    
    for i in xrange(len(frames)):
        command = 'ccdcombine ' + op.join ( frames[i].origloc, frames[i].actionbase["LL"] + frames[i].basename + '_' + frames[i].ifuslot ) + ' ' + frames[i].type   
        run_cure_command( command, 0 )
    
    return command
    

def addphotonnoise(frames, side):
    
    filenames = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in frames]

    gain = filename[0].gain

    command = 'addphotonnoise -g %s' % (gain)

    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('p', side=side) for f in frames] 
    
    return command
    
    
def addfits(filename, opts):

    command = 'addfits %s %s' % (opts,filename)
    
    run_cure_command( command, 0 )
        
    return command 
    
    
def meanlampfits(side, specid, lamp, dest_dir, basename , opts, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in frames]

    mastername = op.join ( dest_dir , basename + '_' + lamp + '_' + specid + '_' + side + '.fits')
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    return command
    
    
def meantracefits(side, specid, dest_dir, basename , opts, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in frames]

    mastername = op.join ( dest_dir , basename + '_' + specid + '_' + side + '.fits')
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    return command
    

def deformer(mastertrace,masterarc,linesfile,wave_range,ref_line,opts):
    
    command = 'deformer %s -s %s -W %s -o \"%s\" -l %s -a %s %s' % (opts,ref_line,wave_range,redux_dir,linesfile,masterarc,mastertrace)  

    run_cure_command( command, 0 )

    return command

def fibextract(frames,base,side,distmodel,fibermodel,opts):

    filenames = [(redux_dir + '/' + sci_dir + '/' + base + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in frames]

    for f in filenames:
    
        command = 'fiberextract %s -x -d %s -f %s %s' %(opts,distmodel,fibermodel,f)
        
        run_cure_command( command, 0 )

    return command
        

def initial_setup ( file_loc_dir = None, redux_dir = None, DIR_DICT = None):
    '''
    Running the initial setup which includes:
    1) Building a standard reduction folder structure
    2) Copying files from specified location into the structure
    3) Creating class variable, VirusFrame, for each file which records
    keywords for other functions as well as time of the observation and basename
    '''
    vframes = [] # will fill this list with VirusFrame class objects for each image
    if redux_dir is None:
        print ( "Please provide a reduction directory \"redux_dir\"." )
        return None
    else:
        if not op.exists ( redux_dir ):
            os.mkdir ( redux_dir )
        else:
            if RESTART_FROM_SCRATCH:
                shutil.rmtree ( redux_dir )
                os.mkdir ( redux_dir )

    if file_loc_dir is None:        
        print ( "Please provide a file location directory \"file_loc_dir\"." )
        print ( "This is a list of the location directories for each file type." )
        return None
        
    if DIR_DICT is None:        
        print ( "Please provide a directory dictionary \"DIR_DICT\"." )
        print ( "The directory dictionary order should match the file location directory." )
        return None
    
    # Loop through the file location directories     
    for i in xrange ( len ( file_loc_dir ) ):
        # If the file_loc_dir[i] is None, then nothing is done
        if file_loc_dir[i] is not None:
            # If the reduction location exists, don't re-make the directory (also, any files in that directory remain)
            if not op.exists ( op.join ( redux_dir, DIR_DICT[i] ) ):
                os.mkdir ( op.join ( redux_dir, DIR_DICT[i] ) )
            # Run through the list of ObsID's, exposures, or other structures to gather all desired images for a particular type (e.g., zro)
            for file_loc in file_loc_dir[i]:
                # Trying to figure out what part of the directory structure was given to properly copy it
                dir1 = op.basename ( file_loc )
                if len ( dir1 ) > 5:
                    fit_path  = "/exp*/lrs2/*.fits"
                    filenames = glob.glob ( file_loc + fit_path )
                elif len ( dir1 ) == 5:
                    if dir1[0:3] == "exp": 
                        fit_path  = "/lrs2/*.fits"                   
                        filenames = glob.glob ( file_loc + fit_path )                    
                    else:      
                        fit_path  = "/*.fits"             
                        filenames = glob.glob ( file_loc + fit_path )                   
                else:               
                    print ( "Did not recognize the " + DIR_DICT[i] + " basename as \"lrs2XXXXXXX\", \"expXX\", or \"lrs2\"." )
                    return None

                # Loop through the retrieved files names to copy to new structure
                # Create a VirusFrame class for each frame that can maintain info for each original frame
                # The VirusFrame is critical for the rest of the reduction pipeline
                # Only gather VirusFrame objects for LL frames as a single VirusFrame serves for all four amps
                if all_copy:
                    #must copy all file to directories first
                    for f in filenames:  
                        shutil.copy ( f, op.join ( redux_dir, DIR_DICT[i] ) )
                    #once files copied it builds the virus frames for each
                    for f in filenames:            
                        temp, temp1, temp2 = op.basename ( f ).split('_')
                        amp                = temp1[3:5]
                        if amp == "LL":
                            a = VirusFrame( op.join( redux_dir, DIR_DICT[i], op.basename ( f ) ) ) 
                            vframes.append(copy.deepcopy(a))
                else:
                    #if not all_copy does not copy files to directories but builds virus frames for each
                    for f in filenames:            
                        temp, temp1, temp2 = op.basename ( f ).split('_')
                        amp                = temp1[3:5]
                        if amp == "LL":
                            a = VirusFrame(  f  ) 
                            vframes.append(copy.deepcopy(a))
                        
    return vframes, first_run

def basicred( file_loc_dir, redux_dir, DIR_DICT, basic = False, dividepf = False,
              normalize = False, masterdark = False, masterarc = False, mastertrace = False):
    '''
    Running the basic reduction which includes:
    1) Overscan subtract and trim all frames
    2) Create blank error frame with readnoise in it
    3) Remove cosmic rays from sci frames (user option)
    4) Create master bias frame for each amp 
    5) Subtract master bias from cmp/flt/sci frames
    6) Create master dark for each amp make scaled frame for each sci exposure time (user option)
    7) Subtract master dark from each science frame (user option)
    8) Ccdcombine all frames which puts units in e-
    9) Add photon noise to combined frames
    10) Divide pixelflat from cmp/flt/sci frames (user option)
    11) Normalize cmps and flts 
    12) Combine cmps and flts into masterarc and mastertrace
    '''

    print ('*************************')
    print ('* BUILDING IMAGE FRAMES *')
    print ('*************************')

    #holds the VIRUS frames for all of the data 
    vframes, first_run = initial_setup ( file_loc_dir, redux_dir, DIR_DICT )

    #first frames to pull basic header information 
    f1 = vframes[0]

    #from the vframes makes lists of files vframes according to type and specid (ucam)
    tframes  = [v for v in vframes if (v.specid == ucam) ] # gives all frames 

    oframes  = [t for t in tframes if t.type != "zro" ] # gives "flt" and "cmp" frames (basically just not "zro")
    lframes  = [t for t in tframes if t.type == "cmp" ] # gives just "cmp" frames
    fframes  = [t for t in tframes if t.type == "flt" ] # gives just "flt" frames
    zframes  = [t for t in tframes if t.type == "zro" ] # gives just "zro" frames

    #Check that data is correct
    if len(lframes) == 0:
        print ("You did not provide the correct comp frames for this data set")
        sys.exit("Either the flt_folder you provided are not flats or these are not for LRS2-"+LRS2_spec)
    if len(fframes) == 0:
        print ("You did not provide the correct flat frames for this data set")
        sys.exit("Either the Hg/Cd/FeAr_folder you provided are not comps or these are not for LRS2-"+LRS2_spec)
    if len(zframes) == 0:
        print ("You did not provide zero frames for this data set")
        sys.exit("Check the data type of the zro_folder you provided to make sure they are zro images")

    # Run basic reduction
    if basic:
        for sp in SPEC:
            trimsec = f1.trimsec["LL"] # Trimsec assumed to be the same for all frames of a given amp
            biassec = f1.biassec["LL"] # Biassec assumed to be the same for all frames of a given amp
            print ('**************************')
            print ('* MAKE ERROR FRAME FOR '+sp+' *')
            print ('**************************')
            mkerrorframe ( tframes, sp )               # for all frames

            print ('***************************')
            print ('* SUBTRACT OVERSCAN FOR '+sp+' *')
            print ('***************************')
            subtractoverscan( biassec, tframes, sp )   # for all frames

            print ('*****************************')
            print ('* EXTRACT DATA REGION FOR '+sp+' *')
            print ('*****************************')
            extractfits ( trimsec, tframes, sp )       # for all frames                   
            
            print ('**************************')
            print ('* BUILD MASTERBIAS FOR '+sp+' *')
            print ('**************************')
            vframesselect  = [v for v in tframes if v.type == "zro" and v.specid == ucam] 
            meanbiasfits   ( sp, ucam, redux_dir, 'masterbias', meanfitsopts, vframesselect ) # meanfits for masterbias for unique specid 

            print ('*****************************')
            print ('* SUBTRACT MASTERBIAS FOR '+sp+' *')
            print ('*****************************')
            masterbiasname = op.join ( redux_dir, 'masterbias' + '_' + ucam + '_' + sp + '.fits' ) 
            oframesselect  = [o for o in oframes if o.specid == ucam]
            subtractbias   ( oframesselect, masterbiasname, sp) # for all frames

                
        print ('*******************************************')
        print ('* COMBINING CCD AMPS - MULTIPYING BY GAIN *')
        print ('*******************************************')

        ccdcombine ( oframes ) # Combine amps and multiply by gain for all frames other than zro's
        for o in oframes:
            o.actionbase["L"] = o.actionbase["LL"]
            o.actionbase["R"] = o.actionbase["RL"]
            if o.clean:
                for sp in SPEC:
                    filename   = op.join ( o.origloc,        o.actionbase[sp] + o.basename + '_' + o.ifuslot + sp + '_' + o.type + '.fits')
                    filename_e = op.join ( o.origloc, 'e.' + o.actionbase[sp] + o.basename + '_' + o.ifuslot + sp + '_' + o.type + '.fits')
                    os.remove ( filename )
                    os.remove ( filename_e )
        for side in SPECBIG:    
            addphotonnoise ( oframes , side )
    
    # Make Master Arc Frames           
    if masterarc:        
        print ('**********************')
        print ('* BUILDING MASTERARC *')
        print ('**********************')
        for side in SPECBIG:
            for lamp in LAMP_DICT.values():
                #Creates a masterarc frame for each arc lamp in LAMP_DICT
                lframesselect = [l for l in lframes if l.specid == ucam]
                #If more than one image for that lamp take median image 
                if len(lframesselect)>1:
                    meanlampfits(side, ucam, lamp, redux_dir, 'masterarc' , arcopts, lframesselect) 
                #If only one frame for the lamp make that frame the master frame for that lamp 
                elif len(lframesselect) == 1:
                    filename = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in lframesselect]
                    mastername = op.join ( redux_dir , 'masterarc' + '_' + lamp + '_' + ucam + '_' + side + '.fits')
                    shutil.copy ( filename[0], mastername )
                    efilename = [op.join ( f.origloc, 'e.' + f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in lframesselect]
                    emastername = op.join ( redux_dir , 'e.masterarc' + '_' + lamp + '_' + ucam + '_' + side + '.fits')
                    shutil.copy ( efilename[0], emastername )
                #If there are zero frames found then the arc lamp information provided is wrong
                else:
                    sys.exit("You did not provide the correct arc lamp data")

            #Combine each arc master frame for each lamp in LAMP_DICT into one masterarc 
            opt = "--file {:s}".format( op.join ( redux_dir, 'masterarc' + '_' + LAMP_DICT[0] + '_' + ucam + '_' + side + '.fits' ) )
            filename = op.join ( redux_dir, 'masterarc' + '_' + LAMP_DICT[1] + '_' + ucam + '_' + side + '.fits' )
            addfits ( filename, opt)
            shutil.copy( op.join ( redux_dir, 'smasterarc' + '_' + LAMP_DICT[1] + '_' + ucam + '_' + side + '.fits' ), 
                         op.join ( redux_dir, 'masterarc' + '_' + ucam + '_' + side + '.fits' ) )
            shutil.copy( op.join ( redux_dir, 'e.smasterarc' + '_' + LAMP_DICT[1] + '_' + ucam + '_' + side + '.fits' ), 
                         op.join ( redux_dir, 'e.masterarc' + '_' + ucam + '_' + side + '.fits' ) )
            os.remove  ( op.join ( redux_dir, 'smasterarc' + '_' + LAMP_DICT[1] + '_' + ucam + '_' + side + '.fits' ) )

        #clean intermediate frames in cmp if clean
        for l in lframes:
            if l.clean:
                for side in SPECBIG:
                    filename   = op.join ( l.origloc,        l.actionbase[side] + l.basename + '_' + l.ifuslot + '_' + l.type + '_' + side + '.fits')
                    filename_e = op.join ( l.origloc, 'e.' + l.actionbase[side] + l.basename + '_' + l.ifuslot + '_' + l.type + '_' + side + '.fits')
                    os.remove ( filename )
                    os.remove ( filename_e )
         
    # Make Master Trace Frames
    if mastertrace:
        print ('************************')
        print ('* BUILDING MASTERTRACE *')
        print ('************************')
        for side in SPECBIG:
            #Build a list of flat frames matching the correct lamp
            fframesselect = [f for f in fframes if f.specid == ucam] 
            #If there is more than one frame take a median from of the ones provided for the mastertrace
            if len ( fframesselect ) > 1:
                meantracefits(side, ucam, redux_dir, 'mastertrace', traceopts, fframesselect)
            #If there was only one frame provided that frame becomes the mastertrace
            else:
                filename = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in fframesselect]
                mastername = op.join ( redux_dir , 'mastertrace' + '_' + ucam + '_' + side + '.fits')
                shutil.copy ( filename[0], mastername )
                efilename = [op.join ( f.origloc, 'e.' + f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in fframesselect]
                emastername = op.join ( redux_dir , 'e.mastertrace' + '_' + ucam + '_' + side + '.fits')
                shutil.copy ( efilename[0], emastername )       
          
    # Run Deformer
    if run_deformer:
        print ('*************************************************************************')
        print ('* RUNNING DEFORMER TO BUILD DISTORTION SOLUTION AND WAVELENGTH SOLUTION *')
        print ('*************************************************************************')
        #check that basic has been run 
        trace_files = glob.glob(op.join(redux_dir,'mastertrace*'))
        arc_files   = glob.glob(op.join(redux_dir,'masterarc*'))
        if len(trace_files) == 0 or len(arc_files) == 0:
            sys.exit("You must run basic reduction before you can run deformer")

        for side in SPECBIG:  
            #copy the lines file used to this directory 
            shutil.copy ( op.join(linesdir,'lines' + '_' + side + '_' + ucam +'.par'), op.join(redux_dir,'lines' + '_' + side + '_' + ucam +'.par' ) )
            #build the names of the files given to deformer 
            mastertrace = op.join ( redux_dir, 'mastertrace' + '_' + ucam + '_' + side + '.fits' )
            masterarc   = op.join ( redux_dir, 'masterarc' + '_' + ucam + '_' + side + '.fits' )
            linefile    = op.join ( redux_dir, 'lines' + '_' + side + '_' + ucam +'.par' )
            deformer ( mastertrace, masterarc, linefile, wave_range, ref_line, deformeropts)

    # Run fiberextract
    if fiberextract:  
        print ('**********************************************************')
        print ('* EXTRACTING SPECTRA IN FLAT FRAMES - Without Resampling *')
        print ('**********************************************************')
        #check that deformer has been run 
        dist_files = glob.glob(op.join(redux_dir,'*.dist'))
        if len(dist_files) == 0:
            sys.exit("You must run deformer before you can run fiber extract")

        base = 'pses'
        for side in SPECBIG:   
            fframesselect = [f for f in fframes if s.specid == ucam]

            distmodel = redux_dir + "/mastertrace_" + ucam + "_" + side + ".dist"
            fibermodel = redux_dir + "/mastertrace_" + ucam + "_" + side + ".fmod"
            fibextract(fframesselect,base,side,distmodel,fibermodel,fibextractopts)

    #CURE saves these files from deformer outside of the redux directory for some reason.
    #This moves them inside of the redux directory.
    left_files = glob.glob('*.log') + glob.glob('*.residuals')
    if len(left_files) > 0:
        for l in left_files:
            os.rename(l, op.join(redux_dir,l))

    return vframes
    
def main():
    frames = basicred( file_loc_dir, redux_dir, DIR_DICT, basic = basic, dividepf = dividepf,
                      normalize = normalize, masterdark = masterdark, masterarc = masterarc, mastertrace = mastertrace )                 
    
if __name__ == '__main__':
    main()  
