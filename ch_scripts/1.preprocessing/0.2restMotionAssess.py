#!/usr/bin/python
# script from Mumford Brain stats to assess Motion in BOLD data
# updated by MVT in Jan 2018 and PAB Feb 13, 2018 and CH May 10, 2018

# DO THIS BEFORE YOU PREPROCESS YOUR DATA!
# this gets Motion information on each subject for each run. 
# makes an output in the SB/PaulCompileTGNG/data/SBxxx/BOLD/motion_assess directories where you can find very useful information! 

import glob
import os
import sys
import subprocess

path = '/danl/Harmon_dynCon/'
 
# get the paths for each BOLD.nii file 
rest_files = glob.glob('%s/7*/Rest/Rest?.nii'%(path))
print (rest_files)


#
# Create a big html file to put all QA info together. 
outhtml = '/danl/Harmon_dynCon/scripts/1.preprocessing/motion.qa.html'

# save lists of subject runs with good motion to keep vs. bad motion to exclude
out_bad_bold_list =  '/danl/Harmon_dynCon/scripts/1.preprocessing/Subject_high_motion_exclude.txt'
out_good_bold_list =  '/danl/Harmon_dynCon/scripts/1.preprocessing/Subject_runs_good_motion_include.txt'

# delete these files so that if you run this script twice, it doesn't keep apending to the old stuff. 
#os.system("rm %s"%(out_bad_bold_list))
#os.system("rm %s"%(outhtml))



# loop through each directory for each bold.nii file - this is just going to show me my BOLD files
for cur_bold in list(rest_files): # Subset rest_files with list indexing to do smaller batches
    print(cur_bold)
    # Store directory name
    cur_dir = os.path.dirname(cur_bold)
     # get current bold file
    cur_bold_no_nii = cur_bold
    
    # Assessing motion using FSL_motion_outliers function
    if os.path.isdir("%s/motion_assess/"%(cur_dir))==False:
      os.system("mkdir %s/motion_assess"%(cur_dir))
      # lots of options: fd = framewise displacement 
      # threshold = setting to a value of 0.9mm which is from the Siegel 2014 paper! 
      # confound.txt file = what you can put into the design matrix 
      # this will generate the column of 1 and 0 for the frames to remove = 1 
      # -p = plot 
      # -v = verbose mode which spits out a lot of information 
      # this will generate a motion_assess folder within the BOLD directory
    os.system("fsl_motion_outliers -i %s -o %s/motion_assess/confound.txt --fd --thresh=0.9 -s -p %s/motion_assess/fd_plot -v > %s/motion_assess/outlier_output.txt"%(cur_bold_no_nii, cur_dir, cur_dir, cur_dir))

    # 4. Put confound info into html file for review later on
    os.system("cat %s/motion_assess/outlier_output.txt >> %s"%(cur_dir, outhtml))
    os.system("echo '<p>=============<p>FD plot %s <br><IMG BORDER=0 SRC=%s/motion_assess/fd_plot.png WIDTH=100%s></BODY></HTML>' >> %s"%(cur_dir, cur_dir,'%', outhtml))

    # 5. Last, if we're planning on modeling out scrubbed volumes later
    #   it is helpful to create an empty file if confound.txt isn't
    #   generated (i.e. for subjects with no scrubbing needed).  It is basically a
    #   place holder to make future scripting easier
    if os.path.isfile("%s/motion_assess/confound.txt"%(cur_dir))==False:
      os.system("touch %s/motion_assess/confound.txt"%(cur_dir))

    # 6. Very last, create a list of subjects who exceed a threshold for
    #  number of scrubbed volumes.  This should be taken seriously.  If
    #  most of your scrubbed data are occurring during task, that's
    #  important to consider (e.g. subject with 20 volumes scrubbed
  

# PAB notes 
# we have 130 TRs per run so i'll keep a run unless more than 40 TRs have > 0.9 mm  FD

    output = subprocess.check_output("grep -o 1 %s/motion_assess/confound.txt | wc -l"%(cur_dir), shell=True)
    num_scrub = [int(s) for s in output.split() if s.isdigit()]
    if num_scrub[0]> 40: 
        with open(out_bad_bold_list, "a") as myfile:
        	myfile.write("%s\n"%(cur_bold))
   
    if num_scrub[0]<= 40: 
        with open(out_good_bold_list, "a") as myfile:
        	myfile.write("%s\n"%(cur_bold))
