#!/bin/bash

flirt -ref $1/../structural/mprage.anat/T1_biascorr_brain.nii.gz\
 -in $1/reg/example_func.nii.gz\
 -omat $1/reg/example_func2highres.mat;

echo warping $1;

#apply warp from FNIRT to preprocessed 4D data
applywarp --ref=$2/standard/MNI152_T1_2mm.nii.gz\
 --in=$1/36par+spikes.feat/stats/res4d.nii.gz\
 --out=$1/36par+spikes.feat/stats/res4d_std.nii.gz\
 --warp=$1/../structural/mprage.anat/T1_to_MNI_nonlin_field.nii.gz\
 --premat=$1/reg/example_func2highres.mat;




