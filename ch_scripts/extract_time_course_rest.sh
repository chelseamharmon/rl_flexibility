#!/bin/bash

for i in /danl/Harmon_dynCon/7*/Learn?_PEpriorD.feat/36par+spikes.feat;
  do
  #extract timeseries (mean or 1st eigenvector, see function) data from each ROI in ~/Harvard-Oxford_ROIs/
  /danl/Harmon_dynCon/rl_flexibility/extract_ROIs.sh $i/stats/res4d_std.nii.gz /danl/Harmon_dynCon/Harvard-Oxford_ROIs/ $i/H-O_rois/;
done



