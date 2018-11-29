#!/bin/bash

sublist=$(<../subjRes4dPath.txt)
for i in $sublist;
  do
  j=${i:20:3}
  k=${i:29:1}

  #extract first volume 
  #fslsplit $i /danl/Harmon_dynCon/scripts/QA/$j #This also works but seperates into all volumes 
  fslroi $i /danl/Harmon_dynCon/scripts/QA/$j$k 0 1  #this just gives us the first volume 
done

fslmerge -t res4d_std_vol1 7*
