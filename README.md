# Network Flexibility and Reinforcement Learning
## Chelsea Harmon, 2018 
### adapted from : Raphael Gerraty, 2015-2016

Descriptions and example scripts for running network preprocessing and analysis functions contained in this repository. See paper for details when it comes out.

### Syncing data 


```.bash
#From original location to server 
rsync -avz -O --omit-dir-times --no-perms  --include="7*" --include="7*/Learn?_PEpriorD.feat" --include="7*/mprage.nii" --include="7*/Rest?.nii" --include="7*/Learn?_PEpriorD.feat/filtered_func_data.nii.gz" --include="7*/Learn?_PEpriorD.feat/reg" --include="7*/Learn?_PEpriorD.feat/mc" --include="7*/Learn?_PEpriorD.feat/mc/*" --include="7*/Learn?_PEpriorD.feat/reg/*" --exclude="*" --exclude="*/*" --exclude="*/*/*" rgerraty@lovelace.psych.columbia.edu:/data/engine/juliet/adoles/ /danl/Harmon_dynCon/

#From server to Habanero 
rsync -avz -O --omit-dir-times --no-perms  --include="7*" --include="7*/Learn?_PEpriorD.feat" --include="7*/Learn?_PEpriorD.feat/*" --include="7*/Learn?_PEpriorD.feat/reg" --include="7*/Learn?_PEpriorD.feat/mc" --include="7*/Learn?_PEpriorD.feat/mc/*" --include="7*/Learn?_PEpriorD.feat/reg/*" --exclude="*" --exclude="*/*" --exclude="*/*/*" charmon@lux.psych.columbia.edu:/danl/Harmon_dynCon/ /rigel/psych/users/cmh2228/dynCon/
```



### Preprocessing 
Prepare folder structure. Run fsl_anat to prepare anatomical scans. Run initial preprocessing of resting state data. 

```.bash
#lux
/danl/Harmon_dynCon/scripts/0.set_up_scripts
./0.1preprepFolderStructure.sh

```
```.bash
#habanero
/rigel/psych/users/cmh2228/dynCon/scripts/
for n in 7* ; do sbatch fslAnat.sh $n ; done 

for n in 7* ; do sbatch restFeat.sh $n ; done

```

```.bash
#lux
/danl/Harmon_dynCon/scripts/1.preprocessing
./restReg.sh

```


### Extended Preprocessing
Because of the known effect of motion on measures of connectivity, we followed up standard preprocessing in FSL with an extended nuisance regression. Affine transformation parameters from motion correction, CSF, white matter, and whole-brain signals are regressed against preprocessed 4D data, along with the squares, derivatives, and squared derivatives of these confounds. See Satterthwaite et al 2013 for details. Bash code:

For Task 
```.bash
#do this on lux 
/danl/Harmon_dynCon/scripts/1.preprocessing
./1.0extract_confts

#do this on habanero
for i in /rigel/psych/users/cmh2228/dynCon/7*/Learn*/filtered_func_data.nii.gz; 
do 
subdir=$(dirname $i); 
sbatch --export=arg1=$i,arg2=$subdir/36par+spikes.txt,arg3=/rigel/psych/users/cmh2228/dynCon/rl_flexibility/conf_reg_design.fsf,arg4=/rigel/psych/app/fsl/data,arg5=$subdir/conf_reg_design.fsf /rigel/psych/users/cmh2228/dynCon/scripts/1.preprocessing/1st_level_conf.sub.sh; 
done
```

For Rest 
```.bash

rsync -avz -O --omit-dir-times --no-perms --include="7*" --include="7*/Rest" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/*" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/logs" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/logs/*" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/mc" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/mc/*" --exclude="*" --exclude="*/*" --exclude="*/*/*" --exclude="*/*/*/*" cmh2228@habanero.rcs.columbia.edu:/rigel/psych/users/cmh2228/dynCon/ /danl/Harmon_dynCon/ 

#must register data before doing: 

#do this on lux 
/danl/Harmon_dynCon/scripts/1.preprocessing
./1.0extract_confts_rest

#do this on habanero
rsync -avz -O --omit-dir-times --no-perms --include="7*" --include="7*/Rest" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/*" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/logs" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/logs/*" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/mc" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/mc/*" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/reg" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/reg/*" --exclude="*" --exclude="*/*" --exclude="*/*/*" --exclude="*/*/*/*" charmon@lux.psych.columbia.edu:/danl/Harmon_dynCon/ /rigel/psych/users/cmh2228/dynCon/

rsync -avz -O --omit-dir-times --no-perms --include="7*" --include="7*/Rest" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/*" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/reg" --include="7*/Rest/Rest?.feat/reg/*" --exclude="*" --exclude="*/*" --exclude="*/*/*" --exclude="*/*/*/*" charmon@lux.psych.columbia.edu:/danl/Harmon_dynCon/ /rigel/psych/users/cmh2228/dynCon/


for i in /rigel/psych/users/cmh2228/dynCon/7*/Rest/Rest?.feat/filtered_func_data.nii.gz; 
do 
subdir=$(dirname $i); 
sbatch --export=arg1=$i,arg2=$subdir/36par+spikes.txt,arg3=/rigel/psych/users/cmh2228/dynCon/rl_flexibility/conf_reg_design.fsf,arg4=/rigel/psych/app/fsl/data,arg5=$subdir/conf_reg_design.fsf /rigel/psych/users/cmh2228/dynCon/scripts/1.preprocessing/1st_level_conf.sub.sh; 
done

```


### Running nonlinear registration with FNIRT
After nuisance regression has been run, the residual timeseries needs to be transformed into standard space (in this case, MNI). Make sure fsl\_anat has been run on each structural image first. The following bash code was used to perform these transformations:

for Task
```.bash
#fnirt has already been run, just applying transformation

for i in /rigel/psych/users/cmh2228/dynCon/7*/Learn?_PEpriorD.feat;
do sbatch --export=arg1=$i,arg2=/rigel/psych/app/fsl/data /rigel/psych/users/cmh2228/dynCon/scripts/non_lin_reg_sub.sh; 
done

```

for Rest
```.bash
#fnirt has already been run, just applying transformation

for i in /rigel/psych/users/cmh2228/dynCon/7*/Rest/Rest?.feat;
do sbatch --export=arg1=$i,arg2=/rigel/psych/app/fsl/data /rigel/psych/users/cmh2228/dynCon/scripts/non_lin_reg_sub.sh; 
done
```

### Syncing data 


```.bash
#Moving newly created extended preprocesses folders (36par+spikes.feat) to lux 

#Practice - dry run displays the files that will be copied 
rsync --dry-run -avz --include="7*" --include="7*/Learn?_PEpriorD.feat" --include="7*/Learn?_PEpriorD.feat/36*.feat" --include="7*/Learn?_PEpriorD.feat/36*.feat/stats" --include="7*/Learn?_PEpriorD.feat/36*.feat/logs/" --include="7*/Learn?_PEpriorD.feat/36*.feat/stats/*" --include="7*/Learn?_PEpriorD.feat/36*.feat/logs/*" --exclude="*" --exclude="*/*" --exclude="*/*/*" --exclude="*/*/*/*" --exclude="*/*/*/*/*" cmh2228@habanero.rcs.columbia.edu:/rigel/psych/users/cmh2228/dynCon/ /danl/Harmon_dynCon/ 

#Task
rsync -avz -O --omit-dir-times --no-perms --include="7*" --include="7*/Learn?_PEpriorD.feat" --include="7*/Learn?_PEpriorD.feat/36*.feat" --include="7*/Learn?_PEpriorD.feat/36*.feat/stats" --include="7*/Learn?_PEpriorD.feat/36*.feat/logs/" --include="7*/Learn?_PEpriorD.feat/36*.feat/stats/*" --include="7*/Learn?_PEpriorD.feat/36*.feat/logs/*" --exclude="*" --exclude="*/*" --exclude="*/*/*" --exclude="*/*/*/*" --exclude="*/*/*/*/*" cmh2228@habanero.rcs.columbia.edu:/rigel/psych/users/cmh2228/dynCon/ /danl/Harmon_dynCon/ 

#Rest 
rsync -avz -O --omit-dir-times --no-perms --include="7*" --include="7*/Rest" --include="7*/Rest/Rest?.feat" --include="7*/Rest/Rest?.feat/36*.feat" --include="7*/Rest/Rest?.feat/36*.feat/*" --include="7*/Rest/Rest?.feat/36*.feat/stats" --include="7*/Rest/Rest?.feat/36*.feat/logs/" --include="7*/Rest/Rest?.feat/36*.feat/stats/*" --include="7*/Rest/Rest?.feat/36*.feat/logs/*" --exclude="*" --exclude="*/*" --exclude="*/*/*" --exclude="*/*/*/*" --exclude="*/*/*/*/*" --exclude="*/*/*/*/*/*" cmh2228@habanero.rcs.columbia.edu:/rigel/psych/users/cmh2228/dynCon/ /danl/Harmon_dynCon/ 

```


### Extracting time courses
Once the preprocessed images have been registered, we extract mean timecourses for each Harvard-Oxford ROI, using the function extract_ROIs.sh. The output of this function is a timecourse for each ROI in the specified input folder, as well as a .txt file containing all of the ROIs. The bash code used to run this function on each learning block for each subject is below:

```.bash
#scp Harvard-Oxford_ROIs folder from rgerraty

for i in /danl/Harmon_dynCon/7*/Learn?_PEprior.feat/36par+spikes.feat/; 
    do 
    #extract timeseries (mean or 1st eigenvector, see function) data from each ROI in ~/Harvard-Oxford_ROIs/ 
    ~/GitHub/rl_flexibility/extract_ROIs.sh $i/stats/res4d_std.nii.gz /danl/Harmon_dynCon/Harvard-Oxford_ROIs/ $i/H-O_rois/;
done

#Also can be run by running script ./extract_time_course.sh
```

for Rest 

```.bash

nohup ./extract_time_course_rest.sh &

```




### Calculate coherence matrices for each time window
Connectivity between pairs of ROIs was measured by average magnitude squared coherence in the .06-.12 Hz band, computed in MATLAB. The code below calls a function for creating a coherence matrices in a specified frequency range for specified time windows (in this case 25 TRs, or 50 s). These are saved as a .mat file for multi-slice community detection. 

for Task 
```.matlab
screen -r 
matlab -nosplash -nojvm 
%Alternatively(and slower) matlab -nosplash -nodesktop

addpath ~/GitHub/rl_flexibility
%read in all subject/run ROI timeseries directories 
[a,b]=system('ls -d /danl/Harmon_dynCon/7*/Learn?_PEpriorD.feat/36par+spikes.feat/H-O_rois');
%do separately above section followed by next section 

c=strread(b,'%s');

for i=1:size(c,1)

    %calculate coherence per time window from concatenated ROI file
    filename=char(strcat(c(i),'/all_rois.txt'))
    %need to specify filename, window length in TR, sampling rate, bandpass 
    conn_cell=coherence_by_block(filename,25,.5,.06,.12);
    save(char(strcat(c(i),'/conn_cells')),'conn_cell')

end

```

Subject 713 has 2 missing ROIs (column 13 & 61) - this code removes those columns in the all_rois.txt file before running community detection
Rerun calculating the coherence matrices having removed these two columns from all 4 Learning blocks for this subject
After rerunning add columns and rows back in 

```.matlab

addpath ~/GitHub/rl_flexibility
%read in all subject/run ROI timeseries directories 
[a,b]=system('ls -d /danl/Harmon_dynCon/713/Learn4_PEpriorD.feat/36par+spikes.feat/H-O_rois');
%do separately above section followed by next section 

c=strread(b,'%s');

incomp_rois=load(char(strcat(c(4),'/all_rois.txt')))
[row, col]=find(~incomp_rois) 
column_missing=unique(col) 

%This is one way to remove missing ROIs for the file that has 0's incomp_rois( :, ~any(incomp_rois,1) ) = []
%Below is to get rid of excluded coher matrices for all runs for this subject 

%Creating new ROI txt files excluding NA columns 

[a,b]=system('ls -d /danl/Harmon_dynCon/713/Learn?_PEpriorD.feat/36par+spikes.feat/H-O_rois');

c=strread(b,'%s');
for i=1:size(c,1)
    incomp_rois=load(char(strcat(c(i),'/all_rois.txt')))
    incomp_rois(:,[13 61]) = []
    save(char(strcat(c(1),'/all_rois_incomp.txt')),'incomp_rois', '-ascii')
end 

%calculating coherence matrices excluding NA rois 

[a,b]=system('ls -d /danl/Harmon_dynCon/713/Learn?_PEpriorD.feat/36par+spikes.feat/H-O_rois');
%do separately above section followed by next section 

c=strread(b,'%s');

for i=1:size(c,1)

    %calculate coherence per time window from concatenated ROI file
    
    filename=char(strcat(c(i),'/all_rois_incomp.txt'))
    %need to specify filename, window length in TR, sampling rate, bandpass 
    conn_cell=coherence_by_block(filename,25,.5,.06,.12);
    save(char(strcat(c(i),'/conn_cells')),'conn_cell')

end




```



for Rest 
```.matlab
screen -r 
matlab -nosplash -nojvm 
%Alternatively(and slower) matlab -nosplash -nodesktop

addpath ~/GitHub/rl_flexibility
%read in all subject/run ROI timeseries directories 
[a,b]=system('ls -d /danl/Harmon_dynCon/7*/Rest/Rest?.feat/36par+spikes.feat/H-O_rois');
%do separately above section followed by next section 

c=strread(b,'%s');

for i=1:size(c,1)

    %calculate coherence per time window from concatenated ROI file
    filename=char(strcat(c(i),'/all_rois.txt'))
    %need to specify filename, window length in TR, sampling rate, bandpass 
    conn_cell=coherence_by_block(filename,25,.5,.06,.12);
    save(char(strcat(c(i),'/conn_cells')),'conn_cell')

end

```



### Run multi-slice community detection and flexibility statistics

Input coherence matrix for each block. Also need number of blocks,
resolution and coupling parameters. In Matlab

Because Subject 713 has 2 missing ROIs - add them back to conn_cells here 

forTask
```.matlab
#Do separately for 70* 71* and 72* by using: 
screen 
control^ a d 
screen -r 
matlab -nosplash -nodisplay 


%need multi-slice, flexibility codes not yet on GitHub for network_diags to run 
addpath ~/Github/rl_flexibility
addpath ~/Github/rl_flexibility/GenLouvain/
addpath ~/Github/rl_flexibility/Bassett_Code/

%read in data
[a,b]=system('ls -d /danl/Harmon_dynCon/72*/Learn?_PEpriorD.feat/36par+spikes.feat/H-O_rois/');

[a,b]=system('ls -d /danl/Harmon_dynCon/713/Learn?_PEpriorD.feat/36par+spikes.feat/H-O_rois/');


%do the above first then do the below 

c=strread(b,'%s');

%concatenate runs for each subject
numruns=4
k=1;
for j=1:size(c,1)/numruns
    c(k)
    conn_cell_cat=[];
    for i=1:numruns 
        load(strcat(char(c(k-1+i)),'/conn_cells'))
        conn_cell_cat=cat(3,conn_cell_cat,conn_cell)
    end

    %network_diags code:
    %runs multi-slice community detection
    %gives flexibility for each run
    %also allegiance matrix (not using yet)
    %need to specify number of blocks, simulations, coupling, resolution
    [a_mat,flex]=network_diags(conn_cell_cat,4,500,1,1.1813)
    save(char(strcat(c(k),'/../../../a_mat')),'a_mat')
    save(char(strcat(c(k),'/../../../flex')),'flex')
    k=k+numruns;
end
```

forRest
```.matlab
screen 
control^ a d 
screen -r 
matlab -nosplash -nodisplay 


%need multi-slice, flexibility codes not yet on GitHub for network_diags to run 
addpath ~/Github/rl_flexibility
addpath ~/Github/rl_flexibility/GenLouvain/
addpath ~/Github/rl_flexibility/Bassett_Code/

%read in data
[a,b]=system('ls -d /danl/Harmon_dynCon/7*/Rest/Rest?.feat/36par+spikes.feat/H-O_rois/');

%do the above first then do the below 

c=strread(b,'%s');

%concatenate runs for each subject
numruns=2
k=1;
for j=1:size(c,1)/numruns
    c(k)
    conn_cell_cat=[];
    for i=1:numruns 
        load(strcat(char(c(k-1+i)),'/conn_cells'))
        conn_cell_cat=cat(3,conn_cell_cat,conn_cell)
    end

    %network_diags code:
    %runs multi-slice community detection
    %gives flexibility for each run
    %also allegiance matrix (not using yet)
    %need to specify number of blocks, simulations, coupling, resolution
    [a_mat,flex]=network_diags(conn_cell_cat,4,500,1,1.1813)
    save(char(strcat(c(k),'/../../../a_mat')),'a_mat')
    save(char(strcat(c(k),'/../../../flex')),'flex')
    k=k+numruns;
end

```

To add back in the missing columns/rows for subject 713 before pulling all subj flexibility stats 
#Should be updated for a_mat.mat & flex.mat 
```
%do this from a script in matlab in the Learn folder seperately to index properly - to add back NA columns before concatenating

load('flex.mat')
flexTEMP=flex
missing=[13 61; 61 13]

    for i=1:length(missing)
        %FOR inserting rows
        len=size(flexTEMP,i);
        A = [flexTEMP(1:missing(i)-1,:)];
        B = NaN(1,4);
        C = [flexTEMP(missing(i):end,:)];
        flexTEMP=[A; B; C];
    end
save('flexTEMP', 'flexTEMP')

load('a_mat.mat')
a_matTEMP=a_mat;
a_matTEST=zeros(110,110,32);
missing=[13 61; 61 13];
for k=1:size(a_mat,3)
    a_matTEMP=a_mat(:,:,k)
    missing=[13 61; 61 13]
    for i=1:length(missing)
            %For inserting rows
            len=size(a_matTEMP,i);
            A = [a_matTEMP(1:missing(i)-1,:)]';
            B = NaN(1,len)';
            C = [a_matTEMP(missing(i):end,:)]';
            a_matTEMP1=[A B C];
            %For inserting columns
            len=size(a_matTEMP1,2)
            A = [a_matTEMP1(1:missing(i)-1,:)];
            B = NaN(1,len);
            C = [a_matTEMP1(missing(i):end,:)];
            a_matTEMP=[A; B; C];
    end
    a_matTEST(:,:,k)=a_matTEMP;
end

```


### Pull flexibility statistics

For plotting and preparing for heirarchical models. Matlab.

forTask
```.matlab

%load data and concatenate flexibility statistics
[a,b]=system('ls -d /danl/Harmon_dynCon/7*/flex.mat');
%do separately above and below 


c=strread(b,'%s');
flex_cat=[];
for j=1:size(c,1)
    load(char(c(j)))
    flex_cat=cat(3,flex_cat,flex)
end
plot(squeeze(mean(flex_cat)))

block=repmat([1:4]',25,1);
sub=repmat([1:25]',1,4)'
sub=sub(:);

%reshape whole-brain average flexibility
meanflex=squeeze(mean(flex_cat));
meanflex=meanflex(:);

%get striatal average flexibility
%check to make sure indices are correct
[trash,roi_names]=system('ls  /danl/Harmon_dynCon/Harvard-Oxford_ROIs/*nii.gz | xargs -n1 basename');
roi_names=strread(roi_names,'%s');
str_ind=[49,51,54,104,106,109];
roi_names(str_ind)
roi_names(hipp_ind)

strflex=squeeze(mean(flex_cat(str_ind,:,:)));
strflex=strflex(:);

plot(squeeze(mean(flex_cat(str_ind,:,:))))

%Hippocampus 
hipp_ind=[9, 69];
left_hipp=[9];
right_hipp=[69];
hippflex=squeeze(mean(flex_cat(hipp_ind,:,:)));
hippflex=hippflex(:);
hippflexR=squeeze(mean(flex_cat(right_hipp,:,:)));
hippflexR=hippflexR(:);
hippflexL=squeeze(mean(flex_cat(right_hipp,:,:)));
hippflexL=hippflexL(:);
plot(squeeze((flex_cat(right_hipp,:,:))))
plot(squeeze((flex_cat(left_hipp,:,:))))

plot(squeeze(mean(flex_cat(hipp_ind,:,:))))

%write out csv for modeling in R
flexdata=[sub block meanflex strflex]
dlmwrite('/danl/Harmon_dynCon/flexdata.csv',flexdata) 

%get flexibility scores for each ROI for each run for whole-brain search
%can prob do this more effeciently but this is easier to see, harder to botch
flex_allrois=[];
for i=1:size(flex_cat,3)
  flex_allrois=[flex_allrois;flex_cat(:,:,i)'];
end

dlmwrite('/danl/Harmon_dynCon/flex_allrois.csv',flex_allrois) 
```

forRest 
```.matlab

%load data and concatenate flexibility statistics
[a,b]=system('ls -d /danl/Harmon_dynCon/7*/Rest/Rest1.feat/flex_rest.mat');
%do separately above and below 


c=strread(b,'%s');
flex_cat_rest=[];
for j=1:size(c,1)
    load(char(c(j)))
    flex_cat_rest=cat(3,flex_cat_rest,flex_rest)
end
plot(squeeze(mean(flex_cat_rest)))

block=repmat([1:4]',25,1);
sub=repmat([1:25]',1,4)'
sub=sub(:);

%reshape whole-brain average flexibility
meanflex=squeeze(mean(flex_cat_rest));
meanflex=meanflex(:);

%get striatal average flexibility
%check to make sure indices are correct
[trash,roi_names]=system('ls  /danl/Harmon_dynCon/Harvard-Oxford_ROIs/*nii.gz | xargs -n1 basename');
roi_names=strread(roi_names,'%s');
str_ind=[49,51,54,104,106,109];
roi_names(str_ind)

strflex=squeeze(mean(flex_cat_rest(str_ind,:,:)));
strflex=strflex(:);

plot(squeeze(mean(flex_cat(str_ind,:,:))))

%write out csv for modeling in R
flexdata=[sub block meanflex strflex]
dlmwrite('/danl/Harmon_dynCon/flexdataRest.csv',flexdata) 

%get flexibility scores for each ROI for each run for whole-brain search
%can prob do this more effeciently but this is easier to see, harder to botch
flex_allrois=[];
for i=1:size(flex_cat_rest,3)
  flex_allrois=[flex_allrois;flex_cat_rest(:,:,i)’];
end

dlmwrite('/danl/Harmon_dynCon/flex_allroisRest.csv',flex_allrois) 


```


### Combine raw reinforcement learning data from individual subjects 

```.matlab
a=dir('/danl/Harmon_dynCon/behavior/7*_tb_l*_svlo*.mat');

a={a.name};

longform=[];
for k = 1:length(a)
	k
	clear subnum chose_right shown_corr trial
	filea=a{(k)};
	load(a{(k)});	
    subnum=str2double(num2str(a{k}(1:3)))';
    chose_right=double(strcmp(resp,'b'))';
    shown_corr=double(shown_corr)';
    no_resp=shown_corr==2;
    chose_right(no_resp)=NaN;
    shown_corr(no_resp)=NaN;
    ntrials=size(resp',1);
    trial=1:ntrials';
    longform=[longform;repmat(subnum,ntrials,1) stim_shown' trial' chose_right shown_corr];
end
header={'sub','stim','trial','choice','fb'};
header2=sprintf('%s,',header{:});header2(end)=[];
dlmwrite('/danl/Harmon_dynCon/behavior/choice_fb_long.csv',...
	header2,'')

dlmwrite('/danl/Harmon_dynCon/behavior/choice_fb_long.csv',...
	longform,'-append','delimiter',',')


#do not run 
!cp /danl/Harmon_dynCon/behavior/choice_fb_long.csv ~/GitHub/rl_flexibility/RL_model/choice_fb_long.csv
```


### Prepare data for stan and run hierarchical bayesian model

```.r

#This code is also found in ~/Columbia/learninglab/dynCon/ch_scripts/RL_model/heirarBeysModel.R

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
dat <- read.csv('~/Columbia/learninglab/dynCon/ch_scripts/RL_model/choice_fb_long.csv')

choices<-unique(na.omit(dat$choice))

dat$chosen<-dat$choice
dat$chosen[dat$choice==choices[1]]=1
dat$chosen[dat$choice==choices[2]]=2
dat$unchosen=abs(dat$chosen-3)
dat$chose_two<-dat$chosen==2

subs = unique(dat$sub);
NS = length(subs);
NStim=length(unique(dat$stim));
MT=max(dat$trial);
NT = array(0,NS);
choice = array(0,c(NS,MT));
unchoice=choice;
choice_two=choice;
rew = choice;
stim=choice;

for (i in 1:NS) {
  NT[i] = nrow(subset(dat,sub==subs[i]));
  
  stim[i,1:NT[i]] = subset(dat,sub==subs[i])$stim;
  
  #choice and reward history
  choice[i,1:NT[i]] = subset(dat,sub==subs[i])$chosen;
  unchoice[i,1:NT[i]] = subset(dat,sub==subs[i])$unchosen;
  
  rew[i,1:NT[i]] = subset(dat,sub==subs[i])$fb;
  
  #based on choosing second option
  choice_two[i,1:NT[i]] = subset(dat,sub==subs[i])$chose_two;
}

choice[is.na(choice)]<--1
unchoice[is.na(unchoice)]<--1
rew[is.na(rew)]<--1
choice_two[is.na(choice_two)]<--1

rl_standata = list(NS=NS, NC=2, NStim=NStim, MT=MT, NT= NT, choice=choice, 
                   stim=stim,choice_two=choice_two,rew=rew );

rl_fit <- stan(file = '~/Columbia/learninglab/dynCon/ch_scripts/RL_model/standard_rl.stan', 
               data = rl_standata, iter = 1250, warmup = 250, chains = 4)

save(rl_fit,file='~/Columbia/learninglab/dynCon/ch_scripts/RL_model/rl_fit')
```


### Extract parameters from model
```.r
#Code can be found in ~/Columbia/learninglab/dynCon/ch_scripts/RL_model/extraxtParamsFromModel.R

library(rstan)
#install.packages('ppcor')
library(ppcor)

fit_rl<-load('~/Columbia/learninglab/dynCon/ch_scripts/RL_model/rl_fit')
fit_extract<-extract(rl_fit,permute=T)

#fit_extract$beta_mean<-fit_extract$b1/fit_extract$b2 #<- dividing is wrong #Group distribution - beta estimate mean for each posterior sample (stan gamma distribution) **Should be multiplied not divided!! because inverse. STAN uses rate - invesre of scale - rather than scale. THe inverse mean formula matters based on inverse   
#beta = inverse temperature 
hist(fit_extract$beta_mean) #model uncertainty about the mean 
#could also calculate SD 
fit_extract$alpha_mean<-fit_extract$a1/(fit_extract$a1+fit_extract$a2)
#alpha learning rate 
hist(fit_extract$alpha_mean)

#variance
hist(fit_extract$b1*fit_extract$b2^2)

#When parameters dont change in the paradigm- inverese temp is similar to percent correct - lower alpha better learning, higher inverse temp, better learning 

#Check with Juliet's paper similar values 
#with unstable estimes - first thing to try is to make hyperparameters - non informed priors - may be worth shrinking the dist. = making priors find more extreme values less plausble - and if the models is still fitting well then maybe thats meaningul -keep priors the same if comparing the two 

#beta how much does value influence decision 



betas<-apply(fit_extract$beta,2,mean)
alphas<-apply(fit_extract$alpha,2,mean)

for (i in 1:25){
  plot(fit_extract$alpha[,i], fit_extract$beta[,i], xlab="alpha", ylab='beta', main=paste('Subj', i))
}


#Check out Gamma distribution calculations for fits - 
#https://en.wikipedia.org/wiki/Gamma_distribution

fit_extract$beta_mean<-fit_extract$b1*fit_extract$b2
fit_extract$alpha_mean<-fit_extract$a1/(fit_extract$a1+fit_extract$a2)

str_acorrel<-matrix(0,4000,6)
str_bcorrel<-matrix(0,4000,6)
str_apbcorrel<-matrix(0,4000,6)
str_bpacorrel<-matrix(0,4000,6)
wb_acorrel<-NULL
wb_bcorrel<-NULL
wb_apbcorrel<-NULL
wb_bpacorrel<-NULL


flexdat<-read.csv("flex_allrois.csv",header=F)
flex_behav<-read.csv("flex_behav.csv",header=T)

View(flex_behav)

flexdat$Sub<-rep(seq(1,25,1),each=4)
flexdat$Block<-rep(seq(1,4,1),times=25)
flexdat$Corr<-flex_behav$correct
flexdat$weights<-flex_behav$weights
flexdat<-melt(flexdat,id=c("Sub","Block","Corr","weights"))
names(flexdat)[c(5,6)]<-c("ROI","flex")
flexdat$ROI<-as.factor(as.numeric(flexdat$ROI))
meanflex_rois<-tapply(flexdat$flex,list(flexdat$Sub,flexdat$ROI),mean)
meanflex<-rowMeans(meanflex_rois)

str_ind=c(49,51,54,104,106,109);

for( i in seq(1,dim(fit_extract$beta)[1],1)){
  
  betas_tmp=fit_extract$beta[i,]
  alphas_tmp=fit_extract$alpha[i,]
  k<-1
  
  wb_acorrel[i]<-cor(meanflex,alphas_tmp)
  wb_bcorrel[i]<-cor(meanflex,betas_tmp)
  wb_apbcorrel[i]<-pcor.test(meanflex,alphas_tmp,betas_tmp)[1] #partial correlations
  wb_bpacorrel[i]<-pcor.test(meanflex,betas_tmp,alphas_tmp)[1] #partial correlations 
  
  for(j in str_ind){
    str_acorrel[i,k]<-cor(meanflex_rois[,k],alphas_tmp)
    str_bcorrel[i,k]<-cor(meanflex_rois[,k],betas_tmp)
    
    str_bpacorrel[i,k]<-pcor.test(meanflex_rois[,k],betas_tmp,alphas_tmp)$estimate
    str_apbcorrel[i,k]<-pcor.test(meanflex_rois[,k],alphas_tmp,betas_tmp)$estimate
    k<-k+1
  }
}

par(mfrow=c(1,2))

plot(rowMeans(str_bcorrel),rowMeans(str_acorrel),
     col=rgb(0,0,0,alpha=0),pch=21,bg=rgb(0,0,0,alpha=.03),
     xlab="Correlation with Beta",ylab="Correlation with Alpha",
     main="Striatum Flexibility",xlim=c(-.1,.5))
abline(v=0,lty=2)
abline(h=0,lty=2)

plot(wb_bcorrel,wb_acorrel,
     col=rgb(0,0,0,alpha=0),pch=21,bg=rgb(0,0,0,alpha=.03),
     xlab="Correlation with Beta",ylab="Correlation with Alpha",
     main="Whole-brain Flexibility",xlim=c(-.1,.5))
abline(v=0,lty=2)
abline(h=0,lty=2)



par(mfrow=c(1,1))
plot(fit_extract$alpha[,1], fit_extract$beta[,1], xlab="alpha", ylab='beta', main='Subj 1')
plot(fit_extract$alpha[,2], fit_extract$beta[,2], xlab="alpha", ylab='beta', main='Subj 2')
plot(fit_extract$alpha[,3], fit_extract$beta[,3], xlab="alpha", ylab='beta', main='Subj 3')

plot(fit_extract$alpha[,1], fit_extract$beta[,1], xlab="alpha", ylab='beta', main='Subj 1')

```

### ML and fully Bayesian hierarchical models 
Estimate the effect of striatal and whole-brain flexibility on reinforcement learning. See models.Rmd and models.pdf for more details. R

```.r
#This code can also be found ~/Columbia/learninglab/dynCon/ch_scripts/RL_model/MLBayesianHeirarchicalModels.R

#prepare data for binomial logistic modelling
library(lme4)
library(reshape2)
library(brms)

flexdat<-read.csv("flex_allrois.csv",header=F)

#read in trial-by-by trial behavioral data 
data<-read.delim('behav_data.tsv',header=1)

data<-read.csv('behav_data.csv',header=1)
data <- data[data$group==7,]

#Only have behavioral data : no brain data 
data <- data[!data$subjectNum==716,]
data <- data[!data$subjectNum==718,]


data$block<-rep(rep(seq(1,4,1),each=30),25)

flexdat$subject<-rep(seq(1,25,1),each=4)
flexdat$Block<-rep(seq(1,4,1),times=25)
flexdat$Corr<-melt(tapply(data$optCor,list(data$block,data$subjectNum),mean,na.rm=1))$value

nacount<-is.na(data$optCor)
weights<-melt(tapply(nacount,list(data$block,data$subjectNum),sum))
flexdat$weights<-30-weights[,3]

str_ind=c(49,51,54,104,106,109);

flexdat<-melt(flexdat,id=c("subject","Block","Corr","weights"))

names(flexdat)[c(5,6)]<-c("ROI","flex")

flexdat$ROI<-as.factor(as.numeric(flexdat$ROI))


flexdat_str<-subset(flexdat,ROI %in% str_ind)
flexdat_str$ROI<-as.factor(flexdat_str$ROI)

flexdat_str_avg<-melt(tapply(flexdat_str$flex,list(flexdat_str$subject,flexdat_str$ROI),mean,na.rm=1))
names(flexdat_str_avg)<-c("subject","ROI","flex")
flexdat_str_avg<-subset(flexdat_str_avg,ROI %in% str_ind)

flexdat_str$meanflex<-rep(flexdat_str_avg$flex,each=4)

flexdat_str_melt<-melt(tapply(flexdat_str$flex,list(flexdat_str$subject,flexdat_str$Block),mean))
names(flexdat_str_melt)<-c("subject","Block","flex")



flexdat_str_melt$correct<-melt(tapply(flexdat_str$Corr,list(flexdat_str$subject,flexdat_str$Block),mean))$value
flexdat_str_melt$weights<-melt(tapply(flexdat_str$weights,list(flexdat_str$subject,flexdat_str$Block),mean))$value



m_str<-glmer(correct~flex+(flex||subject),
             family=binomial,weights=weights,
             data=flexdat_str_melt)

write.csv(flexdat_str_melt,'flex_behav.csv')

#for posterior inference, run bayesian model using brms wrapper for stan
options(mc.cores = parallel::detectCores())
flexdat_str_melt$numcorr<-as.integer(flexdat_str_melt$correct*flexdat_str_melt$weights)
mlearn_stan<-brm(numcorr~flex+(flex|subject),data=flexdat_str_melt,family=binomial)

##
Warning messages:
  1: There were 15 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
2: Examine the pairs() plot to diagnose sampling problems


```

### Whole-brain search for effects of flexibility on learning
R

```.r
#This code can also be found ~/Columbia/learninglab/dynCon/ch_scripts/RL_model/wholeBrainSearchforEffectsofFlexonLearning.R

#load in data
library(lme4)
roi_data<-read.csv('flex_allrois.csv',header=0)
behav<-read.csv('flex_behav.csv')

p<-NULL
b<-NULL


for(i in 1:dim(roi_data)[2]){
  mtmp<-summary(glmer(behav$correct~roi_data[,i]+(roi_data[,i]||behav$subject),family=binomial,weights=behav$weights))
  if(dim(summary(mtmp)[10]$coefficients)[2]>=4){
    p[i]<-summary(mtmp)[10]$coefficients[2,4]
    b[i]<-summary(mtmp)$coefficients[2,1]
  } else{
    p[i]<-NA
    b[i]<-NA
  }
}
write(p,'p_vals_learning_glm.csv')
write(b,'b_vals_learning_glm.csv')

##
Warning message:
In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  Model failed to converge with max|grad| = 0.00105927 (tol = 0.001, component 1)

```

### Combine module allegiance matrices for striatum ROIs across subjects
Matlab

```.matlab
str_ind=[49,51,54,104,106,109];
[a,b]=system('ls -d /data/engine/rgerraty/learn_dyncon/4*/a_mat.mat');
c=strread(b,'%s');

h=1
for s=1:length(c)
  clear a_mat clear a_mat_str
  load(c{s})
  a_mat(a_mat==1)=NaN;
  a_mat_str=a_mat(:,str_ind,:); 
  for i=1:size(a_mat_str,1)
    for j=1:size(a_mat_str,2)
      for k=1:size(a_mat_str,3)
        a_mat_long(h,:)=[a_mat_str(i,j,k),i,j,ceil(k/8),s];
        h=h+1;
      end
    end
  end
end

header={'allegiance','ROI','str_roi','block','sub'};
header2=sprintf('%s,',header{:});header2(end)=[];
dlmwrite('/data/engine/rgerraty/learn_dyncon/alleg_long.csv',...
  header2,'')

dlmwrite('/data/engine/rgerraty/learn_dyncon/alleg_long.csv',...
  a_mat_long,'-append','delimiter',',')

```


### Assessing motion 
Adolescents 
```.matlab
%read in all subject/run mc motion files 
[a,b]=system('ls -d /danl/Harmon_dynCon/7*/Learn?_PEpriorD.feat/mc/prefiltered_func_data_mcf_rel_mean.rms');

%do separately above and below 
c=sort(strread(b,'%s'))
motion = (1:length(c))'

for i=1:length(c)
	load(c{i})
	motion(i) = prefiltered_func_data_mcf_rel_mean
end 

 dlmwrite('/danl/Harmon_dynCon/motion.csv', motion)
```

Adults 
```.matlab
%read in all subject/run mc motion files 
[a,b]=system('ls -d /data/engine/rgerraty/learn_dyncon/4*/Learn?_PEprior.feat/mc/prefiltered_func_data_mcf_rel_mean.rms');

%do separately above and below 
c=sort(strread(b,'%s'))
motion = (1:length(c))'

for i=1:length(c)
	load(c{i})
	motion(i) = prefiltered_func_data_mcf_rel_mean
end 

 dlmwrite('/data/engine/charmon/motion.csv', motion)
```


