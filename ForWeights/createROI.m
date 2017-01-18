
function [maskPath]= createROI(afq_path, dt6_path, outdir)

% This function creates an ROI of the Corpus Callosum, that is only made up
% of voxels that contain fibers of the 8 major callosal tracts/ 
% as input it need the path to:
%         dt6 file created with mrDiffusion's dtiInit.
%             the dt6 file is used to create a rough callosum ROI which we
%             will refine
%         afq file create with afq.
%             the afq contains pointers to the 8 callosel fiber groups that
%             will be used to refine the sub-regions of the corpus
%             callosum.
% The output for this function is the path to the CC roi saved as nifti
%%
%% check the input variable

if notDefined('afq_path')
    error('please provide path to AFQ file')
elseif ~exist(afq_path,'file')
    error('afq path  does not exist')
end

if notDefined('dt6_path')
    error('please provide path to dt6 file')
elseif ~exist(afq_path,'file')
    error('dt6 path does not exist')
end

if notDefined('outdir')
    outdir=fileparts(afq_path);
end
%% step 1: use dt6.mat to get a rough callosum ROI
% get the 'rough' CC ROI using mrDiffusion
dt=dtiLoadDt6(dt6_path); 
ccThresh = 0.5;
[ccCoords_acpc] = dtiFindCallosum(dt.dt6,dt.b0,dt.xformToAcpc,ccThresh);

% dtiFindCallosum provides coordinates in acpc space, we work in image space:
xform=inv(dt.xformToAcpc);
ccCoords =round(mrAnatXformCoords(xform, ccCoords_acpc));

%% step 2: Expand the CC mask
% create matrix with ones in x axis
xones = [ones(size(ccCoords,1),1), zeros(size(ccCoords,1),1), zeros(size(ccCoords,1),1)]; 
% add the matrix to the cc coordinates such that it's wider laterally. In
% our data the mrDiffusion was sslightly tilted so we needed to add more
% voxels on one side compared with the other (to make the mask medial),
% this migh tneed to change in oither data sets: 
ccCoords_expanded = [ccCoords; ccCoords+xones; ccCoords-xones; ccCoords-2*xones; ccCoords-3*xones];
ccCoords_expanded = unique(ccCoords_expanded,'rows');

% turn the coordiantes into a mask 
ccInd=sub2ind(size(b0),ccCoords_expanded(:,1),ccCoords_expanded(:,2),ccCoords_expanded(:,3));
roughmask=zeros(size(b0)); roughmask(ccInd)=1;

%% step 3: load the fiber group
load(afq_path)

load(afq.files.fibers.CC_Occipital_clean{1});fg_a{1}=fg;
load(afq.files.fibers.CC_Temporal_clean{1}); fg_a{2}=fg;
load(afq.files.fibers.CC_Post_Parietal_clean{1});fg_a{3}=fg;
load(afq.files.fibers.CC_Sup_Parietal_clean{1}); fg_a{4}=fg;
load(afq.files.fibers.CC_Motor_clean{1}); fg_a{5}=fg;
load(afq.files.fibers.CC_Sup_Frontal_clean{1}); fg_a{6}=fg;
load(afq.files.fibers.CC_Ant_Frontal_clean{1});fg_a{7}=fg;
load(afq.files.fibers.CC_Orb_Frontal_clean{1}); fg_a{8}=fg;
%
%% step 4: find the CC ROI of each crossing fibers
x=xform(1,4);% take the middle of the brain

CCtractCoords=[];

 % go over the 8 tracts
for ii=1:length(fg_a)  
   
    MidFibCoords=[];

    if isempty(fg_a{ii}.fibers)
        continue 
    else
        
        % go over the single fibers within the tract
        for jj=1:length(fg_a{ii}.fibers)  
            
            % get the coordinatd in image space
            fibCoords= round(mrAnatXformCoords(xform,fg_a{ii}.fibers{jj}));
            
            % take only mid coordinates: these are the coordinate for the
            % medial voxels in which these fibers pass through
            locx=find(fibCoords(:,1)>=x-4 & fibCoords(:,1)<=x+4);
            loc=locx;
            
            % the following variable keep tracks of the voxel coordinates
            % for all of the fibers fromthis tract
            MidFibCoords=[MidFibCoords;fibCoords(loc,:)];  
            
        end
        
        % remove excess coordinates
        MidFibCoords = unique(MidFibCoords,'rows');
        
        % the following variable keep tracks of the voxel coordinates
        % for all of the fibers and all tracts!
        CCtractCoords=[CCtractCoords; MidFibCoords,ii*ones(length(MidFibCoords),1)]; % these are the Actuall coordinates.
    
    end
end

%% step 5: make mask from coordinates

%  initialize masks
maskFib=zeros(size(b0));   CC_ROI=zeros(size(b0));

% switch to image space
ccIndFib=sub2ind(size(b0),CCtractCoords(:,1),CCtractCoords(:,2),CCtractCoords(:,3));
maskFib(ccIndFib)=1;

maskFinal = roughmask & maskFib;

%% step 6: save  CC segmentation nifti

maskPath=fullfile(outdir,'CCseg.nii.gz');
dtiWriteNiftiWrapper(single(maskFinal),dt.xformToAcpc,maskPath);
%
