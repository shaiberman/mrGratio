function [weightsPath]= createWeightedROIs(afq_path, dt6_path, segPath, outdir)

% This function creates a mat file containing the probability of each voxel
% in the CC to "belong" to any of the 8 major Callosal tracts. 
% As input it need the path to:
%         dt6 file created with mrDiffusion's dtiInit.
%             the dt6 file is used to create a rough callosum ROI which we
%             will refine
%         afq file create with afq.
%             the afq contains pointers to the 8 callosel fiber groups that
%             will be used to refine the sub-regions of the corpus
%             callosum.
%         nifti with a CC segmentation created with createROI
% The output is the path to the weights file. 
%         the size of the weights file is the size of our brain (3D), with
%         a 4th dimension that has 8 places, one for each fiber tract. for
%         each vooxel in each tract, we get the percent of fibers in that
%         voxell to belong to the certain tract. For examply, if in voxel
%         in the image space location [x,y,z] has 20 fibers from tract 2,
%         and 80 fibers that belong to tract 3, weill find that:
%         weightsMat_prob(x,y,z,1)=0;
%         weightsMat_prob(x,y,z,2)=0.2;
%         weightsMat_prob(x,y,z,3)=0.8;
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

if notDefined('segPath')
    error('please provide path to CC segmentation nifti')
elseif ~exist(segPath,'file')
    error('CC segmentation nifti does not exist')
end

if notDefined('outdir')
    outdir=fileparts(afq_path);
end
%% step 1: get the xform to image space
dt=dtiLoadDt6(dt6_path);
xform=inv(dt.xformToAcpc);
x=xform(1,4);% the middle of the brain!
%% step 2: load the CC segmentation for the subject to create the mask.

mask=readFileNifti(segPath);mask=logical(mask.data);

% initialize the mask for the weights ffilw
probMask=zeros([size(mask),8]);
voxInd=find(mask);

%% step 3: get the fibers
load(afq_path)
load(afq.files.fibers.CC_Occipital_clean{1});fg_a{1}=fg;
load(afq.files.fibers.CC_Temporal_clean{1}); fg_a{2}=fg;
load(afq.files.fibers.CC_Post_Parietal_clean{1});fg_a{3}=fg;
load(afq.files.fibers.CC_Sup_Parietal_clean{1}); fg_a{4}=fg;
load(afq.files.fibers.CC_Motor_clean{1}); fg_a{5}=fg;
load(afq.files.fibers.CC_Sup_Frontal_clean{1}); fg_a{6}=fg;
load(afq.files.fibers.CC_Ant_Frontal_clean{1});fg_a{7}=fg;
load(afq.files.fibers.CC_Orb_Frontal_clean{1}); fg_a{8}=fg;

%% step 4: find the CC ROI of each crossing fibers

for ii=1:8  %    go over the 8 tracts
    
    MidFibInd=[];
    
    if isempty(fg_a{ii}.fibers)
        continue
    else
        
        for jj=1:length(fg_a{ii}.fibers)  %   go over the single fibers within the tract
            
            % get the coordinatd in image space
            fibCoords= ceil(mrAnatXformCoords(xform,fg_a{ii}.fibers{jj}));
            
            % take only mid coordinates
            locx=find(fibCoords(:,1)>=x-4 & fibCoords(:,1)<=x+4);
            fibCoords=unique(fibCoords(locx,:),'rows');
            
            % make coordinates into indices
            FibInd=sub2ind(size(mask), fibCoords(:,1), fibCoords(:,2), fibCoords(:,3));
            
            % track all indices of the fibers within a tract
            MidFibInd=[MidFibInd;FibInd];
        end
    end
    
    %%  step 5: find the amount of fibers crossing in each voxel from the ROI
        
    probMaskTemp=zeros(size(mask));
    for jj=1:length(voxInd)
        % For each voxel, find the number of fibers that pass through it
        % (AKA how many time do we find the index of the current voxel, in
        % the matrix containing the indices of all of the fibers)
        probMaskTemp(voxInd(jj))=sum(ismember(MidFibInd, voxInd(jj)));
    end
    
    % save the matrix:
    % the size of the matric is the size of our brain, with a 4th dimension
    % that has 8 places, such that for each fiber tract, we get for each
    % voxel the number of fibers fromtha tract that pass through it. 
    probMask(:,:,:,ii)=probMaskTemp;
    
end

%% step 6: make it into percentages

% first find the overall number of fibers that pass through each voxel 
meanFiberNumber=repmat(sum(probMask,4),1,1,1,8); 

% now divide by the overall number of fibers. 
weightsMat_prob = probMask./meanFiberNumber;

% get rid of nans
weightsMat_prob(isnan(probMask))=0;

% Now if voxel in the location [x,y,z] has only fibers from the first
% tract, we will find that weightsMat_prob(x,y,z,1)=1, and also,
% weightsMat_prob(x,y,z,2)=0. 

%% save
weightsPath=fullfile(outdir,'weights.mat');
save(weightsPath,'weightsMat_prob');