%% I. Define paths

% % % the directory where the raw data niftis are located
inDir='';

% % % the directory where the analysis output will be saved
analysisDir='';

% path to a T1w image to that may be used as refernce for alignment of the anlysis.
reffile= fullfile(inDir,'.nii.gz');

%% II. run mrQ
% mrQ will create qT1 maps and MTV maps
% the code is here: https://github.com/mezera/mrQ
% the analysis requires variable flip angle SPGR data, and variable
% inversion time SEIR data (to ermove B1 bias)

% the MTV map can be used to calculate g-ratio, as done in:
% Evaluating g-ratio weighted changes in the corpus callosum as a function of age and sexReview Article
% In Press, Accepted Manuscript, Available online 30 June 2017
% Shai Berman, Kathryn West, Mark D. Does, Jason D. Yeatman, Aviv A. Mezer

mrQoutDir=fullfile(analysisDir,'mrQ');

if exist(reffile,'file')
    mrQ_run(inDir,mrQoutDir,[],[],[],{'ref',reffile});
else
    mrQ_run(inDir,mrQoutDir);
    disp('no ref file was found, so mrQ ran w/o one')
end


%% III. Run DTIinit

% DTIinit will perform the preprocessing for the diffusion data. it will
% create align the different diffusion weighted images, it can perform
% eddy-current correction, and it fits a tensor to the data set with the
% largest number of directions. the code is here:
% https://github.com/vistalab/vistasoft

% set the reference T1 file

mrQstruct_path=fullfile(mrQoutDir, 'mrQ_params.mat');
load(mrQstruct_path);
if ~exist(reffile,'file')
    reffile =fullfile(mrQ.T1w_files,'T1w.nii.gz');
end

% determnine output dir
DTIoutputDir= fullfile(analysisDir,'DTIoutput');

% find the diffusion dataset with the largest number of directions
dwiFile = getDwiFilesStruct(inDir);
if numel(dwiFile) > 1
    msg = 'Multiple diffusion datasets were found! Choosing the scan with the largest number of directions!';
    for j = 1:numel(dwiFile)
        dirs(j) = dwiFile{j}.directions;
    end
    [~, loc] = max(dirs);
    dwiData = dwiFile{loc};
else
    dwiData  = dwiFile{numel(dwiFile)};
end

% These are the parameters we used for the weston havens DWI data
params                      = dtiInitParams;
params.eddyCorrect          = 0; % we used a sequence that minimzes eddy currents
params.noiseCalcMethod      = 'b0';
params.numBootStrapSamples  = 0;
params.clobber              = -1;
params.fitMethod            = 'rt';

% This is where we'll save all the outputs
params.dt6BaseName = fullfile(DTIoutputDir,'DTI',['dt' num2str(dwiData.directions)]);

% Get a path to the dwi file.
dwiFile          = dwiData.nifti;
params.bvecsFile = dwiData.bvec;
params.bvalsFile = dwiData.bval;

% Run dtiInit
fprintf('Running dtiInit: %s \n',dwiFile);
dtPath = dtiInit(dwiFile, reffile,params);

%% IV. Align mrQ maps to DTI
% Even though they are aligned to the same image, diffusion data is
% acquired often with EPI and therefore has warping, and we can use ANts ti
% warp the maps to the diffusion space. This can be avoided when, for
% example, the diffusion data has been unwarped (lik with fsl top-up).

% The first map in the cell array is the t1 - which is the first input
% argument to the register funciton, so we skip that one here.

maps = fieldnames(mrQ.maps); maps = maps(2:end);
otherMapsPath = cell(1,numel(maps));
for ii = 1:numel(maps)
    otherMapsPath{ii} = [mrQ.maps.(maps{ii})];
end
mapsdir = fullfile(params.dt6BaseName,'bin');
b0      = fullfile(mapsdir,'b0.nii.gz');
t1      =  mrQ.maps.T1path;

% Perform the alignment of the maps to the dti data (b0)
disp('Aligning maps to diffusion data...');
AlignedMaps= mrQ_registerMap2DTI(b0,t1,otherMapsPath,mapsdir);


%% V. probablistic tractography with mrTrix
% This setion perfornm probabilistic tractography with the CSD model. the
% code is also distrubuted in vista soft

outDirTrac=fullfile(analysisDir,'mrTrix');
if ~exist(outDirTrac,'dir'),  mkdir(outDirTrac);   end

% load dt file
dt=dtiLoadDt6(dt6Path);

% Lmax defines the number of set of spherical harmonics to fit CSD model.
% In vistasoft-MrTrix bridge, the default value of Lmax is 6.
lmax=6 ;

% Fit the spherical deconvolution (CSD) model to DWI data.
wmMaskFile=dt.files.wmMask;
files = mrtrix_init(dt6Path, lmax, outDirTrac, wmMaskFile);
 %This functionn performs the following operations:
%   1. Convert the raw dwi file into .mif format
%   2. Convert the bvecs, bvals into .b format
%   3. Convert the brain-mask to .mif format 
%   4. Fit DTI and calculate FA and EV images
%   5. Estimate response function for single fibers, in voxels w/ FA > 0.7
%   6. Fit the CSD model. 
%   7. Convert the white-matter mask to .mif format. 

% define parameters for tracking
roi=files.wm;  % the ROI in which  to place the seeds
mask= files.wm;
mode='prob';   % Tracking mode:  {'prob' | 'stream' | 'tensor'} for probabilistic 
               % or  deterministic-SD or determinstic DT tracking. 
nSeeds=500000; % The number of fibers to generate. ~10min For 100k streamlines
bkgrnd=false;  % on unix, whether to perform the operation in another process
verbose=false; %  whether to display standard output to the command window.
clobber=false; % Whether or not to overwrite the fiber group if it was already computed

% perform tractography
[status, results, fg, pbdFile] = mrtrix_track(files, roi, mask, mode,...
    nSeeds, bkgrnd, verbose, clobber);
% Provided a csd estimate, generate estimates of the fibers starting in roi 
% and terminating when they reach the boundary of mask

% convert to a format that is comfortable for AFQ
% fg is in pdb format
% maybe just convert to mat
fgDetails = [pbdFile(1:end-3)];
fgPath = [fileparts(fgDetails),'/Tractography.mat'];
save(fgPath, 'fg');



%% fiber segmentation with AFQ

% After the data was analyzed with dtiinit, mrQ (+alignment to DTI space)
% and mrTrix tractography, we perform fiber tract segments, including the
% 20 mori group fibers as well as the 8 major callosal tracts. 
% we do that with AFQ, which is available here: 
% https://github.com/yeatmanlab/AFQ

% initializt the AFQ atruct using dtpath, we do this for each subject
% separately (so we can do this in parallel, and then we can combine the
% different AFQ structures into one big structure maybe we need to separate
% functions and "find" the correct directories
%  we need to fimd dtPaths and set it up in the correct struct. 
% dtPath is thethe path to the dt6.mat file

% you can define the final output however you find comfortable. most of the
% output will be in the DTIoutput directory anyway, and the final output
% only determines the ouput of the AFQ struct (mat file). I believe it's
% comfotable to place it in the same directory as the tractography from
% whch it originated - 
outDirAFQ = outDirTrac;

sub_dirs{1} = fileparts(char(dtPath));

% Just one group here
sub_group = ones(numel(sub_dirs),1);

% Create afq structure
afq = AFQ_Create('sub_dirs',sub_dirs,'sub_group',sub_group,'clip2rois',0,'outdir',outDirAFQ);
% since we only performed tractography - make sure AFQ knows about it,
% otherwise it will perform determinstic tractography
afq=AFQ_set(afq,'wholebrainfgpath',fgPath); 


% Get and set the fieldnames for those maps we want to analyze
% with AFQ. use any mrQ (the lated one that was loaded)

for jj = 1:numel(alignedMaps)
    for ii=1:length(sub_group)
        %  This needs to be a cell array, s.t. each line (ii) is a
        %  different subjects, each row (jj) is a different map. The
        %  afq.files.images.paths should be organized such that it has a
        %  cell for each field (TV,T1) and within each such cell, there
        %  should be cells containing the paths for the fitting maps of all
        %  subjects, different cell for each subject....
        image{ii} = alignedMaps{jj};
    end
    afq = AFQ_set(afq, 'images', image);
end

% RUN AFQ
disp('Running AFQ...');
afq = AFQ_run(sub_dirs, sub_group, afq);
afq = AFQ_SegmentCallosum(afq);

