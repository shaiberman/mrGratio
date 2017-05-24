close all
clear
currentPath=fileparts(which(mfilename));

%% table of pvalue dependence on age

% load data
load(fullfile(currentPath,'pVAlues.mat'))

% 
for ii=1:length(fgNames)
    fgNames{ii}=regexprep( fgNames{ii},'-','');
    eval([fgNames{ii},'=pa(:,',num2str(ii),')'])   
end

for ii = 1:length(ages)
    ages_str{ii}=['7-',num2str(ages(ii)),' y'];
end
    
eval( [ 'tbl=table( ',fgNames{1},',',fgNames{2},',',fgNames{3},',',fgNames{4},...
        ',',fgNames{5},',',fgNames{6},',',fgNames{7},',',fgNames{8},')'])

tbl.Properties.DimensionNames = {'subRegion' 'age range'};

tbl.Properties.RowNames=ages_str

