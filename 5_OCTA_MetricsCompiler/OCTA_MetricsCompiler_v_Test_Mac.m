% OCTA_MetricsCompiler
% By: Mina Gaffney
% Created on: 3/10/2021
%
% This script was written to Compile the OCTA metrics output from the
% Vessel Loss script. 
% 
% Expected input: VesselLoss_Output folder 
% Output: 1 spreadsheet per subject which contains the descriptive stats for
% the FAZ and PICA from each simulated vessel loss folder. 
%
% Options: None 
% This has some syntax for Mac vs PC - lines 45, 93, 101
% THIS IS ONLY TO BE RUN ON THE TEST DATASET
% JC_12058, 1, 2, 3 % loss with 5 iterations in each folder - TESTING only
% (different machines)

%% Load in files or folders
input_dir = uigetdir();
% Find the first level of folders within the chosen directory:
searchlist_lvl1 = dir(input_dir);

% Loading in and reading the ground truth data (0pct loss)
groundtruth_struct = searchlist_lvl1(contains({searchlist_lvl1.name}, {'AllSubjects_0pct'}));
groundtruth = readcell(fullfile(groundtruth_struct.folder,groundtruth_struct.name));

%Calculate average of ground truth for each metric
groundtruth_avg = mean(cell2mat(groundtruth(2:end,2:end)),1);% Calculating average of the population 
groundtruth_sd = std(cell2mat(groundtruth(2:end,2:end)),0);%Calculating standard deviations of the population, note differences in stdev between MatLab and Excel (N or N-1)

% Fix the DOS era issue with the dir function (loads in the parent
% directories '.' and '..')
searchlist_lvl1 = searchlist_lvl1(~ismember({searchlist_lvl1.name}, {'.', '..'}));

% Ignore Lut or anything that isn't a folder
searchlist_lvl1 = searchlist_lvl1([searchlist_lvl1.isdir]);

% Pre- allocate cell array:
level1_folders = cell(length(searchlist_lvl1),1);

% Loop through level 1 folders in directory:
for iii = 1:length(searchlist_lvl1)
    level1_folders{iii} = searchlist_lvl1(iii).name;
    lvl1_paths{iii,1} = fullfile(input_dir,level1_folders{iii});
      
% Finds the full path for all folders and subfolders in selected folder:
searchlist = genpath(lvl1_paths{iii});
searchlist = strsplit(searchlist, ':'); %: for Mac, ; for PC

% Read in the stuff in the current folder
Current_subfolder = dir(searchlist{1,1}); 
Current_subfolder = Current_subfolder(~ismember({Current_subfolder.name}, {'.', '..'}));
Current_subfolder = Current_subfolder([Current_subfolder.isdir]);
Current_subfolder = struct2cell(Current_subfolder)'; 

% Loop through Percent Loss folders
for iv = 1:size(Current_subfolder,1)
    
Current_lossfolder = dir(fullfile(Current_subfolder{iv,2},Current_subfolder{iv,1})); 
Current_lossfolder = Current_lossfolder(~ismember({Current_lossfolder.name}, {'.', '..'}));
Current_lossfolder = struct2cell(Current_lossfolder)'; 
    
% Load in the RealLoss and the RawData spreadsheets
% RealLoss_dir(iv,:) = Current_lossfolder(~cellfun(@isempty, strfind(Current_lossfolder(:,1), 'RealLoss')),:);
% RealLoss_sheet(:,:,iv) = readcell(fullfile(RealLoss_dir{iv,2},RealLoss_dir{iv,1}));
RawData_dir(iv,:) = Current_lossfolder(~cellfun(@isempty, strfind(Current_lossfolder(:,1), 'Raw_Data')),:);
RawData_sheet(:,:,iv) = readcell(fullfile(RawData_dir{iv,2},RawData_dir{iv,1}));

%% Determine Subject ID 
FileNameSplit = strsplit(RawData_dir{iv,1},'_'); 
SubjectID = [FileNameSplit{1,1} '_' FileNameSplit{1,2}]; 


%% Determine % loss
PercentLoss(iv,1) = str2num(strrep(FileNameSplit{1,3},'pct',''));

end

%% Calculate stats
Metrics_Summary_data = zeros(size(PercentLoss,1),42);

%% Find subject ID within groundtruth 

ID_check(1,:) = groundtruth(contains(groundtruth(:,1),string(SubjectID)),:);
Subject_groundtruth = cell2mat(ID_check(1,2:14));

Subject_Change = groundtruth_avg - Subject_groundtruth; %Calculating the amount of adjustment need for the one subject
Norm_groundtruth = Subject_groundtruth + Subject_Change; %Adjusting the subject's no-loss metrics to the population average

for v = 1:size(PercentLoss,1)

NormData(:,:,v) = cell2mat(RawData_sheet(2:end,2:end,v))+Subject_Change; %normalizing each person's data the average of the population without any simulated loss

normname = RawData_dir{v,1}; %Pulling correct name from the directory of file names
norm_name = strcat(normname(1:end-12),'Normalized_Data.csv'); %Changing name to show z-score calculations
writematrix(NormData(:,:,v), fullfile(RawData_dir{v,2},'/',norm_name)); %Saving z-score matrix to an excel sheet \ for PC and / for Mac

for m = 1:5  %Previous code: m = 1:length(NormData)
    zscore(m,:,v) = (NormData(m,:,v)-groundtruth_avg)./groundtruth_sd; %calculating the zscore for each cell, row by row in the normalized data
end

zname = RawData_dir{v,1}; %Pulling correct name from the directory of file names
zscore_name = strcat(zname(1:end-12),'Zscores.csv'); %Changing name to show z-score calculations
writematrix(zscore(:,:,v), fullfile(RawData_dir{v,2},'/',zscore_name)); %Saving z-score matrix to an excel sheet \ for PC and / for Mac

zlogical = zeros(5,13); %Preallocating the 1st matrix for turning the z-scores that are signficiant to 1 and the rest to 0
for n = 1:13
    for l = 1:5 
        if zscore(l,n,v) < -2 %if the z-score is less than below 2SDs,, make it a 1
            zlogical(l,n) = 1;
        elseif zscore(l,n,v) > 2 %if the z-score is above 2SDs, make it a 1
            zlogical(l,n) = 1;
        end
    end
end

ground_truth_summary(1,1:5) = 0;
ground_truth_summary(1,6:3:44) = Norm_groundtruth(1,:);
ground_truth_summary(1,43:44) = 0; 

% Calculate all descriptive stats by going down columns: 
Metrics_Summary_data(v,1) = PercentLoss(v,1); % Intended Loss
% Metrics_Summary_data(v,2) = mean(cell2mat(RealLoss_sheet(:,:,v)),'all'); %Real Loss
% Metrics_Summary_data(v,3) = std(cell2mat(RealLoss_sheet(:,:,v))); %Real Loss SD
% Metrics_Summary_data(v,4) = min(cell2mat(RealLoss_sheet(:,:,v))); %Calculating min of Real Loss
% Metrics_Summary_data(v,5) = max(cell2mat(RealLoss_sheet(:,:,v))); %Maxs
Metrics_Summary_data(v,6:3:44) = mean(NormData(:,:,v),1); % Means
Metrics_Summary_data(v,7:3:44) = std(NormData(:,:,v),0,1); %SDs
Metrics_Summary_data(v,8:3:44) = sum(zlogical(:,:),1); %summing up z-score
%Metrics_Summary_data(v,6:4:55) = min(cell2mat(NormData(:,:,v)),[],1); % Mins
%Metrics_Summary_data(v,7:4:55) = max(cell2mat(NormData(:,:,v)),[],1); % Max

%% Outputs 
Metrics_Summary_header0 = {' ',' ', ' ',' ',' ',' ',' ',' ',' ',' ',' ','With FAZ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','Without FAZ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',};
Metrics_Summary_header1 = {' ',' ',' ',' ',' ','FAZ Area',' ',' ','Vessel Density',' ',' ', 'Average PICA',' ',' ', 'Summed PICA',' ',' ', 'PICA Regularity (mean/sd)',' ',' ', 'PICA SD',' ',' '...
    'Average PICA',' ',' ','Summed PICA',' ',' ','PICA Regularity (mean/sd)',' ',' ','PICA SD',' ',' ','VCI',' ',' ','VPI',' ',' ','Df',' ',' '};
Metrics_Summary_header2 = {'Intended Loss', 'Real Loss', 'Real Loss SD', 'Real Loss Min','Real Loss Max','Mean', 'SD', 'Zscore', 'Mean', 'SD', 'Zscore', 'Mean', 'SD', 'Zscore','Mean', 'SD', 'Zscore', 'Mean', 'SD', 'Zscore', 'Mean', 'SD', 'Zscore',...
   'Mean', 'SD', 'Zscore', 'Mean', 'SD', 'Zscore', 'Mean', 'SD', 'Zscore', 'Mean', 'SD', 'Zscore', 'Mean', 'SD', 'Zscore', 'Mean', 'SD', 'Zscore', 'Mean', 'SD', 'Zscore'};
Metrics_Summary = cat(1,Metrics_Summary_header0,Metrics_Summary_header1,Metrics_Summary_header2,num2cell(ground_truth_summary),num2cell(Metrics_Summary_data));

OutName = [SubjectID, '_Metrics_Summary.xlsx']; 
writecell(Metrics_Summary, fullfile(lvl1_paths{iii,1},OutName));

clear Subject_groundtruth;
 end
clear RealLoss_sheet RawData_Sheet ground_truth_summary Metrics_Summary_data;
end



