close all 
clear all

%This script creates randomized vessel loss and removes the segments from
%an image
%Inputs:
    %Binary image with all vessel segments present
    %Skeletonized image for segment selection
    %Hardcoded but important: max percent vessel loss

%Outputs:
    %Binary image with randomized vessel loss
    
%Modified 1/29/22 by JC to add in Mac semantics on lines 108 & 111

%Adding support libraries to the Matlab search path (natsortfiles, lib, and Graph2Skel)
mfilefullpath = matlab.desktop.editor.getActiveFilename; 
mfileshortpath = erase(mfilefullpath, 'OCTA_VesselLoss_v5.m'); 
natsort_path = append(mfileshortpath, 'natsortfiles');
lib_path = append(mfileshortpath, 'lib');
addpath(genpath(natsort_path));
addpath(genpath(lib_path));
link_path = append(mfileshortpath, 'Graph2Skel');
addpath(genpath(link_path));

%Batch Upload feature:
%User selects desired folder:
root_path = uigetdir('.','Select directory containing images');
pathparts = strsplit(root_path,'/'); %split up file path parts \ = PC, / = Mac
top_path = fullfile(pathparts{1:2},'\Vessel_Loss_Output'); %Create Vessel Loss folder at the same level as the cropped images folder
mkdir(top_path); %Create Vessel Loss Output folder
addpath(genpath(top_path)); %Add Vessel Loss Output folder to matlab path

%Get percent loss file
[lutf, lutp] = uigetfile(fullfile(root_path,'*.csv'), 'Select percent loss file (csv)');
fid = fopen(fullfile(lutp, lutf), 'r'); %Open csv file
loss_range = cell2mat(textscan(fid, '%f','delimiter', ',')'); %Load in numbers as a float and change the cell array to a matrix
fclose(fid); %Close the csv file

tic
root_dir = dir(root_path); %Get all of the files specified in directory contianing images
root_dir = root_dir(~ismember({root_dir.name}, {'.', '..'})); %remove parent directories
root_dir = struct2cell(root_dir)'; %Change type of data from structure to cell for the image names


%looks for all tifs with binary
skel_dir = root_dir(~cellfun(@isempty, strfind(root_dir(:,1), 'VascMaskedSkel')),:); %Create a directory with only images in the root directory with the word skel in it
unmasked_dir = root_dir(~cellfun(@isempty, strfind(root_dir(:,1), 'VascUnmaskedSkel')),:);
mask_dir = root_dir(~cellfun(@isempty, strfind(root_dir(:,1), 'ArterMask')),:);

%Sort each directory
skel_dir = natsortfiles(skel_dir(:,1)); %Sort the skel directory by name
unmasked_dir = natsortfiles(unmasked_dir(:,1));
mask_dir = natsortfiles(mask_dir(:,1));
SE5 = strel('disk',5); %Sets constant for dilation

%Load in percent loss matrix
for i = 1:size(skel_dir,1) %For the number of images in the binary directory
    fileparts = strsplit(skel_dir{i},'_'); %splits up file name for each image in the skeleton directory
    subID = strcat(fileparts{1},'_',fileparts{2}); %Recreates subject ID from file name
    subj_path = fullfile(top_path,subID); %Creates a path within Vessel Loss Output using subject ID
    mkdir(subj_path); %makes subject folder
    addpath(genpath(subj_path)); %Adds subject folder to matlab path
    
    unmasked_skel = imread(fullfile(root_path,unmasked_dir{i}));
    skel = imread(fullfile(root_path,skel_dir{i})); %Load in the masked skeleton image 
    mask = imread(fullfile(root_path,mask_dir{i}));
    skelsize = size(skel); %Get size of the skeleton image
    skel = skel > 0; %Ensure skeleton is a logical type array
    skelpx = sum(sum(skel)); %Find total path length by getting sum of non-zero pixels

    
    %Obtain locations of segments within image
    [~, N, link] = Skel2Graph3D(skel, 0); %Gets all of the segments (and length of segments) in the skeletonized image
    numsegs = length(link); %Finds the number of segments in the skeletonized image
    
    %Set up code for randomized loss
    iter = 1000; %Set number of iterations for randomization
    
    for a = 1:length(loss_range) %For each percent loss in the csv file
        pctloss = loss_range(a); %Percent loss
        iniloss = pctloss - .002; %Creating .2% below the expected loss so the actual loss doesn't overshoot the intended loss
        filename = num2str(pctloss*100); %Creates folder of percent loss
        pctpath = fullfile(subj_path,filename); %Creates a path with the percent loss folder at the end
        mkdir(pctpath); %Makes the percent loss folder
        numsegs_removed = round(numsegs*iniloss); %Find the first bulk of segments that need to be removed 
        segmentarray = {};
        trueloss = [];
        
        for b = 1:iter
            numsegind = randperm(numsegs); %Randomizes the list of segments in the skeletonized image
            tempmat = zeros(skelsize(1), skelsize(2)); %Create a matrix to store all of the segments lost
            for c = 1:numsegs_removed %For the number of segments needed to be removed
                segmentarray{c} = link(numsegind(c)).point; %Creates a cell array of the segment coordinates
                tempmat(segmentarray{c}) = 1; %the coordinates of the segment array in tempmat become 1
            end
            trueloss(b) = sum(tempmat,'all')/skelpx; %Calculates the true loss of segments removed
            while trueloss(b) < pctloss
                c = c + 1;
                segmentarray{c} = link(numsegind(c)).point;
                tempmat(segmentarray{c}) = 1; %the coordinates of the segment array in tempmat become 1
                trueloss(b) = sum(tempmat,'all')/skelpx;
            end
            loss_skel = unmasked_skel - tempmat; %subtracts the segments designated to be lost from the original skeleton image
            loss_dilated = imdilate(loss_skel,SE5); %dilating the skeleton image by the same filter as used by OCT_processing_v3 (line 110)
            loss_with_mask = loss_dilated | mask; %Add in masked arterioles for more accurate image
            losspic = 255*uint8(loss_with_mask); % Figure 2I
            imwrite(losspic,[pctpath '/' subID '_OD_rloss_' num2str(pctloss*100) 'pct_' num2str(b) 'i.bmp']); %\ = PC, / = Mac
        end
    trueloss = trueloss';
    writematrix(trueloss,[pctpath '/' subID '_OD_rloss_' num2str(pctloss*100) 'pct_RealLoss.csv']); %\ = PC, / = Mac
    end
end
toc

            
            
       
