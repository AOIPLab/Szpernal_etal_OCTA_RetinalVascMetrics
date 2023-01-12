% OCTA_Metrics_v4p3_PC.m
% Created on: 2021?
% Created by: Mina Gaffney and Rachel Linderman
%
% This code was written by Mina Gaffney and edited by Rachel Linderman. It
% is an adaptation of the PICA code from NYEEI written by Toco Chui. The
% purpose of this code is to load in a binary image and to calculate
% parafoveal intercapillary areas for the whole image.
%
% Revision Log:
% OCTA_Metrics_v4p2 updated 1/28/2022 by MG 
% Removed pathparts and top_path
% changed top_dir to look for the LUT in the top folder level. 
% 
% Expected input structure as of 1/28/2022:
% Top folder level (to be selected by user via uigetfile)
    % Contains 
        % Info.csv a 3 column matrix. No Header
            % Format: column 1 contains loss itteration ex: _1i 
            % column 2 contains the eye formatted as OD or OS 
            % column 3 contains the RMF formatted as x.xxxx
        % All subject folders labeled by subject ID ex: JC_#####
            % Subject folders contain percent loss folders labeled by
            % number ex: 1,2,10 etc.
                % Percent loss folders contain: itteations of simulated
                % loss .bmp's 
                % Naming convention: SubjectID_EYE_rloss_#pct_#i.bmp
                % ex: JC_#####_OD_rloss_2pct_1i.bmp
                
%OCTA_Metrics_v4p3 updated 1/28/2022 by JC 
%Added hard code RMF, commented out Bennett's forumala
%Modfied for PC settings - ; : \ / pathparts{#}
                
                
close all;
clear all;

CENTER_OF_OCTA_IND = 1662576;
tic
        
    
mfilefullpath = matlab.desktop.editor.getActiveFilename;
mfileshortpath = erase(mfilefullpath, 'OCTA_Metrics_v4p3_PC.m');
natsort_path = append(mfileshortpath, 'natsortfiles');
lib_path = append(mfileshortpath, 'lib');
addpath(genpath(natsort_path));
addpath(genpath(lib_path));

%Batch Upload feature - load in top level of folder structure
root_path = uigetdir('.','Select folder containing images');
%pathparts = strsplit(root_path,'\'); % MG removed 1/28/2022 
root_dir = dir(root_path);
root_dir = root_dir(~ismember({root_dir.name}, {'.', '..'}));

% Ignore Lut or anything that isn't a folder
root_dir = root_dir([root_dir.isdir]);
rootdir_folders = cell(length(root_dir),1); %pre-allocate cell array
addpath(genpath(root_path));
%top_path = fullfile(pathparts{1:3}); %this needs to change as file path
%changes % MG removed 1/28/2022 This would look for the In.csv sheet
%several levels up
top_dir = dir(root_path); % MG changed 1/28/2022 was: top_dir = dir(top_path);
top_dir = top_dir(~ismember({top_dir.name}, {'.', '..'}));
top_dir = struct2cell(top_dir)';

%Loop through first level of folders in directory
for kk = 1:length(root_dir)
    rootdir_folders{kk} = root_dir(kk).name;
    rootdir_paths{kk,1} = fullfile(root_path,rootdir_folders{kk});

    %Finds the full path for all folders and subfolders in selected folder:
    searchlist = genpath(rootdir_paths{kk});
    searchlist = strsplit(searchlist, ';'); %: for MAC and ; for PC
    
    %Read in files/folders from current folder
    Current_subfolder = dir(searchlist{1,1});
    Current_subfolder = Current_subfolder(~ismember({Current_subfolder.name},{'.','..'}));
    Current_subfolder = Current_subfolder([Current_subfolder.isdir]);
    Current_subfolder = struct2cell(Current_subfolder)';

    %Loop through Percent Loss folders
    for kkk = 1:size(Current_subfolder,1)
        Current_lossfolder = dir(fullfile(Current_subfolder{kkk,2},Current_subfolder{kkk,1}));
        Current_lossfolder = Current_lossfolder(~ismember({Current_lossfolder.name}, {'.','..'}));
        Current_lossfolder = struct2cell(Current_lossfolder)';
        Currentloss_pathparts = strsplit(Current_lossfolder{kkk,2},'\'); %/ for MAC and \ for PC
        
        %Find Axial Length Look Up Table
        batch_dir = top_dir(~cellfun(@isempty, strfind(top_dir(:,1),'Info')),:);
        batch_info = table2cell(readtable(fullfile(batch_dir{1,2},batch_dir{1,1}),'Format','auto'));
        
        %Looks for all bmps
        img_dir = Current_lossfolder(~cellfun(@isempty, strfind(Current_lossfolder(:,1),'bmp')),:);
        img_dir = natsortfiles(img_dir(:,1));
        
        comp_table = cell(0,14);
        
        for i = 1:size(img_dir)
            Subject_ID = batch_info{i,1};
            eye = batch_info{i,2};
            %subjectAL = batch_info{i,3};
                        
            % Pull corresponding scan based on subject ID and eye
            subj_tifs = img_dir(...
                ~cellfun(@isempty, strfind(img_dir(:,1), Subject_ID)),:);
            
            octa_img = char(subj_tifs(...
                ~cellfun(@isempty, strfind(subj_tifs(:,1), eye)), 1));
            
            img = logical(imread(octa_img));
            
            %RMF = (3000*((0.01306*(subjectAL-1.82))/(0.01306*(23.95-1.82))))/1824; %retinal magnification factor(Bennett's formula)
            RMF = batch_info{i,3} %hard code for 2380/1824 by Joe - 1/28/22
            
            %Calculating Vessel Density over whole image
            whitepixels = sum(img,'all');
            [numRows, numCols] = size(img);
            total = numRows*numCols;
            VesDen = (whitepixels/total)*100;
            
            %Calculating Vessel Complexity Index (VCI) and Vessel Perimeter
            %Index (VPI)
            perimeter_img = bwperim(img,8); %creating perimeter matrix
            perimeter_pixels = sum(perimeter_img,'all'); %summing all of the white pixels up in the perimeter image
            VCI = ((perimeter_pixels^2)/(4*pi*whitepixels))/1.6169; %calculating VCI
            VPI = perimeter_pixels/total; %calculating VPI
            
            %Calculating fractal dimension
            [n,r] = boxcount(img);
            f = fit(r',n','power1'); %fit r vs n to get the power function
            coeff_values = coeffvalues(f); %pull the coefficients calculated by the fit function
            Df = abs(coeff_values(2)); %pull the the exponent from the coefficient values to get the fractal dimension
            
            %Calculating PICA over whole image
            clear nonperfused Area Centroidx Centroidy PixelIdxList RingResults NoFAZ large_Areas large_Areasinds
            nonperfused = regionprops(~img,'area','Centroid','PixelIdxList'); %get the area, centroid and pixel list in a structure for all of the regions in the image
            for k = 1:length(nonperfused)
                Area(k) = nonperfused(k,1).Area;
                Centroidx(k) = nonperfused(k,1).Centroid(1);
                Centroidy(k) = nonperfused(k,1).Centroid(2);
                PixelList{k} = nonperfused(k,1).PixelIdxList(:);
            end
            
            % Sorts PICA areas in descending size
            [large_Areas, large_Areasinds]  = sort(Area,'descend');%sorted indexes
            
            %Finding the PICA area which contains the central point of the image (912, 912)
            for k = 1:10
               if any(PixelList{large_Areasinds(k)} == CENTER_OF_OCTA_IND)
                   fazhere = k;
                   break
               end
            end
            
            faz_area = Area(large_Areasinds(k));
            
            %Accounting for ocular magnification for FAZ
            FAZ_Area = (faz_area*RMF^2)/1000000;
            
            %Removing all PICA values <51
            [r2, c2] = find(Area<51);
            for ii = length(r2)
                Area(r2(ii),c2(ii)) = 0;
                Centroidx(r2(ii),c2(ii)) = 0;
                Centroidy(r2(ii),c2(ii)) = 0;
            end
            
            %Saving raw data after removing small PICA values
            Area(Area == 0) = NaN;
            Centroidx(Centroidx == 0) = NaN;
            Centroidy(Centroidy == 0) = NaN;
            
            for iii = 1:length(Area)
                RingResults(iii,:) = [Area(iii),Centroidx(iii),Centroidy(iii)];
                RingResults(iii,1) = RingResults(iii,1)*(RMF^2)/1000000;
            end
            
            averageRingResults = mean(RingResults,'omitnan');
            SDRingResults = std(RingResults,'omitnan');
            sumPICA = sum(RingResults,'omitnan');
            PICAReg = averageRingResults(1)/SDRingResults(1);
            
            %Removing FAZ from PICA measurements
            Area(large_Areasinds(k)) = NaN;
            Centroidx(large_Areasinds(k)) = NaN;
            Centroidy(large_Areasinds(k)) = NaN;
            
            for iii = 1:length(Area)
                NoFAZ(iii,:) = [Area(iii),Centroidx(iii),Centroidy(iii)];
                NoFAZ(iii,1) = NoFAZ(iii,1)*(RMF^2)/1000000;
            end
            
            NoFAZAverage = mean(NoFAZ,'omitnan');
            NoFAZSD = std(NoFAZ,'omitnan');
            NoFAZsumPICA = sum(NoFAZ,'omitnan');
            NoFAZPICAReg = NoFAZAverage(1)/NoFAZSD(1);
            
            %Compiling all of the different metrics into 1 matrix
            output(:,1) = {octa_img};
            output(:,2) = {FAZ_Area}; 
            output(:,3) = {VesDen};
            output(:,4) = {averageRingResults(1)};
            output(:,5) = {sumPICA(1)};
            output(:,6) = {PICAReg(1)};
            output(:,7) = {SDRingResults(1)};
            output(:,8) = {NoFAZAverage(1)};
            output(:,9) = {NoFAZsumPICA(1)};
            output(:,10) = {NoFAZPICAReg};
            output(:,11) = {NoFAZSD(1)};
            output(:,12) = {VCI};
            output(:,13) = {VPI};
            output(:,14) = {Df};
            
            comp_table = [comp_table;output];
            newname2 = strcat(octa_img(1:end-4),'.csv');
            writematrix(RingResults,fullfile(Current_lossfolder{i,2},newname2));
            
        end
        header1 = {'Image','FAZ Area','Vessel Density','Average PICA Area','Summed PICA Area','PICA Regularity',...
            'PICA SD','Average PICA Area','Summed PICA Area','PICA Regularity','PICA SD','VCI','VPI','DF'};
        finaloutput = [header1;comp_table];
        newname = strcat(Currentloss_pathparts{4},'_',Currentloss_pathparts{5},'pct_Raw_Data.csv'); %Need to change {#} based on folder structure (Mac / PC) 4, 5 for PC and 11, 12 for Mac
        writecell(finaloutput,fullfile(Current_lossfolder{1,2},newname));
        
        clear output  
    end
    clear comp_table finaloutput
end
toc
disp('Batch Complete!')
