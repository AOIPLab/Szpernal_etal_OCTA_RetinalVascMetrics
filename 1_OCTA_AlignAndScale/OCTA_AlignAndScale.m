clear all
close all

%Purpose of this code is to align multiple OCTA images to their FAZ
%centers, resize the images to the same scale and crop to a common area.
%The scale is determined by the image with the highest resolution in the
%group (The subject with the shortest axial length will have the highest
%resolution for the image.
%Code was written and verified by R. Linderman - 03/03/2021


%Select directory where images and scaling file are
img_dir = uigetdir(pwd,'Select the folder containing the images');

fNames = read_folder_contents(img_dir,'tif');

[lutf, lutp] = uigetfile(fullfile(img_dir,'*.csv'), 'Select scaling (csv) file');

%Load in data
fid = fopen(fullfile(lutp, lutf), 'r');

lutData = textscan(fid, '%s%s%f%f%f','delimiter', ',');

fclose(fid);

%Find the smallest axial length
AL_min = min(lutData{3});

imgs = cell(length(fNames),1);

num_overlap = length(fNames); %This determines the amount of overlap

for i = 1:length(fNames)
    imgs{i} = im2double(imread(fullfile(img_dir,fNames{i})));
    img_dim = size(imgs{i});
    
    
    [idpiece1, remain] = strtok(fNames{i},'_');
    [idpiece2, ~] = strtok(remain,'_');
    subID = [idpiece1 '_' idpiece2];
    
    LUTindex = find(strcmp(lutData{1},subID));
    
    axialLength = lutData{3}(LUTindex);
    foveal_x = lutData{4}(LUTindex);
    foveal_y = lutData{5}(LUTindex);
    
    linearscale = round((axialLength/AL_min),3);
    real_imgs{i} = imresize(imgs{i},linearscale);
    real_foveal_coord(i,:) = round([foveal_x.*linearscale foveal_y.*linearscale],3);
    real_dim(i,:) = size(real_imgs{i});

    
    dist(1) = real_foveal_coord(i,1);
    dist(2) = real_foveal_coord(i,2);
    dist(3) = real_dim(i,1) - dist(1);
    dist(4) = real_dim(i,2) - dist(2);
    minimum(i) = min(dist);
    clear remain idpiece1 idpiece2;
end

true_min = min(minimum);

for t = 1:length(real_imgs)
    c_min = round(real_foveal_coord(t,1) - true_min + 1);
    c_max = round(real_foveal_coord(t,1) + true_min);

    r_min = round(real_foveal_coord(t,2) - true_min + 1);
    r_max = round(real_foveal_coord(t,2) + true_min);

    crop_img = real_imgs{t}(r_min:r_max,c_min:c_max);
    finalpx = num2str(size(crop_img,1));
    final_img = imresize(crop_img,[304,304]);
    newname = strcat(fNames{t}(1:end-11),finalpx,'_aligned.tif');
    imwrite(final_img,fullfile(lutp,newname));
end
    