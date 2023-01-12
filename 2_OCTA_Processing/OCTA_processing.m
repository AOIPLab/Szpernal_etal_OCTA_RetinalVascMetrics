% OCTA_processing.m
% Adapted from: test_OCTA_foveal_mappingV3_MCW_10102017 by NYEEI - Toco Chui 08/16/2016
% Created on: 2/25/2021
% By: Mina Gaffney 
% 
% Decription: 
%
% Inputs: 
% A 304x304x3 uint8 averaged OCT
%
% Revision log: 
% 

%% Options
flip = 0; % flip OCTA vertically? (1 for true, 0 for false)

%% Constants 
SE2 = strel('disk',2);
SE3 = strel('disk',3);
SE5 = strel('disk',5);
SE6 = strel('disk', 6);

%% Load and read in OCTA image. Expected input 304x304x3 uint8 
[filename, pathname] = uigetfile({'*.tif';'*.bmp';'*.png'}, 'Select a Superficial Layer OCTA Image 304*304 pixels','MultiSelect','on');
cd(pathname);

%% Check if multiple files have been selected: 

if iscell(filename) == 0
    numImg = 1; 
else
    numImg = length(filename); 
end 

averaged = imread(fullfile(pathname,filename));

%% Contrast stretch image 
stretched_average = imadjust(averaged, stretchlim(averaged),[]);

% Colapses image down to 1 layer: 
[sa_x, sa_y, sa_z]=size(stretched_average);
if sa_z>1
    stretched_average = stretched_average(:,:,1);
end

% Flips image if necessary:
if flip==1;
   stretched_average=flipud(stretched_average);
end


%% Background subtraction 
background = imopen(stretched_average,strel('disk',15));
stretched_average = stretched_average - background;

%% Skeletonizing and Thresholding 

% contrast stretch stretched_average again.
sa_stretched_average = imadjust((stretched_average), stretchlim((stretched_average)),[]);

% Upsample to 6X (304x304 to 1824x1824)
sa_stretched_average_6X = imresize(sa_stretched_average,6);

% Thresholding 
sa_sa_6X_thresh = sa_stretched_average_6X;
sa_sa_6X_thresh(sa_sa_6X_thresh < 20) = 0; % drops any pixel less than 20 to 0
sa_sa_6X_thresh(sa_sa_6X_thresh > 35) = 255; % caps out any pixel more than 35 to 255

% Averaging 
sa_sa_6X_thresh = filter2(fspecial('average',10),sa_sa_6X_thresh); % 10 x 10 pix averaging filter

% Separating the background and foreground, finds perimeter of objects, and
% inverts the resulting matrix
separated_sa_sa_6X = adaptivethreshold(sa_sa_6X_thresh,15 ,0.0000000005,0); 
inverted_separated = medfilt2(~(separated_sa_sa_6X - bwperim(separated_sa_sa_6X)),[10 10]);

% Skeletonizing 
skeleton = bwmorph(~inverted_separated,'skel',inf);
skeleton = bwareaopen(skeleton, 10);
skeleton = bwmorph(skeleton,'spur',10);
skeleton = bwmorph(skeleton,'clean');

% Thresholding skeleton 
testmask = sa_sa_6X_thresh < 30; % creates a mask that is equal to 1 whereever sa_stretched_average_6X is less than 30
skeleton = immultiply(~testmask, skeleton); % thresholds skeleton to remove values less than 30

% Separating the background and foreground again
separated_sa_sa_6X_2 = adaptivethreshold(sa_stretched_average_6X,30 ,0.0000000005,0);
separated_sa_sa_6X_2 = bwareaopen(separated_sa_sa_6X_2, 3000);

% Thresholding OCTA
testmask2 = sa_stretched_average_6X < 30; 
thresholded = immultiply(~testmask2,separated_sa_sa_6X_2);

% Creating a dilated skeleton with an r=3 disk
dilated_skel_3 = imdilate(skeleton,SE3);

% Optional display of OCTA with skeleton overlaid on top. Used for
% visualization purposes 
skeletonoverlay = imoverlay(sa_stretched_average_6X, dilated_skel_3, [1 1 0]);
%figure; imshow(skeletonoverlay);

% Creating a skeleton that has been dilated by a r=5 disk 
dilated_skel_5 = imdilate(skeleton,SE5);

% Combining the r=5 dilated skeleton with the thresholded OCTA 
thresholded2 = dilated_skel_5==1 | thresholded==1;

figure; imagesc(imoverlay(sa_stretched_average_6X, imerode(thresholded2,SE2),[0 1 1])); 
colormap gray; axis image; axis on; hold on; truesize
title('Remove Unwanted Pixels - Drag and circle unwanted pixels. Esc to repeat this step. To quit, one click on the FAZ and pres ESC');
xlim([300 1500]);ylim([300 1500]);
h = imfreehand( gca ); setColor(h,'red');
position = wait( h );
BWW = createMask( h );
close(gcf);

cropped_thresholded = logical(thresholded2.*(~BWW));


%% Masking the larger arterioles

bld_vel_mask = imbinarize(sa_stretched_average_6X, 0.7);
bld_vel_mask = bwareaopen(bld_vel_mask, 400);
bld_vel_mask = imdilate(bld_vel_mask, SE6);

RingMask = ones(1824);
RingMask1 = RingMask - bld_vel_mask;

vascAnalyzed_masked = immultiply(RingMask1, cropped_thresholded);
vascAnalyzed_unmasked = immultiply(RingMask, cropped_thresholded);

origvascAnalyzed = 255.*(imresize(vascAnalyzed_masked,1/6,'nearest')>0.5);
origvascAnalyzed2 = 255.*(imresize(vascAnalyzed_unmasked,1/6,'nearest')>0.5);

maskedskelimg = bwmorph(origvascAnalyzed,'skel',Inf);

%% Optional: Display figures 

figure; imshow(cropped_thresholded);
figure; imshow(origvascAnalyzed);
figure; imshow(origvascAnalyzed2);
figure; imshow(maskedskelimg); 

%% Writing outputs to file 

% appedning things to the output filenames
outfilename_croppedthresh = strrep(filename, '.tif', '_CroppedThresholded.tif');
outfilename_masked = strrep(filename, '.tif', '_VascMasked.tif');
outfilename_unmasked = strrep(filename, '.tif', '_Unmasked.tif');
outfilename_maskedskel = strrep(filename, '.tif', '_VascMaskedSkel.tif');

imwrite(cropped_thresholded,fullfile(pathname,outfilename_croppedthresh));
imwrite(origvascAnalyzed,fullfile(pathname,outfilename_masked));
imwrite(origvascAnalyzed2,fullfile(pathname,outfilename_unmasked));
imwrite(maskedskelimg,fullfile(pathname,outfilename_maskedskel));

