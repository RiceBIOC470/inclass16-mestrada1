% Inclass16

%The folder in this repository contains code implementing a Tracking
%algorithm to match cells (or anything else) between successive frames. 
% It is an implemenation of the algorithm described in this paper: 
%
% Sbalzarini IF, Koumoutsakos P (2005) Feature point tracking and trajectory analysis 
% for video imaging in cell biology. J Struct Biol 151:182?195.
%
%The main function for the code is called MatchFrames.m and it takes three
%arguments: 
% 1. A cell array of data called peaks. Each entry of peaks is data for a
% different time point. Each row in this data should be a different object
% (i.e. a cell) and the columns should be x-coordinate, y-coordinate,
% object area, tracking index, fluorescence intensities (could be multiple
% columns). The tracking index can be initialized to -1 in every row. It will
% be filled in by MatchFrames so that its value gives the row where the
% data on the same cell can be found in the next frame. 
%2. a frame number (frame). The function will fill in the 4th column of the
% array in peaks{frame-1} with the row number of the corresponding cell in
% peaks{frame} as described above.
%3. A single parameter for the matching (L). In the current implementation of the algorithm, 
% the meaning of this parameter is that objects further than L pixels apart will never be matched. 

% Continue working with the nfkb movie you worked with in hw4. 

% Part 1. Use the first 2 frames of the movie. Segment them any way you
% like and fill the peaks cell array as described above so that each of the two cells 
% has 6 column matrix with x,y,area,-1,chan1 intensity, chan 2 intensity

image1 = 'nfkb_movie1.tif';
reader1 = bfGetReader(image1);
iplane1 = reader1.getIndex(0,0,0)+1;
frame1 = bfGetPlane(reader1, iplane1);
iplane2 = reader1.getIndex(0,0,1)+1;
frame2 = bfGetPlane(reader1, iplane2);

imshow(frame1, []);
imshow(frame2, []);

lims = [100 2000];
figure;
subplot(1,2,1); imshow(frame1, lims);
subplot(1,2,2); imshow(frame2, lims);
figure; imshowpair(imadjust(frame1), imadjust(frame2));

imwrite(frame1, 'nfkb_frame1.tif', 'tif');
imwrite(frame2, 'nfkb_frame2.tif', 'tif');

% Performed segmentation using Ilastik

mask1 = readIlastikFile('nfkb_frame1_Simple Segmentation.h5');
mask1 = mask1(:,:,1);
imshow(mask1, []);

mask2 = readIlastikFile('nfkb_frame2_Simple Segmentation.h5');
mask2 = mask2(:,:,1);
imshow(mask2, []);

stats1 =regionprops(mask1, 'Area');
figure; hist([stats1.Area]);
xlabel('Cell Area', 'FontSize', 24);
ylabel('Frequency', 'FontSize', 24);

minarea = 50;
mask1 = imfill(mask1, 'holes');
mask1 = bwareaopen(mask1, minarea);
imshow(mask1, []);

mask2 = imfill(mask2, 'holes');
mask2 = bwareaopen(mask2, minarea);
imshow(mask2, []);

figure; 
subplot(1,2,1); imshow(mask1);
subplot(1,2,2); imshow(mask2);
figure; imshowpair(mask1, mask2);

stats1 = regionprops(mask1, frame1, 'Centroid', 'Area', 'MeanIntensity');
stats2 = regionprops(mask2, frame2, 'Centroid', 'Area', 'MeanIntensity');

iplane1 = reader1.getIndex(0,1,0)+1;
frame1_c2 = bfGetPlane(reader1, iplane1);
iplane2 = reader1.getIndex(0,1,1)+1;
frame2_c2 = bfGetPlane(reader1, iplane2);
imwrite(frame1_c2, 'nfkb_frame1_c2.tif', 'tif');
imwrite(frame2_c2, 'nfkb_frame2_c2.tif', 'tif');

% Performed segmentation using Ilastik

mask1_c2 = readIlastikFile('nfkb_frame1_c2_Simple Segmentation.h5');
mask1_c2 = mask1_c2(:,:,1);
imshow(mask1_c2, []);

mask2_c2 = readIlastikFile('nfkb_frame2_c2_Simple Segmentation.h5');
mask2_c2 = mask2_c2(:,:,1);
imshow(mask2_c2, []);

minarea = 50;
mask1_c2 = imfill(mask1_c2, 'holes');
mask1_c2 = bwareaopen(mask1_c2, minarea);
imshow(mask1_c2, []);

mask2_c2 = imfill(mask2_c2, 'holes');
mask2_c2 = bwareaopen(mask2_c2, minarea);
imshow(mask2_c2, []);

stats1_c2 = regionprops(mask1_c2, frame1_c2, 'Centroid', 'Area', 'MeanIntensity');
stats2_c2 = regionprops(mask2_c2, frame2_c2, 'Centroid', 'Area', 'MeanIntensity');

xy1 = cat(1, stats1.Centroid);
xy2 = cat(1, stats2.Centroid);
a1 = cat(1, stats1.Area);
a2 = cat(1, stats2.Area);
mi1 = cat(1, stats1.MeanIntensity);
mi2 = cat(1, stats2.MeanIntensity);
mi1_c2 = cat(1, stats1_c2.MeanIntensity);
mi2_c2 = cat(1, stats2_c2.MeanIntensity);
tmp1 = -1*ones(size(a1));
tmp2 = 1*ones(size(a2));
peaks{1} = [xy1, a1, tmp1, mi1, mi1_c2];
peaks{2} = [xy2, a2, tmp2, mi2, mi2_c2];

% Part 2. Run match frames on this peaks array. ensure that it has filled
% the entries in peaks as described above. 

match_peaks = MatchFrames(peaks, 2, 50);

% Part 3. Display the image from the second frame. For each cell that was
% matched, plot its position in frame 2 with a blue square, its position in
% frame 1 with a red star, and connect these two with a green line.

imshow(frame2, [200, 1000]); hold on;
for ii = 1:size(peaks{1})
    plot(peaks{1}(ii,1), peaks{1}(ii,2), 'r*', 'MarkerSize', 24)
    ind = match_peaks{1}(ii,4);
    if ind > 0
        plot(peaks{2}(ii,1), peaks{2}(ii,2), 'cs', 'MarkerSize', 24)
        plot([peaks{2}(ii,1) peaks{1}(ii,1)], [peaks{2}(ii,2) peaks{1}(ii,2)], 'g');
    end
end
    
    