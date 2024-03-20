% imOverlay
%
% code to overlay one image onto another
close all
clear variables
 
% Base image for displaying on
Im = imread('mri.tif');
subplot(221)
imshow(Im)
 
% Overlay color image
Im_color = reshape(cmaperise(Im),[size(Im),3]);
subplot(222)
imshow(Im_color)
 
% Extract ROI from color
% (I'm doing this manually, but you can get this using imageSegmenter)
% (imageSegmenter has an export for the mask after segmenting)
BW = false(size(Im));
BW(60:80,40:60) = true;
subplot(223)
imshow(BW)
 
% Final plot
subplot(224)
imshow(Im)
hold on
h_overlay = imshow(Im_color);
set(h_overlay,'AlphaData',BW)
