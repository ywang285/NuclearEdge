function [outputImage, number_of_nuclei, afterLoG, binarization] = nuclei_counter(sourceImage, threshold)
% This script accept an input nuclei image taken under microscope,
% then preprocess the image with tophat and bottomhat filter, enhance its 
% contrast by multiplying gamma. Convert it to binary image with a histogram 
% based threshold, count the number of connected components and detect edges
% of nucleus with LoG filter, then plot edges. Or calculate distance transform 
% and detect seed points,use watershed algorithm for segmentation.
% 
% flow chart:
% input image ---> tophat/bottomhat filter ---> image enhancement--->binarization by a histogram based threshold
%  ---> calculate seed points ---> segmentation or plot edges of blobs ---> count connected components
%  
% refrence:
% * Yousef Al-Kofahi, et al. Improved Automatic Detection and Segmentation of Cell Nuclei in Histopathology Images. 
%   IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. 57, NO. 4, APRIL 2010
% * Jos B.T.M, et al, The Watershed Transform: Defnitions, Algorithms and Parallelization Strategies
%   Fundamenta Informaticae 41 (2001) 187{228
% 
% Input
%     sourceImage: 8 bits grayscale image;
%     threshold: threshold [0 1] for binarization.
%     
% Output:
%     outputImage: source image with overlayed contour
%     number_of_nuclei: number of nucleis in image
%     
% If you have any ideas or improvements about this method, welcome to contact me 'zhyh8341@gmail.com'.


% display original image
image = sourceImage;
figure;imshow(image,[min(image(:)) max(image(:))]);
title('input image');

total = numel(image);


% apply top hat and bottom hat filter
se = strel('disk',30);
tophat = imtophat(image,se);
bottomhat = imbothat(image,se);
filterImage = image + (tophat - bottomhat);
se = strel('disk',15);
tophat = imtophat(filterImage,se);
bottomhat = imbothat(filterImage,se);
filterImage = filterImage + (tophat - bottomhat);

% calculate histogram of filtered image
% estimate more than 78.5% area is background (pi/4 = .785)
[counts,x] = imhist(filterImage);
ssum = cumsum(counts);
bg = .215*total;
fg = .99*total;
low = find(ssum>bg, 1, 'first');
high = find(ssum>fg, 1, 'first');
adjustedImage = imadjust(filterImage, [low/255 high/255],[0 1],1.8);


% image binarization, threshold is choosen based on experience
if(nargin < 2)
    matrix = reshape(adjustedImage,total,1);
    matrix = sort(matrix);
    threshold = graythresh(matrix(total*.5:end))/3; % adjust default threshold
end
binarization = im2bw(adjustedImage,threshold);


% open image and then detect edge using laplacian of gaussian
se2 = strel('disk',5);
afterOpening = imopen(binarization,se2);
nsize = 5; sigma = 3;
h = fspecial('log',nsize,sigma);
afterLoG = uint8(imfilter(double(afterOpening)*255,h,'same').*(sigma^2)); 


se2 = strel('disk',5);
afterOpening = imopen(binarization,se2);
number_of_nuclei = bwconncomp(afterOpening);


% % you can either use watershed method to do segmentation
% D = -bwdist(~afterOpening);
% D(~afterOpening) = -Inf;
% L = watershed(D);


outputImage = sourceImage + afterLoG*5;
figure;imshow(outputImage, [min(sourceImage(:)) max(sourceImage(:))]);
title('output image');

    