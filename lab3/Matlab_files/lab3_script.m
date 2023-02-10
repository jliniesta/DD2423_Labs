%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          LAB 3 - SCRIPT 
%
% Javier Lopez Iniesta Diaz del Campo
%         Mathias Näreaho
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close('all')

%% K-means clustering

% Parameters
K = 2;               % number of clusters used
L = 200;              % number of iterations
seed = 14;           % seed used for random initialization
scale_factor = 0.7;  % image downscale factor
image_sigma = 1.0;   % image preblurring scale
image = 'orange';
% image = 'tiger1';
% image = 'tiger2';
% image = 'tiger3';

I = imread(strcat(image,'.jpg'));
I = imresize(I, scale_factor);
Iback = I;
d = 2*ceil(image_sigma*2) + 1;
h = fspecial('gaussian', [d d], image_sigma);
I = imfilter(I, h);

tic
[segm, centers] = kmeans_segm(I, K, L, seed);
toc
Inew = mean_segments(Iback, segm);
I = overlay_bounds(Iback, segm);
imwrite(Inew,strcat('../result/', image,'_k_', num2str(K), '_1.png'))
imwrite(I,strcat('../result/', image,'_k_', num2str(K), '_2.png'))

figure('Name','K-means clustering','NumberTitle','off');
% imshow(Inew);
imshow(Inew);
title(sprintf('K = %d', K), FontSize=16)  
