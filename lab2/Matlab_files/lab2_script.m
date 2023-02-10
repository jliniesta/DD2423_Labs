%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          LAB 2 - SCRIPT 
%
% Javier Lopez Iniesta Diaz del Campo
%         Mathias Näreaho
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close('all')

%% Difference operators

% Load image
tools = few256;

% Difference operators
deltax = [1 0 -1; 2 0 -2; 1 0 -1];
deltay = deltax';

% Discrete derivation approximations
dxtools = conv2(tools, deltax, 'valid');
dytools = conv2(tools, deltay, 'valid');

figure('Name','Difference operators','NumberTitle','off');
subplot(1, 3, 1);
    showgrey(tools);
    title('Original image')
subplot(1, 3, 2);
    showgrey(dxtools);
    title('deltax')
subplot(1, 3, 3);
    showgrey(dytools);
    title('deltay')
 
%% Point-wise thresholding of gradient magnitudes

img = godthem256;

figure('Name','Original image','NumberTitle','off');
showgrey(img);
title('Original image')

% 1) Whithout blurred the image

img_1 = Lv(img);

% Threshold
threshold1 = [50, 100, 130, 200, 300];

figure('Name','Point-wise thresholding of gradient magnitudes','NumberTitle','off');
subplot(2, 6, 1);
    showgrey(img_1);
    title('Original image')

for i=1:length(threshold1)
subplot(2, 6, i+1);
    showgrey((img_1 - threshold1(i)) > 0)
    title(sprintf('Threshold = %d', threshold1(i)))  
end

% 2) Blurring the image
t = 2;
img_filtered = gaussfft(img, t);
img_2 = Lv(img_filtered);

% Threshold
threshold2 = [50, 88, 100, 130, 200];

subplot(2, 6, 7);
    showgrey(img_2);
    title('Blurred image')

for i=1:length(threshold2)
subplot(2, 6, 6+i+1);
    showgrey((img_2 - threshold2(i)) > 0)
    title(sprintf('Threshold = %d', threshold2(i)))  
end

figure('Name','Blurred image','NumberTitle','off');
showgrey(img_filtered);
title('Blurred image')

% figure('Name','Histogram','NumberTitle','off');
% hist(img_1)
% title('Image histogram')
% 
% figure('Name','Histogram after smoothing','NumberTitle','off');
% hist(img_2)
% title('Image histogram after smoothing')

figure('Name','Histograms','NumberTitle','off');
subplot(2, 1, 1);
    hist(img_1)
    title('Image histogram')
subplot(2, 1, 2);
    hist(img_2)
    title('Image histogram after smoothing')

%% Computing differential geometry descriptors

close('all')
clc; 
clear; 

% TEST:
% [x y] = meshgrid(-5:5, -5:5);
% 
% delta_x = [0,   0, 0,    0, 0; 
%            0,   0, 0,    0, 0; 
%            0, 1/2, 0, -1/2, 0;
%            0,   0, 0,    0, 0;
%            0,   0, 0,    0, 0];
% 
% delta_y = delta_x';
% delta_xx = [0, 0,  0, 0, 0; 
%             0, 0,  0, 0, 0; 
%             0, 1, -2, 1, 0;
%             0, 0,  0, 0, 0;
%             0, 0,  0, 0, 0];
% 
% delta_xxx = conv2(delta_xx, delta_x, 'same');
% delta_xxy = conv2(delta_xx, delta_y, 'same');
% 
% test1 = conv2(delta_xxx, x.^3, 'same')
% test2 = conv2(delta_xx,  x.^3, 'same')
% test3 = conv2(delta_xxy, (x.^2).*y, 'same')

house = godthem256;
tools = few256;

scale = [0.0001, 1, 4, 16, 64];

% 0 crossings of Lvv_tilde with different scales

figure('Name','0 crossings of Lvv_tilde with different scales','NumberTitle','off');
subplot(1, 6, 1);
    showgrey(house);
    title('Original image')

for i=1:length(scale)
subplot(1, 6, i+1);
    contour(Lvvtilde(discgaussfft(house, scale(i)), 'same'), [0 0])
    axis('image')
    axis('ij')
    if (i==1)
        title(sprintf('Scale = %1.4f', scale(i)))
    else 
        title(sprintf('Scale = %1.0f', scale(i)))
    end
end

% Sign of Lvvv_tilde with different scales

figure('Name','Sign of Lvvv_tilde with different scales','NumberTitle','off');
subplot(1, 6, 1);
    showgrey(tools);
    title('Original image')

for i=1:length(scale)
subplot(1, 6, i+1);
    showgrey(Lvvvtilde(discgaussfft(tools, scale(i)), 'same') < 0)
    axis('image')
    axis('ij')
    if (i==1)
        title(sprintf('Scale = %1.4f', scale(i)))
    else 
        title(sprintf('Scale = %1.0f', scale(i)))
    end
end

%% Extraction of edge segments

close('all')
clc; 
clear; 

scale     = [0.0001, 1, 4, 16, 64];
threshold = [30,40,50,60,70];

tools = few256;

% figure('Name','Trying different scales for tools image','NumberTitle','off');
% for i=1:length(scale)
% subplot(1, 6, i+1);
%     overlaycurves(tools,extractedge(tools,scale(i),16,'same'))
%     title(sprintf('Scale = %1.4f', scale(i)))  
% end

% figure('Name','Trying different threshold levels for tools image','NumberTitle','off');
% for i=1:length(threshold)
% subplot(1, 6, i+1);
%     overlaycurves(tools,extractedge(tools,4,threshold(i),'same'))
%     title(sprintf('Threshold = %1.0f', threshold(i)))  
% end

house = godthem256;

% figure('Name','Trying different scales for house image','NumberTitle','off');
% for i=1:length(scale)
% subplot(1, 6, i+1);
%     overlaycurves(tools,extractedge(tools,scale(i),16,'same'))
%     title(sprintf('Scale = %1.4f', scale(i)))  
% end

% figure('Name','Trying different threshold levels for house image','NumberTitle','off');
% for i=1:length(threshold)
% subplot(1, 6, i+1);
%     overlaycurves(tools,extractedge(tools,4,threshold(i),'same'))
%     title(sprintf('Threshold = %1.0f', threshold(i)))  
% end

% Best results
figure('Name','Best results','NumberTitle','off');
subplot(1, 2, 1);
    overlaycurves(tools,extractedge(tools,4,50,'same'))
    title(sprintf('Scale = 4; Threshold = 50'))  
subplot(1, 2, 2);
    overlaycurves(house,extractedge(house,4,40,'same'))
    title(sprintf('Scale = 4; Threshold = 40'))     


%% Hough transform

% image = few256;
image = godthem256;
% image = houghtest256;
% image = triangle128;
% image = godthem256;

% image_subsample = binsubsample(image);

scale = 4;
gradmagnthreshold = 60;
nrho = 256;
ntheta = 180;
nlines = 10;
verbose = 0;
     
[linepar, acc] = houghedgeline(image,scale,gradmagnthreshold,nrho,ntheta,nlines,verbose);

